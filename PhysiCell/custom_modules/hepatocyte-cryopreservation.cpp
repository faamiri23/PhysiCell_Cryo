/*
Program name       : hepatocyte-cryopreservation.cpp                                                           
Date Written       : Septemeper 5, 2021                                             
Date Last Modified : June 8, 2023                                             
Authors            : Fatemeh Amiri                                              
Discription        : This file uses rat hepatocyte parameters and                    
                     implement and solve the Monte-Carlo (Gillespie) ice 
					 formation and propagtion model in small tissue using 
					 PhysiCell (Version 1.7.1) [1], with BioFVM [2] to solve 
					 the transport equations.                                                                                                       
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
*/

#include "./hepatocyte-cryopreservation.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <cstdlib>
// Global Variables
#define PI 3.14159265
#define sgn(v) ( ( (v) < 0 ) ? -1 : ( (v) > 0 ) )


Cell_Definition hepatocyte_cell;

std::ofstream ofs;
/////////////////////////////___Custom Cell Lists___///////////////////////////

std::vector<Cell*> neighbors;
std::vector<Cell*> colored_cell_list;

///////////////////////////////////////////////////////////////////////////////

//Create hepatocyte cells type and set the required parameters and variables
void create_hepatocyte_cell_type(void)
{
	hepatocyte_cell = cell_defaults;

	hepatocyte_cell.name = "hepatocyte cell";
	hepatocyte_cell.type = 4;

	// Turn off proliferation;
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live );
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live );

	hepatocyte_cell.phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = 0.0;

	int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model );

	// Turn off secretion
	hepatocyte_cell.phenotype.secretion.uptake_rates[0] *=0;
	hepatocyte_cell.phenotype.secretion.secretion_rates[0] *=0;
	hepatocyte_cell.phenotype.secretion.uptake_rates[1] *=0;
	hepatocyte_cell.phenotype.secretion.secretion_rates[1] *=0;
	hepatocyte_cell.phenotype.secretion.uptake_rates[2] *=0;
	hepatocyte_cell.phenotype.secretion.secretion_rates[2] *=0;
	
	// set apoptosis to survive 10 days (on average)
	hepatocyte_cell.phenotype.death.rates[apoptosis_index] = 0;


	// turn on motility;
	hepatocyte_cell.phenotype.motility.is_motile = false;

    // Set the mechanical parameters
	hepatocyte_cell.phenotype.mechanics.cell_cell_adhesion_strength *= 75;
	hepatocyte_cell.phenotype.mechanics.cell_cell_repulsion_strength *= 75;
    hepatocyte_cell.phenotype.mechanics.relative_maximum_adhesion_distance *= 1;

	// set functions
	// Set the update functions for velocity, Phenotype and custom cell rule
    hepatocyte_cell.functions.update_velocity = standard_update_cell_velocity;

	hepatocyte_cell.functions.update_phenotype = hepatocyte_phenotype_rule;//phenotype_dt
	hepatocyte_cell.functions.custom_cell_rule = hepatocyte_cell_rule;//mechanics_dt 
	
	
	// set custom data values
	hepatocyte_cell.custom_data["elastic coefficient"]=.0000;
	hepatocyte_cell.custom_data[ "initial_volume" ] = 4.9889e+03;
	hepatocyte_cell.custom_data["Lpg"]= 0.1520;
	hepatocyte_cell.custom_data["Elp"]= 341;//kJ/mol
    
	hepatocyte_cell.custom_data["Ps"]= 0.1167;//um/sec
	hepatocyte_cell.custom_data["Lp_salt"]=0.00135;
	hepatocyte_cell.custom_data["Ps_salt"]=0;
    hepatocyte_cell.custom_data["Csalt"]=0.15*1e-15;//mol/kg=mol/L=fmol/mm^3
    
    hepatocyte_cell.custom_data["Tseed"]= 273.15-1.0;//mol/kg=mol/L=fmol/mm^3
    hepatocyte_cell.custom_data["T_f0"]= 273.15-0.52;
    
    hepatocyte_cell.custom_data["omega"]=110e-4;
    hepatocyte_cell.custom_data["kapa"]=14.0e8;

	hepatocyte_cell.custom_data["Area"]=1.4120e+03;
	hepatocyte_cell.custom_data["R"]= 0.0821;
    hepatocyte_cell.custom_data.add_variable("Temperature","unitless" ,271.0);
    hepatocyte_cell.custom_data.add_variable("tau_t","unitless" ,0.0); //dimensionless time
	hepatocyte_cell.custom_data["isosmotic_volume"]=4.9889e+03; 
	hepatocyte_cell.custom_data["Partial_molar_volume"]=0.0368;
	hepatocyte_cell.custom_data["Vb_fraction"]=0.51;

    hepatocyte_cell.custom_data.add_variable("Prev_active_cell_volume", "unitless", 2.4446e+03);
    hepatocyte_cell.custom_data.add_variable("Active_cell_volume", "unitless",2.4446e+03);

	hepatocyte_cell.custom_data.add_variable("Water_volume", "unitless", 2.4446e+03);
	hepatocyte_cell.custom_data.add_variable("Prev_Water_volume", "unitless", 2.4446e+03);
    hepatocyte_cell.custom_data.add_variable("Prev_total_volume","unitless",4.9889e+03);
	
	hepatocyte_cell.custom_data.add_variable("dVw_dt_1", "unitless", 0.0);
    hepatocyte_cell.custom_data.add_variable("dVcs_dt_1", "unitless", 0.0);
	hepatocyte_cell.custom_data.add_variable("dS_dt_1", "unitless", 0.0);

    //Monte-Carlo methods variable
    hepatocyte_cell.custom_data.add_variable("fvecold", "unitless", 0); 
    hepatocyte_cell.custom_data.add_variable("fvec", "unitless", 0); 
    hepatocyte_cell.custom_data.add_variable("knew", "unitless", 0);
    hepatocyte_cell.custom_data.add_variable("Freezing_time", "unitless", 0.0);

	return;
}
void create_cell_types( void )
{
	// use the same random seed so that future experiments have the
	// same initial histogram of oncoprotein, even if threading means
	// that future division and other events are still not identical
	// for all runs
	SeedRandom( parameters.ints("random_seed") );

	// housekeeping

	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );

	// turn the default cycle model to live,
	// so it's easier to turn off proliferation

	cell_defaults.phenotype.cycle.sync_to_cycle_model( live );

	// Make sure we're ready for 2D

	cell_defaults.functions.set_orientation = up_orientation;
	
	cell_defaults.phenotype.motility.restrict_to_2D = false; // 

	// set to no motility for hepatocyte cells
	cell_defaults.phenotype.motility.is_motile = false;

	// use default proliferation and death
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live );
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live );

	int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model );


    // Turn off the default PhysiCell functions that are not needed here
	// set default uptake and secretion
	// oxygen
	cell_defaults.phenotype.secretion.secretion_rates[0] = 0;
	cell_defaults.phenotype.secretion.uptake_rates[0] = 0;
	cell_defaults.phenotype.secretion.saturation_densities[0] = 0;

	// immunostimulatory
	cell_defaults.phenotype.secretion.saturation_densities[1] = 0;

	cell_defaults.functions.update_phenotype = NULL;

	cell_defaults.functions.custom_cell_rule = NULL;

    // Cancer cells are type default 0
	cell_defaults.name = "cancer cell";
	cell_defaults.type = 0;

	// add custom data


		cell_defaults.custom_data.add_variable("elastic coefficient","unitless",0);
		cell_defaults.custom_data.add_variable("initial_volume","unitless",0);
		cell_defaults.custom_data.add_variable("Lp","unitless",0);
		cell_defaults.custom_data.add_variable("Ps","unitless",0);
        
        cell_defaults.custom_data.add_variable("Lpg","unitless",0);
		cell_defaults.custom_data.add_variable("Elp","unitless",0);
        cell_defaults.custom_data.add_variable("Csalt","unitless",0);
        cell_defaults.custom_data.add_variable("Tseed","unitless",0);//mol/kg=mol/L=fmol/mm^3
        cell_defaults.custom_data.add_variable("T_f0","unitless",0);
        cell_defaults.custom_data.add_variable("omega","unitless",0);//13e-4;//5.9e-4;
        cell_defaults.custom_data.add_variable("kapa","unitless",0);//4.8e9;//3.0e9;

		cell_defaults.custom_data.add_variable("Lp_salt","unitless",0);
		cell_defaults.custom_data.add_variable("Ps_salt","unitless",0);

		cell_defaults.custom_data.add_variable("Area","unitless",0);
		cell_defaults.custom_data.add_variable("R","unitless",0);

		cell_defaults.custom_data.add_variable("isosmotic_volume","unitless",0);

		cell_defaults.custom_data.add_variable("Partial_molar_volume","unitless",0);

		cell_defaults.custom_data.add_variable("Vb_fraction","unitless",0);
		cell_defaults.custom_data.add_variable( "kill rate" , "1/min" , 0 ); // how often it tries to kill
		cell_defaults.custom_data.add_variable( "attachment lifetime" , "min" , 0 ); // how long it can stay attached
		cell_defaults.custom_data.add_variable( "attachment rate" , "1/min" ,0 ); // how long it wants to wander before attaching
		cell_defaults.custom_data.add_variable("dVw_dt_1", "unitless", 0.0);
        cell_defaults.custom_data.add_variable("dVcs_dt_1", "unitless", 0.0);
		cell_defaults.custom_data.add_variable("dS_dt_1", "unitless", 0.0);
        cell_defaults.custom_data.add_variable("Temperature", "unitless", 0.0);
        cell_defaults.custom_data.add_variable("tau_t", "unitless", 0.0);
		cell_defaults.custom_data.add_variable("Prev_total_volume","unitless",0);
		cell_defaults.custom_data.add_variable("test_current_voxel","unitless",0);
        
        cell_defaults.custom_data.add_variable("fvecold", "unitless", 0);
        cell_defaults.custom_data.add_variable("fvec", "unitless", 0);
        cell_defaults.custom_data.add_variable("knew", "unitless", 0);
        cell_defaults.custom_data.add_variable("Freezing_time", "unitless", 0.0);// time that the cell freezes
	// create the cell types
           create_hepatocyte_cell_type();
	return;
}

void setup_microenvironment( void )
{
	// set domain parameters


	if( default_microenvironment_options.simulate_2D == true )
	{
		std::cout << "Warning: overriding 2D setting to return to 3D" << std::endl;
		default_microenvironment_options.simulate_2D = false;
	}

	// gradients are needed for this example

	default_microenvironment_options.calculate_gradients = true;

	// add the immunostimulatory factor

	microenvironment.add_density( "blank", "dimensionless" );
	microenvironment.diffusion_coefficients[1] = 0;
	microenvironment.decay_rates[1] = 0;

	microenvironment.add_density( "EG", "dimensionless" );
	microenvironment.diffusion_coefficients[2] = 1100;
	microenvironment.decay_rates[2] = 0;

	microenvironment.add_density( "PBS", "dimensionless" );
	microenvironment.diffusion_coefficients[3] = 1075;
	microenvironment.decay_rates[3] = 0;
	// let BioFVM use oxygen as the default

	default_microenvironment_options.use_oxygen_as_first_field = true;

	// set Dirichlet conditions
	default_microenvironment_options.outer_Dirichlet_conditions = true;
	default_microenvironment_options.Dirichlet_condition_vector[0] = 0;
	default_microenvironment_options.Dirichlet_condition_vector[1] = 0;
	default_microenvironment_options.Dirichlet_condition_vector[2] = 0;


	initialize_microenvironment();

	return;
}

// Create cells hepatocyte positions in the tissue
std::vector<std::vector<double>> create_cells_hepatocyte_positions(
double cell_radius, double sphere_radius, double inner_radius)
{
	std::vector<std::vector<double>> cells;
	int xc=0,yc=0,zc=0;
    double x_spacing= cell_radius*sqrt(3);
	double y_spacing= cell_radius*2;

    int n=ceil(sphere_radius/cell_radius);

	std::vector<double> tempPoint(3,0.0);

tempPoint[0]=0;
tempPoint[1]=0;
tempPoint[2]=0;
if(sqrt(norm_squared(tempPoint))< sphere_radius)
	{
	cells.push_back(tempPoint);
	}

  double teta_0=0.0;

  double dteta=35.0;double teta=0.0;
  double r=0.0;
  for(int j=1;j<3;j++)
	{
    r += 2*cell_radius;
    int k=0, i=0;
    if (j==1)
    {   while (i<7)
        {
            i++;
            teta=teta_0+(i)*dteta; 
            tempPoint[0]= r*cos(teta);
			tempPoint[1]=r*sin(teta);
            if(sqrt(norm_squared(tempPoint))< sphere_radius)
			{cells.push_back(tempPoint);}
        }
    }
    else
    {
     while (i<14)
     {
        i++;
        teta = teta_0+i*0.5*dteta;   
        tempPoint[0]= r*cos(teta);
		tempPoint[1]=r*sin(teta);
        if(sqrt(norm_squared(tempPoint))< sphere_radius)
		{
            ofs.open ("output/temp_points.csv", std::
            ofstream::out | std::ofstream::app);
			ofs <<tempPoint[0]<<","<<tempPoint[1]<<", "
            <<tempPoint[2]<<std::endl;
			ofs.close();
			cells.push_back(tempPoint);

		}
     }

    }

  }

return cells;

}

void setup_tissue( void )
{
   
	double hepatocyte_radius = 10.6;//rat hepatocyte 

	double sphere_radius =50;
	
	double inner_radius =10.0;



	std::cout << "\tPlacing rat hepatocyte cells ... " << std::endl;


   Cell* pCell;
   

std::vector<std::vector<double>> hepatocyte_positions = create_cells_hepatocyte_positions(hepatocyte_radius, sphere_radius,inner_radius);//disc, monolayer
 
for( int i=0; i <hepatocyte_positions.size(); i++ )
   {
       pCell = create_cell(hepatocyte_cell);
          
       pCell->assign_position(hepatocyte_positions[i]); //hepatocyte position

       pCell->set_total_volume(4.9889e+03);//rat hepatocyte
   }

return;
}


std::vector<std::string> follicle_coloring_function( Cell* pCell )
{
	//static int oncoprotein_i = pCell->custom_data.find_variable_index( "oncoprotein" );

	//Set the color to be red for unfrozen cells
	std::vector< std::string > output( 4, "red" );

    // If the cells are frozen change the color to blue
	if( pCell->state.ice_state==1 )
	{
		output[0] = "blue"; // 
		output[1] = "blue"; // 
		output[2] = "blue"; // 
		
		colored_cell_list.push_back(pCell);
		return output;
	}

	return output;
}

//-------------------------------------------------------------------------------------
/* Check for the water volume of each cell (agent), check if the cell is frozen,
update the cell position and volume based on the mount of water lost and freezing time.*/
//-------------------------------------------------------------------------------------
void hepatocyte_cell_rule( Cell* pCell, Phenotype& phenotype, double dt )
{

	std::ofstream ofs;
  //If the cell is not frozent yet
 if ( pCell->state.ice_state == false)
 { 
     // Check the change of water amount of the cell
     WaterTransportModel_NoCPA(pCell,phenotype,dt);
	 // Compute the dimensionless time
     dimensionless_time_NoCPA(pCell,phenotype,dt);
	 // Write some output data 
     ofs.open ("output/rat_hepatocyte_22cells_B400.txt" , std::ofstream::out | std::ofstream::app);
     ofs << PhysiCell_globals.current_time<<", "<< pCell->custom_data["Temperature"]<<",  "<< pCell->custom_data["dVcs_dt_1"]<<", "<< pCell->custom_data["Active_cell_volume"]<<", "<<pCell->phenotype.volume.total<<std::endl;
     ofs.close();
     ofs.open ("output/rat_hepatocyte_Tau_22cells_B400.txt" , std::ofstream::out | std::ofstream::app);
     ofs << pCell->custom_data["Temperature"]<<",  "<<PhysiCell_globals.current_time<<",  "<<pCell->custom_data["tau_t"]<< std::endl;
     ofs.close();
     
    // Check the Monte-Carlo saved results for each cell (agent). The Monte-Carlo 
    //results indicate which cell is frozen.	
     if (pCell->custom_data["Freezing_time"]<= pCell->custom_data["tau_t"] && abs(pCell->custom_data["Freezing_time"]- pCell->custom_data["tau_t"])<0.1&& pCell->custom_data["Active_cell_volume"]>0)  
     {
		 // Update the state of the cell to be frozen
         pCell->state.ice_state =1;
		 // Update the volume of the cell, the volume expand
         pCell->set_total_volume(pCell->custom_data["Prev_total_volume"]/0.9);
         // Write some output data 
		 ofs.open ("output/rathepatocyte_freezing_time_B400.txt" , std::ofstream::out | std::ofstream::app);
         ofs << PhysiCell_globals.current_time<<", "<< pCell<<",  "<< pCell->custom_data["Freezing_time"]<<","<<pCell->custom_data["Temperature"]<< ",  "<<pCell->custom_data["tau_t"]<< ","<<pCell->state.ice_state<< std::endl;
         ofs.close(); 
     }
 }
 //Update the position of the cell as it loses water or freezes
 pCell->update_position(dt); 
	return;
}


void hepatocyte_phenotype_rule( Cell* pCell, Phenotype& phenotype, double dt )
{

	return;
}

// Compute the water change (transport) model for each cell at the given time
void WaterTransportModel_NoCPA(Cell* pCell, Phenotype& phenotype, double dt )
{
	// gets the volume of each voxel, they are all the same
	static double voxel_volume = microenvironment.mesh.voxels[0].volume;
	double TR = 273.15;//refrence Temperature
	double R1 =8.3145e-3; //KJ/molKelvin
	double R2 =0.0821e15; //micm^3 atm /molKelvin
    double B  = 6.6667; //cooling rat
    double Vw = 18e12;//mm^3/mol //Partial_molar_volume of water
    double delta_H = 5.94126e16;//mm3 atm/mol
    double Vs      = 2;
	// get parameters for this cell type
	double Lpg = pCell->custom_data["Lpg"];
	double Elp = pCell->custom_data["Elp"];//activation energy //kJ/mol 
	double Area= pCell->custom_data["Area"]; //gets the surface area of this cell
    double Vb_fraction      = pCell->custom_data["Vb_fraction"]; 
    double Isosmotic_volume = pCell->custom_data["isosmotic_volume"];// gets the  TOTAL isosmotic_volume from the current cell
    double Vb_volume        = Vb_fraction*Isosmotic_volume;
	
    double Vcv  = pCell->custom_data["Prev_active_cell_volume"];
    double Csalt= pCell->custom_data["Csalt"];//mol/kg=mol/L=fmol/mm^3
    double ns   = Csalt*pCell->custom_data["Active_cell_volume"];

   
	double t    = PhysiCell_globals.current_time;

    double Lp   = Lpg*exp(-(Elp/R1)*((1/pCell->custom_data["Temperature"])-(1/TR)));
 
 
    double dVcs_dt_2=-(1/Vw)*(Lp*Area*R2*pCell->custom_data["Temperature"])*(-(delta_H/R2)*((1/TR)-(1/pCell->custom_data["Temperature"]))+log(Vcv/(Vcv+Vw*(Vs*ns))));
    //calculate first step by forward euler
	if(PhysiCell_globals.current_time<0.1)
	{
		//std::cout<< "step 1"<< std::endl;
		//first step calculations
		pCell->custom_data["Active_cell_volume"]=pCell->custom_data["Prev_active_cell_volume"]+dVcs_dt_2*dt;
		pCell->custom_data["dVcs_dt_1"]=dVcs_dt_2;
		pCell->custom_data["Prev_active_cell_volume"]=pCell->custom_data["Active_cell_volume"];
		
	}
    //2nd order Adams-Bashforth
	else 
	{
		pCell->custom_data["Active_cell_volume"]=pCell->custom_data["Prev_active_cell_volume"]+(dt/2)*(3*(dVcs_dt_2)-(pCell->custom_data["dVcs_dt_1"]));

		pCell->custom_data["dVcs_dt_1"]=dVcs_dt_2;

		pCell->custom_data["Prev_active_cell_volume"]=pCell->custom_data["Active_cell_volume"];

	}

    pCell->set_total_volume(pCell->custom_data["Active_cell_volume"]+Vb_volume);

	pCell->custom_data["Prev_total_volume"]=pCell->phenotype.volume.total;

    pCell->custom_data["Temperature"] =272.15000-B*t;

	return;
}

// Compute the dimensionless time for each cell at the given time based on
// the water content and the physics
void dimensionless_time_NoCPA(Cell* pCell, Phenotype& phenotype, double dt )
{

    double Csalt= pCell->custom_data["Csalt"];//mol/kg=mol/L=fmol/mm^3
    double V_cell=pCell->phenotype.volume.total;
    double Vb_fraction=pCell->custom_data["Vb_fraction"]; 
    double Isosmotic_volume = pCell->custom_data["isosmotic_volume"];
    double Vb_volume= Vb_fraction*Isosmotic_volume;
   
    double Vcv_initail=Isosmotic_volume*(1-Vb_fraction);
    double Vcv=pCell->custom_data["Active_cell_volume"];
    double ns = Csalt*Vcv;
    double nw= 55.3419*Vcv;
    double Xw= nw/(nw+2*ns);


    double omega= pCell->custom_data["omega"];
    double kapa = pCell->custom_data["kapa"];
    double Q = 0.609375;
    double phi_m ;
    double eta_w0= 11.1289;
    double eta_0 = 0.0798;
    
    double B = 6.6667;//Cooling rate
 
    
    double R2      =0.0821e15; //micm^3 atm /molKelvin
    double delta_H = 5.94126e16;//mm3 atm/mol//59.4126; //mm3 atm/fmol 
    
    double T_f0  = pCell->custom_data["T_f0"];
    double Tseed = pCell->custom_data["Tseed"];
    double T_f   = 1/((1/T_f0)-((R2/delta_H)*log(Xw)));

    double t= PhysiCell_globals.current_time;
    int    n= 10*t;
    double h= t/n;
    double Int_het=0.0;

    double tau_T =(T_f/T_f0)*(T_f/T_f0)*(T_f/T_f0)*(T_f/T_f0);
    double r3=0.75*3.14*V_cell;

    double r =  std::cbrt(r3);
	
    double Area1=4*3.14*r*r;
    double Area_0=pCell->custom_data["Area"];
    double Area=Area_0;

    
    for (int k=0;k<n;k++)
    {
        double eta = std::pow((((Tseed-B*k*h)/225)-1), 41.0/25);
        
        double J_het_kh    = omega*(eta/eta_0)*sqrt((Tseed-B*k*h)/T_f0)*exp(-kapa/((Tseed-B*k*h)*(Tseed-B*k*h)*(Tseed-B*k*h)*(T_f-(Tseed-B*k*h))*(T_f-(Tseed-B*k*h)))*(tau_T));

        if (k==0 || k==n-1)
        {
          Int_het +=0.5*h*J_het_kh;  

        }
        else
        {
        Int_het+=h*J_het_kh;

        }
    }
   
    pCell->custom_data["tau_t"] = Int_het*Area;
    
    
    return;
}

//essential functions///////////////////////////////////////////////////////////////////////////////
//mechanics functions:

//calculates absolute distance between the centers of two cells uses norm from biofvm
double distance_between(Cell* pCell_1,Cell* pCell_2)
{
	static double distance_between_cells=0.0;
	std::vector<double> displacement_between_cells=pCell_2->position-pCell_1->position;


	distance_between_cells = norm(displacement_between_cells);
	return distance_between_cells;
}
//calculate the neighborh list of each cell (agent)
//function to add a vector list of neighbors to a cells state.neighbors
 std::vector< int> populate_my_neighborhood(Cell* pCell)
{
    std::vector<int> neighbors_index;
    pCell->state.neighbors.resize(0);
	// go through all cells
	for (int i = 0; i < (*all_cells).size(); i++)
	{
		Cell* possible_neighbor=(*all_cells)[i];
		double distance=distance_between(pCell,possible_neighbor);
		// Only consider hepatocyte cells, type 4
		if(possible_neighbor->type==4 && pCell!=possible_neighbor && distance<22)
		{
             
			pCell->state.neighbors.push_back(possible_neighbor);
            neighbors_index.push_back(i);

		}
	}

	return neighbors_index;
}

// Monte-Carlo and Gillespie are stochastics models to determine the 
// probability of each state (number of frozen and unfrozen cells in a tissue) 
//as a function of dimensionless time. 
void Gillespie_Model(void)
{
    std::ofstream ofs;
	//Number of cells in the tissue
    int number_of_cells = (*all_cells).size();
    double r [2]  ; 
    r [0]=0 ,  r [1]=0;
	//Non-dimensional ice propagtion rate alpha=Jp/Ji
    double alpha =10.4;
	//Non-dimensional time 
    double time = 0.0;
    double max_time = 10.0;
	//Non-dimensional time steps
    double dtau = 0.0;
    double s[number_of_cells];

    int  f_sum =0; 
    double a0 = 1;
    int z=0;
    int k=0;
    //Will be used to obtain a seed for the random number engine
    std::random_device rd;  
	//Standard mersenne_twister_engine seeded with rd()
    std::mt19937 gen(rd()); 
	//The uniform distribution on the interval(0,1)
    std::uniform_real_distribution<> dis(0.0, 1.0);
	
    while( time < max_time)
    {
        if (a0>0)
        {
		// Output file: when all the cells are frozen	
        if (f_sum==number_of_cells)
        {
            ofs.open ("output/Gillespie_alpha10_rathepatocyte_22ncell_B400.txt" , std::ofstream::out | std::ofstream::app);
            ofs << time<<std::endl;
            ofs.close();
         }
		    // loop over the number of cells and find the neighbors for each cell
            for (int i=0;i<number_of_cells;i++)
            {
                Cell* pCell = (*all_cells)[i];                         
                
                //
                std::vector<int> neighbor_index = populate_my_neighborhood(pCell);
                z=0;

                for(int j=0; j< neighbor_index.size();j++)
                {
                    Cell* pCell_3 = (*all_cells)[neighbor_index[j]];
                    z+= pCell_3->custom_data["fvecold"];
                    
                }

                pCell->custom_data["knew"] = z;

            }

            a0 = 0;

            for (int i=0;i<number_of_cells;i++)
            {
                Cell* pCell_1=(*all_cells)[i]; 
                a0+=(1-(pCell_1->custom_data["fvec"]))*(1+(pCell_1->custom_data["knew"]*alpha));
                
            }

            // generate two random numbers
            for (int i=0;i<2;i++)
            {
                r[i]=dis(gen);

            }
                                
           
           for (int i=0;i<number_of_cells;i++)
            {
                Cell* pCell =(*all_cells)[i];
 
                if (i ==0)
                {
                    s[0]=(1-pCell->custom_data["fvec"])*(1+pCell->custom_data["knew"]*alpha);
                }
                else
                {
                s[i]=s[i-1]+(1-pCell->custom_data["fvec"])*(1+pCell->custom_data["knew"]*alpha);
                }
                
            }

            
            for (int j=0;j<number_of_cells;j++)
            { 
                Cell* pCell_2 =(*all_cells)[j];
				// Find the next reaction location based on random variable
                if (r[1]*a0<=s[j])
                {
                    pCell_2->custom_data["fvec"]=1;
                  
                    if (a0>0)
                    {
                    pCell_2->custom_data["Freezing_time"]=time-(log(r[0])/a0); 

                    }
                    break;
                    
                    
                }
            }
        // if a0=0, then all the cells are frozen
        if (a0==0)
        {
            return;
        }
        else
        {
			// Estimates the time of next reaction based on random variable
            // The next time step
			dtau =-(log(r[0])/a0);
			// The next time of reaction
            time+=dtau; 
            f_sum =0; 

            for (int i=0;i<number_of_cells;i++)
            {
                Cell* pCell =(*all_cells)[i];
                pCell->custom_data["fvecold"] = pCell->custom_data["fvec"];
                f_sum+=pCell->custom_data["fvec"];
    
            }
        }
    }//end if (a0>0)
    
}//end while
    
    return;
}
// Check the number of frozen in the tissue 
void check_tissue_freezing_state()
{
    int n_of_frozen_cells=0;
    int number_of_cells = (*all_cells).size();
    
    for (int i=0;i<number_of_cells;i++)
    {
         Cell* pCell=(*all_cells)[i];
         if (pCell->state.ice_state==1)
         {
             n_of_frozen_cells++;
         }
    }
    if (n_of_frozen_cells == number_of_cells)
    {
        PhysiCell_settings.enable_SVG_saves = false;
    }
    else
        PhysiCell_settings.enable_SVG_saves = true;   
    return;
}
