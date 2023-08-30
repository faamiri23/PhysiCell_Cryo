/*
###############################################################################                                                                
#Program name       : main file                                               #                           
#Date Written       : Septemeper 5, 2021                                      #
#Date Last Modified : June 8, 2023                                            #
#Programmers name   : Fatemeh Amiri                                           #
#Purpose            : This study applied a lattice-free agent-based model in  #
#combination with stochastic models for ice formation and propagtion          #
#in the dic monolayer rat hepatocyte tissue.                                  #
#Discription:                                                                 #
# We implemented and solved the Monte-Carlo ice formation and propagtion model# 
# in small tissue using PhysiCell (Version 1.7.1) [1],                        #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
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

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <cmath>
#include <omp.h>
#include <fstream>

#include "./core/PhysiCell.h"
#include "./modules/PhysiCell_standard_modules.h"

// custom user modules
#include "./custom_modules/hepatocyte-cryopreservation.h"

using namespace BioFVM;
using namespace PhysiCell;

int main( int argc, char* argv[] )
{
	// find and read(load) the configuration input file
        bool XML_status = false;
	if( argc > 1 )
	{ XML_status = load_PhysiCell_config_file( argv[1] ); }
	else
	{ XML_status = load_PhysiCell_config_file( "./config/PhysiCell_settings.xml" ); }
	if( !XML_status )
	{ exit(-1); }

	// OpenMP setup
	omp_set_num_threads(PhysiCell_settings.omp_num_threads);


	// time setup
	std::string time_units = "sec";

	/* Microenvironment setup */

	setup_microenvironment();//set 

	/* PhysiCell setup */

	// set mechanics voxel size, and match the data structure to BioFVM
	double mechanics_voxel_size = 8.1;
	Cell_Container* cell_container = create_cell_container_for_microenvironment( microenvironment, mechanics_voxel_size );

    // Define and create the cells type
	create_cell_types();

    // Define the Geometry of the tissue and the positions of the cells in the tissue
	setup_tissue();
    // Simulate the probability of ice formation and propagtion
	 Gillespie_Model();
	 
    // Trun off the matlab format saving ouput
	set_save_biofvm_mesh_as_matlab( false);
	set_save_biofvm_data_as_matlab( false );
	set_save_biofvm_cell_data( false );
	set_save_biofvm_cell_data_as_custom_matlab( false );

	// save a simulation snapshot at initial time t=0
	char filename[1024];
	sprintf( filename , "%s/initial" , PhysiCell_settings.folder.c_str() );
	save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time );
	
	// save a quick SVG cross section through z = 0
        PhysiCell_SVG_options.length_bar = 100;

	// for simplicity, set a pathology coloring function
	std::vector<std::string> (*cell_coloring_function)(Cell*) = follicle_coloring_function;

	sprintf( filename , "%s/initial.svg" , PhysiCell_settings.folder.c_str() );
	SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );


	// set the performance timers
    BioFVM::RUNTIME_TIC();
	BioFVM::TIC();


	// The main loop starts here, loop over time
       
	 while( PhysiCell_globals.current_time < PhysiCell_settings.max_time + 0.1*diffusion_dt)
		{
          #pragma omp critical  
                // SVG output is here
				if( PhysiCell_settings.enable_SVG_saves == true )
				{
					
					sprintf( filename , "%s/snapshot%08u.svg" , PhysiCell_settings.folder.c_str() , PhysiCell_globals.SVG_output_index );
					SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );

					PhysiCell_globals.SVG_output_index++;
					PhysiCell_globals.next_SVG_save_time  += PhysiCell_settings.SVG_save_interval;
					
				}
				// Check if the whole tissue is frozen to stop the saving output files
				check_tissue_freezing_state();
			// update the microenvironment
			microenvironment.simulate_diffusion_decay( diffusion_dt);//diffusion_dt );


            //Run PhysiCell (goes through all cells)
			((Cell_Container *)microenvironment.agent_container)->update_all_cells( PhysiCell_globals.current_time );

			PhysiCell_globals.current_time += diffusion_dt;
		}
            

	// save a final simulation snapshot
	sprintf( filename , "%s/final" , PhysiCell_settings.folder.c_str() );
	save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time );

	sprintf( filename , "%s/final.svg" , PhysiCell_settings.folder.c_str() );
	SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );

	return 0;
}
