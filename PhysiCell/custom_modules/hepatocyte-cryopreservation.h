/*
Program name       : hepatocyte-cryopreservation.h                                  
Date Written       : Septemeper 5, 2021                                            
Date Last Modified : June 8, 2023                                             
Authors            : Fatemeh Amiri                                              
Discription        : Header file to predefine and declare functions in 
                     hepatocyte_cryopreservation.cpp and some useful functions.      
*/

#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h"

using namespace BioFVM;
using namespace PhysiCell;

void create_hepatocyte_cell_type(void);
// set the  cell properties, then call the function
// to set up the cells
void create_cell_types( void );
void probability_MonteCarlo( void );

void setup_tissue();

void Gillespie_Model();
void check_tissue_freezing_state();

std::vector< int>  populate_my_neighborhood(Cell* );

// set up the microenvironment to include the immunostimulatory 
//factor
void setup_microenvironment( void );

std::vector<std::string> follicle_coloring_function( Cell* );

void sample_whole_environment(Cell* pCell, Phenotype& phenotype, 
double dt);

void hepatocyte_cell_rule( Cell* pCell, Phenotype& phenotype, double
 dt );

void hepatocyte_phenotype_rule( Cell* pCell, Phenotype& phenotype, 
double dt );

void dimensionless_time_NoCPA(Cell* pCell, Phenotype& phenotype, 
double dt );

std::vector<std::vector<double>> create_cells_hepatocyte_positions(
double cell_radius, double sphere_radius, double inner_radius);

void WaterTransportModel_NoCPA(Cell* pCell, Phenotype& phenotype, 
double dt );

//


		/** The function Fill() assigns val to all of the elements of wa. */
		//@{
		template <class T>
				inline static void Fill ( int N, T val, T *wa )
		{
			register int i;
			for ( i=0; i < N; i++ )
				wa[i] = val;
			return;
		}

		template <class T>
				inline static void Fill3D ( T val, T *wa )
		{
			wa[0] = val;
			wa[1] = val;
			wa[2] = val;
			return;
		}
		//@}
		
		/** The Copy() function copies the elements from wa1 to wa2. */
		//@{
		template <class T>
				inline static void Copy3D ( const T *wa1, T *wa2 )
		{
			wa2[0] = wa1[0];
			wa2[1] = wa1[1];
			wa2[2] = wa1[2];
			return;
		}
		template <class T>
				inline static void Copy ( int N, const T *wa1, T *wa2 )
		{
			register int i;
			for ( i=0; i < N; i++ )
				wa2[i] = wa1[i];
			return;
		}
        template <class T>
                inline static void Copy ( int N, int M, const T *wa1, T *wa2 )
        {
            register int i, j, row;
            for ( i=0; i < N; i++ )
            {
                row = i*M;
                for ( j=0; j<M; j++ )
                    wa2[row+j] = wa1[row+j];
            }
            return;
        }
		template <class T>
				inline static void Copy ( int N, std::vector<T> wa1, T *wa2 )
		{
			register int i;
			for ( i=0; i < N; i++ )
				wa2[i] = wa1[i];
			return;
		}
		//@}
