This code implementes Gillespie (Monte-Carlo) model of ice formation and propagation in hepatocyte disks monolayer of 22 cells tissue, using  PhysiCell http://physicell.org/.
PhysiCell was structured in several critical subdirectories. Here we use the following subdirectories:

1. BioFVM: This includes a working copy of the BioFVM multi-substrate diffusion code. 

2. config: This directory is used for XML configuration files. Any custom configuration files
should be placed in this directory.

3. core: The core library files for PhysiCell are kept in this directory. Users should not modify functions
in the core library.

5. modules: This is where standard code modules palced for distribution with PhysiCell. 
Currently, it includes pathology, MultiCellDS, and SVG (scalable vector graphic) functions. 

4. custom_modules: This directory includes our main code for this project.

6. sample_projects: This directory includes the main.cpp and Makefile.

7. output: Output files will be here.

PhysiCell should successfully compile and run on any C++11 or later compiler that supports OpenMP.
I use g++ on Linux.

#-------------------------------------------------------------------------------------------------
To compile and run the code:
#-------------------------------------------------------------------------------------------------
1) make data-cleanup
2) make reset 
3) make hepatocyte-sample
4) make
5) ./hepatocyte-cryopreservation
#-------------------------------------------------------------------------------------------------

My 3 main code files are: "hepatocyte-cryopreservation.cpp", "hepatocyte-cryopreservation.h" and "main.cpp".  
The "hepatocyte-cryopreservation.cpp" and "hepatocyte-cryopreservation.h" are located in /Codes/PhysiCell/custom_modules
The "main.cpp" is located in /Codes/PhysiCell/sample_projects/hepatocyte
The results are saved in /Codes/PhysiCell/output as .SVG files. In each .SVG file the red cell indicates the unfrozen cell and the blue cell with a red circle inside indicates the frozen cell in the tissue.
