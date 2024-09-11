# VLBI_optical_project
This is the necesary code and data to reproduce the figures and results of the Nature Astro paper 'Detection of an orthogonal alignment between parsec scale AGN
jets and their host galaxies' (D. Fernández Gil, J. A. Hodgson, B. L’Huillier, J. Asorey, C. Saulder, K. Finner, M. J. Jee, D. Parkinson, F. Combes). 

##File descriptions

The 'requirements.txt' lists all neccesary python packages and versions that must be installed.

The 'libraries.py' file imports all necessary python libraries and dependencies.

The 'data_loading_and_xmatching.py' file contains all functions pertaining to loading the catalogues, performing the cross-matches, saving them into csv files and loading them afterwards.

The 'get_catalogue_parameters.py' file contains functions that access the catalogue data and gets their relevantr parameters.

The 'aux_functions.py' file contains auxiliary functions that perform various tasks such as dealing with angles, computing the optical/VLBI angle difference, carrying out a MonteCarlo technique, obtaining statistical siognificances and handling 3D rotations for the Eagle simulation results.

The 'main_functions.py' file contains the root of the data handlin and ganalysis. It creates a combined dataset from the individual catalogues, it prepares the data for the histogram figures, obtains the p-values and plots it.

The 'figures.py' file contains functions that generate every figure of the paper.

The 'pregenerate_neccesary_data.py' script deals with taking the original optical catalogue data and performing a cross-match with the VLBI data, as well as taking a random sample from two of the catalogues which are used in the coverage map of Figure 2. We include this code here for transparency purposes, but the original optical catalogue data is too heavy to upload. We only upload the cross-matched outputs generated from this script in the data files. 

The 'main.py' script uses all other files to load the catalogues and cross-matches and execute every function to generate all the Figures.

The 'data' folder includes the following files:
- Astrogeo_catalogue.csv: The whole VLBI Astrogeo data.
- Astrogeo_{catalogue}.csv: The cross-matched Astrogeo VLBI data with the corresponding optical catalogue.
- {catalogue}_xmatches.csv: The corresponding cross-matched optical catalogue data with the Astrogeo VLBI catalogue.
whole VLBI Astrogeo data, the mentioned optical/VLBI cross-matches, the sampled optical catalogues and the Eagle simulation data.
- SDSS_sampled_catalogue.csv,KIDS_sampled_catalogue.csv: The randomly sampled SDSS and KIDS data used in the coverage map of Figure 2.
- Eagle_sim.csv: The Eagle simulation data.

The 'paper_figures' folder includes all the Figures of the paper. 

## Instructions of use
If the user wants to reproduce the figures, they have to follow this simple steps:
- Save all these files together with the 'data' folder.
- In the same directory where all the files and the data folder is saved, first run 'pip install -r requirements.txt' to install all python dependencies.
- Run the 'main.py' script. One can do so with the following command: 'python3.8 main.py > out.out' to get an output file with the progress.
The figures will be saved in the 'paper_figures' folder. 
