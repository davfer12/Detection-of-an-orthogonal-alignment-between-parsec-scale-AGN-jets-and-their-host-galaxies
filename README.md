# VLBI_optical_project
This is the necessary code and data to reproduce the figures (except Figure 1, which was done by hand) and results of the Nature Astro paper 'Detection of an orthogonal alignment between parsec scale AGN
jets and their host galaxies' (D. Fernández Gil, J. A. Hodgson, B. L’Huillier, J. Asorey, C. Saulder, K. Finner, M. J. Jee, D. Parkinson, F. Combes). 

## File descriptions

The 'requirements.txt' lists all necessary python packages and versions that must be installed. Python 3.0 or newer is required.

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
- Astrogeo_catalogue_no_redshift.csv: The name and coordinates of the VLBI radio sources without a measured redshift.
- Astrogeo_{survey_name}.csv: The cross-matched Astrogeo VLBI data with the corresponding optical survey.
- {survey_name}_xmatches.csv: The corresponding cross-matched optical survey data with the Astrogeo VLBI catalogue.
- SDSS_sampled_catalogue.csv and KIDS_sampled_catalogue.csv: The randomly sampled SDSS and KIDS data used in the coverage map of Figure 2.
- Eagle_sim.csv: The Eagle simulation data.

The 'paper_figures' folder includes all the Figures of the paper. 

## Instructions of use
If the user wants to reproduce the figures, they have to follow these simple steps:
- Save all these files together with the 'data' and 'paper_figures' folders. One can delete the figures inside 'paper_figures' if they wish.
- In the same directory where all the files and folders are saved, first run 'pip install -r requirements.txt' to install all python dependencies.
- Run the 'main.py' script with python 3.0 or newer. One can do so with the following command: 'python3 main.py > out.out' to get an output file with the progress.
The figures will be saved in the 'paper_figures' folder, overwriting them if they were not previously removed from that folder.  
