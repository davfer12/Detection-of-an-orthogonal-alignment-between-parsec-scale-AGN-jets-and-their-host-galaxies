import numpy as np
import data_loading_and_xmatching
import main_functions
import figures

# Set a random seed for reproducibility
np.random.seed(84537)

# Load Astrogeo catalogue and get a filtered version with good entries
Astrogeo_catalogue = data_loading_and_xmatching.load_Astrogeo_catalogue()
Astrogeo_good_catalogue = data_loading_and_xmatching.get_Astrogeo_good_catalogue(Astrogeo_catalogue)

# Load cross-match data for different surveys
Astrogeo_DESI, DESI_xmatches = data_loading_and_xmatching.load_xmatches('DESI')
Astrogeo_SDSS, SDSS_xmatches = data_loading_and_xmatching.load_xmatches('SDSS')
Astrogeo_DES, DES_xmatches = data_loading_and_xmatching.load_xmatches('DES')
Astrogeo_Skymapper, Skymapper_xmatches = data_loading_and_xmatching.load_xmatches('Skymapper')
Astrogeo_KIDS, KIDS_xmatches = data_loading_and_xmatching.load_xmatches('KIDS')

# Load sampled catalogues for SDSS and KIDS surveys
SDSS_sampled, KIDS_sampled = data_loading_and_xmatching.load_sampled_catalogues()

# Combine data from different surveys
combined_data = (Astrogeo_catalogue, Astrogeo_DESI, DESI_xmatches, Astrogeo_SDSS, SDSS_xmatches, 
                  Astrogeo_DES, DES_xmatches, Astrogeo_Skymapper, Skymapper_xmatches, 
                  Astrogeo_KIDS, KIDS_xmatches)
Surveys_Data = main_functions.make_combined_data(combined_data, False)
Surveys_Good_Data = main_functions.make_combined_data(combined_data, True)

# Load Eagle simulation data
Eagle_sim = data_loading_and_xmatching.load_Eagle_sim()

# Create and save figures using the loaded and processed data

# Figure 2: Overview of different catalogues and surveys
Figure_2_data = (Astrogeo_catalogue, DESI_xmatches, SDSS_sampled, DES_xmatches, KIDS_sampled)
figures.Figure_2(Figure_2_data)

# Figure 3: Analysis with DESI data and tolerance parameters
N_bins = 5
max_tol_Z_array = np.array([10, 0.5, 0.1])
max_tol_err = 22.5
Figure_3_data = (Astrogeo_DESI, DESI_xmatches, max_tol_err, max_tol_Z_array, N_bins)
figures.Figure_3(Figure_3_data)

# Figure 4: Analysis with Eagle simulation data
Figure_4_data = Eagle_sim
figures.Figure_4(Figure_4_data)

# Supplementary Figure 1: Astrogeo catalogue pa distribution
Supplementary_Figure_1_data = Astrogeo_catalogue
figures.Supplementary_Figure_1(Supplementary_Figure_1_data)

# Supplementary Figure 2: Comparison of DESI and Skymapper pa distribution
Supplementary_Figure_2_data = (DESI_xmatches, Skymapper_xmatches)
figures.Supplementary_Figure_2(Supplementary_Figure_2_data)

# Supplementary Figure 3: Combined data analysis for various surveys
N_bins = 5
max_tol_Z_array = np.array([10, 0.5, 0.1])
max_tol_err = 22.5
Supplementary_Figure_3_data = (Surveys_Data, Surveys_Good_Data, max_tol_err, max_tol_Z_array, N_bins)
figures.Supplementary_Figure_3(Supplementary_Figure_3_data)

# Supplementary Figure 4: Analysis with SDSS data
N_bins = 5
max_tol_Z_array = np.array([10, 0.5, 0.1])
max_tol_err = 22.5
Supplementary_Figure_4_data = (Astrogeo_SDSS, SDSS_xmatches, max_tol_err, max_tol_Z_array, N_bins)
figures.Supplementary_Figure_4(Supplementary_Figure_4_data)

# Supplementary Figure 5: Analysis with DES data
N_bins = 5
max_tol_Z_array = np.array([10, 0.5, 0.1])
max_tol_err = 22.5
Supplementary_Figure_5_data = (Astrogeo_DES, DES_xmatches, max_tol_err, max_tol_Z_array, N_bins)
figures.Supplementary_Figure_5(Supplementary_Figure_5_data)

# Supplementary Figure 6: Analysis with SkyMapper data
N_bins = 5
max_tol_Z_array = np.array([10, 0.5, 0.1])
max_tol_err = 22.5
Supplementary_Figure_6_data = (Astrogeo_Skymapper, Skymapper_xmatches, max_tol_err, max_tol_Z_array, N_bins)
figures.Supplementary_Figure_6(Supplementary_Figure_6_data)

# Supplementary Figure 7: Analysis with DESI data and tolerance parameters for different bins
N_bins = 10
max_tol_Z_array = np.array([10, 0.5, 0.1])
max_tol_err = 22.5
Supplementary_Figure_7_data = (Astrogeo_DESI, DESI_xmatches, max_tol_err, max_tol_Z_array, N_bins)
figures.Supplementary_Figure_7(Supplementary_Figure_7_data)


# Supplementary Figure 8: Analysis with DESI data including magnitude tolerance
N_bins = 5
max_tol_Z_array = np.array([10, 0.5, 0.1])
max_tol_err = 22.5
max_mag_tol = 14
Supplementary_Figure_8_data = (Astrogeo_DESI, DESI_xmatches, max_tol_err, max_tol_Z_array, max_mag_tol, N_bins)
figures.Supplementary_Figure_8(Supplementary_Figure_8_data)

# Supplementary Figure 9: Analysis of mean and median Delta pa with DESI data for redshift bins
z_tol = 1
Supplementary_Figure_9_data = (Astrogeo_DESI, DESI_xmatches, z_tol)
figures.Supplementary_Figure_9(Supplementary_Figure_9_data)

# Supplementary Figure 10: Scatterplot and histogram of DESI pa vs. VLBI jet pa
z_bound = 0.1
Supplementary_Figure_10_data = (Astrogeo_DESI, DESI_xmatches, z_bound)
figures.Supplementary_Figure_10(Supplementary_Figure_10_data)

# Supplementary Figure 11: Scatterplot and histogram of DESI Delta pa vs. VLBI jet pa
z_bound = 0.1
Supplementary_Figure_11_data = (Astrogeo_DESI, DESI_xmatches, z_bound)
figures.Supplementary_Figure_11(Supplementary_Figure_11_data)





