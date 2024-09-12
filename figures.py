from libraries import *
import data_loading_and_xmatching
import get_catalogue_parameters
import aux_functions
import main_functions

import numpy as np
import matplotlib.pyplot as plt

def Figure_2(data):
    """
    Create a Mollweide projection plot displaying the spatial distribution of different
    astronomical surveys and VLBI data.
    
    Parameters:
    -----------
    data : tuple
        A tuple containing:
            - Astrogeo_catalogue: VLBI data catalogue.
            - DESI_xmatches: DESI cross-match data.
            - SDSS_sampled: SDSS sampled data.
            - DES_xmatches: DES cross-match data.
            - KIDS_sampled: KiDS sampled data.
    """
    Astrogeo_catalogue, DESI_xmatches, SDSS_sampled, DES_xmatches, KIDS_sampled = data
    
    
    # VLBI data
    Astrogeo_RA = np.array(Astrogeo_catalogue.RA_deg)
    Astrogeo_DEC = np.array(Astrogeo_catalogue.DEC_deg)
    Astrogeo_RA, Astrogeo_DEC = aux_functions.process_coordinates(Astrogeo_RA, Astrogeo_DEC)

    # DESI data
    DESI_RA = np.array(DESI_xmatches.RA)
    DESI_DEC = np.array(DESI_xmatches.DEC)
    DESI_RA, DESI_DEC = aux_functions.process_coordinates(DESI_RA, DESI_DEC)
    
    # SDSS data
    SDSS_RA = np.array(SDSS_sampled.ra)
    SDSS_DEC = np.array(SDSS_sampled.dec)
    SDSS_RA, SDSS_DEC = aux_functions.process_coordinates(SDSS_RA, SDSS_DEC)
    
    # DES data
    DES_RA = np.array(DES_xmatches.DES_full_catalogue_Astrogeo_matches_RA)
    DES_DEC = np.array(DES_xmatches.DES_full_catalogue_Astrogeo_matches_DEC)
    DES_RA, DES_DEC = aux_functions.process_coordinates(DES_RA, DES_DEC)
    
    # SkyMapper data (randomly generated)
    size = 20000
    Skymapper_RA = np.random.uniform(-180, 180, size)
    Skymapper_DEC = np.random.uniform(-90, 0, size)
    Skymapper_RA, Skymapper_DEC = aux_functions.process_coordinates(Skymapper_RA, Skymapper_DEC)
    
    # KiDS data
    KIDS_RA = np.array(KIDS_sampled.RAJ2000)
    KIDS_DEC = np.array(KIDS_sampled.DECJ2000)
    KIDS_RA, KIDS_DEC = aux_functions.process_coordinates(KIDS_RA, KIDS_DEC)
    
    # Plotting
    plt.subplot(111, projection='mollweide')
    
    plt.scatter(Skymapper_RA, Skymapper_DEC, c='navajowhite', s=17, marker='s', alpha=0.5, label='SkyMapper')
    plt.scatter(DESI_RA, DESI_DEC, c='salmon', s=17, marker='s', alpha=0.5, label='DESI LS')
    plt.scatter(SDSS_RA, SDSS_DEC, c='red', s=17, marker='s', alpha=0.5, label='SDSS')
    plt.scatter(DES_RA, DES_DEC, c='darkorange', s=17, marker='s', alpha=0.5, label='DES')
    plt.scatter(KIDS_RA, KIDS_DEC, c='maroon', s=17, marker='s', alpha=0.5, label='KiDS')
    plt.scatter(Astrogeo_RA, Astrogeo_DEC, c='cyan', s=1, marker='.', alpha=1, label='VLBI')
    
    lgnd = plt.legend(bbox_to_anchor=(0.85, 1.3), ncol=3, frameon=False)
    for lh in lgnd.legendHandles:
        lh._sizes = [30]  
        lh.set_alpha(1)
    
    plt.xlabel(r'RA ($^{{\circ}}$)')
    plt.ylabel(r'DEC ($^{{\circ}}$)')
    plt.grid(True)
    
    plt.savefig('paper_figures/Figure_2.pdf', dpi=300, bbox_inches='tight')
    plt.show()

    


def Figure_3(data):
    """
    Generate histograms and a plot for the DESI dataset, comparing VLBI and optical
    polarization angles and their differences.
    
    Parameters:
    -----------
    data : tuple
        A tuple containing:
            - Astrogeo_DESI: Astrogeo data for DESI.
            - DESI_xmatches: DESI cross-match data.
            - max_tol_err: Maximum tolerance for error.
            - max_tol_Z_array: Array of maximum tolerance values for redshift.
            - N_bins: Number of bins for histograms.
    """
    
    # Unpack the input data
    Astrogeo_DESI, DESI_xmatches, max_tol_err, max_tol_Z_array, N_bins = data
    
    # Retrieve catalog parameters for DESI
    (DESI_OPTICAL_PA, DESI_OPTICAL_PA_ERR_ORIGINAL, DESI_VLBI_PA, 
     DESI_VLBI_PA_ERR_ORIGINAL, DESI_PA_DIFF, DESI_PA_DIFF_ERR_ORIGINAL, 
     DESI_SOURCE_Z, DESI_catalogue_Astrogeo_B_matches, 
     DESI_catalogue_Astrogeo_TYPE_matches, 
     DESI_catalogue_Astrogeo_MAG_Z_matches) = get_catalogue_parameters.get_DESI_data(Astrogeo_DESI, DESI_xmatches)
    
    # Define good cases based on catalog types
    DESI_good_cases_cut = np.in1d(DESI_catalogue_Astrogeo_TYPE_matches, ['SER', 'EXP', 'DEV'])
    
    # Define survey name and output file name
    Survey_Name = 'DESI LS'
    Png_Name = 'Figure_3.pdf'
    
    # Apply cut based on the magnitude limit
    b_cut = DESI_catalogue_Astrogeo_B_matches > 1.3
    
    # Compute histograms and p-values
    final_p_values_5_bins, final_p_values_2_bins = main_functions.Histograms(
        DESI_VLBI_PA[b_cut],
        DESI_VLBI_PA_ERR_ORIGINAL[b_cut],
        DESI_OPTICAL_PA[b_cut],
        DESI_OPTICAL_PA_ERR_ORIGINAL[b_cut],
        DESI_PA_DIFF[b_cut],
        DESI_PA_DIFF_ERR_ORIGINAL[b_cut],
        DESI_SOURCE_Z[b_cut],
        max_tol_err,
        max_tol_Z_array,
        Survey_Name,
        Png_Name,
        N_bins,
        Combined=False,
        good_cases_cut=DESI_good_cases_cut[b_cut]
    )
    
    
    
def Figure_4(data):
    """
    Create and save plots for the Eagle simulation data, analyzing the position angle
    differences with various spin orientation models.

    Parameters:
    -----------
    data : dict
        A dictionary containing:
            - 'u_nodust': u-band magnitudes of galaxies
            - 'g_nodust': g-band magnitudes
            - 'r_nodust': r-band magnitudes
            - 'i_nodust': i-band magnitudes
            - 'z_nodust': z-band magnitudes
            - 'pa_r': Polarization angles in r-band
            - 'stars_Spin_x': Spin orientation x-component
            - 'stars_Spin_y': Spin orientation y-component
            - 'stars_Spin_z': Spin orientation z-component
    """
    
    # Extract data
    u = data['u_nodust']
    g = data['g_nodust']
    r = data['r_nodust']
    i = data['i_nodust']
    z = data['z_nodust']

    # Separate ellipticals and spirals
    Eagle_sim_ellipticals = data[(u - r) > 2.2]
    Eagle_sim_spirals = data[(u - r) < 2.2]
    
    r_PA = np.array(Eagle_sim_ellipticals['pa_r'])
    star_Spin_X = np.array(Eagle_sim_ellipticals['stars_Spin_x'])
    star_Spin_Y = np.array(Eagle_sim_ellipticals['stars_Spin_y'])
    star_Spin_Z = np.array(Eagle_sim_ellipticals['stars_Spin_z'])
    
    star_Spins = np.array([star_Spin_X, star_Spin_Y, star_Spin_Z]).T

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 10))

    print('Computing with original data...')
    
    # Compute jet PA and delta PA for original data
    jet_PA = np.arctan2(star_Spin_Y, star_Spin_X) * 180 / np.pi
    jet_PA[jet_PA < 0] += 180
    jet_PA -= 90

    delta_PA = np.abs(r_PA - jet_PA)
    delta_PA[delta_PA > 90] = 180 - delta_PA[delta_PA > 90]
    delta_PA = 90 - delta_PA

    # Plot histogram for original data
    bins_array = np.linspace(0, 90, 6)
    bins_center = bins_array[:-1] + 0.5*np.diff(bins_array)
    counts, _ = np.histogram(delta_PA, bins=bins_array)
    
    ax1.hist(delta_PA, bins=bins_array)
    ax1.set_xlabel(r'$\Delta$PA ($^{{\circ}}$)')
    ax1.text(0.5,-0.2, "(a)", size=12, ha="center", transform=ax1.transAxes)
    ax1.set_title(r'Eagle simulation elliptical galaxies with no scatter')
    ax1.text(0.5, -0.2, "(a)", size=12, ha="center", transform=ax1.transAxes)

    print('Computing with uniform spin orientations...')
    
    # Compute delta PA with uniform spin orientations
    delta_PAs = []
    for i in range(1000):
        theta = np.random.uniform(-np.pi, np.pi, len(star_Spins))
        phi = np.random.uniform(0, 2 * np.pi, len(star_Spins))
        star_Spin_X_p = np.cos(theta) * np.cos(phi)
        star_Spin_Y_p = np.cos(theta) * np.sin(phi)
        star_Spin_Z_p = np.sin(theta)
        
        new_jet_PA = np.arctan2(star_Spin_Y_p, star_Spin_X_p) * 180 / np.pi
        new_jet_PA[new_jet_PA < 0] += 180
        new_jet_PA -= 90
        
        delta_PA = np.abs(r_PA - new_jet_PA)
        delta_PA[delta_PA > 90] = 180 - delta_PA[delta_PA > 90]
        delta_PA = 90 - delta_PA
        delta_PAs.append(delta_PA)
        print(f'{i/10:.1f}%', end='\r')

    counts_array = np.array([np.histogram(d, bins=bins_array)[0] for d in delta_PAs])
    counts_avg = np.mean(counts_array, axis=0)
    counts_std = np.std(counts_array, axis=0)

    # Plot histogram for uniform spin orientations
    ax2.bar(bins_array[:-1], counts_avg, bins_array[1] - bins_array[0], align='edge')
    ax2.errorbar(bins_center, counts_avg, yerr=counts_std, fmt='None', capsize=5, ecolor='k')
    ax2.set_xlabel(r'$\Delta$PA ($^{{\circ}}$)')
    ax2.text(0.5,-0.2, "(b)", size=12, ha="center", transform=ax2.transAxes)
    ax2.set_title(r'Eagle simulation elliptical galaxies with uniform spin orientations')
    ax2.set_ylim([0,225])

    epsilon = 0.33
    print(f'Computing with Gaussian scatter with epsilon = {epsilon}...')
    
    # Compute delta PA with Gaussian scatter
    theta_lim = epsilon * 90
    delta_PAs = []
    for i in range(1000):
        cos_theta = np.random.normal(np.cos(theta_lim*np.pi/180),np.cos(theta_lim*np.pi/180)*epsilon,len(star_Spins))
        theta = np.arccos(cos_theta)*180/np.pi
        phi = np.random.uniform(0,360,len(star_Spins))
        star_Spins_p = aux_functions.rot(star_Spins,theta,phi)
        star_Spin_X_p = star_Spins_p[:,0]
        star_Spin_Y_p = star_Spins_p[:,1]

        new_jet_PA = np.arctan2(star_Spin_Y_p,star_Spin_X_p)*180/np.pi
        new_jet_PA[new_jet_PA < 0] = new_jet_PA[new_jet_PA < 0] + 180
        new_jet_PA =  new_jet_PA - 90


        delta_PA = np.abs(r_PA - new_jet_PA)
        delta_PA[delta_PA > 90] = 180 - delta_PA[delta_PA > 90]
        delta_PA = 90 - delta_PA
        delta_PAs.append(delta_PA)
        print('{}%'.format(i/10),end='\r')


    n_bins = 5
    bins_array = np.linspace(0,90,n_bins+1)
    bins_center = 0.5*(bins_array[1:] + bins_array[:-1])
    counts_array = []
    for i in range(1000):
        counts,_ = np.histogram(delta_PAs[i],bins=bins_array)
        counts_array.append(counts)
    counts_array = np.array(counts_array)

    counts_avg = np.mean(counts_array, axis=0)
    counts_std = np.std(counts_array, axis=0)

    # Plot histogram for Gaussian scatter
    ax3.bar(bins_array[:-1], counts_avg, bins_array[1] - bins_array[0], align='edge')
    ax3.errorbar(bins_center, counts_avg, yerr=counts_std, fmt='None', capsize=5, ecolor='k')
    ax3.set_xlabel(r'$\Delta$PA ($^{{\circ}}$)')
    ax3.text(0.5,-0.2, "(c)", size=12, ha="center", transform=ax3.transAxes)
    ax3.set_title(r'Eagle simulation elliptical galaxies with gaussian scatter with $\epsilon$ = {}'.format(epsilon))

    epsilon = 0.33
    print(f'Computing with uniform scatter with epsilon = {epsilon}...')
    
    # Compute delta PA with uniform scatter
    theta_lim = epsilon * 90
    delta_PAs = []
    for i in range(1000):
        cos_theta = np.random.uniform(np.cos(theta_lim*np.pi/180),1,len(star_Spins))
        theta = np.arccos(cos_theta)*180/np.pi
        phi = np.random.uniform(0,360,len(star_Spins))
        star_Spins_p = aux_functions.rot(star_Spins,theta,phi)
        star_Spin_X_p = star_Spins_p[:,0]
        star_Spin_Y_p = star_Spins_p[:,1]

        new_jet_PA = np.arctan2(star_Spin_Y_p,star_Spin_X_p)*180/np.pi
        new_jet_PA[new_jet_PA < 0] = new_jet_PA[new_jet_PA < 0] + 180
        new_jet_PA =  new_jet_PA - 90


        delta_PA = np.abs(r_PA - new_jet_PA)
        delta_PA[delta_PA > 90] = 180 - delta_PA[delta_PA > 90]
        delta_PA = 90 - delta_PA
        delta_PAs.append(delta_PA)
        print('{}%'.format(i/10),end='\r')


    n_bins = 5
    bins_array = np.linspace(0,90,n_bins+1)
    bins_center = 0.5*(bins_array[1:] + bins_array[:-1])
    counts_array = []
    for i in range(1000):
        counts,_ = np.histogram(delta_PAs[i],bins=bins_array)
        counts_array.append(counts)
    counts_array = np.array(counts_array)

    counts_avg = np.mean(counts_array, axis=0)
    counts_std = np.std(counts_array, axis=0)

    # Plot histogram for uniform scatter
    ax4.bar(bins_array[:-1], counts_avg, bins_array[1] - bins_array[0], align='edge')
    ax4.errorbar(bins_center, counts_avg, yerr=counts_std, fmt='None', capsize=5, ecolor='k', label='Counts +/- SD of bin height')
    ax4.set_xlabel(r'$\Delta$PA ($^{{\circ}}$)')
    ax4.text(0.5,-0.2, "(d)", size=12, ha="center", transform=ax4.transAxes)
    ax4.set_title(r'Eagle simulation elliptical galaxies with uniform scatter with $\epsilon$ = {}'.format(epsilon))

    # Adjust layout and save figure
    fig.legend(bbox_to_anchor=(0.9, 1.0), frameon=False, borderaxespad=0., fontsize=15, bbox_transform=fig.transFigure)
    plt.subplots_adjust(hspace=0.4)
    plt.savefig('paper_figures/Figure_4.pdf', dpi=100, bbox_inches='tight')
    plt.show()
    
    
    
    
def Figure_M1(data):
    """
    Generate and save a histogram plot of position angles (PA) from the provided data.

    Parameters:
    -----------
    data : pd.DataFrame
        DataFrame containing the 'pa' column with position angles.
    """
    
    # Extract polarization angles from the data
    pa = np.array(data.pa)
    
    # Create histogram plot
    plt.figure(figsize=(8, 6))
    plt.hist(pa, histtype='step', fill=True, linewidth=3, alpha=0.2, bins=11)
    
    # Set labels and title
    plt.xlabel(r'PA ($^{\circ}$)', fontsize=15)
    plt.ylabel('Counts', fontsize=15)
    plt.title('VLBI Jets', fontsize=20)
    
    # Customize tick parameters
    plt.tick_params(axis='both', which='major', labelsize=15)
    
    # Adjust layout
    plt.gcf().subplots_adjust(bottom=0.15)
    
    # Save and show plot
    plt.savefig('paper_figures/Figure_M1.pdf', dpi=200)
    plt.show()
    
    
def Figure_M2(data):
    """
    Generate and save histograms comparing position angles for DESI and SkyMapper optical cross-matches.

    Parameters:
    -----------
    data : tuple
        A tuple containing:
        - DESI_xmatches: DataFrame with a column 'pos_angle' for DESI cross-matches.
        - Skymapper_xmatches: DataFrame with a column 'PA' for SkyMapper cross-matches.
    """
    
    DESI_xmatches, Skymapper_xmatches = data
    
    # Process DESI cross-matches
    DESI_angles = np.array(DESI_xmatches.pos_angle)
    DESI_angles = np.fmod(DESI_angles + 360, 180)
    DESI_angles[DESI_angles > 90] -= 180
    
    # Process SkyMapper cross-matches
    Skymapper_angles = np.array(Skymapper_xmatches.PA)
    Skymapper_angles = np.fmod(Skymapper_angles + 360, 180)
    Skymapper_angles[Skymapper_angles > 90] -= 180
    
    # Create subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    
    # Plot DESI histogram
    ax1.hist(DESI_angles, histtype='step', fill=True, linewidth=3, alpha=0.2, bins=11)
    ax1.set_xlabel(r'PA ($^{\circ}$)', fontsize=15)
    ax1.set_ylabel('Counts', fontsize=15)
    ax1.tick_params(axis='both', which='major', labelsize=15)
    ax1.text(0.5, -0.2, "(a)", size=12, ha="center", transform=ax1.transAxes)
    ax1.set_title('DESI LS optical cross-matches')
    
    # Plot SkyMapper histogram
    ax2.hist(Skymapper_angles, histtype='step', fill=True, linewidth=3, alpha=0.2, bins=11)
    ax2.set_xlabel(r'PA ($^{\circ}$)', fontsize=15)
    ax2.set_ylabel('Counts', fontsize=15)
    ax2.tick_params(axis='both', which='major', labelsize=15)
    ax2.text(0.5, -0.2, "(b)", size=12, ha="center", transform=ax2.transAxes)
    ax2.set_title('SkyMapper optical cross-matches')
    
    # Adjust layout
    plt.gcf().subplots_adjust(bottom=0.15)
    
    # Save and show plot
    plt.savefig('paper_figures/Figure_M2.pdf', dpi=200)
    plt.show()
    
    
def Figure_M3(data):
    """
    Generate and save histogram plots comparing position angles (PA) between DESI VLBI and optical data.

    Parameters:
    ----------
    data : tuple
        Contains the following elements:
        - Astrogeo_DESI : DataFrame
            Data from the Astrogeo and DESI catalogues.
        - DESI_xmatches : DataFrame
            Cross-match data between DESI and Astrogeo.
        - max_tol_err : float
            Maximum tolerance for errors in position angle.
        - max_tol_Z_array : array
            Array of redshift tolerance values.
        - max_mag_tol : float
            Maximum magnitude tolerance for the dataset.
        - N_bins : int
            Number of bins to use in the histogram.
    """
    # Unpack data
    Astrogeo_DESI, DESI_xmatches, max_tol_err, max_tol_Z_array, max_mag_tol, N_bins = data
    
    # Extract position angles and associated errors from the DESI data
    DESI_OPTICAL_PA, DESI_OPTICAL_PA_ERR_ORIGINAL, DESI_VLBI_PA, DESI_VLBI_PA_ERR_ORIGINAL, \
    DESI_PA_DIFF, DESI_PA_DIFF_ERR_ORIGINAL, DESI_SOURCE_Z, DESI_catalogue_Astrogeo_B_matches, \
    DESI_catalogue_Astrogeo_TYPE_matches, DESI_catalogue_Astrogeo_MAG_Z_matches = \
        get_catalogue_parameters.get_DESI_data(Astrogeo_DESI, DESI_xmatches)
    
    # Filter to only include the specified types of sources
    DESI_good_cases_cut = np.in1d(DESI_catalogue_Astrogeo_TYPE_matches, ['SER', 'EXP', 'DEV'])
    
    # Define survey name and output file name
    Survey_Name = 'DESI LS' + ' MAG_Z < ' + str(max_mag_tol)
    Png_Name = 'Figure_M3.pdf'
    
    # Apply magnitude and B-magnitude cuts
    mag_cut = DESI_catalogue_Astrogeo_MAG_Z_matches < max_mag_tol
    b_cut = DESI_catalogue_Astrogeo_B_matches > 1.3
    
    # Generate histograms comparing VLBI and optical position angles
    final_p_values_5_bins, final_p_values_2_bins = main_functions.Histograms(
        DESI_VLBI_PA[b_cut & mag_cut], 
        DESI_VLBI_PA_ERR_ORIGINAL[b_cut & mag_cut], 
        DESI_OPTICAL_PA[b_cut & mag_cut], 
        DESI_OPTICAL_PA_ERR_ORIGINAL[b_cut & mag_cut], 
        DESI_PA_DIFF[b_cut & mag_cut], 
        DESI_PA_DIFF_ERR_ORIGINAL[b_cut & mag_cut], 
        DESI_SOURCE_Z[b_cut & mag_cut], 
        max_tol_err, 
        max_tol_Z_array, 
        Survey_Name, 
        Png_Name, 
        N_bins, 
        False, 
        good_cases_cut=DESI_good_cases_cut[b_cut & mag_cut]
    )

               
        
def Figure_M4(data):
    """
    Generate and save scatter plots and histograms comparing position angles (PA) between DESI VLBI and optical data.

    Parameters:
    ----------
    data : tuple
        Contains the following elements:
        - Astrogeo_DESI : DataFrame
            Data from the Astrogeo and DESI catalogues.
        - DESI_xmatches : DataFrame
            Cross-match data between DESI and Astrogeo.
        - z_bound : float
            Redshift boundary for filtering the data.
    """
    # Unpack data
    Astrogeo_DESI, DESI_xmatches, z_bound = data
    N_shuffles = 1000  # Number of shuffles for statistical analysis
    folded = True  # Whether to fold the angles

    # Extract position angles and associated errors from the DESI data
    DESI_OPTICAL_PA, DESI_OPTICAL_PA_ERR_ORIGINAL, DESI_VLBI_PA, DESI_VLBI_PA_ERR_ORIGINAL, \
    DESI_PA_DIFF, DESI_PA_DIFF_ERR_ORIGINAL, DESI_SOURCE_Z, DESI_catalogue_Astrogeo_B_matches, \
    DESI_catalogue_Astrogeo_TYPE_matches, DESI_catalogue_Astrogeo_MAG_Z_matches = \
        get_catalogue_parameters.get_DESI_data(Astrogeo_DESI, DESI_xmatches)
    
    # Filter to include only specified types of sources
    DESI_good_cases_cut = np.in1d(DESI_catalogue_Astrogeo_TYPE_matches, ['SER', 'EXP', 'DEV'])
    b_cut = DESI_catalogue_Astrogeo_B_matches > 1.3
    err_cut = DESI_OPTICAL_PA_ERR_ORIGINAL < 360
    z_cut = DESI_SOURCE_Z < z_bound

    # Create subplots
    fig, axes = plt.subplots(2, 2, figsize=(10, 10))

    # Plot 1: Scatter plot of original and modulated PA for all cases
    x_original = DESI_VLBI_PA[b_cut & err_cut]
    y_original = aux_functions.angle_conversion_major_PA(90 + DESI_OPTICAL_PA[b_cut & err_cut])
    x_err = DESI_VLBI_PA_ERR_ORIGINAL[b_cut & err_cut]
    y_err = DESI_OPTICAL_PA_ERR_ORIGINAL[b_cut & err_cut]

    x_mod, y_mod = aux_functions.mod(x_original, y_original, folded)

    ax1 = plt.subplot(221)
    ax1.set_xlim([-180, 180])
    ax1.set_ylim([-270, 270])
    ax1.scatter(x_mod, y_mod, c='b', s=0.1)
    ax1.errorbar(x_mod, y_mod, xerr=x_err, yerr=y_err, fmt='None', ecolor='k', elinewidth=0.5)
    ax1.axvline(x=0, color='k', lw=0.5)
    ax1.axhline(y=0, color='k', lw=0.5)
    ax1.plot([-270, 270], [-270, 270], 'k--', lw=2)
    ax1.set_aspect(1)
    ax1.set_facecolor('white')
    ax1.set_xlabel(r'Jet PA ($^{{\circ}}$)')
    ax1.set_ylabel(r'Closest Optical minor axis PA ($^{{\circ}}$)')
    ax1.set_title('All cases')

    # Plot 2: Histogram of differences for all cases
    ax2 = plt.subplot(222)
    ax2.hist(np.abs(y_mod - x_mod), bins=np.linspace(0, 90, 6), orientation='horizontal')
    ax2.set_ylabel(r'y-x ($^{{\circ}}$)')
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")
    ax2.set_title('All cases')

    # Plot 3: Scatter plot of original and modulated PA for good cases with z < z_bound
    x_original = DESI_VLBI_PA[DESI_good_cases_cut & z_cut & b_cut]
    y_original = aux_functions.angle_conversion_major_PA(90 + DESI_OPTICAL_PA[DESI_good_cases_cut & z_cut & b_cut])
    x_err = DESI_VLBI_PA_ERR_ORIGINAL[DESI_good_cases_cut & z_cut & b_cut]
    y_err = DESI_OPTICAL_PA_ERR_ORIGINAL[DESI_good_cases_cut & z_cut & b_cut]

    x_mod, y_mod = aux_functions.mod(x_original, y_original, folded)

    ax3 = plt.subplot(223, sharex=ax1)
    ax3.set_xlim([-180, 180])
    ax3.set_ylim([-270, 270])
    ax3.scatter(x_mod, y_mod, c='b', s=1)
    ax3.errorbar(x_mod, y_mod, xerr=x_err, yerr=y_err, fmt='None', ecolor='k', elinewidth=0.5, label='PA err.')
    ax3.axvline(x=0, color='k', lw=0.5)
    ax3.axhline(y=0, color='k', lw=0.5)
    ax3.plot([-270, 270], [-270, 270], 'k--', lw=2)
    ax3.set_aspect(1)
    ax3.set_facecolor('white')
    ax3.set_xlabel(r'Jet PA ($^{{\circ}}$)')
    ax3.set_ylabel(r'Closest Optical minor axis PA ($^{{\circ}}$)')
    ax3.set_title('Good cases, z < {}'.format(z_bound))

    # Plot 4: Histogram of differences for good cases with z < z_bound
    ax4 = plt.subplot(224)
    ax4.hist(np.abs(y_mod - x_mod), bins=np.linspace(0, 90, 6), orientation='horizontal')
    ax4.set_ylabel(r'y-x ($^{{\circ}}$)')
    ax4.yaxis.tick_right()
    ax4.yaxis.set_label_position("right")
    ax4.set_title('Good cases, z < {}'.format(z_bound))

    # Add lines to the figure
    line1 = plt.Line2D([0.425, 0.53], [0.82, 0.55], color='k', linewidth=2, transform=fig.transFigure)
    line2 = plt.Line2D([0.425, 0.53], [0.87, 0.865], color='k', linewidth=2, transform=fig.transFigure)
    line3 = plt.Line2D([0.425, 0.53], [0.40, 0.13], color='k', linewidth=2, transform=fig.transFigure)
    line4 = plt.Line2D([0.425, 0.53], [0.45, 0.445], color='k', linewidth=2, transform=fig.transFigure)
    
    fig.add_artist(line1)
    fig.add_artist(line2)
    fig.add_artist(line3)
    fig.add_artist(line4)
    
    # Adjust layout and save figure
    fig.legend(bbox_to_anchor=(0.8, 0.98), frameon=False, borderaxespad=0., fontsize=10, bbox_transform=fig.transFigure)

    # Adjust layout and save figure
    plt.subplots_adjust(wspace=0.1, hspace=0.2)
    fig.suptitle('DESI')
    plt.savefig('paper_figures/Figure_M4.pdf', dpi=100, bbox_inches='tight')
    plt.show()

        

def Figure_M5(data):
    """
    Generate and save histograms and scatter plots comparing position angle differences between DESI VLBI and optical data.

    Parameters:
    ----------
    data : tuple
        Contains the following elements:
        - Astrogeo_DESI : DataFrame
            Data from the Astrogeo and DESI catalogues.
        - DESI_xmatches : DataFrame
            Cross-match data between DESI and Astrogeo.
        - z_bound : float
            Redshift boundary for filtering the data.
    """
    # Unpack data
    Astrogeo_DESI, DESI_xmatches, z_bound = data
    
    # Extract position angles and associated errors from the DESI data
    DESI_OPTICAL_PA, DESI_OPTICAL_PA_ERR_ORIGINAL, DESI_VLBI_PA, DESI_VLBI_PA_ERR_ORIGINAL, \
    DESI_PA_DIFF, DESI_PA_DIFF_ERR_ORIGINAL, DESI_SOURCE_Z, DESI_catalogue_Astrogeo_B_matches, \
    DESI_catalogue_Astrogeo_TYPE_matches, DESI_catalogue_Astrogeo_MAG_Z_matches = \
        get_catalogue_parameters.get_DESI_data(Astrogeo_DESI, DESI_xmatches)
    
    # Filter to include only specified types of sources and redshift boundaries
    DESI_good_cases_cut = np.in1d(DESI_catalogue_Astrogeo_TYPE_matches, ['SER', 'EXP', 'DEV'])
    b_cut = DESI_catalogue_Astrogeo_B_matches > 1.3
    z_cut = DESI_SOURCE_Z < z_bound
    
    # Extract and filter position angles and their errors
    y = DESI_VLBI_PA[DESI_good_cases_cut & b_cut & z_cut]
    x = aux_functions.angle_conversion_major_PA(DESI_PA_DIFF[DESI_good_cases_cut & b_cut & z_cut])
    x_err = DESI_PA_DIFF_ERR_ORIGINAL[DESI_good_cases_cut & b_cut & z_cut]
    y_err = DESI_VLBI_PA_ERR_ORIGINAL[DESI_good_cases_cut & b_cut & z_cut]

    # Calculate the correlation coefficient between x and y
    r = np.corrcoef(x, y)
    rho = r[0, 1]

    # Create subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True, layout='constrained')

    # Plot 1: Histogram of position angle differences
    counts, bins, _ = ax1.hist(x, bins=np.linspace(0, 90, 6), histtype='step', fill=True,
                              color='lightcoral', edgecolor='lightcoral', linewidth=3, alpha=0.5)
    for tk in ax1.get_xticklabels():
        tk.set_visible(False)
    ax1.set_xlim([0, 1.2 * np.max(counts)])
    ax1.set_ylabel('Counts per bin')

    # Plot 2: Scatter plot and 2D histogram of position angle differences vs. VLBI PA
    ax2.scatter(x, y, c='b', s=3)
    _, _, _, im = ax2.hist2d(x, y, bins=[np.linspace(0, 90, 6), np.linspace(-180, 180, 2)],
                             cmap='Blues', vmin=0)
    ax2.errorbar(x, y, xerr=x_err, yerr=y_err, fmt='None', ecolor='k', elinewidth=0.5, label='PA err.')
    ax2.axvline(x=0, color='k', lw=0.5)
    ax2.axhline(y=0, color='k', lw=0.5)
    ax2.set_ylabel('Counts per bin')
    for tk in ax2.get_xticklabels():
        tk.set_visible(True)
    for tk in ax2.get_yticklabels():
        tk.set_visible(True)
    ax2.set_xlim([0, 90])
    ax2.set_ylim([-180, 180])
    fig.colorbar(im, ax=ax2, orientation='vertical', label='Counts per bin')
    ax2.set_ylabel(r'Jet PA ($^{{\circ}}$)')
    ax2.set_xlabel(r'$\Delta$PA ($^{{\circ}}$)')
    
    fig.legend(bbox_to_anchor=(0.8, 1.), frameon=False, borderaxespad=0., fontsize=10, bbox_transform=fig.transFigure)
    
    # Add overall title and save figure
    fig.suptitle(r'DESI good cases, z < {}'.format(z_bound))
    plt.savefig('paper_figures/Figure_M5.pdf', dpi=100, bbox_inches='tight')
    plt.show()

    
            

def Figure_M6(data):
    """
    Generate and save histograms comparing position angles for DESI data with various cuts applied.

    Parameters:
    ----------
    data : tuple
        Contains the following elements:
        - Astrogeo_DESI : DataFrame
            Data from the Astrogeo and DESI catalogues.
        - DESI_xmatches : DataFrame
            Cross-match data between DESI and Astrogeo.
        - max_tol_err : float
            Maximum tolerance for error.
        - max_tol_Z_array : array-like
            Array of maximum tolerance values for redshift.
        - N_bins : int
            Number of bins to use in the histograms.
    """
    # Unpack data
    Astrogeo_DESI, DESI_xmatches, max_tol_err, max_tol_Z_array, N_bins = data
    
    # Extract position angles and associated errors from the DESI data
    DESI_OPTICAL_PA, DESI_OPTICAL_PA_ERR_ORIGINAL, DESI_VLBI_PA, DESI_VLBI_PA_ERR_ORIGINAL, \
    DESI_PA_DIFF, DESI_PA_DIFF_ERR_ORIGINAL, DESI_SOURCE_Z, DESI_catalogue_Astrogeo_B_matches, \
    DESI_catalogue_Astrogeo_TYPE_matches, DESI_catalogue_Astrogeo_MAG_Z_matches = \
        get_catalogue_parameters.get_DESI_data(Astrogeo_DESI, DESI_xmatches)
    
    # Filter to include only specified types of sources
    DESI_good_cases_cut = np.in1d(DESI_catalogue_Astrogeo_TYPE_matches, ['SER', 'EXP', 'DEV'])
    
    # Define survey name and output file name
    Survey_Name = 'DESI LS'
    Png_Name = 'Figure_M6.pdf'
    
    # Apply a brightness cut
    b_cut = DESI_catalogue_Astrogeo_B_matches > 1.3
    
    # Generate histograms with specified parameters
    final_p_values_5_bins, final_p_values_2_bins = main_functions.Histograms(
        DESI_VLBI_PA[b_cut], 
        DESI_VLBI_PA_ERR_ORIGINAL[b_cut], 
        DESI_OPTICAL_PA[b_cut], 
        DESI_OPTICAL_PA_ERR_ORIGINAL[b_cut], 
        DESI_PA_DIFF[b_cut], 
        DESI_PA_DIFF_ERR_ORIGINAL[b_cut], 
        DESI_SOURCE_Z[b_cut], 
        max_tol_err, 
        max_tol_Z_array, 
        Survey_Name, 
        Png_Name, 
        N_bins, 
        False, 
        good_cases_cut=DESI_good_cases_cut[b_cut]
    )

    


def Figure_M7(data):
    """
    Generate and save a histogram and error plot of position angle differences versus redshift bins for DESI data.

    Parameters:
    ----------
    data : tuple
        Contains the following elements:
        - Astrogeo_DESI : DataFrame
            Data from the Astrogeo catalogue.
        - DESI_xmatches : DataFrame
            Cross-match data between DESI and Astrogeo.
        - z_tol : float
            Redshift cutoff value.
    """
    # Unpack data
    Astrogeo_DESI, DESI_xmatches, z_tol = data
    
    # Extract relevant arrays from data
    DESI_SOURCE_Z = np.array(Astrogeo_DESI.Z)
    DESI_catalogue_Astrogeo_TYPE_matches = np.array(DESI_xmatches.TYPE)
    DESI_good_cases_cut = np.in1d(DESI_catalogue_Astrogeo_TYPE_matches, ['SER', 'EXP', 'DEV'])
    DESI_catalogue_Astrogeo_B_matches = np.array(DESI_xmatches.b_axis)
    DESI_catalogue_Astrogeo_THETA_J2000_matches = np.array(DESI_xmatches.pos_angle)
    Astrogeo_DESI_catalogue_JET_PA_matches_good = np.array(Astrogeo_DESI.pa)
    
    # Compute position angle differences
    DESI_PA_DIFF = aux_functions.PA_difference(DESI_catalogue_Astrogeo_THETA_J2000_matches, Astrogeo_DESI_catalogue_JET_PA_matches_good)
    
    # Apply redshift cutoff
    z_cut = DESI_SOURCE_Z < z_tol
    
    # Filter valid data
    x = DESI_SOURCE_Z[DESI_good_cases_cut & (np.isnan(DESI_SOURCE_Z) == False) & z_cut]
    y = DESI_PA_DIFF[DESI_good_cases_cut & (np.isnan(DESI_SOURCE_Z) == False) & z_cut]
    
    # Define bin edges
    n = 10
    n_bins = np.linspace(0, 1, n + 1)

    # Bin the data
    y_bins = np.zeros(len(n_bins) - 1, dtype='object')
    y_means = np.zeros(len(n_bins) - 1)
    for i in range(len(n_bins) - 1):
        y_i = y[(n_bins[i] < x) & (x < n_bins[i + 1])]
        y_bins[i] = y_i
        y_means[i] = np.mean(y_i)

    # Compute statistics for each bin
    mean_stat = binned_statistic(x, y, statistic='mean', bins=n_bins, range=(0, np.nanmax(x)))
    counts_stat = binned_statistic(x, y, statistic='count', bins=n_bins, range=(0, np.nanmax(x)))
    median_stat = binned_statistic(x, y, statistic='median', bins=n_bins, range=(0, np.nanmax(x)))
    std_stat = binned_statistic(x, y, statistic='std', bins=n_bins, range=(0, np.nanmax(x)))
    upp_pctl_stat = binned_statistic(x, y, statistic=lambda y: np.percentile(y, 84), bins=n_bins, range=(0, np.nanmax(x)))
    lower_pctl_stat = binned_statistic(x, y, statistic=lambda y: np.percentile(y, 16), bins=n_bins, range=(0, np.nanmax(x)))

    # Extract statistics for bins
    upp_pctl_bins = upp_pctl_stat.statistic
    lower_pctl_bins = lower_pctl_stat.statistic
    mean_bins = mean_stat.statistic
    median_bins = median_stat.statistic
    std_bins = std_stat.statistic
    bin_centers = (mean_stat.bin_edges[1:] + mean_stat.bin_edges[:-1]) / 2

    # Compute confidence intervals
    conf_intervals_arr = np.zeros((n, 2))
    conf_intervals_arr[:] = [16, 84]
    yerrors = np.array([median_bins - lower_pctl_bins, upp_pctl_bins - median_bins])

    # Create plot
    fig, ax1 = plt.subplots()

    # Plot histogram
    ax1.hist(x, bins=n_bins, alpha=0.2)
    ax1.set_ylabel('Counts')
    ax1.set_xlabel('z bin centers')
   
    
    # Create twin axis for error bar plot
    ax2 = ax1.twinx()
    ax2 = plt.gca()
    plt.errorbar(bin_centers, median_bins, yerr=yerrors, xerr=None, fmt='None', elinewidth=1, capsize=3, ecolor='k', alpha=0.5, label='68% confidence interval')
    ax2.scatter(bin_centers, mean_bins, s=7, c='b', label='mean')
    ax2.scatter(bin_centers, median_bins, s=7, c='r', label='median')
    ax2.set_ylim([0, 90])
    ax2.axhline(y=45, ls='--', color='k', alpha=0.5)
    ax2.set_ylabel(r'$\Delta$PA ($^{{\circ}}$)')
    
    # Add legend and title
    plt.legend(frameon=False, bbox_to_anchor=(0.9, 0.98), ncols=3, fontsize=8)
    plt.title('DESI good cases, z < {}'.format(z_tol))
    
    # Save and show plot
    plt.savefig('paper_figures/Figure_M7.pdf', dpi=100)
    plt.show()

    

def Figure_M8(data):
    """
    Generate and save histograms for combined survey data, including VLBI and optical position angles and their differences.

    Parameters:
    ----------
    data : tuple
        Contains the following elements:
        - Surveys_Data : DataFrame
            Data from various surveys.
        - Surveys_Good_Data : DataFrame
            Data from the good quality surveys.
        - max_tol_err : float
            Maximum tolerance for error.
        - max_tol_Z_array : array-like
            Array of tolerance values for redshift.
        - N_bins : int
            Number of bins for histograms.
    """
    # Unpack data
    Surveys_Data, Surveys_Good_Data, max_tol_err, max_tol_Z_array, N_bins = data
    
    # Extract combined data
    (Combined_VLBI_PA, Combined_VLBI_PA_ERR_ORIGINAL, Combined_Good_VLBI_PA, Combined_Good_VLBI_PA_ERR_ORIGINAL, 
     Combined_OPTICAL_PA, Combined_OPTICAL_PA_ERR_ORIGINAL, Combined_Good_OPTICAL_PA, Combined_Good_OPTICAL_PA_ERR_ORIGINAL, 
     Combined_PA_DIFF, Combined_PA_DIFF_ERR_ORIGINAL, Combined_Good_PA_DIFF, Combined_Good_PA_DIFF_ERR_ORIGINAL, 
     Combined_SOURCE_Z, Combined_Good_SOURCE_Z) = get_catalogue_parameters.get_Combined_data(Surveys_Data, Surveys_Good_Data)
    
    # Define survey name and output file name
    Survey_Name = 'Combined'
    Png_Name = 'Figure_M8.pdf'
    
    # Generate histograms
    final_p_values_5_bins, final_p_values_2_bins = main_functions.Histograms(
        Combined_VLBI_PA,
        Combined_VLBI_PA_ERR_ORIGINAL,
        Combined_OPTICAL_PA,
        Combined_OPTICAL_PA_ERR_ORIGINAL,
        Combined_PA_DIFF,
        Combined_PA_DIFF_ERR_ORIGINAL,
        Combined_SOURCE_Z,
        max_tol_err,
        max_tol_Z_array,
        Survey_Name,
        Png_Name,
        N_bins,
        True,
        VLBI_PA_GOOD=Combined_Good_VLBI_PA,
        VLBI_PA_GOOD_ERR=Combined_Good_VLBI_PA_ERR_ORIGINAL,
        OPTICAL_PA_GOOD=Combined_Good_OPTICAL_PA,
        OPTICAL_PA_GOOD_ERR=Combined_Good_OPTICAL_PA_ERR_ORIGINAL,
        PA_DIFF_GOOD=Combined_Good_PA_DIFF,
        PA_DIFF_GOOD_ERR=Combined_Good_PA_DIFF_ERR_ORIGINAL,
        SOURCE_Z_GOOD=Combined_Good_SOURCE_Z
    )

    
        
    
def Figure_M9(data):
    """
    Generate and save histograms for SDSS r band data, including VLBI and optical position angles and their differences.

    Parameters:
    ----------
    data : tuple
        Contains the following elements:
        - Astrogeo_SDSS : DataFrame
            Data from SDSS.
        - SDSS_xmatches : DataFrame
            Cross-match data for SDSS.
        - max_tol_err : float
            Maximum tolerance for error.
        - max_tol_Z_array : array-like
            Array of tolerance values for redshift.
        - N_bins : int
            Number of bins for histograms.
    """
    # Unpack data
    Astrogeo_SDSS, SDSS_xmatches, max_tol_err, max_tol_Z_array, N_bins = data
    
    # Extract SDSS data
    (SDSS_OPTICAL_PA, SDSS_OPTICAL_PA_ERR_ORIGINAL, SDSS_VLBI_PA, SDSS_VLBI_PA_ERR_ORIGINAL, 
     SDSS_PA_DIFF, SDSS_PA_DIFF_ERR_ORIGINAL, SDSS_SOURCE_Z, SDSS_catalogue_Astrogeo_TYPE_matches, 
     SDSS_catalogue_Astrogeo_B_matches) = get_catalogue_parameters.get_SDSS_data(Astrogeo_SDSS, SDSS_xmatches)
    
    # Apply filters
    SDSS_good_cases_cut = np.array(SDSS_catalogue_Astrogeo_TYPE_matches) == 3
    b_cut = SDSS_catalogue_Astrogeo_B_matches > 1.3

    # Define survey name and output file name
    Survey_Name = 'SDSS r band'
    Png_Name = 'Figure_M9.pdf'

    # Generate histograms
    final_p_values_5_bins, final_p_values_2_bins = main_functions.Histograms(
        SDSS_VLBI_PA[b_cut],
        SDSS_VLBI_PA_ERR_ORIGINAL[b_cut],
        SDSS_OPTICAL_PA[b_cut],
        SDSS_OPTICAL_PA_ERR_ORIGINAL[b_cut],
        SDSS_PA_DIFF[b_cut],
        SDSS_PA_DIFF_ERR_ORIGINAL[b_cut],
        SDSS_SOURCE_Z[b_cut],
        max_tol_err,
        max_tol_Z_array,
        Survey_Name,
        Png_Name,
        N_bins,
        False,
        good_cases_cut=SDSS_good_cases_cut[b_cut]
    )

    


def Figure_M10(data):
    """
    Generate and save histograms for DES data, including VLBI and optical position angles and their differences.

    Parameters:
    ----------
    data : tuple
        Contains the following elements:
        - Astrogeo_DES : DataFrame
            Data from DES.
        - DES_xmatches : DataFrame
            Cross-match data for DES.
        - max_tol_err : float
            Maximum tolerance for error.
        - max_tol_Z_array : array-like
            Array of tolerance values for redshift.
        - N_bins : int
            Number of bins for histograms.
    """
    # Unpack data
    Astrogeo_DES, DES_xmatches, max_tol_err, max_tol_Z_array, N_bins = data
    
    # Extract DES data
    (DES_OPTICAL_PA, DES_OPTICAL_PA_ERR_ORIGINAL, DES_VLBI_PA, DES_VLBI_PA_ERR_ORIGINAL, 
     DES_PA_DIFF, DES_PA_DIFF_ERR_ORIGINAL, DES_SOURCE_Z, 
     DES_catalogue_Astrogeo_EXTENDED_CLASS_COADD_matches, DES_catalogue_Astrogeo_B_matches) = get_catalogue_parameters.get_DES_data(Astrogeo_DES, DES_xmatches)
    
    # Apply filters
    DES_good_cases_cut = np.in1d(DES_catalogue_Astrogeo_EXTENDED_CLASS_COADD_matches, [2, 3])
    b_cut = DES_catalogue_Astrogeo_B_matches > 1.3
    
    # Define survey name and output file name
    Survey_Name = 'DES'
    Png_Name = 'Figure_M10.pdf'    
    
    # Generate histograms
    final_p_values_5_bins, final_p_values_2_bins = main_functions.Histograms(
        DES_VLBI_PA[b_cut],
        DES_VLBI_PA_ERR_ORIGINAL[b_cut],
        DES_OPTICAL_PA[b_cut],
        DES_OPTICAL_PA_ERR_ORIGINAL[b_cut],
        DES_PA_DIFF[b_cut],
        DES_PA_DIFF_ERR_ORIGINAL[b_cut],
        DES_SOURCE_Z[b_cut],
        max_tol_err,
        max_tol_Z_array,
        Survey_Name,
        Png_Name,
        N_bins,
        False,
        good_cases_cut=DES_good_cases_cut[b_cut]
    )

    
    
    
def Figure_M11(data):
    """
    Generate and save histograms for SkyMapper data, including VLBI and optical position angles and their differences.

    Parameters:
    ----------
    data : tuple
        Contains the following elements:
        - Astrogeo_Skymapper : DataFrame
            Data from SkyMapper.
        - Skymapper_xmatches : DataFrame
            Cross-match data for SkyMapper.
        - max_tol_err : float
            Maximum tolerance for error.
        - max_tol_Z_array : array-like
            Array of tolerance values for redshift.
        - N_bins : int
            Number of bins for histograms.
    """
    # Unpack data
    Astrogeo_Skymapper, Skymapper_xmatches, max_tol_err, max_tol_Z_array, N_bins = data
    
    # Extract SkyMapper data
    (Skymapper_OPTICAL_PA, Skymapper_OPTICAL_PA_ERR_ORIGINAL, Skymapper_VLBI_PA, 
     Skymapper_VLBI_PA_ERR_ORIGINAL, Skymapper_PA_DIFF, Skymapper_PA_DIFF_ERR_ORIGINAL, 
     Skymapper_SOURCE_Z, Skymapper_Astrogeo_TYPE_matches, Skymapper_Astrogeo_B_matches) = get_catalogue_parameters.get_Skymapper_data(Astrogeo_Skymapper, Skymapper_xmatches)
    
    # Apply filters
    b_filter = 2.0
    stellar_index_filter = 0.5
    Skymapper_good_cases_cut = (Skymapper_Astrogeo_B_matches > b_filter) & (Skymapper_Astrogeo_TYPE_matches < stellar_index_filter)
    b_cut = Skymapper_Astrogeo_B_matches > 1.3
    
    # Define survey name and output file name
    Survey_Name = 'SkyMapper'
    Png_Name = 'Figure_M11.pdf'
    
    # Generate histograms
    final_p_values_5_bins, final_p_values_2_bins = main_functions.Histograms(
        Skymapper_VLBI_PA[b_cut],
        Skymapper_VLBI_PA_ERR_ORIGINAL[b_cut],
        Skymapper_OPTICAL_PA[b_cut],
        Skymapper_OPTICAL_PA_ERR_ORIGINAL[b_cut],
        Skymapper_PA_DIFF[b_cut],
        Skymapper_PA_DIFF_ERR_ORIGINAL[b_cut],
        Skymapper_SOURCE_Z[b_cut],
        max_tol_err,
        max_tol_Z_array,
        Survey_Name,
        Png_Name,
        N_bins,
        False,
        good_cases_cut=Skymapper_good_cases_cut[b_cut]
    )

    
    
    
