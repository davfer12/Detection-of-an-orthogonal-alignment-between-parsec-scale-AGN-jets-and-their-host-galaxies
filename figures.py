from libraries import *
import data_loading_and_xmatching
import aux_functions
import main_functions

def Figure_2(data):
    Astrogeo_catalogue,DESI_xmatches,SDSS_sampled,DES_xmatches,KIDS_sampled = data
    
  
    Astrogeo_RA_good = np.array(Astrogeo_catalogue.RA_deg)
    Astrogeo_DEC_good = np.array(Astrogeo_catalogue.DEC_deg)
    Astrogeo_RA_good[Astrogeo_RA_good > 180] = Astrogeo_RA_good[Astrogeo_RA_good > 180] -360
    Astrogeo_DEC_good[Astrogeo_DEC_good > 180] = Astrogeo_DEC_good[Astrogeo_DEC_good > 180] -360
    Astrogeo_RA_good = Astrogeo_RA_good*np.pi/180
    Astrogeo_DEC_good = Astrogeo_DEC_good*np.pi/180
    
    
    DESI_RA,DESI_DEC = np.array(DESI_xmatches.RA),np.array(DESI_xmatches.DEC)
    DESI_RA[DESI_RA > 180] = DESI_RA[DESI_RA > 180] - 360
    DESI_RA = DESI_RA*np.pi/180
    DESI_DEC[DESI_DEC > 180] = DESI_DEC[DESI_DEC > 180] - 360
    DESI_DEC = DESI_DEC*np.pi/180
    
    
    SDSS_RA = np.array(SDSS_sampled.ra) #[SDSS_random_index]
    SDSS_DEC = np.array(SDSS_sampled.dec) #[SDSS_random_index]
    SDSS_RA[SDSS_RA > 180] = SDSS_RA[SDSS_RA > 180] - 360
    SDSS_RA = SDSS_RA*np.pi/180
    SDSS_DEC[SDSS_DEC > 180] = SDSS_DEC[SDSS_DEC > 180] - 360
    SDSS_DEC = SDSS_DEC*np.pi/180
    
    DES_RA  = DES_xmatches.DES_full_catalogue_Astrogeo_matches_RA
    DES_DEC = DES_xmatches.DES_full_catalogue_Astrogeo_matches_DEC
    DES_RA[DES_RA > 180] = DES_RA[DES_RA > 180] - 360
    DES_RA = DES_RA*np.pi/180
    DES_DEC[DES_DEC > 180] = DES_DEC[DES_DEC > 180] - 360
    DES_DEC = DES_DEC*np.pi/180

    size = 20000
    Skymapper_RA = np.random.uniform(-180,180,size)
    Skymapper_DEC = np.random.uniform(-90,0,size)

    
    Skymapper_RA[Skymapper_RA > 180] = Skymapper_RA[Skymapper_RA > 180] - 360
    Skymapper_DEC[Skymapper_DEC > 180] = Skymapper_DEC[Skymapper_DEC > 180] - 360
    Skymapper_RA = Skymapper_RA*np.pi/180
    Skymapper_DEC = Skymapper_DEC*np.pi/180
    
    KIDS_RA = np.array(KIDS_sampled.RAJ2000)
    KIDS_DEC = np.array(KIDS_sampled.DECJ2000)
    KIDS_RA[KIDS_RA > 180] = KIDS_RA[KIDS_RA > 180] - 360
    KIDS_DEC[KIDS_DEC > 180] = KIDS_DEC[KIDS_DEC > 180] - 360
    KIDS_RA = KIDS_RA*np.pi/180
    KIDS_DEC = KIDS_DEC*np.pi/180
    
    
    plt.figure()
    plt.subplot(111,projection='mollweide')
    plt.scatter(Skymapper_RA,Skymapper_DEC,c='navajowhite',s=17,marker='s',alpha=0.5,label='SkyMapper')
    plt.scatter(DESI_RA,DESI_DEC,c='salmon',s=17,marker='s',alpha=0.5,label='DESI LS')
    plt.scatter(SDSS_RA,SDSS_DEC,c='red',s=17,marker='s',alpha=0.5,label='SDSS')
    plt.scatter(DES_RA,DES_DEC,c='darkorange',s=17,marker='s',alpha=0.5,label='DES')
    plt.scatter(KIDS_RA,KIDS_DEC,c='maroon',s=17,marker='s',alpha=0.5,label='KiDS')
    plt.scatter(Astrogeo_RA_good,Astrogeo_DEC_good,c='cyan',s=0.5,marker='.',alpha=1,label='VLBI')
    lgnd = plt.legend(bbox_to_anchor=(0.85,1.3),ncol=3)
    for lh in lgnd.legendHandles:
        lh._sizes = [30]  
        lh.set_alpha(1)   
    plt.xlabel('RA')
    plt.ylabel('DEC')
    ax = plt.gca()
    ax.set_xticklabels([])
    plt.grid(True)
    plt.savefig('Paper_images/Footprints_HD.png',dpi=300,bbox_inches='tight')
    plt.show()
    


def Figure_3(data):
    Astrogeo_DESI,DESI_xmatches,max_tol_err,max_tol_Z_array,N_bins = data     
    
    DESI_OPTICAL_PA,DESI_OPTICAL_PA_ERR_ORIGINAL,DESI_VLBI_PA,DESI_VLBI_PA_ERR_ORIGINAL,DESI_PA_DIFF,DESI_PA_DIFF_ERR_ORIGINAL,DESI_SOURCE_Z,DESI_catalogue_Astrogeo_B_matches,DESI_catalogue_Astrogeo_TYPE_matches,DESI_catalogue_Astrogeo_MAG_Z_matches  = data_loading_and_xmatching.get_DESI_data(Astrogeo_DESI,DESI_xmatches)
 
    DESI_good_cases_cut = np.in1d(DESI_catalogue_Astrogeo_TYPE_matches,['SER','EXP','DEV'])
    
    Survey_Name = 'DESI LS'
    Png_Name = 'DESI_Final'
    
    
    b_cut = DESI_catalogue_Astrogeo_B_matches > 1.3
    final_p_values_5_bins,final_p_values_2_bins = main_functions.Histograms(DESI_VLBI_PA[b_cut],DESI_VLBI_PA_ERR_ORIGINAL[b_cut],DESI_OPTICAL_PA[b_cut],DESI_OPTICAL_PA_ERR_ORIGINAL[b_cut],DESI_PA_DIFF[b_cut],DESI_PA_DIFF_ERR_ORIGINAL[b_cut],DESI_SOURCE_Z[b_cut],max_tol_err,max_tol_Z_array,Survey_Name,Png_Name,N_bins,False,good_cases_cut=DESI_good_cases_cut[b_cut])
    
    
    
def Figure_4(data):
    Eagle_sim = data
    u = Eagle_sim['u_nodust']   
    g = Eagle_sim['g_nodust']
    r = Eagle_sim['r_nodust']
    i = Eagle_sim['i_nodust']
    z = Eagle_sim['z_nodust']
    Eagle_sim_ellipticals = Eagle_sim[(u-r) > 2.2]
    Eagle_sim_spirals = Eagle_sim[(u-r) < 2.2]
    
    r_PA = np.array(Eagle_sim_ellipticals['pa_r'])
    star_Spin_X = np.array(Eagle_sim_ellipticals['stars_Spin_x'])
    star_Spin_Y = np.array(Eagle_sim_ellipticals['stars_Spin_y'])
    star_Spin_Z = np.array(Eagle_sim_ellipticals['stars_Spin_z'])

    star_Spins = np.array([star_Spin_X,star_Spin_Y,star_Spin_Z]).T

    fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(15,10))

    print('Computing with original data...')

    jet_PA = np.arctan(star_Spin_Y/star_Spin_X)*180/np.pi

    jet_PA[jet_PA < 0] = jet_PA[jet_PA < 0] + 180
    jet_PA =  jet_PA - 90

    delta_PA = np.abs(r_PA - jet_PA)
    delta_PA[delta_PA > 90] = 180 - delta_PA[delta_PA > 90]
    delta_PA = 90 - delta_PA



    n_bins = 5
    bins_array = np.linspace(0,90,n_bins+1)
    bins_center = 0.5*(bins_array[1:] + bins_array[:-1])
    counts,bins = np.histogram(delta_PA, bins=bins_array)

    ax1.hist(delta_PA, bins=bins_array)
    ax1.set_xlabel(r'$\Delta$PA')
    ax1.text(0.5,-0.2, "(a)", size=12, ha="center", transform=ax1.transAxes)
    ax1.set_title(r'Eagle simulation elliptical galaxies with no scatter')   
    print('Computing with uniform spin orientations...')

    delta_PAs = []
    for i in range(1000):
        theta = np.random.uniform(-np.pi,np.pi,len(star_Spins))
        phi = np.random.uniform(0,2*np.pi,len(star_Spins))
        r = 5
        star_Spin_X_p = r*np.cos(theta)*np.cos(phi)
        star_Spin_Y_p = r*np.cos(theta)*np.sin(phi)
        star_Spin_Z_p = r*np.sin(theta)

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





    ax2.bar(bins_array[:-1], counts_avg, bins_array[1] - bins_array[0], align='edge')
    ax2.errorbar(bins_center, counts_avg, yerr=counts_std, fmt='None', ecolor='k')
    ax2.set_xlabel(r'$\Delta$PA')
    ax2.text(0.5,-0.2, "(b)", size=12, ha="center", transform=ax2.transAxes)
    ax2.set_title(r'Eagle simulation elliptical galaxies with uniform spin orientations')   




    epsilon = 0.33
    print(r'Computing with gaussian scatter with epsilon = {}...'.format(epsilon))
    theta_lim = epsilon*90
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


    ax3.bar(bins_array[:-1], counts_avg, bins_array[1] - bins_array[0], align='edge')
    ax3.errorbar(bins_center, counts_avg, yerr=counts_std, fmt='None', ecolor='k')
    ax3.set_xlabel(r'$\Delta$PA')
    ax3.text(0.5,-0.2, "(c)", size=12, ha="center", transform=ax3.transAxes)
    ax3.set_title(r'Eagle simulation elliptical galaxies with gaussian scatter with $\epsilon$ = {}'.format(epsilon))



    epsilon = 0.33
    print(r'Computing with uniform scatter with epsilon = {}...'.format(epsilon))
    theta_lim = epsilon*90
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


    ax4.bar(bins_array[:-1], counts_avg, bins_array[1] - bins_array[0], align='edge')
    ax4.errorbar(bins_center, counts_avg, yerr=counts_std, fmt='None', ecolor='k')
    ax4.set_xlabel(r'$\Delta$PA')
    ax4.text(0.5,-0.2, "(d)", size=12, ha="center", transform=ax4.transAxes)
    ax4.set_title(r'Eagle simulation elliptical galaxies with uniform scatter with $\epsilon$ = {}'.format(epsilon))


    plt.subplots_adjust(hspace=0.4)
    plt.savefig('Paper_images/Eagle_sim_random_ellipticals.png',dpi=100, bbox_inches='tight')
    plt.show()
    
    
    
    
def Figure_M1(data):
    Astrogeo_catalogue = data
    plt.hist(np.array(Astrogeo_catalogue.pa),histtype = 'step', fill = True, linewidth=3, alpha = 0.2)
    plt.xlabel(r'PA ($^{\circ}$)',fontsize=15)
    plt.ylabel('Counts',fontsize=15)
    plt.gca().tick_params(axis='both', which='major', labelsize=15)
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.title('VLBI Jets',fontsize=20)
    plt.savefig('Paper_images/jet_pa.png',dpi=200)
    plt.show()
    
    
def Figure_M2(data):
    DESI_xmatches,Skymapper_xmatches = data
    
    DESI_catalogue_Astrogeo_THETA_J2000_matches = np.array(DESI_xmatches.pos_angle)
    DESI_catalogue_Astrogeo_THETA_J2000_matches = np.fmod(DESI_catalogue_Astrogeo_THETA_J2000_matches+360,180)
    DESI_catalogue_Astrogeo_THETA_J2000_matches[DESI_catalogue_Astrogeo_THETA_J2000_matches > 90] = DESI_catalogue_Astrogeo_THETA_J2000_matches[DESI_catalogue_Astrogeo_THETA_J2000_matches > 90] - 180
    
    Skymapper_Astrogeo_THETA_J2000_matches = np.array(Skymapper_xmatches.PA)
    Skymapper_Astrogeo_THETA_J2000_matches = np.fmod(Skymapper_Astrogeo_THETA_J2000_matches+360,180)
    Skymapper_Astrogeo_THETA_J2000_matches[Skymapper_Astrogeo_THETA_J2000_matches > 90] = Skymapper_Astrogeo_THETA_J2000_matches[Skymapper_Astrogeo_THETA_J2000_matches > 90] - 180 
    
    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(13,5))
    
    ax1.hist(np.array(DESI_catalogue_Astrogeo_THETA_J2000_matches),histtype = 'step', fill = True, linewidth=3, alpha = 0.2)
    ax1.set_xlabel(r'PA ($^{\circ}$)',fontsize=15)
    ax1.set_ylabel('Counts',fontsize=15)
    ax1.tick_params(axis='both', which='major', labelsize=15)
    ax1.text(0.5,-0.2, "(a)", size=12, ha="center", 
         transform=ax1.transAxes)
    ax1.set_title('DESI LS optical cross-matches')
    
    
    ax2.hist(np.array(Skymapper_Astrogeo_THETA_J2000_matches),histtype = 'step', fill = True, linewidth=3, alpha = 0.2)
    ax2.set_xlabel(r'PA ($^{\circ}$)',fontsize=15)
    ax2.set_ylabel('Counts',fontsize=15)
    ax2.tick_params(axis='both', which='major', labelsize=15)
    ax2.text(0.5,-0.2, "(b)", size=12, ha="center", 
         transform=ax2.transAxes)
    ax2.set_title('Skymapper optical cross-matches')
    
    
  
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.savefig('Paper_images/desi_pa.png',dpi=200)
    plt.show()
    
    
def Figure_M3(data):
    Astrogeo_DESI,DESI_xmatches,max_tol_err,max_tol_Z_array,max_mag_tol,N_bins = data     
    DESI_OPTICAL_PA,DESI_OPTICAL_PA_ERR_ORIGINAL,DESI_VLBI_PA,DESI_VLBI_PA_ERR_ORIGINAL,DESI_PA_DIFF,DESI_PA_DIFF_ERR_ORIGINAL,DESI_SOURCE_Z,DESI_catalogue_Astrogeo_B_matches,DESI_catalogue_Astrogeo_TYPE_matches,DESI_catalogue_Astrogeo_MAG_Z_matches = data_loading_and_xmatching.get_DESI_data(Astrogeo_DESI,DESI_xmatches)
    
    DESI_good_cases_cut = np.in1d(DESI_catalogue_Astrogeo_TYPE_matches,['SER','EXP','DEV'])
    
    Survey_Name = 'DESI LS' + ' MAG_Z < ' + str(max_mag_tol)
    Png_Name = 'DESI_Final' + ' MAG_Z_' + str(max_mag_tol)  
    
    
    mag_cut = DESI_catalogue_Astrogeo_MAG_Z_matches < max_mag_tol
    b_cut = DESI_catalogue_Astrogeo_B_matches > 1.3
    final_p_values_5_bins,final_p_values_2_bins = main_functions.Histograms(DESI_VLBI_PA[b_cut & mag_cut],DESI_VLBI_PA_ERR_ORIGINAL[b_cut & mag_cut],DESI_OPTICAL_PA[b_cut & mag_cut],DESI_OPTICAL_PA_ERR_ORIGINAL[b_cut & mag_cut],DESI_PA_DIFF[b_cut & mag_cut],DESI_PA_DIFF_ERR_ORIGINAL[b_cut & mag_cut],DESI_SOURCE_Z[b_cut & mag_cut],max_tol_err,max_tol_Z_array,Survey_Name,Png_Name,N_bins,False,good_cases_cut=DESI_good_cases_cut[b_cut & mag_cut])
                                                                                   


def Figure_M4(data):
    Astrogeo_DESI,DESI_xmatches,max_tol_err,max_tol_Z_array,N_bins = data 
    
    DESI_OPTICAL_PA,DESI_OPTICAL_PA_ERR_ORIGINAL,DESI_VLBI_PA,DESI_VLBI_PA_ERR_ORIGINAL,DESI_PA_DIFF,DESI_PA_DIFF_ERR_ORIGINAL,DESI_SOURCE_Z,DESI_catalogue_Astrogeo_B_matches,DESI_catalogue_Astrogeo_TYPE_matches,DESI_catalogue_Astrogeo_MAG_Z_matches = data_loading_and_xmatching.get_DESI_data(Astrogeo_DESI,DESI_xmatches)
    
    DESI_good_cases_cut = np.in1d(DESI_catalogue_Astrogeo_TYPE_matches,['SER','EXP','DEV'])
    
    Survey_Name = 'DESI LS'
    Png_Name = 'DESI_Final' +'{}_bins'.format(N_bins) 
    
    b_cut = DESI_catalogue_Astrogeo_B_matches > 1.3
    final_p_values_5_bins,final_p_values_2_bins = main_functions.Histograms(DESI_VLBI_PA[b_cut],DESI_VLBI_PA_ERR_ORIGINAL[b_cut],DESI_OPTICAL_PA[b_cut],DESI_OPTICAL_PA_ERR_ORIGINAL[b_cut],DESI_PA_DIFF[b_cut],DESI_PA_DIFF_ERR_ORIGINAL[b_cut],DESI_SOURCE_Z[b_cut],max_tol_err,max_tol_Z_array,Survey_Name,Png_Name,N_bins,False,good_cases_cut=DESI_good_cases_cut[b_cut])
    


def Figure_M5(data):
    Astrogeo_DESI,DESI_xmatches,z_tol = data
    DESI_SOURCE_Z = np.array(Astrogeo_DESI.Z)
    DESI_catalogue_Astrogeo_TYPE_matches = np.array(DESI_xmatches.TYPE)
    DESI_good_cases_cut = np.in1d(DESI_catalogue_Astrogeo_TYPE_matches,['SER','EXP','DEV'])
    DESI_catalogue_Astrogeo_B_matches = np.array(DESI_xmatches.b_axis)
    DESI_catalogue_Astrogeo_THETA_J2000_matches = np.array(DESI_xmatches.pos_angle)
    Astrogeo_DESI_catalogue_JET_PA_matches_good = np.array(Astrogeo_DESI.pa)
    DESI_PA_DIFF = aux_functions.PA_difference(DESI_catalogue_Astrogeo_THETA_J2000_matches,Astrogeo_DESI_catalogue_JET_PA_matches_good)
    b_cut = DESI_catalogue_Astrogeo_B_matches > 1.3
    z_cut = DESI_SOURCE_Z < z_tol
    x = DESI_SOURCE_Z[DESI_good_cases_cut & (np.isnan(DESI_SOURCE_Z) == False) & z_cut & b_cut]
    y = DESI_PA_DIFF[DESI_good_cases_cut & (np.isnan(DESI_SOURCE_Z) == False) & z_cut & b_cut]

    n = 10
    n_bins = np.linspace(0,1,n+1)

    y_bins = np.zeros(len(n_bins)-1,dtype='object')
    y_means = np.zeros(len(n_bins)-1)
    for i in range (len(n_bins) - 1):
        y_i = y[(n_bins[i] < x) & (x < n_bins[i+1])]
        y_bins[i] = y_i
        y_means[i] = np.mean(y_i)

    mean_stat = binned_statistic(x, y, statistic='mean', bins=n_bins, range=(0, np.nanmax(x)))
    counts_stat = binned_statistic(x, y, statistic='count', bins=n_bins, range=(0, np.nanmax(x)))
    median_stat = binned_statistic(x, y, statistic='median', bins=n_bins, range=(0, np.nanmax(x))) 
    std_stat = binned_statistic(x, y, statistic='std', bins=n_bins, range=(0, np.nanmax(x)))
    upp_pctl_stat = binned_statistic(x, y, statistic=lambda y: np.percentile(y, 84), bins=n_bins, range=(0, np.nanmax(x)))
    lower_pctl_stat = binned_statistic(x, y, statistic=lambda y: np.percentile(y, 16), bins=n_bins, range=(0, np.nanmax(x)))

    upp_pctl_bins = upp_pctl_stat.statistic
    lower_pctl_bins = lower_pctl_stat.statistic



    mean_bins = mean_stat.statistic
    median_bins = median_stat.statistic
    std_bins = std_stat.statistic
    bin_centers = (mean_stat.bin_edges[1:]+mean_stat.bin_edges[:-1])/2

    conf_intervals_arr = np.zeros((n,2))
    conf_intervals_arr[:] = [16,84]


    yerrors = np.array([median_bins - lower_pctl_bins,upp_pctl_bins - median_bins])

    fig, ax1 = plt.subplots()

    ax1.hist(x, bins=n_bins,alpha=0.2)
    ax1.set_ylabel('Counts')
    ax1.set_xlabel('z bin centers')
    ax2 = ax1.twinx()
    ax2 = plt.gca()
    plt.errorbar(bin_centers,median_bins,yerr=yerrors,xerr=None,fmt='None',elinewidth=1,capsize=3,ecolor='k',alpha=0.5)
    ax2.scatter(bin_centers,mean_bins,s=7,c='b',label='mean')
    ax2.scatter(bin_centers,median_bins,s=7,c='r',label='median')
    ax2.set_ylim([0,90])
    ax2.axhline(y=45,ls='--',color='k',alpha=0.5)
    ax2.set_ylabel(r'$\Delta$PA ($^{{\circ}}$)')
    plt.legend()
    plt.title('DESI b > 1.3", good cases, z < {}'.format(z_tol))
    plt.savefig('Paper_images/DESI_Final_median_delta_redshift.png',dpi=100)
    plt.show()
    

def Figure_M6(data):
    Surveys_Data,Surveys_Good_Data,max_tol_err,max_tol_Z_array,N_bins = data   
    Combined_VLBI_PA,Combined_VLBI_PA_ERR_ORIGINAL,Combined_Good_VLBI_PA,Combined_Good_VLBI_PA_ERR_ORIGINAL,Combined_OPTICAL_PA,Combined_OPTICAL_PA_ERR_ORIGINAL,Combined_Good_OPTICAL_PA,Combined_Good_OPTICAL_PA_ERR_ORIGINAL,Combined_PA_DIFF,Combined_PA_DIFF_ERR_ORIGINAL,Combined_Good_PA_DIFF,Combined_Good_PA_DIFF_ERR_ORIGINAL,Combined_SOURCE_Z,Combined_Good_SOURCE_Z = data_loading_and_xmatching.get_Combined_data(Surveys_Data,Surveys_Good_Data)
    
    Survey_Name = 'Combined'
    Png_Name = 'Combined_Final'
    
    final_p_values_5_bins,final_p_values_2_bins = main_functions.Histograms(Combined_VLBI_PA,Combined_VLBI_PA_ERR_ORIGINAL,Combined_OPTICAL_PA,Combined_OPTICAL_PA_ERR_ORIGINAL,Combined_PA_DIFF,Combined_PA_DIFF_ERR_ORIGINAL,Combined_SOURCE_Z,max_tol_err,max_tol_Z_array,Survey_Name,Png_Name,N_bins,True,VLBI_PA_GOOD=Combined_Good_VLBI_PA,VLBI_PA_GOOD_ERR=Combined_Good_VLBI_PA_ERR_ORIGINAL,OPTICAL_PA_GOOD=Combined_Good_OPTICAL_PA,OPTICAL_PA_GOOD_ERR=Combined_Good_OPTICAL_PA_ERR_ORIGINAL,PA_DIFF_GOOD=Combined_Good_PA_DIFF,PA_DIFF_GOOD_ERR=Combined_Good_PA_DIFF_ERR_ORIGINAL,SOURCE_Z_GOOD =Combined_Good_SOURCE_Z)
    
        
    
def Figure_M7(data):
    Astrogeo_SDSS,SDSS_xmatches,max_tol_err,max_tol_Z_array,N_bins = data
    
    SDSS_OPTICAL_PA,SDSS_OPTICAL_PA_ERR_ORIGINAL,SDSS_VLBI_PA,SDSS_VLBI_PA_ERR_ORIGINAL,SDSS_PA_DIFF,SDSS_PA_DIFF_ERR_ORIGINAL,SDSS_SOURCE_Z,SDSS_catalogue_Astrogeo_TYPE_matches,SDSS_catalogue_Astrogeo_B_matches = data_loading_and_xmatching.get_SDSS_data(Astrogeo_SDSS,SDSS_xmatches)
    
    
    SDSS_good_cases_cut = np.array(SDSS_catalogue_Astrogeo_TYPE_matches) == 3
    b_cut = SDSS_catalogue_Astrogeo_B_matches > 1.3

    Survey_Name = 'SDSS r band'
    Png_Name = 'SDSS_Final'

    final_p_values_5_bins,final_p_values_2_bins = main_functions.Histograms(SDSS_VLBI_PA[b_cut],SDSS_VLBI_PA_ERR_ORIGINAL[b_cut],SDSS_OPTICAL_PA[b_cut],SDSS_OPTICAL_PA_ERR_ORIGINAL[b_cut],SDSS_PA_DIFF[b_cut],SDSS_PA_DIFF_ERR_ORIGINAL[b_cut],SDSS_SOURCE_Z[b_cut],max_tol_err,max_tol_Z_array,Survey_Name,Png_Name,N_bins,False,good_cases_cut=SDSS_good_cases_cut[b_cut])
    


def Figure_M8(data):
    Astrogeo_DES,DES_xmatches,max_tol_err,max_tol_Z_array,N_bins = data
    DES_OPTICAL_PA,DES_OPTICAL_PA_ERR_ORIGINAL,DES_VLBI_PA,DES_VLBI_PA_ERR_ORIGINAL,DES_PA_DIFF,DES_PA_DIFF_ERR_ORIGINAL,DES_SOURCE_Z,DES_catalogue_Astrogeo_EXTENDED_CLASS_COADD_matches,DES_catalogue_Astrogeo_B_matches = data_loading_and_xmatching.get_DES_data(Astrogeo_DES,DES_xmatches)
    
    DES_good_cases_cut = np.in1d(DES_catalogue_Astrogeo_EXTENDED_CLASS_COADD_matches,[2,3])
    b_cut = DES_catalogue_Astrogeo_B_matches > 1.3
    
    Survey_Name = 'DES'
    Png_Name = 'DES_Final'    
        
    final_p_values_5_bins,final_p_values_2_bins = main_functions.Histograms(DES_VLBI_PA[b_cut],DES_VLBI_PA_ERR_ORIGINAL[b_cut],DES_OPTICAL_PA[b_cut],DES_OPTICAL_PA_ERR_ORIGINAL[b_cut],DES_PA_DIFF[b_cut],DES_PA_DIFF_ERR_ORIGINAL[b_cut],DES_SOURCE_Z[b_cut],max_tol_err,max_tol_Z_array,Survey_Name,Png_Name,N_bins,False,good_cases_cut=DES_good_cases_cut[b_cut])
    
    
    
def Figure_M9(data):
    Astrogeo_Skymapper,Skymapper_xmatches,max_tol_err,max_tol_Z_array,N_bins = data
    
    Skymapper_OPTICAL_PA,Skymapper_OPTICAL_PA_ERR_ORIGINAL,Skymapper_VLBI_PA,Skymapper_VLBI_PA_ERR_ORIGINAL,Skymapper_PA_DIFF,Skymapper_PA_DIFF_ERR_ORIGINAL,Skymapper_SOURCE_Z,Skymapper_catalogue_Astrogeo_TYPE_matches,Skymapper_catalogue_Astrogeo_B_matches = data_loading_and_xmatching.get_Skymapper_data(Astrogeo_Skymapper,Skymapper_xmatches)
    
    b_filter = 2.
    stellar_index_filter = 0.5
    Skymapper_good_cases_cut = (Skymapper_Astrogeo_B_matches > b_filter) & (Skymapper_Astrogeo_TYPE_matches < stellar_index_filter)
    b_cut = Skymapper_catalogue_Astrogeo_B_matches > 1.3
    
    Survey_Name = 'SkyMapper'
    Png_Name = 'Skymapper_Final'
    
    final_p_values_5_bins,final_p_values_2_bins = main_functions.Histograms(Skymapper_VLBI_PA[b_cut],Skymapper_VLBI_PA_ERR_ORIGINAL[b_cut],Skymapper_OPTICAL_PA[b_cut],Skymapper_OPTICAL_PA_ERR_ORIGINAL[b_cut],Skymapper_PA_DIFF[b_cut],Skymapper_PA_DIFF_ERR_ORIGINAL[b_cut],Skymapper_SOURCE_Z[b_cut],max_tol_err,max_tol_Z_array,Survey_Name,Png_Name,N_bins,False,good_cases_cut=Skymapper_good_cases_cut[b_cut])
    
    
    