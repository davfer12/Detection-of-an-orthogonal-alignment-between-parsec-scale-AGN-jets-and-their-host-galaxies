from libraries import *
import aux_functions

def make_combined_data(data,Good):
    
    print('Generating combined dataset...')
    Astrogeo_catalogue,Astrogeo_DESI,DESI_xmatches,Astrogeo_SDSS,SDSS_xmatches,Astrogeo_DES,DES_xmatches,Astrogeo_Skymapper,Skymapper_xmatches,Astrogeo_KIDS,KIDS_xmatches = data
    
    n = len(Astrogeo_catalogue) 
    
    VLBI_indexes = np.zeros(n)
    VLBI_Jet_PAs = np.zeros(n)
    VLBI_Jet_PA_errors = np.zeros(n)
    N_matches = np.zeros(n)

    SDSS_major_PAs = np.zeros(n)
    SDSS_major_PA_errors = np.zeros(n)
    SDSS_PA_diffs = np.zeros(n)
    SDSS_PA_diff_errors = np.zeros(n)
    SDSS_TYPE = np.zeros(n)
    SDSS_BA = np.zeros(n)
    SDSS_BA_errors = np.zeros(n)

    DESI_major_PAs = np.zeros(n)
    DESI_major_PA_errors = np.zeros(n)
    DESI_PA_diffs = np.zeros(n)
    DESI_PA_diff_errors = np.zeros(n)
    DESI_TYPE = np.zeros(n,dtype=object)
    DESI_BA = np.zeros(n)
    DESI_BA_errors = np.zeros(n)


    skymapper_major_PAs = np.zeros(n)
    skymapper_major_PA_errors = np.zeros(n)
    skymapper_PA_diffs = np.zeros(n)
    skymapper_PA_diff_errors = np.zeros(n)
    skymapper_TYPE = np.zeros(n)
    skymapper_A = np.zeros(n)
    skymapper_B = np.zeros(n)
    skymapper_BA = np.zeros(n)


    KIDS_major_PAs = np.zeros(n)
    KIDS_major_PA_errors = np.zeros(n)
    KIDS_PA_diffs = np.zeros(n)
    KIDS_PA_diff_errors = np.zeros(n)
    KIDS_TYPE = np.zeros(n)
    KIDS_A = np.zeros(n)
    KIDS_A_errors = np.zeros(n)
    KIDS_B = np.zeros(n)
    KIDS_B_errors = np.zeros(n)
    KIDS_BA = np.zeros(n)
    KIDS_BA_errors = np.zeros(n)

    DES_major_THETA_J2000_PAs = np.zeros(n)
    DES_major_THETA_IMAGE_PA_errors = np.zeros(n)
    DES_major_THETAWIN_IMAGE_G_PAs = np.zeros(n)
    DES_major_THETAWIN_IMAGE_G_PA_errors = np.zeros(n)
    DES_THETA_J2000_PA_diffs = np.zeros(n)
    DES_THETA_IMAGE_PA_diff_errors = np.zeros(n)
    DES_THETAWIN_IMAGE_G_PA_diffs = np.zeros(n)
    DES_THETAWIN_IMAGE_G_PA_diff_errors = np.zeros(n)
    DES_TYPE = np.zeros(n)
    DES_EXTENDED_CLASS_COADD = np.zeros(n)
    DES_A = np.zeros(n)
    DES_A_errors = np.zeros(n)
    DES_B = np.zeros(n)
    DES_B_errors = np.zeros(n)
    DES_BA = np.zeros(n)
    DES_BA_errors = np.zeros(n)

    Weighted_average_Optical_Major_PAs = np.zeros(n) 
    Weighted_average_Optical_Major_PA_errors = np.zeros(n)

    Weighted_average_Jet_Minor_diffs = np.zeros(n) 
    Weighted_average_Jet_Minor_diffs_errors = np.zeros(n)

    Surveys_present = np.zeros(n)


    b_filter = 1.3
    stellar_index_filter = 0.5
    DES_stellar_index_filter = 0.01
    
    
    if Good:
        i_SDSS_good = (np.array(SDSS_xmatches.type_r) == 3) & (np.array(SDSS_xmatches.modelAB_r)*np.array(SDSS_xmatches.devrad_r) > b_filter)
        i_DESI_good = np.in1d(np.array(DESI_xmatches.TYPE),['SER','EXP','DEV']) & (np.array(DESI_xmatches.b_axis) > b_filter)
        i_skymapper_good = (np.array(Skymapper_xmatches.ClassStar) < stellar_index_filter) & (0.5*np.array(Skymapper_xmatches.b) > b_filter) 
        i_KIDS_good = (np.array(KIDS_xmatches.CLASS_STAR) < stellar_index_filter) & (0.2*np.array(KIDS_xmatches.B_IMAGE) > b_filter)
        i_DES_good = np.in1d(np.array(DES_xmatches.DES_full_catalogue_Astrogeo_EXTENDED_CLASS_COADD_matches),[2,3]) & (0.26*np.array(DES_xmatches.DES_full_catalogue_Astrogeo_B_matches) > b_filter)
   

        SDSS_filter = i_SDSS_good 
        DESI_filter = i_DESI_good 
        skymapper_filter = i_skymapper_good 
        KIDS_filter = i_KIDS_good 
        DES_filter = i_DES_good 
    
    else:
        SDSS_filter = np.array(SDSS_xmatches.modelAB_r)*np.array(SDSS_xmatches.devrad_r) > b_filter
        DESI_filter = np.array(DESI_xmatches.b_axis) > b_filter
        skymapper_filter = 0.5*np.array(Skymapper_xmatches.b) > b_filter
        KIDS_filter = 0.2*np.array(KIDS_xmatches.B_IMAGE) > b_filter
        DES_filter = 0.26*np.array(DES_xmatches.DES_full_catalogue_Astrogeo_B_matches) > b_filter


    for i in np.arange(n):
        temporal_survey_names = []
        temporal_survey_present = 0
        VLBI_indexes[i] = i
        VLBI_Jet_PAs[i] = Astrogeo_catalogue.pa[i]
        VLBI_Jet_PA_errors[i] = Astrogeo_catalogue.pa_err[i]
        Optical_PAs = []
        Optical_PA_errors = []
        Jet_Minor_diffs = []
        Jet_Minor_diffs_errors = []
        SDSS_index = np.where(np.array(Astrogeo_SDSS.index)[SDSS_filter] == i)[0]
        if len(SDSS_index) == 1:
            temporal_survey_names.append('SDSS')
            temporal_survey_present += 2**4
            
            
            SDSS_catalogue_Astrogeo_THETA_J2000_matches = np.fmod(np.array(SDSS_xmatches.modelPhi_r)+360,180)
            SDSS_catalogue_Astrogeo_THETA_J2000_matches[SDSS_catalogue_Astrogeo_THETA_J2000_matches > 90] = SDSS_catalogue_Astrogeo_THETA_J2000_matches[SDSS_catalogue_Astrogeo_THETA_J2000_matches > 90] - 180  
            SDSS_major_PAs[i] = SDSS_catalogue_Astrogeo_THETA_J2000_matches[SDSS_filter][SDSS_index[0]]
            Optical_PAs.append(SDSS_major_PAs[i])
            
            SDSS_catalogue_Astrogeo_THETA_J2000_ERR_ORIGINAL_matches = np.array(SDSS_xmatches.modelPhiErr)
            SDSS_major_PA_errors[i] = SDSS_catalogue_Astrogeo_THETA_J2000_ERR_ORIGINAL_matches[SDSS_filter][SDSS_index[0]]
            Optical_PA_errors.append(SDSS_major_PA_errors[i])
            
            Astrogeo_SDSS_catalogue_JET_PA_matches_good = np.array(Astrogeo_SDSS.pa)
            Astrogeo_SDSS_catalogue_JET_PA_ERR_ORIGINAL_matches_good = np.array(Astrogeo_SDSS.pa_err)
            
            SDSS_catalogue_Astrogeo_PA_DIFFERENCE_matches = aux_functions.PA_difference(SDSS_catalogue_Astrogeo_THETA_J2000_matches,Astrogeo_SDSS_catalogue_JET_PA_matches_good)
            SDSS_catalogue_Astrogeo_PA_DIFFERENCE_ERR_ORIGINAL_matches = np.sqrt(Astrogeo_SDSS_catalogue_JET_PA_ERR_ORIGINAL_matches_good**2+SDSS_catalogue_Astrogeo_THETA_J2000_ERR_ORIGINAL_matches**2)  
            
            SDSS_PA_diffs[i] = SDSS_catalogue_Astrogeo_PA_DIFFERENCE_matches[SDSS_filter][SDSS_index[0]]
            Jet_Minor_diffs.append(SDSS_catalogue_Astrogeo_PA_DIFFERENCE_matches[SDSS_filter][SDSS_index[0]])
            
            SDSS_PA_diff_errors[i] = SDSS_catalogue_Astrogeo_PA_DIFFERENCE_ERR_ORIGINAL_matches[SDSS_filter][SDSS_index[0]]
            Jet_Minor_diffs_errors.append(SDSS_PA_diff_errors[i])
            
            SDSS_TYPE[i] = np.array(SDSS_xmatches.type_r)[SDSS_filter][SDSS_index[0]]
            SDSS_BA[i] = np.array(SDSS_xmatches.modelAB_r)[SDSS_filter][SDSS_index[0]]
            SDSS_BA_errors[i] = np.array(SDSS_xmatches.modelABErr_r)[SDSS_filter][SDSS_index[0]]   
        else:
            SDSS_major_PAs[i] = np.nan
            SDSS_major_PA_errors[i] = np.nan
            SDSS_PA_diffs[i] = np.nan
            SDSS_PA_diff_errors[i] = np.nan
            SDSS_TYPE[i] = np.nan
            SDSS_BA[i] = np.nan
            SDSS_BA_errors[i] = np.nan

        DESI_index = np.where(np.array(Astrogeo_DESI.index)[DESI_filter] == i)[0]
        if len(DESI_index) == 1:
            temporal_survey_names.append('DESI')
            temporal_survey_present += 2**3
            
            DESI_catalogue_Astrogeo_THETA_J2000_matches = np.fmod(np.array(DESI_xmatches.pos_angle)+360,180)
            DESI_catalogue_Astrogeo_THETA_J2000_matches[DESI_catalogue_Astrogeo_THETA_J2000_matches > 90] = DESI_catalogue_Astrogeo_THETA_J2000_matches[DESI_catalogue_Astrogeo_THETA_J2000_matches > 90] - 180
            
            DESI_major_PAs[i] = np.array(DESI_catalogue_Astrogeo_THETA_J2000_matches)[DESI_filter][DESI_index[0]]
            Optical_PAs.append(DESI_major_PAs[i])
            
            DESI_catalogue_Astrogeo_THETA_J2000_ERR_ORIGINAL_matches = np.array(DESI_xmatches.pos_angle_err)
            DESI_error = DESI_catalogue_Astrogeo_THETA_J2000_ERR_ORIGINAL_matches[DESI_filter][DESI_index[0]]
            if np.isnan(DESI_error):
                DESI_error = 360.
            DESI_major_PA_errors[i] = DESI_error
            Optical_PA_errors.append(DESI_major_PA_errors[i])
            
            Astrogeo_DESI_catalogue_JET_PA_matches_good = np.array(Astrogeo_DESI.pa)
            Astrogeo_DESI_catalogue_JET_PA_ERR_ORIGINAL_matches_good = np.array(Astrogeo_DESI.pa_err)
            
            DESI_catalogue_Astrogeo_PA_DIFFERENCE_matches = aux_functions.PA_difference(DESI_catalogue_Astrogeo_THETA_J2000_matches,Astrogeo_DESI_catalogue_JET_PA_matches_good)   
            DESI_catalogue_Astrogeo_PA_DIFFERENCE_ERR_ORIGINAL_matches = np.sqrt(Astrogeo_DESI_catalogue_JET_PA_ERR_ORIGINAL_matches_good**2+DESI_catalogue_Astrogeo_THETA_J2000_ERR_ORIGINAL_matches**2)
            
            DESI_PA_diffs[i] = DESI_catalogue_Astrogeo_PA_DIFFERENCE_matches[DESI_filter][DESI_index[0]]
            DESI_PA_diff_errors[i] = DESI_catalogue_Astrogeo_PA_DIFFERENCE_ERR_ORIGINAL_matches[DESI_filter][DESI_index[0]]
            Jet_Minor_diffs.append(DESI_PA_diffs[i])
            Jet_Minor_diffs_errors.append(DESI_PA_diff_errors[i])
            
            DESI_TYPE[i] = np.array(DESI_xmatches.TYPE)[DESI_filter][DESI_index[0]]
            DESI_BA[i] = np.array(DESI_xmatches.BA_ratio)[DESI_filter][DESI_index[0]]
            DESI_BA_errors[i] = np.array(DESI_xmatches.BA_ratio_err)[DESI_filter][DESI_index[0]]
        else:
            DESI_major_PAs[i] = np.nan
            DESI_major_PA_errors[i] = np.nan
            DESI_PA_diffs[i] = np.nan
            DESI_PA_diff_errors[i] = np.nan
            DESI_TYPE[i] = np.nan
            DESI_BA[i] = np.nan
            DESI_BA_errors[i] = np.nan

        skymapper_index = np.where(np.array(Astrogeo_Skymapper.index)[skymapper_filter] == i)[0] 
        if len(skymapper_index) == 1:
            temporal_survey_names.append('Skymapper')
            temporal_survey_present += 2**2
            
            skymapper_Astrogeo_THETA_J2000_matches = np.fmod(np.array(Skymapper_xmatches.PA)+360,180)
            skymapper_Astrogeo_THETA_J2000_matches[skymapper_Astrogeo_THETA_J2000_matches > 90] = skymapper_Astrogeo_THETA_J2000_matches[skymapper_Astrogeo_THETA_J2000_matches > 90] - 180 
            
            skymapper_major_PAs[i] = skymapper_Astrogeo_THETA_J2000_matches[skymapper_filter][skymapper_index[0]]
            Optical_PAs.append(skymapper_major_PAs[i])
            
            skymapper_Astrogeo_THETA_J2000_ERR_ORIGINAL_matches = np.array(Skymapper_xmatches.e_pa)
            skymapper_error = skymapper_Astrogeo_THETA_J2000_ERR_ORIGINAL_matches[skymapper_filter][skymapper_index[0]]
            skymapper_major_PA_errors[i] = skymapper_error
            Optical_PA_errors.append(skymapper_major_PA_errors[i])
            
                      
            Astrogeo_skymapper_JET_PA_good_matches = np.array(Astrogeo_Skymapper.pa)
            Astrogeo_skymapper_JET_PA_ERR_ORIGINAL_good_matches = np.array(Astrogeo_Skymapper.pa_err)
            
            skymapper_Astrogeo_PA_DIFFERENCE_matches = aux_functions.PA_difference(skymapper_Astrogeo_THETA_J2000_matches,Astrogeo_skymapper_JET_PA_good_matches)
            skymapper_Astrogeo_PA_DIFFERENCE_ERR_ORIGINAL_matches = np.sqrt(Astrogeo_skymapper_JET_PA_ERR_ORIGINAL_good_matches**2+skymapper_Astrogeo_THETA_J2000_ERR_ORIGINAL_matches**2)
            
            skymapper_PA_diffs[i] = np.array(skymapper_Astrogeo_PA_DIFFERENCE_matches)[skymapper_filter][skymapper_index[0]]
            Jet_Minor_diffs.append(skymapper_PA_diffs[i])  
                       
            skymapper_PA_diff_errors[i] = skymapper_Astrogeo_PA_DIFFERENCE_ERR_ORIGINAL_matches[skymapper_filter][skymapper_index[0]]
            Jet_Minor_diffs_errors.append(skymapper_PA_diff_errors[i])
            
            skymapper_A[i] = 0.5*np.array(Skymapper_xmatches.a)[skymapper_filter][skymapper_index[0]]
            skymapper_B[i] = 0.5*np.array(Skymapper_xmatches.b)[skymapper_filter][skymapper_index[0]]
            skymapper_BA[i] = skymapper_B[i]/skymapper_A[i]
            skymapper_TYPE[i] = np.array(Skymapper_xmatches.ClassStar)[skymapper_filter][skymapper_index[0]]

        else:
            skymapper_major_PAs[i] = np.nan
            skymapper_major_PA_errors[i] = np.nan
            skymapper_PA_diffs[i] = np.nan
            skymapper_PA_diff_errors[i] = np.nan
            skymapper_A[i] = np.nan
            skymapper_B[i] = np.nan
            skymapper_BA[i] = np.nan
            skymapper_TYPE[i] = np.nan

        KIDS_index = np.where(np.array(Astrogeo_KIDS.index)[KIDS_filter] == i)[0] 
        if len(KIDS_index) == 1:
            temporal_survey_names.append('KIDS')
            temporal_survey_present += 2**1
            
            KIDS_catalogue_Astrogeo_THETA_J2000_matches = np.fmod(-np.array(data.THETA_J2000)+360,180)
            KIDS_catalogue_Astrogeo_THETA_J2000_matches[KIDS_catalogue_Astrogeo_THETA_J2000_matches > 90] = KIDS_catalogue_Astrogeo_THETA_J2000_matches[KIDS_catalogue_Astrogeo_THETA_J2000_matches > 90] - 180 
            
            KIDS_major_PAs[i] = np.array(KIDS_catalogue_Astrogeo_THETA_J2000_matches)[KIDS_filter][KIDS_index[0]]
            Optical_PAs.append(KIDS_major_PAs[i])
            
            KIDS_Astrogeo_THETA_J2000_ERR_ORIGINAL_matches = np.array(data.ERRTHETA_J2000)
            KIDS_error = KIDS_Astrogeo_THETA_J2000_ERR_ORIGINAL_matches[KIDS_filter][KIDS_index[0]]
            KIDS_major_PA_errors[i] = KIDS_error
            Optical_PA_errors.append(KIDS_major_PA_errors[i])
            
            Astrogeo_KIDS_JET_PA_good_matches = np.array(Astrogeo_KIDS.pa)
            Astrogeo_KIDS_JET_PA_ERR_ORIGINAL_good_matches = np.array(Astrogeo_KIDS.pa_err)
            
            KIDS_Astrogeo_PA_DIFFERENCE_matches = aux_functions.PA_difference(KIDS_Astrogeo_THETA_J2000_matches,Astrogeo_KIDS_JET_PA_good_matches)
            KIDS_Astrogeo_PA_DIFFERENCE_ERR_ORIGINAL_matches = np.sqrt(Astrogeo_KIDS_JET_PA_ERR_ORIGINAL_good_matches**2+KIDS_Astrogeo_THETA_J2000_ERR_ORIGINAL_matches**2)
            
            KIDS_PA_diffs[i] = KIDS_Astrogeo_PA_DIFFERENCE_matches[KIDS_filter][KIDS_index[0]]
            Jet_Minor_diffs.append(KIDS_PA_diffs[i])
            
            KIDS_PA_diff_errors[i] = KIDS_Astrogeo_PA_DIFFERENCE_ERR_ORIGINAL_matches[KIDS_filter][KIDS_index[0]]
            Jet_Minor_diffs_errors.append(KIDS_PA_diff_errors[i])
            
            KIDS_TYPE[i] = np.array(KIDS_catalogue_Astrogeo_CLASS_STAR_matches)[KIDS_filter][KIDS_index[0]]
            KIDS_A[i] = np.array(KIDS_xmatches.A_IMAGE)[KIDS_filter][KIDS_index[0]]
            KIDS_A_errors[i] = np.array(KIDS_xmatches.ERRA_IMAGE)[KIDS_filter][KIDS_index[0]]
            KIDS_B[i] = np.array(KIDS_xmatches.B_IMAGE)[KIDS_filter][KIDS_index[0]]
            KIDS_B_errors[i] = np.array(KIDS_xmatches.ERRB_IMAGE)[KIDS_filter][KIDS_index[0]]
            KIDS_BA[i] = KIDS_B[i]/KIDS_A[i]
            KIDS_BA_errors[i] = np.sqrt((KIDS_B_errors[i]/KIDS_A[i])**2 + (KIDS_A_errors[i]*KIDS_B[i]/KIDS_A[i]**2)**2)
        else:
            KIDS_major_PAs[i] = np.nan
            KIDS_major_PA_errors[i] = np.nan
            KIDS_PA_diffs[i] = np.nan
            KIDS_PA_diff_errors[i] = np.nan
            KIDS_TYPE[i] = np.nan
            KIDS_A[i] = np.nan
            KIDS_A_errors[i] = np.nan
            KIDS_B[i] = np.nan
            KIDS_B_errors[i] = np.nan
            KIDS_BA[i] = np.nan
            KIDS_BA_errors[i] = np.nan

        #Apparently the DES index to Astrogeo are not unique! So I will get more than one!
        DES_index = np.where(np.array(DES_xmatches.DES_full_catalogue_Astrogeo_matches_indexes)[DES_filter] == i)[0]
        if len(DES_index) == 1:
            temporal_survey_names.append('DES')
            temporal_survey_present += 2**0
            
            
            DES_major_THETA_J2000_PAs[i] = np.array(DES_xmatches.DES_full_catalogue_Astrogeo_THETA_J2000_matches)[DES_filter][DES_index[0]]
            Optical_PAs.append(DES_major_THETA_J2000_PAs[i])
            DES_THETA_IMAGE_error = np.array(DES_xmatches.tot_err)[DES_filter][DES_index[0]]
            DES_major_THETA_IMAGE_PA_errors[i] = DES_THETA_IMAGE_error
            Optical_PA_errors.append(DES_major_THETA_IMAGE_PA_errors[i])
            DES_major_THETAWIN_IMAGE_G_PAs[i] = np.array(DES_xmatches.DES_full_catalogue_Astrogeo_THETAWIN_IMAGE_G_matches)[DES_filter][DES_index[0]]
            DES_THETAWIN_IMAGE_G_error = np.array(DES_xmatches.DES_full_catalogue_Astrogeo_ERRTHETAWIN_IMAGE_G_matches)[DES_filter][DES_index[0]]
            DES_major_THETAWIN_IMAGE_G_PA_errors[i] = DES_THETAWIN_IMAGE_G_error
            DES_THETA_J2000_PA_diffs[i] = np.array(DES_xmatches.DES_full_catalogue_Astrogeo_PA_DIFFERENCE_matches)[DES_filter][DES_index[0]]
            Jet_Minor_diffs.append(np.array(DES_xmatches.DES_full_catalogue_Astrogeo_PA_DIFFERENCE_matches)[DES_filter][DES_index[0]])
            DES_THETA_IMAGE_PA_diff_errors[i] = np.sqrt(Astrogeo_catalogue.tot_err[i]**2 + DES_THETA_IMAGE_error**2)
            Jet_Minor_diffs_errors.append(np.sqrt(Astrogeo_catalogue.tot_err[i]**2 + DES_THETA_IMAGE_error**2))
            DES_THETAWIN_IMAGE_G_PA_diffs[i] = np.array(DES_xmatches.DES_full_catalogue_Astrogeo_PA_DIFFERENCE_THETAWIN_IMAGE_G_matches)[DES_filter][DES_index[0]]
            DES_THETAWIN_IMAGE_G_PA_diff_errors[i] = np.sqrt(Astrogeo_catalogue.tot_err[i]**2 + DES_THETAWIN_IMAGE_G_error**2)
            DES_TYPE[i] = np.array(DES_xmatches.DES_full_catalogue_Astrogeo_TYPE_matches)[DES_filter][DES_index[0]]
            DES_EXTENDED_CLASS_COADD[i] = np.array(DES_xmatches.DES_full_catalogue_Astrogeo_EXTENDED_CLASS_COADD_matches)[DES_filter][DES_index[0]]
            DES_A[i] = np.array(DES_xmatches.DES_full_catalogue_Astrogeo_A_matches)[DES_filter][DES_index[0]]
            DES_A_errors[i] = np.array(DES_xmatches.DES_full_catalogue_Astrogeo_A_ERR_matches)[DES_filter][DES_index[0]]
            DES_B[i] = np.array(DES_xmatches.DES_full_catalogue_Astrogeo_B_matches)[DES_filter][DES_index[0]]
            DES_B_errors[i] = np.array(DES_xmatches.DES_full_catalogue_Astrogeo_B_ERR_matches)[DES_filter][DES_index[0]]
            DES_BA[i] = np.array(DES_xmatches.DES_full_catalogue_Astrogeo_BA_matches)[DES_filter][DES_index[0]]
            DES_BA_errors[i] =  np.array(DES_xmatches.DES_full_catalogue_Astrogeo_BA_ERR_matches)[DES_filter][DES_index[0]]             
        else:
            DES_major_THETA_J2000_PAs[i] = np.nan
            DES_major_THETA_IMAGE_PA_errors[i] = np.nan
            DES_major_THETAWIN_IMAGE_G_PAs[i] = np.nan
            DES_major_THETAWIN_IMAGE_G_PA_errors[i] = np.nan
            DES_THETA_J2000_PA_diffs[i] = np.nan
            DES_THETA_IMAGE_PA_diff_errors[i] = np.nan
            DES_THETAWIN_IMAGE_G_PA_diffs[i] = np.nan
            DES_THETAWIN_IMAGE_G_PA_diff_errors[i] = np.nan
            DES_TYPE[i] = np.nan
            DES_EXTENDED_CLASS_COADD[i] = np.nan
            DES_A[i] = np.nan
            DES_A_errors[i] = np.nan
            DES_B[i] = np.nan
            DES_B_errors[i] = np.nan
            DES_BA[i] = np.nan
            DES_BA_errors[i] = np.nan

        N_matches[i] = len(temporal_survey_names)
        Surveys_present[i] = temporal_survey_present
        if len(Jet_Minor_diffs) > 0:
            Jet_Minor_diffs = np.array(Jet_Minor_diffs)
            Jet_Minor_diffs_errors = np.array(Jet_Minor_diffs_errors)
            Weights = Jet_Minor_diffs_errors**(-2)/np.sum(Jet_Minor_diffs_errors**(-2))
            Weighted_average_Optical_PAs = np.sum(Weights*Optical_PAs)
            Weighted_average_Optical_PA_errors = np.sum(Weights*Optical_PA_errors)
            Weighted_average_Optical_Major_PAs[i] = Weighted_average_Optical_PAs
            Weighted_average_Optical_Major_PA_errors[i] = Weighted_average_Optical_PA_errors
            Weighted_average_Jet_Minor_diff = np.sum(Weights*Jet_Minor_diffs)
            Weighted_average_Jet_Minor_diff_error = np.sum(Jet_Minor_diffs_errors**(-2))**(-1/2)
            Weighted_average_Jet_Minor_diffs[i] = Weighted_average_Jet_Minor_diff
            Weighted_average_Jet_Minor_diffs_errors[i] = Weighted_average_Jet_Minor_diff_error
        else:
            Weighted_average_Optical_Major_PAs[i] = np.nan
            Weighted_average_Optical_Major_PA_errors[i] = np.nan
            Weighted_average_Jet_Minor_diffs[i] = np.nan
            Weighted_average_Jet_Minor_diffs_errors[i] = np.nan
        print('{:.0f}%'.format(100*i/n), end='\r')
    
    data = {
        'VLBI_indexes':VLBI_indexes,
        'VLBI_Jet_PAs':VLBI_Jet_PAs,
        'VLBI_Jet_PA_errors':VLBI_Jet_PA_errors,
        'N_matches':N_matches,
        'Surveys_present':Surveys_present,
        'SDSS_major_PAs':SDSS_major_PAs,
        'SDSS_major_PA_errors':SDSS_major_PA_errors,
        'SDSS_PA_diffs':SDSS_PA_diffs,
        'SDSS_PA_diff_errors':SDSS_PA_diff_errors,
        'SDSS_TYPE':SDSS_TYPE,
        'SDSS_BA':SDSS_BA,
        'SDSS_BA_errors':SDSS_BA_errors,
        'DESI_major_PAs':DESI_major_PAs,
        'DESI_major_PA_errors':DESI_major_PA_errors,
        'DESI_PA_diffs':DESI_PA_diffs,
        'DESI_PA_diff_errors':DESI_PA_diff_errors,
        'DESI_TYPE':DESI_TYPE,
        'DESI_BA':DESI_BA,
        'DESI_BA_errors':DESI_BA_errors,
        'skymapper_major_PAs':skymapper_major_PAs,
        'skymapper_major_PA_errors':skymapper_major_PA_errors,
        'skymapper_PA_diffs':skymapper_PA_diffs,
        'skymapper_PA_diff_errors':skymapper_PA_diff_errors,
        'skymapper_TYPE':skymapper_TYPE,
        'skymapper_A':skymapper_A,
        'skymapper_B':skymapper_B,
        'skymapper_BA':skymapper_BA,
        'KIDS_major_PAs':KIDS_major_PAs,
        'KIDS_major_PA_errors':KIDS_major_PA_errors,
        'KIDS_PA_diffs':KIDS_PA_diffs,
        'KIDS_PA_diff_errors':KIDS_PA_diff_errors,
        'KIDS_TYPE':KIDS_TYPE,
        'KIDS_A':KIDS_A,
        'KIDS_A_errors':KIDS_A_errors,
        'KIDS_B':KIDS_B,
        'KIDS_B_errors':KIDS_B_errors,
        'KIDS_BA':KIDS_BA,
        'KIDS_BA_errors':KIDS_BA_errors,
        'DES_major_THETA_J2000_PAs':DES_major_THETA_J2000_PAs,
        'DES_major_THETA_IMAGE_PA_errors':DES_major_THETA_IMAGE_PA_errors,
        'DES_major_THETAWIN_IMAGE_G_PAs':DES_major_THETAWIN_IMAGE_G_PAs,
        'DES_major_THETAWIN_IMAGE_G_PA_errors':DES_major_THETAWIN_IMAGE_G_PA_errors,
        'DES_THETA_J2000_PA_diffs':DES_THETA_J2000_PA_diffs,
        'DES_THETA_IMAGE_PA_diff_errors':DES_THETA_IMAGE_PA_diff_errors,
        'DES_THETAWIN_IMAGE_G_PA_diffs':DES_THETAWIN_IMAGE_G_PA_diffs,
        'DES_THETAWIN_IMAGE_G_PA_diff_errors':DES_THETAWIN_IMAGE_G_PA_diff_errors,
        'DES_TYPE':DES_TYPE,
        'DES_EXTENDED_CLASS_COADD':DES_EXTENDED_CLASS_COADD,
        'DES_A':DES_A,
        'DES_A_errors':DES_A_errors,
        'DES_B':DES_B,
        'DES_B_errors':DES_B_errors,
        'DES_BA':DES_BA,
        'DES_BA_errors':DES_BA_errors,
        'Weighted_average_Optical_Major_PAs':Weighted_average_Optical_Major_PAs,
        'Weighted_average_Optical_Major_PA_errors':Weighted_average_Optical_Major_PA_errors,
        'Weighted_average_Jet_Minor_diffs':Weighted_average_Jet_Minor_diffs,
        'Weighted_average_Jet_Minor_diffs_errors':Weighted_average_Jet_Minor_diffs_errors
    }

    Surveys_Data = pd.DataFrame(data)
    if Good:
        Z_Surveys_Good_Data = np.array(Astrogeo_catalogue.Z)
        Z_Surveys_Good_Data[np.isnan(Surveys_Data.Weighted_average_Jet_Minor_diffs)] = np.nan
        Z = Z_Surveys_Good_Data
    else:
        Z = np.array(Astrogeo_catalogue.Z) 
    
    Surveys_Data['Z'] = Z
    return Surveys_Data
    #Surveys_Data.to_csv('Surveys_Data.csv',index=False)

    
    
def divide_data_into_panels(VLBI_PA_ERR,OPTICAL_PA_ERR,PA_DIFF,PA_DIFF_ERR,SOURCE_Z,max_tol_Z_array,max_tol_err,Combined,good_cases_cut,PA_DIFF_GOOD,PA_DIFF_GOOD_ERR,SOURCE_Z_GOOD):
    

    
    pa_jet_angle_cut = VLBI_PA_ERR<max_tol_err
    pa_optical_angle_cut = OPTICAL_PA_ERR<max_tol_err
    pa_diff_angle_cut = pa_jet_angle_cut & pa_optical_angle_cut

    
    PA_diffs_complete,PA_diffs_err_complete = PA_DIFF,PA_DIFF_ERR
    PA_diffs_angle_cut,PA_diffs_angle_cut_err = PA_DIFF[pa_diff_angle_cut],PA_DIFF_ERR[pa_diff_angle_cut]
    if Combined:
        PA_diffs_good_cases,PA_diffs_good_cases_err = PA_DIFF_GOOD,PA_DIFF_GOOD_ERR
    else:
        PA_diffs_good_cases,PA_diffs_good_cases_err = PA_DIFF[good_cases_cut],PA_DIFF_ERR[good_cases_cut]
    
    max_size = len(PA_diffs_complete)
    size_1 = max_size
    size_2 = len(PA_diffs_angle_cut)
    size_3 = len(PA_diffs_good_cases)
    
    n_rows = len(max_tol_Z_array) + 1

  
    data_array = np.empty((n_rows,3,max_size))
    data_err_array = np.empty((n_rows,3,max_size))
    
    data_array[:,:,:] = np.nan
    data_err_array[:,:,:] = np.nan

    
    
    data_array[0,0,:size_1],data_array[0,1,:size_2],data_array[0,2,:size_3] = PA_diffs_complete,PA_diffs_angle_cut,PA_diffs_good_cases
    data_err_array[0,0,:size_1],data_err_array[0,1,:size_2],data_err_array[0,2,:size_3] = PA_diffs_err_complete,PA_diffs_angle_cut_err,PA_diffs_good_cases_err
    
    
    
    
    for i,max_tol_Z in enumerate(max_tol_Z_array):
        Z_filter_i = (SOURCE_Z < max_tol_Z) & (SOURCE_Z > 0)
        jet_err_Z_filter_i = Z_filter_i & (pa_jet_angle_cut)
        optical_err_Z_filter_i = Z_filter_i & (pa_optical_angle_cut)
        err_Z_filter_i = Z_filter_i & (pa_diff_angle_cut)
        if Combined:
            good_Z_filter_i = (SOURCE_Z_GOOD < max_tol_Z) & (SOURCE_Z_GOOD > 0)
            PA_diffs_good_cases_Z_i_cut,PA_diffs_good_cases_Z_i_cut_err = PA_DIFF_GOOD[good_Z_filter_i],PA_DIFF_GOOD_ERR[good_Z_filter_i]
        else:
            good_Z_filter_i = Z_filter_i & (good_cases_cut)
            PA_diffs_good_cases_Z_i_cut,PA_diffs_good_cases_Z_i_cut_err = PA_DIFF[good_Z_filter_i],PA_DIFF_ERR[good_Z_filter_i]
            
        PA_diffs_Z_i_cut,PA_diffs_Z_i_cut_err = PA_DIFF[Z_filter_i],PA_DIFF_ERR[Z_filter_i]
        PA_diffs_angle_Z_i_cut,PA_diffs_angle_Z_i_cut_err = PA_DIFF[err_Z_filter_i],PA_DIFF_ERR[err_Z_filter_i]
        
        
        size_1 = len(PA_diffs_Z_i_cut)
        size_2 = len(PA_diffs_angle_Z_i_cut)
        size_3 = len(PA_diffs_good_cases_Z_i_cut)
        
        
        data_array[i+1,0,:size_1],data_array[i+1,1,:size_2],data_array[i+1,2,:size_3] = PA_diffs_Z_i_cut,PA_diffs_angle_Z_i_cut,PA_diffs_good_cases_Z_i_cut
        data_err_array[i+1,0,:size_1],data_err_array[i+1,1,:size_2],data_err_array[i+1,2,:size_3] = PA_diffs_Z_i_cut_err,PA_diffs_angle_Z_i_cut_err,PA_diffs_good_cases_Z_i_cut_err
        
    return data_array,data_err_array,n_rows


   
def Histograms_data(VLBI_PA_ERR,OPTICAL_PA_ERR,PA_DIFF,PA_DIFF_ERR,SOURCE_Z,max_tol_Z_array,max_tol_err,N_bins,Combined,good_cases_cut,PA_DIFF_GOOD,PA_DIFF_GOOD_ERR,SOURCE_Z_GOOD):

    
    data_array,data_err_array,n_rows = divide_data_into_panels(VLBI_PA_ERR,OPTICAL_PA_ERR,PA_DIFF,PA_DIFF_ERR,SOURCE_Z,max_tol_Z_array,max_tol_err,Combined,good_cases_cut,PA_DIFF_GOOD,PA_DIFF_GOOD_ERR,SOURCE_Z_GOOD)
    
    all_counts = np.zeros((n_rows,3,N_bins))
    
    N_draws = 10000
    N_bins_2 = 20
    
    for i in range(n_rows):
        for j in range(3):
            
            counts,bins = np.histogram(data_array[i,j],bins=np.linspace(0,90,N_bins+1))
            all_counts[i,j] = counts    
         
    return all_counts 



def Histograms_plots(data_array,data_err_array,max_tol_Z_array,Survey_Name,Png_Name,counts_all_all,counts_all_mean,counts_all_std,p_values,counts_all_all_5_bins,counts_all_mean_5_bins,counts_all_std_5_bins,p_values_5_bins,N_bins,N_draws,N_bins_2,N_bins_3,n_rows):
    
    fig,axes = plt.subplots(n_rows,3,sharex=True,figsize=(30,30))
    
    pad = 0
    cols = ['Any',r'$\sigma_{{PA}} < 22.5^{{\circ}}$','Good cases']
    rows = ['Any ','Only with z'] + ['Only with z < {}'.format(max_tol_Z_array[i]) for i in range(1,len(max_tol_Z_array))]
    
    
    for ax, col in zip(axes[0], cols):
        ax.annotate(col, xy=(0.5, 1.1), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                fontsize=30, ha='center', va='baseline')
    
    for ax, row in zip(axes[:,2], rows):
        ax.annotate(row, xy=(1.1, 0.5), xytext=(0,pad),
                xycoords='axes fraction', textcoords='offset points',
                fontsize=30, ha='right', va='center',rotation = 90)
        
        

    N_draws = 10000
    N_bins_2 = 20
    


    p_values_blue_5_bins = np.zeros((n_rows,3))
    p_values_blue_2_bins = np.zeros((n_rows,3))
    
    final_p_values_5_bins = np.zeros((n_rows,3))
    final_p_values_2_bins = np.zeros((n_rows,3))
    
    labels = np.array([['a)','b)','c)'],['d)','e)','f)'],['g)','h)','i)'],['j)','k)','l)']])
    label_seps = np.array([[-0.1,-0.1,-0.1],[-0.1,-0.1,-0.1],[-0.1,-0.1,-0.1],[-0.3,-0.3,-0.3]])
    percentage = 0
    for i in range(n_rows):
        for j in range(3):
         
    
            Original_counts_fig,Mean_counts_fig,Y_errs_fig,Bins_centers_fig,Bins_array_fig,Cov_matrix_blue = aux_functions.MonteCarlo(data_array[i,j],data_err_array[i,j],N_draws,N_bins,N_bins_2)
            bin_2_Original_counts_fig,bin_2_Mean_counts_fig,bin_2_Y_errs_fig,bin_2_Bins_centers_fig,bin_2_Bins_array_fig,bin_2_Cov_matrix_blue = aux_functions.MonteCarlo(data_array[i,j],data_err_array[i,j],N_draws,2,N_bins_2)
            
            try:
                popt_fig,pcov_fig = curve_fit(aux_functions.linear, Bins_centers_fig, Original_counts_fig, p0=None, sigma=np.mean(Y_errs_fig,axis=0), absolute_sigma=True, check_finite=True, bounds=(-np.inf, np.inf))
                perr_fig = np.sqrt(np.diag(pcov_fig))
                x_line_fig = np.linspace(0,90,100)
                Line_fig = x_line_fig*popt_fig[0]+popt_fig[1]
                significance_fig = popt_fig[0]/abs(perr_fig[0])
            except Exception:
                Line_fig = np.repeat(0,len(x_line_fig))

            blue_counts_2_bins,blue_bins_2_bins = np.histogram(data_array[i,j],bins=np.linspace(0,90,3))
            bin_2_sign_fig5 = aux_functions.bin_N_significance(data_array[i,j],data_err_array[i,j],2)
            
            p_values_red_2_bins = p_values[i,j,:]
            p_value_blue_2_bins = (1/2)*np.sum((bin_2_Original_counts_fig - counts_all_mean[i,j,:])**2/(counts_all_std[i,j,:]**2))

            final_p_value_2_bins = len(p_values_red_2_bins[p_values_red_2_bins > p_value_blue_2_bins])/len(p_values_red_2_bins)

            p_values_red_5_bins = p_values_5_bins[i,j,:]
            p_value_blue_5_bins = (1/N_bins)*np.sum((Original_counts_fig - counts_all_mean_5_bins[i,j,:])**2/(counts_all_std_5_bins[i,j,:]**2))

            final_p_value_5_bins = len(p_values_red_5_bins[p_values_red_5_bins > p_value_blue_5_bins])/len(p_values_red_5_bins)
            
            final_p_values_5_bins[i,j] = final_p_value_5_bins
            final_p_values_2_bins[i,j] = final_p_value_2_bins
            
            p_values_blue_5_bins[i,j] = p_value_blue_5_bins
            p_values_blue_2_bins[i,j] = p_value_blue_2_bins

            fill_x = np.linspace(0,90,N_bins+1)
            fill_y = np.array(list(counts_all_mean_5_bins[i,j])+[counts_all_mean_5_bins[i,j,-1]])
            fill_err = np.array(list(counts_all_std_5_bins[i,j])+[counts_all_std_5_bins[i,j,-1]])

            bin_2_fill_x = np.linspace(0,90,3)
            bin_2_fill_y = np.array(list(counts_all_mean[i,j])+[counts_all_mean[i,j,-1]])
            bin_2_fill_err = np.array(list(counts_all_std[i,j])+[counts_all_std[i,j,-1]])

        

           
            if final_p_value_5_bins < 0.0001:
                string_5_bins = r'$p_{{{:.0f}}} < 0.0001$'.format(N_bins)
               
            elif 0.01 > final_p_value_5_bins:
                string_5_bins = r'$p_{{{:.0f}}}$ = {:.3f}'.format(N_bins,final_p_value_5_bins)
            else:
                string_5_bins = r'$p_{{{:.0f}}}$ = {:.2f}'.format(N_bins,final_p_value_5_bins)

            
            if final_p_value_2_bins < 0.0001:
                string_2_bins = r'$p_2 < 0.0001$'
 
            elif 0.01 > final_p_value_2_bins:
                string_2_bins = r'$p_2$ = {:.3f}'.format(final_p_value_2_bins)
                
            else:
                string_2_bins = r'$p_2$ = {:.2f}'.format(final_p_value_2_bins)

            min_value = np.min(Original_counts_fig)
            min_pos = np.argmin(Original_counts_fig)
            min_value_error = Y_errs_fig[1][min_pos]
            max_value = np.max(Original_counts_fig)
            max_pos = np.argmax(Original_counts_fig)
            max_value_error = Y_errs_fig[0][max_pos]

            
            
            axes[i,j].text(0.1,0.05,string_5_bins,fontsize=30,transform=axes[i,j].transAxes,color='white')
            axes[i,j].text(0.5,0.05,string_2_bins,fontsize=30,transform=axes[i,j].transAxes,color='white')
            
            axes[i,j].hist(data_array[i,j],bins=Bins_array_fig,histtype = 'step', fill = True, color = 'dodgerblue', edgecolor='dodgerblue', linewidth=3, alpha = 0.7)
            axes[i,j].hist(data_array[i,j],bins=bin_2_Bins_array_fig,histtype = 'step', fill = True, color = 'lightcoral', edgecolor='lightcoral', linewidth=3, alpha=0.5)
            
            axes[i,j].errorbar(Bins_centers_fig,Original_counts_fig,yerr=Y_errs_fig,xerr=None,fmt='None',ecolor='k',elinewidth=2)
            axes[i,j].errorbar(bin_2_Bins_centers_fig,bin_2_Original_counts_fig,yerr=bin_2_Y_errs_fig,xerr=None,fmt='None',ecolor='k',elinewidth=2)
    
            axes[i,j].fill_between(fill_x,fill_y-fill_err,fill_y+fill_err, alpha=0.8,color='midnightblue',step='post')
            axes[i,j].fill_between(bin_2_fill_x,bin_2_fill_y-bin_2_fill_err,bin_2_fill_y+bin_2_fill_err, alpha=0.5,color='darkred',step='post')
            axes[i,j].text(0.5,label_seps[i,j], labels[i,j], size=25, ha="center", 
         transform=axes[i,j].transAxes)
            
            kstest_res_fig = kstest(data_array[i,j][np.isnan(data_array[i,j]) == False], 'uniform',args=(0,90))
            
            if kstest_res_fig.pvalue < 0.0001:
                kstest_label = r'$p_{{ks}}$ < 0.0001'
            elif kstest_res_fig.pvalue < 0.01:
                kstest_label = r'$p_{{ks}}$ = {:.3f}'.format(kstest_res_fig.pvalue)
            else: 
                kstest_label = r'$p_{{ks}}$ = {:.2f}'.format(kstest_res_fig.pvalue)
            
            axes[i,j].text(0.25,0.2,r'$p_{{ks}}$ = {:.3f}'.format(kstest_res_fig.pvalue),fontsize=30,transform=axes[i,j].transAxes,color='white')
           
        
            top_ylim = axes[i,j].get_ylim()[1]
            axes[i,j].set_ylim([0, top_ylim*1.1])
            axes[i,j].text(0.4,0.9,'N = {}'.format(len(data_array[i,j][np.isnan(data_array[i,j]) == False])),fontsize=30,transform=axes[i,j].transAxes,color='black')
                      
            
            if i == 3:
                string = r'$\Delta$PA ($^{{\circ}}$)'
                axes[i,j].set_xlabel(string,fontsize=30)
            axes[i,j].tick_params(axis='both', which='major', labelsize=30)
             
            percentage += 1./(n_rows*3)
            print('{:.0f}'.format(100*percentage)+'%', end='\r') 
    p_values_blue_2_bins = np.zeros((4,3))

    fig.suptitle(Survey_Name,fontsize=50)
    plt.subplots_adjust(top=0.9)
    plt.savefig('Paper_images/'+Png_Name,dpi=100,bbox_inches = 'tight',
    pad_inches = 0)

    plt.show()
    
    return final_p_values_5_bins,final_p_values_2_bins




def Histograms(VLBI_PA,VLBI_PA_ERR,OPTICAL_PA,OPTICAL_PA_ERR,PA_DIFF,PA_DIFF_ERR,SOURCE_Z,max_tol_err,max_tol_Z_array,Survey_Name,Png_Name,N_bins,Combined,good_cases_cut=None,VLBI_PA_GOOD=np.zeros(5),VLBI_PA_GOOD_ERR=np.zeros(5),OPTICAL_PA_GOOD=np.zeros(5),OPTICAL_PA_GOOD_ERR=np.zeros(5),PA_DIFF_GOOD=np.zeros(5),PA_DIFF_GOOD_ERR=np.zeros(5),SOURCE_Z_GOOD=np.zeros(5)):
    
    N_bins_shuffled = 2
    N_shuffles = 1000 
    
    if good_cases_cut is None:
        good_cases_cut = np.zeros(len(PA_DIFF))
    
    VLBI_indexes = np.arange(0,len(VLBI_PA))
    VLBI_indexes_copy = VLBI_indexes.copy()
    OPTICAL_indexes = np.arange(0,len(OPTICAL_PA))
    OPTICAL_indexes_copy = OPTICAL_indexes.copy()
    
    VLBI_GOOD_indexes = np.arange(0,len(VLBI_PA_GOOD))
    VLBI_GOOD_indexes_copy = VLBI_GOOD_indexes.copy()
    OPTICAL_GOOD_indexes = np.arange(0,len(OPTICAL_PA_GOOD))
    OPTICAL_GOOD_indexes_copy = OPTICAL_GOOD_indexes.copy()
    
    
    n_rows = len(max_tol_Z_array) + 1
    
    counts_all_all_2_bins = np.zeros((n_rows,3,2,N_shuffles))
    counts_all_all_5_bins = np.zeros((n_rows,3,N_bins,N_shuffles))
   
    p_values_all_2_bins = np.zeros((n_rows,3,N_shuffles))
    p_values_all_5_bins = np.zeros((n_rows,3,N_shuffles))
    
    
    bin_2_data_significance,bin_2_data_errors = aux_functions.bin_N_significance(PA_DIFF,PA_DIFF_ERR,2)
    bin_5_data_significance,bin_5_data_errors = aux_functions.bin_N_significance(PA_DIFF,PA_DIFF_ERR,N_bins)
    
    print('Computing angle shuffling for the null hypotheses...')
    
    for i in range(N_shuffles):
        print('{:.2f}%'.format(100*i/N_shuffles), end='\r')
        random.shuffle(VLBI_indexes_copy)
        random.shuffle(OPTICAL_indexes_copy)
        random.shuffle(VLBI_GOOD_indexes_copy)
        random.shuffle(OPTICAL_GOOD_indexes_copy)
        
        VLBI_PA_shuffled,VLBI_PA_ERR_shuffled = VLBI_PA[VLBI_indexes_copy],VLBI_PA_ERR[VLBI_indexes_copy]
        OPTICAL_PA_shuffled,OPTICAL_PA_ERR_shuffled,SOURCE_Z_shuffled = OPTICAL_PA[OPTICAL_indexes_copy],OPTICAL_PA_ERR[OPTICAL_indexes_copy],SOURCE_Z[OPTICAL_indexes_copy]
        PA_DIFF_shuffled = aux_functions.PA_difference(OPTICAL_PA_shuffled,VLBI_PA_shuffled)
        PA_DIFF_ERR_shuffled = np.sqrt(VLBI_PA_ERR_shuffled**2+OPTICAL_PA_ERR_shuffled**2)
        
        OPTICAL_PA_GOOD_shuffled,OPTICAL_PA_GOOD_ERR_shuffled = OPTICAL_PA_GOOD[OPTICAL_GOOD_indexes_copy],OPTICAL_PA_GOOD_ERR[OPTICAL_GOOD_indexes_copy]
        VLBI_PA_GOOD_shuffled,VLBI_PA_GOOD_ERR_shuffled = VLBI_PA_GOOD[VLBI_GOOD_indexes_copy],VLBI_PA_GOOD_ERR[VLBI_GOOD_indexes_copy]
            
        SOURCE_Z_GOOD_shuffled = SOURCE_Z_GOOD[OPTICAL_GOOD_indexes_copy]
        PA_DIFF_GOOD_shuffled = aux_functions.PA_difference(OPTICAL_PA_GOOD_shuffled,VLBI_PA_GOOD_shuffled) 
        PA_DIFF_GOOD_ERR_shuffled = np.sqrt(VLBI_PA_GOOD_ERR_shuffled**2+OPTICAL_PA_GOOD_ERR_shuffled**2)
        
        
        good_cases_cut_shuffled = good_cases_cut[OPTICAL_indexes_copy]
         
      
        counts_all_2_bins = Histograms_data(VLBI_PA_ERR_shuffled,OPTICAL_PA_ERR_shuffled,PA_DIFF_shuffled,PA_DIFF_ERR_shuffled,SOURCE_Z_shuffled,max_tol_Z_array,max_tol_err,2,Combined,good_cases_cut_shuffled,PA_DIFF_GOOD_shuffled,PA_DIFF_GOOD_ERR_shuffled,SOURCE_Z_GOOD_shuffled)
        counts_all_5_bins = Histograms_data(VLBI_PA_ERR_shuffled,OPTICAL_PA_ERR_shuffled,PA_DIFF_shuffled,PA_DIFF_ERR_shuffled,SOURCE_Z_shuffled,max_tol_Z_array,max_tol_err,N_bins,Combined,good_cases_cut_shuffled,PA_DIFF_GOOD_shuffled,PA_DIFF_GOOD_ERR_shuffled,SOURCE_Z_GOOD_shuffled)
        
        
   
        counts_all_all_2_bins[:,:,:,i] = counts_all_2_bins
        counts_all_all_5_bins[:,:,:,i] = counts_all_5_bins
  
    counts_all_mean_2_bins = np.mean(counts_all_all_2_bins,axis=3)
    counts_all_std_2_bins = np.std(counts_all_all_2_bins,axis=3)
    counts_all_mean_5_bins = np.mean(counts_all_all_5_bins,axis=3)
    counts_all_std_5_bins = np.std(counts_all_all_5_bins,axis=3)
  
    for j in range(N_shuffles):
               
        p_values_all_2_bins[:,:,j] = (1./2)*np.sum(np.abs(counts_all_all_2_bins[:,:,:,j]-counts_all_mean_2_bins)**2/(counts_all_std_2_bins**2),axis=2)

        p_values_all_5_bins[:,:,j] = (1./N_bins)*np.sum(np.abs(counts_all_all_5_bins[:,:,:,j]-counts_all_mean_5_bins)**2/(counts_all_std_5_bins**2),axis=2)

    
    data_array,data_err_array,n_rows = divide_data_into_panels(VLBI_PA_ERR,OPTICAL_PA_ERR,PA_DIFF,PA_DIFF_ERR,SOURCE_Z,max_tol_Z_array,max_tol_err,Combined,good_cases_cut,PA_DIFF_GOOD,PA_DIFF_GOOD_ERR,SOURCE_Z_GOOD)
        


        

    N_draws = 10000
    N_bins_2 = 20
    N_bins_3 = 100
   
    print('Plotting histogram...')
    final_p_values_5_bins,final_p_values_2_bins = Histograms_plots(data_array,data_err_array,max_tol_Z_array,Survey_Name,Png_Name,counts_all_all_2_bins,counts_all_mean_2_bins,counts_all_std_2_bins,p_values_all_2_bins,counts_all_all_5_bins,counts_all_mean_5_bins,counts_all_std_5_bins,p_values_all_5_bins,N_bins,N_draws,N_bins_2,N_bins_3,n_rows)
  
    
    return final_p_values_5_bins,final_p_values_2_bins