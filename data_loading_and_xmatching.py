from libraries import *
import aux_functions

def load_Astrogeo_catalogue():
    print('Loading Astrogeo...')
    Astrogeo_catalogue = pd.read_csv('Source_coords_csv.csv')
    return Astrogeo_catalogue

def get_Astrogeo_good_catalogue(Astrogeo_catalogue):
    Astrogeo_RA = np.array(Astrogeo_catalogue.RA_deg)
    Astrogeo_DEC = np.array(Astrogeo_catalogue.DEC_deg)
    RA_bad = (((Astrogeo_RA > 0) & (Astrogeo_RA < 360)) == False) & (np.isnan(Astrogeo_RA) == True)
    DEC_bad = (((Astrogeo_DEC > (-90)) & (Astrogeo_DEC < 90)) == False) & (np.isnan(Astrogeo_DEC) == True)
    RA_good = (((Astrogeo_RA > 0) & (Astrogeo_RA < 360)) == True) & (np.isnan(Astrogeo_RA) == False)
    DEC_good = (((Astrogeo_DEC > (-90)) & (Astrogeo_DEC < 90)) == True) & (np.isnan(Astrogeo_DEC) == False)
    Astrogeo_good = Astrogeo_catalogue[RA_good & DEC_good]
    return Astrogeo_good
    

def load_catalogues():
    print('Loading DESI...')
    DESI_catalogue = pd.read_csv('VLBI_source_matches_dr9-10_combi_complete.csv',header=0)
    print('Loading SDSS...')
    SDSS_catalogue = pd.read_csv('SDSS_DR17_Gal_QSO_jasorey_complete.csv',header=0)
    print('Loading DES...')
    DES_catalogue = pd.read_csv('DES_full_catalogue_Astrogeo.csv',header=0)
    print('Loading Skymapper...')
    Skymapper_catalogue = pd.read_csv('skymapper.csv',header=0)
    print('Loading KIDS...')
    with fits.open('KiDS_DR4.1_ugriZYJHKs_SOM_gold_WL_cat.fits',memmap=True) as KIDS_fits:
        KIDS_catalogue = KIDS_fits[1].data
    return DESI_catalogue,SDSS_catalogue,DES_catalogue,Skymapper_catalogue,KIDS_catalogue

def load_sampled_catalogues(size,SDSS_catalogue,KIDS_catalogue):
    print('Loading sample of SDSS...')
    SDSS_index = np.arange(0,len(SDSS_catalogue))
    SDSS_random_index = np.array(random.sample(list(SDSS_index), k=size))
    SDSS_sampled_catalogue = SDSS_catalogue.iloc[SDSS_random_index]
    print('Loading sample of KIDS...')
    KIDS_index = np.arange(0,len(KIDS_catalogue))
    KIDS_random_index = np.array(random.sample(list(KIDS_index), k=size))
    KIDS_sampled_catalogue = KIDS_catalogue[KIDS_random_index]
    return SDSS_sampled_catalogue,KIDS_sampled_catalogue


def perform_xmatch(Astrogeo_good,catalogue,catalogue_name,catalogue_RA_name,catalogue_DEC_name):
    print('Performing xmatch on '+catalogue_name+'...')
    Astrogeo_RA_good = np.array(Astrogeo_good.RA_deg)
    Astrogeo_DEC_good = np.array(Astrogeo_good.DEC_deg)
    Astrogeo_coords_good = SkyCoord(Astrogeo_RA_good,Astrogeo_DEC_good,frame='icrs',unit='deg')
    
    catalogue_RA = catalogue[catalogue_RA_name]
    catalogue_DEC = catalogue[catalogue_DEC_name]
    catalogue_coords = SkyCoord(catalogue_RA,catalogue_DEC,frame='icrs',unit='deg')
    
    max_sep = 1 * u.arcsec
    idx, d2d, d3d = Astrogeo_coords_good.match_to_catalog_sky(catalogue_coords)
    sep_constraint = d2d < max_sep
    Astrogeo_catalogue_matches = Astrogeo_good[sep_constraint]
    if catalogue_name == 'KIDS':
        catalogue_Astrogeo_matches = catalogue[idx[sep_constraint]]
    else:
        catalogue_Astrogeo_matches = catalogue.iloc[idx[sep_constraint]]
    
    return Astrogeo_catalogue_matches,catalogue_Astrogeo_matches

def load_Eagle_sim():
    print('Loading Eagle simulation...')
    with fits.open('eagle_new.fits') as Eagle_fits:
        Eagle_sim = Eagle_fits[1].data
    Eagle_sim = Table(Eagle_sim)
    return Eagle_sim



def get_DESI_data(Astrogeo_DESI,DESI_xmatches):
    DESI_catalogue_Astrogeo_B_matches = np.array(DESI_xmatches.b_axis)
    DESI_catalogue_Astrogeo_THETA_J2000_matches = np.array(DESI_xmatches.pos_angle)
    DESI_catalogue_Astrogeo_THETA_J2000_ERR_ORIGINAL_matches = np.array(DESI_xmatches.pos_angle_err)
    DESI_catalogue_Astrogeo_TYPE_matches = np.array(DESI_xmatches.TYPE)
    DESI_catalogue_Astrogeo_MAG_Z_matches = np.array(DESI_xmatches.mag_z)
    Astrogeo_DESI_catalogue_JET_PA_matches_good = np.array(Astrogeo_DESI.pa)
    Astrogeo_DESI_catalogue_JET_PA_ERR_ORIGINAL_matches_good = np.array(Astrogeo_DESI.pa_err)
    Astrogeo_DESI_catalogue_Z_good_matches = np.array(Astrogeo_DESI.Z)
    
    DESI_catalogue_Astrogeo_THETA_J2000_matches = np.fmod(DESI_catalogue_Astrogeo_THETA_J2000_matches+360,180)
    DESI_catalogue_Astrogeo_THETA_J2000_matches[DESI_catalogue_Astrogeo_THETA_J2000_matches > 90] = DESI_catalogue_Astrogeo_THETA_J2000_matches[DESI_catalogue_Astrogeo_THETA_J2000_matches > 90] - 180      
    DESI_catalogue_Astrogeo_PA_DIFFERENCE_matches = aux_functions.PA_difference(DESI_catalogue_Astrogeo_THETA_J2000_matches,Astrogeo_DESI_catalogue_JET_PA_matches_good)   
    DESI_catalogue_Astrogeo_PA_DIFFERENCE_ERR_ORIGINAL_matches = np.sqrt(Astrogeo_DESI_catalogue_JET_PA_ERR_ORIGINAL_matches_good**2+DESI_catalogue_Astrogeo_THETA_J2000_ERR_ORIGINAL_matches**2)
   
    DESI_SOURCE_Z = np.array(Astrogeo_DESI_catalogue_Z_good_matches)
    DESI_VLBI_PA = np.array(Astrogeo_DESI_catalogue_JET_PA_matches_good)
    DESI_VLBI_PA_ERR_ORIGINAL = np.array(Astrogeo_DESI_catalogue_JET_PA_ERR_ORIGINAL_matches_good)
    DESI_OPTICAL_PA = np.array(DESI_catalogue_Astrogeo_THETA_J2000_matches)
    DESI_OPTICAL_PA_ERR_ORIGINAL = np.array(DESI_catalogue_Astrogeo_THETA_J2000_ERR_ORIGINAL_matches)
    DESI_PA_DIFF = np.array(DESI_catalogue_Astrogeo_PA_DIFFERENCE_matches)
    DESI_PA_DIFF_ERR_ORIGINAL = np.array(DESI_catalogue_Astrogeo_PA_DIFFERENCE_ERR_ORIGINAL_matches)
    
    return DESI_OPTICAL_PA,DESI_OPTICAL_PA_ERR_ORIGINAL,DESI_VLBI_PA,DESI_VLBI_PA_ERR_ORIGINAL,DESI_PA_DIFF,DESI_PA_DIFF_ERR_ORIGINAL,DESI_SOURCE_Z,DESI_catalogue_Astrogeo_B_matches,DESI_catalogue_Astrogeo_TYPE_matches,DESI_catalogue_Astrogeo_MAG_Z_matches



def get_SDSS_data(Astrogeo_SDSS,SDSS_xmatches):
  
    SDSS_catalogue_Astrogeo_THETA_J2000_matches = np.array(SDSS_xmatches.modelPhi_r)
    SDSS_catalogue_Astrogeo_THETA_J2000_ERR_ORIGINAL_matches = np.array(SDSS_xmatches.modelPhiErr)
    Astrogeo_SDSS_catalogue_JET_PA_matches_good = np.array(Astrogeo_SDSS.pa)
    Astrogeo_SDSS_catalogue_JET_PA_ERR_ORIGINAL_matches_good = np.array(Astrogeo_SDSS.pa_err)
    Astrogeo_SDSS_catalogue_Z_matches_good = np.array(Astrogeo_SDSS.Z)
    SDSS_catalogue_Astrogeo_TYPE_matches = np.array(SDSS_xmatches.type_r)
    SDSS_catalogue_Astrogeo_B_matches = np.array(SDSS_xmatches.modelAB_r)*np.array(SDSS_xmatches.devrad_r)
    
    SDSS_catalogue_Astrogeo_THETA_J2000_matches = np.fmod(SDSS_catalogue_Astrogeo_THETA_J2000_matches+360,180)
    SDSS_catalogue_Astrogeo_THETA_J2000_matches[SDSS_catalogue_Astrogeo_THETA_J2000_matches > 90] = SDSS_catalogue_Astrogeo_THETA_J2000_matches[SDSS_catalogue_Astrogeo_THETA_J2000_matches > 90] - 180    
    SDSS_catalogue_Astrogeo_PA_DIFFERENCE_matches = aux_functions.PA_difference(SDSS_catalogue_Astrogeo_THETA_J2000_matches,Astrogeo_SDSS_catalogue_JET_PA_matches_good)
    SDSS_catalogue_Astrogeo_PA_DIFFERENCE_ERR_ORIGINAL_matches = np.sqrt(Astrogeo_SDSS_catalogue_JET_PA_ERR_ORIGINAL_matches_good**2+SDSS_catalogue_Astrogeo_THETA_J2000_ERR_ORIGINAL_matches**2)    
    
    SDSS_SOURCE_Z = np.array(Astrogeo_SDSS_catalogue_Z_matches_good)
    SDSS_VLBI_PA = np.array(Astrogeo_SDSS_catalogue_JET_PA_matches_good)
    SDSS_VLBI_PA_ERR_ORIGINAL = np.array(Astrogeo_SDSS_catalogue_JET_PA_ERR_ORIGINAL_matches_good)
    SDSS_OPTICAL_PA = np.array(SDSS_catalogue_Astrogeo_THETA_J2000_matches)
    SDSS_OPTICAL_PA_ERR_ORIGINAL = np.array(SDSS_catalogue_Astrogeo_THETA_J2000_ERR_ORIGINAL_matches)
    SDSS_PA_DIFF = np.array(SDSS_catalogue_Astrogeo_PA_DIFFERENCE_matches)
    SDSS_PA_DIFF_ERR_ORIGINAL = np.array(SDSS_catalogue_Astrogeo_PA_DIFFERENCE_ERR_ORIGINAL_matches)
   
    return SDSS_OPTICAL_PA,SDSS_OPTICAL_PA_ERR_ORIGINAL,SDSS_VLBI_PA,SDSS_VLBI_PA_ERR_ORIGINAL, SDSS_PA_DIFF,SDSS_PA_DIFF_ERR_ORIGINAL,SDSS_SOURCE_Z,SDSS_catalogue_Astrogeo_TYPE_matches,SDSS_catalogue_Astrogeo_B_matches


def get_DES_data(Astrogeo_DES,DES_xmatches):

    DES_OPTICAL_PA = np.array(DES_xmatches.DES_full_catalogue_Astrogeo_THETA_J2000_matches)
    DES_OPTICAL_PA_ERR_ORIGINAL = np.array(DES_xmatches.DES_full_catalogue_Astrogeo_ERRTHETA_IMAGE_matches)
    DES_VLBI_PA = np.array(Astrogeo_DES.pa)
    DES_VLBI_PA_ERR_ORIGINAL = np.array(Astrogeo_DES.pa_err)
    DES_PA_DIFF = np.array(DES_xmatches.DES_full_catalogue_Astrogeo_PA_DIFFERENCE_matches)
    DES_PA_DIFF_ERR_ORIGINAL = np.sqrt(DES_VLBI_PA_ERR_ORIGINAL**2+DES_OPTICAL_PA_ERR_ORIGINAL**2)  
    DES_SOURCE_Z = np.array(DES_xmatches.Astrogeo_DES_full_catalogue_Z_good_matches)
    DES_catalogue_Astrogeo_B_matches = 0.26*np.array(DES_xmatches.DES_full_catalogue_Astrogeo_B_matches)
    DES_catalogue_Astrogeo_EXTENDED_CLASS_COADD_matches = np.array(DES_xmatches.DES_full_catalogue_Astrogeo_EXTENDED_CLASS_COADD_matches)
    
    return DES_OPTICAL_PA,DES_OPTICAL_PA_ERR_ORIGINAL,DES_VLBI_PA,DES_VLBI_PA_ERR_ORIGINAL,DES_PA_DIFF,DES_PA_DIFF_ERR_ORIGINAL,DES_SOURCE_Z,DES_catalogue_Astrogeo_EXTENDED_CLASS_COADD_matches,DES_catalogue_Astrogeo_B_matches


def get_Skymapper_data(Astrogeo_Skymapper,Skymapper_xmatches):
    
    skymapper_Astrogeo_B_matches = 0.5*np.array(Skymapper_xmatches.b)
    skymapper_Astrogeo_THETA_J2000_matches = np.array(Skymapper_xmatches.PA)
    skymapper_Astrogeo_THETA_J2000_ERR_ORIGINAL_matches = np.array(Skymapper_xmatches.e_pa)
    Astrogeo_skymapper_JET_PA_good_matches = np.array(Astrogeo_Skymapper.pa)
    Astrogeo_skymapper_JET_PA_ERR_ORIGINAL_good_matches = np.array(Astrogeo_Skymapper.pa_err)
    skymapper_Astrogeo_TYPE_matches = np.array(Skymapper_xmatches.ClassStar)
    Astrogeo_skymapper_Z_good_matches = np.array(Astrogeo_Skymapper.Z)

    skymapper_Astrogeo_THETA_J2000_matches = np.fmod(skymapper_Astrogeo_THETA_J2000_matches+360,180)
    skymapper_Astrogeo_THETA_J2000_matches[skymapper_Astrogeo_THETA_J2000_matches > 90] = skymapper_Astrogeo_THETA_J2000_matches[skymapper_Astrogeo_THETA_J2000_matches > 90] - 180 
    skymapper_Astrogeo_PA_DIFFERENCE_matches = aux_functions.PA_difference(skymapper_Astrogeo_THETA_J2000_matches,Astrogeo_skymapper_JET_PA_good_matches)
    skymapper_Astrogeo_PA_DIFFERENCE_ERR_ORIGINAL_matches = np.sqrt(Astrogeo_skymapper_JET_PA_ERR_ORIGINAL_good_matches**2+skymapper_Astrogeo_THETA_J2000_ERR_ORIGINAL_matches**2)
        
    skymapper_SOURCE_Z = np.array(Astrogeo_skymapper_Z_good_matches)
    skymapper_VLBI_PA = np.array(Astrogeo_skymapper_JET_PA_good_matches)
    skymapper_VLBI_PA_ERR_ORIGINAL = np.array(Astrogeo_skymapper_JET_PA_ERR_ORIGINAL_good_matches)
    skymapper_OPTICAL_PA = np.array(skymapper_Astrogeo_THETA_J2000_matches)
    skymapper_OPTICAL_PA_ERR_ORIGINAL = np.array(skymapper_Astrogeo_THETA_J2000_ERR_ORIGINAL_matches)
    skymapper_PA_DIFF = np.array(skymapper_Astrogeo_PA_DIFFERENCE_matches)
    skymapper_PA_DIFF_ERR_ORIGINAL = np.array(skymapper_Astrogeo_PA_DIFFERENCE_ERR_ORIGINAL_matches)
   
    
    return skymapper_OPTICAL_PA,skymapper_OPTICAL_PA_ERR_ORIGINAL,skymapper_VLBI_PA,skymapper_VLBI_PA_ERR_ORIGINAL,skymapper_PA_DIFF,skymapper_PA_DIFF_ERR_ORIGINAL,skymapper_SOURCE_Z,skymapper_Astrogeo_TYPE_matches,skymapper_Astrogeo_B_matches



def get_Combined_data(Surveys_Data,Surveys_Good_Data):
    
    
    Combined_VLBI_PA,Combined_VLBI_PA_ERR_ORIGINAL,Combined_Good_VLBI_PA,Combined_Good_VLBI_PA_ERR_ORIGINAL,Combined_OPTICAL_PA,Combined_OPTICAL_PA_ERR_ORIGINAL,Combined_Good_OPTICAL_PA,Combined_Good_OPTICAL_PA_ERR_ORIGINAL,Combined_PA_DIFF,Combined_PA_DIFF_ERR_ORIGINAL,Combined_Good_PA_DIFF,Combined_Good_PA_DIFF_ERR_ORIGINAL,Combined_SOURCE_Z,Combined_Good_SOURCE_Z = np.array(Surveys_Data.VLBI_Jet_PAs)[pd.isnull(Surveys_Data.Weighted_average_Jet_Minor_diffs) == False],np.array(Surveys_Data.VLBI_Jet_PA_errors)[pd.isnull(Surveys_Data.Weighted_average_Jet_Minor_diffs) == False],np.array(Surveys_Good_Data.VLBI_Jet_PAs)[pd.isnull(Surveys_Good_Data.Weighted_average_Jet_Minor_diffs) == False],np.array(Surveys_Good_Data.VLBI_Jet_PA_errors)[pd.isnull(Surveys_Good_Data.Weighted_average_Jet_Minor_diffs) == False],np.array(Surveys_Data.Weighted_average_Optical_Major_PAs)[pd.isnull(Surveys_Data.Weighted_average_Jet_Minor_diffs) == False],np.array(Surveys_Data.Weighted_average_Optical_Major_PA_errors)[pd.isnull(Surveys_Data.Weighted_average_Jet_Minor_diffs) == False],np.array(Surveys_Good_Data.Weighted_average_Optical_Major_PAs)[pd.isnull(Surveys_Good_Data.Weighted_average_Jet_Minor_diffs) == False],np.array(Surveys_Good_Data.Weighted_average_Optical_Major_PA_errors)[pd.isnull(Surveys_Good_Data.Weighted_average_Jet_Minor_diffs) == False],np.array(Surveys_Data.Weighted_average_Jet_Minor_diffs)[pd.isnull(Surveys_Data.Weighted_average_Jet_Minor_diffs) == False],np.array(Surveys_Data.Weighted_average_Jet_Minor_diffs_errors)[pd.isnull(Surveys_Data.Weighted_average_Jet_Minor_diffs) == False],np.array(Surveys_Good_Data.Weighted_average_Jet_Minor_diffs)[pd.isnull(Surveys_Good_Data.Weighted_average_Jet_Minor_diffs) == False],np.array(Surveys_Good_Data.Weighted_average_Jet_Minor_diffs_errors)[pd.isnull(Surveys_Good_Data.Weighted_average_Jet_Minor_diffs) == False],np.array(Surveys_Data.Z)[pd.isnull(Surveys_Data.Weighted_average_Jet_Minor_diffs) == False],np.array(Surveys_Good_Data.Z)[pd.isnull(Surveys_Good_Data.Weighted_average_Jet_Minor_diffs) == False]
    return Combined_VLBI_PA,Combined_VLBI_PA_ERR_ORIGINAL,Combined_Good_VLBI_PA,Combined_Good_VLBI_PA_ERR_ORIGINAL,Combined_OPTICAL_PA,Combined_OPTICAL_PA_ERR_ORIGINAL,Combined_Good_OPTICAL_PA,Combined_Good_OPTICAL_PA_ERR_ORIGINAL,Combined_PA_DIFF,Combined_PA_DIFF_ERR_ORIGINAL,Combined_Good_PA_DIFF,Combined_Good_PA_DIFF_ERR_ORIGINAL,Combined_SOURCE_Z,Combined_Good_SOURCE_Z 