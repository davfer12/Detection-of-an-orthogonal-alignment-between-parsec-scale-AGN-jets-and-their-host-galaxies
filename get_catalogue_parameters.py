from libraries import *
import aux_functions

def get_DESI_data(Astrogeo_DESI, DESI_xmatches):
    """
    Processes and extracts various data arrays from the DESI cross-matched results.

    This function extracts and processes relevant data from the `DESI_xmatches` object 
    and the `Astrogeo_DESI` DataFrame. It performs calculations and returns arrays for 
    different parameters including positional angles, differences, and errors.

    Args:
        Astrogeo_DESI (pd.DataFrame): The filtered Astrogeo DESI catalogue containing data on positional angles and redshifts.
        DESI_xmatches (pd.DataFrame): The DESI cross-matched results containing various parameters.

    Returns:
        tuple: A tuple containing:
            - np.ndarray: Optical positional angles from DESI catalogue.
            - np.ndarray: Errors on optical positional angles.
            - np.ndarray: VLBI positional angles from Astrogeo catalogue.
            - np.ndarray: Errors on VLBI positional angles.
            - np.ndarray: Difference between optical and VLBI positional angles.
            - np.ndarray: Errors on the positional angle differences.
            - np.ndarray: Redshifts from the Astrogeo DESI catalogue.
            - np.ndarray: B-axis values from DESI cross-matches.
            - np.ndarray: Types from DESI cross-matches.
            - np.ndarray: Magnitudes (z-band) from DESI cross-matches.
    """
    
    # Extract arrays from DESI cross-matched results
    DESI_catalogue_Astrogeo_B_matches = np.array(DESI_xmatches.b_axis)
    DESI_catalogue_Astrogeo_THETA_J2000_matches = np.array(DESI_xmatches.pos_angle)
    DESI_catalogue_Astrogeo_THETA_J2000_ERR_ORIGINAL_matches = np.array(DESI_xmatches.pos_angle_err)
    DESI_catalogue_Astrogeo_TYPE_matches = np.array(DESI_xmatches.TYPE)
    DESI_catalogue_Astrogeo_MAG_Z_matches = np.array(DESI_xmatches.mag_z)
    
    # Extract arrays from Astrogeo DESI catalogue
    Astrogeo_DESI_catalogue_JET_PA_matches_good = np.array(Astrogeo_DESI.pa)
    Astrogeo_DESI_catalogue_JET_PA_ERR_ORIGINAL_matches_good = np.array(Astrogeo_DESI.pa_err)
    Astrogeo_DESI_catalogue_Z_good_matches = np.array(Astrogeo_DESI.Z)
    
    # Adjust the positional angles to the range [0, 180) degrees
    DESI_catalogue_Astrogeo_THETA_J2000_matches = np.fmod(DESI_catalogue_Astrogeo_THETA_J2000_matches + 360, 180)
    DESI_catalogue_Astrogeo_THETA_J2000_matches[DESI_catalogue_Astrogeo_THETA_J2000_matches > 90] -= 180
    
    # Compute the positional angle differences and errors
    DESI_catalogue_Astrogeo_PA_DIFFERENCE_matches = aux_functions.PA_difference(
        DESI_catalogue_Astrogeo_THETA_J2000_matches, Astrogeo_DESI_catalogue_JET_PA_matches_good
    )
    DESI_catalogue_Astrogeo_PA_DIFFERENCE_ERR_ORIGINAL_matches = np.sqrt(
        Astrogeo_DESI_catalogue_JET_PA_ERR_ORIGINAL_matches_good**2 +
        DESI_catalogue_Astrogeo_THETA_J2000_ERR_ORIGINAL_matches**2
    )
    
    # Prepare arrays for output
    DESI_SOURCE_Z = np.array(Astrogeo_DESI_catalogue_Z_good_matches)
    DESI_VLBI_PA = np.array(Astrogeo_DESI_catalogue_JET_PA_matches_good)
    DESI_VLBI_PA_ERR_ORIGINAL = np.array(Astrogeo_DESI_catalogue_JET_PA_ERR_ORIGINAL_matches_good)
    DESI_OPTICAL_PA = np.array(DESI_catalogue_Astrogeo_THETA_J2000_matches)
    DESI_OPTICAL_PA_ERR_ORIGINAL = np.array(DESI_catalogue_Astrogeo_THETA_J2000_ERR_ORIGINAL_matches)
    DESI_PA_DIFF = np.array(DESI_catalogue_Astrogeo_PA_DIFFERENCE_matches)
    DESI_PA_DIFF_ERR_ORIGINAL = np.array(DESI_catalogue_Astrogeo_PA_DIFFERENCE_ERR_ORIGINAL_matches)
    
    # Return the processed data as a tuple
    return (
        DESI_OPTICAL_PA, DESI_OPTICAL_PA_ERR_ORIGINAL,
        DESI_VLBI_PA, DESI_VLBI_PA_ERR_ORIGINAL,
        DESI_PA_DIFF, DESI_PA_DIFF_ERR_ORIGINAL,
        DESI_SOURCE_Z, DESI_catalogue_Astrogeo_B_matches,
        DESI_catalogue_Astrogeo_TYPE_matches, DESI_catalogue_Astrogeo_MAG_Z_matches
    )




def get_SDSS_data(Astrogeo_SDSS, SDSS_xmatches):
    """
    Processes and extracts various data arrays from the SDSS cross-matched results.

    This function extracts and processes relevant data from the `SDSS_xmatches` object 
    and the `Astrogeo_SDSS` DataFrame. It performs calculations to obtain arrays for 
    positional angles, differences, errors, and other parameters.

    Args:
        Astrogeo_SDSS (pd.DataFrame): The filtered Astrogeo SDSS catalogue containing data on positional angles and redshifts.
        SDSS_xmatches (pd.DataFrame): The SDSS cross-matched results containing various parameters.

    Returns:
        tuple: A tuple containing:
            - np.ndarray: Optical positional angles from the SDSS catalogue.
            - np.ndarray: Errors on optical positional angles.
            - np.ndarray: VLBI positional angles from the Astrogeo catalogue.
            - np.ndarray: Errors on VLBI positional angles.
            - np.ndarray: Difference between optical and VLBI positional angles.
            - np.ndarray: Errors on the positional angle differences.
            - np.ndarray: Redshifts from the Astrogeo SDSS catalogue.
            - np.ndarray: Types from SDSS cross-matches.
            - np.ndarray: B-axis values from SDSS cross-matches, scaled by deviance radius.
    """
    
    # Extract arrays from SDSS cross-matched results
    SDSS_catalogue_Astrogeo_THETA_J2000_matches = np.array(SDSS_xmatches.modelPhi_r)
    SDSS_catalogue_Astrogeo_THETA_J2000_ERR_ORIGINAL_matches = np.array(SDSS_xmatches.modelPhiErr)
    SDSS_catalogue_Astrogeo_TYPE_matches = np.array(SDSS_xmatches.type_r)
    SDSS_catalogue_Astrogeo_B_matches = np.array(SDSS_xmatches.modelAB_r) * np.array(SDSS_xmatches.devrad_r)
    
    # Extract arrays from Astrogeo SDSS catalogue
    Astrogeo_SDSS_catalogue_JET_PA_matches_good = np.array(Astrogeo_SDSS.pa)
    Astrogeo_SDSS_catalogue_JET_PA_ERR_ORIGINAL_matches_good = np.array(Astrogeo_SDSS.pa_err)
    Astrogeo_SDSS_catalogue_Z_matches_good = np.array(Astrogeo_SDSS.Z)
    
    # Adjust the positional angles to the range [0, 180) degrees
    SDSS_catalogue_Astrogeo_THETA_J2000_matches = np.fmod(SDSS_catalogue_Astrogeo_THETA_J2000_matches + 360, 180)
    SDSS_catalogue_Astrogeo_THETA_J2000_matches[SDSS_catalogue_Astrogeo_THETA_J2000_matches > 90] -= 180
    
    # Compute the positional angle differences and errors
    SDSS_catalogue_Astrogeo_PA_DIFFERENCE_matches = aux_functions.PA_difference(
        SDSS_catalogue_Astrogeo_THETA_J2000_matches, Astrogeo_SDSS_catalogue_JET_PA_matches_good
    )
    SDSS_catalogue_Astrogeo_PA_DIFFERENCE_ERR_ORIGINAL_matches = np.sqrt(
        Astrogeo_SDSS_catalogue_JET_PA_ERR_ORIGINAL_matches_good**2 +
        SDSS_catalogue_Astrogeo_THETA_J2000_ERR_ORIGINAL_matches**2
    )
    
    # Prepare arrays for output
    SDSS_SOURCE_Z = np.array(Astrogeo_SDSS_catalogue_Z_matches_good)
    SDSS_VLBI_PA = np.array(Astrogeo_SDSS_catalogue_JET_PA_matches_good)
    SDSS_VLBI_PA_ERR_ORIGINAL = np.array(Astrogeo_SDSS_catalogue_JET_PA_ERR_ORIGINAL_matches_good)
    SDSS_OPTICAL_PA = np.array(SDSS_catalogue_Astrogeo_THETA_J2000_matches)
    SDSS_OPTICAL_PA_ERR_ORIGINAL = np.array(SDSS_catalogue_Astrogeo_THETA_J2000_ERR_ORIGINAL_matches)
    SDSS_PA_DIFF = np.array(SDSS_catalogue_Astrogeo_PA_DIFFERENCE_matches)
    SDSS_PA_DIFF_ERR_ORIGINAL = np.array(SDSS_catalogue_Astrogeo_PA_DIFFERENCE_ERR_ORIGINAL_matches)
    
    # Return the processed data as a tuple
    return (
        SDSS_OPTICAL_PA, SDSS_OPTICAL_PA_ERR_ORIGINAL,
        SDSS_VLBI_PA, SDSS_VLBI_PA_ERR_ORIGINAL,
        SDSS_PA_DIFF, SDSS_PA_DIFF_ERR_ORIGINAL,
        SDSS_SOURCE_Z, SDSS_catalogue_Astrogeo_TYPE_matches,
        SDSS_catalogue_Astrogeo_B_matches
    )


def get_Skymapper_data(Astrogeo_Skymapper, Skymapper_xmatches):
    """
    Processes and extracts various data arrays from the Skymapper cross-matched results.

    This function extracts and processes relevant data from the `Skymapper_xmatches` object 
    and the `Astrogeo_Skymapper` DataFrame. It performs calculations to obtain arrays for 
    positional angles, differences, errors, and other parameters.

    Args:
        Astrogeo_Skymapper (pd.DataFrame): The filtered Astrogeo Skymapper catalogue containing data on positional angles and redshifts.
        Skymapper_xmatches (pd.DataFrame): The Skymapper cross-matched results containing various parameters.

    Returns:
        tuple: A tuple containing:
            - np.ndarray: Optical positional angles from the Skymapper catalogue.
            - np.ndarray: Errors on optical positional angles.
            - np.ndarray: VLBI positional angles from the Astrogeo catalogue.
            - np.ndarray: Errors on VLBI positional angles.
            - np.ndarray: Difference between optical and VLBI positional angles.
            - np.ndarray: Errors on the positional angle differences.
            - np.ndarray: Redshifts from the Astrogeo Skymapper catalogue.
            - np.ndarray: Types from Skymapper cross-matches.
            - np.ndarray: B-axis values from Skymapper cross-matches, scaled by 0.5.
    """
    
    # Extract arrays from Skymapper cross-matched results
    skymapper_Astrogeo_B_matches = 0.5 * np.array(Skymapper_xmatches.b)
    skymapper_Astrogeo_THETA_J2000_matches = np.array(Skymapper_xmatches.PA)
    skymapper_Astrogeo_THETA_J2000_ERR_ORIGINAL_matches = np.array(Skymapper_xmatches.e_pa)
    skymapper_Astrogeo_TYPE_matches = np.array(Skymapper_xmatches.ClassStar)
    
    # Extract arrays from Astrogeo Skymapper catalogue
    Astrogeo_skymapper_JET_PA_good_matches = np.array(Astrogeo_Skymapper.pa)
    Astrogeo_skymapper_JET_PA_ERR_ORIGINAL_good_matches = np.array(Astrogeo_Skymapper.pa_err)
    Astrogeo_skymapper_Z_good_matches = np.array(Astrogeo_Skymapper.Z)
    
    # Adjust the positional angles to the range [0, 180) degrees
    skymapper_Astrogeo_THETA_J2000_matches = np.fmod(skymapper_Astrogeo_THETA_J2000_matches + 360, 180)
    skymapper_Astrogeo_THETA_J2000_matches[skymapper_Astrogeo_THETA_J2000_matches > 90] -= 180
    
    # Compute the positional angle differences and errors
    skymapper_Astrogeo_PA_DIFFERENCE_matches = aux_functions.PA_difference(
        skymapper_Astrogeo_THETA_J2000_matches, Astrogeo_skymapper_JET_PA_good_matches
    )
    skymapper_Astrogeo_PA_DIFFERENCE_ERR_ORIGINAL_matches = np.sqrt(
        Astrogeo_skymapper_JET_PA_ERR_ORIGINAL_good_matches**2 +
        skymapper_Astrogeo_THETA_J2000_ERR_ORIGINAL_matches**2
    )
    
    # Prepare arrays for output
    skymapper_SOURCE_Z = np.array(Astrogeo_skymapper_Z_good_matches)
    skymapper_VLBI_PA = np.array(Astrogeo_skymapper_JET_PA_good_matches)
    skymapper_VLBI_PA_ERR_ORIGINAL = np.array(Astrogeo_skymapper_JET_PA_ERR_ORIGINAL_good_matches)
    skymapper_OPTICAL_PA = np.array(skymapper_Astrogeo_THETA_J2000_matches)
    skymapper_OPTICAL_PA_ERR_ORIGINAL = np.array(skymapper_Astrogeo_THETA_J2000_ERR_ORIGINAL_matches)
    skymapper_PA_DIFF = np.array(skymapper_Astrogeo_PA_DIFFERENCE_matches)
    skymapper_PA_DIFF_ERR_ORIGINAL = np.array(skymapper_Astrogeo_PA_DIFFERENCE_ERR_ORIGINAL_matches)
    
    # Return the processed data as a tuple
    return (
        skymapper_OPTICAL_PA, skymapper_OPTICAL_PA_ERR_ORIGINAL,
        skymapper_VLBI_PA, skymapper_VLBI_PA_ERR_ORIGINAL,
        skymapper_PA_DIFF, skymapper_PA_DIFF_ERR_ORIGINAL,
        skymapper_SOURCE_Z, skymapper_Astrogeo_TYPE_matches,
        skymapper_Astrogeo_B_matches
    )

def get_DES_data(Astrogeo_DES, DES_xmatches):
    """
    Extracts and processes data from Astrogeo and DES catalogues, including position angles (PAs),
    PA differences, redshifts, and classification data.

    Parameters:
    ----------
    Astrogeo_DES : pd.DataFrame
        DataFrame containing VLBI (Very Long Baseline Interferometry) data, including position angles (`pa`)
        and associated errors (`pa_err`).

    DES_xmatches : pd.DataFrame
        DataFrame containing cross-matched DES (Dark Energy Survey) data, including optical position angles,
        errors, redshift (`Z`), and other catalogue properties.

    Returns:
    -------
    tuple : 
        A tuple containing arrays with the following data:
        - Optical position angles from DES.
        - Errors in DES optical position angles.
        - VLBI position angles from Astrogeo.
        - Errors in Astrogeo VLBI position angles.
        - PA differences between DES and VLBI.
        - Errors in the PA difference.
        - Redshift values for matched sources.
        - Classification of sources as point or extended.
        - Optical flux densities of sources (scaled by 0.26).
    """
    
    # Extract optical position angle from DES matches
    DES_OPTICAL_PA = np.array(DES_xmatches.DES_full_catalogue_Astrogeo_THETA_J2000_matches)
    
    # Extract errors for DES optical position angles
    DES_OPTICAL_PA_ERR_ORIGINAL = np.array(DES_xmatches.DES_full_catalogue_Astrogeo_ERRTHETA_IMAGE_matches)
    
    # Extract VLBI position angle from Astrogeo data
    DES_VLBI_PA = np.array(Astrogeo_DES.pa)
    
    # Extract errors for VLBI position angles
    DES_VLBI_PA_ERR_ORIGINAL = np.array(Astrogeo_DES.pa_err)
    
    # Calculate position angle difference between DES optical and VLBI
    DES_PA_DIFF = np.array(DES_xmatches.DES_full_catalogue_Astrogeo_PA_DIFFERENCE_matches)
    
    # Calculate combined error in position angle difference (sum of squared errors)
    DES_PA_DIFF_ERR_ORIGINAL = np.sqrt(DES_VLBI_PA_ERR_ORIGINAL**2 + DES_OPTICAL_PA_ERR_ORIGINAL**2)
    
    # Extract redshift values for matched sources
    DES_SOURCE_Z = np.array(DES_xmatches.Astrogeo_DES_full_catalogue_Z_good_matches)
    
    # Extract flux densities of sources and scale them by a factor of 0.26
    DES_catalogue_Astrogeo_B_matches = 0.26 * np.array(DES_xmatches.DES_full_catalogue_Astrogeo_B_matches)
    
    # Extract classification of sources as point or extended
    DES_catalogue_Astrogeo_EXTENDED_CLASS_COADD_matches = np.array(DES_xmatches.DES_full_catalogue_Astrogeo_EXTENDED_CLASS_COADD_matches)
    
    return (DES_OPTICAL_PA, DES_OPTICAL_PA_ERR_ORIGINAL, DES_VLBI_PA, DES_VLBI_PA_ERR_ORIGINAL, 
            DES_PA_DIFF, DES_PA_DIFF_ERR_ORIGINAL, DES_SOURCE_Z, 
            DES_catalogue_Astrogeo_EXTENDED_CLASS_COADD_matches, DES_catalogue_Astrogeo_B_matches)


def get_Skymapper_data(Astrogeo_Skymapper, Skymapper_xmatches):
    """
    Extracts and processes SkyMapper and Astrogeo data, including position angles (PAs), PA differences, 
    redshifts, and other properties from the SkyMapper-Astrogeo cross-matched catalogues.

    Parameters:
    ----------
    Astrogeo_Skymapper : pd.DataFrame
        DataFrame containing VLBI (Very Long Baseline Interferometry) data from the Astrogeo catalogue, 
        including position angles (`pa`), errors (`pa_err`), and redshift (`Z`).

    Skymapper_xmatches : pd.DataFrame
        DataFrame containing cross-matched SkyMapper data, including optical position angles, errors, 
        and other catalogue properties such as flux densities and object classification.

    Returns:
    -------
    tuple : 
        A tuple containing NumPy arrays with the following data:
        - Optical position angles from SkyMapper.
        - Errors in the optical position angles from SkyMapper.
        - VLBI position angles from Astrogeo.
        - Errors in the VLBI position angles.
        - PA differences between SkyMapper optical and Astrogeo VLBI.
        - Combined errors in the PA difference (sum of squared errors).
        - Redshift values for the matched sources.
        - Object classification from SkyMapper (e.g., point or extended sources).
        - Optical flux densities of the sources, scaled by a factor of 0.5.
    """
    
    # Scale the optical flux density by 0.5
    skymapper_Astrogeo_B_matches = 0.5 * np.array(Skymapper_xmatches.b)
    
    # Extract SkyMapper optical position angles and errors
    skymapper_Astrogeo_THETA_J2000_matches = np.array(Skymapper_xmatches.PA)
    skymapper_Astrogeo_THETA_J2000_ERR_ORIGINAL_matches = np.array(Skymapper_xmatches.e_pa)
    
    # Extract VLBI position angles and errors from Astrogeo
    Astrogeo_skymapper_JET_PA_good_matches = np.array(Astrogeo_Skymapper.pa)
    Astrogeo_skymapper_JET_PA_ERR_ORIGINAL_good_matches = np.array(Astrogeo_Skymapper.pa_err)
    
    # Extract SkyMapper object classification and redshift values
    skymapper_Astrogeo_TYPE_matches = np.array(Skymapper_xmatches.ClassStar)
    Astrogeo_skymapper_Z_good_matches = np.array(Astrogeo_Skymapper.Z)
    
    # Adjust the optical position angles to fit within [0, 180]
    skymapper_Astrogeo_THETA_J2000_matches = np.fmod(skymapper_Astrogeo_THETA_J2000_matches + 360, 180)
    skymapper_Astrogeo_THETA_J2000_matches[skymapper_Astrogeo_THETA_J2000_matches > 90] -= 180
    
    # Calculate PA difference between SkyMapper and VLBI
    skymapper_Astrogeo_PA_DIFFERENCE_matches = aux_functions.PA_difference(skymapper_Astrogeo_THETA_J2000_matches, Astrogeo_skymapper_JET_PA_good_matches)
    
    # Compute the combined error in the PA difference
    skymapper_Astrogeo_PA_DIFFERENCE_ERR_ORIGINAL_matches = np.sqrt(
        Astrogeo_skymapper_JET_PA_ERR_ORIGINAL_good_matches**2 + 
        skymapper_Astrogeo_THETA_J2000_ERR_ORIGINAL_matches**2
    )
    
    # Extract relevant data for return
    skymapper_SOURCE_Z = np.array(Astrogeo_skymapper_Z_good_matches)
    skymapper_VLBI_PA = np.array(Astrogeo_skymapper_JET_PA_good_matches)
    skymapper_VLBI_PA_ERR_ORIGINAL = np.array(Astrogeo_skymapper_JET_PA_ERR_ORIGINAL_good_matches)
    skymapper_OPTICAL_PA = np.array(skymapper_Astrogeo_THETA_J2000_matches)
    skymapper_OPTICAL_PA_ERR_ORIGINAL = np.array(skymapper_Astrogeo_THETA_J2000_ERR_ORIGINAL_matches)
    skymapper_PA_DIFF = np.array(skymapper_Astrogeo_PA_DIFFERENCE_matches)
    skymapper_PA_DIFF_ERR_ORIGINAL = np.array(skymapper_Astrogeo_PA_DIFFERENCE_ERR_ORIGINAL_matches)
    
    return (skymapper_OPTICAL_PA, skymapper_OPTICAL_PA_ERR_ORIGINAL, skymapper_VLBI_PA, skymapper_VLBI_PA_ERR_ORIGINAL, 
            skymapper_PA_DIFF, skymapper_PA_DIFF_ERR_ORIGINAL, skymapper_SOURCE_Z, 
            skymapper_Astrogeo_TYPE_matches, skymapper_Astrogeo_B_matches)

def get_Combined_data(Surveys_Data, Surveys_Good_Data):
    """
    Extracts and combines data from two survey datasets based on the availability of weighted average jet minor differences.

    This function extracts arrays of various parameters from the `Surveys_Data` and `Surveys_Good_Data` DataFrames.
    It filters the data to exclude entries where the weighted average jet minor differences are missing and combines 
    the relevant data into arrays for further analysis.

    Args:
        Surveys_Data (pd.DataFrame): The DataFrame containing data from the surveys, including positional angles and errors.
        Surveys_Good_Data (pd.DataFrame): The DataFrame containing additional good data from the surveys.

    Returns:
        tuple: A tuple containing:
            - np.ndarray: VLBI positional angles from the Surveys_Data.
            - np.ndarray: Errors on VLBI positional angles from the Surveys_Data.
            - np.ndarray: VLBI positional angles from the Surveys_Good_Data.
            - np.ndarray: Errors on VLBI positional angles from the Surveys_Good_Data.
            - np.ndarray: Optical positional angles from the Surveys_Data.
            - np.ndarray: Errors on optical positional angles from the Surveys_Data.
            - np.ndarray: Optical positional angles from the Surveys_Good_Data.
            - np.ndarray: Errors on optical positional angles from the Surveys_Good_Data.
            - np.ndarray: Jet minor differences from the Surveys_Data.
            - np.ndarray: Errors on jet minor differences from the Surveys_Data.
            - np.ndarray: Jet minor differences from the Surveys_Good_Data.
            - np.ndarray: Errors on jet minor differences from the Surveys_Good_Data.
            - np.ndarray: Redshifts from the Surveys_Data.
            - np.ndarray: Redshifts from the Surveys_Good_Data.
    """
    
    # Filter and extract arrays based on the availability of weighted average jet minor differences
    valid_indices_data = pd.notnull(Surveys_Data.Weighted_average_Jet_Minor_diffs)
    valid_indices_good_data = pd.notnull(Surveys_Good_Data.Weighted_average_Jet_Minor_diffs)
    
    Combined_VLBI_PA = np.array(Surveys_Data.VLBI_Jet_PAs)[valid_indices_data]
    Combined_VLBI_PA_ERR_ORIGINAL = np.array(Surveys_Data.VLBI_Jet_PA_errors)[valid_indices_data]
    Combined_Good_VLBI_PA = np.array(Surveys_Good_Data.VLBI_Jet_PAs)[valid_indices_good_data]
    Combined_Good_VLBI_PA_ERR_ORIGINAL = np.array(Surveys_Good_Data.VLBI_Jet_PA_errors)[valid_indices_good_data]
    
    Combined_OPTICAL_PA = np.array(Surveys_Data.Weighted_average_Optical_Major_PAs)[valid_indices_data]
    Combined_OPTICAL_PA_ERR_ORIGINAL = np.array(Surveys_Data.Weighted_average_Optical_Major_PA_errors)[valid_indices_data]
    Combined_Good_OPTICAL_PA = np.array(Surveys_Good_Data.Weighted_average_Optical_Major_PAs)[valid_indices_good_data]
    Combined_Good_OPTICAL_PA_ERR_ORIGINAL = np.array(Surveys_Good_Data.Weighted_average_Optical_Major_PA_errors)[valid_indices_good_data]
    
    Combined_PA_DIFF = np.array(Surveys_Data.Weighted_average_Jet_Minor_diffs)[valid_indices_data]
    Combined_PA_DIFF_ERR_ORIGINAL = np.array(Surveys_Data.Weighted_average_Jet_Minor_diffs_errors)[valid_indices_data]
    Combined_Good_PA_DIFF = np.array(Surveys_Good_Data.Weighted_average_Jet_Minor_diffs)[valid_indices_good_data]
    Combined_Good_PA_DIFF_ERR_ORIGINAL = np.array(Surveys_Good_Data.Weighted_average_Jet_Minor_diffs_errors)[valid_indices_good_data]
    
    Combined_SOURCE_Z = np.array(Surveys_Data.Z)[valid_indices_data]
    Combined_Good_SOURCE_Z = np.array(Surveys_Good_Data.Z)[valid_indices_good_data]
    
    return (
        Combined_VLBI_PA, Combined_VLBI_PA_ERR_ORIGINAL,
        Combined_Good_VLBI_PA, Combined_Good_VLBI_PA_ERR_ORIGINAL,
        Combined_OPTICAL_PA, Combined_OPTICAL_PA_ERR_ORIGINAL,
        Combined_Good_OPTICAL_PA, Combined_Good_OPTICAL_PA_ERR_ORIGINAL,
        Combined_PA_DIFF, Combined_PA_DIFF_ERR_ORIGINAL,
        Combined_Good_PA_DIFF, Combined_Good_PA_DIFF_ERR_ORIGINAL,
        Combined_SOURCE_Z, Combined_Good_SOURCE_Z
    ) 