from libraries import *
import aux_functions

def load_Astrogeo_catalogue():
    """
    Loads the Astrogeo catalogue from a CSV file.

    This function reads a CSV file called 'Source_coords_csv.csv' which contains 
    the Astrogeo catalogue and loads it into a pandas DataFrame.

    Returns:
        pd.DataFrame: A pandas DataFrame containing the data from the Astrogeo catalogue.
    """
    # Print a message indicating the loading process has started
    print('Loading Astrogeo...')
    
    # Load the CSV file into a pandas DataFrame
    Astrogeo_catalogue = pd.read_csv('data/Astrogeo_catalogue.csv')
    
    # Return the loaded DataFrame
    return Astrogeo_catalogue

def get_Astrogeo_good_catalogue(Astrogeo_catalogue):
    """
    Filters the Astrogeo catalogue to return only entries with valid RA and DEC values.

    This function checks the Right Ascension (RA) and Declination (DEC) values in the 
    Astrogeo catalogue and filters out rows with invalid or NaN values. The valid 
    range for RA is between 0 and 360 degrees, and for DEC it is between -90 and 90 degrees.

    Args:
        Astrogeo_catalogue (pd.DataFrame): The full Astrogeo catalogue loaded from a CSV file.

    Returns:
        pd.DataFrame: A pandas DataFrame containing only the entries with valid RA and DEC values.
    """
    # Convert the RA and DEC columns to numpy arrays for easier filtering
    Astrogeo_RA = np.array(Astrogeo_catalogue.RA_deg)
    Astrogeo_DEC = np.array(Astrogeo_catalogue.DEC_deg)
    
    # Identify bad RA values: values outside the 0-360 range or NaN values
    RA_bad = (((Astrogeo_RA > 0) & (Astrogeo_RA < 360)) == False) & (np.isnan(Astrogeo_RA) == True)
    
    # Identify bad DEC values: values outside the -90 to 90 range or NaN values
    DEC_bad = (((Astrogeo_DEC > (-90)) & (Astrogeo_DEC < 90)) == False) & (np.isnan(Astrogeo_DEC) == True)
    
    # Identify good RA values: values within the 0-360 range and not NaN
    RA_good = (((Astrogeo_RA > 0) & (Astrogeo_RA < 360)) == True) & (np.isnan(Astrogeo_RA) == False)
    
    # Identify good DEC values: values within the -90 to 90 range and not NaN
    DEC_good = (((Astrogeo_DEC > (-90)) & (Astrogeo_DEC < 90)) == True) & (np.isnan(Astrogeo_DEC) == False)
    
    # Filter the catalogue to include only rows where both RA and DEC are good
    Astrogeo_good = Astrogeo_catalogue[RA_good & DEC_good]
    
    # Return the filtered catalogue
    return Astrogeo_good

    

def load_catalogues(catalogue):
    """
    Loads the specified astronomical catalogue based on the input parameter.

    Depending on the catalogue name provided, this function loads the corresponding 
    CSV or FITS file into memory. If 'All' is provided as the catalogue name, it 
    will load multiple catalogues simultaneously and return them as a tuple.

    Args:
        catalogue (str): The name of the catalogue to load. Possible values are:
                         'DESI', 'SDSS', 'DES', 'Skymapper', 'KIDS', or 'All'.

    Returns:
        pd.DataFrame or tuple: A pandas DataFrame for the specific catalogue, or a tuple 
        of DataFrames and FITS data if 'All' is selected.
    """
    if catalogue == 'DESI':
        # Loading the DESI catalogue from a CSV file
        print('Loading DESI...')
        DESI_catalogue = pd.read_csv('data/DESI_catalogue.csv', header=0)
        return DESI_catalogue
    
    if catalogue == 'SDSS':
        # Loading the SDSS catalogue from a CSV file
        print('Loading SDSS...')
        SDSS_catalogue = pd.read_csv('data/SDSS_catalogue.csv', header=0)
        return SDSS_catalogue
    
    if catalogue == 'DES':
        # Loading the DES catalogue from a CSV file
        print('Loading DES...')
        DES_catalogue = pd.read_csv('data/DES_full_catalogue_Astrogeo.csv', header=0)
        return DES_catalogue
    
    if catalogue == 'Skymapper':
        # Loading the Skymapper catalogue from a CSV file
        print('Loading Skymapper...')
        Skymapper_catalogue = pd.read_csv('data/Skymapper_catalogue.csv', header=0)
        return Skymapper_catalogue
    
    if catalogue == 'KIDS':
        # Loading the KIDS catalogue from a FITS file
        print('Loading KIDS...')
        with fits.open('data/KiDS_catalogue.fits', memmap=True) as KIDS_fits:
            KIDS_catalogue = KIDS_fits[1].data  # Access the data from the first extension
        return KIDS_catalogue
    
    if catalogue == 'All':
        # Loading all the catalogues and returning them as a tuple
        print('Loading DESI...')
        DESI_catalogue = pd.read_csv('data/DESI_catalogue.csv', header=0)
        
        print('Loading SDSS...')
        SDSS_catalogue = pd.read_csv('data/SDSS_catalogue.csv', header=0)
        
        print('Loading DES...')
        DES_catalogue = pd.read_csv('data/DES_catalogue.csv', header=0)
        
        print('Loading Skymapper...')
        Skymapper_catalogue = pd.read_csv('data/Skymapper_catalogue.csv', header=0)
        
        print('Loading KIDS...')
        with fits.open('data/KiDS_catalogue.fits', memmap=True) as KIDS_fits:
            KIDS_catalogue = KIDS_fits[1].data
        
        # Return all loaded catalogues as a tuple
        return DESI_catalogue, SDSS_catalogue, DES_catalogue, Skymapper_catalogue, KIDS_catalogue

def generate_sampled_catalogues(size, SDSS_catalogue, KIDS_catalogue):
    """
    Generates a random sample from the SDSS and KIDS catalogues, to latyer use in Figure 2.

    This function takes the full SDSS and KIDS catalogues, and generates random samples 
    of the specified size from each. The samples are selected without replacement.

    Args:
        size (int): The number of random samples to extract from each catalogue.
        SDSS_catalogue (pd.DataFrame): The full SDSS catalogue in pandas DataFrame format.
        KIDS_catalogue (array or DataFrame): The full KIDS catalogue, typically an array from a FITS file.

    Saves:
        - pd.DataFrame: Randomly sampled rows from the SDSS catalogue.
        - np.ndarray or pd.DataFrame: Randomly sampled rows from the KIDS catalogue.
    """
    
    # Load a random sample from the SDSS catalogue
    print('Generating sample of SDSS...')
    
    # Generate an array of indices for the SDSS catalogue
    SDSS_index = np.arange(0, len(SDSS_catalogue))
    
    # Randomly select 'size' number of indices from the SDSS index array
    SDSS_random_index = np.array(random.sample(list(SDSS_index), k=size))
    
    # Use the randomly selected indices to get a sample from the SDSS catalogue
    SDSS_sampled_catalogue = SDSS_catalogue.iloc[SDSS_random_index]
    
    # Load a random sample from the KIDS catalogue
    print('Generating sample of KIDS...')
    
    # Generate an array of indices for the KIDS catalogue
    KIDS_index = np.arange(0, len(KIDS_catalogue))
    
    # Randomly select 'size' number of indices from the KIDS index array
    KIDS_random_index = np.array(random.sample(list(KIDS_index), k=size))
    
    # Use the randomly selected indices to get a sample from the KIDS catalogue
    KIDS_sampled_catalogue = Table(KIDS_catalogue[KIDS_random_index])
    
    # Save the sampled catalogues as a pd.DataFrame
    SDSS_sampled_catalogue.to_csv('data/SDSS_sampled_catalogue.csv', index=False) 
    KIDS_sampled_catalogue.write('data/KIDS_sampled_catalogue.csv', delimiter=',', format='ascii') 


def perform_xmatch(Astrogeo_catalogue, catalogue, catalogue_name, catalogue_RA_name, catalogue_DEC_name):
    """
    Performs a cross-match between the Astrogeo catalogue and a given external catalogue.

    This function matches the coordinates of the `Astrogeo_good` catalogue with those of the 
    external catalogue specified by `catalogue_name`. It uses a maximum separation of 1 arcsecond 
    to find the closest matches.

    Args:
        Astrogeo_good (pd.DataFrame): The filtered Astrogeo catalogue with valid coordinates.
        catalogue (pd.DataFrame or array-like): The external catalogue to be matched with the Astrogeo catalogue.
        catalogue_name (str): The name of the external catalogue (e.g., 'KIDS').
        catalogue_RA_name (str): The name of the column in `catalogue` that contains RA values.
        catalogue_DEC_name (str): The name of the column in `catalogue` that contains DEC values.

    Saves:
        - pd.DataFrame: The subset of `Astrogeo_good` that has matches in the external catalogue.
        - pd.DataFrame or array-like: The subset of the external catalogue that matches the `Astrogeo_good` entries.
    """
    print('Performing xmatch on ' + catalogue_name + '...')
    
    # Extract RA and DEC from the Astrogeo catalogue and convert to SkyCoord objects
    Astrogeo_RA_good = np.array(Astrogeo_good.RA_deg)
    Astrogeo_DEC_good = np.array(Astrogeo_good.DEC_deg)
    Astrogeo_coords_good = SkyCoord(Astrogeo_RA_good, Astrogeo_DEC_good, frame='icrs', unit='deg')
    
    # Extract RA and DEC from the external catalogue and convert to SkyCoord objects
    catalogue_RA = catalogue[catalogue_RA_name]
    catalogue_DEC = catalogue[catalogue_DEC_name]
    catalogue_coords = SkyCoord(catalogue_RA, catalogue_DEC, frame='icrs', unit='deg')
    
    # Define the maximum separation for matches
    max_sep = 1 * u.arcsec
    
    # Perform the match and obtain indices and separation distances
    idx, d2d, d3d = Astrogeo_coords_good.match_to_catalog_sky(catalogue_coords)
    
    # Apply the separation constraint to filter matches
    sep_constraint = d2d < max_sep
    
    # Get the matched rows from the Astrogeo catalogue
    Astrogeo_catalogue_matches = Astrogeo_good[sep_constraint]
    
    # Get the matched rows from the external catalogue
    if catalogue_name == 'KIDS':
        catalogue_Astrogeo_matches = Table(catalogue[idx[sep_constraint]])
        catalogue_Astrogeo_matches.write('data/'+catalogue_name+'_xmatches'+'.csv', delimiter=',', format='ascii')
    else:
        catalogue_Astrogeo_matches = catalogue.iloc[idx[sep_constraint]]
        catalogue_Astrogeo_matches.to_csv('data/'+catalogue_name+'_xmatches'+'.csv',index=False)
    Astrogeo_catalogue_matches.to_csv('data/Astrogeo_'+catalogue_name+'.csv',index=False)
    
    
def load_sampled_catalogues():
    """
    Load the SDSS and KIDS sampled catalogues from CSV files.

    Returns:
    --------
    SDSS_sampled_catalogue : pandas.DataFrame
        Data from the SDSS sampled catalogue.
    KIDS_sampled_catalogue : pandas.DataFrame
        Data from the KIDS sampled catalogue.
    """
    SDSS_sampled_catalogue = pd.read_csv('data/SDSS_sampled_catalogue.csv')
    KIDS_sampled_catalogue = pd.read_csv('data/KIDS_sampled_catalogue.csv')
    
    return SDSS_sampled_catalogue, KIDS_sampled_catalogue
    
def load_xmatches(catalogue_name):
    """
    Load Astrogeo and cross-matches catalogues based on the provided catalogue name.

    Parameters:
    -----------
    catalogue_name : str
        The name of the catalogue to load (e.g., 'SDSS', 'DESI').

    Returns:
    --------
    Astrogeo_catalogue : pandas.DataFrame
        Data from the Astrogeo catalogue corresponding to the given catalogue name.
    catalogue_xmatches : pandas.DataFrame
        Cross-matches data for the given catalogue.
    """
    print('Loading '+catalogue_name+'...')
    Astrogeo_catalogue = pd.read_csv(f'data/Astrogeo_{catalogue_name}.csv')
    catalogue_xmatches = pd.read_csv(f'data/{catalogue_name}_xmatches.csv')
    
    return Astrogeo_catalogue, catalogue_xmatches

def load_Eagle_sim():
    """
    Loads the Eagle simulation data from a FITS file.

    This function reads the Eagle simulation data from the FITS file 'eagle_new.fits' 
    and converts it into an Astropy Table for easier manipulation and analysis.

    Returns:
        astropy.table.Table: The Eagle simulation data as an Astropy Table.
    """
    print('Loading Eagle simulation...')
    
    # Open the FITS file and access the data from the first extension
    with fits.open('data/Eagle_sim.fits') as Eagle_fits:
        Eagle_sim = Eagle_fits[1].data
    
    # Convert the data to an Astropy Table
    Eagle_sim = Table(Eagle_sim)
    
    # Return the Astropy Table containing the Eagle simulation data
    return Eagle_sim
