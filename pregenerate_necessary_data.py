
#imports
import numpy as np
import data_loading_and_xmatching


#Load the Astrogeo VLBI catalogue
Astrogeo_catalogue = data_loading_and_xmatching.load_Astrogeo_catalogue()
Astrogeo_good_catalogue = data_loading_and_xmatching.get_Astrogeo_good_catalogue(Astrogeo_catalogue)

#Load the whole optical catalogues
DESI_catalogue,SDSS_catalogue,DES_catalogue,Skymapper_catalogue,KIDS_catalogue = data_loading_and_xmatching.load_catalogues('All')ç


#Perform cross-matches bewteen the VLBI and each optical survey
data_loading_and_xmatching.perform_xmatch(Astrogeo_good_catalogue,DESI_catalogue,'DESI','RA','DEC')
data_loading_and_xmatching.perform_xmatch(Astrogeo_good_catalogue,SDSS_catalogue,'SDSS','ra','dec')
data_loading_and_xmatching.perform_xmatch(Astrogeo_good_catalogue,DES_catalogue,'DES','DES_full_catalogue_Astrogeo_matches_RA','DES_full_catalogue_Astrogeo_matches_DEC')
data_loading_and_xmatching.perform_xmatch(Astrogeo_good_catalogue,Skymapper_catalogue,'Skymapper','ra','dec')
data_loading_and_xmatching.perform_xmatch(Astrogeo_good_catalogue,KIDS_catalogue,'KIDS','RAJ2000','DECJ2000')

#Create a random sample of SDSS and KIDS for Figure 2
size = 20000
data_loading_and_xmatching.generate_sampled_catalogues(size,SDSS_catalogue,KIDS_catalogue)