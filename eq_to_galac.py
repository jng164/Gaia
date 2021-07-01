import pandas as pd
import csv
import numpy as np
import math as m
from astropy.coordinates import SkyCoord
import astropy.units as u
from multiprocessing import Pool

def equator_to_galac(row):

	coord = SkyCoord(ra = row['ra']*u.degree,dec = row['dec']*u.degree,
                  pm_ra_cosdec = row['pmra']*u.mas/u.yr, 
      
                  pm_dec = row['pmdec']*u.mas/u.yr)
	gala = coord.galactic

	return (gala.pm_l_cosb.value, gala.pm_b.value)

def apply(df):
	df[['pm_l','pm_b']] = df.apply(equator_to_galac,axis = 1,result_type = "expand")
	
	return df


def parallelize_dataframe(df,func, n_cores):
	df_split = np.array_split(df,n_cores)
	pool = Pool(n_cores)
	df = pd.concat(pool.map(func,df_split))
	pool.close()
	pool.join()
	return df

data = pd.read_csv("/home/hpc/atena/jnogue/Documents/Gaia/Data_error_5/data_error_5.csv")
# Preview the first 5 lines of the loaded data 
#test = data.head(20)
test = data
n_cores = 6

	
new_df = parallelize_dataframe(test,apply,n_cores)

new_df.to_csv('/home/hpc/atena/jnogue/Documents/Gaia/Data_error_5/new_data_error_5.csv', index=False)

#data[['pm_l','pm_b']] = data.apply(equator_to_galac,axis = 1,result_type = "expand")
#data.to_csv('test_tot.csv', index=False)
