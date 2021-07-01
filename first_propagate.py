import pandas as pd
import csv
import numpy as np
import math as m
from astropy.coordinates import SkyCoord
import astropy.units as u
from multiprocessing import Pool


def check_poles(row,t):

	
	new_l = (row['l'] + t*row['pm_l']/(3600*1000))%360
	new_b = row['b'] + t*row['pm_b']/(3600*1000)
	#new_l = (row['l'] + t*row['pm_l']/(3600*1000))%360
	#new_b = row['b'] + t*row['pm_b']/(3600*1000)
	
	if abs(new_b) >= 90:
		#change to equatorials and propagate
		coord_gal = SkyCoord(l = row['l']*u.degree,b = row['b']*u.degree, frame='galactic')
		

		ra = coord_gal.icrs.ra.value
		dec = coord_gal.icrs.dec.value
		new_ra = (ra+ t*row['pmra']/(3600*1000))%360
		new_dec = dec + t*row['pmdec']/(3600*1000)

		coord = SkyCoord(ra = new_ra*u.degree,dec = new_dec*u.degree)
		gala = coord.galactic
		new_b = gala.b.value

		
	return new_l,new_b

def apply(df):
	df[['new_l','new_b']] = df.apply(check_poles,axis = 1,result_type = "expand",args=(t,))
	#df[['new_l','new_b']] = df.apply(check_poles,axis = 1,args=(t,))
	return df


def parallelize_dataframe(df,func, n_cores):
	df_split = np.array_split(df,n_cores)
	pool = Pool(n_cores)
	df = pd.concat(pool.map(func,df_split))
	pool.close()
	pool.join()
	return df

# Preview the first 5 lines of the loaded data 
#test = data.head(50)


"""Prova per mirar si calcula 5E4 a partir de 4E4"""

deltat = -1E4
n_cores = 6
data_0 = pd.read_csv("/home/hpc/atena/jnogue/Documents/Gaia/Data_error_5/new_data_error_5.csv")
t = deltat
new_df = parallelize_dataframe(data_0,apply,n_cores)
new_df.to_csv('/home/hpc/atena/jnogue/Documents/Gaia/Data_error_5/-100000.csv', index=False)