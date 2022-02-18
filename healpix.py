import healpy
from astropy_healpix import HEALPix
import pandas as pd
import numpy as np
from astropy . coordinates import SkyCoord
import astropy . units as u
from multiprocessing import Pool
import matplotlib . pyplot as plt

def propagate (row ,hpx ,t):
	coord = SkyCoord (l=row['l']*u.deg , b=row['b']*u.deg ,pm_l_cosb =row['pm_l']*u.mas/u.yr ,pm_b = row['pm_b']*u.mas /u.yr , distance =200* u.kpc , frame =’galactic ’)
	coord = coord . apply_space_motion (dt=t*u.yr)
	new_b = coord . galactic .b. value
	new_l = coord . galactic .l. value
	pm_l = coord . galactic . pm_l_cosb . value
	pm_b = coord . galactic . pm_b . value
	# Assigns to each pair of coordinates a healpix
	h8 = hpx . lonlat_to_healpix ( new_l *u.deg , new_b *u.deg )
	return new_l ,new_b , pm_l , pm_b ,h8

def apply (df):
# modified apply function where a healpix8 column is added to the file
	df[['l','b','pm_l','pm_b ','healpix8 ']] = df. apply (propagate , axis = 1 , result_type = " expand ",args =( hpx ,t))
	return df

def parallelize_dataframe (df ,func , n_cores ):
	#parallelize dataframe for multiprocessing
	df_split = np.array_split (df , n_cores )
	pool = Pool (n_cores)
	df = pd.concat (pool.map (func , df_split ) )
	pool.close ()
	pool.join ()
	return df


# Healpix characteristics , more level equals more
resolution
level =7
hpx = HEALPix (nside=4096/2**(12 - level ) , order ='nested')
num_hp = int( healpy . nside2npix ( nside =4096/2**(12 - level ) ) )
area_hp =(4* np.pi/ num_hp ) #in sr
nside =4096/2**(12 - level )

# Here begins the propagation : Reads the csv file for 100000 years , does 91 steps of 10000 years and computes the csv file for 1 milion years .
t0 = 1E5
deltat = 1E4
n_cores = 6
data_0 = pd. read_csv (" path_to_file /100000. csv ")

for i in range (1 ,91) :
	t_past = int(t0 + (i -1) * deltat )
	print ( t_past )
	t = t0 + i* deltat
	filename = str(int(t) )
	print ( filename )
	new_df = parallelize_dataframe (data_0 ,apply , n_cores )
	# HEALPIX BEGINS
	healpixRVals = dict ()
	healpixGVals = dict ()
	healpixBVals = dict ()
	healpixCount = dict ()
	for index , row in new_df . iterrows () :
		if index % 10000 == 0:
			print (str( index ) + " of " + str(len( new_df ) ) )
		healpix = row ['healpix8 ']

		if healpix in healpixRVals :
			healpixRVals [ healpix ] += row ['phot_rp_mean_flux']
		else :
			healpixRVals [ healpix ] = row ['phot_rp_mean_flux']

		if healpix in healpixGVals :
			healpixGVals [ healpix ] += row ['phot_g_mean_flux']
		else :
			healpixGVals [ healpix ] = row ['phot_g_mean_flux']
		if healpix in healpixBVals :
			healpixBVals [ healpix ] += row ['phot_bp_mean_flux']
		else :
			healpixBVals [ healpix ] = row ['phot_bp_mean_flux']
		if healpix in healpixCount :
			healpixCount [ healpix ] += 1
		else :
			healpixCount [ healpix ] = 1


	r_list = []
	b_list = []
	g_list = []
	count_list = []
	for healpix_num in range ( num_hp ) :
		if healpix_num in healpixRVals and not np. isnan (healpixRVals [ healpix_num ]) :
			r_val = healpixRVals [ healpix_num ]
		else :
			r_val = 0

		if healpix_num in healpixGVals and not np. isnan (healpixGVals [ healpix_num ]) :
			g_val = healpixGVals [ healpix_num ]
		else :
			g_val = 0

		if healpix_num in healpixBVals and not np. isnan (healpixBVals [ healpix_num ]) :
			b_val = healpixBVals [ healpix_num ]
		else :
			b_val = 0

		if healpix_num in healpixCount and not np. isnan (healpixCount [ healpix_num ]) :
			count_val = healpixCount [ healpix_num ]
		else :
			count_val = 0
		r_list . append ( r_val )
		g_list . append ( g_val )
		b_list . append ( b_val )
		count_list . append ( count_val)

	fonts =12
	fig = plt . figure ( figsize =(8 ,5) )
	ax = healpy . projaxes . HpxMollweideAxes (fig,[0.1 ,0.1 ,0.9 ,0.9] , rot =(0 ,0.0 ,0.0) ,coord =['G'])
	arrR =ax. projmap (np. array ( r_list ) ,nest =True , coord =['G'])
	arrG =ax. projmap (np. array ( g_list ) ,nest =True , coord =['G'])
	arrB =ax. projmap (np. array ( b_list ) ,nest =True , coord =['G'])
	arrCount =ax. projmap (np. array ( count_list ) ,nest =True , coord =['G'])
	# Transformations explained in equation (6)
	arrR = 24.7479 -2.5* np.log ( arrR )
	arrG = 25.6085 -2.5* np.log ( arrG )
	arrB = 25.3385 -2.5* np.log ( arrB )

	# Reescalte to obtain a value of magnitude between 0 and 1 (RGB scale )
	arrR = -(arrR -np.min ( arrR ) *np. ones (np. shape ( arrR ) ) ) /( np.max ( arrR ) -np.min( arrR ) ) +np. ones (np. shape ( arrR ) )
	arrG = -(arrG -np.min ( arrG ) *np. ones (np. shape ( arrG ) ) ) /( np.max ( arrG ) -np.min( arrG ) ) +np. ones (np. shape ( arrR ) )
	arrB = -(arrB -np.min ( arrB ) *np. ones (np. shape ( arrB ) ) ) /( np.max ( arrB ) -np.min( arrB ) ) +np. ones (np. shape ( arrR ) )

	rgbArr = np. stack (( arrR ,arrG , arrB ) , axis =2 , out = None )
	## Color Map
	plt . figure ( figsize = (10 ,6) )
	fig = plt . imshow (rgbArr , extent =(180 , -180 ,90 , -90) , aspect ='auto ')
	ax = plt .gca ()
	ax. set_title (str (int (t) ) +' years ',fontsize = fonts +5)
	ax. invert_yaxis ()
	ax. set_xlabel (r'$\ell$ (deg)',fontsize = fonts )
	ax. set_ylabel (r'$b$ (deg)',fontsize = fonts )

	plt . savefig (' path_to_file / healpix_color_hp7 '+str(int(t) ) +'. png ',format ='png',dpi =200)
	plt . close ()
	## Density Map
	plt . figure ( figsize = (10 ,6) )
	fig = plt . imshow (np.log ( arrCount ) ,'jet',extent=(180 , -180 ,90 , -90) , aspect ='auto')
	ax = plt .gca ()
	ax. set_title (str (int (t) ) +'years' ,fontsize = fonts +5)
	ax. invert_yaxis ()
	ax. set_xlabel (r'$ \ ell$ (deg)',fontsize = fonts )
	ax. set_ylabel (r'$b$ (deg)',fontsize = fonts )
	plt . savefig ('path_to_file / healpix_count_hp7 '+str(int(t) ) +'. png ',format ='png',dpi =200)
	plt . close ()
	# Obtain final csv
	new_df . to_csv ('path_to_file /'+ filename +'.csv', index =False )
	print (" Fitxer : "+ filename +'.csv fet ')
	data_0 = new_df
