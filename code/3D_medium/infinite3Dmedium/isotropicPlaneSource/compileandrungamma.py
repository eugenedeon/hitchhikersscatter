# Compile and run simulations: 

# 3D INFINITE MEDIUM ISOTROPIC PLANE SOURCE

import sys
import os
import math

if not os.path.exists('data'):
	os.mkdir('data')

# classical - various phase functions
#os.system('g++ infinite3D_isotropicplane_isotropicscatter_gamma2.cpp -o infinite3D_isotropicplane_isotropicscatter_gamma2 -O3 -I ../../../include/')
#os.system('g++ infinite3D_isotropicplane_isotropicscatter_gamma2C_tracklength.cpp -o infinite3D_isotropicplane_isotropicscatter_gamma2C_tracklength -O3 -I ../../../include/')
#os.system('g++ infinite3D_isotropicplane_isotropicscatter_gamma2U_tracklength.cpp -o infinite3D_isotropicplane_isotropicscatter_gamma2U_tracklength -O3 -I ../../../include/')
os.system('g++ infinite3D_isotropicplane_HG_exponential.cpp -o infinite3D_isotropicplane_HG_exponential -O3 -I ../../../include/')

#cs = [ 0.01, 0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99, 0.999 ]
cs = [ 0.9 ]

mfp = 2

for g in [0.0]:
	for c in cs:
		expected_num_samples_inv = math.sqrt( ( 1.0 - c ) / ( c * c ) )
		numsamples = min( 10000000, int( 10000000.0 * expected_num_samples_inv ) )
		print 'num samples: ' + str(numsamples)
		maxr = 25.0 * mfp
		dr = 0.2 * mfp
		du = 2.0 / 4.0
		numorders = 10
		nummoments = 8
		#filename = 'data/infinite3D_isotropicplane_isotropicscatter_gamma2C_tracklength_c' + str(c) + '_mut' + str(mut) + '.txt'
		#print 'computing: ' + filename
		#os.system( './infinite3D_isotropicplane_isotropicscatter_gamma2C_tracklength ' + str(c) + ' ' + str(mut) + ' ' + str(maxr) + ' ' + str(dr) + ' ' + str(du) + ' ' + str( numsamples ) + ' ' + str( numorders ) + ' ' + str( nummoments ) + ' > ' + filename )
		#filename = 'data/infinite3D_isotropicplane_isotropicscatter_gamma2U_tracklength_c' + str(c) + '_mut' + str(mut) + '.txt'
		#print 'computing: ' + filename
		#os.system( './infinite3D_isotropicplane_isotropicscatter_gamma2U_tracklength ' + str(c) + ' ' + str(mut) + ' ' + str(maxr) + ' ' + str(dr) + ' ' + str(du) + ' ' + str( numsamples ) + ' ' + str( numorders ) + ' ' + str( nummoments ) + ' > ' + filename )
		filename = 'data/infinite3D_isotropicplane_HG_exponential_c' + str(c) + '_mfp' + str(mfp) + '_g' + str(g) + '.txt'
		print 'computing: ' + filename
		os.system( './infinite3D_isotropicplane_HG_exponential ' + str(c) + ' ' + str(mfp) + ' ' + str(maxr) + ' ' + str(dr) + ' ' + str(du) + ' ' + str( numsamples ) + ' ' + str( numorders ) + ' ' + str( nummoments ) + ' ' + str(g) + ' > ' + filename )

