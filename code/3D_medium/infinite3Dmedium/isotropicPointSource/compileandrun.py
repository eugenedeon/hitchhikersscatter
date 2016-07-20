# Compile and run simulations: 

# 3D INFINITE MEDIUM ISOTROPIC POINT SOURCE

import sys
import os
import math

if not os.path.exists('data'):
	os.mkdir('data')

# classical - various phase functions
os.system('g++ infinite3D_isotropicpoint_isotropicscatter_exponential.cpp -o infinite3D_isotropicpoint_isotropicscatter_exponential -O3 -I ../../../include/')
os.system('g++ infinite3D_isotropicpoint_linanisoscatter_exponential.cpp -o infinite3D_isotropicpoint_linanisoscatter_exponential -O3 -I ../../../include/')

cs = [ 0.01, 0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99, 0.999 ]

# ISOTROPIC SCATTERING

for mut in [1.0,3.0]:
	for c in cs:
		expected_num_samples_inv = math.sqrt( ( 1.0 - c ) / ( c * c ) )
		numsamples = min( 10000000, int( 10000000.0 * expected_num_samples_inv ) )
		print 'num samples: ' + str(numsamples)
		maxr = 25.0 / mut
		dr = 0.2 / mut
		du = 2.0 / 30.0
		numorders = 50
		nummoments = 5
		filename = 'data/inf3d_isotropicpoint_isotropicscatter_c' + str(c) + '_mut' + str(mut) + '.txt'
		print 'computing: ' + filename
		os.system( './infinite3d_isotropicpoint_isotropicscatter_exponential ' + str(c) + ' ' + str(mut) + ' ' + str(maxr) + ' ' + str(dr) + ' ' + str(du) + ' ' + str( numsamples ) + ' ' + str( numorders ) + ' ' + str( nummoments ) + ' > ' + filename )

# LINEAR ANISOTROPIC SCATTERING

for mut in [1.0,3.0]:
	for c in cs:
		for b in [-0.9,0.7]:
			expected_num_samples_inv = math.sqrt( ( 1.0 - c ) / ( c * c ) )
			numsamples = min( 10000000, int( 10000000.0 * expected_num_samples_inv ) )
			print 'num samples: ' + str(numsamples)
			maxr = 25.0 / mut
			dr = 0.2 / mut
			du = 2.0 / 30.0
			numorders = 10
			nummoments = 8
			filename = 'data/inf3d_isotropicpoint_linanisoscatter_c' + str(c) + '_mut' + str(mut) + '_b' + str(b) + '.txt'
			print 'computing: ' + filename
			os.system( './infinite3D_isotropicpoint_linanisoscatter_exponential ' + str(c) + ' ' + str(mut) + ' ' + str(maxr) + ' ' + str(dr) + ' ' + str(du) + ' ' + str( numsamples ) + ' ' + str( numorders ) + ' ' + str( nummoments ) + ' ' + str(b) + ' > ' + filename )

			