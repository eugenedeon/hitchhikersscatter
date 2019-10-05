# Compile and run simulations: 

# 3D INFINITE MEDIUM ISOTROPIC POINT SOURCE - power law scattering

import sys
import os
import math

if not os.path.exists('data3'):
	os.mkdir('data3')

# power law uncorrelated
os.system('g++ infinite3D_isotropicpoint_isotropicscatter_powerLawU.cpp -o infinite3D_isotropicpoint_isotropicscatter_powerLawU -O3 -I ../../../include/')
os.system('g++ infinite3D_isotropicpoint_isotropicscatter_powerLawC.cpp -o infinite3D_isotropicpoint_isotropicscatter_powerLawC -O3 -I ../../../include/')

cs = [ 0.5, 0.7, 0.9, 0.99, 0.999 ]

# ISOTROPIC SCATTERING

mfp = 1.0
#for a in [1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 10.0, 20.0, 50.0, 100.0]:
for a in [20.0, 50.0, 100.0]:
	for c in cs:
		expected_num_samples_inv = math.sqrt( ( 1.0 - c ) / ( c * c ) )
		numsamples = 4 * min( 10000000, int( 10000000.0 * expected_num_samples_inv ) )
		print 'num samples: ' + str(numsamples)
		maxr = 50.0 * mfp * math.sqrt( 1 / (1 - c ) ) / 1.5
		dr = maxr / 400.0
		du = 2.0 / 30.0
		numorders = 2
		nummoments = 3
		filename = 'data3/inf3d_isotropicpoint_isotropicscatter_powerlawC' + str(a) + '_c' + str(c) + '_mfp' + str(mfp) + '.txt'
		print 'computing: ' + filename
		os.system( './infinite3D_isotropicpoint_isotropicscatter_powerLawC ' + str(c) + ' ' + str(mfp) + ' ' + str(maxr) + ' ' + str(dr) + ' ' + str(du) + ' ' + str( numsamples ) + ' ' + str( numorders ) + ' ' + str( nummoments ) + ' ' + str(a) + ' > ' + filename )
		filename = 'data3/inf3d_isotropicpoint_isotropicscatter_powerlawU' + str(a) + '_c' + str(c) + '_mfp' + str(mfp) + '.txt'
		print 'computing: ' + filename
		os.system( './infinite3D_isotropicpoint_isotropicscatter_powerLawU ' + str(c) + ' ' + str(mfp) + ' ' + str(maxr) + ' ' + str(dr) + ' ' + str(du) + ' ' + str( numsamples ) + ' ' + str( numorders ) + ' ' + str( nummoments ) + ' ' + str(a) + ' > ' + filename )

