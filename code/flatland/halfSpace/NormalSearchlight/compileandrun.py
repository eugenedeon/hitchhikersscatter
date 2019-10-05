# Compile and run simulations: 

# FLATLAND HALF SPACE SEARCHLIGHT

import sys
import os
import math

# classical - various phase functions
os.system('g++ halfSpaceFlatland_normalSearchlight_isotropicscatter_exponential.cpp -o halfSpaceFlatland_normalSearchlight_isotropicscatter_exponential -O3 -I ../../../include/')

cs = [ 0.01, 0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99, 0.999 ]

# ISOTROPIC SCATTERING

for mut in [1.0,3.0]:
	for c in cs:
		expected_num_samples_inv = math.sqrt( ( 1.0 - c ) / ( c * c ) )
		numsamples = min( 50000000, int( 50000000.0 * expected_num_samples_inv ) )
		print 'num samples: ' + str(numsamples)
		maxr = 25.0 / mut
		dr = 0.2 / mut
		du = 1.0 / 30.0
		filename = 'data/halfSpaceFlatland_normalSearchlight_isotropicscatter_exponential_c' + str(c) + '_mut' + str(mut) + '.txt'
		print 'computing: ' + filename
		os.system( './halfSpaceFlatland_normalSearchlight_isotropicscatter_exponential ' + str(c) + ' ' + str(mut) + ' ' + str(maxr) + ' ' + str(dr) + ' ' + str(du) + ' ' + str( numsamples ) + ' > ' + filename )