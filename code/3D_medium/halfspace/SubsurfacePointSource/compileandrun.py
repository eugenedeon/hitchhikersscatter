# Compile and run simulations: 

# Half space subsurface point source isotropic scattering vacuum boundary

import sys
import os
import math

if not os.path.exists('data'):
	os.mkdir('data')

# classical
os.system('g++ HalfspaceSubsurfacePointExponential.cpp -o HalfspaceSubsurfacePointExponential -O3 -I ../../../include/')

dr = 0.1
maxr = 10.0
maxz = 25.0
dz = 0.1
numorders = 50
nummoments = 5

cs = [ 0.01, 0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99, 0.999 ]

# ISOTROPIC SCATTERING

# exponential
for depth in [1.0, 0.5, 5.0]:
	for c in [ 0.5, 0.9 ]:
		expected_num_samples_inv = math.sqrt( ( 1.0 - c ) / ( c * c ) )
		numsamples = min( 10000000, int( 10000000.0 * expected_num_samples_inv ) )
		print 'num samples: ' + str(numsamples)
		mut = 1
		filename = 'data/subsurfacePoint_exp_c' + str(c) + '_mut' + str(mut) + '_depth' + str(depth) + '.txt'
		print 'computing: ' + filename
		os.system( './HalfspaceSubsurfacePointExponential ' + str(c) + ' ' + str(mut) + ' ' + str(depth)+ ' ' + str(dr) + ' ' + str(maxr) + ' ' + str(dz) + ' ' + str(maxz) + ' ' + str( numsamples ) + ' > ' + filename )

