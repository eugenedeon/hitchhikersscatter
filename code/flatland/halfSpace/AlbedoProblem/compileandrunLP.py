# Compile and run simulations: 

# Half space albedo problem isotropic scattering vacuum boundary

import sys
import os
import math

if not os.path.exists('data'):
	os.mkdir('data')


os.system('g++ HalfspaceAlbedoProblemDeltaLP.cpp -o HalfspaceAlbedoProblemDeltaLP -O3 -I ../../../include/')

du = 0.01
maxz = 25.0
dz = 0.1
numorders = 50
nummoments = 5

cs = [ 0.01, 0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99, 0.999 ]

# ISOTROPIC SCATTERING

# Delta illumination
for mui in [0.25, 0.5, 1.0]:
	for c in [0.5, 0.7, 0.9, 0.99, 0.999]:
		pa = 0.5
		sigratio = 6.5
		lambsum = 6.0
		expected_num_samples_inv = math.sqrt( ( 1.0 - c ) / ( c * c ) )
		numsamples = min( 10000000, int( 10000000.0 * expected_num_samples_inv ) )
		print 'num samples: ' + str(numsamples)
		filename = 'data/HalfspaceAlbedoProblemDeltaLP_c' + str(c) + '_mui' + str(mui) + '_pa' + str(pa) + '_sigratio' + str(sigratio) + '_lambsum' + str(lambsum) + '.txt'
		print 'computing: ' + filename
		os.system( './HalfspaceAlbedoProblemDeltaLP ' + str(c) + ' ' + str(pa) + ' ' + str(mui)+ ' ' + str(du) + ' ' + str(dz) + ' ' + str(maxz) + ' ' + str( numsamples ) + ' ' + str( numorders ) + ' ' + str( nummoments ) + ' ' + str(sigratio) + ' ' + str(lambsum) + ' > ' + filename )
