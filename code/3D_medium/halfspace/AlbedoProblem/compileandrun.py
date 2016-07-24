# Compile and run simulations: 

# Half space albedo problem isotropic scattering vacuum boundary

import sys
import os
import math

if not os.path.exists('data'):
	os.mkdir('data')

# classical - various phase functions
os.system('g++ HalfspaceAlbedoProblemDelta.cpp -o HalfspaceAlbedoProblemDelta -O3 -I ../../../include/')
os.system('g++ HalfspaceAlbedoProblemWhiteSky.cpp -o HalfspaceAlbedoProblemWhiteSky -O3 -I ../../../include/')

du = 0.01
maxz = 25.0
dz = 0.1
numorders = 50
nummoments = 5

cs = [ 0.01, 0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99, 0.999 ]

# ISOTROPIC SCATTERING

# Delta illumination
for mui in [1.0, 0.5, 0.25, 0.1]:
	for c in cs:
		expected_num_samples_inv = math.sqrt( ( 1.0 - c ) / ( c * c ) )
		numsamples = min( 10000000, int( 10000000.0 * expected_num_samples_inv ) )
		print 'num samples: ' + str(numsamples)
		mut = 1
		filename = 'data/albedoproblem_delta_c' + str(c) + '_mut' + str(mut) + '_mui' + str(mui) + '.txt'
		print 'computing: ' + filename
		os.system( './HalfspaceAlbedoProblemDelta ' + str(c) + ' ' + str(mut) + ' ' + str(mui)+ ' ' + str(du) + ' ' + str(dz) + ' ' + str(maxz) + ' ' + str( numsamples ) + ' ' + str( numorders ) + ' ' + str( nummoments ) + ' > ' + filename )

# Whitesky illumination
for c in cs:
	expected_num_samples_inv = math.sqrt( ( 1.0 - c ) / ( c * c ) )
	numsamples = min( 10000000, int( 10000000.0 * expected_num_samples_inv ) )
	print 'num samples: ' + str(numsamples)
	mut = 1
	filename = 'data/albedoproblem_whitesky_c' + str(c) + '_mut' + str(mut) + '.txt'
	print 'computing: ' + filename
	os.system( './HalfspaceAlbedoProblemWhiteSky ' + str(c) + ' ' + str(mut) + ' ' + str(du) + ' ' + str(dz) + ' ' + str(maxz) + ' ' + str( numsamples ) + ' ' + str( numorders ) + ' ' + str( nummoments ) + ' > ' + filename )
