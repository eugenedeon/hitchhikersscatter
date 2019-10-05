# Compile and run simulations: 

# 4D albedo problem isotropic scattering vacuum boundary

import sys
import os
import math

if not os.path.exists('data'):
	os.mkdir('data')

# classical - various phase functions
os.system('g++ 4DAlbedoProblemWhiteSky.cpp -o 4DAlbedoProblemWhiteSky -O3 -I ../../../include/')

du = 0.01
maxz = 25.0
dz = 0.1
numorders = 50
nummoments = 5

cs = [ 0.01, 0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99, 0.999 ]

# ISOTROPIC SCATTERING

# whitesky illumination
for c in cs:
	expected_num_samples_inv = math.sqrt( ( 1.0 - c ) / ( c * c ) )
	numsamples = min( 10000000, int( 10000000.0 * expected_num_samples_inv ) )
	print 'num samples: ' + str(numsamples)
	mut = 1
	filename = 'data/albedoproblem_whitesky_c' + str(c) + '_mut' + str(mut)  + '.txt'
	print 'computing: ' + filename
	os.system( './4DAlbedoProblemWhiteSky ' + str(c) + ' ' + str(mut) + ' ' + str(du) + ' ' + str(dz) + ' ' + str(maxz) + ' ' + str( numsamples ) + ' ' + str( numorders ) + ' ' + str( nummoments ) + ' > ' + filename )
