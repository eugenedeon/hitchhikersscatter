# Compile and run simulations: 

# Half space albedo problem isotropic scattering vacuum boundary

import sys
import os
import math

if not os.path.exists('data'):
	os.mkdir('data')

# compile
os.system('g++ HalfspaceAlbedoProblemEscapeProbability.cpp -o HalfspaceAlbedoProblemEscapeProbability -O3 -I ../../../include/')


numz = 200
dz = 0.1

cs = [ 0.01, 0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99, 0.999 ]

# ISOTROPIC SCATTERING

# escape probability
for c in cs:
	expected_num_samples_inv = math.sqrt( ( 1.0 - c ) / ( c * c ) )
	numsamples = min( 10000000, int( 10000000.0 * expected_num_samples_inv ) ) / 20
	print 'num samples: ' + str(numsamples)
	filename = 'data/albedoproblem_escape_c' + str(c) + '.txt'
	print 'computing: ' + filename
	os.system( './HalfspaceAlbedoProblemEscapeProbability ' + str(c) + ' ' + str(dz) + ' ' + str(numz) + ' ' + str( numsamples ) + ' > ' + filename )
