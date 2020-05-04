# Compile and run simulations: 

# 3D INFINITE MEDIUM ISOTROPIC POINT SOURCE

import sys
import os
import math

if not os.path.exists('MCdata'):
	os.mkdir('MCdata')

# besselK flights - Isotropic and Rayleigh scattering
os.system('g++ infinite3D_isotropicpoint_isotropicscatter_besselK.cpp -o infinite3D_isotropicpoint_isotropicscatter_besselK -O3 -I ../../../include/')
os.system('g++ infinite3D_isotropicpoint_rayleighscatter_besselK.cpp -o infinite3D_isotropicpoint_rayleighscatter_besselK -O3 -I ../../../include/')

cs = [ 0.01, 0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99, 0.999 ]
numorders = 25
nummoments = 10

for c in cs: 			# single-scattering albedo
    nu_0 = 1/math.sqrt(1 - math.pow(c,2.4429445001914587 + 0.5786368322364553/c - 0.021581332427913873*c))
    expected_num_samples_inv = 1.0 - c
    numsamples = min( 10000000, int( 10000000.0 * expected_num_samples_inv ) )
    print 'num samples: ' + str(numsamples)
    maxr = nu_0 * 25.0
    dr = maxr / 500.0
    du = 2.0 / 30.0
    maxt = maxr
    dt = 0.1

    filename = 'MCdata/inf3d_isotropicpoint_isotropicscatter_besselK_c' + str(c) + '.txt'
    print 'computing: ' + filename
    os.system( './infinite3D_isotropicpoint_isotropicscatter_besselK ' + str(c) + ' ' + str(maxr) + ' ' + str(dr) + ' ' + str(maxt) + ' ' + str(dt) + ' ' + str(du) + ' ' + str( numsamples ) + ' ' + str( numorders ) + ' ' + str( nummoments ) + ' > ' + filename )

    filename = 'MCdata/inf3d_isotropicpoint_rayleighscatter_besselK_c' + str(c) + '.txt'
    print 'computing: ' + filename
    os.system( './infinite3d_isotropicpoint_rayleighscatter_besselK ' + str(c) + ' ' + str(maxr) + ' ' + str(dr) + ' ' + str(maxt) + ' ' + str(dt) + ' ' + str(du) + ' ' + str( numsamples ) + ' ' + str( numorders ) + ' ' + str( nummoments ) + ' > ' + filename )
