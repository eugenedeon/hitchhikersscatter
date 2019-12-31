# Compile and run simulations: 

# 3D INFINITE MEDIUM ISOTROPIC POINT SOURCE

import sys
import os
import math

if not os.path.exists('MCdata'):
	os.mkdir('MCdata')

# classical - various phase functions
os.system('g++ infinite3D_isotropicpoint_isotropicscatter_exponential.cpp -o infinite3D_isotropicpoint_isotropicscatter_exponential -O3 -I ../../../include/')
os.system('g++ infinite3D_isotropicpoint_rayleighscatter_exponential.cpp -o infinite3D_isotropicpoint_rayleighscatter_exponential -O3 -I ../../../include/')
os.system('g++ infinite3D_isotropicpoint_LSscatter_exponential.cpp -o infinite3D_isotropicpoint_LSscatter_exponential -O3 -I ../../../include/')
os.system('g++ infinite3D_isotropicpoint_linanisoscatter_exponential.cpp -o infinite3D_isotropicpoint_linanisoscatter_exponential -O3 -I ../../../include/')
os.system('g++ infinite3D_isotropicpoint_HG_exponential.cpp -o infinite3D_isotropicpoint_HG_exponential -O3 -I ../../../include/')

cs = [ 0.01, 0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99, 0.999 ]
numorders = 25
nummoments = 10

# ISOTROPIC, RAYLEIGH & LAMBERT SPHERE SCATTERING

for mfp in [1.0, 0.3]: 		# mean-free path
	for c in cs: 			# single-scattering albedo
		nu_0 = 1/math.sqrt(1 - math.pow(c,2.4429445001914587 + 0.5786368322364553/c - 0.021581332427913873*c))
		expected_num_samples_inv = 1.0 - c
		numsamples = min( 10000000, int( 10000000.0 * expected_num_samples_inv ) )
		print 'num samples: ' + str(numsamples)
		maxr = nu_0 * 25.0 * mfp
		dr = maxr / 500.0
		du = 2.0 / 30.0
		
		filename = 'MCdata/inf3d_isotropicpoint_isotropicscatter_c' + str(c) + '_mfp' + str(mfp) + '.txt'
		print 'computing: ' + filename
		os.system( './infinite3d_isotropicpoint_isotropicscatter_exponential ' + str(c) + ' ' + str(mfp) + ' ' + str(maxr) + ' ' + str(dr) + ' ' + str(du) + ' ' + str( numsamples ) + ' ' + str( numorders ) + ' ' + str( nummoments ) + ' > ' + filename )

		filename = 'MCdata/inf3d_isotropicpoint_rayleighscatter_c' + str(c) + '_mfp' + str(mfp) + '.txt'
		print 'computing: ' + filename
		os.system( './infinite3d_isotropicpoint_rayleighscatter_exponential ' + str(c) + ' ' + str(mfp) + ' ' + str(maxr) + ' ' + str(dr) + ' ' + str(du) + ' ' + str( numsamples ) + ' ' + str( numorders ) + ' ' + str( nummoments ) + ' > ' + filename )

		filename = 'MCdata/inf3d_isotropicpoint_LSscatter_c' + str(c) + '_mfp' + str(mfp) + '.txt'
		print 'computing: ' + filename
		os.system( './infinite3d_isotropicpoint_LSscatter_exponential ' + str(c) + ' ' + str(mfp) + ' ' + str(maxr) + ' ' + str(dr) + ' ' + str(du) + ' ' + str( numsamples ) + ' ' + str( numorders ) + ' ' + str( nummoments ) + ' > ' + filename )

# LINEAR ANISOTROPIC SCATTERING

for mfp in [1.0, 0.3]: 			# mean-free path
	for c in cs: 				# single-scattering albedo
		for b in [-0.9,0.7]:	# anisotropy parameter b

			g = b / 3.0
			nu_0 = 1/math.sqrt(1 - math.pow(c,2.4429445001914587 + 0.5786368322364553/c - 0.021581332427913873*c))
			expected_num_samples_inv = 1.0 - c
			numsamples = min( 10000000, int( 10000000.0 * expected_num_samples_inv ) )
			print 'num samples: ' + str(numsamples)
			maxr = nu_0 * ( 1.0 - g ) * 25.0 * mfp
			dr = maxr / 500.0
			du = 2.0 / 30.0
			numorders = 50
			nummoments = 5

			filename = 'MCdata/inf3d_isotropicpoint_linanisoscatter_c' + str(c) + '_mfp' + str(mfp) + '_b' + str(b) + '.txt'
			print 'computing: ' + filename
			os.system( './infinite3D_isotropicpoint_linanisoscatter_exponential ' + str(c) + ' ' + str(mfp) + ' ' + str(maxr) + ' ' + str(dr) + ' ' + str(du) + ' ' + str( numsamples ) + ' ' + str( numorders ) + ' ' + str( nummoments ) + ' ' + str(b) + ' > ' + filename )

# HENYEY GREENSTEIN (HG)

for mfp in [1.0, 0.3]: 					# mean-free path
	for c in cs: 						# single-scattering albedo
		for g in [-0.5, 0.3, 0.5, 0.7, 0.8, 0.9]:	# mean cosine g

			nu_0 = 1/math.sqrt(1 - math.pow(c,2.4429445001914587 + 0.5786368322364553/c - 0.021581332427913873*c))
			expected_num_samples_inv = 1.0 - c
			numsamples = min( 10000000, int( 10000000.0 * expected_num_samples_inv ) )
			print 'num samples: ' + str(numsamples)
			maxr = nu_0 * ( 1.0 - g ) * 25.0 * mfp
			dr = maxr / 500.0
			du = 2.0 / 30.0
			numorders = 50
			nummoments = 5

			filename = 'MCdata/inf3d_isotropicpoint_HG_c' + str(c) + '_mfp' + str(mfp) + '_g' + str(g) + '.txt'
			print 'computing: ' + filename
			os.system( './infinite3D_isotropicpoint_HG_exponential ' + str(c) + ' ' + str(mfp) + ' ' + str(maxr) + ' ' + str(dr) + ' ' + str(du) + ' ' + str( numsamples ) + ' ' + str( numorders ) + ' ' + str( nummoments ) + ' ' + str(g) + ' > ' + filename )

			