import sys
import os

# internal distributions
os.system('g++ infiniteRod_isotropicpoint_isotropicscatter_exponential.cpp -o infiniteRod_isotropicpoint_isotropicscatter_exponential -O3 -I ../../../include/')
os.system('g++ infiniteRod_isotropicpoint_anisotropicscatter_exponential.cpp -o infiniteRod_isotropicpoint_anisotropicscatter_exponential -O3 -I ../../../include/')

# fluence parameters
maxx = 10.0
dx = 20.0 / 101

cs = [ 0.1, 0.3, 0.5, 0.7, 0.9 ]
muts = [1.0, 3.0]

numsamples = 10000000
numorders = 20
nummoments = 10

if 1:
  for mut in muts:
      for c in cs:
          filename = 'data/infrod_isotropicpoint_isotropicscatter_exp_c' + str(c) + '_mut' + str(mut) + '.txt'
          print 'computing: ' + filename
          os.system( './infiniteRod_isotropicpoint_isotropicscatter_exponential ' + str(c) + ' ' + str(mut) + ' ' + str(maxx) + ' ' + str(dx) + ' ' + str( numsamples ) + ' ' + str( numorders ) + ' ' + str( nummoments ) + ' > ' + filename )

if 1:
	for g in [-0.9, -0.3, 0.3, 0.9]:
	  for mut in muts:
	      for c in cs:
	          filename = 'data/infrod_isotropicpoint_anisotropicscatter_exp_c' + str(c) + '_mut' + str(mut) + '_g' + str(g) + '.txt'
	          print 'computing: ' + filename
	          os.system( './infiniteRod_isotropicpoint_anisotropicscatter_exponential ' + str(c) + ' ' + str(mut) + ' ' + str(maxx) + ' ' + str(dx) + ' ' + str( numsamples ) + ' ' + str( numorders ) + ' ' + str( nummoments ) + ' ' + str( g ) + ' > ' + filename )
