import sys
import os

os.mkdir('data')

# internal distributions
os.system('g++ halfRod_albedoproblem_isotropicscatter_exponential.cpp -o halfRod_albedoproblem_isotropicscatter_exponential -O3 -I ../../../include/')

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
          filename = 'data/halfrod_albedoproblem_isotropicscatter_exp_c' + str(c) + '_mut' + str(mut) + '.txt'
          print 'computing: ' + filename
          os.system( './halfRod_albedoproblem_isotropicscatter_exponential ' + str(c) + ' ' + str(mut) + ' ' + str(maxx) + ' ' + str(dx) + ' ' + str( numsamples ) + ' ' + str( numorders ) + ' ' + str( nummoments ) + ' > ' + filename )

numsamples = 1000000
if 1:
  for mut in muts:
      for c in [0.95,0.98,0.99,0.999]:
          filename = 'data/halfrod_albedoproblem_isotropicscatter_exp_c' + str(c) + '_mut' + str(mut) + '.txt'
          print 'computing: ' + filename
          os.system( './halfRod_albedoproblem_isotropicscatter_exponential ' + str(c) + ' ' + str(mut) + ' ' + str(maxx) + ' ' + str(dx) + ' ' + str( numsamples ) + ' ' + str( numorders ) + ' ' + str( nummoments ) + ' > ' + filename )
