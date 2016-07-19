import sys
import os

if not os.path.exists('data'):
  os.mkdir('data')

# internal distributions
os.system('g++ halfRod_albedoproblem_anisotropicscatter_exponential.cpp -o halfRod_albedoproblem_anisotropicscatter_exponential -O3 -I ../../../include/')

maxx = 10.0
dx = 20.0 / 101

cs = [ 0.1, 0.3, 0.5, 0.7, 0.9 ]
muts = [1.0, 3.0]
gs = [-0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9]

numsamples = 10000000
numorders = 20
nummoments = 10

if 1:
  for g in gs:
    for c in cs:
      for mut in muts:
        filename = 'data/halfrod_albedoproblem_anisotropicscatter_exp_c' + str(c) + '_mut' + str(mut) + '_g' + str(g) + '.txt'
        print 'computing: ' + filename
        os.system( './halfRod_albedoproblem_anisotropicscatter_exponential ' + str(c) + ' ' + str(mut) + ' ' + str(maxx) + ' ' + str(dx) + ' ' + str( numsamples ) + ' ' + str( numorders ) + ' ' + str( nummoments ) + ' ' + str( g ) + ' > ' + filename )

numsamples = 1000000
if 1:
  for g in gs:
    for c in [0.95,0.98,0.99,0.999]:
      for mut in muts:
        filename = 'data/halfrod_albedoproblem_anisotropicscatter_exp_c' + str(c) + '_mut' + str(mut) + '_g' + str(g) + '.txt'
        print 'computing: ' + filename
        os.system( './halfRod_albedoproblem_anisotropicscatter_exponential ' + str(c) + ' ' + str(mut) + ' ' + str(maxx) + ' ' + str(dx) + ' ' + str( numsamples ) + ' ' + str( numorders ) + ' ' + str( nummoments ) + ' ' + str( g ) + ' > ' + filename )