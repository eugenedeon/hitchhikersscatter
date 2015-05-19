import sys
import os

# internal distributions
os.system('g++ halfRod_albedoproblem_anisotropicscatter_exponential.cpp -o halfRod_albedoproblem_anisotropicscatter_exponential -O3 -I ../../../include/')

# variation in g
numsamples = 1000000
c = 0.95
mut = 1
maxx = 10
dx = 1
numorders = 2
nummoments = 2
if 1:
  g = -0.9
  while g <= 0.9:    
    filename = 'data/halfrod_albedoproblem_anisotropicscatter_exp_c' + str(c) + '_mut' + str(mut) + '_g' + str(g) + '.txt'
    print 'computing: ' + filename
    os.system( './halfRod_albedoproblem_anisotropicscatter_exponential ' + str(c) + ' ' + str(mut) + ' ' + str(maxx) + ' ' + str(dx) + ' ' + str( numsamples ) + ' ' + str( numorders ) + ' ' + str( nummoments ) + ' ' + str( g ) + ' > ' + filename )
    g += 0.1

numsamples = 10000000
if 1:
  dx = 0.1
  g = 0.7
  c = 0.7
  numorders = 10
  nummoments = 10
  filename = 'data/halfrod_albedoproblem_anisotropicscatter_exp_c' + str(c) + '_mut' + str(mut) + '_g' + str(g) + '.txt'
  print 'computing: ' + filename
  os.system( './halfRod_albedoproblem_anisotropicscatter_exponential ' + str(c) + ' ' + str(mut) + ' ' + str(maxx) + ' ' + str(dx) + ' ' + str( numsamples ) + ' ' + str( numorders ) + ' ' + str( nummoments ) + ' ' + str( g ) + ' > ' + filename )
    
if 1:
  dx = 0.1
  g = 0.7
  c = 0.7
  mut = 0.4
  numorders = 10
  nummoments = 10
  filename = 'data/halfrod_albedoproblem_anisotropicscatter_exp_c' + str(c) + '_mut' + str(mut) + '_g' + str(g) + '.txt'
  print 'computing: ' + filename
  os.system( './halfRod_albedoproblem_anisotropicscatter_exponential ' + str(c) + ' ' + str(mut) + ' ' + str(maxx) + ' ' + str(dx) + ' ' + str( numsamples ) + ' ' + str( numorders ) + ' ' + str( nummoments ) + ' ' + str( g ) + ' > ' + filename )
    