import sys
import os

if not os.path.exists('data'):
  os.mkdir('data')

# internal distributions
# os.system('g++ halfRod_pointsource_isotropicscatter_exponential.cpp -o halfRod_pointsource_isotropicscatter_exponential -O3 -I ../../../include/')
os.system('g++ halfRod_pointsource_isotropicscatter_K2.cpp -o halfRod_pointsource_isotropicscatter_K2 -O3 -I ../../../include/')

# fluence parameters
maxx = 10.0
dx = 20.0 / 101

cs = [ 0.1, 0.3, 0.5, 0.7, 0.9 ]
x0s = [0.2, 1.0, 3.0]

numsamples = 10000000
numorders = 20
nummoments = 10

if 1:
  for x0 in x0s:
      for c in cs:
          filename = 'data/halfrod_pointsource_isotropicscatter_K2_c' + str(c) + '_x0_' + str(x0) + '.txt'
          print 'computing: ' + filename
          os.system( './halfRod_pointsource_isotropicscatter_K2 ' + str(c) + ' ' + str(x0) + ' ' + str(maxx) + ' ' + str(dx) + ' ' + str( numsamples ) + ' ' + str( numorders ) + ' ' + str( nummoments ) + ' > ' + filename )

numsamples = 1000000
if 1:
  for x0 in x0s:
      for c in [0.95,0.98,0.99,0.999]:
          filename = 'data/halfrod_pointsource_isotropicscatter_K2_c' + str(c) + 'x0_' + str(x0) + '.txt'
          print 'computing: ' + filename
          os.system( './halfRod_pointsource_isotropicscatter_K2 ' + str(c) + ' ' + str(x0) + ' ' + str(maxx) + ' ' + str(dx) + ' ' + str( numsamples ) + ' ' + str( numorders ) + ' ' + str( nummoments ) + ' > ' + filename )
