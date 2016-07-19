import sys
import os

os.mkdir('data')

# compile MC simulator
os.system('g++ finiteRod_albedoproblem_isotropicscatter_exponential.cpp -o finiteRod_albedoproblem_isotropicscatter_exponential -O3 -I ../../../include/')
os.system('g++ finiteRod_albedoproblem_anisotropicscatter_exponential.cpp -o finiteRod_albedoproblem_anisotropicscatter_exponential -O3 -I ../../../include/')

# R/T length variation
length = 0.1
g = 0.7
if 1:
  while length < 100:

    numsamples = 10000000
    numorders = 2
    nummoments = 2
    dx = 20.0 / 101
    mut = 1
    c = 0.7

    filename = 'data/finiteRod_albedoproblem_isotropicscatter_exp_c' + str(c) + '_mut' + str(mut) + '_length' + str(length) + '.txt'
    print 'computing: ' + filename
    os.system( './finiteRod_albedoproblem_isotropicscatter_exponential ' + str(c) + ' ' + str(mut) + ' ' + str(length) + ' ' + str(dx) + ' ' + str( numsamples ) + ' ' + str( numorders ) + ' ' + str( nummoments ) + ' > ' + filename )

    filename = 'data/finiteRod_albedoproblem_anisotropicscatter_exp_c' + str(c) + '_mut' + str(mut) + '_length' + str(length) + '_g' + str(g) + '.txt'
    print 'computing: ' + filename
    os.system( './finiteRod_albedoproblem_anisotropicscatter_exponential ' + str(c) + ' ' + str(mut) + ' ' + str(length) + ' ' + str(dx) + ' ' + str( numsamples ) + ' ' + str( numorders ) + ' ' + str( nummoments ) + ' ' + str(g) + ' > ' + filename )

    length *= 1.5

# R/T albedo variation
length = 2.6
c = 0.999

if 1:
  while c > 0.2:

    numsamples = 10000000
    numorders = 2
    nummoments = 2
    dx = length / 101
    mut = 1

    g = 0.7

    filename = 'data/finiteRod_albedoproblem_isotropicscatter_exp_c' + str(c) + '_mut' + str(mut) + '_length' + str(length) + '.txt'
    print 'computing: ' + filename
    os.system( './finiteRod_albedoproblem_isotropicscatter_exponential ' + str(c) + ' ' + str(mut) + ' ' + str(length) + ' ' + str(dx) + ' ' + str( numsamples ) + ' ' + str( numorders ) + ' ' + str( nummoments ) + ' > ' + filename )

    filename = 'data/finiteRod_albedoproblem_anisotropicscatter_exp_c' + str(c) + '_mut' + str(mut) + '_length' + str(length) + '_g' + str(g) + '.txt'
    print 'computing: ' + filename
    os.system( './finiteRod_albedoproblem_anisotropicscatter_exponential ' + str(c) + ' ' + str(mut) + ' ' + str(length) + ' ' + str(dx) + ' ' + str( numsamples ) + ' ' + str( numorders ) + ' ' + str( nummoments ) + ' ' + str(g) + ' > ' + filename )

    c *= 0.9

# R/T g variation

g = -0.9

if 1:
  while g <= 0.91:

    length = 1.2
    c = 0.8

    numsamples = 10000000
    numorders = 2
    nummoments = 2
    dx = length / 101
    mut = 1

    filename = 'data/finiteRod_albedoproblem_anisotropicscatter_exp_c' + str(c) + '_mut' + str(mut) + '_length' + str(length) + '_g' + str(g) + '.txt'
    print 'computing: ' + filename
    os.system( './finiteRod_albedoproblem_anisotropicscatter_exponential ' + str(c) + ' ' + str(mut) + ' ' + str(length) + ' ' + str(dx) + ' ' + str( numsamples ) + ' ' + str( numorders ) + ' ' + str( nummoments ) + ' ' + str(g) + ' > ' + filename )

    g += 0.1

numsamples = 10000000
numorders = 10
nummoments = 10

for c in [0.3, 0.7, 0.9]:

    for length in [0.3, 1.0, 3.0]:
        dx = length / 101
        for mut in [1.0, 3.0]:

            filename = 'data/finiteRod_albedoproblem_isotropicscatter_exp_c' + str(c) + '_mut' + str(mut) + '_length' + str(length) + '.txt'
            print 'computing: ' + filename
            os.system( './finiteRod_albedoproblem_isotropicscatter_exponential ' + str(c) + ' ' + str(mut) + ' ' + str(length) + ' ' + str(dx) + ' ' + str( numsamples ) + ' ' + str( numorders ) + ' ' + str( nummoments ) + ' > ' + filename )

            dx = length / 51

            for g in [-0.5, 0.3, 0.7]:

                filename = 'data/batch_finiteRod_albedoproblem_anisotropicscatter_exp_c' + str(c) + '_mut' + str(mut) + '_length' + str(length) + '_g' + str(g) + '.txt'
                print 'computing: ' + filename
                os.system( './finiteRod_albedoproblem_anisotropicscatter_exponential ' + str(c) + ' ' + str(mut) + ' ' + str(length) + ' ' + str(dx) + ' ' + str( numsamples ) + ' ' + str( numorders ) + ' ' + str( nummoments ) + ' ' + str(g) + ' > ' + filename )

