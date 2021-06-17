# hitchhikersscatter

This repo contains C++ and Mathematica code for Monte Carlo and exact/approximate analytical solutions to multiple scattering problems as described in the book:

A Hitchhiker's guide to Multiple Scattering
(c) 2016 Eugene d'Eon, www.eugenedeon.com/hitchhikers

This is a collection of benchmark solutions for non-trivial multiple scattering processes under the umbrella of scalar radiative transfer (and some random-media non-Beerian generalizations).  This codebase includes code for (not limited to):
* Phase functions
* NDFs for rough surfaces
* Green's functions for an isotropic point source in infinite dD spaces under a variety of scattering and absorption conditions
* The albedo problem in a 1D rod, 2D Flatland, 3D and higher dimensional spaces
* Diffusion approximations
* H-functions and solutions to Fredholm integral equations

## Monte Carlo code
Each scattering benchmark problem has a python script that compiles and generates MC data to be loaded in the corresponding Mathematica notebook.  
The submitted code assumes Mac OS g++ with drand48().  It should be easily portable to other environments by changing the RandomReal() function in random.h.

## Deterministic code
The deterministic exact and approximate solutions of various transport theory problems are provided as Mathematica notebooks that load the corresponding Monte Carlo data to generate plots that ensure the two are in agreement.  Each Mathematica notebook is also saved as a .pdf for convenience (these may get out of sync now and then).
