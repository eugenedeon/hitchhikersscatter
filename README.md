# hitchhikersscatter

This repo contains C++ and Mathematica code for Monte Carlo and exact/approximate analytical solutions to multiple scattering problems as described in the book:

A Hitchhiker's guide to Multiple Scattering
(c) 2016 Eugene d'Eon, www.eugenedeon.com/hitchhikers

We cover a variety of topics such as
* phase functions
* NDFs for rough surfaces
* The albedo problem in a 1D rod, 2D Flatland, 3D and higher dimensional spaces
* Diffusion approximations
* H-functions and solutions to Fredholm integral equations

## Monte Carlo code
Each problem has a python script that compiles and generates MC data to be loaded in the correspoding Mathematica notebook.  The assumption is Mac OS g++ with drand48().  It should be easily portable to other environments by changing the RandomReal() function in random.h.

## Deterministic code
The deterministic exact and approximate solutions of various transport theory problems are provided as Mathematica notebooks that load the corresponding Monte Carlo data to generate plots that ensure the two are in agreement.  Each Mathematica notebook is also saved as a .pdf for convenience (these may get out of sync now and then).
