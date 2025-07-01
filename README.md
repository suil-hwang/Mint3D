This repository contains a reference implementation of the 🍀[paper]{https://www.cs.utexas.edu/~josh/papers/Mint3d_final.pdf}🍀: 
# Mint: Discretely Integrable Moments for Symmetric Frame Fields
to be published in SGP 2025 proceedings. 

![Mint3D teaser](/teaser.jpg "A menagerie")


Please go to: https://github.com/the13fools/Mint3D to find the latest version of this code.  

What is contained in this repo?  

- 🤝 Symmetric Frame Field Solver 
    - This app implements a solver which is by default configured to recover approximately orthonormal, integrable, frame fields subject to symmetric integrability and boundary alignment constraints.  
- 🤝 Symmetric Frame Field Visualizer 
    - This app provides a front end for (pre)visualizing symmetric frame fields by plotting field streamlines, best fit global parameterizations, and other quality metrics.
- 🤝 Reference CubeCover Implementation 
    - We include an independent implementation of cubecover in this code release with a minimal amount of bells and whistles.  This implementation is callable from the visualizer, not as a stand alone app in this implementation.

--------

# How to Guide

Build the code using cmake as described below.  Then run by calling (from the build directory)

```./bin/00_basic_test``` 
to check if the build worked. 

If so can call 
``` ./bin/01_Mint3D -h ```
to get the commandline options for the main app.  If you want to directly launch the viewer, can also do this from commandline
```./bin/02_SymFrameView -h```   

Can try to run: 
```./bin/01_Mint3D -m ../tools/shared/cube_lil.mesh```
Can launch parameteriation from within the app to see the fields viewer.   

To view the results from an earlier solve, can load and replay the results by using the -d flag, e.g. 

```./bin/01_Mint3D -d ../../results/orthog/4_9_7_21__cube_lil_ss_1.00e+00_is_1.00e-08_o_1.00e-01_u_1.00e-05_vsc_1.00e-09_plane_1e-6_0.00e+00_bsdf```

*PLEASE LET THE SOLVER RUN TO COMPLETION BEFORE USING FIELDS FOR APPLICATIONS UNLESS YOU ARE **SURE** THEY ARE SUFFICIENTLY CONVERGED FOR YOUR NEEDS*

We include the specific meshes that we tested our solver on in this paper in the /meshes/ directory.


--------

This example repository includes make files potentially with extra dependencies, and this will be addressed in a patch post publication.  

We are including a repository with some scripts we used but these are mostly for reference in case anyone needs them, however this is not a batteries included repo at the moment.  Scripts will likely need to be tweaked to work on a different machine, and might be removed in a subsequent patch.  



--------

This repository includes a lot of gpt generated code, but its a cyborg project with lots of other snippets mixed in.

We wrote this code with the help of contemporary language models circa the big recent boom, and it's perhaps an artifact which captures what was possible to do with the tools available as this project was being developed, as has always been the case.  



## Compile

Compile this project using the standard cmake routine:

    mkdir build
    cd build
    cmake ..
    make

This should find and build the dependencies and create a `bin/##` binary.

## Run

From within the `build` directory just issue:

    ./example



# MAC BUILD DETAILS 

fresh install details: 
brew install cmake eigen suitesparse llvm 

there might be a compilation issue where the following file path needs to have the "register" keyword removed 
Mint3D/build/_deps/comiso-src/ext/gmm-4.2/include/gmm

------------

This might not be an issue if the openmp dependency is removed, but in case you need it: 

Thank you eerii .  Anonomous internet comments really came through here.  

https://github.com/glfw/glfw/issues/1743#issuecomment-1229177189

For now until we improve the cmake file we added the following to our ~/.zshrc

```
export CC=$(which clang)
export CXX=$(which g++)
export OpenMP_ROOT=$(brew --prefix)/opt/libomp
```



# Linux Build Details 

If you are running the fields viewer and it's not opening stuff, it's probably because the stack space is getting filled up.  Try running the following command before calling the viewer (note this is automatically called when the viewer is launched from the mint3d solver gui)

```
ulimit -s 512000
```

Addendum:

In the spirit of Buckminster Fuller, who changed how humanity understood parameterization in the 20th century - may we use these new tools at our disposal be good stewards of our lands, and to improve spaceship earth for future generations. ðŸ?€


