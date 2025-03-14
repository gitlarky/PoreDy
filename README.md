PoreDy
=

PoreDy: A Porous Media Evaluation & Design Software Featuring Coupled Pore-Lab Scale Dynamics Solver

## Announcement:

PoreDy is now open sourced under Apache 2.0 license. The governance of further development is transfered from Dr. Zhenyu Xu to the Composite Laboratory of University of Wisconsin Milwaukee under Professor Krishna M Pillai's supervision. The updating on this repository would not be guaranted since Release1.2 on 2020-04-15, but you can always find the newest progress in Prof. Pillai's fork page.

You are encouraged to use this code in any academic or industrial work, my research works and their simulation-setting files are archived in the test folder. Any form of citations are welcomed but not required.

## Description

PoreDy is a tool originally designed for simulating porous media drying problems, which is capable of couple solving the pore-scale porous media interior and lab-scale outside field simultaneously. Right now, PoreDy is mainly used in three categories of problems: general fluid flow, porous media interior transportation, transportation between porous media and external flow field. The results can be output in VTK format, which can be post-processed with ParaView.

PoreDy is written in C++ and Python with OOP technology, which make it easier to extend to include more physical models, algorithms and schemes. Open source linear library Eigen is used for the linear equations solver, which keeps the cost low and is the only module currently capable of parallel computing. 

#### Major Features for Pore-Scale Porous Media Interior Solver:
1. Implicit time marching solver for 2D/3D pore network evaporation model with invasion-percolation rule and film effect considered
2. Complex shape pore network is acceptable, 2D regular/irregular pore network generator is included

#### Major Features for Lab-Scale Outside Field Solver:
1. Implicit 2D/3D steady flow finite volume solver for Navier-Stokes equations
2. Implicit 2D/3D time marching solver for convection and diffusion evaporation model
3. Complex geometry is acceptable, 2D meshing tool is included
4. Optional Boundary Conditions: Inlet, Outlet, Periodic, Static Wall, Sliding Wall, Porous Media Opening
5. Optional Algorithms for Flow Solver: SIMPLE, SIMPLER
6. Optional Algorithm for Evaporation Solver: Convection-Diffusion Operator Splitting, Diffusion-Convection Operator Splitting, Non Operator Splitting
7. Optional Schemes: Hybrid, Hayase QUICK


## Release Note

#### Release 1.2 (Apr 15th, 2020)
- [x] The feature of receiving individually assigned throat pore-network added
- [x] The specific pore-network generator added
- [x] Bugs on uniform calculation segment fault fixed

#### Release 1.0 (Dec 24th, 2015)
- [x] Release downloadable published
- [x] Couple solver is optimized
- [x] Cases setting files for testing Version 1.0 are added
- [x] The software is tested on clusters of UWM HPC center

#### Release 0.6 (May 4th, 2015)
- [x] Release downloadable published
- [x] The film effect is added and tested
- [x] The pore network interior solver and outside solver is coupled together
- [x] All solvers are extended to 3D
- [x] A 2D meshing tool is added and tested
- [x] VTK format outputting module for pore network is added and tested; MATLAB figure outputting module for pore network is dropped

#### Release 0.4 (Nov 13th, 2014)
- [x] 2D/3D invasion-percolation pore network model is added and tested
- [x] MATLAB figure outputting module for pore network is added and tested
- [x] 2D convection and diffusion evaporation solver is added and tested
- [x] Eigen is used to replace our home made linear solver

#### Release 0.2 (Oct 5th, 2013)
- [x] 2D Navier-Stokes Solver is developed and tested
- [x] VTK format outputting module for outside field is developed and tested
