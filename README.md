# PoreDy

PoreDy: A Porous Media Simulation & Design Software Coupling Pore-Level Dynamics with Overflow Fluid Mechanics

# Description

PoreDy is a CFD tool originally designed for simulating porous media drying problems, which is capable of couple solving the pore-scale porous media interior and lab-scale outside field simultaneously. Right now, PoreDy is mainly used in 3 categories of problems: general fluid flow, porous media interior transportation, and transportation between porous media and external flow field. The results can be output in VTK format, which can be post-processed with ParaView.

PoreDy is written in C++ with OOP technology, which make it easier to extend to include more physical models, algorithms and schemes. Free linear library Eigen is used for the linear equations solver, which keeps the cost low and is the only module currently capable of parallel computing.

## Major Features for Pore-Scale Porous Media Interior:
	Implicit time marching solver for 2D/3D pore network evaporation model with invasion-percolation rule and film effect considered
	Complex shape pore network is acceptable, 2D regular/irregular pore network generator is included
## Major Features for Lab-Scale Outside Field:
	Implicit 2D/3D steady flow finite volume solver for Navier-Stokes equations
	Implicit 2D/3D time marching solver for convection and diffusion evaporation model
	Complex geometry is acceptable, 2D meshing tool is included
	Optional Boundary Conditions: Inlet, Outlet, Periodic, Static Wall, Sliding Wall, Porous Media Opening
	Optional Algorithms for Flow Solver: SIMPLE, SIMPLER
	Optional Algorithm for Evaporation Solver: Convection-Diffusion Operator Splitting, Diffusion-Convection Operator Splitting, Non Operator Splitting
	Optional Schemes: Hybrid, Hayase QUICK

# Release Note

## Release 1.0 (Dec 24th, 2015)
	Couple solver is optimized
## Release 0.8 (Aug 30th, 2015)
	The film effect is added and tested
	The software is tested on clusters of UWM HPC center
## Release 0.6 (May 4th, 2015)
	The pore network interior solver and outside solver is coupled together
	All solvers are extended to 3D
	A 2D meshing tool is added and tested
	VTK format outputting module for pore network is added and tested; MATLAB figure outputting module for pore network is dropped
## Release 0.4 (Nov 13th, 2014)
	2D/3D invasion-percolation pore network model is added and tested
	MATLAB figure outputting module for pore network is added and tested
	2D convection and diffusion evaporation solver is added and tested
	Eigen is used to replace our home made linear solver
## Release 0.2 (Oct 5th, 2013)
	2D Navier-Stokes Solver is developed and tested
	VTK format outputting module for outside field is developed and tested


