# Hydro-PED: Hydro Poro-elastic damage simulation

Hydro-PED models nucleation and propagation of damage zones and seismicity patterns induced by
wellbore fluid injection. The model formulation of Hydro-PED accounts for the following general aspects of brittle rock deformation: (1) Nonlinear elasticity that connects the effective elastic moduli to a damage variable and loading conditions; (2) Evolution of the damage variable as a function of the ongoing deformation and gradual conversion of elastic strain to permanent inelastic deformation during material degradation; (3) Macroscopic brittle instability at a critical level of damage and related rapid conversion of elastic strain to permanent inelastic strain; (4) Coupling between deformation and porous fluid flow through poro-elastic constitutive relationships incorporating damage rheology with Biot's poroelasticity.

**This subversion contains only the mechanics module. An asynchronous offload to two Xeon-Phi KNC co-processors is implemented.**

A Doxygen documentation of the code is available at ([https://harellevin.github.io/hydroped](https://harellevin.github.io/hydroped)).

# Instructions
## Requirements

1. Intel Fortran compiler with compatibility to Xeon-Phi KNC offload.
 

## Compilation

 1. Edit make.inc as needed.
 2. make
 3. ./poro.exe

The tool *iaja.f* creates CSR coordinates according to nodes.dat and elements.dat.
The tool *dat2tec.f90* converts *.nodes.dat and *.tetra.dat output files to tecplot files.

# Citation

If you use the code for science or any form of scientific and technical dissemination activity, we kindly ask to cite the code using the following references:

 - Shalev, E., Lyakhovsky, V.: Modeling reservoir stimulation induced by wellbore fluid injection. In: Thirty Eighth Workshop on Geothermal Reservoir Engineering, Stanford University Stanford, California (2013)
 - Levin, H., Oren, G., Shalev, E., Lyakhovsky, V.: Acceleration of Hydro Poro-elastic Damage Simulation in a Shared-Memory Environment, ([https://bit.ly/2LKUmC5](https://www.researchgate.net/publication/333651797_Acceleration_of_Hydro_Poro-elastic_Damage_Simulation_in_a_Shared-Memory_Environment)) (2019)
