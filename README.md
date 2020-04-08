# Hydro-PED: Hydro Poro-elastic damage simulation

Hydro-PED models nucleation and propagation of damage zones and seismicity patterns induced by
wellbore fluid injection. The model formulation of Hydro-PED accounts for the following general aspects of brittle rock deformation: (1) Nonlinear elasticity that connects the effective elastic moduli to a damage variable and loading conditions; (2) Evolution of the damage variable as a function of the ongoing deformation and gradual conversion of elastic strain to permanent inelastic deformation during material degradation; (3) Macro-
scopic brittle instability at a critical level of damage and related rapid conversion of elastic strain to permanent inelastic strain; (4) Coupling between deformation and porous fluid flow through poro-elastic constitutive relationships incorporating damage rheology with Biot's poroelasticity.

A Doxygen documentation of the code is available at ([https://harellevin.github.io/hydroped](https://harellevin.github.io/hydroped)).

# Instructions
## Requirements

1. C++ and Fortran compiler (as GNU, Intel or XL).
1. BLAS implementation (OpenBLAS, MKL, ESSL, etc.).
1. LLNL's Trilinos ([link](https://trilinos.github.io/)), built with one the following packages:
	 - Belos  
	 - TeuchosCore 
	 - TeuchosComm 
	 - TeuchosRemainder 
	 - TeuchosParameterList 
	 - TeuchosNumerics 
 1. And one of the following Trilinos package combinations:
 
 |Templated Version|Basic Version  |
 |--|--|
 | Tpetra | Epetra|
 | BelosTpetra | BelosEpetra |
 | Kokkos | |
 

## Compilation

 1. Edit the places marked on make.inc.AIX or make.inc.Linux (depends on your OS).
 2. make
 3. ./poro.exe

The tool *iaja.f* creates CSR coordinates according to nodes.dat and elements.dat.
The tool *dat2tec.f90* converts *.nodes.dat and *.tetra.dat output files to tecplot files.

# Citation

If you use the code for science or any form of scientific and technical dissemination activity, we kindly ask to cite the code using the following references:

 - Shalev, E., Lyakhovsky, V.: Modeling reservoir stimulation induced by wellbore fluid injection. In: Thirty Eighth Workshop on Geothermal Reservoir Engineering, Stanford University Stanford, California (2013)
Additionally, please cite
```
@inproceedings{DBLP:conf/parco/LevinOSL19.bib,
  author    = {Harel Levin, 
               Gal Oren, 
               Eyal Shalev and 
               Vladimir Lyakhovsky},
  editor    = {Ian Foster and
               Gerhard R. Joubert and
               Luděk Kučera and
               Wolfgang E. Nagel and
               Frans Peters},
  title     = {Acceleration of Hydro Poro-Elastic Damage Simulation in a Shared-Memory Environment},
  booktitle = {Parallel Computing: Technology Trends, Proceedings of the International
               Conference on Parallel Computing, Prague, Czech Republic, 10-13 September 2019},
  series    = {Advances in Parallel Computing},
  volume    = {36},
  pages     = {341--353},
  publisher = {{IOS} Press},
  year      = {2020},
  doi       = {10.3233/APC200075},
  biburl    = {https://dblp.org/rec/conf/parco/LevinOSL19.bib},
  bibsource = {dblp computer science bibliography, https://dblp.org}
}
```
