[![license](https://img.shields.io/github/license/flavianowilliams/HICOLM?style=plastic)](https://github.com/flavianowilliams/HICOLM/blob/master/LICENSE)
[![latest](https://img.shields.io/github/v/release/flavianowilliams/HICOLM?color=blue&include_prereleases&style=plastic)](https://github.com/flavianowilliams/HICOLM/tree/v2.0.0)
![size](https://img.shields.io/github/repo-size/flavianowilliams/HICOLM?color=yellow&style=plastic)
![releasedate](https://img.shields.io/github/release-date-pre/flavianowilliams/HICOLM?color=brown&style=plastic)

# HICOLM: Multi-Methods for Molecules and Condensed Systems

Hicolm is an open source package focused in molecular modelling of condensed systems.

<p align="center">
    <img width=500 height=auto src=docs/pictures/infrared.png>
</p>

Below is a list of publications involving the HICOLM program:
<p>F. W. Fernandes, HICOLM: High-Performance Platform of Physical Simulations by Using Low Computational Cost Methods. <a href="https://seer.ufrgs.br/rita/article/view/RITA_VOL26_NR3_90">doi: 10.22456/2175-2745.92486.</a></p>

# Installation

To install just run "as a superuser" the script install.sh in the SRC directory,

```
$ sudo ./install.sh
```

Notice: gfortran library is necessary to compile the Hicolm program.

# Highlighted topics

* Velocity Verlet integrator.
* Emulation of NVT, NPT and NVE ensemble by Berendsen or Nosè-Hoover barostat and thermostat.
* Amber Protocol available in the Force Field section. <sup>:new:</sup>
* Optimization method by steepest descent algorithm. <sup>:new:</sup>

# Usage

After install the user must create in the working directory the following input files:

* Hicolm.in, with relevant informations about the structure, MD parameters and the force field:

```
 # direct molecular dynamics
 
 @MDRUNNING

# System information

&STRUCT    
...
&END

# MD parameters

&MD
...
&END

# Optimization parameters

&OPT
...
&END

# Force Field information

&FORCE
...
&END
```

The detailed description of each directive can be seen in <a href="https://seer.ufrgs.br/rita/article/view/RITA_VOL26_NR3_90">doi: 10.22456/2175-2745.92486.</a></p>

* Hicolm.sys, with coordinates and other information about each molecular type. Below is an example of H_2O condensed system, the first line shows the quantity of molecular types, the second shows the molecular name, quantity of molecules and quantity of atoms of each molecule, respectively. Next shows types of each site followed by the partial charges. After that, the atomic number followed by atomic coordinates of each atom are represented in the next lines.

```
    1
H2O         100    3
 OW HW HW
-0.8340 0.4170 0.4170
    8      5.29770100     12.22297800     19.92931800
    1      4.72481900     13.03438700     19.81360300
    1      5.30418400     11.42751000     20.53525100
    ...
```

After that, run the command `hicolm` in your working directory. Four output files will be create:

1. **Hicolm.out**: Contains relevant informations about the simulation;
2. **Hicolm.md**: With values of energy, force, velocity and positions of each atom at each MD cycle;
3. **Hicolm.XSF**: The latest force, velocity and positions of each atom in XSF format (see [XCrysDen format](http://www.xcrysden.org/doc/XSF.html));
4. **Hicolm.AXSF**: Conjunction of frames containing force, velocity and positions in AXSF format.

For more information, you can run the examples in the examples directory, or send an issue in [Issues page](https://github.com/flavianowilliams/HICOLM/issues).