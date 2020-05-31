[![license](https://img.shields.io/github/license/flavianowilliams/HICOLM?style=plastic)](https://github.com/flavianowilliams/HICOLM/blob/master/LICENSE)
[![latest](https://img.shields.io/github/v/release/flavianowilliams/HICOLM?style=plastic)](https://github.com/flavianowilliams/HICOLM/releases/latest)
![downloads](https://img.shields.io/github/downloads/flavianowilliams/HICOLM/total?style=plastic)
![size](https://img.shields.io/github/repo-size/flavianowilliams/HICOLM?color=yellow&style=plastic)
![releasedate](https://img.shields.io/github/release-date-pre/flavianowilliams/HICOLM?color=brown&style=plastic)

# HICOLM: Multi-Methods for Molecules and Condensed Systems

Hicolm is an open source package focused in molecular modelling of condensed systems.

<p align="center">
    <img width=500 height=auto src=DOCS/pictures/input_file.png>
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
* Optimization method by steepest descent algorithm. <sup>:new2:</sup>

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

* Hicolm.str, with coordinates and species identification of each atom. Below is an example of C2H4 molecule, the first collumn in the atomic number, from the second to fourth collumns are the atomic coordinates. Next are the specie identification and the last one is the freeze option of each atom (0 freeze, 1 free).

```
    6    0.00632000    0.00000000    0.00000000        1      1
    6    1.34000000    0.00000000    0.00000000        2      1
    1   -0.53868000   -0.94397000    0.00000000        3      1
    1   -0.53868000    0.94397000    0.00000000        4      1
    1    1.88500000    0.94397000    0.00000000        5      1
    1    1.88500000   -0.94397000    0.00000000        6      1
```

After that, run the command `hicolm` in your working directory. Four output files will be create:

1. **Hicolm.out**: Contains relevant informations about the simulation;
2. **Hicolm.md**: With values of energy, force, velocity and positions of each atom at each MD cycle;
3. **Hicolm.XSF**: The latest force, velocity and positions of each atom in XSF format (see [XCrysDen format](http://www.xcrysden.org/doc/XSF.html));
4. **Hicolm.AXSF**: Conjunction of frames containing force, velocity and positions in AXSF format.

For more information, you can run the examples in the EXAMPLE directory, or send an issue in [Issues page](https://github.com/flavianowilliams/HICOLM/issues).