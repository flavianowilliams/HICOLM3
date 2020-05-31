[![license](https://img.shields.io/github/license/flavianowilliams/HICOLM?style=plastic)](https://github.com/flavianowilliams/HICOLM/blob/master/LICENSE)
[![latest](https://img.shields.io/github/v/release/flavianowilliams/HICOLM?style=plastic)](https://github.com/flavianowilliams/HICOLM/releases/latest)
![downloads](https://img.shields.io/github/downloads/flavianowilliams/HICOLM/total?style=plastic)
![size](https://img.shields.io/github/repo-size/flavianowilliams/HICOLM?color=yellow&style=plastic)
![releasedate](https://img.shields.io/github/release-date-pre/flavianowilliams/HICOLM?color=brown&style=plastic)

# HICOLM Multi-Methods for Molecules and Condensed Systems

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

Notice: gfortran libraries is necessary to compile the Hicolm program.

# What's new?<sup>:star2:</sup>

* Amber Protocol available in Force Field section.
* Optimization method by steepest descent algorithm.

# Usage

After install the user must create in working directory the following input files:

* Hicolm.in, with relevant informations of structure, MD parameters and force field.

```
@MDRUNNING # direct molecular dynamics

&STRUCT
...
&END

&MD
&END

&OPT
&END

&FORCE
&END
```

* Hicolm.str, with coordinates and species identification of each atom.

After that, run the command `hicolm` in your working directory. Four output files will be create:

1. **Hicolm.out**: Which contains relevant information about the simulation;
2. **Hicolm.md**: With values of energy, force, velocity and positions of each atom at each MD cycle;
3. **Hicolm.XSF**: The latest force, velocity and positions of each atom in XSF format (see [XCrysDen format](http://www.xcrysden.org/doc/XSF.html));
4. **Hicolm.AXSF**: Conjunction of frames containing force, velocity and positions in AXSF format.

For more information, you can run some examples in EXAMPLE directory, or send an issue in [Issues page](https://github.com/flavianowilliams/HICOLM/issues).