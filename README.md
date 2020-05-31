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

* Hicolm.in

* Hicolm.str

After that, run the command `hicolm` in your working directory. Four output files will be create:

1. **Hicolm.out**: Which contains relevant informations about the simulation;
2. **Hicolm.md**: With informations about energy, force, velocity and positions of each atom at each MD cycle;
3. **Hicolm.XSF**: The latest force, velocity and positions of each atom;
4. **Hicolm.AXSF**: Conjunction of frames containing force, velocity and positions.