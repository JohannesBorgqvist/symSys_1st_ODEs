# symSys_1st_ODEs
Welcome to the git repositry "*symSys_1st_ODEs*". This is a Python and LaTeX based project which conducts automated symbolic calculations of the symmetries of systems of first order ODEs using sympy.

## Description of the project
This repositry provides an automated framework for calulculating the symmetry generators for a system of ODEs. This repositry is linked to the article (**Reference to future article**). 
![Hydons_ODEs](ODE_sys.jpg)

## Pre-requisites to run the scripts
The scripts have been developed on a laptop with the following operating system:<br>
**Ubuntu 20.10**<br>
with the following properties:<br>
**Kernel: Linux 5.8.0-53-generic**<br>
**Architecture: x86-64**.<br>
The programming has been done in Python, and as an output the script generates a LaTeX-file with the calculated generator. The versions of packages required for the programming scripts are the following:<br>
1. *python*, version 3.8.3,
2. *pandas*, version 1.0.5,
3. *notebook*, version 6.0.3,
4. *anaconda*, version 2020.07 (*not required but convenient*). <br>

The script also generates a LaTeX report with *pdflatex*. The version of *pdflatex* that has been used for this project is the following:<br>
**pdfTeX 3.14159265-2.6-1.40.21 (TeX Live 2020/Debian)**,<br>
**kpathsea version 6.3.2**,<br>
**Copyright 2020 Han The Thanh (pdfTeX) et al.**<br>

## Setting up the project using anaconda
The easiest way to get the scripts to work properly is to install [anaconda](https://docs.anaconda.com/anaconda/install/). When anaconda is installed, the required dependencies in order to run the script is accessed by creating an anconda environment with the provided yml-file. This is done with the following command:<br>
*conda env create -f symSys_1st_ODEs.yml*.<br>
After this, the environment is activated with the command:<br>
*conda activate symSys_1st_ODEs*,<br>
and it is deactivated using<br>
*conda deactivate*.<br>
