# symSys_1st_ODEs
**Author**: Johannes Borgqvist<br>
**Date**: 2021-06-02<br>
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

To run all scripts and generate a report, a bash script called "*run\_all.sh*" has been created. It might be neccessary to get permission to run the file first which can be done by typing "*chmod +x run\_all.sh*" and then the script is executed with the command "*./run\_all.sh*".
## Setting up the project using anaconda
The easiest way to get the scripts to work properly is to install [anaconda](https://docs.anaconda.com/anaconda/install/). When anaconda is installed,create an environment for the project

conda create -n *sym\_sys\_1st\_ODEs* python sympy pandas scipy
conda install -n *sym\_sys\_1st\_ODEs* -c conda-forge pickle5
conda install -n *sym\_sys\_1st\_ODEs* -c anaconda csvkit
conda install -n *sym\_sys\_1st\_ODEs* -c anaconda more-itertools
conda install -n *sym\_sys\_1st\_ODEs* -c conda-forge multiprocess

After this, the environment is activated with the command:<br>
conda activate *sym\_sys\_1st\_ODEs*,<br>
and it is deactivated using<br>
conda deactivate.<br>
