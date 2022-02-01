# DESCRIPTIONS OF THE SCRIPTS
**Author**: Johannes Borgqvist<br>
**Date**: 2021-12-16<br>
In order to calculate the symmetries of a given model, you can execute the script called *launch\_symmetry\_calculations.py*. This script takes two compulsory extra arguments and there is a third optional argument. For example,  you can run the following command in order to calculate the symmetries of Hydon's model using a set polynomial ans\"atze of degree 2:<br>

*python launch\_symmetry\_calculations.py hydons_model 2*.<br>

Here, the first input argument is the name of an xlsx-file stored in the Input folder and the second input argument is the degree of the polynomial ans\"atze implying that in this case second order polynomials are used. This script is automatically set to terminate after two hours, and in order to change this termination time you can pass a third input argument being the time before termination *in seconds*. Thus, if you want the script to terminate after 4 hours instead you can alter the above command to the following:<br>

*python launch\_symmetry\_calculations.py hydons_model 2 14400*.<br>

This script conducts all calculations by using the other scripts in this folder. Here comes a brief description of these scripts. 

1. *conduct\_symmetry\_calculations.py*: this scripts takes in the model name and the degree in the tangential ans\"atze from the script *launch\_symmetry\_calculations.py* and then it conducts all the symmetry calculations by using the subsequent two scripts in this list,
2. *read\_and\_write\_data.py*: this script reads the xlsx files in the Input folder and it writes all the output files in the Output folder,
3. *symmetry\_toolbox\_first\_order\_ODEs.py*: this script is the main script containing all symmetry based calculations relying on [sympy](https://www.sympy.org/en/index.html),
4. *generate\_report.py*: this script generates a pdf-report using LaTeX that is stored in the Output folder and the content of this report is based on the content of the Input folder. 
