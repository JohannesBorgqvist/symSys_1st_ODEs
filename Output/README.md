# DESCRIPTIONS OF THE OUTPUT FILES
**Author**: Johannes Borgqvist<br>
**Date**: 2021-12-16<br>
All output files are stored in folders that are named after the name of each model defined by the input files. In each folder, numerous symbolic Matrices are stored as "pickle" files which can be [opened and manipulated afterwards](https://www.datacamp.com/community/tutorials/pickle-python-tutorial). Due to some difficulties of storing undefined symbolic functions, I converted all the symbolic functions to strings. Thus, if you want to obtain the symbolic objects after you have opened these files using pickle, you can convert the strings to symbolic objects using [sympify](https://www.tutorialspoint.com/sympy/sympy_sympify_function.htm). 

Otherwise, the results from the symmetry calculations of each model are stored in a tex-file called *out.tex*. In order to generate a pdf-report via LaTeX which accesses all of these files named *out.tex* you can execute the script with relative path "*../Code/generate\_report.py*" which generates a file called "*summary\_report.tex*" in the current directory. This file can then be converted to a pdf using *pdflatex* for example. 
