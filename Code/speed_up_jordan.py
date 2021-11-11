# Import some useful packages
import matplotlib.pyplot as plt
import seaborn
import timeit



# Benchmark jordan form
time_jordan = timeit.timeit(
    stmt = 'f(np.random.random())',
    setup='# Import sympy as well;\
           import sympy as sp;\
           from sympy.utilities.autowrap import autowrap;\
           import numpy as np;\
           # Create a symbol;\
           t = sp.Symbol("t");\
           # Define the dimension of the matrix at hand;\
           dimension = 17;\
           # Create a zero matrix;\
           B = sp.zeros(dimension);\
           # Add all non-zero elements;\
           B[0,13] = 1;\
           B[1,15] = 1/(2*t);\
           B[3,16] = 1;\
           B[6,15] = 2;\
           B[7,16] = 8;\
           B[12,16] = 4*t;\
           B[15,15] = 4/t;\
           hello = B.jordan_form();\
           P = hello[0];\
           J = hello[1];\
           J = t*J;\
           exp_J = J.exp();\
           f1 = autowrap(hello, backend="cython", args=(t,));\
           #f2 = autowrap(exp_J, backend="cython", args=(t,));\' 
          )
