# =================================================================================
# =================================================================================
# Script:"print_equations_phase_plane_lin_sym.py"
# Date: 2022-06-29
# Implemented by: Johannes Borgqvist
# Description:
# Here, we implement a script for just printing the equations stemming from
# plugging in the polynomial ansatze into the linearised symmetry condition.
# =================================================================================
# =================================================================================
# Import Libraries
# =================================================================================
# =================================================================================
# Import sympy
from sympy import *
# Translate a string to symbolic expression
from sympy.parsing.sympy_parser import parse_expr
# Finding monomials
from sympy.polys.monomials import itermonomials, monomial_count
# Ordering monomials 
from sympy.polys.orderings import monomial_key
# For printing to a file
import sys
# Import mathematical function
import math as m
# To extract all subsets of the coefficient vector
#from itertools import chain, combinations
import itertools
# Import numpy as well
import numpy as np
# Import the multiprocessing package for parallelisation
import multiprocessing as mp # For parallelisation
# Import time for timing
import time
from datetime import datetime
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# FUNCTION 1: "create_tangent_ansatze"
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# The function takes two inputs:
#1. The states in the list states, 
#2. The degree of the polynomial ("degree_polynomial").
# It returns three outputs:
# 1.The unknown coefficients in the tangential ansatze,
# 2. A list of all tangents ("eta_list"),
# 3. The list of all monomials.
# The function uses various functionalities in sympy and the polynomials are generated using the built in functions in the polys.monomial part of the sympy library. To automate the generation of the sympy functions and symbol, the "exec" command is used at multiple places in the function.
def create_tangential_ansatze(states,degree_polynomial):
    #----------------------------------------------------------------------------------
    # INITIALISATION: ALLOCATE MEMORY FOR OUR OUTPUT
    #----------------------------------------------------------------------------------
    # Calculate the number of variables and number of states
    num_of_states = len(states)
    # Allocate our to lists which will be the output:
    eta_list = [] # The function with the tangents        
    c = [] # The list of the coefficients
    polynomial_list = [] # The function with the tangents        
    #----------------------------------------------------------------------------------
    # STEP 1: GENERATE POLYNOMIALS IN THE TEMPORARY VARIABLES
    #----------------------------------------------------------------------------------
    # Generate all monomials for our polynomials
    M = list(itermonomials(states, degree_polynomial))
    # Sort the list
    M = sorted(M, key=monomial_key('lex', states))
    # Calculate the number of terms in each of the polynomials
    num_of_monomials = monomial_count(num_of_states,degree_polynomial)
    #----------------------------------------------------------------------------------
    # STEP 2: DEFINE THE POLYNOMIAL ANSATZE
    #----------------------------------------------------------------------------------
    # Allocate our common polynomial
    P_0 = symbols('P_0', cls=Function)
    # Loop over the states and allocate our polynomial
    for state_index in range(2*num_of_states):
        # Initialise it to zero
        P_0 = 0    
        # Loop over the number of terms
        for monomial_index,monomial in enumerate(M):
            # Define the current index of our coefficient
            index = monomial_index + state_index*num_of_monomials
            # Allocate a coefficient
            exec("c_%d = symbols(\'c_%d\') "%(index,index)) 
            # Add this coefficient to the set of unknowns
            exec("c.append(c_%d)"%(index))
            # Add this coefficent to our polynomial at hand
            P_0+=c[-1]*M[monomial_index]
        polynomial_list.append(P_0)
    # When we are done we add P_0 to our list of tangential ansatze
    eta_list.append(polynomial_list[0]/polynomial_list[1])
    eta_list.append(polynomial_list[2]/polynomial_list[3])
    # Return the output    
    return c, eta_list, M, polynomial_list
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# FUNCTION 2: "lin_sym_phase_plane"
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# This function sets up the linearised symmetry condition for the phase plane. 
# It takes three inputs:
# 1. The states,
# 2. The reaction terms,
# 3. The tangents which are rational functions,
# 4. The polynomials making up the tangents.
# It returns three outputs:
# 1. The linearised symmetry condition in a variable called lin_sym,
# 2. The system of equation stemming from the monomials being linearly independent,
# 3. The monomials that each equation stems from.
def  lin_sym_phase_plane(states,omega,polynomial_list):
    # Allocate memory for our outputs
    lin_sym = [] # The linearised symmetry condition
    eq_sys = [] # The system of equations
    monomial_sys = [] # The monomials each equation stems from
    #----------------------------------------------------------------
    # STEP 1: Calculate the partial derivatives of the rational functions
    #----------------------------------------------------------------
    # Here, we only calculate the nominators of the of our lovely derivatives
    der_1 = Derivative(polynomial_list[0],states[0]).doit()*polynomial_list[1]-Derivative(polynomial_list[1],states[0]).doit()*polynomial_list[0]
    der_2 = Derivative(polynomial_list[0],states[1]).doit()*polynomial_list[1]-Derivative(polynomial_list[1],states[1]).doit()*polynomial_list[0]
    der_3 = Derivative(polynomial_list[2],states[0]).doit()*polynomial_list[3]-Derivative(polynomial_list[3],states[0]).doit()*polynomial_list[2]
    der_4 = Derivative(polynomial_list[2],states[1]).doit()*polynomial_list[3]-Derivative(polynomial_list[3],states[1]).doit()*polynomial_list[2]
    #----------------------------------------------------------------
    # STEP 2: SETUP THE LINEARISED SYMMETRY CONDITION
    #----------------------------------------------------------------
    # The linearised symmetry condition is hard coded
    temp_eq = 0
    # The LHS
    temp_eq += (omega[0]**2)*der_3*(polynomial_list[1]**2)
    temp_eq += (omega[0]*omega[1])*((der_4*(polynomial_list[1]**2))-(der_1*(polynomial_list[3]**2)))
    temp_eq += -(omega[1]**2)*der_2*(polynomial_list[3]**2)
    # The RHS
    temp_eq += -(Derivative(omega[1],states[0]).doit()*omega[0]- omega[1]*Derivative(omega[0],states[0]).doit())*polynomial_list[0]*polynomial_list[1]*(polynomial_list[3]**2)
    temp_eq += -(Derivative(omega[1],states[1]).doit()*omega[0]- omega[1]*Derivative(omega[0],states[1]).doit())*polynomial_list[2]*polynomial_list[3]*(polynomial_list[1]**2)
    # Simplify our linearised symmetry condition by multiplying with the product of the squares of the 
    # polynomials in the denominators of the rational functions in the tangential ansatze
    temp_eq = expand(temp_eq)
    # Append our lovely equation
    lin_sym.append(temp_eq)
    #----------------------------------------------------------------
    # STEP 3: EXTRACT THE EQUATIONS STEMMING FROM THE LIN SYM
    #----------------------------------------------------------------    
    # Generate all monomials of a high degree
    M = list(itermonomials(states, 20))
    # Sort the list
    M = sorted(M, key=monomial_key('lex', states))
    # Loop through the monomials, and extract all coefficients that are non-zero
    for m_index,monomial in enumerate(M):
        # Calculate the coefficient at hand
        temp_eq = lin_sym[0].coeff(monomial)
        # Set all higher order monomials to zero
        for state_index,state in enumerate(states):
            temp_eq = temp_eq.subs(state,0)
        # If we have a non-zero equation we append it to our list of equations
        if temp_eq != 0:
            eq_sys.append(temp_eq)
            monomial_sys.append(monomial)
    # Return the output 
    return lin_sym,eq_sys,monomial_sys
# =================================================================================
# =================================================================================
# The LV-model
# =================================================================================
# =================================================================================
print("=======================================================================")
print("\tLV-model")
print("=======================================================================\n")
# Dependent variables
u,v = symbols('u v')
# Define our parameter a in the LV model
a = symbols('a')
# Define the states
states_LV = [u, v]
# Define our reaction terms
omega_LV = [u*(1-v), a*v*(u-1)]
# Create polynomial ansatze
c_LV, eta_LV, M_LV, poly_LV = create_tangential_ansatze(states_LV,2)
# Plug in the tangential ansatze and calculate the resulting equations
lin_sym_LV,eq_sys_LV,monomials_LV = lin_sym_phase_plane(states_LV,omega_LV,poly_LV)
# Write the solutions to file
#Open text file
text_file = open("../Notes/rational_ansatze_lin_sym/Input/generators_LV.tex", "w")
# Write the equations
text_file.write("The LV-model in the phase plane is given by:\n\\begin{equation}\n")
text_file.write("\\dv{%s}{%s}=\\dfrac{%s}{%s}\n\\label{eq:LV}\n"%(latex(states_LV[1]),latex(states_LV[0]),latex(omega_LV[1]),latex(omega_LV[0])))    
text_file.write("\\end{equation}\n")
# Write the monomials
text_file.write("The monomials are:\n%s\n"%(latex(Matrix(len(M_LV),1,M_LV),mode='equation').replace("\\begin{equation}","\\begin{equation}\nM=").replace("\\end{equation}",".\n\\end{equation}\n")))
# Write the coefficient
text_file.write("The unknown coefficient in our polynomial ansatze are:\n%s\n"%(latex(Matrix(len(c_LV),1,c_LV),mode='equation').replace("\\begin{equation}","\\begin{equation}\n\\mathbf{c}=").replace("\\end{equation}",".\n\\end{equation}\n")))
# Write the polynomials
text_file.write("Now, we use rational ansatze of the type $\\eta_1=P_1/P_2$ and $\\eta_2=P_3/P_4$. Our four polynomials are:\n\n\\begin{align*}\n")
text_file.write("P_1&=%s,\\\\\n"%(latex(poly_LV[0])))
text_file.write("P_2&=%s,\\\\\n"%(latex(poly_LV[1])))
text_file.write("P_3&=%s,\\\\\n"%(latex(poly_LV[2])))
text_file.write("P_4&=%s.\\\\\n"%(latex(poly_LV[3])))
text_file.write("\\end{align*}\n")
# Write the align environment
text_file.write("The system of equations resulting from plugging in these ansatze into the linearised symmetry condition contains %d equations.\n"%(len(eq_sys_LV)))
#text_file.write("\\begin{align*}\n")
# Loop over the output and write them all to our file
#for eq_temp in eq_sys_LV:
    #write string to file
    #text_file.write("%s&=0,\\\\\n"%(latex(eq_temp)))
#text_file.write("\\end{align*}\n")    
print("=======================================================================\n")
# Close file
text_file.close()
