# =================================================================================
# =================================================================================
# Script:"phase_plane.py"
# Date: 2022-06-30
# Implemented by: Johannes Borgqvist
# Description:
# The script calculates the infinitesimal generators of the Lie group restricted
# to the phase plane. Here, we focus on the Lotka-Volterra (LV)-model with
# competition. 
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
# =================================================================================
# =================================================================================
# FUNCTIONS
# =================================================================================
# =================================================================================
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# FUNCTION 1: "create_tangent_ansatze"
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# The function takes three inputs:
#1. The states in the list states, 
#2. The degree of the polynomial ("degree_polynomial"),
# 3. The list of reaction terms ("omega_list")
# It returns three outputs:
# 1.The unknown coefficients in the tangential ansatze,
# 2. A list of all tangents ("eta_list"),
# 3. The list of all monomials.
# The function uses various functionalities in sympy and the polynomials are generated using the built in functions in the polys.monomial part of the sympy library. To automate the generation of the sympy functions and symbol, the "exec" command is used at multiple places in the function.
def create_tangential_ansatze(states,degree_polynomial,omega_list):
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
    M_high = list(itermonomials(states, degree_polynomial+1))
    # Sort the list
    M_high = sorted(M_high, key=monomial_key('lex', states))
    # Calculate the number of terms in each of the polynomials
    num_of_monomials_high = monomial_count(num_of_states,degree_polynomial+1)
    # Generate all monomials for our polynomials
    M_low = list(itermonomials(states, degree_polynomial-1))
    # Sort the list
    M_low = sorted(M_low, key=monomial_key('lex', states))
    # Calculate the number of terms in each of the polynomials
    num_of_monomials_low = monomial_count(num_of_states,degree_polynomial-1)    
    #----------------------------------------------------------------------------------
    # STEP 2: DEFINE THE POLYNOMIAL ANSATZE
    #----------------------------------------------------------------------------------
    # Allocate our common polynomial
    P_0 = symbols('P_0', cls=Function)
    # Loop over the states and allocate our polynomial
    for state_index in range(4):
        # Initialise it to zero
        P_0 = 0    
        # Loop over the number of terms
        for monomial_index,monomial in enumerate(M_low):
            # Define the current index of our coefficient
            index = monomial_index + state_index*num_of_monomials_low
            # Allocate a coefficient
            exec("c_%d = symbols(\'c_%d\') "%(index,index)) 
            # Add this coefficient to the set of unknowns
            exec("c.append(c_%d)"%(index))
            # Add this coefficent to our polynomial at hand
            P_0+=c[-1]*M_low[monomial_index]
        polynomial_list.append(P_0)
    # Loop over the states and allocate our polynomial
    for state_index in range(2):
        # Initialise it to zero
        P_0 = 0    
        # Loop over the number of terms
        for monomial_index,monomial in enumerate(M_high):
            # Define the current index of our coefficient
            index = (4*num_of_monomials_low) + monomial_index + state_index*num_of_monomials_high
            # Allocate a coefficient
            exec("c_%d = symbols(\'c_%d\') "%(index,index)) 
            # Add this coefficient to the set of unknowns
            exec("c.append(c_%d)"%(index))
            # Add this coefficent to our polynomial at hand
            P_0+=c[-1]*M_high[monomial_index]
        polynomial_list.append(P_0)        
    # Define some nice functionalities
    g1 = omega_list[0]*polynomial_list[0] + omega_list[1]*polynomial_list[1]
    g2 = omega_list[0]*polynomial_list[2] + omega_list[1]*polynomial_list[3]    
    # When we are done we add P_0 to our list of tangential ansatze
    eta_list.append((g1/g2)*polynomial_list[4])
    eta_list.append((g1/g2)*polynomial_list[5])
    # Return the output    
    return c, eta_list, M_high, polynomial_list
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
    # Extract all lower order polynomials
    P1 = polynomial_list[0]
    P2 = polynomial_list[1]
    P3 = polynomial_list[2]
    P4 = polynomial_list[3]
    # Extract all higher order polynomials
    P5 = polynomial_list[4]
    P6 = polynomial_list[5]
    # Extract the reaction terms as well
    omega_1 = omega[0]
    omega_2 = omega[1]
    # Extract the states as well
    u = states[0]
    v = states[1]
    #----------------------------------------------------------------
    # STEP 1: Calculate the partial derivatives of the rational functions
    #----------------------------------------------------------------
    # Define our glue functions as linear combinations of the reaction terms
    g1 = P1*omega_1 + P2*omega_2
    g2 = P3*omega_1 + P4*omega_2    
    # Calculate the partial derivatives of these functions
    du_g1 = (P1*Derivative(omega_1,u).doit() + omega_1*Derivative(P1,u).doit()) + (P2*Derivative(omega_2,u).doit() + omega_2*Derivative(P2,u).doit())
    dv_g1 = (P1*Derivative(omega_1,v).doit() + omega_1*Derivative(P1,v).doit()) + (P2*Derivative(omega_2,v).doit() + omega_2*Derivative(P2,v).doit())
    du_g2 = (P3*Derivative(omega_1,u).doit() + omega_1*Derivative(P3,u).doit()) + (P4*Derivative(omega_2,u).doit() + omega_2*Derivative(P4,u).doit())
    dv_g2 = (P3*Derivative(omega_1,v).doit() + omega_1*Derivative(P3,v).doit()) + (P4*Derivative(omega_2,v).doit() + omega_2*Derivative(P4,v).doit())        
    # Calculate our partial derivatives of the numerators of the tangents
    du_n1 = (du_g1*P5+g1*Derivative(P5,u).doit())*g2 - (g1*P5)*du_g2
    dv_n1 = (dv_g1*P5+g1*Derivative(P5,v).doit())*g2 - (g1*P5)*dv_g2
    du_n2 = (du_g1*P6+g1*Derivative(P6,u).doit())*g2 - (g1*P6)*du_g2
    dv_n2 = (dv_g1*P6+g1*Derivative(P6,v).doit())*g2 - (g1*P6)*dv_g2        
    #----------------------------------------------------------------
    # STEP 2: SETUP THE LINEARISED SYMMETRY CONDITION
    #----------------------------------------------------------------
    # The linearised symmetry condition is hard coded
    temp_eq = 0
    # The LHS
    temp_eq += (omega_1**2)*du_n2
    temp_eq += (omega_1*omega_2)*(dv_n2-du_n1)
    temp_eq += -(omega_2**2)*dv_n1
    # The RHS
    temp_eq += -(Derivative(omega_2,u).doit()*omega_1- omega_2*Derivative(omega_1,u).doit())*g1*g2*P5
    temp_eq += -(Derivative(omega_2,v).doit()*omega_1- omega_2*Derivative(omega_1,v).doit())*g1*g2*P6
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
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# FUNCTION 3: "extract_generators"
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# This function sets up the linearised symmetry condition for the phase plane. 
# It takes four inputs:
# 1. c the coefficients in the initial ansatze,
# 2. LHS being the coefficients we have solved for,
# 3. RHS being the value of these coefficients in terms of other coefficients,
# 4. eta being the list of the initial tangents.
# It returns a single output being the list of all generators as a tuple of the two 
# infinitesimals namely (eta_1,eta_2).
def extract_generators(c,LHS,RHS,eta):
    # Allocate memory of our output
    generator_list = []    
    #----------------------------------------------------------------
    # STEP 1: FIND THE CONSTANTS THAT APPEAR IN THE FINAL SOLUTION
    #----------------------------------------------------------------
    # See which constants that exists among the constants
    constants_final = []
    # Loop over all algebraic expressions
    for RHS_index,RHS_temp in enumerate(RHS):
            # Loop over all candidate constants
            for c_index,c_temp in enumerate(c):
                # See if the constant exists in the final solution
                if RHS_temp.coeff(c_temp) != 0:
                    # Add the constant
                    constants_final.append(c_temp)
    # Only save the unique ones
    constants_final = list(set(constants_final))
    #----------------------------------------------------------------
    # STEP 2: SUBSTITUTE CONSTANTS INTO THE TANGENTS
    #----------------------------------------------------------------
    # Allocate memory for the final tangents
    eta_final = eta
    # Loop through the tangents and substitute the values for the coefficients
    for eta_index,eta_temp in enumerate(eta_final):
    # Loop through the coefficients as well
        for c_index,LHS_temp in enumerate(LHS):
            eta_final[eta_index] = eta_final[eta_index].subs(LHS_temp,RHS[c_index])
        # Simplify in the end to clean it up
        eta_final[eta_index] = expand(eta_final[eta_index])
    #----------------------------------------------------------------
    # STEP 3: EXTRACT GENERATORS
    #----------------------------------------------------------------    
    # ALlocate memory for the trivial generator
    trivial_generator = eta_final
    trivial_generator_returned = []
    # Loop over all our coefficients and save our generators
    for constant_index,constant in enumerate(constants_final):
        # Append the current generator
        generator_list.append((factor(eta_final[0].coeff(constant)),factor(eta_final[1].coeff(constant))))
        # Substract the current generator from the trivial generator
        trivial_generator[0] = expand(trivial_generator[0] - eta_final[0].coeff(constant)*constant)
        trivial_generator[1] = expand(trivial_generator[1] - eta_final[1].coeff(constant)*constant)
    # Lastly, append the trivial generator as well
    if simplify(trivial_generator[0]) != 0 and simplify(trivial_generator[1])!= 0:
        trivial_generator_returned.append((simplify(trivial_generator[0]),simplify(trivial_generator[1])))
    # Return our generators
    return generator_list, trivial_generator_returned
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# FUNCTION 4: "findsubsets"
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
def findsubsets(s, n):
    return list(itertools.combinations(s, n))
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# FUNCTION 5: "filtration_tangential_ansatze"
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# This function takes five inputs:
# 1. coeff_vec being the list of all coefficient vectors,
# 2. c_subsets beint the list of all non-zero subsets of the coefficient vectors,
# 3. poly being the list of the polynomials that constitute our tangential ansatze,
# 4. c being the list of the original unknown coefficients in the tangential ansatze,
# 5. states containg the states in our model.
# This function returns two outputs:
# 1. coeff_vec_reduced being the reduced coefficient list,
# 2. c_subsets_reduced being the reduced list of subsets.
def filtration_tangential_ansatze(coeff_vec,c_subsets,poly,c,states):
    # Define the subsets that we want to return in the end
    coeff_vec_reduced = coeff_vec
    c_subsets_reduced = c_subsets
    # Allocate memory of the indices we want to remove in the end
    remove_indices = []
    # Loop through our coefficient candidates and remove the ones resulting in unreasonable ansatze
    for coeff_index, coeff_vec_temp in enumerate(coeff_vec_reduced):
        # Extract our two polynomials in the denominator
        P1 = poly[0]
        P2 = poly[1]
        P3 = poly[2]
        P4 = poly[3]
        P5 = poly[4]
        P6 = poly[5]        
        # Loop through the coefficients in c and substitute our candidate values
        for c_index, coeff_temp in enumerate(c):
            P1 = P1.subs(coeff_temp,coeff_vec_temp[c_index])
            P2 = P2.subs(coeff_temp,coeff_vec_temp[c_index])
            P3 = P3.subs(coeff_temp,coeff_vec_temp[c_index])
            P4 = P4.subs(coeff_temp,coeff_vec_temp[c_index])
            P5 = P5.subs(coeff_temp,coeff_vec_temp[c_index])
            P6 = P6.subs(coeff_temp,coeff_vec_temp[c_index])            
        # Loop through our states and substitute them by zero
        for state in states:
            P1 = P1.subs(state,0)
            P2 = P2.subs(state,0)
            P3 = P3.subs(state,0)
            P4 = P4.subs(state,0)
            P5 = P5.subs(state,0)
            P6 = P6.subs(state,0)                    
        # Remove indices in the following cases
        if P3==0 and P4==0: # Denominator is zero=>things blow up!
            remove_indices.append(coeff_index)
        elif P5==0 and P6==0:# Both nominators are zero=> Trivial generator!
            remove_indices.append(coeff_index)
        elif P1 ==0 and P2 == 0: # Both nominators are zero=> Trivial generator!
            remove_indices.append(coeff_index)            
    # Sort our remove indices
    remove_indices.sort(reverse=True)
    # Lastly, remove the unneccessary coefficients
    for remove_index in remove_indices:
        del coeff_vec_reduced[remove_index]
        del c_subsets_reduced[remove_index]    
    # Return our reduced lists
    return coeff_vec_reduced, c_subsets_reduced
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# FUNCTION 6: "solve_algebraic_system"
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# This function solves an individual algebraic system stemming from plugging in a 
# coefficient vector into the original massive system of equations. 
# It takes four inputs:
# 1. coeff_vec the coefficients we are plugging in,
# 2. c the original coefficients,
# 3. c_subset send in the correspond,
# 3. eta the tangents of the model in question,
# 4. eq_sys being the full system of equations that we want to reduce,
# 5. The polynomials in the ansatze in the list poly,
# 6. The states of the system for fancy printing.
# It returns a single output which is a text string of the calculated solutions. 
def solve_algebraic_system(coeff_vec,c,c_subset,eta,eq_sys,poly,states,current_iteration,final_iteration,omega):
    # Prompt to user
    print("\tIteration %d out of %d started," %(current_iteration,final_iteration))
    # Allocate memory for our text string 
    generator_string = ""
    # We start by creating a temporary tangents and temporary polynomials
    eta_temp = [tangent for tangent in eta]
    poly_temp = [poly_temp for poly_temp in poly]
    # Substitute our candidate coefficient vector into our tangential ansätze
    for coeff_index,coeff_temp in enumerate(coeff_vec):
        # Substitute all zeroes here...
        if coeff_temp == 0:
            eta_temp[0] = eta_temp[0].subs(c[coeff_index],coeff_temp)
            eta_temp[1] = eta_temp[1].subs(c[coeff_index],coeff_temp)
            poly_temp[0] = poly_temp[0].subs(c[coeff_index],coeff_temp)
            poly_temp[1] = poly_temp[1].subs(c[coeff_index],coeff_temp)
            poly_temp[2] = poly_temp[2].subs(c[coeff_index],coeff_temp)
            poly_temp[3] = poly_temp[3].subs(c[coeff_index],coeff_temp)
            poly_temp[4] = poly_temp[4].subs(c[coeff_index],coeff_temp)
            poly_temp[5] = poly_temp[5].subs(c[coeff_index],coeff_temp)            
    # Extract the reduced system of equations
    eq_sys_reduced = []
    # Loop through equations in the original system of equations
    for eq_temp in eq_sys:
        # Extract an individual equations
        temp_eq = eq_temp
        # Loop through all coefficients and substitute these values
        for coeff_index,coeff_temp in enumerate(coeff_vec):
            temp_eq = temp_eq.subs(c[coeff_index],coeff_temp)
        # Save all non-trivial equations after the substitutions
        if temp_eq != 0:
            eq_sys_reduced.append(temp_eq)
    # Solve the system if the reduced system contains equations
    if len(eq_sys_reduced)!=0:
        solutions = solve(eq_sys_reduced,c_subset,set=True)      
        # Allocate memory for the generators
        generators = []
        # Loop through our solutions and substitute the left hand side in
        # solutions[0] with the different right hand sides in solutions[1]
        for RHS_index,RHS in enumerate(solutions[1]):
            # Extract all polynomials
            P1 = poly_temp[0]
            P2 = poly_temp[1]
            P3 = poly_temp[2]
            P4 = poly_temp[3]
            P5 = poly_temp[4]
            P6 = poly_temp[5]                        
            # Substitute coefficients into polynomials
            for c_index,coeff_temp in enumerate(RHS):
                P1 = P1.subs(solutions[0][c_index],coeff_temp)
                P2 = P2.subs(solutions[0][c_index],coeff_temp)
                P3 = P3.subs(solutions[0][c_index],coeff_temp)
                P4 = P4.subs(solutions[0][c_index],coeff_temp)
                P5 = P5.subs(solutions[0][c_index],coeff_temp)
                P6 = P6.subs(solutions[0][c_index],coeff_temp)                                       
            # Save all tangents that we have created that do not blow up
            if (P3!=0 or P4!=0) and (P1!=0 or P2!=0) and (P5!=0 or P6!=0):# Denominator non-zero
                # Define our rational functions g1 and g2
                g1 = omega[0]*P1+omega[1]*P2
                g2 = omega[0]*P3+omega[1]*P4
                # Append our generators in the end
                generators.append([(g1/g2)*P5, (g1/g2)*P6])
        # If we found some generators we save them as a text string
        if len(generators)!=0:
            # Loop through all generators
            for generator in generators:
                # We neglect all trivial ones where both tangents are zero
                if generator[0]!=0 or generator[1]!=0:
                    # u-directional
                    if generator[0]!=0 and generator[1]==0:
                        generator_string += "\\begin{align*}\n"
                        generator_string += "X_{}&=" + latex(generator[0]) + "\\partial_{" + latex(states[0]) + "},\\\\\n"
                        generator_string += "\\end{align*}\n"
                    # v-directional                        
                    elif generator[0]==0 and generator[1]!=0:
                        generator_string += "\\begin{align*}\n"
                        generator_string += "X_{}&=" + latex(generator[1]) + "\\partial_{" + latex(states[1]) + "},\\\\\n"
                        generator_string += "\\end{align*}\n"
                    # Tangents in both directions
                    else:
                        generator_string += "\\begin{align*}\n"
                        generator_string += "X_{}&=" + latex(generator[0]) + "\\partial_{" + latex(states[0]) + "}+" + latex(generator[1]) + "\\partial_{" + latex(states[1]) + "},\\\\\n"
                        generator_string += "\\end{align*}\n"                                
    # Prompt to user again                    
    print("\t\tIteration %d out of %d finished," %(current_iteration,final_iteration))
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("\t\t\tCurrent Time =", current_time)
    # Return our generator string
    return generator_string
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
a, b, p = symbols('a b p')
# Define the states
states_LV = [u, v]
# Define our reaction terms
#omega_LV = [u*(1-v), a*v*(u-1)] # Standard Lotka-Volterra
omega_LV = [u*(1-u-a*v), p*v*(1-v-b*u)] # Lotka-Volterra with competition
# Create polynomial ansatze
print("\tCreate tangential ansatze:")
c_LV, eta_LV, M_LV, poly_LV = create_tangential_ansatze(states_LV,2,omega_LV)
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("\t\t\tDone!\tCurrent Time =", current_time)
# Plug in the tangential ansatze and calculate the resulting equations
print("\tSetup linearised symmetry condition:")
lin_sym_LV,eq_sys_LV,monomials_LV = lin_sym_phase_plane(states_LV,omega_LV,poly_LV)
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("\t\t\tDone!\tCurrent Time =", current_time)
# Find all tuples of subsets of a given size
subset_size = 10
print("\tFind all subsets of size %d:"%(subset_size))
c_subsets_LV = findsubsets(c_LV, subset_size)
# Convert all elements to lists
c_subsets_LV = [list(subset_as_tuple) for subset_as_tuple in c_subsets_LV]
print("\t\tWe have %d such subsets in total..."%(len(c_subsets_LV)))
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("\t\t\tDone!\tCurrent Time =", current_time)
# Filtration step
print("\tOrganise these subsets in terms of coefficient vectors:")
coeff_vec_LV = []
remove_indices = []
# Set the remaining coefficients to zero if they are not in the subsets
for index,c_subset_LV in enumerate(c_subsets_LV):
    c_list = list(np.zeros((len(c_LV)), dtype=int))
    for const_index,const in enumerate(c_LV):
        for sub_index,c_temp in enumerate(c_subset_LV):
            if const == c_temp:
                c_list[const_index] = const
                break
    # Extract our polynomials
    P1 = poly_LV[0]
    P2 = poly_LV[1]
    P3 = poly_LV[2]
    P4 = poly_LV[3]
    P5 = poly_LV[4]
    P6 = poly_LV[5]        
    # Loop through the coefficients in c and substitute our candidate values
    for c_index, coeff_temp in enumerate(c_LV):
        P1 = P1.subs(coeff_temp,c_list[c_index])
        P2 = P2.subs(coeff_temp,c_list[c_index])
        P3 = P3.subs(coeff_temp,c_list[c_index])
        P4 = P4.subs(coeff_temp,c_list[c_index])
        P5 = P5.subs(coeff_temp,c_list[c_index])
        P6 = P6.subs(coeff_temp,c_list[c_index])            
    # Loop through our states and substitute them by zero
    for state in states_LV:
        P1 = P1.subs(state,0)
        P2 = P2.subs(state,0)
        P3 = P3.subs(state,0)
        P4 = P4.subs(state,0)
        P5 = P5.subs(state,0)
        P6 = P6.subs(state,0)                    
    # Remove indices in the following cases
    if P3==0 and P4==0: # Denominator is zero=>things blow up!
        remove_indices.append(index)
    elif P5==0 and P6==0:# Both nominators are zero=> Trivial generator!
        remove_indices.append(index)
    elif P1 ==0 and P2 == 0: # Both nominators are zero=> Trivial generator!
        remove_indices.append(index)
    else: # The candidate coefficient vector is perfectly fine so let's keep it! 
        # Append our list
        coeff_vec_LV.append(c_list)
# Sort our remove indices
remove_indices.sort(reverse=True)
# Lastly, remove the unneccessary coefficients
for remove_index in remove_indices:
    del c_subsets_LV[remove_index]    
# Prompt to the user
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("\t\t\tDone!\tCurrent Time =", current_time)
# Filtration step
#print("\tFiltration step:")
#print("\t\tNumber of coefficients before filtration:\t\t%d"%(len(coeff_vec_LV)))
#print("\t\tNumber of coefficient subsets before filtration:\t%d"%(len(c_subsets_LV)))
# Do the filtration step
#coeff_vec_LV, c_subsets_LV = filtration_tangential_ansatze(coeff_vec_LV,c_subsets_LV,poly_LV,c_LV,states_LV)
print("\t\tNumber of coefficients after filtration:\t\t%d"%(len(coeff_vec_LV)))
print("\t\tNumber of coefficient subsets after filtration:\t\t%d"%(len(c_subsets_LV)))
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("\t\t\tDone!\tCurrent Time =", current_time)
# Solve all subsystems in parallel
print("\tSolve algebraic systems in parallel.")
print("\t\tSolving sub-systems in parallel:")
# We do the fitting in parallel (to use all cores, use mp.cpu_count())
pool = mp.Pool(mp.cpu_count())
# Solve the equations in parallel instead
generator_strings = pool.starmap(solve_algebraic_system,[(coeff_vec_temp,c_LV,c_subsets_LV[index],eta_LV,eq_sys_LV,poly_LV,states_LV,index+1,len(coeff_vec_LV),omega_LV) for index,coeff_vec_temp in enumerate(coeff_vec_LV)]) 
# Close the pool
pool.close()
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("\t\t\tDone!\tCurrent Time =", current_time)
# Write the solutions to file
#Open text file
print("\tWrite the answer to a tex-file.")
text_file = open("../Output/generators_LV_competitive.tex", "w")
#text_file.write(solution_string)
# Loop over the output and write them all to our file
for generator_string in generator_strings:
    #write string to file
    text_file.write(generator_string)
print("=======================================================================\n")
# Close file
text_file.close()
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("\t\t\tDone!\tCurrent Time =", current_time)
# =================================================================================
# =================================================================================
# The BZ-model
# =================================================================================
# =================================================================================


# =================================================================================
# =================================================================================
# The Brusselator
# =================================================================================
# =================================================================================
