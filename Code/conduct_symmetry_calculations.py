#=================================================================================
#=================================================================================
# Script:"conduct_symmetry_calculations"
# Date: 2021-06-01
# Implemented by: Johannes Borgqvist
# Description:
# The script conducts all the symmetry calculations using functions that are
# stored in the script "symmetry_toolbox_first_order_ODEs.py". This script
# contains merely one function called "calculate_symmetries_ODEs" which reads a
# model from a given data file and attempts to calculate the symmetries.
#=================================================================================
#=================================================================================
# Import Libraries
#=================================================================================
#=================================================================================
import read_and_write_data # Home-made
import symmetry_toolbox_first_order_ODEs # Home-made
# For printing with latex
from sympy import latex
# To time each part of the program
import time
#=================================================================================
#=================================================================================
# The Function
#=================================================================================
#=================================================================================
# The script contains only one function called "calculate_symmetries_ODEs".
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# FUNCTION 1: "calculate_symmetries_ODEs"
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# The function takes two inputs:
# 1. The string "file_name" with the name of the model file,
# 2. The degree of the polynomials in the tangent ansätze called "tangent degree".
# The function stores the calculated generator and all introduced variables in a folder named after the string "file_name" in the output folder with path "../Output/".
def calculate_symmetries_ODEs(file_name,tangent_degree):
    t0_total = time.time()
    # Define the file name of the model
    #file_name = "hydons_model"
    # Print to the user that we are creating ansätze
    print("\t\t\tStep 1 out of 6: Loading the model...")
    # Read the model given by the file name
    model_reading_successful, variables, parameters, reaction_terms, omega_list = read_and_write_data.read_input_model(file_name)
    if model_reading_successful:
        print("\t\t\t\tThe model was successfully loaded!\n\t\t\t\tDone!\n")
        # Calculate the number of variables and number of states
        num_of_variables = 1 # Number of variables 
        num_of_states = len(variables)-num_of_variables # number of states
        # Print to the user that we are creating ansätze
        print("\t\t\tStep 2 out of 6: Creating tangent ansätze...")
        # Calculate our new tangents and the variables
        t0 = time.time()
        x, c, eta_list = symmetry_toolbox_first_order_ODEs.create_tangent_ansatze(num_of_variables,num_of_states,tangent_degree)
        t1 = time.time()
        print("\t\t\tTime elapsed:\t%0.5f\tseconds"%(t1-t0))
        # Print that this is done
        print("\t\t\t\tDone!")
        # Print to the user that we are calculating the linearised symmetry conditions
        print("\t\t\tStep 3 out of 6: Calculating the linearised symmetry conditions...")
        # Calculate the linearised symmetry conditions
        t0 = time.time()
        lin_sym_list = symmetry_toolbox_first_order_ODEs.lin_sym_cond(x,eta_list,omega_list)
        t1 = time.time()
        print("\t\t\tTime elapsed:\t%0.5f\tseconds"%(t1-t0))        
        # Print that this is done
        print("\t\t\t\tDone!")
        # Print to the user that we are calculating the determining equations
        print("\t\t\tStep 4 out of 6: Deriving the determining equations...")
        t0 = time.time()
        # Define the degree we want on our monomial
        degree_monomial = 7
        # Step 4.3: Calculate our so called "determining equations"
        det_eq, monomials, lin_sym_eq_number = symmetry_toolbox_first_order_ODEs.determining_equations(x,lin_sym_list,degree_monomial)
        t1 = time.time()
        print("\t\t\tTime elapsed:\t%0.5f\tseconds"%(t1-t0))        
        # Print that this is done
        print("\t\t\t\tDone!")
        # Print to the user that we are solving the determining equations
        print("\t\t\tStep 5 out of 6: Solving the determining equations...")
        t0 = time.time()
        X = symmetry_toolbox_first_order_ODEs.solve_determining_equations(x,eta_list,c,det_eq,variables)
        t1 = time.time()
        print("\t\t\tTime elapsed:\t%0.5f\tseconds"%(t1-t0))        
        # Change the name of the variables
        for i in range(len(variables)):
            from_str = str(latex(x[i]))
            to_str = str(latex(variables[i]))
            X = X.replace(from_str,to_str)
        # Print that this is done
        print("\t\t\t\tDone!")
        # Print to the user that we are solving the determining equations
        print("\t\t\tStep 6 out of 6: Saving the data...")    
        # Print that we are saving the data
        read_and_write_data.write_output_generator(tangent_degree,file_name,variables,x,X,reaction_terms,omega_list)
        print("\t\t\t\tDone!\n")
        t1_total = time.time()
        print("\t\t\tTotal time elapsed:\t%0.5f\tseconds\n"%(t1_total-t0_total))
        print("\t\t\tThe calculations are finished.")
        
    else:# We could not read the model: => Print an error message!
        print("\t\t\tError when reading the model!\n\t\t\tAborting...\n\n")
