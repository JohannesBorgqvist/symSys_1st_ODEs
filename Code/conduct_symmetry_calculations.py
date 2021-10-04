# =================================================================================
# =================================================================================
# Script:"conduct_symmetry_calculations"
# Date: 2021-06-01
# Implemented by: Johannes Borgqvist
# Description:
# The script conducts all the symmetry calculations using functions that are
# stored in the script "symmetry_toolbox_first_order_ODEs.py". This script
# contains merely one function called "calculate_symmetries_ODEs" which reads a
# model from a given data file and attempts to calculate the symmetries.
# =================================================================================
# =================================================================================
# Import Libraries
# =================================================================================
# =================================================================================
import read_and_write_data  # Home-made
import symmetry_toolbox_first_order_ODEs  # Home-made
import sympy
# For printing with latex
from sympy import latex
# To time each part of the program
import time
# Import matrices
from sympy import Matrix

# =================================================================================
# =================================================================================
# The Function
# =================================================================================
# =================================================================================
# The script contains only one function called "calculate_symmetries_ODEs".
# ----------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
# FUNCTION 1: "calculate_symmetries_ODEs"
# ----------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
# The function takes two inputs:
# 1. The string "file_name" with the name of the model file,
# 2. The degree of the polynomials in the tangent ansätze called "tangent degree".
# The function stores the calculated generator and all introduced variables in a folder named after the string "file_name" in the output folder with path "../Output/".


def calculate_symmetries_ODEs(file_name, tangent_degree):
    t0_total = time.time()
    # Define the file name of the model
    #file_name = "hydons_model"
    # Print to the user that we are creating ansätze
    print("\t\t\tStep 1 out of 6: Loading the model...")
    # Read the model given by the file name
    model_reading_successful, variables, parameters, reaction_terms, omega_list = read_and_write_data.read_input_model(
        file_name)
    if model_reading_successful:
        print("\t\t\t\tThe model was successfully loaded!\n\t\t\t\tDone!")
        # Calculate the number of variables and number of states
        num_of_variables = 1  # Number of variables
        num_of_states = len(variables)-num_of_variables  # number of states
        # Print to the user that we are creating ansätze
        print("\t\t\tStep 2 out of 6: Creating tangent ansätze...")
        # Calculate our new tangents and the variables
        x, c, eta_list = symmetry_toolbox_first_order_ODEs.create_tangent_ansatze(
            num_of_variables, num_of_states, tangent_degree)
        #print("# Defining the tangents")
        eta_temp = eta_list
        #print("The tangents are:<br>")
        #print("\\begin{align*}")
        #tangent_counter = 1
        #for eta in eta_temp:
            #for index in range(len(x)):
                #eta = eta.subs(x[index],variables[index])
            #if tangent_counter == 1:
                #print("\\xi&=%s\\\\"%(str(latex(eta))))
            #else:
                #print("\\eta_%d&=%s\\\\"%(tangent_counter-1,str(latex(eta))))                
                #print(eta.atoms(sympy.Mul))
            #print(eta.args)
            #tangent_counter+=1
        #print("\\end{align*}")
        #print("The unknown coefficients:<br>")
        #print("\\begin{equation}\n\\mathbf{c}=%s\n\\end{equation}"%(latex(Matrix(len(c),1,c))))
        # Print that this is done
        print("\t\t\t\tDone!")
        # Print to the user that we are calculating the linearised symmetry conditions
        print("\t\t\tStep 3 out of 6: Calculating the linearised symmetry conditions...")
        # Calculate the linearised symmetry conditions
        lin_sym_list = symmetry_toolbox_first_order_ODEs.lin_sym_cond(
            x, eta_list, omega_list)
        #print("# The linearised symmetry conditions")
        #lin_sym_counter = 1
        #for lin_sym in lin_sym_list:
            #print("Linearised symmetry condition\t%d:<br>"%(lin_sym_counter))
            #print("\\begin{align*}")
            #temp_str = "0&="
            #help_counter = 1
            #for index in range(len(lin_sym.args)):
                #term = lin_sym.args[index]
                #for index in range(len(x)):
                    #term = term.subs(x[index],variables[index]) 
                #temp_str += str("\\left(" + latex(term) + "\\right)")
                #if help_counter<6:
                    #temp_str += "+"
                #else:
                    #temp_str += "\\\\"
                    #print("%s"%(temp_str))
                    #help_counter = 1
                    #temp_str = "&+"
                #help_counter += 1
                
            #print("\\end{align*}")
            #lin_sym_counter+=1        
        # Print that this is done
        print("\t\t\t\tDone!")
        # Print to the user that we are calculating the determining equations
        print("\t\t\tStep 4 out of 6: Deriving the determining equations...")
        # Define the degree we want on our monomial
        degree_monomial = 7
        # Step 4.3: Calculate our so called "determining equations"
        det_eq, monomials, lin_sym_eq_number = symmetry_toolbox_first_order_ODEs.determining_equations(
            x, lin_sym_list, degree_monomial)
        
        #print(lin_sym_eq_number)
        #print("# The determining equations")

        #print("Determining equations from linearised symmetry condition\t%d:<br>"%(1))
        #print("\\begin{align*}")
        #condition_counter = 2
        #for det_eq_index in range(len(det_eq)):
            #mon = monomials[det_eq_index]
            #temp_mon = ""
            #for var_index in range(len(variables)-1):
                #power = mon[var_index]
                #variable = variables[var_index+1]
                #temp_mon += latex(variable) + "^{" + latex(power) + "}"
            #temp_det = det_eq[det_eq_index]
            #for index in range(len(x)):
                #temp_det = temp_det.subs(x[index],variables[index])            
            #temp_str = temp_mon + ": 0&=" + latex(temp_det) + "\\\\"                        
            #print("%s"%(temp_str))
            #temp_str = ""
            #if det_eq_index >0:
                #if det_eq_index==(len(det_eq)-1):
                    #print("\\end{align*}")                    
                #elif lin_sym_eq_number[det_eq_index]<lin_sym_eq_number[det_eq_index+1]:
                    #print("\\end{align*}")
                    #print("Determining equations from linearised symmetry condition\t%d:<br>"%(condition_counter))
                    #condition_counter += 1
                    #print("\\begin{align*}")
                
                

            
        #print("\\end{align*}")
        #lin_sym_counter+=1
        
        # Print that this is done
        print("\t\t\t\tDone!")
        # Print to the user that we are solving the determining equations
        print("\t\t\tStep 5 out of 6: Solving the determining equations...")
        X = symmetry_toolbox_first_order_ODEs.solve_determining_equations(
            x, eta_list, c, det_eq, variables,omega_list)
        # Change the name of the variables
        for i in range(len(variables)):
            from_str = str(latex(x[i]))
            to_str = str(latex(variables[i]))
            X = X.replace(from_str, to_str)
        # Print that this is done
        print("\t\t\t\tDone!")
        # Print to the user that we are solving the determining equations
        print("\t\t\tStep 6 out of 6: Saving the data...")
        # Print that we are saving the data
        read_and_write_data.write_output_generator(
            tangent_degree, file_name, variables, x, X, reaction_terms, omega_list)
        print("\t\t\t\tDone!\n")
        print("\t\t\tThe calculations are finished.")

    else:  # We could not read the model: => Print an error message!
        print("\t\t\tError when reading the model!\n\t\t\tAborting...\n\n")
