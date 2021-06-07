#=================================================================================
#=================================================================================
# Script:"read_and_write_data"
# Date: 2021-05-18
# Implemented by: Johannes Borgqvist
# Description:
# The script contains two functions which
# reads a certain model as input from
# a pre-defined file and then it writes
# the results of the symmetry calculations
# to a folder which it also creates. It
# also contains a help function which
# translates the names of all variables
# from the ones used during the symmetry
# calculations (e.g. (x_0,x_1,...)) to
# the variable names provided in the
# data file. 
#=================================================================================
#=================================================================================
# Import Libraries
#=================================================================================
#=================================================================================
# For reading the data frame
# using pandas
import pandas as pd
# For symbolic calculations
from sympy import *
# To create a data folder
import os
# To get the date and time
import datetime
# To read and write binary files
# in order to store all equations
# from all runs
import pickle
#=================================================================================
#=================================================================================
# The Functions
#=================================================================================
#=================================================================================
# The following functions are defined subsequently:
# 1. *read_input_model* (reads the model and creates output folder)
# 2. *write_output_generator* (write the output of the symmetry calculations to files)
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
# FUNCTION 1: "read_input_model"
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
# The function takes the name of the model file (which is a string) as input and
# then it looks for an xlsx-files with that name stored in the folder "../Input".
# It has four outputs: a logical variable telling whether or not the reading of the model was successful, the list "variables" which contains the original variables
# as defined in the model file, the list "parameters" containing the parameters
#  of the model and the list "omega_list" which contains the reaction terms
# defined in the model file but translated to the generic variables (x0,x1,x2,...)
# where by variables we mean both the independent and dependent variables. The
# script does also create an output directory in the folder "../Output" which
# is later used for storing the data files.
def read_input_model(file_name):
    # Allocate memory for the lists that we will return if we are successful
    variables = [] # All the dependent and independent variables
    parameters = [] # All the parameters of the model
    reaction_terms = [] # All the reaction terms
    omega_list = [] # All the reaction terms written in generic variables (i.e. (x0,x1,x2,...))
    # Define a logical variable telling us if we could read the data file
    model_reading_successful = False
    # Check if the file exists
    file_str = "../Input/" + file_name + ".xlsx"
    print(file_str)
    #  Read the model
    if os.path.isfile(file_str):
        # Store the model in a data frame df
        df = pd.read_excel(file_str)
        # Check the categories in our model file
        # iterating the columns
        list_of_columns = [col for col in df.columns]
        # Columns in a correctly written data file
        correct_columns = ["Model name", "Variable", "States", "Parameters", "Reaction terms"]
        # See if the data is correctly provided
        if list_of_columns == correct_columns:
            # Read the name of the model
            model_name = str(df["Model name"][0])  
            # Check if a name has been provided.
            if len(model_name)== 0 or model_name == "nan": # No name given
                print("\t\tERROR: No model name has been given. Add a name of the model to the model file. ")
            else:# A name has been given
                # VARIABLES
                # Save the independent variable
                if str(df["Variable"][0]) != "nan":
                    variables.append(df["Variable"][0])
                # Save the dependent variables
                for var in df["States"]:
                    if str(var) != "nan": 
                        variables.append(var)
                # REACTION TERMS
                reaction_terms = [r for r in df["Reaction terms"] if str(r)!="nan"]
                # PARAMETERS
                parameters = [r for r in df["Parameters"] if str(r)!="nan"]         
                # Create a of the reaction terms list which the program can process
                omega_list = reaction_terms.copy()
                # Loop over the reaction terms and replace with the appropriate variable names
                for var_index in range(len(variables)):
                    # Loop over the reaction terms
                    for reac_index in range(len(omega_list)):
                        # Do the substitution
                        exec("omega_list[reac_index] = omega_list[reac_index].replace(variables[var_index],'x_%d')"%(int(var_index)))
                # Render these lists as symbolic objects
                reaction_terms = [sympify(r) for r in reaction_terms]
                omega_list = [sympify(r) for r in omega_list]
                variables = [sympify(r) for r in variables]
                if len(parameters) != 0:
                    parameters = [sympify(r) for r in parameters]                
                # Check that we have the correct number of ODEs and that
                # all parameters and states that are introduced are the only
                # thing that occurs in the reaction terms
                list_var_par = variables + parameters # Parameters and variables
                list_reaction = [] # Reaction terms
                # Loop over reaction terms and find all symbols in them
                for i in range(len(reaction_terms)):
                    list_reaction += [a for a in reaction_terms[i].atoms(Symbol)]
                    list_reaction += [a for a in reaction_terms[i].atoms(Function)]                     
                print("Variables+Parameters:")
                print(list_var_par)
                print("Reaction terms:")
                print(list_reaction)                
                # See if both these guys contains the same symbols and that the correct number of ODEs is provided.                    
                if (set(list_var_par)==set(list_reaction)) and (len(reaction_terms) == (len(variables)-1) ):
                    model_reading_successful = True                       
                else:
                    print("\t\tERROR: The reaction terms and/or the states are given with the wrong format. Please check that the number of ODEs matches the number of states and that only defined symbols (parameters and variables) occurs in the reaction terms.")                        
        else: # The wrong columns have been provided
            print("\t\tERROR: Please re-write the model file with the following columns (CASE SENSITIVE):")
            print(correct_columns)
    else: # The model not correctly provided
        print("\tERROR:File does not exist!")
        print("\nNote that only the name of the file (e.g. model1) that is stored in the folder ../Input should be provided and no path should be provided. Thus both ../Input/model1.xlsx and model1.xlsx will give the wrong answer.")
    # Return the output
    return model_reading_successful, variables, parameters, reaction_terms, omega_list
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
# FUNCTION 2: "write_output_generator"
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
# The function takes various lists of equations, tangents, the variables,
# the parameters, the algebraic equations, the solution of the algebraic equations and it saves all these lists in the pickle-format (i.e. binary files).
# The files are stored in a subsub-folder of the folder "../Output" where the first folder is determined by the name in the input file and the second folder is a unique folder named after the time and date that the script was launched. Also, the generatoris written into the Markdown-file "out.md" stored in the sub folder of Output named after the name of the model which is defined in the model file. 
def write_output_generator(tangent_degree,folder_name,variables,parameters,reaction_terms,x,omega_list,eta_list,lin_sym_list,det_eq,monomials,lin_sym_eq_number,X,c_original,c,eq_alg,sol_alg):
    # Begin by creating our dear output folder
    path = "../Output/" + folder_name
    os.makedirs(path, exist_ok=True)
    # Create a sub folder in which we store all lists.
    # Each of these folders is named after the current date and time
    date_str = datetime.datetime.now().strftime("%I_%M%p_%d_%B-%Y")
    sub_folder = path + "/" + date_str
    os.makedirs(sub_folder, exist_ok=True)
    #-----------------------------------------------------------------------
    # Part 1: Save all lists using pickle 
    #-----------------------------------------------------------------------
    # Save the states and the independent variable
    temp_str = sub_folder + "/variables.pickle"
    with open(temp_str, 'wb') as outf:
        outf.write(pickle.dumps(variables))
    # Save the parameters
    temp_str = sub_folder + "/parameters.pickle"
    with open(temp_str, 'wb') as outf:
        outf.write(pickle.dumps(parameters))
    # Save the reaction terms
    temp_str = sub_folder + "/reaction_terms.pickle"
    with open(temp_str, 'wb') as outf:
        outf.write(pickle.dumps(reaction_terms))        
    # Save the re-named variables
    temp_str = sub_folder + "/x.pickle"
    with open(temp_str, 'wb') as outf:
        outf.write(pickle.dumps(x))
    # Save the re-named reaction terms
    temp_str = sub_folder + "/omega_list.pickle"
    with open(temp_str, 'wb') as outf:
        outf.write(pickle.dumps(omega_list))
    # Save the tangents
    eta_list = [str(eta_list[i]) for i in range(len(eta_list))]
    temp_str = sub_folder + "/eta_list.pickle"
    with open(temp_str, 'wb') as outf:
        outf.write(pickle.dumps(eta_list))
    # Save the lists of the linearised symmetry condition
    lin_sym_list = [str(lin_sym_list[i]) for i in range(len(lin_sym_list))]
    temp_str = sub_folder + "/lin_sym_list.pickle"
    with open(temp_str, 'wb') as outf:
        outf.write(pickle.dumps(lin_sym_list))        
    # Save the determining equations
    det_eq = [str(det_eq[i]) for i in range(len(det_eq))]
    temp_str = sub_folder + "/det_eq.pickle"
    with open(temp_str, 'wb') as outf:
        outf.write(pickle.dumps(det_eq))
    # Save the monomials
    monomials = [str(monomials[i]) for i in range(len(monomials))]
    temp_str = sub_folder + "/monomials.pickle"
    with open(temp_str, 'wb') as outf:
        outf.write(pickle.dumps(monomials))
    # Save the list connecting the determining equations
    # to the corresponding linearised symmetry condition
    temp_str = sub_folder + "/lin_sym_eq_number.pickle"
    with open(temp_str, 'wb') as outf:
        outf.write(pickle.dumps(lin_sym_eq_number))
    # Save the original coordinates
    c_original = [str(c_original[i]) for i in range(len(c_original))]
    temp_str = sub_folder + "/c_original.pickle"
    with open(temp_str, 'wb') as outf:
        outf.write(pickle.dumps(c_original))        
    # Save the coordinates before the algebraic substitution
    c = [str(c[i]) for i in range(len(c))]
    temp_str = sub_folder + "/c.pickle"
    with open(temp_str, 'wb') as outf:
        outf.write(pickle.dumps(c))
    # Save the algebraic equations
    eq_alg = [str(eq_alg[i]) for i in range(len(eq_alg))]    
    temp_str = sub_folder + "/eq_alg.pickle"
    with open(temp_str, 'wb') as outf:
        outf.write(pickle.dumps(eq_alg))
    # Save the solutions to the algebraic equations
    sol_alg = [str(sol_alg[i]) for i in range(len(sol_alg))]
    temp_str = sub_folder + "/sol_alg.pickle"
    with open(temp_str, 'wb') as outf:
        outf.write(pickle.dumps(sol_alg))        
    #-----------------------------------------------------------------------
    # Part 2: Write the output to a Markdown file
    #-----------------------------------------------------------------------
    # Define the file which we write to
    file_path = path + "/out.tex"
    # Importantly append to the file
    f = open(file_path, "a")
    # Write the date as the title
    date_str = date_str.replace("_","\\_")
    f.write("\\section*{Run %s}\n" % (date_str))
    # Print the degree of the polynomial in the ansätze
    f.write("\\textbf{Degree in tangential ansätze:}\t%d\\\\\n"%(int(tangent_degree)))
    # Print the system of ODEs
    f.write("The reaction terms are:\n\n")
    f.write("\\begin{align*}\n")
    for i in range(len(omega_list)):
        f.write('\\dv{%s}{%s}&=%s\\\\\n'%(latex(variables[i+1]),latex(variables[0]),latex(reaction_terms[i])))
    f.write("\\end{align*}\n\n")
    # Change name of the generator
    for i in range(len(x)):
        exec("X = X.replace('x_{%d}',str(latex(variables[%d])))"%(int(i),int(i)))
    # Print the generator
    f.write("The calculated generator is:\n\n")
    f.write("%s\n"%(X))
    # Close the file
    f.close()    
