#=
================================================================================
=================================================================================
 Script:"conduct_symmetry_calculations"
 Date: 2021-07-07
 Implemented by: Johannes Borgqvist
 Description:
 The script conducts all the symmetry calculations using functions that are
 stored in the script "symmetry_toolbox_first_order_ODEs.jl". This script
 contains merely one function called "calculate_symmetries_ODEs" which reads a
 model from a given data file and attempts to calculate the symmetries.
=================================================================================
=================================================================================
 Import Libraries
=================================================================================
=================================================================================
=#
include("read_and_write_data.jl")
#=
=================================================================================
=================================================================================
 Run the code
=================================================================================
=================================================================================
=#
# Give the file name
file_name = "../Input/ODE_models.xlsx"
# Give the name of the sheet
sheet_number = 1
# Allocate some symbolic variables
variables_symbolic = []
parameters_symbolic = []
reaction_terms_symbolic = []
omega_list = []
#model_name = ""
# Define our lovely model
x_sym,model_name = read_input_model(file_name,sheet_number,variables_symbolic, parameters_symbolic,reaction_terms_symbolic,omega_list)
# Print the variables
print("\t\tThe variables are:\n")
for (index, value) in enumerate(variables_symbolic)
    print("\t\t\t$value\n")
end
# Print the parameters
print("\t\tThe parameters are:\n")
for (index, value) in enumerate(parameters_symbolic)
    print("\t\t\t$value\n")
end
# Print the reaction_terms
print("\t\tThe reaction terms are:\n")
for (index, value) in enumerate(reaction_terms_symbolic)
    print("\t\t\t$value\n")
end
# Print the general variables
print("\t\tThe general variables are:\n")
for (index, value) in enumerate(x)
    print("\t\t\t$value\n")
end
# Print the general reaction terms
print("\t\tThe general reaction terms are:\n")
for (index, value) in enumerate(omega_list)
    print("\t\t\t$value\n")
end
# Print the model name
print("\t\tThe name of the model is:\t$model_name")



