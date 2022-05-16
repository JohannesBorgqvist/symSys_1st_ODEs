# =================================================================================
# =================================================================================
# Script:"LV_lin_sym_in_Julia.jl"
# Date: 2022-05-13
# Implemented by: Johannes Borgqvist
# Description:
# The script sets up the linearised symmetry conditions as well as various
# tangential ansatze and then it tries to solve the resulting algebraic equations.
# We have previously made all these scripts in Python, but for performance reasons
# we will now try to implement this in Julia.
# =================================================================================
# =================================================================================
# Install and add packages
# =================================================================================
# =================================================================================
#using Pkg
#Pkg.add("Symbolics")
using Symbolics
#Pkg.add("Latexify")
using Latexify
# =================================================================================
# =================================================================================
# The script starts
# =================================================================================
# =================================================================================
# Define our symbolic variables
@variables t u(t) v(t)
# Save all variables
variables = [t, u, v]
# Allocate memory for our tangents
tangents = []
# Loop over our three tangents
for (i,variable) in enumerate(variables)
    if i ==1
        eval(Meta.parse("@variables xi(t)"))
    else
        temp_str_1 = string.("@variables eta_",string.(i-1),"(")
        temp_str_2 = string.("eta_",string.(i-1))
        # Loop over all variables
        for (sub_index,sub_variable) in enumerate(variables)
            # We want to define things as functions of u rather than u(t)
            variable_str = replace(string.(sub_variable),"(t)" => "")
            # Parenthesis at final iteration and comma otherwise
            if sub_index == (length(variables))
                temp_str_1 = string.(temp_str_1,variable_str,")")
            else
                temp_str_1 = string.(temp_str_1,variable_str,",")                
            end
        end
        # Allocate the tangent as a symbolic function
        eval(Meta.parse(temp_str_1))
        # Append our lovely tangent to a list
        eval(Meta.parse(string.("append!(tangents,",temp_str_2,")")))
    end
end
# Print the states
states = [u, v]
println("The states are:")
for state in states
    println(string.(state))
end
# Print the tangents
println("The tangents are:\n")
for tangent in tangents
    println(tangent)
end



