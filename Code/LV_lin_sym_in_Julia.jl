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
#Pkg.add("TypedPolynomials")
#using TypedPolynomials
#Pkg.add("DynamicPolynomials")
using DynamicPolynomials
#Pkg.add("MultivariateBases")
using MultivariateBases
#Pkg.add("SymPoly")
#using SymPoly
#Pkg.add("SymbolicUtils")
#using SymbolicUtils
#Pkg.add("Symbolics")
using Symbolics
#Pkg.add("Latexify")
using Latexify
#Pkg.add("Combinatorics")
#Pkg.add("BenchmarkTools")
#using Combinatorics, BenchmarkTools
#using Combinatorics
# =================================================================================
# =================================================================================
# Functions
# =================================================================================
# =================================================================================
# Function 1: allocate_monomials
# This functions takes two inputs:
# 1. An integer named degree which determines the degree of the polynomial,
# 2. A vector of variables containing our symbolic dependent and independent variables.
# The function returns one output:
# 1. A vector basis_vec_states containing all the basis monomials in terms of the states.
function allocate_monomials(degree::Int64,variables::Vector{Num})
    # The number of variables
    num_variables = length(variables)
    # Define a polynomial array with our variables    
    eval(Meta.parse(string.("@polyvar x[1:",string.(num_variables),"]")))
    # Define a list of all our monomials in our multivariate polynomial
    basis_vec = MultivariateBases.monomials(x, 0:degree)
    # Extract the exponents for each monomial as a tuple
    exponents_of_basis = [exponents(monomial) for monomial in basis_vec]
    # Allocate memory for basis functions in terms of the states
    basis_vec_states = []
    # Loop over all monomials in terms of their exponents
    for (outer_index,exponent) in enumerate(exponents_of_basis)
        # Allocate a temporary monomial
        temp_monomial = 1
        # Loop over all powers in our current exponent and construct our new basis monomial
        for (inner_index,power) in enumerate(exponent)
            temp_monomial *= variables[inner_index]^power
        end
        # Append our new basis function to our list of basis functions in terms of the states
        append!(basis_vec_states,temp_monomial)
    end
    return basis_vec_states
end
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
# Create the tangential ansatze
monomials = allocate_monomials(2,variables)


