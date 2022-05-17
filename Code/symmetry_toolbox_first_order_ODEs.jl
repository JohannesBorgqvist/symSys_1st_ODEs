# =================================================================================
# =================================================================================
# Script:"symmetry_toolbox_first_order_ODEs.jl"
# Date: 2022-05-17
# Implemented by: Johannes Borgqvist
# Description:
# This script contains the functions regarding symbolic manipulations in julia.
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
    # Return out new basis vector 
    return basis_vec_states
end
# Function 2: create_tangential_ansatze
# This functions takes two inputs:
# 1. An integer named degree which determines the degree of the polynomial,
# 2. A vector of variables containing our symbolic dependent and independent variables.
# The function returns two outputs:
# 1. A vector monomials containing all our monomials in terms of the states,
# 2. A vector tangential_ansatze containing our allocated tangents,
# 3. A vector unknown_coefficients containing the unknown coefficients in the tangential anstaze.
function create_tangential_ansatze(degree::Int64,variables::Vector{Num})
    # Allocate memory for our tangential ansatze
    tangential_ansatze = []
    # Calculate our monomials
    monomials = allocate_monomials(degree,variables)
    # Calculate the number of monomials
    num_monomials = length(monomials)        
    # Allocate memory for our unknown coefficients
    unknown_coefficients = [eval(Meta.parse(string.("@variables c_",string(i)))) for i in 1:(3*length(monomials))]
    # Now, we generate a few polynomials that will be used to create the tangential ansatze in the end
    polynomials = []
    # Loop through the number of states
    for pol_index in 1:length(variables)
        # Allocate a temporary polynomial
        @variables temp_pol
        temp_pol = 0
        # Loop through all monomials and add them to our polynomial
        for (i,monomial) in enumerate(monomials)
            coefficient_index = i+((pol_index-1)*num_monomials)
            temp_term = unknown_coefficients[coefficient_index]*monomial
            temp_pol += temp_term[1]
        end
        # Append the polynomial to our vector of polynomials
        append!(polynomials,temp_pol)
    end
    # Loop over the tangents and save them in our lovely list of tangents
    for (i,variable) in enumerate(variables)
        if i ==1
            @variables xi(t)
            append!(tangential_ansatze,xi)
        else
            append!(tangential_ansatze,polynomials[i]*exp(polynomials[1]))            
        end
    end
    # Return our monomials, tangential ansatze and unknown coefficients
    return monomials, tangential_ansatze, unknown_coefficients 
end
