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
# The script starts
# =================================================================================
# =================================================================================
# Include our home-made symmetry toolbox
include("symmetry_toolbox_first_order_ODEs.jl")
# Define our symbolic variables
@variables t u(t) v(t)
# Save all variables
variables = [t, u, v]
# Print the states
states = [u, v]
println("The states are:")
for state in states
    println(string.(state))
end
# Create the tangential ansatze
degree = 2
monomials, tangential_ansatze, unknown_coefficients = create_tangential_ansatze(degree,variables)
print("\t\tMONOMIALS:\n",monomials,"\n\n")
print("\t\tTANGENTIAL ANSATZE:\n",tangential_ansatze,"\n\n")
print(tangential_ansatze[1])
print("\t\tUNKNOWN COEFFICIENTS:\n",unknown_coefficients,"\n\n")


