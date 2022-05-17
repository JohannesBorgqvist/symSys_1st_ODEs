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
# Define our growth rate a
@variables a
# Define our reaction terms
omega_list = [u*(1-v), a*v*(u-1)]
# Create the tangential ansatze
degree = 2
monomials, tangential_ansatze, unknown_coefficients, polynomials = create_tangential_ansatze(degree,variables)
# Formulate the linearised symmetry condition
lin_syms = lin_sym_conds_fibre_preserving(degree,variables,tangential_ansatze,omega_list,polynomials)
# Write things to a latex document
str_1 = latexify(lin_syms[1])
str_1 = replace(str_1, "\\begin{equation}" => "\\begin{align*}")
str_1 = replace(str_1, "\n\\end{equation}" => "&=0\\\\")
str_2 = latexify(lin_syms[2])
str_2 = replace(str_2,"\\begin{equation}\n"=>"")
str_2 = replace(str_2, "\n\\end{equation}" => "&=0\\\\\n\\end{align*}")
tot_str = string(str_1,str_2)
write("latex_output_julia/lin_syms.tex",tot_str)



