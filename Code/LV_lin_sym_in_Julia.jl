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
using Pkg
#Pkg.add("Symbolics")
using Symbolics
# =================================================================================
# =================================================================================
# The script starts
# =================================================================================
# =================================================================================
print("\n\t\tHello world!\n\n")


