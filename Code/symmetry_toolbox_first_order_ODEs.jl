#=
================================================================================
=================================================================================
 Script:"symmetry_toolbox_first_order_ODEs"
 Date: 2021-07-12
 Implemented by: Johannes Borgqvist
 Description:
 The script contains all the function that are necessary in order to
 calculate a symmetry generator for a given system of ODEs. There is a function
 called "create_tangent_ansatze" which sets up the so called projective tangent
 ansätze that are polynomials in the states where each coefficient is a function
 of the independent variable. The function "lin_sym_cond" formulates the
 linearised symmetry condition, and given these the function "determining_equations" formulates the determining equations. Lastly, the major function which does
 most of the work is called "solve_determining_equations"which solves the
 determining equations for the given system of ODEs. 
=================================================================================
=================================================================================
 Import Libraries
=================================================================================
=================================================================================
=#
using Symbolics # For doing the calculations
using ModelingToolkit # For some of the variables
using DynamicPolynomials # For generating polynomials
using SymbolicUtils # To do symbolic calculations
using SymPoly # To convert a polynomic variable to a symbolic one
#=
=================================================================================
=================================================================================
 Functions
=================================================================================
=================================================================================
=#
#=
 Function 1: "create_tangent_ansatze"
The function takes the number of variables,
the symbolic array with the variables and
the degree of the mutlivariate polynomial that
will constitute the ansatz. It then returns a
list of the polynomial ansätze eta_list as well
as a list of the coefficients in the ansätze
denoted by c. 
=#
function create_tangent_ansatze(x,num_of_variables,degree_polynomial)
    # Allocate two lists
    c = []
    eta_list = []
    # Calculate the number of terms
    num_terms = factorial(degree_polynomial+(size(x,1)-num_of_variables),degree_polynomial)
    print("\n\t\tThe number of terms is:\t$num_terms\n")
    @polyvar y[1:size(x,1)]
    # Print the variables
    print("\t\tOur variables are:\n")
    for (index, value) in enumerate(y)
        print("\t\t\t$value\n")
    end
    # Generate some monomials
    list_mon = monomials(y, 1:degree_polynomial)
    # Print
    print("\t\tThe list of polynomials are:\t$list_mon\n")
    # Create a polynomial
    eq = x[1]+x[2]+x[3]
    print("\t\tCreated an arbitrary polynomial:\t$eq\n")
    # Create polynomial
    #dict = Dict(x[1] => x[1], x[2] => y[2], x[3] => y[3])
    #p = poly(eq, x=>y)
    print("\t\tConverted it to a symbolic polynomial:\t$p")
    # Return the output 
    return c, eta_list
end
