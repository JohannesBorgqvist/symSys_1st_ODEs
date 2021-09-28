from sympy import *
#from sympy.core.function import AppliedUndef
from sympy.core.function import *


def identify_basis_functions(expr,coeff_list):
    # Define the arbitrary functions
    arbitrary_functions = list(expr.atoms(AppliedUndef))
    # Allocate memory an empty list for the basis functions
    basis_functions = []
    # Define a list with all arguments
    arguments = list(expr.args)
    # Loop through all arguments
    for argument in arguments:
        # Loop through all coefficients 
        for coefficient in coeff_list:
            # Save all basis functions provided that
            # the exist as a linear combination of the
            # coefficients in coeff_list
            if argument.coeff(coefficient) != 0: 
                basis_functions.append(argument.coeff(coefficient))
    # Lastly, we only want unique basis functions so we take the
    # set of these functions
    basis_functions = list(set(basis_functions).difference(set(arbitrary_functions)))
    # We also need to remove doublettes. Start with the constants
    constant_basis_list = []
    for base_index in range(len(basis_functions)):
        if len(basis_functions[base_index].atoms(Symbol))==0:
            constant_basis_list.append(base_index)
    # Make sure that the constant correspond to the basis element 1
    basis_functions[constant_basis_list[0]] = basis_functions[constant_basis_list[0]]/basis_functions[constant_basis_list[0]]
    # Remove the zeroth element
    del constant_basis_list[0]
    # Sort in reverse order
    constant_basis_list.sort(reverse=True)
    # Remove all these constants from the basis list
    for index in constant_basis_list:
        del basis_functions[index]
        
    # Lastly, save only unique values
    unique_values = []
    # Save the first element
    unique_values.append(basis_functions[0])
    # Loop over all basis functions and save the unique ones
    for base in basis_functions:
        # Skip the first one
        if base != basis_functions[0]:
            # Allocate a temp sum
            temp_sum = 0
            # Loop over all unique values and see if we append
            # the value or not
            for unique_value in unique_values:
                # We only care about non-constant functions
                if len(unique_value.atoms(Symbol))!=0:
                    if len(simplify(base/unique_value).atoms(Symbol))==0:
                        # Add term to the temp sum
                        temp_sum += 1
            # If no terms were added then we have a unique value
            if temp_sum == 0:
                divide_by_this = factor_list(base)[0]
                unique_values.append(simplify(base/divide_by_this))
    # Finally, we assign the unique values to the basis functions
    basis_functions = unique_values        
    # Return the output
    return arbitrary_functions, basis_functions

def solve_algebraic_equation(expr,coeff_list):
    # Find all basis functions and arbitrary functions in the expression at hand
    arbitrary_functions,basis_functions = identify_basis_functions(expr,coeff_list)
    # Allocate memory for the LHS being the solutions to the equations
    # stemming from the expression at hand and the RHS being the
    # corresponding expressions for the solutions
    LHS = []
    RHS = []
    # If we have an arbitrary function then we solve the equation
    # for that, otherwise we solve for each coefficient
    if len(arbitrary_functions)!=0: # Arbitrary functions
        # Solve the equation for the arbitrary function
        sol_temp = solve(expr,arbitrary_functions[0])
        # Save the solution and the equation
        LHS.append(arbitrary_functions[0])
        RHS.append(sol_temp[0])        
    else: # No arbitrary functions
        # Create a temporary sum in order to define the
        # equation stemming from the constant if such
        # a constant exist among the basis functions
        temp_sum = 0
        # Loop through the basis functions and save
        # the solution and its corresponding expression
        for func_temp in basis_functions:
            # We ignore the constant basis function
            if func_temp != 1:
                # We define the equation which is the constant
                # in front of the basis function
                eq_temp = expr.coeff(func_temp)
                # Find the coefficient which we solve for
                LHS_temp = 0
                for koeff in coeff_list:
                    if eq_temp.coeff(koeff)!=0:
                        LHS_temp = eq_temp.coeff(koeff)*koeff
                        break
                # Solve the equation for this coefficient
                RHS_temp = solve(eq_temp,LHS_temp)[0]
                # Append both the LHS and the RHS
                LHS.append(LHS_temp)
                RHS.append(RHS_temp)                
                # Increase the temporary sum
                temp_sum += eq_temp*func_temp
        # Define the zeroth equation
        eq_zero = expr - temp_sum
        # Find the coefficient which we solve for
        LHS_temp = 0
        for koeff in coeff_list:
            if eq_zero.coeff(koeff)!=0:
                LHS_temp = eq_zero.coeff(koeff)*koeff
                break
        # Solve the equation for this coefficient
        RHS_temp = solve(eq_zero,LHS_temp)[0]
        # Append both the LHS and the RHS
        LHS.append(LHS_temp)
        RHS.append(RHS_temp)                        
    # Lastly, return the LHS and the RHS
    return LHS, RHS




# Create arbitrary constants C's
C_list = []
for i in range(1,9):
    # Create a symbol for the current constant
    exec("C%d=Symbol('C%d')"%(int(i), int(i)))
    exec("C_list.append(C%d)"%(int(i)))

print("\n\t\tC_list:\t%s"%(str(C_list)))
# Create arbitrary constants k's
for i in range(1,3):
    # Create a symbol for the current constant
    exec("k%d=Symbol('k%d')"%(int(i),int(i)))
# Create symbol which is the time
t = Symbol('t')
# Create an arbitrary function
f = symbols('f',cls=Function)
# Create the first expression
e_1 = C1 + C2*t + C3*(t**5) + C4 * exp(k1*t) + C5*sin(k2*t) + C6 * f(t) - 2*C7 * (t**5) - 3*C8
# Create the second expression
e_2 = C1 + C2*t + C3*(t**5) + C4 * exp(k1*t) + C5*sin(k2*t)- 2*C7 * (t**5) - 3*C8
# Calculate the basis functions and the arbitrary expressions
a_f_1, b_f_1 = identify_basis_functions(e_1,C_list)
a_f_2, b_f_2 = identify_basis_functions(e_2,C_list)
LHS_list_1, RHS_list_1 = solve_algebraic_equation(e_1,C_list)
LHS_list_2, RHS_list_2 = solve_algebraic_equation(e_2,C_list)
# Print
print("\n\t\tExpression 1:\n\t\t\te1\t=\t%s"%(str(e_1)))
print("\t\tBasis functions 1:\t%s"%(str(b_f_1)))
print("\t\tArbitrary functions 1:\t%s"%(str(a_f_1)))
print("\t\tSolutions:")
for index in range(len(LHS_list_1)):
    print("\t\t\t%s\t=\t%s"%(LHS_list_1[index],RHS_list_1[index]))
print("\n\t\tExpression 2:\n\t\t\te2\t=\t%s"%(str(e_2)))
print("\t\tBasis functions 2:\t%s"%(str(b_f_2)))
print("\t\tArbitrary functions 2:\t%s"%(str(a_f_2)))
print("\t\tSolutions:")
for index in range(len(LHS_list_2)):
    print("\t\t\t%s\t=\t%s"%(LHS_list_2[index],RHS_list_2[index]))


