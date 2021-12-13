#=================================================================================
#=================================================================================
# Script:"symmetry_toolbox_first_order_ODEs"
# Date: 2021-09-28
# Implemented by: Johannes Borgqvist
# Description:
# The script contains all the function that are necessary in order to
# calculate a symmetry generator for a given system of ODEs. There is a function
# called "create_tangent_ansatze" which sets up the so called projective tangent
# ansätze that are polynomials in the states where each coefficient is a function
# of the independent variable. The function "lin_sym_cond" formulates the
# linearised symmetry condition, and given these the function "determining_equations" formulates the determining equations. Lastly, the major function which does
# most of the work is called "solve_determining_equations"which solves the
# determining equations for the given system of ODEs. 
#=================================================================================
#=================================================================================
# Import Libraries
#=================================================================================
#=================================================================================
# Import SymPy
from sympy import *
from sympy.core.function import *
# Other functionalities in sympy:
# Finding monomials
from sympy.polys.monomials import itermonomials, monomial_count
# Ordering monomials 
from sympy.polys.orderings import monomial_key
# To solve ODE systems
from sympy import symbols, Eq, Function
# Import all matrix stuff
from sympy.matrices import Matrix, eye, zeros, ones, diag, GramSchmidt
# Import math for combinatorials
from math import *
import math
# To do the fancy iterations
import itertools
# To manipulate strings
import string
# To time each part of the program
import time
# To be able to break the code after a certain amount of time has passed
import signal
#=================================================================================
#=================================================================================
# The Functions
#=================================================================================
#=================================================================================
#The following six functions are defined subsequently:
#1. *create_tangent_ansatze* (defines variables and creates polynomial tangent ansätze),
#2. *Lie_generator* (Defines the partial differential operator $X$ being the infinitesimal generator of the Lie group),
#3. *total_derivative* (Defines the total derivative $D_t$),
#4. *lin_sym_cond* (Generates the linearised symmetry conditions),
#5. *determining_equations* (Generates the determining equations),
#6. *identify_basis_function* (The function changes the name of the arbitrary integration constants)
#7.  *post_processing_algebraic_solutions* (Help function 1 for function 9)
#8. *extract_equation_linear_independence* (Help function 2 for function 9)
#9. *solve_algebraic_equations* (Solve an algebraic equation one at a time)
# 10. *"integration_by_parts"* (Integration by parts on the non-homogenous source terms)
# 11. *"solve_determining_equations"* (the main function solving the determining equations)
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# FUNCTION 1: "create_tangent_ansatze"
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# The function takes three inputs:
#1. The number of variables ("num_of_variables"),
#2. The number of states ("num_of_states"), 
#2. The degree of the polynomial ("degree_polynomial").
# It returns two outputs:
# 1.The variables (independent and dependent) in the symmetry framework,
# 2. A list of all tangents ("eta_list").
# The function uses various functionalities in sympy and the polynomials are generated using the built in functions in the polys.monomial part of the sympy library. To automate the generation of the sympy functions and symbol, the "exec" command is used at multiple places in the function.
def create_tangent_ansatze(num_of_variables,num_of_states,degree_polynomial):
    # Allocate our to lists which will be the output:
    x = [] # The independent and dependent variables
    states = [] # The independent and dependent variables
    eta_list = [] # The function with the tangents        
    c = [] # The list of the coefficients
#----------------------------------------------------------------------------------
    # STEP 1: GENERATE VARIABLES    
#----------------------------------------------------------------------------------
    # Add variables to the list
    for index in range(num_of_variables+num_of_states):
        # Allocate variable
        exec("x_%d = Symbol(\'x_%d\') "%(index,index))  
        # Add variable to the list    
        exec("x.append(x_%d)"%(index)) 
        # Add state to list
        if not index < (num_of_variables):
            # Add variable to the list    
            exec("states.append(x_%d)"%(index)) 
#----------------------------------------------------------------------------------
    # STEP 2: GENERATE FOUR POLYNOMIAL VECTORS, ONE FOR EACH VARIABLE
#----------------------------------------------------------------------------------
    # Generate all monomials for our polynomials
    M = list(itermonomials(states, degree_polynomial))
    # Sort the list
    #M = sorted(M, key=monomial_key('grlex', states))
    M = sorted(M, key=monomial_key('lex', states))
    #M = sorted(M, key=monomial_key('grevlex', states))
    #M = M.reverse()
    # Calculate the number of terms in each of the polynomials
    num_of_terms = monomial_count(num_of_states,degree_polynomial)
#----------------------------------------------------------------------------------
    # STEP 3: DEFINE THE TANGENTIAL ANSÄTZE
#----------------------------------------------------------------------------------
    # Loop over the variables and create our tangent eta for each variable
    # and add it to the eta_list. Lastly, reset it to zero and start over 
    # again. 
    for iteration in range(num_of_variables+num_of_states):
        # Create our tangent vector
        exec("eta_%d = symbols(\'eta_%d\', cls=Function) "%(iteration,iteration))  
        exec("eta_%d = 0"%(iteration))
        # Loop over the desired terms
        for index in range(num_of_terms):            
            # Calculate a temporary index
            temp_index = index
            if index != 0:
                temp_index = (len(range(num_of_terms)) - 1)- (index-1)             
            # Allocate coefficient
            exec("c_%d_%d = symbols(\'c_%d_%d\', cls=Function) "%(iteration,temp_index,iteration,temp_index))  
            # Save the coefficient in the list
            exec("c.append(c_%d_%d)"%(iteration,temp_index))              
            # Print variable    
            exec("eta_%d+=c_%d_%d(x_0)*M[%d]" % (iteration,iteration,temp_index,index))         
        # Add our tangent to the list of tangents "eta_list"
        exec("eta_list.append(eta_%d)"% (iteration))
    # Return the output    
    return x, c, eta_list, M
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# FUNCTION 2: "Lie_generator"
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# The function takes three inputs:
# 1. The variables x,
# 2. The tangents in the eta_list,
# 3. The function f which is to be differentiated.
# The function returns the output which is the "derivative" corresponding to letting the Lie generator X act on the function f
def Lie_generator(x,eta_list,f):
    # Define the derivative which we will return
    der = symbols('der', cls=Function)
    der = 0
    # Loop over all variables and add terms succesively
    for i in range(len(x)):
        der += eta_list[i]*diff(f,x[i])
    return der # Return the derivative
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# FUNCTION 3: "total_derivative"
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# The function takes three inputs:
# 1. The variables x,
# 2. The omega list,  
# 3. The function f which is to be differentiated.
# The function returns the output which is the "derivative" corresponding to taking the total derivative D_t on the function f
def total_derivative(x,omega_list,f):
    # Calculate the number of states
    num_of_states = len(omega_list)
    # Calcuate the number of variables using this value
    num_of_variables = len(x) - num_of_states
    # Define the derivative which we will return
    der = symbols('der', cls=Function)
    der = 0 # Initialise it
    # Loop over all variables and add terms succesively
    for i in range(len(x)):
        if i < num_of_variables: # INDEPENDENT VARIABLES
            der += diff(f,x[i])
        else:# DEPENDENT VARIABLES (I.E. THE STATES)
            der += omega_list[i-num_of_variables]*diff(f,x[i])
    return der # Return the derivative
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# FUNCTION 4: "lin_sym_cond"
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# The function takes three inputs:
# 1. The variables x,
# 2. The tangents in the eta_list,
# 3. The omega list,  
# The function returns the list "lin_sym_list" which is a list which is as long as the number of states and it contains one equation corresponding to one linearised symmetry condition (one for each state).
def lin_sym_cond(x,eta_list,omega_list):
    #------------------------------------------------------------------------------
    # STEP 1: CALCULATE THE LINEARISED SYMMETRY CONDITIONS
    #------------------------------------------------------------------------------
    # Define list of the symmetry conditions
    lin_sym_list = []
    # Calculate the number of states
    num_of_states = len(omega_list)   
    # Calcuate the number of variables using this value
    num_of_variables = len(x) - num_of_states    
    # Define the derivative which we will return
    sym_temp = symbols('sym_temp', cls=Function)
    sym_temp = 0
    # Loop over the number of states and calculate 
    # the linearised symmetry conditions
    for i in range(num_of_states):
        # Term 1: X(omega)
        sym_temp += Lie_generator(x,eta_list,omega_list[i])
        # Term 2: D_t eta
        sym_temp -= total_derivative(x,omega_list,eta_list[i+num_of_variables])
        # Term 3: omega * D_t eta_0
        sym_temp += omega_list[i]*total_derivative(x,omega_list,eta_list[num_of_variables-1])
        # Expand the term so that we do not have any factors
        sym_temp = expand(sym_temp)
        # Add this to the list
        lin_sym_list.append(sym_temp)
        # Reset the symmetry condition
        sym_temp = 0
    #------------------------------------------------------------------------------
    # STEP 2: SIMPLIFY THE LINEARISED SYMMETRY CONDITIONS
    #------------------------------------------------------------------------------
    # Define a help counter
    help_counter = 0
    # Loop over the symmetry conditions and simplify them in the following ways:
    # 1. Write them as a fraction using the function "simplify",
    # 2. Extract only the numerator using the function "fraction",
    # 3. Cancel all terms, expand them by performing multiplications and organise them in terms
    # of the various monomials with respect to the states in the vector x. 
    for lin_sym in lin_sym_list:
        # Write the linearised symmetry condition as a fraction 
        # with a common denominator
        tidy_eq = simplify(lin_sym)  
        # Extract the nominator
        tidy_eq, denom = fraction(tidy_eq)
        # Cancel terms, expand by performing multiplications and organise in
        # terms of monomials
        lin_sym_list[help_counter] = powsimp(expand(tidy_eq))
        # Increment the help_counter
        help_counter += 1
    return lin_sym_list # Return the equations
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# FUNCTION 5: "determining_equations"
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# The function takes three inputs: 
# 1. "x" the list of variables,
# 2. "lin_sym_list" the list of the linearised symmetry conditions in polynomial form,
#3. A list called "degree" which contains the highest and lowest degrees of all the states.
# The function returns a (huge) list of the determining equations as output. 
def determining_equations(x,lin_sym_list,degree_monomial):
    # In essence, the linear symmetry condition in the case of a polynomial ansatz
    # corresponds to find the roots of a multivariate polynomial. What we do here
    # is that we have a polynomial in the states of the models which is zero and
    # we force each of the coefficients in this polynomial to be zero in order for
    # this to hold for all values of the states. It is setting the coefficients
    # of the monomials in this multivariate polynomial that corresponds to the
    # determining equations. Implementation-wise, we 
    # define a list of the powers of a given monomial.
    # The monomial then combines the powers in all these
    # lists to give us all the monomials in the multivariate
    # polynomial.
    degree_list = list(range(degree_monomial + 1))
    # Calculate the number of states
    num_of_states = len(x) - 1
    # Generate a list which has three lists of powers
    # of the monomials. First allocate memory:
    degree = []    
    # Define the degrees for our states in the model
    # by appending the degree list to it (one for each state). 
    for i in range(num_of_states):
        degree.append(degree_list)
    # Using the degree list, the monomials will be sorted
    # by combining the powers of each state in all possible
    # ways...
    # Allocate memory for the determining equations
    det_eq = []        
    # This is the KEY STEP, and it is solved thanks to "itertools".
    # We generate all powers in the respective monomials in our 
    # multivariate polynomial using the product iterator which 
    # calcluates all the combinations in the "degree" list. 
    powers = list(itertools.product(*degree))        
    # Calculate the number of states
    num_of_states = len(lin_sym_list)
    # Calculate the number of variables
    num_of_variables = len(x) - num_of_states
    # Allocate some memory for the monomials
    # That is, we wish to know which monomial
    # that each determining equation comes from
    monomials = []
    # State which symmetry condition that the determining
    # equation comes from
    lin_sym_eq_number = []
    # To keep track of this one we define a help index 
    help_index = -1
    # Loop over the linearised symmetry conditions
    for lin_eq in lin_sym_list:        
        # Increase the help_index so that we can
        # save the determining equation
        help_index += 1
        # Loop over all powers
        for power in powers:        
            # We initiate the monomial to our symmetry condition
            # of interest
            mon = lin_eq
            # Loop over all states
            for i in range(num_of_states):            
                # Collect the monomials with respect to
                # the variable of interest
                mon = collect(mon, x[num_of_variables+i])
                # Extract the power of interest
                mon = mon.coeff(x[num_of_variables+i],power[i])   
            # Save the determining equation if it is non-zero
            # Save the monomial i.e. the tuple
            # Save the number of the symmetry condition
            if mon != 0:
                # DETERMINING EQUATION 
                det_eq.append(mon)
                # THE MONOMIAL
                monomials.append(power)
                # THE PRECISE EQUATION AMONG
                # THE LINEARISED SYMMETRY CONDITIONS
                # THAT THE EQUATION STEMS FROM
                lin_sym_eq_number.append(help_index)
    # Return the determining equations
    return det_eq, monomials, lin_sym_eq_number
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# FUNCTION 6: "identify_basis_functions"
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# The function takes three inputs which are the following:
# 1. An algebraic equation being a linear combination of functions stored in the variable "alg_eq",
# 2. A list of all arbitrary integration coefficients stored in "const_list",
# 3. The list of variables called x (it is really only the independent variable x[0] that is of interest here).
# The script returns two outputs:
# 1. A list of arbitrary functions stored in "arbitrary_functions",
# 2. A list of all basis functions stored in "basis_functions".
def identify_basis_functions(alg_eq,const_list,x):
    # Define the arbitrary functions
    arbitrary_functions = list(alg_eq.atoms(AppliedUndef))
    # Allocate memory an empty list for the basis functions
    basis_functions = []
    # ---------------------------------------------------
    # Step 1: Extract all arguments in the algebraic equation
    # ---------------------------------------------------    
    # Define a list with all arguments (or terms) in the
    # algebraic equation
    arguments = list(alg_eq.args)
    # Loop through all arguments
    for argument in arguments:
        # Loop through all coefficients 
        for constant in const_list:
            # Save all basis functions provided that
            # they exist as a linear combination of the
            # arbitrary integration constants in coeff_list
            if argument.coeff(constant) != 0: 
                basis_functions.append(argument.coeff(constant))
    # ---------------------------------------------------
    # Step 2: Loop through the basis functions and normalise
    # them in order to remove constants or factors. We
    # only want to keep the important part of the basis
    # function so that if we obtain "-3t**2" we only keep
    # "t**2" as a basis function. There are three keys to
    # this: (1) fraction to get the nominator and denominator,
    # (2) factor_list to get all factors in the product,
    # (3) Derivative to just extract the functions of t.
    # ---------------------------------------------------    
    # Loop over the basis functions
    for base_index in range(len(basis_functions)):
        # Extract the numerator and the denominator from the basis function
        nom, denom = fraction(basis_functions[base_index])
        # The basis function will be rational in the sense that it contains
        # a contribution from the numerator num and denominator denom above.
        # These will be constructed as a product of all the factors in num
        # and denom respectively. 
        # Initiate these functions with 1
        f_nom = 1 # Nominator
        f_denom = 1 # Denominator
        # Start with the nominator: If we have
        # any factors in this list which are functions
        # of the independent variable x[0], we multiply
        # the current function with this factor
        if len(factor_list(nom)[1])!=0:
            # Loop over the factors and multiply
            for factor_tuple in factor_list(nom)[1]:
                # If it is a function of the independent variable x[0]
                # then we extract this function and multiply the contribution
                # to the basis function from the nominator
                if Derivative(factor_tuple[0],x[0]).doit()!=0: # We have a function!
                    f_nom = f_nom*((factor_tuple[0])**(factor_tuple[1]))
        # Do the same steps with the denominator: If we have
        # any factors in this list which are functions
        # of the independent variable x[0], we multiply
        # the current function with this factor
        if len(factor_list(denom)[1])!=0:
            # Loop over the factors and multiply
            for factor_tuple in factor_list(denom)[1]:
                # If it is a function of the independent variable x[0]
                # then we extract this function and multiply the contribution
                # to the basis function from the nominator
                if Derivative(factor_tuple[0],x[0]).doit()!=0: # We have a function!
                    f_denom = f_denom*((factor_tuple[0])**(factor_tuple[1]))
        # Lastly, update the basis function at hand
        basis_functions[base_index] = ((f_nom)/(f_denom))
    # ---------------------------------------------------
    # Step 3: Only return the unique basis functions
    # ---------------------------------------------------    
    # Lastly, we only want unique basis functions so we take the
    # set of these functions    
    basis_functions = list(set(basis_functions).difference(set(arbitrary_functions)))
    # Return the output
    return arbitrary_functions, basis_functions
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# FUNCTION 7: "post_processing_algebraic_solutions"
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# The function takes four inputs which are the following:
# 1. A list of the left hand sides of the solved equations called LHS_list,
# 2. A list of the right hand sides of the solved equations called RHS_list,
# 3. A list of arbitrary integration constants called "const_list"
# The script returns the following outputs:
# 1. A list of the updated left hand sides called LHS_list_updated,
# 2. A list of the updated right hand sides called RHS_list_updated.
#----------------------------------------------------------------------------------
def post_processing_algebraic_solutions(LHS_list,RHS_list,const_list):
    # Allocate memory for the output
    LHS_list_updated = []
    RHS_list_updated = []
    # Allocate memory for the elements being removed which are
    # equations of the type "0=0"
    indices_remove_solutions = [i for i in range(len(LHS_list)) if (LHS_list[i]==0 and RHS_list[i]==0)]
    # Sort these indices in reverse
    indices_remove_solutions.sort(reverse=True)
    # Remove all useless equations with no information in them
    for index in indices_remove_solutions:
        del LHS_list[index]
        del RHS_list[index]
    # Allocate memory for the constant that exist in our equations
    const_mat = []
    # Loop through all constants and save the ones that exist in our
    # equations
    for constant in const_list:
        # Loop through the solutions as well
        for i in range(len(LHS_list)):
            # Save the ones that exist in our equations
            if LHS_list[i].coeff(constant)!=0 or RHS_list[i].coeff(constant)!=0:
                const_mat.append(constant)
    # Lastly, we are merely intested in the unique constants
    const_mat = list(set(const_mat))
    # Create a solution list which we will later convert to a matrix
    sol_list = [LHS_list[i]-RHS_list[i] for i in range(len(LHS_list))]
    # Now we will make a matrix out of the solutions in sol_list
    A = []
    # Loop through the solutions
    for rI in range(len(sol_list)):
        # Loop through the coefficients
        for cI in range(len(const_mat)):
            # Save the coefficients of each integration constant
            A.append(sol_list[rI].coeff(const_mat[cI]))
    # Now, we convert this into a matrix
    A = Matrix(len(sol_list),len(const_mat),A)
    # Convert the list of constants into a matrix too
    const_mat = Matrix(len(const_mat),1,const_mat)
    # Next, we row-reduce the matrix A to remove linearly dependent solutions
    A, pivot_columns = A.rref()
    # Remove all zero rows of A:
    # Caclulate the non-zero rows
    rows_nonzero = [i for i in range(len(sol_list)) if any(A[i, j] != 0 for j in range(len(const_mat)))]
    # Update A
    A = A[rows_nonzero, :]
    # Now, we define an updated RHS list
    RHS_list_updated = A*const_mat
    # The left hand side consists of the pivot constants
    LHS_list_updated = [const_mat[pC] for pC in list(pivot_columns)]
    # Update the RHS
    RHS_list_updated = list(-(RHS_list_updated-Matrix(len(LHS_list_updated),1,LHS_list_updated)))
    # Return the output
    return LHS_list_updated, RHS_list_updated
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# FUNCTION 8: "extract_equation_linear_independence"
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# The function takes three inputs which are the following:
# 1. An expression expr which defines an equation of the type "expr=0" with
# linearly independent basis functions,
# 2. A basis function basis_function which is not zero and is a function of
# the independent variable x[0],
# 3. A list of the variables x where the independent variable x[0] is of
# interest.
# The script returns the following output:
# 1. The extracted equation called extracted_equation.
#----------------------------------------------------------------------------------
def extract_equation_linear_independence(expr,basis_function,x):
    # Allicate memory for the equation we return at the end
    extracted_equation = 0
    # Loop through the terms in the expression at hand
    # and add the constants in front of each basis function
    for term in expr.args:
        # Add the term if it contains our basis function at hand
        if Derivative(simplify(term/basis_function),x[0]).doit()==0:
            extracted_equation += simplify(term/basis_function)
    # Return the extracted equation
    return extracted_equation
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# FUNCTION 9: "solve_algebraic_equation"
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# The function takes three inputs which are the following:
# 1. An arbitrary expression being a linear combination of functions stored in the variable "expr",
# 2. A list of all arbitrary integration coefficients stored in "coeff_list",
# 3. A list of the variables x (where the independent variable x[0] is used to identify the basis functions). 
# The script returns two outputs:
# 1. A list of the LHSs being the solutions stored in "LHS",
# 2. A list of the RHSs being the expressions stored in "RHS".
# The function uses the previous function (i.e. Function 6) called "identify_basis_functions" in
# order to find the basis functions and the arbitrary functions in the expression at hand.
def solve_algebraic_equation(expr,coeff_list,x):
    eq_str = ""
    # Extractthe nominator
    expr, denom = fraction(expr)
    # In case it is we only keep the numerator
    if denom != 1:
        expr = nom
    # Then we expand our lovely expression as far as we can
    expr = expand(expr)
    # Find all basis functions and arbitrary functions in the expression at hand
    arbitrary_functions,basis_functions = identify_basis_functions(expr,coeff_list,x)
    # Allocate memory for the LHS being the solutions to the equations
    # stemming from the expression at hand and the RHS being the
    # corresponding expressions for the solutions
    LHS = []
    RHS = []
    # Equations that we do not need in most cases
    LHS_before = []
    RHS_before = []
    # Allocate an empty equation string as well
    eq_Str = ""
    # If we have an arbitrary function then we solve the equation
    # for that, otherwise we solve for each coefficient
    if len(arbitrary_functions)!=0: # Arbitrary functions
        # We have two cases here, one difficult and one ease choice. The latter
        # one is when we just have an arbitrary function which we directly
        # solve the expression for. The harder case is when we have an unknown
        # integral which we must solve for. In this case, it is necessary
        # to do some work in order just to extract this arbitrary function in
        # sympy. We use an if-statement where the first statement is the difficult
        # case and the second statement is the easy case
        if expr.coeff(arbitrary_functions[0]) == 0: # WE HAVE AN INTEGRAL TERM
            # Define a list of all arguments in which the unknown integral
            # term resides
            arguments = expr.args
            # Loop over all term and find the integral term which is characterised
            # by the fact that it has an undefined function in it
            for argument in arguments:
                # We have an unknown function
                if len(argument.atoms(AppliedUndef)) != 0:
                    unknown_term = argument
                    break
            # Define a list of all factors in this integral term
            factors = factor_list(unknown_term)
            # Extract the unknown integral term which is hidden in
            # the back of the factor list
            unknown_function = factors[1][len(factors[1])-1][0]
            # Calculate the factor in front of this term
            unknown_factor = simplify(unknown_term/unknown_function)
            # Now we can solve the equation at hand
            sol_temp = solve(expr,unknown_function)
            # Save the left hand side
            LHS.append(unknown_function)
            # Save the right hand side
            RHS.append(sol_temp[0])
        else: # WE DO NOT HAVE AN INTEGRAL TERM
            # Solve the equation for the arbitrary function
            sol_temp = solve(expr,arbitrary_functions[0])
            # Save the solution and the equation
            LHS.append(arbitrary_functions[0])
            #RHS.append(sol_temp[0])
            RHS.append(expand(sol_temp[0]))        
    else: # No arbitrary functions
        # Create a temporary sum in order to define the
        # equation stemming from the constant if such
        # a constant exist among the basis functions
        temp_sum = 0
        # Loop through the basis functions and save
        # the solution and its corresponding expression
        for basis_function in basis_functions:
            # We ignore the constant basis function
            if basis_function != 1:
                #==============================================================================
                # We define the equation which is the constant
                # in front of the basis function
                #eq_temp = expr.coeff(func_temp)
                eq_temp = extract_equation_linear_independence(expr,basis_function,x)
                # Printing statements
                eq_str += "Equation for the basis function $" + latex(basis_function) + "$:\n"
                eq_str += "$$" + latex(eq_temp) + "$$\n"
                # Find the coefficient which we solve for
                LHS_temp = 0
                for koeff in coeff_list:
                    if eq_temp.coeff(koeff)!=0:
                        LHS_temp = koeff
                        break
                eq_str += "This equation was solved for: $" + latex(LHS_temp) + "$ which gave the solution:\n"
                # Solve the equation for this coefficient
                #RHS_temp = solve(eq_temp,LHS_temp)[0]
                RHS_temp = expand(solve(eq_temp,LHS_temp)[0])
                eq_str += "$$" + latex(RHS_temp) + "$$\n" 
                # Append both the LHS and the RHS
                LHS.append(LHS_temp)
                RHS.append(RHS_temp)                
                # Increase the temporary sum
                temp_sum += eq_temp*basis_function
        # Define the zeroth equation
        eq_zero = expand(expr - temp_sum)
        # Find the coefficient which we solve for
        LHS_temp = 0
        for koeff in coeff_list:
            if eq_zero.coeff(koeff)!=0:
                LHS_temp = koeff
                break
        # Solve the equation for this coefficient
        RHS_temp = expand(solve(eq_zero,LHS_temp)[0])        
        # Append both the LHS and the RHS
        LHS.append(LHS_temp)
        RHS.append(RHS_temp)
        # Return the solutions before for clarity's sake
        LHS_before = LHS
        RHS_before = RHS
        # Lastly, we do a quick post-processing of the solutions we have
        # which makes sure that we have linearly independent solutions only
        LHS,RHS = post_processing_algebraic_solutions(LHS,RHS,coeff_list)                        
    # Lastly, return the LHS and the RHS
    return LHS, RHS, basis_functions, LHS_before, RHS_before, eq_str
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# FUNCTION 10: "integration_by_parts"
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# This function conduct the integration by parts necessary to simplify all
# non-homogeneous terms in case there are non-pivot columns in the ODE system.
# The inputs are:
# 1. The "function_list" being a list of all involved arbitrary unknown functions,
# 2. The "constant_list" being a list of all arbitrary integration constants
# 3. The "integrand_list" being a list of all integrands that are of the type
# "df/dt*exp(k*t)",
# 4. The variable "variable" which is the current integration variable,
# 4. The dummy variable "s" which is used for the integration.
# stemming from the integration by parts where the arbitrary functions are evaluated
# at zero,
# 3. The "integrand_list" being a list of all integrands,
# 4. The "variable" which in this setting is just x[0] being the independent variable.
# The output is:
# 1. The "integrallist" being a list of all integrals.
def integration_by_parts(function_list,constant_list,integrand_list,variable,dummy):
    # At the end, we want to return the correct integrals
    integral_list = []
    # Loop through the integrands
    for integrand in integrand_list:
        # Allocate a temporary integral
        temp_int = 0
        # Loop through the functions and find the coefficients in the integrand
        for func_index in range(len(function_list)):
            # Extract the function at hand
            func = function_list[func_index]
            # If we have no coefficients, we'll just move on
            if integrand.coeff(Derivative(func(variable),variable)) == 0:
                continue
            else:# Otherwise, we extract the other function in the integrand
                # Extract the coefficient which is the "other function" in the integrand
                other_function = integrand.coeff(Derivative(func(variable),variable))
            # Now, we have two options:
            # 1. The coefficient is a constant (e.g. 1, 5 or k),
            # 2. The coefficient is a function meaning that we have to use integration
            # by parts. 
            # Alternative 1: Constant
            if len(list(other_function.atoms(Function)))==0:
                # Just solve the integral
                temp_int += (other_function*func(variable)) - (other_function*constant_list[func_index])
                # If we have a polynomial of t we need to add an extra term
                if degree(other_function,variable)!=0:
                    # The integral in the integration by parts
                    new_integrand = (func(variable)*Derivative(other_function,variable).doit()).subs(variable,dummy)
                    # Add the integral term as well
                    temp_int += -Integral(new_integrand,(dummy,0,variable))
            else:# Alternative 2: Integration by parts
                # The boundary term in the integration by parts
                temp_int += (other_function*func(variable)) - (other_function.subs(variable,0)*constant_list[func_index]) 
                # The integral in the integraiton by parts
                new_integrand = (func(variable)*Derivative(other_function,variable).doit()).subs(variable,dummy)
                temp_int += -Integral(new_integrand,(dummy,0,variable))
        # Evaluate the integrand if possible
        if temp_int != 0:
            temp_int = temp_int.doit()
        # Add the evaluated integrals to our integral list
        integral_list.append(temp_int)
    # Return the integral list
    return integral_list
#----------------------------------------------------------------------------------
# FUNCTION 11: "solve_determining_equations"
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# The function takes four inputs which are the following:
# 1. The variables (dependent and independent) in a vector "x",
# 2. The tangent vector called "eta_list",
# 3. The coefficients in the tangents called "c"
# 4. The list of the determining equations called "det_eq".
# The script return the following output:
# 1. The calculated generator of the symmetry which is stored in a string called "X".
def solve_determining_equations(x,eta_list,c,det_eq,variables,omega_list,M):
    alg_str = ""
    #------------------------------------------------------------------------------
    # STEP 0 of 6: ENABLE FOR EASY CALCULATIONS OF TANGENTS WITHOUT SUBSTITUTIONS
    #-----------------------------------------------------------------------------
    # The whole aim of this function is to determine the unknown coefficients stored
    # in the list c which is provided as an input. Later in this script, the vector
    # c_mat (technically stored as a Matrix) will contain the calculated constants
    # after solving the involved matrix equations. Initially, the tangents corresponding
    # to the final solution were retrieved by substituting the value of the elements
    # in the list c by the corresponding element in c_mat. This substitution was
    # conducted using the "subs" function in sympy. However, it turns out that when
    # the expressions in c_mat contains integrals of arbitrary functions where the
    # integrals contain a dummy variable, then it is not possible for "subs" to
    # conduct the substitution. Therefore another approach is implemented where
    # each coefficient in c has a matching monomial in the list "monomial_list"
    # and a matching tangent in the list "tangent_indicator". In this way, it is
    # just a matter of looping through one list (as these three lists have the
    # same dimensions) and then reconstruct the tangents from there. This part
    # of the code sets up these lists which will enable an easy assembly of the
    # tangents in the very end of the code after all calculations are done.
    # THE CODE BEGINS
    # Create a list of the monomials in the tangential ansätze of a certain order
    monomial_list = list(M)
    # Allocate memory for a list indicating which tangent the corresponding monomial
    # belongs to
    tangent_indicators = []
    # Fill this list with values. These are indices indicating which tangent the
    # particular monomial and the particular coefficient belong to.
    for index in range(len(c)):
        # Since we always use the same number of unknown coefficients in
        # the ansätze, it is possible to assign an integer to each tangent
        # based on the number of monomials we have in the ansätze
        tangent_indicators.append(math.floor(index/len(monomial_list)))        
    # Calculate the number of tangents
    num_tangents = len(eta_list)
    # Add the monomial_list to itself "num_tangents-1" times so that each
    # coefficient in c matches a monomial. 
    for index in range(num_tangents-1):
        monomial_list += monomial_list
    # Allocate a vector for the tangents which we will calculate
    eta_list_final = []
    # Loop over the tangents in the ansätze, and initialise this
    # list with zeroes. 
    for index in range(len(eta_list)):
        eta_list_final.append(0)
    #for index in range(len(c)):
        #eta_list_new[tangent_indicator[index]] += c[index](x[0])*monomial_list[index]
    #------------------------------------------------------------------------------
    # STEP 1 of 6: FORMULATE THE DETERMINING EQUATIONS ON MATRIX FORM A*dc/dt=B*c
    #-----------------------------------------------------------------------------
    # Allocate memory for the two matrices. 
    # We allocate them as lists because we can
    # easily generate matrices from a list where
    # we provide the list of course but before 
    # that we just state how many rows and columns
    # we have. 
    A_mat = []
    B_mat = []
    help_counter = 1
    A_row = []
    B_row = []    
    # Loop over all ODEs and save the coefficients
    for isolated_eq in det_eq:
        # Expand the equation at hand
        eq_temp = expand(isolated_eq)
        # Loop over all coefficients in the tangential ansätze
        # and define our matrices
        for koeff in c:
            # Save the derivative of the A coefficient
            A_mat.append(-eq_temp.coeff(Derivative(koeff(x[0]),x[0])))
            # Save the coefficient of the coefficient in the tangential
            # ansätze
            B_mat.append(eq_temp.coeff(koeff(x[0])))
    # Finally we generate the matrices from our lists A_mat and B_mat
    num_of_rows = len(det_eq) # Number of rows
    num_of_cols = len(c) # Number of columns
    A = Matrix(num_of_rows,num_of_cols,A_mat) # Matrix A
    B = Matrix(num_of_rows,num_of_cols,B_mat) # Matrix B
    #------------------------------------------------------------------------------
    # STEP 2 of 6: TRANSFER THE SYSTEM TO ONE MATRIX M=[-A|B] AND REDUCE THE NUMBER OF 
    # EQUATIONS BY FINDING A BASIS FOR THE COLUMN SPACE OF M^T. FINISH BY SPLITTING
    # UP THE MATRIX INTO ITS COMPONENTS PART, I.E. A AND B RESPECTIVELY
    #-----------------------------------------------------------------------------
    # Find the dimensions of the matrix A
    num_of_eq, n = A.shape        
    # FIRST ROW OF M
    # Begin by giving the matrix A the value of the first row 
    # in A
    M = -(A.row(0))
    # Then we add the zeroth row of B in the columns to the left
    for j in range(len(c)):
        M = M.col_insert(j+len(c), Matrix([B[0,j]]))
    # THE REMAINING ROWS OF M
    for i in range(num_of_eq):
        if i > 0:
            temp_row = -(A.row(i))
            # Then we add the zeroth row of B in the columns to the left
            for j in range(len(c)):
                temp_row = temp_row.col_insert(j+len(c), Matrix([B[i,j]]))
            # Insert the temporary row
            M = M.row_insert(i, temp_row)
    # Let us calculate a basis for "col(M^T)". 
    eq_sys = M.T # Take the transpose of the matrix
    reduced_sys = eq_sys.columnspace()
    M_tilde = reduced_sys[0].T
    for i in range(len(reduced_sys)):
        if i > 0:
            M_tilde = M_tilde.row_insert(i, reduced_sys[i].T)          
    # Row-reduce the expanded matrix
    M_tilde = M_tilde.rref()[0]
    # Split the matrix into its component parts
    A = M_tilde[:,0:(len(c))]
    B = -M_tilde[:,len(c):(2*len(c))]
    #------------------------------------------------------------------------------
    # STEP 3 of 6: FIND THE PURELY ALGEBRAIC EQUATIONS WHICH CORRESPONDS TO THE ZERO ROWS IN THE MATRIX A. THE SAME ROWS IN THE MATRIX B WILL BE FORMULATED INTO A NEW MATRIX B_algebraic WHICH CONTAINS ONLY THE ALGEBRAIC EQUATIONS
    #-----------------------------------------------------------------------------
    # Remove the zero rows from A:
    # Calculate the dimensions 
    m, n = A.shape
    # Caclulate the non-zero rows
    rows_nonzero = [i for i in range(m) if any(A[i, j] != 0 for j in range(n))]
    # Caclulate the zero rows
    rows_zero = [i for i in range(m) if i not in rows_nonzero]
    # Update A
    A = A[rows_nonzero, :]
    # Extract zero bit from B which constitutes our algebraic equations
    B_algebraic = B[rows_zero,:]
    # Update B
    B = B[rows_nonzero,:]
    #------------------------------------------------------------------------------
    # STEP 4 of 6: REMOVE (POTENTIAL) EXTRA COLUMNS 
    #-----------------------------------------------------------------------------
    # Calculate rows and columns of matrices
    r_algebraic, cols = B_algebraic.shape # Dimensions of B_algebraic
    rows, cols = A.shape # Dimensions of A
    # Create a set of all the columns
    col_set = set(range(cols))
    # Find the pivot columns of A
    pivot_columns = A.rref()[1]
    # Calculate the non-pivot columns
    non_pivot_columns = list(col_set.difference(set(pivot_columns)))
    # We assume that we have a homogeneous system
    non_homogeneous = False
    # Check if we have non-pivot columns, that is we have
    # an underdetermined system with extra columns (i.e. more columns than equations).
    # In this case, we move these column to the right hand side and we then treat
    # these as source terms if you will, so that we get a quadratic system with some
    # extra inhomogeneities. 
    if len(non_pivot_columns)!=0: # Yes, we do have non-pivot coefficients
        # Indicate that we have a non_homogeneous system
        non_homogeneous = True
        # Define a list of all non-pivot elements.
        # We will refer to this as a non-homogeneous
        # source term
        source_ODE = []
        source_ODE_der = []
        source_algebraic = []        
        # DEFINE SOURCE TERM FOR ODEs
        for r in range(rows):
            # Define two temporary source terms
            temp_source = 0
            temp_source_der = 0
            # Loop over the non-pivot columns and perform
            # the matrix multiplication in order to define
            # the source term
            for pC in non_pivot_columns:
                temp_source_der += -A[r,pC]*Derivative(c[pC](x[0]),x[0])
                temp_source += B[r,pC]*c[pC](x[0])
            # Append the temporary source
            source_ODE.append(temp_source)
            source_ODE_der.append(temp_source_der)            
        # DEFINE SOURCE TERM FOR ALGEBRAIC EQUATIONS
        for r in range(r_algebraic):
            # Define a temporary source term
            temp_source_alg = 0
            # Loop over the non-pivot columns and perform
            # the matrix multiplication
            for pC in non_pivot_columns:
                temp_source_alg += B_algebraic[r,pC]*c[pC](x[0])
            # Append the non-homogeneous source 
            source_algebraic.append(temp_source_alg)
        # Collect the source terms
        source_ODE = Matrix(len(source_ODE),1,source_ODE)
        source_ODE_der = Matrix(len(source_ODE_der),1,source_ODE_der)        
        source_alg = Matrix(len(source_algebraic),1,source_algebraic)
    # Save the unknown arbitrary functions being the non-pivot elements
    non_pivot_functions = [c[non_pivot_column] for non_pivot_column in non_pivot_columns]
    # Save the monomials of the unknown non-pivot functions as well
    non_pivot_monomials = [monomial_list[non_pivot_column] for non_pivot_column in non_pivot_columns]
    # Save the tangent indicator of the unknown non-pivot functions as well
    non_pivot_tangent_indicators = [tangent_indicators[non_pivot_column] for non_pivot_column in non_pivot_columns]    
    # Remove all non-pivot tangential coefficient
    non_pivot_columns.sort(reverse=True)
    # Loop over all non_pivot_columns and remove them from our three
    # matrices A, B and B_algebraic and from the coefficient vector c
    for index in range(len(non_pivot_columns)):
        # Remove column matrices
        A.col_del(non_pivot_columns[index])
        B.col_del(non_pivot_columns[index])
        B_algebraic.col_del(non_pivot_columns[index])
        # Remove the non-pivot parameters from the coefficient vector c,
        # the monomial_list and the tangent_indicators
        del c[non_pivot_columns[index]]
        del monomial_list[non_pivot_columns[index]]
        del tangent_indicators[non_pivot_columns[index]]
    # Make a copy of the original coefficients
    c_original = c
    #------------------------------------------------------------------------------
    # STEP 5 of 6: SOLVE THE ODE SYSTEM PROVIDED THAT IT IS QUADRATIC AND
    # SUBSTITUTE THE ALGEBRAIC EQUATIONS
    #-----------------------------------------------------------------------------
    # Calculate the dimensions of B
    num_rows, num_cols = B.shape
    r_algebraic, num_cols = B_algebraic.shape    
    # Check if B is quadratic which it must be in order for the script to
    # find the solution
    if num_rows == num_cols:
        # Check that the matrix A is an identity matrix
        if A != eye(num_rows):#If not we stop the script
            X = "\\Huge\\textsf{$A$ is quadratic but not an identity matrix!}\normalsize\\[2cm]" 
        else: # A is an identity matrix
            print("\t\t\t\tBeginning to solve ODE system...")
            X = ""
            #----------------------------------------------
            # PART 1: Allocate arbitrary coefficients
            #----------------------------------------------
            # Allocate memory for a constant list
            constant_list = []
            # Loop over the number of desired constants
            for i in range(len(c)):
                # Create a symbol for the current constant
                exec("C%d=Symbol('C%d')"%(int(i+1),int(i+1)))
                # Append the constant
                exec("constant_list.append(C%d)"%(int(i+1)))
            # Define a coefficient matrix
            c_mat = Matrix(len(constant_list),1,constant_list)
            # Allocate memory for the constants stemming from the inhomogeneities
            const_inhomo = []
            # Loop over the number of desired constants
            for i in range(len(c),len(c)+len(non_pivot_functions)):
                # Create a symbol for the current constant
                exec("C%d=Symbol('C%d')"%(int(i+1),int(i+1)))
                # Append the constant
                exec("const_inhomo.append(C%d)"%(int(i+1)))
            #----------------------------------------------
            # PART 2: Solve ODE system
            #----------------------------------------------
            t_start = time.perf_counter()
            print("\t\t\t\t\tCalculating Jordan form...")
            # Simplify each element if possible
            B_list = [simplify(cancel(B[i,j])) for i in range(num_rows) for j in range(num_cols)]
            B = Matrix(num_rows,num_cols,B_list)
            # Calculate Jordan form
            P,J = B.jordan_form()
            t_end = time.perf_counter()
            print("\t\t\t\t\t\tIt took %d seconds..."%(round(t_end-t_start)))
            # Calculate the inverse matrix of P once and once only
            t_start = time.perf_counter()
            print("\t\t\t\t\tCalculating inv(P)...")
            P_inv = P.inv()
            t_end = time.perf_counter()
            print("\t\t\t\t\t\tIt took %d seconds..."%(round(t_end-t_start)))            
            # Re-define our matrix J by scaling it by the
            # independent variable x[0]
            t_start = time.perf_counter()
            print("\t\t\t\t\tCalculating exp(Jt)...")            
            J = x[0]*J
            # We also calculate the inverse matrix
            #J_inv = -x[0]*J            
            # Re-define J by taking the matrix exponential
            J = J.exp()
            t_end = time.perf_counter()
            print("\t\t\t\t\t\tIt took %d seconds..."%(round(t_end-t_start)))            
            #J = expMt(J,x)
            # Re-define J_inv by taking the matrix exponential
            #J_inv = J_inv.exp()
            t_start = time.perf_counter()
            print("\t\t\t\t\tCalculating inv(exp(Jt))...")
            J_inv = J.inv()
            t_end = time.perf_counter()
            print("\t\t\t\t\t\tIt took %d seconds..."%(round(t_end-t_start)))            
            # Calculate matrix factor 1
            mat_factor_1 = (P*J*P_inv)
            # Calculate matrix factor 2
            mat_factor_2 = (P*J_inv*P_inv)            
            # Calculate the homogeneous solution
            homo_sol = mat_factor_1*c_mat
            # Add the particular solution if such a solution exists
            if non_homogeneous: # non-homogeneous ODE
                # Begin by defining the integrand of the particular
                # solution as the exponential matrix times the source
                # term of the ODE
                part_sol = expand(mat_factor_2*source_ODE)                
                part_sol_der = expand(mat_factor_2*source_ODE_der)
                # Introduce a dummy variable with which we conduct the integration
                s = Symbol('s')              
                # Convert the derivative part to a list
                part_sol_der = list(part_sol_der)
                # Simplify the derivative terms by using integration by parts
                part_sol_der = integration_by_parts(non_pivot_functions,const_inhomo,part_sol_der,x[0],s)
                # Re-define these terms as a Matrix
                part_sol_der = Matrix(len(part_sol_der),1,part_sol_der)
                # Define the dimensions of the particular solution at
                # hand
                m,n = part_sol.shape
                # Cast the particular solution to te matrix format
                part_sol = Matrix(part_sol)
                # Loop over the rows and columns and integrate each
                # element of the vector with the particular solution
                for rI in range(m):
                    for cI in range(n):
                        # Replace the current variable x[0] with the dummy variable s
                        part_sol[rI,cI] = part_sol[rI,cI].subs(x[0],s)
                        # Define each matrix element as an integral 
                        # with respect to the zeroth variable 
                        part_sol[rI,cI] = Integral(part_sol[rI,cI],(s,0,x[0]))
                        # If possible, try to evaluate the integral at hand
                        part_sol[rI,cI] = part_sol[rI,cI].doit()
                # Lastly, multiply with J again to construct the particular solution
                part_sol = simplify(mat_factor_1*(part_sol + part_sol_der))
            else:# homogeneous ODE
                part_sol = zeros(len(c_mat),1) # We just add zeroes
            # Construct the solution by adding the particular
            # and homogeneous solution
            c_mat = expand(homo_sol + part_sol)
            print("\t\t\t\tODE system is solved...")
            #----------------------------------------------
            # PART 3: Solve Algebraic equations
            #----------------------------------------------
            # Derive the algebraic equations where each
            # element equals zero
            if non_homogeneous:# The non-homogeneous case
                c_alg = expand(B_algebraic*c_mat + source_alg)
            else:
                c_alg = expand(B_algebraic*c_mat)
            # Define a list of remaining constants
            const_remaining = []
            # Define the rows in c_alg
            alg_row,alg_col = c_alg.shape            
            # Define the number of unknown coefficients in this
            # system of equations
            for constant in constant_list:
                # Loop over equations
                for eq_index in range(alg_row):
                    # Lastly, add the constant if it exist
                    if c_alg[eq_index,0].coeff(constant) != 0:
                        const_remaining.append(constant)
                        break
            # Add the inhomogeneous constants as well
            const_remaining = const_remaining + const_inhomo
            # Convert to standard matrices
            c_mat = Matrix(c_mat)
            c_alg = Matrix(c_alg)
            # Loop through the algebraic equations and solve them
            for eq_temp_index in range(len(c_alg)):
                # Extract the current equation
                #eq_temp =expand(cancel(expand(c_alg[eq_temp_index])))
                eq_temp = expand(c_alg[eq_temp_index])
                c_alg[eq_temp_index] = eq_temp
                alg_str += "\nEquation:$" + latex(eq_temp) + "=0$\n"
                # Solve the corresponding equations given by the current
                # algebraic equations
                LHS_list,RHS_list,basis_functions,LHS_before,RHS_before,eq_str  = solve_algebraic_equation(eq_temp,const_remaining,x)
                alg_str += "Basis functions:\n$$" + latex(basis_functions) + "$$\n"
                #alg_str += eq_str
                #alg_str += "Solutions \\textit{before} processing:\n"
                #alg_str += "\\begin{align*}\n"
                #for temp_index in range(len(LHS_before)):
                #    alg_str += latex(LHS_before[temp_index]) + "&=" + latex(RHS_before[temp_index]) + "\\\\\n"
                #alg_str += "\\end{align*}\n"
                alg_str += "Solutions \\textit{after} processing:\n"
                alg_str += "\\begin{align*}\n"
                for temp_index in range(len(LHS_list)):
                    alg_str += latex(LHS_list[temp_index]) + "&=" + latex(RHS_list[temp_index]) + "\\\\\n"
                alg_str += "\\end{align*}\n"
                # Substitute the solution of the algebraic equation
                # into the solution of the ODE for the tangential coefficients
                for sub_index in range(len(c_mat)):
                    for index in range(len(LHS_list)):                        
                        c_mat[sub_index] = c_mat[sub_index].subs(LHS_list[index],RHS_list[index])
                        #c_mat[sub_index] = cancel(expand(c_mat[sub_index]))
                        #c_mat[sub_index] = expand(cancel(expand(c_mat[sub_index])))
                        c_mat[sub_index] = expand(c_mat[sub_index])
                # Current ODE solutions
                alg_str += "Current ODE solutions:\n"
                alg_str += "\\begin{equation*}\n\\mathbf{c}=" + latex(c_mat) + "\n\\end{equation*}\n" 
                # Substitute the solution of the current algebraic equation into the remaining
                # algebraic equations
                # Find the next index
                next_index = eq_temp_index + 1
                # If we are not at the second last algebraic equation, we shall also substitute
                # the current algebraic solution into the remaining algebraic equations
                if next_index != (len(c_alg) - 1):
                    # Loop over the remaining algebraic equations and substitute the current
                    # conditions
                    for sub_index in range(next_index,len(c_alg)):
                        for index in range(len(LHS_list)):
                            c_alg[sub_index] = c_alg[sub_index].subs(LHS_list[index],RHS_list[index])
                            #c_alg[sub_index] = expand(cancel(expand(c_alg[sub_index])))
                            c_alg[sub_index] = expand(c_alg[sub_index])
                # Current algebraic solutions
                alg_str += "Current algebraic equations solutions:\n"
                alg_str += "\\begin{equation*}\n" + latex(c_alg) + "=" + latex(zeros(len(c_alg),1)) + "\n\\end{equation*}\n"                             
            alg_str += "\\huge\\textit{Solutions after all is done:}\\normalsize\n"
            alg_str += "\\begin{equation*}\n\\mathbf{c}=" + latex(c_mat) + "\n\\end{equation*}\n" 
            #----------------------------------------------
            # PART 4: Substitute the solution into the tangents
            # and each sub-generator
            #----------------------------------------------
            # Add the inhomogenous constants to the list of constants
            constant_list = constant_list + const_inhomo
            # See which constants that exists among the constants
            constants_in_final_solution = []
            # Loop over all algebraic expressions
            for alg_index in range(len(c_mat)):
                # Loop over all candidate constants
                for constant_index in range(len(constant_list)):
                    # See if the constant exists in the final solution
                    if c_mat[alg_index].coeff(constant_list[constant_index]) != 0:
                        # Add the constant
                        constants_in_final_solution.append(constant_list[constant_index])
            # Only save the unique ones
            constants_in_final_solution = list(set(constants_in_final_solution))
            # Loop over the constants and assemble the final tangents
            for index in range(len(c)):
                eta_list_final[tangent_indicators[index]] += expand(c_mat[index]*monomial_list[index])
            # Also, loop over the non-pivot arbitrary functions and add these to the same tangents
            for index in range(len(non_pivot_functions)):
                eta_list_final[non_pivot_tangent_indicators[index]] += non_pivot_functions[index](x[0])*non_pivot_monomials[index]
            # Finally, just loop through the tangents and simplify
            alg_str += "Tangents before any manipulation:\n\\begin{align*}\n"
            for index in range(len(eta_list)):
                eta_list_final[index] = expand(eta_list_final[index])
                alg_str += "\\eta_{" + str(index) + "}&=" + latex(eta_list_final[index]) + "\\\\\n"
            alg_str += "\\end{align*}\n"
            
            # In the non-homogeneous case, we need to save a non-homogeneous tangent as well
            if non_homogeneous:
                # Allocate memory
                non_homo_tangent = []
                # Initiate with zeroes 
                for tangent_number in range(len(eta_list)):
                    non_homo_tangent.append(0)                
            # Split the tangents into its component parts
            tangent_component = []                
            # Loop over all arbitrary integration coefficients and save the tangents
            for constant in constants_in_final_solution:
                # Allocate a list for the temporary generator
                temp_eta = []
                # Loop over all tangents and add the generator
                # corresponding to the current arbitrary coefficient
                for tangent_number in range(len(eta_list)):
                    # Save the sub-generator
                    temp_eta.append(eta_list_final[tangent_number].coeff(constant))
                    # The non-homogeneous case: Save all terms involving the coefficients so that we can substract it from the original tangent later in order to calculate the extra terms.
                    if non_homogeneous:
                        # Add the term in order to calculate the non-homogeneous part
                        non_homo_tangent[tangent_number] += eta_list_final[tangent_number].coeff(constant)*constant
                # Append the extracted generator to the list of tangents
                tangent_component.append(temp_eta)
            # Calculate the non-homogeneous tangent if such a tangent exist
            if non_homogeneous:
                # Loop over the non-homogenous terms and calculate the part of the tangent that is not associated with an arbitrary integration constant
                for non_homo_index in range(len(non_homo_tangent)):
                    non_homo_tangent[non_homo_index] = expand(simplify(eta_list_final[non_homo_index] - non_homo_tangent[non_homo_index]))
                # Append the non-homogeneous tangent to the list of tangents
                tangent_component.append(non_homo_tangent)

            alg_str += "Generators before some are removed:\n"
            for tangent_index in range(len(tangent_component)):
                alg_str += "$$" + latex(tangent_component[tangent_index]) + "$$\n"
            #----------------------------------------------
            # PART 5: Check if the generators satisfy
            # the linearised symmetry conditions
            #----------------------------------------------
            # Check the linearised symmetry conditions
            # Allocate a vector of the linearised symmetry conditions
            lin_sym_index = []
            lin_sym_failure = []
            #print("# The component tangents and their linearised symmetry conditions")
            alg_str += "The tangents and their symmetry conditions:\n"
            # Loop over all sub tangents and check the linearised symmetry conditions
            for tangent_index in range(len(tangent_component)):
                alg_str += "$$\\eta_{" + str(tangent_index) + "}=" + latex(tangent_component[tangent_index]) + "$$\n"
                # Calculate the symmetry conditions for the tangent at hand
                temp_list = lin_sym_cond(x, tangent_component[tangent_index], omega_list)
                temp_list = [expand(expr) for expr in temp_list]
                alg_str += "$$" + latex(temp_list) + "$$\n"
                # Loop over all tangents in the generator at hand
                for sub_index in range(len(tangent_component[tangent_index])):
                    # Extract a tangent
                    tangent_temp = tangent_component[tangent_index][sub_index]
                    # Substitute to the original variables
                    for temp_index in range(len(x)):
                        tangent_temp = tangent_temp.subs(x[temp_index],variables[temp_index])
                # We assume that the tangent at hand is a symmetry. This would mean
                # that all symmetry conditions equal zero. 
                not_a_symmetry = False
                # Loop over the symmetry conditions
                for term in temp_list:
                    # If one of them is non-zero then
                    # we do not have a symmetry and
                    # we stop looping!
                    if term != 0:
                        not_a_symmetry = True
                        break
                # The generator in question was miscalculated, and
                # we remove it from the list of symmetries
                if not_a_symmetry:
                    lin_sym_index.append(tangent_index)
                    lin_sym_failure.append(temp_list)
            # Sort the list in reverse
            lin_sym_index.sort(reverse=True)
            if len(lin_sym_index)>0:
                # We had failed symmetries which shall be printed, hey?
                not_a_symmetry = True
            # Remove all tangents that were miscalculated
            for index in lin_sym_index:
                del tangent_component[index]
            alg_str += "Generators after some were removed:\n"
            for tangent_index in range(len(tangent_component)):
                alg_str += "$$" + latex(tangent_component[tangent_index]) + "$$\n"                
            #----------------------------------------------
            # PART 6: Printing the generator
            #----------------------------------------------
            # Substitute from the pseudo to the new variables   
            for tangent_index in range(len(tangent_component)):
                for sub_index in range(len(tangent_component[tangent_index])):
                    for index in range(len(x)):
                        tangent_component[tangent_index][sub_index] = tangent_component[tangent_index][sub_index].subs(x[index],variables[index])
            # Change the name of the arbitrary functions to something with f
            arb_func_counter = 1
            arbitrary_functions = []
            for non_pivot_function in non_pivot_functions:
                # Introduce an arbitrary function
                exec("f_%d = symbols(\'f_%d\', cls=Function) "%(arb_func_counter,arb_func_counter))
                # Save the same arbitrary function
                exec("arbitrary_functions.append(f_%d) "%(arb_func_counter))
                # Increment the counter
                arb_func_counter += 1
            # Define a counter for the generators
            generator_counter = 0
            # Initialise the output string
            X = ""
            #-----------------------------------------------------
            # Some stuff for pretty printing
            #-----------------------------------------------------
            # Define a term counter and a line counter
            term_counter = 1
            # Define two tolerances for the number of terms we tolerate
            # and the number of lines we tolerate
            term_tolerance = 6
            # Loop over the component tangent vectors
            for tangent in tangent_component:
                # Allocate memory for a list with all terms
                term_list = []
                # Reset the term_counter
                term_counter = 1
                # Increment the generator_counter
                generator_counter += 1
                # Loop over each coordinate in the tangent vector   
                for tangent_part in range(len(tangent)):
                    # Add all non-trivial tangents to the generator    
                    if tangent[tangent_part]!=0:
                        # Replace the name of all arbitrary functions
                        for arb_func_ind in range(len(non_pivot_functions)):                        
                            tangent[tangent_part] = tangent[tangent_part].subs(non_pivot_functions[arb_func_ind](variables[0]),arbitrary_functions[arb_func_ind](variables[0]))
                            tangent[tangent_part] = tangent[tangent_part].subs(non_pivot_functions[arb_func_ind](s),arbitrary_functions[arb_func_ind](s))                            

                        # Save a temporary string for the tangents
                        temp_str = latex(tangent[tangent_part])
                        # Chop the tangent into its pieces
                        chopped_generator = temp_str.split("+")
                        alg_str += "Generator:$$" + temp_str + "$$\n"
                        # Calculate the arguments
                        tangent_arguments = tangent[tangent_part].args
                        # Case one, we only have a translation generator
                        # which is identified in the following manner
                        #if len(tangent_arguments) == 0 and len(chopped_generator) == 1:
                        if len(chopped_generator) == 1:
                            # We stay with the chopped_generator which is just the
                            # translation operator for instance
                            chopped_generator = chopped_generator
                        else:
                            chopped_generator = [latex(argument) for argument in tangent_arguments]
                        alg_str += "Component parts of generator:$$" + str(chopped_generator) + "$$\n"                             
                        # Loop over the pieces of the tangent and add them
                        for term_index in range(len(chopped_generator)):
                            # If we only have one generator we add it directly
                            if len(chopped_generator)==1:
                                term_list.append("\\left(" + chopped_generator[term_index] + " \\right)" + "\\partial " + str(latex(variables[tangent_part])))
                            # For the first element we add a "\left(" before
                            elif term_index == 0:
                                term_list.append("\\left(" + chopped_generator[term_index])
                            # The last element we add a "\right)" after
                            elif term_index == (len(chopped_generator)-1):
                                term_list.append(chopped_generator[term_index] + " \\right)" + "\\partial " + str(latex(variables[tangent_part])))
                                # Also add an extra element so that we can recognise these elements
                                term_list.append(1)                                
                            # The other terms in the middle of the tangent is just about adding
                            # the darn element
                            else:
                                term_list.append(chopped_generator[term_index])
                # Back to the tangent! Let's loop over it and add the generator at hand
                # Start the alignment
                X += "\\begin{align*}\nX_{" + str(generator_counter) + "}=&"
                # Loop over the terms
                for term_index in range(len(term_list)):
                    # If the element is one we move on
                    if term_list[term_index] == 1:
                        continue
                    else: # Else, we mean business!
                        # Add the term at hand
                        X += term_list[term_index]
                        # Increment the term_counter
                        term_counter += 1
                        # Check the cases when we should add a plus sign and so on
                        if term_index == (len(term_list)-1): # The very last element
                            # The very last generator ends with a point
                            if generator_counter == len(tangent_component):
                                X += ".\\\\\n"
                            else: # The other generators ends with a comma
                                X += ",\\\\\n"                            
                        else: # The other elements
                                # If we have too many terms we break
                                # and add a new line
                                if term_counter > term_tolerance:
                                    # Re-set the term_counter
                                    term_counter = 1
                                    # We are at the end of a tangent coordinate
                                    if term_list[term_index + 1]==1:
                                        # We are at the end of the tangent and
                                        # therefore we should close the generator
                                        if (term_index + 1) == (len(term_list)-1):
                                            X+= "\\\\\n"
                                        else: # Otherwise, we continue adding terms on the next row
                                            X+= "\\\\\n&+"
                                    else: # We break in the middle of a tangent
                                        X+= "\\right.\\\\\n&+\\left."
                                else: # We just add a plus sign
                                    X += "+"
                # We finish up with our tangent
                X += "\\end{align*}\n\n"
                # Correct the ending of the alignment
                X = X.replace("+\\end{align*}","\n\\end{align*}")
                # Correct the minus signs as well
                X = X.replace("+-","-")                
        # Add a string in case we have a non-homogeneous function
        if non_homogeneous:
            # Add all the arbitrary functions to the output
            temp_str =  "\n\n\\noindent Some of the generators might contain the following arbitrary functions:\n"
            temp_str += "\\begin{align*}\n"
            for arbitrary_function in arbitrary_functions:
                temp_str += "&" + latex(arbitrary_function) + "\\\\\n"
            temp_str += "\\end{align*}\n\n"
            # Add these to the generator as well
            X = X + temp_str
        # Add a string in case we have non-generators             
        if not_a_symmetry:
            # Add a string saying that something went wrong and that the list is not complete
            temp_str = "\\noindent\\huge\\textbf{WARNING}:\\\\\n"
            temp_str += "\\noindent\\Large\\textit{Some of the calculated generators did not satisfy the linearised symmetry conditions. Thus, the presented list here is not complete and consists exclusively of the calculated generators that satisfy the linearised symmetry conditions.}\\normalsize\\\\[2cm]\n"
            # Add these to the generator as well
            X = X + temp_str
    else:
        # Return that the matrix is not quadratic
        X = "\\Huge\\textsf{Not a quadratic matrix}\\normalsize\\\\[2cm]" 
    #X += alg_str
    # Return the solved coefficients and the generator
    return X, c_mat, c_alg, c_original, eta_list_final 


    
