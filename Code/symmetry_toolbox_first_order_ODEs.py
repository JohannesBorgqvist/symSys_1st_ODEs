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
#from symengine import *
#import sympy
from sympy import *
from sympy.core.function import *
#from sympy import simplify
#from sympy import fraction
#from sympy import powsimp
#from sympy import cancel
#from sympy import collect
#from sympy import Matrix
#from sympy import solve
#from sympy import latex
#from sympy.core.function import *
# Other functionalities in sympy:
# Finding monomials
from sympy.polys.monomials import itermonomials, monomial_count
# Ordering monomials 
from sympy.polys.orderings import monomial_key
# Import some variables just to test
from sympy.abc import x, y
# To solve ODE systems
from sympy import symbols, Eq, Function
# Import all matrix stuff
from sympy.matrices import Matrix, eye, zeros, ones, diag, GramSchmidt
# For printing with latex
from sympy import latex
# Import math for combinatorials
from math import *
# To do the fancy iterations
import itertools
# To manipulate string
import string
# To time each part of the program
import time
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
#6. *sub_const* (The function changes the name of the arbitrary integration constants)
#7. *solve_linear_sys_ODEs* (Help function 1 for function 9)
#8. *solve_system_with_algebraic_equations* (Help function 2 for function 9)
#9. *solve_determining_equations* (Solves the determining equations).
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
    return x, c, eta_list
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
        lin_sym_list[help_counter] = powsimp(expand(cancel(tidy_eq)))
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
# The function takes two inputs which are the following:
# 1. An arbitrary expression being a linear combination of functions stored in the variable "expr",
# 2. A list of all arbitrary integration coefficients stored in "coeff_list",
# The script returns two outputs:
# 1. A list of arbitrary functions stored in "arbitrary_functions",
# 2. A list of all basis functions stored in "basis_functions".
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
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# FUNCTION 7: "solve_algebraic_equation"
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# The function takes two inputs which are the following:
# 1. An arbitrary expression being a linear combination of functions stored in the variable "expr",
# 2. A list of all arbitrary integration coefficients stored in "coeff_list",
# The script returns two outputs:
# 1. A list of the LHSs being the solutions stored in "LHS",
# 2. A list of the RHSs being the expressions stored in "RHS".
# The function uses the previous function (i.e. Function 6) called "identify_basis_functions" in
# order to find the basis functions and the arbitrary functions in the expression at hand.
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
                        #LHS_temp = eq_temp.coeff(koeff)*koeff
                        LHS_temp = koeff
                        break
                # Solve the equation for this coefficient
                RHS_temp = solve(eq_temp,LHS_temp)[0]
                #print("\t\t\tSolve\t%s\tfor\t%s\tyields\t%s"%(str(eq_temp),str(LHS_temp),str(RHS_temp)))
                # Append both the LHS and the RHS
                LHS.append(LHS_temp)
                RHS.append(RHS_temp)                
                # Increase the temporary sum
                temp_sum += eq_temp*func_temp
        # Define the zeroth equation
        eq_zero = simplify(cancel(expand(expr - temp_sum)))
        # Find the coefficient which we solve for
        LHS_temp = 0
        for koeff in coeff_list:
            if eq_zero.coeff(koeff)!=0:
                #LHS_temp = eq_zero.coeff(koeff)*koeff
                LHS_temp = koeff
                break
        # Solve the equation for this coefficient
        RHS_temp = solve(eq_zero,LHS_temp)[0]
        #print("\t\t\tSolve\t%s\tfor\t%s\tyields\t%s"%(str(eq_zero),str(LHS_temp),str(RHS_temp)))
        # Append both the LHS and the RHS
        LHS.append(LHS_temp)
        RHS.append(RHS_temp)                        
    # Lastly, return the LHS and the RHS
    return LHS, RHS
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# FUNCTION 8: "integration_by_parts"
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
# 1. The "integral_list" being a list of all integrals.
def integration_by_parts(function_list,constant_list,integrand_list,variable,s):
    # At the end, we want to return the correct integrals
    integral_list = []
    # Loop through the integrands
    for integrand in integrand_list:
        # Allocate a temporary integral
        temp_int = 0
        # Split this into further parts and solve these
        for sub_integrand in integrand.args:
            # Now, firstly we isolate the factors
            factors_in_subintegrand = sub_integrand.factor().args
            # The derivative term
            derivative_term = 0
            # Loop over our factors in the integrand 
            # and see if any of them contain a derivative 
            for factor in factors_in_subintegrand:
                # Find the arbitrary functions
                a_f = list(factor.atoms(AppliedUndef))
                if len(a_f)>0:
                    derivative_term = factor
            # The other func is the coefficient of the derivative term
            other_func = sub_integrand.coeff(derivative_term)
            # Calculate the primitive function
            primitive_function = Integral(derivative_term,variable).doit()
            # Let's add the three terms to our temporary integral
            temp_int += primitive_function*other_func
            #temp_int += -primitive_function*other_func
            # Loop over integration constants as well
            for f_i in range(len(function_list)):
                if function_list[f_i](variable)==primitive_function:
                    temp_int += -constant_list[f_i]*other_func.subs(variable,0)
                    #temp_int += constant_list[f_i]*other_func.subs(variable,0)
            # Lastly, we add the integral term:
            # Define the new integrand
            new_integrand = primitive_function*Derivative(other_func,variable).doit()
            # Change integration to a dummy variable s
            new_integrand = new_integrand.subs(variable,s)
            # Lastly, add the integral term as well
            temp_int += -Integral(new_integrand,(s,0,variable))
            #temp_int += Integral(new_integrand,(variable,0,variable))
        # Make sure that the integral is evaluated 
        if temp_int != 0:
            temp_int = temp_int.doit()
        # Add the evaluated integrals to our integral list
        integral_list.append(temp_int)
    # Return the integral list
    return integral_list
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# FUNCTION 9: "solve_determining_equations"
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# The function takes four inputs which are the following:
# 1. The variables (dependent and independent) in a vector "x",
# 2. The tangent vector called "eta_list",
# 3. The coefficients in the tangents called "c"
# 4. The list of the determining equations called "det_eq".
# The script return the following output:
# 1. The calculated generator of the symmetry which is stored in a string called "X".
def solve_determining_equations(x,eta_list,c,det_eq,variables,omega_list):
    #print("# Solving the determining equations")
    #------------------------------------------------------------------------------
    # STEP 1 of 6: FORMULATE THE DETERMINING EQUATIONS ON MATRIX FORM A*dc/dt=B*c
    #-----------------------------------------------------------------------------
    #print("## Step 0 of 6: Defining the matrices from the determining equations")
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
        #print("### Determining equation %d out of %d"%(help_counter,len(det_eq)))
        #help_counter+=1
        #eq_temp_copy = eq_temp
        #for index in range(len(x)):
            #eq_temp_copy = eq_temp_copy.subs(x[index],variables[index])
        #print("\\begin{equation}\n%s=0\n\\end{equation}"%(latex(eq_temp_copy)))
        
        # Loop over all coefficients in the tangential ansätze
        # and define our matrices
        for koeff in c:
            # Save the derivative of the A coefficient
            A_mat.append(-eq_temp.coeff(Derivative(koeff(x[0]),x[0])))
            #A_row.append(-eq_temp.coeff(Derivative(koeff(x[0]),x[0])))
            # Save the coefficient of the coefficient in the tangential
            # ansätze
            B_mat.append(eq_temp.coeff(koeff(x[0])))
            #B_row.append(eq_temp.coeff(koeff(x[0])))
        #A_row = Matrix(1,len(c),A_row)
        #B_row = Matrix(1,len(c),B_row)            
        #for index in range(len(x)):
            #A_row = A_row.subs(x[index],variables[index])
            #B_row = B_row.subs(x[index],variables[index])
        #print("\\begin{equation}\n%s\\dfrac{\\mathrm{d}}{\\mathrm{d}t}%s=%s%s\n\\end{equation}"%(latex(A_row),latex(Matrix(len(c),1,c)),latex(B_row),latex(Matrix(len(c),1,c))))
        #A_row = []
        #B_row = []
    # Finally we generate the matrices from our lists A_mat and B_mat
    num_of_rows = len(det_eq) # Number of rows
    num_of_cols = len(c) # Number of columns
    A = Matrix(num_of_rows,num_of_cols,A_mat) # Matrix A
    B = Matrix(num_of_rows,num_of_cols,B_mat) # Matrix B
    #print("## Step 1 of 6: the initial matrices")
    #print("Dimension of matrices:\t%dX%d<br>"%(int(num_of_rows),int(num_of_cols)))
    #print("Matrix A<br>")
    #A_temp = A
    #print("Matrix B<br>")
    #B_temp = B
    #for index in range(len(x)):
        #A_temp = A_temp.subs(x[index],variables[index])
        #B_temp = B_temp.subs(x[index],variables[index])
    #print("\\begin{equation}\nA=%s\n\\end{equation}"%(latex(A_temp)))
    #print("\\begin{equation}\nB=%s\n\\end{equation}"%(latex(B_temp)))    
    
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
    #A = -M_tilde[:,0:(len(c))]
    #B = M_tilde[:,len(c):(2*len(c))]
    A = M_tilde[:,0:(len(c))]
    B = -M_tilde[:,len(c):(2*len(c))]    
    #print("## Step 2 of 6: the reduced based on col(M^T) where M=[-A|B]<br>")
    #r_temp, c_temp = A.shape
    #print("Dimension of matrices:\t%dX%d<br>"%(int(r_temp),int(c_temp)))
    #A_temp = A
    #B_temp = B
    #for index in range(len(x)):
        #A_temp = A_temp.subs(x[index],variables[index])
        #B_temp = B_temp.subs(x[index],variables[index])
    #print("Matrix A<br>")        
    #print("\\begin{equation}\nA=%s\n\\end{equation}"%(latex(A_temp)))
    #print("Matrix B<br>")
    #print("\\begin{equation}\nB=%s\n\\end{equation}"%(latex(B_temp)))     
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
    #print("## Step 3 of 6: Splitting up to A, B and B_algebraic")
    #r_temp, c_temp = A.shape
    #r_temp_2, c_temp_2 = B_algebraic.shape    
    #print("Dimension of matrices A and B:\t%dX%d<br>"%(int(r_temp),int(c_temp)))
    #print("Dimension of matrices B_algebraic:\t%dX%d<br>"%(int(r_temp_2),int(c_temp_2)))    
    #A_temp = A
    #B_temp = B
    #B_alg_temp = B_algebraic
    #for index in range(len(x)):
    #    A_temp = A_temp.subs(x[index],variables[index])
    #    B_temp = B_temp.subs(x[index],variables[index])
    #    B_alg_temp = B_alg_temp.subs(x[index],variables[index])        
    #print("Matrix A<br>")        
    #print("\\begin{equation}\nA=%s\n\\end{equation}"%(latex(A_temp)))
    #print("Matrix B<br>")
    #print("\\begin{equation}\nB=%s\n\\end{equation}"%(latex(B_temp)))
    #print("Matrix B_algebraic<br>")
    #print("\\begin{equation}\nB_{\\textrm{algebraic}}=%s\n\\end{equation}"%(latex(B_alg_temp)))
    #print("Coefficient matrix c:")
    #print("\\begin{equation}\n\\mathbf{c}=%s\n\\end{equation}"%(latex(Matrix(len(c),1,c))))
    #------------------------------------------------------------------------------
    # STEP 4 of 6: REMOVE (POTENTIAL) EXTRA COLUMNS 
    #-----------------------------------------------------------------------------
    # Calculate rows and columns of matrices
    r_algebraic, cols = B_algebraic.shape # Dimensions of B_algebraic
    rows, cols = A.shape # Dimensions of A
    #print("Dimensions of A:\t%dX%d\n"%(rows,cols))
    #print("Dimensions of B_algebraic:\t%dX%d\n"%(r_algebraic,cols))
    # Create a set of all the columns
    col_set = set(range(cols))
    # Find the pivot columns of A
    pivot_columns = A.rref()[1]
    # Calculate the non-pivot columns
    non_pivot_columns = list(col_set.difference(set(pivot_columns)))
    # Substitute the value zero for all non-pivot columns in the coefficient
    # vector
    # Loop over tangents
    #for tangent_number in range(len(eta_list)):
        # Loop through all coefficients and replace them in all tangents
        #for index in non_pivot_columns:            
            # Substitute coefficient in the tangents
            #eta_list[tangent_number] = eta_list[tangent_number].subs(c[index](x[0]),0)
    # We assume that we have a homogeneous system
    non_homogeneous = False
    # Check if we have non-pivot columns, that is we have
    # an underdetermined system with extra columns (i.e. more columns than equations).
    # In this case, we move these column to the right hand side and we then treat
    # these as source terms if you will, so that we get a quadratic system with some
    # extra inhomogeneities. 
    if len(non_pivot_columns)!=0: # Yes, we do have non-pivot coefficients
        # Print the non-pivot columns
        #print("Non-Pivot columns")
        #print(non_pivot_columns)
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
        # Print the matrix
        #print("The non-homogeneous source term:<br>")
        #print(latex(Matrix(len(source_ODE),1,source_ODE),mode='equation*'))
        source_ODE = Matrix(len(source_ODE),1,source_ODE)
        source_ODE_der = Matrix(len(source_ODE_der),1,source_ODE_der)        
        #print("Dimensions:<br>")
        #print(source_ODE.shape)
        #print("The algebraic source term:<br>")
        source_alg = Matrix(len(source_algebraic),1,source_algebraic)
        #print(latex(source_alg,mode='equation*'))
        #print("Dimensions:<br>")
        #print(source_alg.shape)
    # Save the unknown arbitrary functions being the non-pivot elements
    non_pivot_functions = [c[non_pivot_column] for non_pivot_column in non_pivot_columns]
    # Remove all non-pivot tangential coefficient
    non_pivot_columns.sort(reverse=True)
    # Loop over all non_pivot_columns and remove them from our three
    # matrices A, B and B_algebraic and from the coefficient vector c
    for index in range(len(non_pivot_columns)):
        # Remove column matrices
        A.col_del(non_pivot_columns[index])
        B.col_del(non_pivot_columns[index])
        B_algebraic.col_del(non_pivot_columns[index])
        # Remove the non-pivot parameters from the coefficient vector
        del c[non_pivot_columns[index]]
    #print("## Step 4 of 6: Removing potential extra pivot columns")
    r_temp, c_temp = A.shape
    r_temp_2, c_temp_2 = B_algebraic.shape    
    #print("Dimension of matrices A and B:\t%dX%d<br>"%(int(r_temp),int(c_temp)))
    #print("Dimension of matrices B_algebraic:\t%dX%d<br>"%(int(r_temp_2),int(c_temp_2)))    
    #A_temp = A
    #B_temp = B
    #B_alg_temp = B_algebraic
    #for index in range(len(x)):
        #A_temp = A_temp.subs(x[index],variables[index])
        #B_temp = B_temp.subs(x[index],variables[index])
        #B_alg_temp = B_alg_temp.subs(x[index],variables[index])        
    #print("Matrix A<br>")        
    #print("\\begin{equation}\nA=%s\n\\end{equation}"%(latex(A_temp)))
    #print("Matrix B<br>")
    #print("\\begin{equation}\nB=%s\n\\end{equation}"%(latex(B_temp)))
    #------------------------------------------------------------------------------
    # STEP 5 of 6: SOLVE THE ODE SYSTEM PROVIDED THAT IT IS QUADRATIC AND
    # SUBSTITUTE THE ALGEBRAIC EQUATIONS
    #-----------------------------------------------------------------------------
    num_rows, num_cols = B.shape
    #print("Dimensions of B:\t%dX%d<br>"%(num_rows,num_cols))
    if num_rows == num_cols:
        # Check that the matrix A is an identity matrix
        #if A != -eye(num_rows):
        #if abs(A) != eye(num_rows):
        if A != eye(num_rows):
            X = "\\Huge\\textsf{$A$ is quadratic but not an identity matrix!}\normalsize\\[2cm]" 
        else: # A is an identity matrix
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
            #print("Coefficients:<br>")
            #print(latex(c_mat,mode='equation'))
            #print("Number of unknowns:\t%d<br>"%(int(len(c_mat))))
            # PRINT TO THE USER
            #print("## Step 5 of 6: Solving the ODE system")
            #r_temp, c_temp = B.shape
            #r_temp_2, c_temp_2 = B_algebraic.shape    
            #print("Dimension of the matrix B:\t%dX%d<br>"%(int(r_temp),int(c_temp)))
            #print("Dimension of the matrix B_algebraic:\t%dX%d<br>"%(int(r_temp_2),int(c_temp_2)))    
            #B_temp = B
            #B_alg_temp = B_algebraic
            #for index in range(len(x)):
                #A_temp = A_temp.subs(x[index],variables[index])
                #B_temp = B_temp.subs(x[index],variables[index])
                #B_alg_temp = B_alg_temp.subs(x[index],variables[index])        
            #print("ODE system:<br>")        
            #print("\\begin{equation}\n\dfrac{\mathrm{d}}{\mathrm{d}t}%s=%s%s\n\\end{equation}"%(latex(Matrix(len(c),1,c)),latex(B_temp),latex(Matrix(len(c),1,c))))     
            #----------------------------------------------
            # PART 2: Solve ODE system
            #----------------------------------------------
            #print("Solve the ODE system:<br>")
            #print("Initial conditions for $\\mathbf{c}$ denoted by $\mathbf{c}_0$ in terms of arbitrary integration constants:")
            #print("\\begin{equation}\n\mathbf{c}_{0}=%s=%s\n\\end{equation}"%(latex(Matrix(len(c),1,c)),latex(c_mat)))
            # Calculate the Jordan form of the matrix B
            P,J = B.jordan_form()
            #P,J = (-B).jordan_form()
            #print("Jordan form:<br>")
            #J_temp = J
            #P_temp = P
            #for index in range(len(x)):
                #J_temp = J_temp.subs(x[index],variables[index])
                #P_temp = P_temp.subs(x[index],variables[index])
            #print("\\begin{equation}\nP=%s\n\\end{equation}"%(latex(P_temp)))
            #print("\\begin{equation}\nJ=%s\n\\end{equation}"%(latex(J_temp)))                            
            # Re-define our matrix J by scaling it by the
            # independent variable x[0]
            J = x[0]*J
            # Re-define J by taking the matrix exponential
            J = J.exp()
            #J_temp = J
            #for index in range(len(x)):
                #J_temp = J_temp.subs(x[index],variables[index])
            #print("Exponential form:<br>")
            #print("\\begin{equation}\n\\exp\\left(J\\cdot t\\right)=%s\n\\end{equation}"%(latex(J_temp)))
            # Calculate the homogeneous solution
            homo_sol = (P*J*P.inv())*c_mat
            #homo_sol_temp = homo_sol
            #for index in range(len(x)):
                #homo_sol_temp = homo_sol_temp.subs(x[index],variables[index])            
            #print("Solution to the ODE system:")
            #print("\\begin{equation}\nP\\exp\\left(J\\cdot t\\right)P^{-1}\mathbf{c}_{0}=%s\n\\end{equation}"%(latex(homo_sol_temp)))                 
            # Add the particular solution
            if non_homogeneous: # non-homogeneous ODE
                # Begin by defining the integrand of the particular
                # solution as the exponential matrix times the source
                # term of the ODE
                part_sol = (P*J.inv()*P.inv())*source_ODE
                part_sol_der = (P*J.inv()*P.inv())*source_ODE_der
                part_sol_der_before = part_sol_der
                # Introduce a dummy variable with which we conduct the integration
                s = Symbol('s')
                # Convert the derivative part to a list
                part_sol_der = list(part_sol_der)
                # Simplify the derivative terms using integration by parts
                part_sol_der = integration_by_parts(non_pivot_functions,const_inhomo,part_sol_der,x[0],s)
                # Re-define it as a Matrix
                part_sol_der = Matrix(len(part_sol_der),1,part_sol_der)
                # Define the dimensions of the particular solution at
                # hand
                m,n = part_sol.shape
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
                        # Replace the current variable x[0] with the dummy variable s
                        part_sol_der_before[rI,cI] = part_sol_der_before[rI,cI].subs(x[0],s)
                        # Define each matrix element as an integral 
                        # with respect to the zeroth variable 
                        part_sol_der_before[rI,cI] = Integral(part_sol_der_before[rI,cI],(s,0,x[0]))
                        # If possible, try to evaluate the integral at hand
                        part_sol_der_before[rI,cI] = part_sol_der_before[rI,cI].doit()
                # Print before and after
                print("# Comparison with and without integration by parts")
                part_sol_der_print = part_sol_der
                part_sol_der_before_print = part_sol_der_before
                for rI in range(m):
                    for cI in range(n):
                        for index in range(len(x)):
                            part_sol_der_print[rI,cI].subs(x[index],variables[index])
                            part_sol_der_before_print[rI,cI] = part_sol_der_before_print[rI,cI].subs(x[index],variables[index])
                print("Solution without integration by parts:<br>")
                print("\\begin{equation}\nP^{-1}\\exp(Jt)P\\int^{t}_{0}P^{-1}\\exp(-Js)P\\dfrac{\\mathrm{d}c}{\\mathrm{d}s}\\mathrm{d}s=%s\n\\end{equation}"%(latex(part_sol_der_before_print)))
                print("Solution with integration by parts:<br>")
                print("\\begin{equation}\nP^{-1}\\exp(Jt)P\\int^{t}_{0}P^{-1}\\exp(-Js)P\\dfrac{\\mathrm{d}c}{\\mathrm{d}s}\\mathrm{d}s=%s\n\\end{equation}"%(latex(part_sol_der_print)))
                # Lastly, multiply with J again to construct the particular solution
                part_sol = simplify(cancel(expand((P*J*P.inv())*(part_sol + part_sol_der))))
            else:# homogeneous ODE
                part_sol = zeros(len(c_mat),1) # We just add zeroes
            # Construct the solution by adding the particular
            # and homogeneous solution
            c_mat = homo_sol + part_sol
            #print("\n\n\n")
            #print("Coefficients after ODE:<br>")
            #print(latex(c_mat,mode='equation'))
            #print("Number of unknowns:\t%d<br>"%(int(len(c_mat))))
            #----------------------------------------------
            # PART 3: Solve Algebraic equations
            #----------------------------------------------
            # Derive the algebraic equations where each
            # element equals zero
            #print("## Step 6 of 6: Solving the algebraic system<br>")
            #m, n = B_algebraic.shape
            #print("Number of algebraic equations:\t%d\n"%(m))
            #print("Matrix B_algebraic<br>")
            #print("\\begin{equation}\nB_{\\textrm{algebraic}}=%s\n\\end{equation}"%(latex(B_alg_temp)))
            #print("Algebraic equations:<br>")
            #print("\\begin{equation}\n%s%s=%s\n\\end{equation}"%(latex(B_alg_temp),latex(Matrix(len(c),1,c)),latex(zeros(len(c),1))))
            if non_homogeneous:
                c_alg = B_algebraic*c_mat + source_alg
            else:
                c_alg = B_algebraic*c_mat
            #c_alg_temp = c_alg
            #for index in range(len(x)):
                #c_alg_temp = c_alg_temp.subs(x[index],variables[index])                        
            #print("Algebraic equations after substitution of the solution to the ODE system:<br>")
            #print("\\begin{equation}\n%s%s=%s=%s\n\\end{equation}"%(str(latex(B_alg_temp)),str(latex(homo_sol_temp)),str(latex(c_alg_temp)),latex(zeros(m,1))))
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
            #print("<br>*Solving the algebraic equations:*<br>")
            # Loop through the algebraic equations and solve them
            for eq_temp in list(c_alg):
                # Solve the corresponding equations given by the current
                # algebraic equations
                #eq_print = eq_temp
                #for index in range(len(x)):
                    #eq_print = eq_print.subs(x[index],variables[index])
                #print("Equation:<br>\n\\begin{equation}\n%s=0\n\\end{equation}"%(latex(eq_print)))
                LHS_list,RHS_list = solve_algebraic_equation(eq_temp,const_remaining)
                #print("Solutions:<br>\n\\begin{align*}")
                #for index in range(len(LHS_list)):
                    #LHS_print = LHS_list[index]
                    #RHS_print = RHS_list[index]
                    #for index in range(len(x)):
                        #LHS_print = LHS_print.subs(x[index],variables[index])
                        #RHS_print = RHS_print.subs(x[index],variables[index])                                            
                    #print("%s&=%s\\\\"%(latex(LHS_print),latex(RHS_print)))
                #print("\\end{align*}")
                # Substitute the solution of the algebraic equation
                # into the solution of the ODE for the tangential coefficients
                for sub_index in range(len(c_mat)):
                    for index in range(len(LHS_list)):                    
                        c_mat[sub_index] = c_mat[sub_index].subs(LHS_list[index],RHS_list[index])
                        c_mat[sub_index] = expand(cancel(c_mat[sub_index]))
            # Create a readable version in original variables
            c_mat_temp = c_mat
            # Substitute to the original variables
            for index in range(len(x)):
                c_mat_temp = c_mat_temp.subs(x[index],variables[index])            
            #print("Solution *after* algebraic substitution:<br>")
            #print("\\begin{equation*}\n\mathbf{c}=%s\n\\end{equation*}"%(latex(c_mat_temp)))                
            #----------------------------------------------
            # PART 4: Substitute the solution into the tangents
            # and each sub-generator
            #----------------------------------------------
            print("\n\n# The very final step")
            print("The very final step: substituting the solution into the tangents and print the results:<br>")
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
            print("Arbitrary integration constants in the final solution:<br>")
            print(latex(Matrix(len(constants_in_final_solution),1,constants_in_final_solution),mode='equation*'))
            print("Number of generators which are divided based on the number of constants:\t%d<br>\n"%(len(constants_in_final_solution)))
            # Loop over tangents and substitute the calculated coefficients into these expressions
            for tangent_number in range(len(eta_list)):
                # Expand the tangents
                eta_list[tangent_number] = eta_list[tangent_number].expand()
                # Loop through all coefficients and replace them in all tangents
                for index in range(len(c)):            
                    # Substitute coefficient in the tangents
                    eta_list[tangent_number] = eta_list[tangent_number].subs(c[index](x[0]),c_mat[index])                    
                    eta_list[tangent_number] = expand(cancel(eta_list[tangent_number]))            
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
                    temp_eta.append(eta_list[tangent_number].coeff(constant))
                    # The non-homogeneous case: Save all terms involving the coefficients so that we can substract it from the original tangent later in order to calculate the extra terms.
                    if non_homogeneous:
                        # Add the term in order to calculate the non-homogeneous part
                        non_homo_tangent[tangent_number] += eta_list[tangent_number].coeff(constant)*constant
                # Append the extracted generator to the list of tangents
                tangent_component.append(temp_eta)
            # Calculate the non-homogeneous tangent if such a tangent exist
            if non_homogeneous:
                # Loop over the non-homogenous terms and calculate the part of the tangent that is not associated with an arbitrary integration constant
                for non_homo_index in range(len(non_homo_tangent)):
                    non_homo_tangent[non_homo_index] = expand(eta_list[non_homo_index]) - expand(non_homo_tangent[non_homo_index])
                    non_homo_tangent[non_homo_index] = simplify(cancel(non_homo_tangent[non_homo_index]))
                # Append the non-homogeneous tangent to the list of tangents
                tangent_component.append(non_homo_tangent)
            #----------------------------------------------
            # PART 5: Check if the generators satisfy
            # the linearised symmetry conditions
            #----------------------------------------------
            #print("Number of component tangents before removing:\t%d<br>"%(int(len(tangent_component))))
            # Check the linearised symmetry conditions
            # Allocate a vector of the linearised symmetry conditions
            lin_sym_index = []
            # Loop over all sub tangents and check the linearised symmetry conditions
            for tangent_index in range(len(tangent_component)):
                # Calculate the symmetry conditions for the tangent at hand
                temp_list = lin_sym_cond(x, tangent_component[tangent_index], omega_list)
                # Print for the viewers pleasure
                print("Generator %d out of %d:<br>"%(tangent_index+1,int(len(tangent_component))))
                # Loop over all tangents in the generator at hand
                print("\n\\begin{align*}")
                for sub_index in range(len(tangent_component[tangent_index])):
                    # Print for the viewers pleasure
                    tangent_temp = tangent_component[tangent_index][sub_index]
                    # Substitute to the original variables
                    for temp_index in range(len(x)):
                        tangent_temp = tangent_temp.subs(x[temp_index],variables[temp_index])
                    if sub_index == 0:
                        print("\\xi&=%s\\\\"%(latex(tangent_temp)))
                    else:
                        print("\\eta_%d&=%s\\\\"%(sub_index,latex(tangent_temp)))                        
                    #print(latex(tangent_component[tangent_index],mode='equation*'))
                print("\\end{align*}\n")
                print("Checking the %d linearised symmetry conditions of generator $X_%d$:<br>"%(int(len(temp_list)),tangent_index+1,))
                print("Lin syms<br>")
                print(latex(temp_list,mode='equation*'))
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
            # Sort the list in reverse
            lin_sym_index.sort(reverse=True)
            # Remove all tangents that were miscalculated
            for index in lin_sym_index:
                del tangent_component[index]
            print("Number of component tangents after removing:\t%d"%(int(len(tangent_component))))
            #----------------------------------------------
            # PART 6: Printing the generator
            #----------------------------------------------
            # Substitute from the pseudo to the new variables   
            for tangent_index in range(len(tangent_component)):
                for sub_index in range(len(tangent_component[tangent_index])):
                    for index in range(len(x)):
                        tangent_component[tangent_index][sub_index] = tangent_component[tangent_index][sub_index].subs(x[index],variables[index])
            # Define a counter for the generators
            generator_counter = 0
            # Define one generator per constant in the final solution
            #for constant in constants_in_final_solution:
            for tangent in tangent_component:
                # Increment the counter for the generator
                generator_counter += 1
                # Set the plus indicator to false
                plus_indicator = False
                # It is our first generator we begin align
                if generator_counter == 1:
                    # Define our generator
                    X = "\\begin{align*}\nX_{" + str(generator_counter) + "}&="
                else:
                    # Define our generator
                    X += ",\\\\\nX_{" + str(generator_counter) + "}&="                    
                for tangent_part in range(len(tangent)):
                    # Add all non-trivial tangents to the generator    
                    if tangent[tangent_part]!=0:
                        if plus_indicator:
                            X += "+\\left( " + str(latex(tangent[tangent_part])) + " \\right)" + "\\partial " + str(latex(variables[tangent_part])) + ""
                        else:
                            X += "\\left( " + str(latex(tangent[tangent_part])) + " \\right)" + "\\partial " + str(latex(variables[tangent_part])) + ""
                            plus_indicator = True
                # Add a line break
                X += ",\\\\\n"
            # Add the final touch to the generator
            X += "\\end{align*}"
            X = X.replace(",\\\\\n\\end{align*}", ".\\\\\n\\end{align*}")
    else:
        # Return that the matrix
        X = "\\Huge\\textsf{Not a quadratic matrix}\\normalsize\\\\[2cm]" 

    print("<br>The final generators are given by:<br>")
    print("%s"%(X))
    # Return the solved coefficients and the generator
    return X 

