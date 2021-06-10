#=================================================================================
#=================================================================================
# Script:"symmetry_toolbox_first_order_ODEs"
# Date: 2021-06-01
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
# Import some variables just to test
from sympy.abc import x, y
# To solve ODE systems
from sympy import symbols, Eq, Function
# Import math for combinatorials
from math import *
# To do the fancy iterations
import itertools
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
# FUNCTION 6: "sub_const"
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# The function takes three inputs:
# 1. A string called "expr" which contains the expression where we wish to do the substitution,
# 2. A character with the name of the arbitrary integration constant that we want the expression to have. In this application the character char = 'K' is implemented,
# 3. An integer "coefficent_counter" which corresponds to the index of the arbitrary integration coefficient.
# In the application, a certain number of integrations are performed which introduced one arbitrary integration constant. In sympy, the standard name of this integration constant is C1 (if two integrations are performed, the next one is called C2). Since multiple single integrations are performed succesively, we wish to keep track of all integration coefficients (so that not all of them are named C1) and therefore we have an integer called "coefficient_counter" corresponding to the index of the latest integration constant introduced. This function renames expression given the value of the character "char" and the integer "coefficient_counter".
def sub_const(expr, char,coefficient_counter):
    # We start with looking for the first arbitrary
    # integration constant called C1.
    i = 1
    # Rename all arbitrary integration cosntants with the
    # name C. For example, rename C1 to K"coefficient_counter".
    while True:
        old = symbols('C{}'.format(i)) # The old integration constant
        new = symbols(char + '_{}'.format(coefficient_counter + (i-1)))# The new integration constant
        new_expr = expr.subs(old, new)# The new expression
        # If we have replaced all constants we exit
        if new_expr == expr:
            return expr
        # Assign the new expression to the old one
        # before we continue re-naming the integration coeffficients,
        expr = new_expr
        # Increment our to counters
        i += 1 # For the coefficient C
        coefficient_counter += 1 # For our new coefficients
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# FUNCTION 7: "solve_linear_sys_ODEs"
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# The function is a help function to Function 9 called "solve_determining_equations". It solves the remaining equations in the ODE-system for the coefficients in the tangential ansätze. If one has algebraic equations in addition to the differential equations, the function converts the matrix system into one list containing a system of differential equations and another list containing the algebraic equations. Then the algebraic equations are substituted into the system of differential equations. Hopefully, one will end up with a situation where the system of differential equations contains single ODEs which can be integrated directly. Then, the solution to this equation is substituted into the remaining system and hopefully this will trickle down so that it becomes a matter of solving single first order ODEs.
def solve_linear_sys_ODEs(eq_sys,x,c_original,c,c_reduced,coefficient_counter,variables):
    #------------------------------------------------------------------------------
    # STEP 1: WRITE THE SYSTEM ON MATRIX FORM
    #------------------------------------------------------------------------------
    # Allocate memory for our matrices
    A_mat = []
    B_mat = []
    # Loop over all determining equations
    for isolated_eq in eq_sys:
        # Expand the equation if needed
        temp_eq = expand(isolated_eq)
        # Loop over the coefficients and construct the lists
        # which contains the elements of the matrices
        for coeff in c_reduced:
            A_mat.append(-temp_eq.coeff(Derivative(coeff(x[0]),x[0])))
            B_mat.append(temp_eq.coeff(coeff(x[0])))
    # Dimensions of our matrices
    num_of_rows = len(eq_sys) # Number of rows
    num_of_cols = len(c_reduced) # Number of columns
    # Construct the matrices from our lists    
    A = Matrix(num_of_rows,num_of_cols,A_mat) # Matrix A
    B = Matrix(num_of_rows,num_of_cols,B_mat) # Matrix B
    #------------------------------------------------------------------------------
    # STEP 2: MERGE THE TWO MATRICES INTO ONE MATRIX M=[-A|B]
    #------------------------------------------------------------------------------
    # Find the dimensions of the matrix A
    num_of_eq, n = A.shape        
    # FIRST ROW OF M
    # Begin by giving the matrix A the value of the first row 
    # in A
    M = -(A.row(0))
    # Then we add the zeroth row of B in the columns to the left
    for j in range(len(c_reduced)):
        M = M.col_insert(j+len(c_reduced), Matrix([B[0,j]]))
    # THE REMAINING ROWS OF M
    for i in range(num_of_eq):
        if i > 0:
            temp_row = -(A.row(i))
            # Then we add the zeroth row of B in the columns to the left
            for j in range(len(c_reduced)):
                temp_row = temp_row.col_insert(j+len(c_reduced), Matrix([B[i,j]]))
            # Insert the temporary row
            M = M.row_insert(i, temp_row)
    #------------------------------------------------------------------------------
    # STEP 3: REDUCE THE SYSTEM OF ODEs BY CALCULATING A BASIS FOR col(M^t)
    #------------------------------------------------------------------------------
    # Row-reduce the expanded matrix M=[-A|B]
    M= M.rref(simplify=True)[0]
    pivot_columns = M.rref(simplify=True)[1]
    #------------------------------------------------------------------------------
    # STEP 4: EXTRACT THE PURELY ALGEBRAIC EQUATIONS FROM THE SYSTEM OF ODEs
    #------------------------------------------------------------------------------    # Split the matrix into its component parts
    A = -M[:,0:(len(c_reduced))]
    B = M[:,len(c_reduced):(2*len(c_reduced))]
    # Calculate the dimensions of A
    num_of_rows, num_of_cols = A.shape
    # Extract the zero rows from $A$
    non_zero_rows = [i for i in range(num_of_rows) for j in range(num_of_cols) if A[i,j]!=0] # The non-zero rows
    zero_rows = [i for i in range(num_of_rows) if i not in non_zero_rows] # The zero rows
    # Define an algebraic matrix based on the zero rows in A
    B_algebraic = B[zero_rows,:]
    # Re-define A and B
    A = A[non_zero_rows,:]
    B = B[non_zero_rows,:]
    # Define the algebraic equations
    eq_alg = []
    # Define the number of algebraic equations
    num_of_alg, cols = B_algebraic.shape
    pivot_columns = B_algebraic.rref(simplify=True)[1]        
    # Add the equations hey?
    for rI in range(num_of_alg):
        # Define a temporary sum
        temp_sum = 0
        # Loop over the columns and perform
        # the matrix multiplication
        for cI in range(cols):
            temp_sum += B_algebraic[rI,cI]*c_reduced[cI](x[0])
        # Append the equation
        eq_alg.append(temp_sum)
    # Solve each equation
    print(latex(B_algebraic))
    print(latex(eq_alg))
    print(pivot_columns)
    print(latex(Matrix(len(c_reduced),1,c_reduced)))
    eq_alg = [solve(eq_alg[i],c_reduced[i](x[0])) for i in pivot_columns]
    #------------------------------------------------------------------------------
    # STEP 5: IF THE SYSTEM IS DIAGONISABLE IT WILL BE SOLVED WITH LINEAR ALGEBRA
    #------------------------------------------------------------------------------
    # If A is diagonizable we will try to
    # solve the system at hand
    if A.is_diagonalizable():
        #-------------------------------------------------------
        # MATRIX CALCULATIONS
        #-------------------------------------------------------
        # Allocate memory for our solutions
        solutions = zeros(len(c_reduced),1)
        # Diagonalise matrix inv(A)*B
        M = A.inv()*B
        # CASE 1: M is diagonizable
        if M.is_diagonalizable():
            (P, D) = M.diagonalize()
            # Define a new matrix exp_D
            exp_D = D
            # Calculate the rows in this matrix
            rows, cols = exp_D.shape
            # Loop over the diagonal, take the
            # exponent of that value times the time
            # and multiply with an arbitrary coefficient
            for rI in range(rows):
                # For some reason, we need to define the name
                # C1 in order to do a substitution
                exec("C%d = symbols('C%d')"%(coefficient_counter,coefficient_counter))
                # Change the matrix at hand 
                exec("exp_D[rI,rI] = C%d*exp(exp_D[rI,rI]*x[0])"%(coefficient_counter))
            # Calculate the solution
            num_of_eqs, cols = M.shape
            # Loop over columns and construct the solutions
            for c in range(cols):
                solutions += exp_D[c,c]*P[:,c]
            # Substitute these value in to each and every one of the equations
            for i in range(len(c_reduced)):
                # Loop over all coefficients in c and do the substitutions
                for index in range(len(c)):
                    # If the element contains our coefficient, replace it
                    # by the substitution
                    if c[index]==c_reduced[i]:
                        c[index] = simplify(solutions[i,0])
                    else:
                        # Neglect all coefficients that are integer, which most often mean that they take the value zero. In this case we just move on.
                        if type(c[index]) != int:
                            # Check if the current coefficient is a function
                            list_func = [a.func for a in c[index].atoms(function.AppliedUndef)]
                            # MOST IMPORTANTLY: Check if the current coefficient is multiplied by some strange factor. In this case, we shall only do the substitution directly. Otherwise, we should just move on.
                            list_mul = [a.func for a in c[index].atoms(Mul)]
                            # Also, we need to make sure that we can differentiate between symbols which the integration constants, e.g. C1, are and our lovely coefficients ,e.g. c_0_1, are not.
                            list_symb = [a.func for a in c[index].atoms(Symbol)] 
                            # Apparently there is a symbolic zero which can cause problems which is stored as a symbolic number. Apparently, there is a type called exactly number which the symbolic zero is stored as but not the coefficients. This will allow us to differentiate between the two. 
                            list_number = [a.func for a in c[index].atoms(Number)] 
                            # CASE 1: We have a function multiplied with various parameters, a so called "mul"=> Substitute the solution directly! 
                            if len(list_mul) != 0:
                                # Substitute the solution into the coefficient of interest
                                temp_expression = c[index]
                                c[index] = temp_expression.subs(c_reduced[i](x[0]),simplify(solutions[i,0]))          
            # The indices that we want to remove from the vector c
            indices_to_remove = [i for i in range(len(c_reduced)) for j in range(len(c)) if c[j]==c_reduced[i]]
            # Sort the indices in reversed order
            indices_to_remove.sort(reverse=True)
            # Loop over indices and remove them
            for i in indices_to_remove:
                del c_reduced[i]
            # Just the equation system to zero
            eq_sys = []
            eq_sys_new = eq_sys.copy()
            eq_alg = []
            return eq_alg, eq_sys_new
    #------------------------------------------------------------------------------
    # STEP 6: THE SYSTEM IS NOT DIAGONISABLE=> SOLVE EACH EQUATION INDIVIDUALLY
    #------------------------------------------------------------------------------
    # A is not diagonisable=> we solve each equation invidually
    # Make a copy of the original system
    eq_sys_modified = eq_sys.copy()
    # Substitute these value in to each and every one of the equations
    for i in pivot_columns:
        # Loop over all coefficients in c and do the substitutions
        for index in range(len(c)):
            # If the element contains our coefficient, replace it
            # by the substitution
            if c[index]==c_reduced[i]:
                c[index] = eq_alg[i][0]
            else:
                # Neglect all coefficients that are integer, which most often mean that they take the value zero. In this case we just move on.
                if type(c[index]) != int:
                    # Check if the current coefficient is a function
                    list_func = [a.func for a in c[index].atoms(function.AppliedUndef)]
                    # MOST IMPORTANTLY: Check if the current coefficient is multiplied by some strange factor. In this case, we shall only do the substitution directly. Otherwise, we should just move on.
                    list_mul = [a.func for a in c[index].atoms(Mul)]
                    # Also, we need to make sure that we can differentiate between symbols which the integration constants, e.g. C1, are and our lovely coefficients ,e.g. c_0_1, are not.
                    list_symb = [a.func for a in c[index].atoms(Symbol)] 
                    # Apparently there is a symbolic zero which can cause problems which is stored as a symbolic number. Apparently, there is a type called exactly number which the symbolic zero is stored as but not the coefficients. This will allow us to differentiate between the two. 
                    list_number = [a.func for a in c[index].atoms(Number)] 
                    # CASE 1: We have a function multiplied with various parameters, a so called "mul"=> Substitute the solution directly! 
                    if len(list_mul) != 0:
                        # Substitute the solution into the coefficient of interest
                        temp_expression = c[index]
                        c[index] = temp_expression.subs(c_reduced[i](x[0]),eq_alg[i][0])          
        # Loop over all equations and do the substitutions
        for eq_num in range(len(eq_sys)):
            eq_temp =  eq_sys_modified[eq_num].subs(c_reduced[i](x[0]),eq_alg[i][0])
            eq_sys_modified[eq_num] = simplify(eq_temp.doit())
    # Decrease the c_reduced vector:
    # List of indices to remove
    indices_to_remove = list(pivot_columns).copy()
    # Reverse the order of them
    indices_to_remove.sort(reverse=True)
    # Loop over equations and remove them
    for i in indices_to_remove:
        del c_reduced[i]
    # Let's see if we can solve the remaining system of equations
    # We have a variable saying that the system contains solvable equations
    solvable_equations_exists = True
    # A counter of the coefficients so that we can rename the automatically generated
    # Save all algebraic equations that are introduced
    eq_alg = []
    # Define the indices of the algebraic equations: 
    eq_alg_ind = [] # Needed in order to remove them from the equation system eq_sys
    # Update the system of equations
    eq_sys = eq_sys_modified
    # Iteration counter...
    iteration_counter = 0
    # Iteration_max
    iter_max = 3
    # We solve all solvable equations, until there exist no more of them!
    while solvable_equations_exists:
        # We do not want to be trapped here forever...
        if iteration_counter > iter_max:
            break
        # Start by checking that we have any equations left:
        if len(eq_sys) == 0: # All equations are solved?
            # In case the answer is yes: Time to get out!
            solvable_equations_exists = False
            # Go to the next iteration
            continue
        #------------------------------------------------------------------
        # REMOVE ALGEBRAIC EQUATIONS
        #------------------------------------------------------------------
        # Indices zero equations
        indices_zero_equations = []
        # Make sure that all derivatives are evaluated before we continue
        for index in range(len(eq_sys)):
            num, denom = fraction(simplify(eq_sys[index].doit()))
            eq_sys[index] = num
            # If we have a zero equation we will have to remove it...
            if eq_sys[index] == 0:
                indices_zero_equations.append(index)
            # Sort the indices in reverse order
            indices_zero_equations.sort(reverse=True)
            # Remove these equations
            for i in indices_zero_equations:
                del eq_sys[i]
        # Re-set indices
        indices_zero_equations = []
        # Before any equations are solved, we need to be sure
        # that no algebraic equations exist among our equations.
        # If so the program will crash when it tries to use
        # "dsolve" to solve the differential equations
        # Define a help_counter to keep track where the algebraic equation
        # is located in the list eq_sys
        help_counter = 0
        # Define a logical variable telling us if we have a differential equation or not
        differential_equation = False
        # Loop over the equations 
        for isolated_eq in eq_sys:
            # Find all coefficients in our current equation
            list_func = [a.func for a in isolated_eq.atoms(function.AppliedUndef)]
            # Loop over all coefficients
            for coeff in list_func:
                # Loop over all terms in the equation
                for term in isolated_eq.args:
                    # Check if we have the derivative, and 
                    # if so it is a differential equation 
                    # meaning that we abort
                    if term.has(Derivative(coeff(x[0]),x[0])):
                        differential_equation = True
                        break
                    # If we have a differential equation break the loop
                    if differential_equation:
                        break
                # Algebraic equation
                if not differential_equation:
                    eq_alg.append(isolated_eq)
                    eq_alg_ind.append(help_counter)
                else: # Differential equation
                    differential_equation = False
                help_counter += 1
            # Remove algebraic equations from the differential equations:
            # Sort the indices of the algebraic equations in descending order
            eq_alg_ind.sort(reverse=True)
            # Loop over these indices and remove them
            for index in eq_alg_ind:
                del eq_sys[index]
                # Very important detail, namely to reset the indices of the algebraic equations
            eq_alg_ind = []  
            #----------------------------------------------------------------------
            # SOLVE DIFFERENTIAL EQUATIONS
            #----------------------------------------------------------------------
            # Loop over all equations and save the solvable ones
            # Allocate memory for all solvable equations
            solutions_solvable_equations = []
            # Save a help counter
            help_counter = 0
            # Make a list where we save the indices where the solvable equations occurred
            indices_solvable_equations = []
            # Loop over all equations and save the ones with a single function
            for eq_temp in eq_sys:
                # Save all functions in a list
                list_of_functions = [a.func for a in eq_temp.atoms(function.AppliedUndef)]
                # See if we have merely one function
                if len(list_of_functions) == 1:
                    # Extract the equation at hand, and solve it
                    eq_save = dsolve(Eq(eq_temp,0))
                    # Extract the RHS of the solved equation
                    RHS = eq_save.rhs
                    # Introduce a temporary character needed
                    # for the substitution
                    temp_char = 'K'
                    # Use the home-built function "sub_const"
                    # which renames the newly introduced coefficients
                    # to the correct number given by the coefficient
                    # counter
                    RHS = sub_const(RHS, temp_char,coefficient_counter)           
                    # Lastly re-define our newly solved equation
                    eq_save = Eq(eq_save.lhs,RHS)
                    # Save the solution to our newly solved equation
                    # that is re-named appropriately as well
                    solutions_solvable_equations.append(eq_save)
                    # Increase the coefficient_counter
                    coefficient_counter += 1
                    # Save where in the equation system it occurs
                    indices_solvable_equations.append(help_counter)
                    # We break as soon as we find a solution
                    break                
                # Increment the help_counter
                help_counter += 1
            # If the number of solvable equations are zero, we are finished
            if len(solutions_solvable_equations) == 0:
                # Time to get out
                solvable_equations_exists = False
                continue
            else: # We have solvable equations
                #------------------------------------------------------------------
                # UPDATE COEFFICIENT VECTORS
                #------------------------------------------------------------------
                # Sort the indices of the solvable equations in descending
                # order (largest to smallest)
                indices_solvable_equations.sort(reverse=True)
                # Loop over these indices and remove them from the equations
                # as we have already solved them! 
                for index in indices_solvable_equations:
                    del eq_sys[index]
                # Make a vector where the solved coefficient occurs in the vector
                # c_reduced
                indices_solved_coefficients = []            
                # Loop over all solutions, and then we will substitute these solutions in three steps. Firstly, we will substitute them into the equation system. Secondly, we will substitute them in the coefficient vector c. Thirdly, wewill remove these coefficients in the reduced coefficient vector c_reduced. 
                for solution in solutions_solvable_equations:
                    # STEP 1: Loop over the solutions and substitute these solutions into all remaining equations Loop over the remaining equations 
                    for eq_number in range(len(eq_sys)):
                        # Save the current equation where the substitution is done
                        eq_temp = eq_sys[eq_number].subs(solution.lhs,solution.rhs)
                        num, denom = fraction(simplify(eq_temp.doit()))            
                        eq_temp = num
                        # Evaluate derivatives and simplify
                        eq_sys[eq_number] = eq_temp
                    # STEP 2: Loop over the coefficient vector c and substitute our newly find coefficients into the coefficient vector c.
                        # Loop over the coefficients
                    for c_index in range(len(c)):
                        # Neglect all coefficients that are integer, which most often mean that they take the value zero. In this case we just move on.
                        if type(c[c_index]) != int:
                            # If we have the solution, let's replace it!
                            if c[c_index] == solution.lhs:
                                # Replace the coefficient at hand
                                c[c_index] = solution.rhs
                            # Check if the current coefficient is a function
                            list_func = [a.func for a in c[c_index].atoms(function.AppliedUndef)]
                            # MOST IMPORTANTLY: Check if the current coefficient is multiplied by some strange factor. In this case, we shall only do the substitution directly. Otherwise, we should just move on.
                            list_mul = [a.func for a in c[c_index].atoms(Mul)]
                            # Also, we need to make sure that we can differentiate between symbols which the integration constants, e.g. C1, are and our lovely coefficients ,e.g. c_0_1, are not.
                            list_symb = [a.func for a in c[c_index].atoms(Symbol)] 
                            # Apparently there is a symbolic zero which can cause problems which is stored as a symbolic number. Apparently, there is a type called exactly number which the symbolic zero is stored as but not the coefficients. This will allow us to differentiate between the two. 
                            list_number = [a.func for a in c[c_index].atoms(Number)]  
                            # CASE 1: We have a function multiplied with various parameters, a so  called "mul"=> Substitute the solution directly! 
                            if len(list_mul) != 0:
                                # Substitute the solution into the coefficient of interest
                                temp_expression = c[c_index]
                                c[c_index] = temp_expression.subs(solution.lhs,solution.rhs)
                            else: # We do not have an "mul"-type expression. 
                                # Could be the desired coefficients, but we also have these sneaky integration coefficients, e.g. C1, which are symbols. We need to just ignore the integration coefficients. So we only do this for non-symbols, non-function and non-muls.
                                if len(list_func) == 0 and len(list_mul) == 0 and len(list_symb) == 0 and len(list_number) == 0: # We do not have an integration constant
                                    temp_expression = c[c_index]
                                    if temp_expression(x[0]) == solution.lhs:
                                        # Replace the coefficient at hand
                                        c[c_index] = solution.rhs
                    # STEP 3: Loop over the reduced coefficient vector c_reduced and find out where our coefficient of interest is stored. These indices will then be sorted in descending order so that we can remove these indices from our vector c_reduced.
                    for i in range(len(c_reduced)):
                        # We have our coefficient at the current index!
                        temp_expression = c_reduced[i]
                        if temp_expression(x[0]) == solution.lhs:
                            # Save our lovely index
                            indices_solved_coefficients.append(i)
                            # Only remove the unique indices so we do not remove a coefficient twice...
                        indices_solved_coefficients = list(set(indices_solved_coefficients))
                        # Now, we will sort the indices of the coefficient that are to be removed from c_reduced in descending order again. 
                    indices_solved_coefficients.sort(reverse=True)
                    # Loop over these indices and remove them
                    for index in indices_solved_coefficients:
                        del c_reduced[index]
                        # Re-set the solved indices
                    indices_solved_coefficients = []            
                    # Now we are back for another iteration in our while loop...
    #------------------------------------------------------------------------------
    # STEP 7: HANDLE THE ALGEBRAIC EQUATIONS
    #------------------------------------------------------------------------------
    # Introduce a help_counter
    help_counter = 0
    # Save the index for all zero equations
    null_equation_indices = []
    # Extract an algebraic equation
    for eq_curr in eq_alg:
        # We introduce a temporary equation
        eq_temp = eq_curr
        # Do we have a trivial equation? If yes, save the index
        if eq_curr == 0:
            null_equation_indices.append(help_counter)
        else:
            # See the number of functions we have
            #list_func = isolated_eq.atoms(function.AppliedUndef)
            list_func = [a for a in eq_curr.atoms(function.AppliedUndef)]
            # Find the position of the coordinate in list_func 
            # in the vector c_original
            index_coefficient = [i for i in range(len(c_original)) for j in range(len(list_func)) if c_original[i](x[0])==list(list_func)[j]]
            if len(list_func) != 0:
                # Substitute the expression for the coefficient and simplify the value
                num,denom = fraction(simplify(eq_temp.subs(list(list_func)[0],c[index_coefficient[0]])))
                eq_temp = num
            else:
                num,denom = fraction(simplify(eq_temp))
                eq_temp = num
                # Save the new algebraic equation after the substitution
                eq_alg[help_counter] = eq_temp
                # Update the help_counter
                help_counter+=1
    # Remove all zero equations from the list of equations
    # Sort the indices in descending order again. 
    null_equation_indices.sort(reverse=True)
    # Loop over these indices and remove them
    for index in null_equation_indices:
        del eq_alg[index]
    eq_sys_new = eq_sys.copy()
    return eq_alg, eq_sys_new
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# FUNCTION 8: "solve_system_with_algebraic_equations"
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# The function is a help function to Function 9 called "solve_determining_equations". It solves the remaining equations of the system of equation if one is in the case where one has algebraic equations in addition to the differential equations. In this case, the function converts the matrix system into one list containing a system of differential equations and another list containing the algebraic equations. Then the algebraic equations are substituted into the system of differential equations. Hopefully, one will end up with a situation where the system of differential equations contains single ODEs which can be integrated directly. Then the solution to this equation is substituted into the remaining system and hopefully this will trickle down so that it becomes a matter of solving single first order ODEs. 
def solve_system_with_algebraic_equations(x,eta_list,c_original,c,c_reduced,A,B,B_algebraic,variables):
    # Calculate where the zero rows of A are.
    # Calculate the dimension of the new A
    m,n = A.shape
    # If we already removed all rows we exit directly!
    if m == 0 or n == 0:
        eq_alg_original = []
        sol_alg = []
        return X,c,eq_alg_original,sol_alg    
    # Caclulate the non-zero rows in the new A
    rows_nonzero = [i for i in range(m) if any(A[i, j] != 0 for j in range(n))]
    cols_nonzero = [j for j in range(n) if any(A[i, j] != 0 for i in range(m))]
    # Caclulate the zero rows in the new A
    rows_zero = [i for i in range(m) if i not in rows_nonzero]
    # Insert new rows in the B_extract matrix
    num_rows_extract, cols = B_algebraic.shape
    counter = 0
    # Loop over the zero rows in A and transfer the corresponding 
    # rows in B to B_algebraic
    for index in rows_zero:
        B_algebraic = B_algebraic.row_insert(num_rows_extract+counter, B[index,:])
        counter += 1
    # Now, we have to tidy up B_extract again
    pivot_columns = B_algebraic.rref()[1]
    # And the row_reduced matrix. 
    B_algebraic = B_algebraic.rref()[0]    
    #------------------------------------------------------------------------------
    # STEP 1: HANDLE THE NON-ZERO PIVOT ELEMENTS
    #------------------------------------------------------------------------------
    # This is done by transferring the system from a matrix system to a system of equations. Then, each algebraic equation is solved and the corresponding solution is substituted into the system of equations which shrinks. Also, the coefficient vector c is updated with coefficient that are solved along the way. 
    #----------------------------------------------------
    # DIFFERENTIAL EQUATIONS
    #----------------------------------------------------
    # Define a list for the equations in the matrix A
    eq_sys = []
    # Loop over the rows to extract each equation
    for i in rows_nonzero:
        # Define the equation in a temporary variable
        temp_eq = 0
        # Loop over the columns and perform the matrix multiplication
        for j in range(n):
            # Extract the coefficient
            coeff = c_reduced[j]
            # Perform the matrix multiplication for the derivatives
            # of the coefficients (i.e the matrix A)
            temp_eq += A[i,j]*Derivative(coeff(x[0]),x[0])
            # Perform the matrix multiplication for the coefficients
            # themselves (i.e. the matrix)
            temp_eq += -( B[i,j] * coeff(x[0]) )
        # Save only the non-trivial equations
        if temp_eq != 0:
            # Save the equation in our vector
            eq_sys.append(temp_eq)
    #----------------------------------------------------
    # ALGEBRAIC EQUATIONS
    #----------------------------------------------------
    # Now, we do the same for the algebraic equations
    row_alg, n = B_algebraic.shape
    # Define a list for the algebraic equations
    eq_alg = []
    # Loop over the rows to extract each equation
    for i in range(row_alg):
        # Define the equation in a temporary variable
        temp_eq = 0
        # Loop over the columns and perform the matrix multiplication
        for j in range(n):
            # Extract the coefficient
            coeff = c_reduced[j]
            temp_eq += B_algebraic[i,j] * coeff(x[0]) 
        # Save only the non-trivial equations
        if temp_eq != 0:
            # Save the equation in our vector
            eq_alg.append(temp_eq)
    #----------------------------------------------------
    # SOLVE ALGEBRAIC EQUATIONS, UPDATE C AND REMOVE
    # COEFFICIENTS FROM THE SYSTEM OF DIFFERENTIAL 
    # EQUATIONS THAT WE HAVE LEFT CALLED eq_sys
    #----------------------------------------------------
    # Loop over all algebraic equations 
    for i in range(len(eq_alg)):
        # Extract the coefficient
        coeff = c_reduced[pivot_columns[i]]
        # Define the equation 
        eq = Eq(-eq_alg[i]+coeff(x[0]),coeff(x[0]))
        # Solve the equation
        sol = solve(eq,coeff(x[0]))
        # Loop through all non-algebraic equations and conduct
        # the substitution at hand
        for j in range(len(eq_sys)):
            # Extract an equation
            eq_temp = eq_sys[j]
            # Simplify and evaluate
            eq_temp = simplify(eq_temp.doit())
            # Do the appropriate substitution
            eq_sys[j] = eq_temp.subs(coeff(x[0]),sol[0])
        # Loop through the vector c (the original coefficient vector)
        # and do the same substitution
        for index in range(len(c)):
            # If the element contains our coefficient, replace it
            # by the substitution
            if c[index]==coeff:
                c[index] = sol[0]
    # Save all indices of equations that are potentially zero
    indices_null_equations = []
    # Add all indices of trivial or null equations
    for eq_index in range(len(eq_sys)):
        # Extract the equation
        eq_temp = eq_sys[eq_index]
        # The equation is zero?
        if eq_temp == 0:
            indices_null_equations.append(eq_index) # Save the index!
    # Now, we will sort the indices of the null equations in descending order
    # so that we can remove them from the system of equations
    indices_null_equations.sort(reverse=True)
    # Loop over these indices and remove them
    for index in indices_null_equations:
        del eq_sys[index]
    for index in pivot_columns[::-1]:
        # Remove these elements in c_reduced
        del c_reduced[index]
    #------------------------------------------------------------------------------
    # STEP 2: SOLVE EACH DIFFERENTIAL EQUATION THAT CONTAINS MERELY ONE FUNCTION
    #------------------------------------------------------------------------------
    # Due to the previous steps of the algorithm, the system is almost guarantueed to have ODEs with only on function that can be solved by caring out a single integration. If this solution is substituted into the remaining equations, more such equations will be created. This step is carried out until no more equations with just one unknown occurs. If some algebraic equations are created they will be saved in a list called "eq_alg" and they will be handled subsequently in the very last step.
    # We have a variable saying that the system contains solvable equations
    solvable_equations_exists = True
    # A counter of the coefficients so that we can rename the automatically generated coefficient which is always called C1 if just one integration is performed
    coefficient_counter = 1
    # Save all algebraic equations that are introduced
    eq_alg = []
    # Define the indices of the algebraic equations: 
    eq_alg_ind = [] # Needed in order to remove them from eq_sys
    # We solve all solvable equations, until there exist no more of them!
    while solvable_equations_exists:
        # Start by checking that we have any equations left:
        if len(eq_sys) == 0: # All equations are solved?
            # In case the answer is yes: Time to get out!
            solvable_equations_exists = False
            # Go to the next iteration
            continue
        #----------------------------------------------------------------------
        # REMOVE ALGEBRAIC EQUATIONS
        #----------------------------------------------------------------------
        # Indices for the zero equations
        indices_zero_equations = []
        # Make sure that all derivatives are evaluated before we continue
        for index in range(len(eq_sys)):
            num, denom = fraction(simplify(eq_sys[index].doit()))
            eq_sys[index] = num
            # If we have a zero equation we will have to remove it...
            if eq_sys[index] == 0:
                indices_zero_equations.append(index)
        # Sort the indices in reverse order
        indices_zero_equations.sort(reverse=True)
        # Remove these equations
        for i in indices_zero_equations:
            del eq_sys[i]
        # Re-set these equations
        indices_zero_equations = []            
        # Before any equations are solved, we need to be sure
        # that no algebraic equations exist among our equations.
        # If so the program will crash when it tries to use
        # "dsolve" to solve the differential equations
        # Define a help_counter to keep track where the algebraic equation
        # is located in the list eq_sys
        help_counter = 0
        # Define a logical variable telling us if we have a differential equation
        # or not
        differential_equation = False
        # Loop over the equations 
        for isolated_eq in eq_sys:
            # Find all coefficients in our current equation
            list_func = [a.func for a in isolated_eq.atoms(function.AppliedUndef)]
            # Loop over all coefficients
            for coeff in list_func:
                # Loop over all terms in the equation
                for term in isolated_eq.args:
                    # Check if we have the derivative, and 
                    # if so it is a differential equation 
                    # meaning that we abort
                    if term.has(Derivative(coeff(x[0]),x[0])):
                        differential_equation = True
                        break
                # If we have a differential equation break the loop
                if differential_equation:
                    break
            # Algebraic equation
            if not differential_equation:
                eq_alg.append(isolated_eq)
                eq_alg_ind.append(help_counter)
            else: # Differential equation
                differential_equation = False
            help_counter += 1
        # Remove algebraic equations from the differential equations:
        # Sort the indices of the algebraic equations in descending order
        eq_alg_ind.sort(reverse=True)
        # Loop over these indices and remove them
        for index in eq_alg_ind:
            del eq_sys[index]
        # Very important detail, namely to reset the indices of the algebraic equations
        eq_alg_ind = []  
        #----------------------------------------------------------------------
        # SOLVE DIFFERENTIAL EQUATIONS
        #----------------------------------------------------------------------
        # Loop over all equations and save the solvable ones
        # Allocate memory for all solvable equations
        solutions_solvable_equations = []
        # Save a help counter
        help_counter = 0
        # Make a list where we save the indices where the solvable equations occurred
        indices_solvable_equations = []
        # Loop over all equations and save the ones with a single function
        for eq_temp in eq_sys:
            # Save all functions in a list
            list_of_functions = [a.func for a in eq_temp.atoms(function.AppliedUndef)]
            # See if we have merely one function
            if len(list_of_functions) == 1:
                # Extract the equation at hand, and solve it
                eq_save = dsolve(Eq(eq_temp,0))
                # Extract the RHS of the solved equation
                RHS = eq_save.rhs
                # Introduce a temporary character needed
                # for the substitution
                temp_char = 'K'
                # Use the home-built function "sub_const"
                # which renames the newly introduced coefficients
                # to the correct number given by the coefficient
                # counter
                RHS = sub_const(RHS, temp_char,coefficient_counter)           
                # Lastly re-define our newly solved equation
                eq_save = Eq(eq_save.lhs,RHS)
                # Save the solution to our newly solved equation
                # that is re-named appropriately as well
                solutions_solvable_equations.append(eq_save)
                # Increase the coefficient_counter
                coefficient_counter += 1
                # Save where in the equation system it occurs
                indices_solvable_equations.append(help_counter)
                # We break as soon as we find a solution
                break                
            # Increment the help_counter
            help_counter += 1
        # If the number of solvable equations are zero, we are finished
        if len(solutions_solvable_equations) == 0:
            # Time to get out
            solvable_equations_exists = False
            continue
        else: # We have solvable equations
            #----------------------------------------------------------------------
            # UPDATE COEFFICIENT VECTORS
            #----------------------------------------------------------------------
            # Sort the indices of the solvable equations in descending
            # order (largest to smallest)
            indices_solvable_equations.sort(reverse=True)
            # Loop over these indices and remove them from the equations
            # as we have already solved them! 
            for index in indices_solvable_equations:
                del eq_sys[index]
            # Make a vector where the solved coefficient occurs in the vector
            # c_reduced
            indices_solved_coefficients = []            
            # Loop over all solutions, and then we will substitute these solutions in three steps. Firstly, we will substitute them into the equation system. Secondly, we will substitute them in the coefficient vector c. Thirdly, we will remove these coefficients in the reduced coefficient vector c_reduced. 
            for solution in solutions_solvable_equations:
                # STEP 1: Loop over the solutions and substitute these solutions
                # into all remaining equations
                # Loop over the remaining equations 
                for eq_number in range(len(eq_sys)):
                    # Save the current equation where the substitution is done
                    eq_temp = eq_sys[eq_number].subs(solution.lhs,solution.rhs)
                    num, denom = fraction(simplify(eq_temp.doit()))                
                    eq_temp = num
                    # Evaluate derivatives and simplify
                    eq_sys[eq_number] = eq_temp
                # STEP 2: Loop over the coefficient vector c and substitute our newly find coefficients into the coefficient vector c.Loop over the coefficients
                for c_index in range(len(c)):
                    # Neglect all coefficients that are integer, which most often mean that they take the value zero. In this case we just move on.
                    if type(c[c_index]) != int:
                        # If we have the solution, let's replace it!
                        if c[c_index] == solution.lhs:
                            # Replace the coefficient at hand
                            c[c_index] = solution.rhs
                        # Check if the current coefficient is a function
                        list_func = [a.func for a in c[c_index].atoms(function.AppliedUndef)]
                        # MOST IMPORTANTLY: Check if the current coefficient is multiplied by some strange factor. In this case, we shall only do the substitution directly. Otherwise, we should just move on.
                        list_mul = [a.func for a in c[c_index].atoms(Mul)]
                        # Also, we need to make sure that we can differentiate between symbols which the integration constants, e.g. C1, are and our lovely coefficients ,e.g. c_0_1, are not.
                        list_symb = [a.func for a in c[c_index].atoms(Symbol)] 
                        # Apparently there is a symbolic zero which can cause problems which is stored as a symbolic number. Apparently, there is a type called exactly number which the symbolic zero is stored as but not the coefficients. This will allow us to differentiate between the two. 
                        list_number = [a.func for a in c[c_index].atoms(Number)]  
                        # CASE 1: We have a function multiplied with various parameters, a so called "mul"=> Substitute the solution directly! 
                        if len(list_mul) != 0:
                            # Substitute the solution into the coefficient of interest
                            temp_expression = c[c_index]
                            c[c_index] = temp_expression.subs(solution.lhs,solution.rhs)
                        else: # We do not have an "mul"-type expression. 
                            # Could be the desired coefficients, but we also have these sneaky integration coefficients, e.g. C1, which are symbols. We need to just ignore the integration coefficients. So we only do this for non-symbols, non-function and non-muls.
                            if len(list_func) == 0 and len(list_mul) == 0 and len(list_symb) == 0 and len(list_number) == 0: # We do not have an integration constant
                                temp_expression = c[c_index]
                                if temp_expression(x[0]) == solution.lhs:
                                    # Replace the coefficient at hand
                                    c[c_index] = solution.rhs
                # STEP 3: Loop over the reduced coefficient vector c_reduced and find out where our coefficient of interest is stored. These indices will then be sorted in descending order so that we can remove these indices from our vector c_reduced.
                for i in range(len(c_reduced)):
                    # We have our coefficient at the current index!
                    temp_expression = c_reduced[i]
                    if temp_expression(x[0]) == solution.lhs:
                        # Save our lovely index
                        indices_solved_coefficients.append(i)
            # Only remove the unique indices so we do not remove a coefficient twice...
            indices_solved_coefficients = list(set(indices_solved_coefficients))
            # Now, we will sort the indices of the coefficient that are to be removed from c_reduced in descending order again. 
            indices_solved_coefficients.sort(reverse=True)
            # Loop over these indices and remove them
            for index in indices_solved_coefficients:
                del c_reduced[index]
            # Re-set the solved indices
            indices_solved_coefficients = []            
            # Now we are back for another iteration in our while loop...
    #------------------------------------------------------------------------------
    # STEP 3: SOLVE REMAINING DIFFERENTIAL EQUATIONS
    #------------------------------------------------------------------------------
    # If we still have differential equations we will have to solve them.
    # For this purpose we call our help function "solve_linear_sys_ODEs"
    # We call this "max_iter" times in the hope that we can solve the
    # remaining system
    max_iter = 5
    # Add a counter which helps us interrupt the while-loop
    iter_counter = 1
    # Loop until all equations are solved
    while len(eq_sys) != 0:
        # Try to solve the remaining system 
        eq_alg_extra, eq_sys_new = solve_linear_sys_ODEs(eq_sys,x,c_original,c,c_reduced,coefficient_counter,variables)
        eq_sys = eq_sys_new
        # If any extra algebraic equations were generated,
        # add them to our list of algebraic equations
        eq_alg += eq_alg_extra
        # Increase the iterations 
        iter_counter += 1
        # If we have done this more than the maximum number
        # of iterations, we simply give up
        if iter_counter>max_iter:
            break
    #If the function is successful it will find the remaining coefficients
    # in the ansätze for the tangents which currently reside in the vector c. 
    #------------------------------------------------------------------------------
    # STEP 4: CONDUCT SUBSTITUTIONS IN THE ALGEBRAIC EQUATIONS
    #------------------------------------------------------------------------------
    # Some of the algebraic equations contain most probably an expression for the constants in the tangent ansätze. So what we do is that we subsitute the value of the coefficients in the tangent ansätze into the algebraic equations. 
    # Introduce a help_counter
    help_counter = 0
    # Save the index for all zero equations
    null_equation_indices = []
    # Extract an algebraic equation
    for eq_curr in eq_alg:
        # We introduce a temporary equation
        eq_temp = eq_curr
        # Do we have a trivial equation? If yes, save the index
        if eq_curr == 0:
            null_equation_indices.append(help_counter)
        else:
            # See the number of functions we have
            #list_func = isolated_eq.atoms(function.AppliedUndef)
            list_func = [a for a in eq_curr.atoms(function.AppliedUndef)]
            # Find the position of the coordinate in list_func 
            # in the vector c_original
            index_coefficient = [i for i in range(len(c_original)) for j in range(len(list_func)) if c_original[i](x[0])==list(list_func)[j]]
            if len(list_func) != 0:
                # Substitute the expression for the coefficient and simplify the value
                num,denom = fraction(simplify(eq_temp.subs(list(list_func)[0],c[index_coefficient[0]])))
                eq_temp = num
            else:
                num,denom = fraction(simplify(eq_temp))
                eq_temp = num
            # Save the new algebraic equation after the substitution
            eq_alg[help_counter] = eq_temp
        # Update the help_counter
        help_counter+=1
    # Remove all zero equations from the list of equations
    # Sort the indices in descending order again. 
    null_equation_indices.sort(reverse=True)
    # Loop over these indices and remove them
    for index in null_equation_indices:
        del eq_alg[index]
    #------------------------------------------------------------------------------
    # STEP 5: SOLVE THE ALGEBRAIC EQUATIONS
    #------------------------------------------------------------------------------
    # As each individual integration was performed, there were most likely multiple arbitrary integration constants introduced called K1,...Kp where p is the number of integation that were carried out. We also have a system of algebraic equations containing these integration constants. So this part of the function solves the last algebraic equations. 
    # List of arbitrary integration constants (e.g. K1,K2,...)
    integration_const_list = []
    ## Coefficients in front of these
    integration_const_coeff_list = []
    # Save the RHS in each equation
    RHS_list = []
    # Loop over all algebraic equations
    for eq_curr in eq_alg:
        # Loop over all coefficients we have introduced
        for coeff_number in range(coefficient_counter):
            # Let's skip the very first coefficient
            if coeff_number > 0:
                # Loop over all terms in the equation
                for term in eq_curr.args:
                    # Allocate an arbitrary integration constant
                    exec("K_%d=symbols('K_%d')"%(int(coeff_number),int(coeff_number)))    
                    # If the integration constant exists we save it in a list
                    # of all arbitrary integration constants involved in
                    # the algebraic equations
                    exec("if eq_curr.coeff(K_%d)!=0:integration_const_list.append(K_%d)"%(int(coeff_number),int(coeff_number)))
    # Make sure that the integration constants only occur once:            
    integration_const_list = list(set(integration_const_list))                                    
    # Once again, loop over all algebraic equations
    for eq_curr in eq_alg:
        # Define a temporary sum
        temp_sum = 0
        # Define a help_counter
        help_counter = 0
        # Now given our constants we can also find their coefficients
        for integration_const in integration_const_list:
            # Extract the current coefficient
            coeff_curr = eq_curr.coeff(integration_const)
            # Append the coefficients
            integration_const_coeff_list.append(coeff_curr)
            # Increase the sum
            temp_sum += integration_const_list[help_counter] * coeff_curr
            # Increase the help_counter
            help_counter +=1
        # Now, we can calculate the RHS
        RHS_list.append(cancel(expand(temp_sum-eq_curr)))
    # Formulate the equations for the arbitrary integration coefficients
    # as a matrix equation
    integration_const_mat = Matrix(len(eq_alg),len(integration_const_list),integration_const_coeff_list)    
    coeff_mat = Matrix(len(integration_const_list),1,integration_const_list)
    RHS_mat = Matrix(len(eq_alg),1,RHS_list)
    # Define an expanded matrix M=[A|RHS]
    M = integration_const_mat
    # Add an additional column with the RHS
    M= M.col_insert(len(integration_const_list),RHS_mat[:,:])
    # Row-reduce the expanded matrix
    M= M.rref()[0]
    pivot_columns = M.rref()[1]
    # Go back to the original form
    integration_const_mat  = M[:,0:(len(integration_const_list))]
    RHS_mat = M[:,len(integration_const_list)]  
    # Define new algebraic equations
    eq_alg_original = eq_alg.copy()
    eq_alg = []
    # Define the new RHS 
    RHS_alg = []
    # Calculate the dimensions of our system
    rows, cols = integration_const_mat.shape
    # Loop over the rows
    for row in range(rows):
        # We use a temporary variable here to define
        # the equation
        eq_temp = 0
        # Loop over the columns
        for col in range(cols):
            # Perform the matrix multiplication
            eq_temp += integration_const_mat[row,col]*coeff_mat[col,0]
        # Add the equation to the list if it is non-trivial
        if eq_temp != 0:
            eq_alg.append(eq_temp)
            RHS_alg.append(RHS_mat[row,0])
    # See if we made it here...
    # Solve equations for our coefficients
    sol_alg = []
    help_iterator = 0
    for eq in eq_alg:
        # Append the solution i.e. the value of the coefficient
        sol_alg.append(solve(Eq(eq,RHS_alg[help_iterator]),coeff_mat[pivot_columns[help_iterator],0]))
        # Append the coefficient that takes on the solution
        sol_alg[help_iterator].append(coeff_mat[pivot_columns[help_iterator],0])
        # Check that we do not have any functions in the right hand side
        list_func = sol_alg[help_iterator][0].atoms(Function)
        # In this case, we take the easy way out and add the value 0
        # as a solution or rather as the value to our coefficient
        if len(list_func) != 0:
            sol_alg[help_iterator][0] = 0
        help_iterator += 1
    # See if we managed this...
    for c_i in range(len(sol_alg)):
        # Coefficient to be substituted
        coeff = sol_alg[c_i][1]
        # Solution that is the substitute
        solution = sol_alg[c_i][0]
        # Loop over the coefficients and substitute our solved coefficients
        for i in range(len(c)):
            list_symb = [a for a in c[i].atoms(Symbol)]
            if type(c[i]) != int and len(list_symb)!=0:
                # Substitute the value in each place of the coefficient vector
                c[i] = c[i].subs(coeff,solution)
    # Since we use dsolve it saves each solution as a list with
    # with the arguments lhs and rhs. So we only extract the solution
    sol_alg_temp = [sol_alg[i][0] for i in range(len(sol_alg))]
    #------------------------------------------------------------------------------
    # STEP 8: DEFINE OUR GENERATOR BY INSERTING THE SOLVED COEFFICIENTS INTO THE TANGENTS
    #-----------------------------------------------------------------------------
    # Define our generator
    X = "\\begin{equation}X="
    # Define a help_counter
    help_counter = 0
    # Loop over tangents
    for tangent_number in range(len(eta_list)):
        # Loop through all coefficients and replace them in all tangents
        for index in range(len(c)):
            # We can only substitute the values of c if the current element
            # in c is a symbol or a number (and not an unknown function)
            list_symb = [a.func for a in c[index].atoms(Symbol)] 
            list_number = [a.func for a in c[index].atoms(Number)]  
            # Do all substitutions 
            if len(list_symb)!= 0 or len(list_number) != 0:
                eta_list[tangent_number] = eta_list[tangent_number].subs(c_original[index](x[0]),c[index])
        if eta_list[tangent_number]!=0:
            X += "\\left(" + str(latex(eta_list[tangent_number])) + "\\right)" + "\partial x_{" + str(tangent_number) + "}"
            if tangent_number < (len(eta_list)-1):
                X+="+"                       
    # Add the final touch to the generator
    X += "\\end{equation}"
    X = X.replace("+\\end{equation}", ".\\end{equation}")
    # Return the solved coefficients and the generator
    return X,c,eq_alg_original,sol_alg    
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
# The script returns two outputs:
# 1. the coefficients in the tangent ansätze stored in a vector "c",
# 2. The calculated generator of the symmetry which is stored in a string called "X".
# It is common that a number of algebraic equations is returned as well. So the function attempts to solve these equations as well. In the end the function returns, the generator before the algebraic equations are solved, the generator after the algebraic equations are solved, the coefficient vectors before and after the the algebraic equations are solved, the algebraic equations and the solutions to the algebraic equation that the script finds.
def solve_determining_equations(x,eta_list,c,det_eq,variables):
    #------------------------------------------------------------------------------
    # STEP 1: FORMULATE THE DETERMINING EQUATIONS ON MATRIX FORM A*dc/dt=B*c
    #-----------------------------------------------------------------------------
    # Allocate memory for the two matrices. 
    # We allocate them as lists because we can
    # easily generate matrices from a list where
    # we provide the list of course but before 
    # that we just state how many rows and columns
    # we have. 
    A_mat = []
    B_mat = []
    # Loop over the the coefficient vector c and 
    # save the coefficents in front of the c's 
    # and the derivatives of c in our two matrices.
    no_occurrences_der = 0
    no_occurrences = 0
    help_iterator = 0    
    # Loop over all determining equations
    for isolated_eq in det_eq:
        # Loop over all coefficients in the tangential ansätze
        # and define our matrices
        for coeff in c:
            A_mat.append(-isolated_eq.coeff(Derivative(coeff(x[0]),x[0])))
            B_mat.append(isolated_eq.coeff(coeff(x[0])))
    # Finally we generate the matrices from our lists A_mat and B_mat
    num_of_rows = len(det_eq) # Number of rows
    num_of_cols = len(c) # Number of columns
    A = Matrix(num_of_rows,num_of_cols,A_mat) # Matrix A
    B = Matrix(num_of_rows,num_of_cols,B_mat) # Matrix B
    #------------------------------------------------------------------------------
    # Copy the coefficient vectors
    #-----------------------------------------------------------------------------
    # Copy the coefficient vector twice
    c_reduced = c.copy() # Will be reduced when we remove parameters
    c_original = c.copy() # Contains the original coordinates    
    #------------------------------------------------------------------------------
    # STEP 2: TRANSFER THE SYSTEM TO ONE MATRIX M=[A|-B] AND REDUCE THE NUMBER OF 
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
    for j in range(len(c_reduced)):
        M = M.col_insert(j+len(c_reduced), Matrix([B[0,j]]))
    # THE REMAINING ROWS OF M
    for i in range(num_of_eq):
        if i > 0:
            temp_row = -(A.row(i))
            # Then we add the zeroth row of B in the columns to the left
            for j in range(len(c_reduced)):
                temp_row = temp_row.col_insert(j+len(c_reduced), Matrix([B[i,j]]))
            # Insert the temporary row
            M = M.row_insert(i, temp_row)
    # Let us calculate a basis for "col(M^T)". 
    eq_sys = M.T # Take the transpose of the matrix
    reduced_sys = eq_sys.columnspace()
    M_tilde = reduced_sys[0].T
    for i in range(len(reduced_sys)):
        if i > 0:
            M_tilde = M_tilde.row_insert(i, reduced_sys[i].T)          
    # Split the matrix into its component parts
    A = -M_tilde[:,0:(len(c_reduced))]
    B = M_tilde[:,len(c_reduced):(2*len(c_reduced))]
    #------------------------------------------------------------------------------
    # STEP 3: FIND THE PURELY ALGEBRAIC EQUATIONS WHICH CORRESPONDS TO THE ZERO ROWS IN THE MATRIX A. THE SAME ROWS IN THE MATRIX B WILL BE FORMULATED INTO A NEW MATRIX B_algebraic WHICH CONTAINS ONLY THE ALGEBRAIC EQUATIONS
    #-----------------------------------------------------------------------------
    # Remove the zero rows from A:
    # Calculate the dimensions 
    m, n = A.shape
    # Caclulate the non-zero rows
    rows_nonzero = [i for i in range(m) if any(A[i, j] != 0 for j in range(n))]
    #cols_nonzero = [j for j in range(n) if any(A[i, j] != 0 for i in range(m))]
    # Caclulate the zero rows
    rows_zero = [i for i in range(m) if i not in rows_nonzero]
    # Update A
    A = A[rows_nonzero, :]
    # Extract zero bit from B which constitutes our algebraic equations
    B_algebraic = B[rows_zero,:]
    # Update B
    B = B[rows_nonzero,:]
    #------------------------------------------------------------------------------
    # STEP 4: FIND THE PIVOT ELEMENTS OF THE NEW MATRIX B_algebraic AND CLASSIFY THEM INTO TWO CATEGORIES: THE PIVOT ELEMENTS WHICH ARE ZERO AND THE NON-ZERO PIVOTELEMENTS. PUT ALL COEFFICIENTS OF THE ZERO PIVOT ELEMENTS TO ZERO AND REMOVE THE CORRESPONDING COLUMNS IN THE MATRICES A, B AND B_algebraic. 
    #-----------------------------------------------------------------------------
    # Row-reduce matrix M_tilde and calculate the pivot_columns
    pivot_columns = B_algebraic.rref()[1]
    # And the row_reduced matrix. 
    B_algebraic = B_algebraic.rref()[0]
    # Loop over the matrix and find out which pivot elements that correspond 
    # to a null row (except for the pivot element itself that is...)
    num_of_rows, num_of_cols = B_algebraic.shape    
    # Given the pivot columns, calculate the pivot rows
    pivot_rows = [i for i in range(num_of_rows) for j in pivot_columns if B_algebraic[i,j]==1 ]
    # CLASSIFY THE PIVOT ELEMENTS IN B_extract AS
    # EITHER ZERO PIVOT ELEMENTS OR NON-ZERO PIVOT
    # ELEMENTS:
    # Find the zero pivot columns
    pivot_columns_zero = []
    # Find the non-zero pivot columns
    pivot_columns_non_zero = []
    # Save pivot rows
    pivot_rows_zero = []
    # Find the zero-pivot element if you
    # know what I mean...
    # Use a help index
    help_index = 0
    # Loop over pivot rows
    for i in pivot_rows:
        # Allocate a row_sum
        row_sum = 0
        # Loop over all columns
        for j in range(num_of_cols):
            # Sum all elements except for the pivot element
            # itself
            if j != pivot_columns[help_index]:
                row_sum += B_algebraic[i,j]
        # If this sum is zero we have a so called zero pivot element
        if row_sum == 0:
            pivot_rows_zero.append(i) # Save this index
        # Increment the help index
        help_index += 1
    # Find the corresponding column for these rows.
    # In other words, find the Pivot columns that 
    # correspond to rows with only one non-zero element.
    for i in pivot_rows_zero:
        for j in pivot_columns:
            if B_algebraic[i,j]==1:
                pivot_columns_zero.append(j)
    # Find the other pivot elements
    pivot_columns_non_zero = list(set(pivot_columns).difference(pivot_columns_zero))
    # Firstly, set these elements to zero in c.
    # The beautiful programming trick here is to
    # do the looping in reverse order.
    for index in pivot_columns_zero[::-1]:
        # Set all these elements in c to zero.
        c[index] = 0
        # Remove these element in c_reduced
        del c_reduced[index]
        # Remove the columns
        A.col_del(index)
        B.col_del(index)
        B_algebraic.col_del(index)
    # Remove zero-rows in B_algebraic
    for i in pivot_rows_zero[::-1]:
        B_algebraic.row_del(i)
    #------------------------------------------------------------------------------
    # STEP 5: SOLVE THE REMAINING SYSTEM IF WE HAVE ALGEBRAIC EQUATIONS
    #-----------------------------------------------------------------------------
    # Solve the remaining part in the case where we have algebraic equations     
    X,c,eq_alg_original,sol_alg = solve_system_with_algebraic_equations(x,eta_list,c_original,c,c_reduced,A,B,B_algebraic,variables)
    #------------------------------------------------------------------------------
    # STEP 6: Return the output
    #-----------------------------------------------------------------------------
    # Return the solved coefficients and the generator
    return X,c_original,c,eq_alg_original,sol_alg

