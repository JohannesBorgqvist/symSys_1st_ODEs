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
    # STEP 1 of 7: FORMULATE THE DETERMINING EQUATIONS ON MATRIX FORM A*dc/dt=B*c
    #-----------------------------------------------------------------------------
    # Allocate memory for the two matrices. 
    # We allocate them as lists because we can
    # easily generate matrices from a list where
    # we provide the list of course but before 
    # that we just state how many rows and columns
    # we have. 
    A_mat = []
    B_mat = []
    # Loop over all ODEs and save the coefficients
    for isolated_eq in eq_sys:
        # Expand the equation at hand
        eq_temp = expand(isolated_eq)
        # Loop over all coefficients in the tangential ansätze
        # and define our matrices
        for coeff in c_reduced:
            # Save the derivative of the A coefficient
            A_mat.append(-eq_temp.coeff(Derivative(coeff(x[0]),x[0])))
            # Save the coefficient of the coefficient in the tangential
            # ansätze
            B_mat.append(eq_temp.coeff(coeff(x[0])))
    # Finally we generate the matrices from our lists A_mat and B_mat
    num_of_rows = len(eq_sys) # Number of rows
    num_of_cols = len(c_reduced) # Number of columns
    A = Matrix(num_of_rows,num_of_cols,A_mat) # Matrix A
    B = Matrix(num_of_rows,num_of_cols,B_mat) # Matrix B
    print("### The original matrix formulation.")
    print("Dimensions of A")
    print(A.shape)
    print("\\begin{equation*}")
    print("A=%s"%(str(latex(A))))
    print("\\end{equation*}")
    print("Dimensions of B")
    print(B.shape)
    print("\\begin{equation*}")
    print("B=%s"%(str(latex(B))))
    print("\\end{equation*}")
    print("The coefficient vector c:")
    print("\\begin{equation*}")
    print("\mathbf{c}=%s"%(str(latex(Matrix(len(c),1,c)))))
    print("\\end{equation*}")
    print("The coefficient vector c_reduced:")
    print("\\begin{equation*}")
    print("\mathbf{c}_{\mathrm{reduced}}=%s"%(str(latex(Matrix(len(c_reduced),1,c_reduced)))))
    print("\\end{equation*}")
    print("%d unknowns remain"%(int(len(c_reduced))))
    #------------------------------------------------------------------------------
    # STEP 2 of 7: TRANSFER THE SYSTEM TO ONE MATRIX M=[A|-B] AND REDUCE THE NUMBER OF 
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
    A = M_tilde[:,0:(len(c_reduced))]
    B = -M_tilde[:,len(c_reduced):(2*len(c_reduced))]
    print("### Reduction based on column space.")
    print("Dimensions of A")
    print(A.shape)
    print("\\begin{equation*}")
    print("A=%s"%(str(latex(A))))
    print("\\end{equation*}")
    print("Dimensions of B")
    print(B.shape)
    print("\\begin{equation*}")
    print("B=%s"%(str(latex(B))))
    print("\\end{equation*}")
    print("The coefficient vector c:")
    print("\\begin{equation*}")
    print("\mathbf{c}=%s"%(str(latex(Matrix(len(c),1,c)))))
    print("\\end{equation*}")
    print("The coefficient vector c_reduced:")
    print("\\begin{equation*}")
    print("\mathbf{c}_{\mathrm{reduced}}=%s"%(str(latex(Matrix(len(c_reduced),1,c_reduced)))))
    print("\\end{equation*}")
    print("%d unknowns remain"%(int(len(c_reduced))))    
    # Row-reduce the expanded matrix
    M_tilde = M_tilde.rref()[0]
    # Split the matrix into its component parts
    A = M_tilde[:,0:(len(c_reduced))]
    B = -M_tilde[:,len(c_reduced):(2*len(c_reduced))]
    print("### Row-reduction after column space.")
    print("Dimensions of A")
    print(A.shape)
    print("\\begin{equation*}")
    print("A=%s"%(str(latex(A))))
    print("\\end{equation*}")
    print("Dimensions of B")
    print(B.shape)
    print("\\begin{equation*}")
    print("B=%s"%(str(latex(B))))
    print("\\end{equation*}")
    print("The coefficient vector c:")
    print("\\begin{equation*}")
    print("\mathbf{c}=%s"%(str(latex(Matrix(len(c),1,c)))))
    print("\\end{equation*}")
    print("The coefficient vector c_reduced:")
    print("\\begin{equation*}")
    print("\mathbf{c}_{\mathrm{reduced}}=%s"%(str(latex(Matrix(len(c_reduced),1,c_reduced)))))
    print("\\end{equation*}")
    print("%d unknowns remain"%(int(len(c_reduced))))    
    #------------------------------------------------------------------------------
    # STEP 3 of 7: FIND THE PURELY ALGEBRAIC EQUATIONS WHICH CORRESPONDS TO THE ZERO ROWS IN THE MATRIX A. THE SAME ROWS IN THE MATRIX B WILL BE FORMULATED INTO A NEW MATRIX B_algebraic WHICH CONTAINS ONLY THE ALGEBRAIC EQUATIONS
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
    print("### Defining algebraic equations.")
    print("Dimensions of A")
    print(A.shape)
    print("\\begin{equation*}")
    print("A=%s"%(str(latex(A))))
    print("\\end{equation*}")
    print("Dimensions of B")
    print(B.shape)
    print("\\begin{equation*}")
    print("B=%s"%(str(latex(B))))
    print("\\end{equation*}")
    print("Dimensions of B_algebraic")
    print(B_algebraic.shape)
    print("\\begin{equation*}")
    print("B_{\mathrm{algebraic}}=%s"%(str(latex(B_algebraic))))
    print("\\end{equation*}")    
    print("The coefficient vector c:")
    print("\\begin{equation*}")
    print("\mathbf{c}=%s"%(str(latex(Matrix(len(c),1,c)))))
    print("\\end{equation*}")
    print("The coefficient vector c_reduced:")
    print("\\begin{equation*}")
    print("\mathbf{c}_{\mathrm{reduced}}=%s"%(str(latex(Matrix(len(c_reduced),1,c_reduced)))))
    print("\\end{equation*}")
    print("%d unknowns remain"%(int(len(c_reduced))))
    # PRINT THE EIGEN VECTORS OF B
    print("Eigen values of B")
    print(latex(B.eigenvals()))  #returns eigenvalues and their algebraic multiplicity
    print("Eigen vectors of B")
    eigen_vectors_B = B.eigenvects()[0]
    eigen_vectors_B = eigen_vectors_B[2]
    print(latex(eigen_vectors_B))  #returns eigenvalues, eigenvects
    print("Number of eigen vectors")
    print(len(eigen_vectors_B))
    
    print("Is B diagonisible?")
    print(B.is_diagonalizable())
    print("Calculate rank of the matrix B:")
    print(B.rank())
    # Try to solve the
    mat_temp = B
    num_rows, num_cols = B.shape
    v_temp = Matrix(len(list(eigen_vectors_B[0])),1,list(eigen_vectors_B[0]))
    print("Eigen vector of interest:")
    print(latex(v_temp,mode='equation'))
    mat_temp = mat_temp.col_insert(num_cols,v_temp)
    print("Original matrix B:")
    print(latex(B,mode='equation'))
    print("Expanded matrix [B|v_1]:")
    print(latex(mat_temp,mode='equation'))
    print("Row-reduced expanded matrix")
    print(latex(mat_temp.rref()[0],mode='equation'))
    print("A solution using the LU solver")
    soln = B.LUsolve(v_temp)
    print(latex(soln,mode='equation'))
    

    temp_list = [1, -1, 1, 3]
    temp_mat = Matrix(2,2,temp_list)
    print("Example matrix")
    print(latex(temp_mat,mode='equation'))
    print("Eigen values of example matrix")
    print(latex(temp_mat.eigenvals()))
    print("Eigen vectors of example Matrix")
    eigen_vectors_temp = temp_mat.eigenvects()[0]
    eigen_vectors_temp = eigen_vectors_temp[2]
    print(latex(eigen_vectors_temp))  #returns eigenvalues, eigenvects
    print("Number of eigen vectors")
    print(len(eigen_vectors_temp))    
    #------------------------------------------------------------------------------
    # STEP 7 of 7: RETURN THE SYSTEM OF EQUATIONS AND THE ALGEBRAIC EQUATIONS
    #------------------------------------------------------------------------------
    # For some reason the system of ODEs does not
    # seem to update when we exit the function.
    # So we make a copy and return that copy which
    # works    
    #eq_sys_new = eq_sys.copy()
    eq_sys_new = []
    # The same goes for the coefficient_counter which
    # is not updated. So we return a copy of it
    coefficient_counter_new = coefficient_counter
    # Algebraic equations
    eq_alg = []
    # Return the algebraic equations and the ODEs
    return eq_alg, eq_sys_new,coefficient_counter_new    
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# FUNCTION 8: "solve_algebraic_equations"
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
def solve_algebraic_equations(x,c,c_original,eq_alg,coefficient_counter,eta_list):
    # Before anything is done, let us make a copy of the original algebraic
    # equations
    eq_alg_original = eq_alg.copy()
    # The reasoning behind this is that if we fail to solve these algebraic
    # equations it is of interest to save the equations so that they can be
    # analysed afterwards.
    #------------------------------------------------------------------------------
    # STEP 1 OF 3: CONDUCT SUBSTITUTIONS IN THE ALGEBRAIC EQUATIONS
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
                print("list_func")
                print(list_func)
                print("c")
                print(c)
                print("Current coefficient")
                print(c[index_coefficient[0]])
                print("Index coefficient")
                print(index_coefficient)
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
    # STEP 2 OF 3: SOLVE THE ALGEBRAIC EQUATIONS
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
    # STEP 3 OF 3: DEFINE OUR GENERATOR BY INSERTING THE SOLVED COEFFICIENTS INTO THE TANGENTS
    #-----------------------------------------------------------------------------
    # Define our generator
    X = "\\begin{align*}\nX&="
    # Define a help_counter
    help_counter = 0
    if help_counter == 1:
        # Loop over tangents
        for tangent_number in range(len(eta_list)):
            # Loop through all coefficients and replace them in all tangents
            for index in range(len(c)):            
                # Substitute coefficient in the tangents
                eta_list[tangent_number] = eta_list[tangent_number].subs(c_original[index](x[0]),c[index])
                # Add all non-trivial tangents to the generator    
                if eta_list[tangent_number]!=0:
                    X += "\\left(" + str(latex(eta_list[tangent_number])) + "\\right)" + "\partial x_{" + str(tangent_number) + "}\\\\"
                if tangent_number < (len(eta_list)-1):
                    X+="&+"                       
        # Add the final touch to the generator
        X += "\\end{align*}"
        X = X.replace("+\\end{align*}", ".\\end{align*}")
    else:
        X += "\\\\\n\\end{align*}"
        sol_alg = []
    # Return the solved coefficients and the generator
    return X,eq_alg_original,sol_alg    
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
    # To simplify it we have only two steps:
    # 1. Solve the ODEs,
    # 2. Solve the algebraic equations.
    # Both these are done in two functions called "solve_linear_sys_ODEs"
    # and "solve_algebraic_equations". To set this up we need to define
    # a few variables first.
    eq_sys = det_eq # System of ODEs
    eq_alg = []# System of algebraic equations
    # Coefficient vectors: we have three of them.
    # 1. c which is the coefficient vector that is updated
    # as each coefficient is solved.
    # 2. c_reduced which is reduced throughout the calculations
    # where each solved coefficient is removed.
    # 3. c_original which remains the same throughout the iterations
    # and is there in order to access the original coefficient which
    # us needed for certain substitutions.
    c_reduced = c.copy() # The reduced vector
    c_original = c.copy() # The original parameters
    # Throughout the calculations a number of integrations are performed
    # which entails the introduction of arbitrary integration coefficients.
    # In order to keep track of these we need a "coefficient counter" which
    # tells us what index the newest integration constant should have.
    coefficient_counter = 1
    #----------------------------------------------------------------
    #----------------------------------------------------------------
    # Solve ODEs until no ODE remain
    #----------------------------------------------------------------
    #----------------------------------------------------------------    
    # We call this "max_iter" times in the hope that we can solve the
    # remaining system
    max_iter = 1
    # Add a counter which helps us interrupt the while-loop
    iter_counter = 1
    # Loop until all equations are solved
    while len(eq_sys) != 0:
        # Try to solve the remaining system with the solution algorithm
        print("---------------------------------------")
        print("## ITER %d OUT OF %d" %(int(iter_counter),int(max_iter)))
        print("---------------------------------------")        
        eq_alg_extra, eq_sys_new, coefficient_counter_new = solve_linear_sys_ODEs(eq_sys,x,c_original,c,c_reduced,coefficient_counter,variables)
        # Update the coefficient counter
        coefficient_counter += coefficient_counter_new
        # Update the system of ODEs
        eq_sys = eq_sys_new
        # Update the system of algebraic equations
        eq_alg += eq_alg_extra
        # Increase the iterations 
        iter_counter += 1
        # If we have done this more than the maximum number
        # of iterations, we simply give up
        if iter_counter>max_iter:
            break
    #----------------------------------------------------------------
    #----------------------------------------------------------------
    # Solve algebraic equations
    #----------------------------------------------------------------
    #----------------------------------------------------------------
    print("---------------------------------------")
    print("## SOLVE ALGEBRAIC EQUATIONS")
    print("---------------------------------------")        
    X,eq_alg_original,sol_alg = solve_algebraic_equations(x,c,c_original,eq_alg,coefficient_counter,eta_list)    
    # Return the solved coefficients and the generator
    return X,c_original,c,eq_alg_original,sol_alg

