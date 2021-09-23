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
#from symengine import *
#import sympy
from sympy import *
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
# FUNCTION 6: "solve_determining_equations"
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
def solve_determining_equations(x,eta_list,c,det_eq,variables,omega_list):
    print("# Solving the determining equations")
    #------------------------------------------------------------------------------
    # STEP 1 of 6: FORMULATE THE DETERMINING EQUATIONS ON MATRIX FORM A*dc/dt=B*c
    #-----------------------------------------------------------------------------
    print("## Step 0 of 6: Defining the matrices from the determining equations")
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
        print("### Determining equation %d out of %d"%(help_counter,len(det_eq)))
        help_counter+=1
        eq_temp_copy = eq_temp
        for index in range(len(x)):
            eq_temp_copy = eq_temp_copy.subs(x[index],variables[index])
        print("\\begin{equation}\n%s=0\n\\end{equation}"%(latex(eq_temp_copy)))
        
        # Loop over all coefficients in the tangential ansätze
        # and define our matrices
        for koeff in c:
            # Save the derivative of the A coefficient
            A_mat.append(-eq_temp.coeff(Derivative(koeff(x[0]),x[0])))
            A_row.append(-eq_temp.coeff(Derivative(koeff(x[0]),x[0])))
            # Save the coefficient of the coefficient in the tangential
            # ansätze
            B_mat.append(eq_temp.coeff(koeff(x[0])))
            B_row.append(eq_temp.coeff(koeff(x[0])))
        A_row = Matrix(1,len(c),A_row)
        B_row = Matrix(1,len(c),B_row)            
        for index in range(len(x)):
            A_row = A_row.subs(x[index],variables[index])
            B_row = B_row.subs(x[index],variables[index])
        print("\\begin{equation}\n%s\\dfrac{\\mathrm{d}}{\\mathrm{d}t}%s=%s%s\n\\end{equation}"%(latex(A_row),latex(Matrix(len(c),1,c)),latex(B_row),latex(Matrix(len(c),1,c))))
        A_row = []
        B_row = []
    # Finally we generate the matrices from our lists A_mat and B_mat
    num_of_rows = len(det_eq) # Number of rows
    num_of_cols = len(c) # Number of columns
    A = Matrix(num_of_rows,num_of_cols,A_mat) # Matrix A
    B = Matrix(num_of_rows,num_of_cols,B_mat) # Matrix B
    print("## Step 1 of 6: the initial matrices")
    print("Dimension of matrices:\t%dX%d<br>"%(int(num_of_rows),int(num_of_cols)))
    print("Matrix A<br>")
    A_temp = A
    print("Matrix B<br>")
    B_temp = B
    for index in range(len(x)):
        A_temp = A_temp.subs(x[index],variables[index])
        B_temp = B_temp.subs(x[index],variables[index])
    print("\\begin{equation}\nA=%s\n\\end{equation}"%(latex(A_temp)))
    print("\\begin{equation}\nB=%s\n\\end{equation}"%(latex(B_temp)))    
    
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
    print("## Step 2 of 6: the reduced based on col(M^T) where M=[-A|B]<br>")
    r_temp, c_temp = A.shape
    print("Dimension of matrices:\t%dX%d<br>"%(int(r_temp),int(c_temp)))
    A_temp = A
    B_temp = B
    for index in range(len(x)):
        A_temp = A_temp.subs(x[index],variables[index])
        B_temp = B_temp.subs(x[index],variables[index])
    print("Matrix A<br>")        
    print("\\begin{equation}\nA=%s\n\\end{equation}"%(latex(A_temp)))
    print("Matrix B<br>")
    print("\\begin{equation}\nB=%s\n\\end{equation}"%(latex(B_temp)))     
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
    print("## Step 3 of 6: Splitting up to A, B and B_algebraic")
    r_temp, c_temp = A.shape
    r_temp_2, c_temp_2 = B_algebraic.shape    
    print("Dimension of matrices A and B:\t%dX%d<br>"%(int(r_temp),int(c_temp)))
    print("Dimension of matrices B_algebraic:\t%dX%d<br>"%(int(r_temp_2),int(c_temp_2)))    
    A_temp = A
    B_temp = B
    B_alg_temp = B_algebraic
    for index in range(len(x)):
        A_temp = A_temp.subs(x[index],variables[index])
        B_temp = B_temp.subs(x[index],variables[index])
        B_alg_temp = B_alg_temp.subs(x[index],variables[index])        
    print("Matrix A<br>")        
    print("\\begin{equation}\nA=%s\n\\end{equation}"%(latex(A_temp)))
    print("Matrix B<br>")
    print("\\begin{equation}\nB=%s\n\\end{equation}"%(latex(B_temp)))
    print("Matrix B_algebraic<br>")
    print("\\begin{equation}\nB_{\\textrm{algebraic}}=%s\n\\end{equation}"%(latex(B_alg_temp)))
    print("Coefficient matrix c:")
    print("\\begin{equation}\n\\mathbf{c}=%s\n\\end{equation}"%(latex(Matrix(len(c),1,c))))
    #------------------------------------------------------------------------------
    # STEP 4 of 6: REMOVE (POTENTIAL) EXTRA COLUMNS 
    #-----------------------------------------------------------------------------
    # Calculate rows and columns of matrices
    r_algebraic, cols = B_algebraic.shape # Dimensions of B_algebraic
    rows, cols = A.shape # Dimensions of A
    print("Dimensions of A:\t%dX%d\n"%(rows,cols))
    print("Dimensions of B_algebraic:\t%dX%d\n"%(r_algebraic,cols))
    # Create a set of all the columns
    col_set = set(range(cols))
    # Find the pivot columns of A
    pivot_columns = A.rref()[1]
    # Calculate the non-pivot columns
    non_pivot_columns = list(col_set.difference(set(pivot_columns)))
    # Substitute the value zero for all non-pivot columns in the coefficient
    # vector
    # Loop over tangents
    for tangent_number in range(len(eta_list)):
        # Loop through all coefficients and replace them in all tangents
        for index in non_pivot_columns:            
            # Substitute coefficient in the tangents
            eta_list[tangent_number] = eta_list[tangent_number].subs(c[index](x[0]),0)
    # We assume that we have a homogeneous system
    non_homogeneous = False
    # Check if we have non-pivot columns, that is we have
    # an underdetermined system with extra columns (i.e. more columns than equations).
    # In this case, we move these column to the right hand side and we then treat
    # these as source terms if you will, so that we get a quadratic system with some
    # extra inhomogeneities. 
    if len(non_pivot_columns)!=0: # Yes, we do have non-pivot coefficients
        # Print the non-pivot columns
        print("Non-Pivot columns")
        print(non_pivot_columns)
        # Indicate that we have a non_homogeneous system
        non_homogeneous = True
        # Define a list of all non-pivot elements.
        # We will refer to this as a non-homogeneous
        # source term
        source_ODE = []
        source_algebraic = []        
        # DEFINE SOURCE TERM FOR ODEs
        for r in range(rows):
            # Define a temporary source term
            temp_source = 0
            # Loop over the non-pivot columns and perform
            # the matrix multiplication in order to define
            # the source term
            for pC in non_pivot_columns:
                temp_source += -A[r,pC]*Derivative(c[pC](x[0]),x[0])
                temp_source += B[r,pC]*c[pC](x[0])
            # Append the temporary source
            source_ODE.append(temp_source)
        # DEFINE SOURCE TERM FOR ALGEBRAIC EQUATIONS
        for r in range(r_algebraic):
            # Define a temporary source term
            temp_source_alg = 0
            # Loop over the non-pivot columns and perform
            # the matrix multiplication
            for pC in non_pivot_columns:
                temp_source_alg += -B_algebraic[r,pC]*c[pC](x[0])
            # Append the non-homogeneous source 
            source_algebraic.append(temp_source_alg)
        # Print the matrix
        print("The non-homogeneous source term:<br>")
        print(latex(Matrix(len(source_ODE),1,source_ODE),mode='equation*'))
        source_ODE = Matrix(len(source_ODE),1,source_ODE)
        print("Dimensions:<br>")
        print(source_ODE.shape)
        print("The algebraic source term:<br>")
        source_alg = Matrix(len(source_algebraic),1,source_algebraic)
        print(latex(source_alg,mode='equation*'))
        print("Dimensions:<br>")
        print(source_alg.shape)
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
    print("## Step 4 of 6: Removing potential extra pivot columns")
    r_temp, c_temp = A.shape
    r_temp_2, c_temp_2 = B_algebraic.shape    
    print("Dimension of matrices A and B:\t%dX%d<br>"%(int(r_temp),int(c_temp)))
    print("Dimension of matrices B_algebraic:\t%dX%d<br>"%(int(r_temp_2),int(c_temp_2)))    
    A_temp = A
    B_temp = B
    B_alg_temp = B_algebraic
    for index in range(len(x)):
        A_temp = A_temp.subs(x[index],variables[index])
        B_temp = B_temp.subs(x[index],variables[index])
        B_alg_temp = B_alg_temp.subs(x[index],variables[index])        
    print("Matrix A<br>")        
    print("\\begin{equation}\nA=%s\n\\end{equation}"%(latex(A_temp)))
    print("Matrix B<br>")
    print("\\begin{equation}\nB=%s\n\\end{equation}"%(latex(B_temp)))
    #------------------------------------------------------------------------------
    # STEP 5 of 6: SOLVE THE ODE SYSTEM PROVIDED THAT IT IS QUADRATIC AND
    # SUBSTITUTE THE ALGEBRAIC EQUATIONS
    #-----------------------------------------------------------------------------
    num_rows, num_cols = B.shape
    print("Dimensions of B:\t%dX%d<br>"%(num_rows,num_cols))
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
            print("Coefficients:<br>")
            print(latex(c_mat,mode='equation'))
            print("Number of unknowns:\t%d<br>"%(int(len(c_mat))))
            # PRINT TO THE USER
            print("## Step 5 of 6: Solving the ODE system")
            r_temp, c_temp = B.shape
            r_temp_2, c_temp_2 = B_algebraic.shape    
            print("Dimension of the matrix B:\t%dX%d<br>"%(int(r_temp),int(c_temp)))
            print("Dimension of the matrix B_algebraic:\t%dX%d<br>"%(int(r_temp_2),int(c_temp_2)))    
            B_temp = B
            B_alg_temp = B_algebraic
            for index in range(len(x)):
                A_temp = A_temp.subs(x[index],variables[index])
                B_temp = B_temp.subs(x[index],variables[index])
                B_alg_temp = B_alg_temp.subs(x[index],variables[index])        
            print("ODE system:<br>")        
            print("\\begin{equation}\n\dfrac{\mathrm{d}}{\mathrm{d}t}%s=%s%s\n\\end{equation}"%(latex(Matrix(len(c),1,c)),latex(B_temp),latex(Matrix(len(c),1,c))))     
            #----------------------------------------------
            # PART 2: Solve ODE system
            #----------------------------------------------
            print("Solve the ODE system:<br>")
            print("Initial conditions for $\\mathbf{c}$ denoted by $\mathbf{c}_0$ in terms of arbitrary integration constants:")
            print("\\begin{equation}\n\mathbf{c}_{0}=%s=%s\n\\end{equation}"%(latex(Matrix(len(c),1,c)),latex(c_mat)))
            # Calculate the Jordan form of the matrix B
            P,J = B.jordan_form()
            #P,J = (-B).jordan_form()
            print("Jordan form:<br>")
            J_temp = J
            P_temp = P
            for index in range(len(x)):
                J_temp = J_temp.subs(x[index],variables[index])
                P_temp = P_temp.subs(x[index],variables[index])
            print("\\begin{equation}\nP=%s\n\\end{equation}"%(latex(P_temp)))
            print("\\begin{equation}\nJ=%s\n\\end{equation}"%(latex(J_temp)))                            
            # Re-define our matrix J by scaling it by the
            # independent variable x[0]
            J = x[0]*J
            # Re-define J by taking the matrix exponential
            J = J.exp()
            J_temp = J
            for index in range(len(x)):
                J_temp = J_temp.subs(x[index],variables[index])
            print("Exponential form:<br>")
            print("\\begin{equation}\n\\exp\\left(J\\cdot t\\right)=%s\n\\end{equation}"%(latex(J_temp)))
            # Calculate the homogeneous solution
            homo_sol = (P*J*P.inv())*c_mat
            homo_sol_temp = homo_sol
            for index in range(len(x)):
                homo_sol_temp = homo_sol_temp.subs(x[index],variables[index])            
            print("Solution to the ODE system:")
            print("\\begin{equation}\nP\\exp\\left(J\\cdot t\\right)P^{-1}\mathbf{c}_{0}=%s\n\\end{equation}"%(latex(homo_sol_temp)))                 
            # Add the particular solution
            if non_homogeneous: # non-homogeneous ODE
                # Begin by defining the integrand of the particular
                # solution as the exponential matrix times the source
                # term of the ODE
                part_sol = (P*J.inv()*P.inv())*source_ODE
                # Define the dimensions of the particular solution at
                # hand
                m,n = part_sol.shape
                # Loop over the rows and columns and integrate each
                # element of the vector with the particular solution
                for rI in range(m):
                    for cI in range(n):
                        # Define each matrix element as an integral 
                        # with respect to the zeroth variable 
                        part_sol[rI,cI] = Integral(part_sol[rI,cI],x[0])
                        # If possible, try to evaluate the integral at hand
                        part_sol[rI,cI] = part_sol[rI,cI].doit()
                # Lastly, multiply with J again to construct the particular solution
                part_sol = (P*J*P.inv())*part_sol
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
            print("## Step 6 of 6: Solving the algebraic system<br>")
            m, n = B_algebraic.shape
            print("Number of algebraic equations:\t%d\n"%(m))
            print("Matrix B_algebraic<br>")
            print("\\begin{equation}\nB_{\\textrm{algebraic}}=%s\n\\end{equation}"%(latex(B_alg_temp)))
            print("Algebraic equations:<br>")
            print("\\begin{equation}\n%s%s=%s\n\\end{equation}"%(latex(B_alg_temp),latex(Matrix(len(c),1,c)),latex(zeros(len(c),1))))
            c_alg = B_algebraic*homo_sol
            c_alg_temp = c_alg
            for index in range(len(x)):
                c_alg_temp = c_alg_temp.subs(x[index],variables[index])                        
            print("Algebraic equations after substitution of the solution to the ODE system:<br>")
            print("\\begin{equation}\n%s%s=%s=%s\n\\end{equation}"%(str(latex(B_alg_temp)),str(latex(homo_sol_temp)),str(latex(c_alg_temp)),latex(zeros(m,1))))




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
            # Define a vector based on this
            c_r = Matrix(len(const_remaining),1,const_remaining)
            # Now, we rewrite c_alg as a matrix muliplication
            # according to c_alg = alg_mat*c_r
            alg_mat_list = []
            # Loop over the algebraic equations 
            for eq_index in range(alg_row):
                for const in const_remaining:
                    # Extract temporary equation
                    eq_temp = c_alg[eq_index,0]
                    # Append the coefficient
                    alg_mat_list.append(eq_temp.coeff(const))
            # Now we define a matrix, hey?
            alg_mat = Matrix(alg_row,len(const_remaining),alg_mat_list)
            # Consistency check to see if the two matrices are equal
            #if c_alg == alg_mat*c_r:
                #print("The algebraic equations:<br>")
                #print("\\begin{equation}\n%s=%s%s\n\\end{equation}"%(str(latex(c_alg)),str(latex(alg_mat)),str(latex(c_r))))

            # Define a right hand side
            if non_homogeneous:
                RHS = -B_algebraic*part_sol + source_alg # The non-homogeneities
            else:
                RHS = zeros(r_algebraic,1) # Zeros as the right hand side
            # Now, we define an extended matrix, row reduce it and split
            # it back up to its component parts!
            # Create the extendent matrix by adding the RHS
            # to alg_mat
            extended_mat = alg_mat
            rows_alg, cols_alg = alg_mat.shape
            extended_mat = extended_mat.col_insert(cols_alg,RHS)
            # Row reduce this matrix
            extended_mat = extended_mat.rref()[0]
            pivot_elements_2 = extended_mat.rref()[1]
            # Split it back up into its component parts
            alg_mat = extended_mat[:,0:cols_alg]
            RHS = extended_mat[:,-1]
            # Update algebraic equations
            c_alg = alg_mat*c_r - RHS            
            # Loop through the algebraic equations and substitute
            # the value in the solution to the ODE
            for index in range(len(list(c_alg))):
                # Solve the algebraic equation for the constant corresponding to the pivot element in B_algebraic
                sol_temp = solve(c_alg[index],const_remaining[pivot_elements_2[index]])
                # Save temporary variables that we can print in a nice fashion
                eq_temp_temp = c_alg[index]
                const_temp = const_remaining[pivot_elements_2[index]]
                sol_temp_temp = sol_temp[0]
                # Substitute to the original variables
                for temp_index in range(len(x)):
                    eq_temp_temp = eq_temp_temp.subs(x[temp_index],variables[temp_index])
                    const_temp = const_temp.subs(x[temp_index],variables[temp_index])
                    sol_temp_temp = sol_temp_temp.subs(x[temp_index],variables[temp_index])
                # Print for the viewer's pleasure 
                print("Equation:\t$%s=0$\t,\tSolution:\t$%s=%s$<br>"%(str(latex(eq_temp_temp)),str(latex(const_temp)),str(latex(sol_temp_temp))))          
                # Substitute the solution of the algebraic equation
                # into the solution of the ODE for the tangential coefficients
                for sub_index in range(len(c_mat)):
                    c_mat[sub_index] = c_mat[sub_index].subs(const_remaining[pivot_elements_2[index]],sol_temp[0])
                    c_mat[sub_index] = expand(cancel(c_mat[sub_index]))
            # Create a readable version in original variables
            c_mat_temp = c_mat
            # Substitute to the original variables
            for index in range(len(x)):
                c_mat_temp = c_mat_temp.subs(x[index],variables[index])            
            print("Solution *after* algebraic substitution:<br>")
            print("\\begin{equation*}\n\mathbf{c}=%s\n\\end{equation*}"%(latex(c_mat_temp)))
            #----------------------------------------------
            # PART 4: Substitute the solution into the tangents
            # and each sub-generator
            #----------------------------------------------
            print("\n\n# The very final step")
            print("The very final step: substituting the solution into the tangents and print the results:<br>")
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
            print("Number of component tangents before removing:\t%d<br>"%(int(len(tangent_component))))
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

