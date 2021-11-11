

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# FUNCTION 7: "solve_determining_equations"
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
    # Loop over all ODEs and save the coefficients
    for isolated_eq in det_eq:
        # Expand the equation at hand
        eq_temp = expand(isolated_eq)
        # Loop over all coefficients in the tangential ansätze
        # and define our matrices
        for coeff in c:
            # Save the derivative of the A coefficient
            A_mat.append(-eq_temp.coeff(Derivative(coeff(x[0]),x[0])))
            # Save the coefficient of the coefficient in the tangential
            # ansätze
            B_mat.append(eq_temp.coeff(coeff(x[0])))
    # Finally we generate the matrices from our lists A_mat and B_mat
    num_of_rows = len(det_eq) # Number of rows
    num_of_cols = len(c) # Number of columns
    A = Matrix(num_of_rows,num_of_cols,A_mat) # Matrix A
    B = Matrix(num_of_rows,num_of_cols,B_mat) # Matrix B
    #------------------------------------------------------------------------------
    # STEP 2 of 6: TRANSFER THE SYSTEM TO ONE MATRIX M=[A|-B] AND REDUCE THE NUMBER OF 
    # EQUATIONS BY FINDING A BASIS FOR THE COLUMN SPACE OF M^T. FINISH BY SPLITTING
    # UP THE MATRIX INTO ITS COMPONENTS PART, I.E. A AND B RESPECTIVELY
    #-----------------------------------------------------------------------------
    t0_cs = time.time()
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
    # Split the matrix into its component parts
    A = M_tilde[:,0:(len(c))]
    B = -M_tilde[:,len(c):(2*len(c))]    
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
    #rows_nonzero = [i for i in range(m) if any(A[i, j] != 0 for j in range(n))]
    rows_nonzero = [i for i in range(m) if A[i, j] != 0 for j in range(n)]
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
    # Find the dimensions of the matrix A
    m, n = A.shape
    # Create a set of all the columns
    col_set = set(range(n))
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
    do_non_homo = False
    # Check if we have Pivot columns
    if len(non_pivot_columns)!=0 and do_non_homo:
        # Print the non-pivot columns
        print("Non-Pivot columns")
        print(non_pivot_columns)
        # Indicate that we have a non_homogeneous system
        non_homogeneous = True
        # Define a list of all non-pivot elements.
        # We will refer to this as a non-homogeneous
        # source term
        source_term = []
        source_algebraic = []
        # Introduce a new symbol 
        coefficient_number = int(len(x))
        exec("x_%d = Symbol(\'x_%d\') "%(coefficient_number,coefficient_number))
        temp_variable = []
        exec("temp_variable.append(x_%d)"%(coefficient_number))
        # Allocate for the start time
        tau_0 =  Symbol('tau_0')
        temp_variable.append(tau_0)
        print("Extra symbolic variable:")
        print(temp_variable)
        print("A the non-pivot part")
        print(latex(A[:,non_pivot_columns[0]],mode='equation*'))
        print("B the non-pivot part")
        print(latex(B[:,non_pivot_columns[0]],mode='equation*'))
        print("B_algebraic the non-pivot part")
        print(latex(B_algebraic[:,non_pivot_columns[0]],mode='equation*'))        
        print("The coefficients")
        print(latex(Matrix(len(c),1,c)))
        # Calculate the rows and columns in A
        r_algebraic, cols = B_algebraic.shape
        rows, cols = A.shape
        print("Dimensions of A:\t%dX%d\n"%(rows,cols))
        print("Dimensions of B_algebraic:\t%dX%d\n"%(r_algebraic,cols))        
        # DEFINE SOURCE TERM FOR ODEs
        for r in range(rows):
            # Define a temporary source term
            temp_source = 0
            # Loop over the pivot columns
            for pC in non_pivot_columns:
                temp_source += -A[r,pC]*Derivative(c[pC](x[0]),x[0])
                temp_source += B[r,pC]*c[pC](x[0])
            # Append the temporary source
            source_term.append(temp_source)
        # DEFINE SOURCE TERM FOR ALGEBRAIC EQUATIONS
        for r in range(r_algebraic):
            # Define a temporary source term
            temp_source_alg = 0
            # Loop over the pivot columns and perform
            # the matrix multiplication
            for pC in non_pivot_columns:
                temp_source_alg += -B_algebraic[r,pC]*c[pC](x[0])
            # Append the non-homogeneous source 
            source_algebraic.append(temp_source_alg)
        # Print the matrix
        print("The non-homogeneous source term:<br>")
        print(latex(Matrix(len(source_term),1,source_term),mode='equation*'))
        source_matrix = Matrix(len(source_term),1,source_term)
        print("Dimensions:<br>")
        print(source_matrix.shape)
        print("The algebraic source term:<br>")
        source_alg_matrix = Matrix(len(source_algebraic),1,source_algebraic)
        print(latex(source_alg_matrix,mode='equation*'))
        print("Dimensions:<br>")
        print(source_alg_matrix.shape)
    # Remove all non-pivot tangential coefficient
    non_pivot_columns.sort(reverse=True)
    # Loop over all non_pivot_columns and remove them from our three
    # matrices A, B and B_algebraic and from the coefficient vector c
    for index in range(len(non_pivot_columns)):
        # Remove column matrices
        A.col_del(non_pivot_columns[index])
        B.col_del(non_pivot_columns[index])
        B_algebraic.col_del(non_pivot_columns[index])
        # Remove parameter coefficient vector
        del c[non_pivot_columns[index]]
    #print("The condition number of $B$: \n")
    #print(B.condition_number())
    #------------------------------------------------------------------------------
    # STEP 5 of 6: SOLVE THE ODE SYSTEM PROVIDED THAT IT IS QUADRATIC AND
    # SUBSTITUTE THE ALGEBRAIC EQUATIONS
    #-----------------------------------------------------------------------------
    num_rows, num_cols = B.shape
    if num_rows == num_cols:
        # Check that the matrix A is an identity matrix
        if A != eye(num_rows):
            X = "\\Huge\\textsf{$A$ is quadratic but not an identity matrix!}\normalsize\\[2cm]" 
        else: # A is an identity matrix
            #----------------------------------------------
            # PART 1: SOLVE SYSTEM OF ODEs
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
            # Calculate the Jordan form of the matrix B
            P,J = B.jordan_form()
            # Re-define our matrix J by scaling it by the
            # independent variable x[0]
            J = x[0]*J
            # Re-define J by taking the matrix exponential
            J = J.exp()
            jordan_copy = J.inv()
            # Solve the system of ODEs: the homogenous solution
            c_mat = (P*J*P.inv())*c_mat
            if non_homogeneous:
                non_homo = (P*jordan_copy*P.inv())*source_matrix
                #non_homo = (P*J*P.inv())*source_matrix
                m,n = non_homo.shape               
                for rI in range(m):
                    for cI in range(n):
                        #non_homo[rI,cI] = integrate(non_homo[rI,cI],x[0])
                        non_homo[rI,cI] = Integral(non_homo[rI,cI],x[0])
                        non_homo[rI,cI] = non_homo[rI,cI].doit()
                # Multiply with J again
                non_homo = (P*J*P.inv())*non_homo
                # Save the homogeneous solution as well
                homo = c_mat
                # Add the non-homogeneous part to the total solution
                c_mat += non_homo
                        
            #----------------------------------------------
            # PART 2: SUBSTITUE ALGEBRAIC EQUATIONS
            #----------------------------------------------            
            # Calculate the pivot elements of the matrix
            # B_algebraic
            B_algebraic = B_algebraic.rref()[0]
            pivot_elements = B_algebraic.rref()[1]
            # Define the algebraic equations for the
            # coefficients
            if non_homogeneous:
                c_alg = B_algebraic*homo
            else:
                c_alg = B_algebraic*c_mat
            # Define a list of remaining constants
            const_remaining = []
            # Define the rows in c_alg
            alg_row,alg_col = c_alg.shape
            # Extract remaining coefficients in our algebraic
            # equation corresponding to c_alg = 0
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
            # Print
            const_remaining_mat = Matrix(len(const_remaining),1,const_remaining)
            print("Consistency check to see that our matrix definition is correct:")
            print("c_alg-alg_mat\*const_remaining_mat=?")
            print(latex(c_alg-(alg_mat*const_remaining_mat),mode='equation*'))
            # Print the solution
            print("Solution:<br>")
            print(latex(c_mat,mode='equation*'))
            print("Algebraic matrix:<br>")
            print(latex(B_algebraic,mode='equation*'))
            print("Algebraic equations: c_alg=B_algebraic*c_mat")
            print(latex(c_alg,mode='equation*'))
            # Non-homogeneos crap
            if non_homogeneous:
                print("Non_homo<br>")
                print(latex(non_homo,mode='equation'))
                print("Dimensions:<br>")
                print(non_homo.shape)
                print("-B_algebraic*non_homo<br>")
                print(latex(-B_algebraic*non_homo,mode='equation'))
                print("Dimensions:<br>")
                print((-B_algebraic*non_homo).shape)                
                print("source alg mat<br>")
                print(latex(source_alg_matrix,mode='equation'))
                print("Dimensions:<br>")
                print(source_alg_matrix.shape)
                # Define the right hand side
                RHS = -B_algebraic*non_homo + source_alg_matrix
                # Create the extendent matrix by adding the RHS
                # to alg_mat
                extended_mat = alg_mat
                rows_alg, cols_alg = alg_mat.shape
                extended_mat = extended_mat.col_insert(cols_alg,RHS)
                print("alg_mat<br>\n")
                print(latex(alg_mat,mode='equation*'))
                print("Dimensions:<br>")
                print(alg_mat.shape)                
                print("RHS<br>\n")
                print(latex(RHS,mode='equation*'))
                print("Dimensions:<br>")
                print(RHS.shape)                                
                print("extended_mat<br>\n")
                print(latex(extended_mat,mode='equation*'))
                print("Dimensions:<br>")
                print(extended_mat.shape)                                
                # Row reduce the extended matrix
                extended_mat = extended_mat.rref()[0]
                pivot_elements_2 = extended_mat.rref()[1]
                print("extended_mat after row reduction<br>\n")
                print(latex(extended_mat,mode='equation*'))
                print("The pivot elements:<br>")
                print(pivot_elements_2)
                # Extract the matrix and the non-homogeneous part
                print("Alg mat after row reduction<br>")
                alg_mat = extended_mat[:,0:cols_alg]
                print(latex(alg_mat,mode='equation*'))
                print("Dimensions:<br>")
                print(alg_mat.shape)
                print("The right hand side")
                RHS = extended_mat[:,-1]
                print(latex(RHS,mode='equation*'))
                print("Dimensions:<br>")
                print(RHS.shape)                
            else:
                # Row reduce the same matrix
                alg_mat = alg_mat.rref()[0]
                pivot_elements_2 = alg_mat.rref()[1]            
            # Update the algebraic equations
            if non_homogeneous:
                c_alg = alg_mat*const_remaining_mat - RHS
            else:
                c_alg = alg_mat*const_remaining_mat
            # Convert the matrix with the ODE coefficients into
            # a list
            c_list = list(c_mat)

            print("Remaining constants:")
            print(latex(const_remaining_mat,mode='equation*'))
            print("Dimensions")
            print(const_remaining_mat.shape)
            print("Algebraic equations")
            print("\\begin{equation*}")
            print("%s=%s%s"%(str(latex(c_alg)),str(latex(alg_mat)),str(latex(const_remaining_mat))))
            print("\\end{equation*}")
            #print(latex(c_alg,mode='equation*'))
            print("Dimensions")
            print(c_alg.shape)
            print("Pivot elements")
            print(pivot_elements_2)
            print("SOLVING EQUATIONS, HEY?")
            # Loop through the algebraic equations and substitute
            # the value in the solution to the ODE
            for index in range(len(list(c_alg))):
                # Solve the algebraic equation for the constant corresponding
                # to the pivot element in B_algebraic
                #sol_temp = solve(c_alg[index],constant_list[pivot_elements[index]])
                #print("Index\t%d\n"%(int(index)))
                #print(c_alg[index])
                sol_temp = solve(c_alg[index],const_remaining[pivot_elements_2[index]])
                print("Equation:\t$%s=0$\t,\tSolution:\t$%s=%s$<br>"%(str(latex(c_alg[index])),str(latex(const_remaining[pivot_elements_2[index]])),str(latex(sol_temp))))                
                #print(sol_temp)
                #print("\n")
                # Substitute the solution of the algebraic equation
                # into the solution of the ODE for the tangential coefficients
                for sub_index in range(len(c_list)):
                    #c_list[sub_index] = c_list[sub_index].subs(constant_list[pivot_elements[index]],sol_temp[0])
                    c_list[sub_index] = c_list[sub_index].subs(const_remaining[pivot_elements_2[index]],sol_temp[0])
                    c_list[sub_index] = expand(cancel(c_list[sub_index]))
                    #for state_index in range(len(x)):
                        #c_list[sub_index] = c_list[sub_index].subs(x[state_index],variables[state_index])
            print("Solution before algebraic substitution:<br>")
            print(latex(c_mat))
            print("Solution after algebraic substitution:<br>")
            print(latex(Matrix(len(c_list),1,c_list),mode='equation*'))
            # See which constants that exists among the constants
            constants_in_final_solution = []
            # Loop over all algebraic expressions
            for alg_index in range(len(c_list)):
                # Loop over all candidate constants
                for constant_index in range(len(constant_list)):
                    #print("Expression:\t%s,\tConstant:\t%s,\tCoeff:\t%s"%(str(c_list[alg_index]),str(constant_list[constant_index]),str(c_list[alg_index].coeff(constant_list[constant_index]))))
                    # See if the constant exists in the final solution
                    if c_list[alg_index].coeff(constant_list[constant_index]) != 0:
                        # Add the constant
                        constants_in_final_solution.append(constant_list[constant_index])
                        # Move on
                        #continue
            # Only save the unique ones
            constants_in_final_solution = list(set(constants_in_final_solution))
            print("constants_in_final_solution")
            print(latex(Matrix(len(constants_in_final_solution),1,constants_in_final_solution),mode='equation*'))
            print("Number of arbitrary constants:\t%d\n"%(len(constants_in_final_solution)))
            # Loop over tangents and substitute the calculated coefficients
            # into these expressions
            for tangent_number in range(len(eta_list)):
                # Expand the tangents
                eta_list[tangent_number] = eta_list[tangent_number].expand()
                # Loop through all coefficients and replace them in all tangents
                for index in range(len(c)):            
                    # Substitute coefficient in the tangents
                    eta_list[tangent_number] = eta_list[tangent_number].subs(c[index](x[0]),c_list[index])                    
                    eta_list[tangent_number] = expand(cancel(eta_list[tangent_number]))
                #print("Tangent %d"%(tangent_number))
                #print(latex(eta_list[tangent_number]))
            print("The tangents then:<br>")
            print("\\begin{align*}")
            counter = 0
            for eta in eta_list:
                print("eta_%d&=%s\\\\"%(int(counter),str(latex(eta))))
                counter += 1
            print("\\end{align*}")
            # Make space for the non-homogoneous tangent as well
            non_homo_tangent = []
            # Add zeros to the non-homogeneous tangent at first
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
                    # Add the term in order to calculate the non-homogeneous part
                    non_homo_tangent[tangent_number] += eta_list[tangent_number].coeff(constant)*constant
                # Append the extracted generator to the list of tangents
                tangent_component.append(temp_eta)
            # Print component tangents
            print("Tangent component:")
            print(latex(tangent_component))
            # Print constants in final solution
            print("Constants in final solution:")
            print(latex(constants_in_final_solution))
            print("Non-homo tangent")
            print(non_homo_tangent)
            print("eta")
            print(eta_list)
            # Finally calculate the non-homogeneous tangent
            for non_homo_index in range(len(non_homo_tangent)):
                non_homo_tangent[non_homo_index] = expand(eta_list[non_homo_index]) - expand(non_homo_tangent[non_homo_index])
                non_homo_tangent[non_homo_index] = simplify(cancel(non_homo_tangent[non_homo_index]))
            # Now, we append the non-homogeneous tangent
            tangent_component.append(non_homo_tangent)
            print("Number of component tangents:\t%d"%(int(len(tangent_component))))
            # For clarity, let's print the non-homogeneous tangent
            print("Non-homogeneous tangent:")
            temp_str = "$$X="
            for non_homo_index in range(len(non_homo_tangent)):
                temp_str += str(latex(non_homo_tangent[non_homo_index])) + "\partial "+ str(latex(x[non_homo_index]))
                if non_homo_index != 0 or non_homo_index != len(non_homo_tangent)-1:
                    temp_str += "+"
            temp_str += "$$"
            print("%s"%(temp_str))
            print("Linearised symmetry condition:")
            temp_list = lin_sym_cond(x, non_homo_tangent, omega_list)
            temp_list = [simplify(cancel(expand(ls))) for ls in temp_list]
            print("\\begin{align*}")
            for lin_sym in temp_list:
                print("%s&=0\\\\"%(str(latex(lin_sym))))
            print("\\end{align*}")
            # Check the linearised symmetry conditions
            # Allocate a vector of the linearised symmetry conditions
            lin_sym_index = []
            # Loop over all sub tangents and check the linearised symmetry conditions
            for tangent_index in range(len(tangent_component)):
                # Calculate the symmetry conditions for the tangent at hand
                temp_list = lin_sym_cond(x, tangent_component[tangent_index], omega_list)
                print("Tangent:<br>")
                print(latex(tangent_component[tangent_index],mode='equation*'))
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
            # Add the final touch to the generator
            X += "\\end{align*}"
            X = X.replace("\\\\\n\\end{align*}", ".\\\\\n\\end{align*}")
    else:
        # Return that the matrix
        X = "\\Huge\\textsf{Not a quadratic matrix}\\normalsize\\\\[2cm]" 
    # Return the solved coefficients and the generator
    return X














### VERY LAST BIT HEY?

#----------------------------------------------
            # PART 4: Substitute the solution into the tangents
            # and each sub-generator
            #----------------------------------------------
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
            print("constants_in_final_solution")
            print(latex(Matrix(len(constants_in_final_solution),1,constants_in_final_solution),mode='equation*'))
            # Loop over tangents and substitute the calculated coefficients
            # into these expressions
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
                    # The non-homogeneous case: Save all terms involving the coefficients
                    # so that we can substract it from the original tangent later in order
                    # to calculate the extra terms.
                    if non_homogeneous:
                        # Add the term in order to calculate the non-homogeneous part
                        non_homo_tangent[tangent_number] += eta_list[tangent_number].coeff(constant)*constant
                # Append the extracted generator to the list of tangents
                tangent_component.append(temp_eta)
            # Calculate the non-homogeneous tangent if such a tangent exist
            if non_homogeneous:
                # Loop over the non-homogenous terms and calculate the part of the tangent
                # that is not associated with an arbitrary integration constant
                for non_homo_index in range(len(non_homo_tangent)):
                    non_homo_tangent[non_homo_index] = expand(eta_list[non_homo_index]) - expand(non_homo_tangent[non_homo_index])
                    non_homo_tangent[non_homo_index] = simplify(cancel(non_homo_tangent[non_homo_index]))
                # Append the non-homogeneous tangent to the list of tangents
                tangent_component.append(non_homo_tangent)
            #----------------------------------------------
            # PART 5: Check if the generators satisfy
            # the linearised symmetry conditions
            #----------------------------------------------
            print("Number of component tangents before removing:\t%d"%(int(len(tangent_component))))
            # Check the linearised symmetry conditions
            # Allocate a vector of the linearised symmetry conditions
            lin_sym_index = []
            # Loop over all sub tangents and check the linearised symmetry conditions
            for tangent_index in range(len(tangent_component)):
                # Calculate the symmetry conditions for the tangent at hand
                temp_list = lin_sym_cond(x, tangent_component[tangent_index], omega_list)
                print("Tangent:<br>")
                print(latex(tangent_component[tangent_index],mode='equation*'))
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
            # Add the final touch to the generator
            X += "\\end{align*}"
            X = X.replace("\\\\\n\\end{align*}", ".\\\\\n\\end{align*}")




# ALGEBRAIC BIT FOR LATERS
            #----------------------------------------------
            # PART 2: Solve algebraic equations
            #----------------------------------------------
            # Row reduce the same matrix
            B_algebraic = B_algebraic.rref()[0]
            pivot_alg = B_algebraic.rref()[1]
            print("Pivot elements:\t%s<br>"%(str(pivot_alg)))
            c_alg = B_algebraic*c_mat
            print("Algebraic equations:<br>")
            print(latex(c_alg,mode='equation'))
            # Loop over the algebraic equations
            for index in range(len(list(c_alg))):
                # Non-homogeneous equation
                if non_homogeneous:
                    # Solve the equation for the pivot coefficient
                    sol_temp = solve(c_alg[index]-source_alg[index],c_mat[pivot_alg[index]])                    
                else: # Homogeneous equation
                    sol_temp = solve(c_alg[index],c_mat[pivot_alg[index]])
                # Loop over all coefficients and substitute the value
                for sub_index in range(len(c_mat)):
                    c_mat[sub_index] = c_mat[sub_index].subs(c_mat[pivot_alg[index]],sol_temp[0])
                    c_mat[sub_index] = expand(cancel(c_mat[sub_index]))
            print("Coefficients:<br>")
            print(latex(c_mat,mode='equation'))
    #-----------------------------------------------------------------------------
        #-----------------------------------------------------------------------------
            #-----------------------------------------------------------------------------

            # Loop over tangents and substitute the calculated coefficients into these expressions
            for tangent_number in range(len(eta_list)):
                # Expand the tangents
                eta_list[tangent_number] = eta_list[tangent_number].expand()
                # Loop through all coefficients and replace them in all tangents
                for index in range(len(c)):            
                    # Substitute coefficient in the tangents
                    #eta_list[tangent_number] = substitute_tangents(eta_list[tangent_number],c[index](x[0]),c_mat[index])
                    #if eta_list[tangent_number].coeff(c[index](x[0]))!=0:
                    temp_expression = eta_list[tangent_number]
                    temp_expression = temp_expression.subs(c[index](x[0]),cancel(expand(simplify(c_mat[index]))))                    
                    eta_list[tangent_number] = expand(cancel(temp_expression))


















#==================================================================
# FIXING THE PRINTING OF THE GENERATOR
#==================================================================
                # Increment the counter for the generator
                generator_counter += 1
                # Set the plus indicator to false
                plus_indicator = False
                # Start
                X += "\\begin{align*}\nX_{" + str(generator_counter) + "}&="
                # Re-set the term_counter
                term_counter = 1
                # Loop over each coordinate in the tangent vector   
                for tangent_part in range(len(tangent)):
                    # Add all non-trivial tangents to the generator    
                    if tangent[tangent_part]!=0:
                        # Replace the name of all arbitrary functions
                        for arb_func_ind in range(len(non_pivot_functions)):
                            tangent[tangent_part] = tangent[tangent_part].subs(non_pivot_functions[arb_func_ind],arbitrary_functions[arb_func_ind])
                        # Write the tangent to the string
                        if plus_indicator:
                            X += "+\\left( "                                
                        else:
                            #X += "\\left( " + str(latex(tangent[tangent_part])) + " \\right)" + "\\partial " + str(latex(variables[tangent_part])) + ""
                            X += "\\left( "                                
                            plus_indicator = True
                        # Chop the tangent into its pieces
                        chopped_generator = latex(tangent[tangent_part]).split("+")
                        # Loop over the pieces of the tangent
                        for term_index in range(len(chopped_generator)):
                            # Add the terms one by one
                            if term_index == (len(chopped_generator)-1):
                                # The very last generator or the other ones
                                X += chopped_generator[term_index] + " \\right)" + "\\partial " + str(latex(variables[tangent_part])) + ""
                                # Reset the term counter if this happens?
                                #term_counter = 1
                            else: # Just add the term
                                X += chopped_generator[term_index] + "+"                                    
                            # Increment the term counter
                            term_counter += 1
                            # Add a line break if we have too many rows
                            if term_counter > term_tolerance:
                                # Re-set the term_counter
                                term_counter =1
                                # Two cases if we end a line or not
                                if term_index == (len(chopped_generator)-1):
                                    # Add a new row
                                    X += "\\\\\n"
                                else: # Not the last term then we need to add a \left...    
                                    # Add a new row
                                    X += "\\right.\\\\\n\\left.&+"
                # End the alignment                    
                X += "\\end{align*}\n"                   
                            
        # Add a line break
        #X += ",\\\\\n"
        # Re-set the term counter when this is done
        #term_counter = 1
        # Add a new line break
        #line_counter += 1        
        # Add the final touch to the generator
        #if line_counter != 1:
            #X += "\\end{align*}"
        #X = X.replace(",\\\\\n\\end{align*}", ".\\\\\n\\end{align*}")
        #









                        # Start the alignment
                X += "\\begin{align*}\nX_{" + str(generator_counter) + "}&="
                # Set the plus indicator to false
                plus_indicator = False
                # Re-set the term_counter
                term_counter = 1                
                # Loop over each coordinate in the tangent vector   
                for tangent_part in range(len(tangent)):
                    # Add all non-trivial tangents to the generator    
                    if tangent[tangent_part]!=0:
                        # Write the tangent to the string
                        if plus_indicator:
                            X += "+\\left( "                                
                        else:
                            X += "\\left( "                                
                        # Chop the tangent into its pieces
                        chopped_generator = latex(tangent[tangent_part]).split("+")
                        # Loop over the pieces of the tangent and add them
                        for term_index in range(len(chopped_generator)):
                            # Add the next term
                            X += chopped_generator[term_index]
                            # Add the terms one by one
                            if term_index == (len(chopped_generator)-1):
                                # The very last generator or the other ones
                                X += " \\right)" + "\\partial " + str(latex(variables[tangent_part])) + "+"
                            else: # Just add the term
                                X += "+"
                            # Increment the term counter
                            term_counter += 1
                            print("Term counter\t=\t%d"%(term_counter))
                            # Add a line break if we have too many rows
                            if term_counter > term_tolerance:
                                # Re-set the term_counter
                                term_counter =1
                                if term_index == (len(chopped_generator)-1):                          
                                    X += "\\\\"
                                elif term_index<(len(chopped_generator)-1):
                                    X += "\\right.\\\\\n&\\left."
                # End the alignment                    
                X += "\n\\end{align*}\n"












                        if len(basis_functions[base_index].atoms(Symbol))==0:
            constant_basis_list.append(base_index)
        else:
            # Make sure we do not have a negative sign in front
            # of the basis element for all polynomial bases
            num, denom = fraction(basis_functions[base_index])
            # We have a polynomial type thingy (monomial or exponential).
            # Then, we can calculate the sign through factor_list. If we have
            # A rational thingy (e.g. 1/x[0]) then factor_list does not work.
            # This is why we calculate the denominator before hand. 
            if denom == 1:
                # Calculate the factors in this polynomial which will give
                # us the sign
                factors = factor_list(basis_functions[base_index])
                #sign_factor = simplify(basis_functions[base_index]/Abs(basis_functions[base_index]))
                #basis_functions[base_index] = sign_factor * basis_functions[base_index]
                # Change the sign of the negative basis elements
                #if factors[0]<0:
                    #basis_functions[base_index] = -1 * basis_functions[base_index]
                    #if factors[0]==-1:
                    #    basis_functions[base_index] = -1 * basis_functions[base_index]
                    #else:
                basis_functions[base_index] = ((1)/(factors[0])) * basis_functions[base_index]
            else:
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







    #==============================================================================
    # We define the equation which is the constant
    # in front of the basis function
    eq_temp = expr.coeff(func_temp)
    # Ok, let's remove all terms which has a dependence
    # on the independent variable.
    # Allocate memory for the dependent terms which we will remove
    dependent_terms = 0
    # Loop over the terms in our equation at hand
    for term in eq_temp.args:
        # Save the one with a dependence on the
        # independent variable
        if Derivative(term,x[0]).doit()!=0: # We have a function!
            dependent_terms += term
            # Now, we update our equation and remove the dependent terms
            eq_temp = eq_temp - dependent_terms
            #==============================================================================
