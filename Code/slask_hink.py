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
