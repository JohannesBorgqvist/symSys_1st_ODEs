#=
=================================================================================
=================================================================================
 Script:"read_and_write_data"
 Date: 2021-07-06
 Implemented by: Johannes Borgqvist
 Description:
 The script contains two functions which
 reads a certain model as input from
 a pre-defined file and then it writes
 the results of the symmetry calculations
 to a folder which it also creates. It
 also contains a help function which
 translates the names of all variables
 from the ones used during the symmetry
 calculations (e.g. (x_0,x_1,...)) to
 the variable names provided in the
 data file. 
=================================================================================
=================================================================================
 Import Libraries
=================================================================================
=================================================================================
=#
using DataFrames, XLSX # For reading the models from a data file
using Symbolics # For doing the calculations
using ModelingToolkit # For some of the variables
#=
=================================================================================
=================================================================================
 Functions
=================================================================================
=================================================================================
=#
#=
 The function takes the name of the model file (which is a string) as input and
 then it looks for an xlsx-files with that name stored in the folder "../Input".
 It has four outputs: a logical variable telling whether or not the reading of the model was successful, the list "variables" which contains the original variables
 as defined in the model file, the list "parameters" containing the parameters
  of the model and the list "omega_list" which contains the reaction terms
 defined in the model file but translated to the generic variables (x0,x1,x2,...)
 where by variables we mean both the independent and dependent variables. The
 script does also create an output directory in the folder "../Output" which
 is later used for storing the data files.
=#
function read_input_model(file_name,sheet_number,variables_symbolic, parameters_symbolic,reaction_terms_symbolic,omega_list)
    # Read the file
    xf = XLSX.readxlsx(file_str)
    # Extract the first sheet name
    sheet_names = XLSX.sheetnames(xf)
    sheet_name = sheet_names[sheet_number,1]
    # Load the sheet as a data frame instead
    df = DataFrame(XLSX.readtable(file_str, sheet_name)...)
    # Calculate the dimensions of the data frame
    num_rows = size(df,1)
    num_cols = size(df,2)
    # Save the names of the columns
    names_df = names(df)
    # Allocate a name for the model
    model_name = "Error in input file"
    # Allcate some memory for our symbolic variables
    x_sym = []        
    # Check if the format is correct of the data file at hand
    if num_cols == 5 && names_df[1] == "Model name" && names_df[2] == "Variable" && names_df[3] == "States" && names_df[4] == "Parameters" && names_df[5] == "Reaction terms"
        # Since the model was correctly defined, let's return the model name
        model_name = df[1,1]
        # Define all symbolic parameters
        # Save the variable
        variable = df[:,2]    
        # Remove alla missing values
        deleteat!(variable, findall(x->typeof(x)==Missing,variable))
        # Loop over the variables and save them
        for (index, value) in enumerate(variable)
            # Define the independent variable
            @eval @variables $(Symbol(string(value)))
            @eval hej = $(Symbol(string(value)))
            # Save the independent variable at hand
            @eval append!(variables_symbolic,[$(Symbol(string(value)))])
        end
        # Save all states
        states = df[:,3]
        # Delete all missing values
        deleteat!(states, findall(x->typeof(x)==Missing,states))
        for (index, value) in enumerate(states)
            # Define the dependent variable
            @eval @variables $(Symbol(string(value)))
            # Save the dependent variable (or the state) at hand
            @eval append!(variables_symbolic,[$(Symbol(string(value)))])
        end    
        # Save all parameters
        parameters = df[:,4]
        deleteat!(parameters, findall(x->typeof(x)==Missing,parameters))
        for (index, value) in enumerate(parameters)
            # Define a variable
            @eval @variables $(Symbol(string(value)))
            # Save the parameter at hand
            @eval append!(parameters_symbolic,[$(Symbol(string(value)))])
        end    
        # Save all reaction terms
        reaction_terms = df[:,5]
        deleteat!(reaction_terms, findall(x->typeof(x)==Missing,reaction_terms))
        for (index, value) in enumerate(variables_symbolic)
            temp_type = typeof(value)
        end
        for (index, value) in enumerate(reaction_terms)
            # Lastly, define a symbolic variable containing the reaction term
            @eval $(Symbol("r"*string(index)))=eval(Meta.parse($(string(value))))
            # Save the reaction term as well
            @eval append!(reaction_terms_symbolic,[$(Symbol(string("r",string(index))))])
        end    
        # Calculate the symbolic variables
        num_var = size(variables_symbolic,1)
        # Create a vector of symbolic variables
        @eval @variables x[1:$(num_var)]
        # Loop over all reaction terms
        for (index, reaction_term) in enumerate(reaction_terms_symbolic)
            # Loop over all variables and conduct the substitution
            for (sub_index, variable) in enumerate(variables_symbolic)
                reaction_term = substitute.(reaction_term, (Dict(variable => x[sub_index]),))
            end
            # Substitute the variables and save the new reaction term
            append!(omega_list,reaction_term)
        end
    else
        # Print an error message
        print("\n\tERROR: INCORRECT MODEL FILE\n")
        print("\tThe data file must have five columns with the following names:\n\n")
        print("\t\t1. Model name\n\t\t2. Variable\n\t\t3. States\n\t\t4. Parameters\n\t\t5. Reaction terms\n\n")
        print("Aborting...\n")
    end
    # Return the output 
    return x_sym,model_name
end
