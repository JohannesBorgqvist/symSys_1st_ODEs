#=================================================================================
#=================================================================================
# Script:"launch_symmetry_calculations"
# Date: 2021-06-01
# Implemented by: Johannes Borgqvist
# Description:
# The script launches the symmetry calculations
# for a given model file and a given degree of
# the tangent vector. 
#=================================================================================
#=================================================================================
# Import Libraries
#=================================================================================
#=================================================================================
import conduct_symmetry_calculations # Home-made
import read_and_write_data # Home-made
# To read the input arguments defined by the user
import sys
#=================================================
# FUNCTIONS
#=================================================
# FUNCTION 1: "is_integer"    
# The function excepts a string and returns a
# logical variable. If the provided string is
# an integer it will return "True" else it will
# return "False". 
def is_integer(n):
    try:
        float(n)
    except ValueError:
        return False
    else:
        return float(n).is_integer()
#=================================================
# HANDLE INPUT ARGUMENTS BY USER
#=================================================
# Assign values to the above variables by going
# through the inputs provided by the user:
# Check the number of inputs
if (len(sys.argv)<3):# Too few
    print("\tToo few arguments are provided!\n")
    print("Execute the script as follows:")
    print("python3 launch_symmetry_calculations file_name tangent_degree")
    print("\tAborting execution...\n")
    sys.exit(0) # Exit the program
if (len(sys.argv)>3):# Too many
    print("\tToo many arguments are provided!\n")
    print("python3 launch_symmetry_calculations file_name tangent_degree")
    print("\tAborting execution...\n")
    sys.exit(0) # Exit the program
elif ((len(sys.argv)==3)):# Correct number of inputs
    file_name = str(sys.argv[1])
    tangent_degree = int(sys.argv[2])
    # Is integer provided for the tangential ansätze?
    if (is_integer(sys.argv[2])==False): # No integer provided?
        print("\n\tERROR: The second input must be a positive integer!\n")
        print("\tThe second input corresponds to the degree of the polynomials in the tangentials ansätze.\n")
    else: # An integer was provided
        if tangent_degree < 1: # Provided a silly number
            print("\n\tERROR: The second input must be a positive integer!\n")
            print("\t You have provided a degree of the polynomials in the tangential ansätze which is smaller than 1.\n")
        else:# Ok, now we are ready to conduct the simulations
            # Prompt to the user
            print("==========================================================================================================")            
            print("INITIATING CALCULATIONS WITH INPUT:\n\tmodel\t\t=\t%s,\n\ttangent degree\t=\t%d.\n\n"%(file_name,tangent_degree))
            # Conduct simulations!
            conduct_symmetry_calculations.calculate_symmetries_ODEs(file_name,tangent_degree)
            print("==========================================================================================================")            




