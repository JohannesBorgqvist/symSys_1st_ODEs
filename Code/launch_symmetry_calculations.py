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
#import conduct_symmetry_calculations_fake # Home-made
import conduct_symmetry_calculations # Home-made
import read_and_write_data # Home-made
# To read the input arguments defined by the user
import sys
# Import the multiprocessing package to handle threads
# which is important in order to terminate script in time
from multiprocessing import Process
# Import the time functionality in order to terminate the script
# after a given number of seconds
import time
# Import the math library to use the rounding function floor
import math
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
# Initially, we assume that we shall not proceed
# with the calculations. Then if all inputs are
# provided in the correct format we will change
# the value of this logical variable to True
# and proceed with the calculations
initiate_calculations = False
# We also have a logical variable stating whether
# or not the script was finished on time
finished_on_time = False
# Assign values to the above variables by going
# through the inputs provided by the user:
# Check the number of inputs
if (len(sys.argv)<3):# Too few
    print("\tToo few arguments are provided!\n")
    print("Execute the script as follows:")
    print("python3 launch_symmetry_calculations file_name tangent_degree")
    print("\tAborting execution...\n")
    sys.exit(0) # Exit the program
if (len(sys.argv)>4):# Too many
    print("\tToo many arguments are provided!\n")
    print("python3 launch_symmetry_calculations file_name tangent_degree")
    print("\tAborting execution...\n")
    sys.exit(0) # Exit the program
elif ((len(sys.argv)==3)):# Correct number of inputs
    # Extract the file name and the tangent degree
    file_name = str(sys.argv[1])
    tangent_degree = int(sys.argv[2])
    # Setting the maximum time to two hours
    maximum_time = 7200
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
            # Calculate the time in hours minutes and so on
            hours = math.floor(maximum_time/3600)
            minutes = maximum_time-(hours*3600)
            minutes = math.floor(((minutes)/(60)))
            seconds = math.floor(maximum_time-(hours*3600) - (minutes*60))
            # Prompt to the user
            print("INITIATING CALCULATIONS WITH INPUT:\n\tmodel\t\t=\t%s,\n\ttangent degree\t=\t%d.\n\tThe computations will be stopped after %d hours, %d minutes and %d seconds.\n\n"%(file_name,tangent_degree,hours,minutes,seconds))
            # Ok, then we can conduct the simulations
            initiate_calculations = True
            # Conduct simulations!
            #conduct_symmetry_calculations.calculate_symmetries_ODEs(file_name,tangent_degree)
elif ((len(sys.argv)==4)):# Correct number of inputs, added a time limit
    # Extract the file name, the tangent degree and the maximum time
    file_name = str(sys.argv[1])
    tangent_degree = int(sys.argv[2])
    maximum_time = int(sys.argv[3])
    # Is integer provided for the tangential ansätze?
    if (is_integer(sys.argv[2])==False): # No integer provided?
        print("\n\tERROR: The second input must be a positive integer!\n")
        print("\tThe second input corresponds to the degree of the polynomials in the tangentials ansätze.\n")
    else: # An integer was provided
        if tangent_degree < 1: # Provided a silly number
            print("\n\tERROR: The second input must be a positive integer!\n")
            print("\t You have provided a degree of the polynomials in the tangential ansätze which is smaller than 1.\n")
        else:
            # Is integer provided for the time limit?
            if (is_integer(sys.argv[3])==False): # No integer provided?
                print("\n\tERROR: The third input must be a positive integer!\n")
                print("\tThe third input corresponds to the maximum time limit for the calculations in seconds.\n")
            else: # An integer was provided
                if maximum_time < 1: # Provided a silly number
                    print("\n\tERROR: The third input must be a positive integer! This is the maximum allowed computation time in seconds. \n")
                    print("\t You have provided a maximum computation time which is smaller than 1 second!\n")
                else:# Ok, now we are ready to conduct the simulations
                    # Prompt to the user
                    print("==========================================================================================================")            
                    # Calculate the time in hours minutes and so on
                    hours = math.floor(maximum_time/3600)
                    minutes = maximum_time-(hours*3600)
                    minutes = math.floor(((minutes)/(60)))
                    seconds = math.floor(maximum_time-(hours*3600) - (minutes*60))
                    # Prompt to the user
                    print("INITIATING CALCULATIONS WITH INPUT:\n\tmodel\t\t=\t%s,\n\ttangent degree\t=\t%d.\n\tThe computations will be stopped after %d hours, %d minutes and %d seconds.\n\n"%(file_name,tangent_degree,hours,minutes,seconds))
                    # Ok, then we can conduct the simulations
                    initiate_calculations = True
                    # Conduct simulations!
                    #conduct_symmetry_calculations.calculate_symmetries_ODEs(file_name,tangent_degree)      

# Okay, we are ready to rock and roll! If the logical variable is true then we will go ahead and
# execute the script. If we are finished on time, the script was succesfully executed and if we
# were not finished on time we will write this to the output file
if __name__ == '__main__':
    # Conduct symmetry calculations
    if initiate_calculations:
        # Define the main thread 
        program = Process(target=conduct_symmetry_calculations.calculate_symmetries_ODEs, args=(file_name,tangent_degree))
        # Start the process
        program.start()
        # End the program after a given amount of time
        program.join(timeout=maximum_time)
        # Terminate the original thread
        program.terminate()
        print("==========================================================================================================")      
