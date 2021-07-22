#=================================================================================
#=================================================================================
# Script:"conduct_symmetry_calculations_in_paralell"
# Date: 2021-07-05
# Implemented by: Johannes Borgqvist
# Description:
# The script conducts all the symmetry calculations using functions that are
# stored in the script "symmetry_toolbox_first_order_ODEs.py". This script
# contains merely one function called "calculate_symmetries_ODEs" which reads a
# model from a given data file and attempts to calculate the symmetries.
# This particular script launches all these simulations in paralell using
# paralellisation on CPUs. 
#=================================================================================
#=================================================================================
# Import Libraries
#=================================================================================
#=================================================================================
import read_and_write_data # Home-made
import conduct_symmetry_calculations # Home-made
#from symengine import *
import multiprocessing as mp # For parallelisation over the CPUs
#from multiprocessing import Process # Import processes
import time
from timeit import default_timer as timer
#=================================================================================
#=================================================================================
# The Function
#=================================================================================
#=================================================================================
# The script contains only one function called "calculate_symmetries_ODEs".
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# FUNCTION 1: "calculate_symmetries_ODEs_in_paralell"
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# The function takes two inputs:
# 1. The string "file_name" with the name of the model file,
# 2. The degree of the polynomials in the tangent ansätze called "tangent degree".
# The function stores the calculated generator and all introduced variables in a folder named after the string "file_name" in the output folder with path "../Output/".
def calculate_symmetries_ODEs_in_paralell(i):
    # Define the name of the files and the
    # degree in the tangent ansätze
    if i == 1:
        file_name = "Lotka_Volterra"
        tangent_degree = 1
    elif i == 2:
        file_name = "Lotka_Volterra"
        tangent_degree = 2
    elif i == 3:
        file_name = "Lotka_Volterra"
        tangent_degree = 3
    elif i == 4:
        file_name = "Lotka_Volterra_realistic"
        tangent_degree = 1
    elif i == 5:
        file_name = "Lotka_Volterra_realistic"
        tangent_degree = 2
    elif i == 6:
        file_name = "Lotka_Volterra_realistic"
        tangent_degree = 3        
    
    """
    if i == 1:
        file_name = "hydons_model"
        tangent_degree = 1
    elif i == 2:
        file_name = "DBH_model"
        tangent_degree = 1
    elif i == 3:
        file_name = "DBH_model"
        tangent_degree = 2
    elif i == 4:
        file_name = "DBH_model"
        tangent_degree = 3        
    """

    """elif i == 4:
        file_name = "BZ_model"
        tangent_degree = 3
    elif i == 5:
        file_name = "Goodwin_1"
        tangent_degree = 3
    elif i == 6:
        file_name = "Goodwin_2"
        tangent_degree = 3
    elif i == 7:
        file_name = "Lactose_operon"
        tangent_degree = 3
    elif i == 8:
        file_name = "Lorenz_model"
        tangent_degree = 4
    elif i == 9:
        file_name = "Lotka_Volterra"
        tangent_degree = 3
    elif i == 10:
        file_name = "Lotka_Volterra_realistic"
        tangent_degree = 3
    elif i == 11:
        file_name = "SIR"
        tangent_degree = 3
    elif i == 12:
        file_name = "Brusselator"
        tangent_degree = 3
    elif i == 13:
        file_name = "AIDS_epidemic"
        tangent_degree = 3"""        
    # Conduct the symmetry calculations 
    conduct_symmetry_calculations.calculate_symmetries_ODEs(file_name,tangent_degree)
    return i
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# RUN THE SIMULATIONS IN PARALELL 
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#Define the number of CPUs, as the number of available minus one.
#num_of_CPUs = mp.cpu_count()-1
num_of_CPUs = 6
# Launch symmetry calculations in parallel:
# Set up a pool of workers
#pool = mp.Pool(num_of_CPUs)
#----------------------------------------------------------------------------------
# ATTEMPT USING PROCESSES
#----------------------------------------------------------------------------------
#if __name__ == "__main__":  # confirms that the code is under main function
#    procs = []
#    proc = Process(target=calculate_symmetries_ODEs_in_paralell)  # instantiating without any argument
#    procs.append(proc)
#    proc.start()

    # instantiating process with arguments
#    for i in range(1,14):
        # print(name)
#        proc = Process(target=calculate_symmetries_ODEs_in_paralell, args=(i,))
#        procs.append(proc)
#        proc.start()

    # complete the processes
#    for proc in procs:
#        proc.join()
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# ATTEMPT USING A POOL
#----------------------------------------------------------------------------------    
def main():

    start = timer()

    print(f'starting computations on {num_of_CPUs} cores')

    #values = range(1,14)
    values = range(1,7)

    with mp.Pool(num_of_CPUs) as pool:
        res = pool.map(calculate_symmetries_ODEs_in_paralell, values)
        print(res)

    end = timer()
    print(f'elapsed time: {end - start}')

if __name__ == '__main__':
    main()
