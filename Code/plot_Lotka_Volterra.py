#=================================================================================
#=================================================================================
# Script:"plot_Lotka_Volterra"
# Date: 2022-06-01
# Implemented by: Johannes Borgqvist
# Description:
# The script solves the Lotka-Volterra system and it then plots the solutions
# individually and the phase portrait. We also implement some candidate symmetries
# to compare the effect of these symmetries on the solution curves with just
# changing the symmetry. 
#=================================================================================
#=================================================================================
#=================================================================================
#=================================================================================
# Import Libraries
#=================================================================================
#=================================================================================
#!python
from numpy import *
import pylab as p
from scipy import integrate
#=================================================================================
#=================================================================================
# The Functions
#=================================================================================
#=================================================================================
# Function 1: dX_dt
# This is a help function defining the ODE system we want to solve.
def dX_dt(X, t=0,a=1):
    """ Return the growth rate of fox and rabbit populations. """
    return array([ X[0]*(1-X[1]) ,
                   a*X[1]*(X[0]-1)])
# Function 2: scaling_transformation
# This functions scales trajectories in the state space
def scaling_transformation(rabbits,foxes,epsilon):
    # Defined the transformed rabbits
    rabbits_transformed = asarray([exp(-epsilon)*rabbit for rabbit in rabbits])
    # Define the transformed foxes
    foxes_transformed = asarray([exp(epsilon)*fox for fox in foxes])
    # Return the output
    return rabbits_transformed, foxes_transformed
# Function 3: rotate_solutions
# This functions rotates trajectories in the state space
def rotate_solutions(rabbits,foxes,epsilon):
    # Defined the transformed rabbits
    rabbits_transformed = asarray([rabbit*cos(epsilon)-foxes[index]*sin(epsilon) for index,rabbit in enumerate(rabbits)])
    # Define the transformed foxes
    foxes_transformed = asarray([rabbits[index]*sin(epsilon)+fox*cos(epsilon) for index,fox in enumerate(foxes)])
    # Return the output
    return rabbits_transformed, foxes_transformed    
#=================================================================================
#=================================================================================
# Solving the ODEs
#=================================================================================
#=================================================================================
t = linspace(0, 10, 200)              # time
X0 = array([1, 0.10])                     # initials conditions: 10 rabbits and 5 foxes
#X0_2 = array([20, 10])                     # initials conditions: 20 rabbits and 10 foxes
X1, infodict = integrate.odeint(dX_dt, X0, t, args = (1,),full_output=True)
infodict['message']                     # >>> 'Integration successful.'

#!python
X2, infodict = integrate.odeint(dX_dt, X0, t, args = (2,),full_output=True)
infodict['message']                     # >>> 'Integration successful.'

#!python
rabbits_1, foxes_1 = X1.T
rabbits_2, foxes_2 = X2.T
#=================================================================================
#=================================================================================
# Plotting the solutions
#=================================================================================
#=================================================================================
#Define the first figure
f1, axs = p.subplots(1, 2, constrained_layout=True, figsize=(10, 4))
axs[0].plot(t, rabbits_1, '-', label="Rabits population 1: $\\alpha=1$" ,color=(0/256,68/256,27/256))
axs[0].plot(t, foxes_1  , '-', label='Foxes population 1: $\\alpha=1$',color=(153/256,216/256,201/256))
axs[0].plot(t, rabbits_2, '-', label='Rabits population 2: $\\alpha=2$',color=(77/256,0/256,75/256))
axs[0].plot(t, foxes_2  , '-', label='Foxes population 2: $\\alpha=2$',color=(140/256,150/256,198/256))
axs[0].grid()
axs[0].legend(loc='best',prop={"size":16})
axs[1].plot(rabbits_1, foxes_1,'-', label='Population 1: $\\alpha=1$',color=(0/256,68/256,27/256))
axs[1].plot(rabbits_2,foxes_2, '-', label='Population 2: $\\alpha=2$',color=(77/256,0/256,75/256))
axs[1].grid()
axs[1].legend(loc='best',prop={"size":16})
# Set fontsize of labels
axs[0].set_xlabel(xlabel='Rabbits, $u(t)$',fontsize=25)
axs[0].set_ylabel(ylabel='Foxes, $v(t)$',fontsize=25)
axs[1].set_xlabel(xlabel='Time, $t$',fontsize=25)
axs[1].set_ylabel(ylabel='Population size',fontsize=25)
# Title and saving the figure
f1.suptitle('Evolution of fox and rabbit populations',fontsize=30,weight='bold')
f1.savefig('../Figures/rabbits_and_foxes_dynamics.png')

#=================================================================================
#=================================================================================
# Transforming the solutions
#=================================================================================
#=================================================================================
# Convert everything to numpy arrays
rabbits_1 = asarray(list(rabbits_1))
rabbits_2 = asarray(list(rabbits_2))
foxes_1 = asarray(list(foxes_1))
foxes_2 = asarray(list(foxes_2))
# Scale the rabbits and foxes
rabbits_scaled, foxes_scaled = scaling_transformation(rabbits_1,foxes_1,0.5)
# Rotate the rabbits and foxes
rabbits_rotated, foxes_rotated = rotate_solutions(rabbits_1,foxes_1,pi/3)
#=================================================================================
#=================================================================================
# Plotting the transformed solutions
#=================================================================================
#=================================================================================
#Define the first figure
f1, axs = p.subplots(1, 3, constrained_layout=True, figsize=(10, 4),sharex=True, sharey=True)
# The first subplot
axs[0].plot(rabbits_1, foxes_1,'-', label='Population 1: $\\alpha=1$',color=(0/256,68/256,27/256))
axs[0].plot(rabbits_2,foxes_2, '-', label='Population 2: $\\alpha=2$',color=(77/256,0/256,75/256))
axs[0].grid()
axs[0].legend(loc='best',prop={"size":16})
# The second subplot
axs[1].plot(rabbits_1, foxes_1,'-', label='Population 1: $\\alpha=1$',color=(0/256,68/256,27/256))
axs[1].plot(rabbits_scaled,foxes_scaled, '-', label='Scaled population: $\\epsilon=0.5$',color=(77/256,0/256,75/256))
axs[1].grid()
axs[1].legend(loc='best',prop={"size":16})
# The third subplot
axs[2].plot(rabbits_1, foxes_1,'-', label='Population 1: $\\alpha=1$',color=(0/256,68/256,27/256))
axs[2].plot(rabbits_rotated,foxes_rotated, '-', label='Rotated population: $\\epsilon=\\pi/3$',color=(77/256,0/256,75/256))
axs[2].grid()
axs[2].legend(loc='best',prop={"size":16})
# Set fontsize of labels
axs[0].set_xlabel(xlabel='Rabbits, $u(t)$',fontsize=25)
axs[0].set_ylabel(ylabel='Foxes, $v(t)$',fontsize=25)
axs[0].set_title(label='Correct answer',fontsize=30,style='italic')
axs[1].set_xlabel(xlabel='Rabbits, $u(t)$',fontsize=25)
axs[1].set_ylabel(ylabel='Foxes, $v(t)$',fontsize=25)
axs[1].set_title(label='Scalings: $\\Gamma_{\\epsilon}:(u,v)\\mapsto(e^{-\\epsilon}u,e^{\\epsilon}v)$',fontsize=30,style='italic')
axs[2].set_xlabel(xlabel='Rabbits, $u(t)$',fontsize=25)
axs[2].set_ylabel(ylabel='Foxes, $v(t)$',fontsize=25)
axs[2].set_title(label='Anti-clockwise rotations',fontsize=30,style='italic')
# Title and saving the figure
f1.suptitle('Evolution of fox and rabbit populations',fontsize=30,weight='bold')
p.show()
f1.savefig('../Figures/rabbits_and_foxes_transformations.png')
