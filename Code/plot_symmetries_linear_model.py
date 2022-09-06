#=================================================================================
#=================================================================================
# Script:"plot_symmetries_linear_model"
# Date: 2022-06-28
# Implemented by: Johannes Borgqvist
# Description:
# The script visualises the symmetries of the linear model
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
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.optimize import fsolve
#=================================================================================
#=================================================================================
# The Functions
#=================================================================================
#=================================================================================
# Function 1: dX_dt
# This is a help function defining the ODE system we want to solve.
def dX_dt(X, t=0,*parameters):
    # Extract the parameters
    a, b, c, d = parameters
    # Return the dynamics of the linear system
    return array([ a*X[0]+b*X[1],
                   c*X[0]+d*X[1]])
# Function 2: The symmetry in the r-direction with constant theta
def Gamma_r(p_hat,*arguments):
    # Define our transformed point that we want to find
    u_hat, v_hat = p_hat
    # Extract our arguments
    u,v,epsilon = arguments
    # Define the r-coordinates
    r = arctan(v/u)
    r_hat = arctan(v_hat/u_hat)
    # Define the s-coordinates
    s = (1/2)*log(v**2+u**2)
    s_hat = (1/2)*log(v_hat**2+u_hat**2)
    # Finally, return our equations that we wish to solve
    return (r-r_hat, s+epsilon-s_hat)
# Function 2: The symmetry in the theta-direction with constant r
def Gamma_theta(p_hat,*arguments):
    # Define our transformed point that we want to find
    u_hat, v_hat = p_hat
    # Extract our arguments
    u,v,epsilon,a,b,c,d = arguments
    # Define the r-coordinates
    r = sqrt(v**2+u**2)
    r_hat = sqrt(v_hat**2+u_hat**2)
    # Define our angles as well
    theta = arctan(v/u)
    theta_hat = arctan(v_hat/u_hat)
    # Define the trace, determinant and discriminant
    trace = a + d
    determinant = a*d - b*c
    discriminant = sqrt(trace**2-4*determinant)
    # Define the s-coordinates, but we do it in pieces
    s = ((trace)/(discriminant))*arctanh(((a+2*b*tan(theta)-d)/(discriminant)))
    s += -(1/2)*log(abs((d-a)*sin(2*theta)+(b+c)*cos(2*theta) - b + c))
    s_hat = ((trace)/(discriminant))*arctanh(((a+2*b*tan(theta_hat)-d)/(discriminant)))
    s_hat += -(1/2)*log(abs((d-a)*sin(2*theta_hat)+(b+c)*cos(2*theta_hat) - b + c))  
    # Finally, return our equations that we wish to solve
    return (r-r_hat, s+epsilon-s_hat)
# Function 3: Gamma_r_ODE
def Gamma_r_ODE(X, epsilon=0,*parameters):
    # Extract the parameters
    a, b, c, d = parameters
    # Return the dynamics of the linear system
    return array([ X[0],X[1]])
                   
# Function 4: Gamma_theta_ODE
def Gamma_theta_ODE(X, epsilon=0,*parameters):
    # Extract the parameters
    a, b, c, d = parameters
    # Return the dynamics of the linear system
    return array([ -((c*X[0]**2+(d-a)*X[0]*X[1]-b*X[1]**2)/(a*X[0]**2+(b+c)*X[0]*X[1]+d*X[1]**2))*X[1],
                   ((c*X[0]**2+(d-a)*X[0]*X[1]-b*X[1]**2)/(a*X[0]**2+(b+c)*X[0]*X[1]+d*X[1]**2))*X[0]])
# Function 5: dX_dt_os
def dX_dt_os(X, t=0,*parameters):
    # Extract the parameters
    omega = parameters[0]
    # Define the radius
    r = sqrt(X[0]**2+X[1]**2)    
    # Return the dynamics of the linear system
    return array([ X[0]*(1-r)-omega*X[1],
                   X[1]*(1-r)+omega*X[0]])
# Function 6: Gamma_r_os
def Gamma_r_os(X, epsilon=0,*parameters):
    # Extract the parameters
    omega = parameters[0]
    # Define the radius
    r = sqrt(X[0]**2+X[1]**2)
    # Return the dynamics of the biological oscillator
    return array([(1-r)*X[0],(1-r)*X[1]])
                   
# Function 7: Gamma_theta_os
def Gamma_theta_os(X, epsilon=0,*parameters):
    # Extract the parameters
    omega = parameters[0]
    # Return the dynamics of the biological oscillator
    return array([ -(X[1]/omega),(X[0]/omega)])
# Function 8: plot_LaTeX_2D
# This functions enable us to reproduce our plots using pgfplots in LaTeX
def plot_LaTeX_2D(t,y,file_str,plot_str,legend_str):
    # Open a file with the append option
    # so that we can write to the same
    # file multiple times
    f = open(file_str, "a")
    # Create a temporary string which
    # is the one that does the plotting.
    # Here we incorporate the input plot_str
    # which contains the color, and the markers
    # of the plot at hand
    if len(legend_str)==0:
        temp_str = "\\addplot[\nforget plot,\n" + plot_str+ "\n]\n"
    else:
        temp_str = "\\addplot[\n" + plot_str+ "\n]\n"
    # Add the coordinates
    temp_str += "coordinates {%\n"
    # Loop over the input files and add
    # them to the file
    for i in range(len(t)):
        temp_str += "(" + str(t[i]) + "," + str(y[i]) + ")\n"
    # The plotting is done, let's close the shop    
    temp_str += "};\n"
    # Add a legend if one is provided
    if len(legend_str) > 0:
        temp_str += "\\addlegendentry{" + legend_str + "}\n"
    # Finally, we write the huge string
    # we have created
    f.write("%s"%(temp_str))
    # Close the file
    f.close()
#=================================================================================
#=================================================================================
# Transforming the linear model with our symmetries
#=================================================================================
#=================================================================================
# ---------------------------------------------------------------------------
# Overall properties (a,b,c,d)=(-1,0,0,-2) Stable node
# ---------------------------------------------------------------------------
# Define the parameters for the stable point
a_stable = -1
b_stable = 0
c_stable = 0
d_stable = -2
# Transformation parameters
epsilon_r=0.33
epsilon_theta = pi/6
# Define our initial condition
u0 = 1
v0 = 3
# Define the time vector
t = linspace(0, 2, 500)              # Time
X0 = array([u0, v0])                  # ICs
X1, infodict = integrate.odeint(dX_dt, X0, t, args = (a_stable,b_stable,c_stable,d_stable),full_output=True)
infodict['message']                     # >>> 'Integration successful.'
# Extract the original solution with the defined parameters and initial conditions
u_1, v_1 = X1.T
# The indices we wish to transform
lin_indices = arange(0,len(t),15)
# ---------------------------------------------------------------------------
# Stable node (a,b,c,d)=(-1,0,0,-2) radial symmetry
# ---------------------------------------------------------------------------
# Plot the action of the symmetry
epsilon_vec = arange(0,epsilon_r-0.005,epsilon_r/50)
# Allocate our empty lists
Gamma_trans_r_u_1 = []
Gamma_trans_r_v_1 = []
# Loop over our indices and find the transformations
for lin_index in lin_indices:
    # Get our transformations
    X_r, infodict = integrate.odeint(Gamma_r_ODE, array([u_1[lin_index],v_1[lin_index]]), epsilon_vec, args = (a_stable,b_stable,c_stable,d_stable),full_output=True)
    infodict['message']                     # >>> 'Integration successful.'
    # Extract transformations
    Gamma_trans_r_u_temp, Gamma_trans_r_v_temp = X_r.T
    # Save our transformations
    Gamma_trans_r_u_1.append(Gamma_trans_r_u_temp)
    Gamma_trans_r_v_1.append(Gamma_trans_r_v_temp)    
# Solve to get the transformed angle solution
X1_r, infodict = integrate.odeint(dX_dt, array([Gamma_trans_r_u_1[0][-1], Gamma_trans_r_v_1[0][-1]]), t, args = (a_stable,b_stable,c_stable,d_stable),full_output=True)
u_r_1, v_r_1 = X1_r.T
# Allocate our empty lists
Gamma_trans_r_u_2 = []
Gamma_trans_r_v_2 = []
# Loop over our indices and find the transformations
for lin_index in lin_indices:
    # Get our transformations
    X_r, infodict = integrate.odeint(Gamma_r_ODE, array([u_r_1[lin_index],v_r_1[lin_index]]), epsilon_vec, args = (a_stable,b_stable,c_stable,d_stable),full_output=True)
    infodict['message']                     # >>> 'Integration successful.'
    # Extract transformations
    Gamma_trans_r_u_temp, Gamma_trans_r_v_temp = X_r.T
    # Save our transformations
    Gamma_trans_r_u_2.append(Gamma_trans_r_u_temp)
    Gamma_trans_r_v_2.append(Gamma_trans_r_v_temp)    
# Solve to get the transformed angle solution
X2_r, infodict = integrate.odeint(dX_dt, array([Gamma_trans_r_u_2[0][-1], Gamma_trans_r_v_2[0][-1]]), t, args = (a_stable,b_stable,c_stable,d_stable),full_output=True)
u_r_2, v_r_2 = X2_r.T
# ---------------------------------------------------------------------------
# Stable node (a,b,c,d)=(-1,0,0,-2) angle symmetry
# ---------------------------------------------------------------------------
# Plot the action of the symmetry
epsilon_vec = arange(0,epsilon_theta-0.005,epsilon_theta/50)
# Allocate our empty lists
Gamma_trans_theta_u_1 = []
Gamma_trans_theta_v_1 = []
# Loop over our indices and find the transformations
for lin_index in lin_indices:
    # Get our transformations
    X_theta, infodict = integrate.odeint(Gamma_theta_ODE, array([u_1[lin_index],v_1[lin_index]]), epsilon_vec, args = (a_stable,b_stable,c_stable,d_stable),full_output=True)
    infodict['message']                     # >>> 'Integration successful.'
    # Extract transformations
    Gamma_trans_theta_u_temp, Gamma_trans_theta_v_temp = X_theta.T
    # Save our transformations
    Gamma_trans_theta_u_1.append(Gamma_trans_theta_u_temp)
    Gamma_trans_theta_v_1.append(Gamma_trans_theta_v_temp)    
# Solve to get the transformed angle solution
X1_angle, infodict = integrate.odeint(dX_dt, array([Gamma_trans_theta_u_1[0][-1], Gamma_trans_theta_v_1[0][-1]]), t, args = (a_stable,b_stable,c_stable,d_stable),full_output=True)
u_angle_1, v_angle_1 = X1_angle.T
# Allocate our empty lists
Gamma_trans_theta_u_2 = []
Gamma_trans_theta_v_2 = []
# Loop over our indices and find the transformations
for lin_index in lin_indices:
    # Get our transformations
    X_theta, infodict = integrate.odeint(Gamma_theta_ODE, array([u_angle_1[lin_index],v_angle_1[lin_index]]), epsilon_vec, args = (a_stable,b_stable,c_stable,d_stable),full_output=True)
    infodict['message']                     # >>> 'Integration successful.'
    # Extract transformations
    Gamma_trans_theta_u_temp, Gamma_trans_theta_v_temp = X_theta.T
    # Save our transformations
    Gamma_trans_theta_u_2.append(Gamma_trans_theta_u_temp)
    Gamma_trans_theta_v_2.append(Gamma_trans_theta_v_temp)    
# Solve to get the transformed angle solution
X2_angle, infodict = integrate.odeint(dX_dt, array([Gamma_trans_theta_u_2[0][-1], Gamma_trans_theta_v_2[0][-1]]), t, args = (a_stable,b_stable,c_stable,d_stable),full_output=True)
u_angle_2, v_angle_2 = X2_angle.T
# ---------------------------------------------------------------------------
# Overall properties (a,b,c,d)=(1,0,0,-2) Saddle point
# ---------------------------------------------------------------------------
# Define the parameters for the stable point
a_saddle = 1
b_saddle = 0
c_saddle = 0
d_saddle = -2
# Solve the system for the saddle parameters
X2, infodict = integrate.odeint(dX_dt, X0, t, args = (a_saddle,b_saddle,c_saddle,d_saddle),full_output=True)
infodict['message']                     # >>> 'Integration successful.'
# Extract the original solution with the defined parameters and initial conditions
u_2, v_2 = X2.T
# ---------------------------------------------------------------------------
# Saddle point (a,b,c,d)=(1,0,0,-2) radial symmetry
# ---------------------------------------------------------------------------
# Plot the action of the symmetry
epsilon_vec = arange(0,epsilon_r-0.005,epsilon_r/50)
# Allocate our empty lists
Gamma_saddle_r_u_1 = []
Gamma_saddle_r_v_1 = []
# Loop over our indices and find the transformations
for lin_index in lin_indices:
    # Get our transformations
    X_r, infodict = integrate.odeint(Gamma_r_ODE, array([u_2[lin_index],v_2[lin_index]]), epsilon_vec, args = (a_saddle,b_saddle,c_saddle,d_saddle),full_output=True)
    infodict['message']                     # >>> 'Integration successful.'
    # Extract transformations
    Gamma_saddle_r_u_temp, Gamma_saddle_r_v_temp = X_r.T
    # Save our transformations
    Gamma_saddle_r_u_1.append(Gamma_saddle_r_u_temp)
    Gamma_saddle_r_v_1.append(Gamma_saddle_r_v_temp)    
# Solve to get the transformed angle solution
X1_r, infodict = integrate.odeint(dX_dt, array([Gamma_saddle_r_u_1[0][-1], Gamma_saddle_r_v_1[0][-1]]), t, args = (a_saddle,b_saddle,c_saddle,d_saddle),full_output=True)
u_saddle_r_1, v_saddle_r_1 = X1_r.T
# Allocate our empty lists
Gamma_saddle_r_u_2 = []
Gamma_saddle_r_v_2 = []
# Loop over our indices and find the transformations
for lin_index in lin_indices:
    # Get our transformations
    X_r, infodict = integrate.odeint(Gamma_r_ODE, array([u_saddle_r_1[lin_index],v_saddle_r_1[lin_index]]), epsilon_vec, args = (a_saddle,b_saddle,c_saddle,d_saddle),full_output=True)
    infodict['message']                     # >>> 'Integration successful.'
    # Extract transformations
    Gamma_trans_r_u_temp, Gamma_trans_r_v_temp = X_r.T
    # Save our transformations
    Gamma_saddle_r_u_2.append(Gamma_trans_r_u_temp)
    Gamma_saddle_r_v_2.append(Gamma_trans_r_v_temp)    
# Solve to get the transformed angle solution
X2_r, infodict = integrate.odeint(dX_dt, array([Gamma_saddle_r_u_2[0][-1], Gamma_saddle_r_v_2[0][-1]]), t, args = (a_saddle,b_saddle,c_saddle,d_saddle),full_output=True)
u_saddle_r_2, v_saddle_r_2 = X2_r.T
# ---------------------------------------------------------------------------
# Stable node (a,b,c,d)=(-1,0,0,-2) angle symmetry
# ---------------------------------------------------------------------------
# Plot the action of the symmetry
epsilon_vec = arange(0,epsilon_theta-0.005,epsilon_theta/50)
# Allocate our empty lists
Gamma_saddle_theta_u_1 = []
Gamma_saddle_theta_v_1 = []
# Loop over our indices and find the transformations
for lin_index in lin_indices:
    # Get our transformations
    X_theta, infodict = integrate.odeint(Gamma_theta_ODE, array([u_2[lin_index],v_2[lin_index]]), epsilon_vec, args = (a_saddle,b_saddle,c_saddle,d_saddle),full_output=True)
    infodict['message']                     # >>> 'Integration successful.'
    # Extract transformations
    Gamma_trans_theta_u_temp, Gamma_trans_theta_v_temp = X_theta.T
    # Save our transformations
    Gamma_saddle_theta_u_1.append(Gamma_trans_theta_u_temp)
    Gamma_saddle_theta_v_1.append(Gamma_trans_theta_v_temp)    
# Define a temporary time vector since we need to solve this ODE for a longer time
t_trans_1 = linspace(0, 2.75, 500)              # Time
# Solve to get the transformed angle solution
X1_angle, infodict = integrate.odeint(dX_dt, array([Gamma_saddle_theta_u_1[0][-1], Gamma_saddle_theta_v_1[0][-1]]), t_trans_1, args = (a_saddle,b_saddle,c_saddle,d_saddle),full_output=True)
# Extract the solution
u_saddle_angle_1, v_saddle_angle_1 = X1_angle.T
# Allocate our empty lists
Gamma_saddle_theta_u_2 = []
Gamma_saddle_theta_v_2 = []
# Loop over our indices and find the transformations
for lin_index in lin_indices:
    # Get our transformations
    X_theta, infodict = integrate.odeint(Gamma_theta_ODE, array([u_saddle_angle_1[lin_index],v_saddle_angle_1[lin_index]]), epsilon_vec, args = (a_saddle,b_saddle,c_saddle,d_saddle),full_output=True)
    infodict['message']                     # >>> 'Integration successful.'
    # Extract transformations
    Gamma_trans_theta_u_temp, Gamma_trans_theta_v_temp = X_theta.T
    # Save our transformations
    Gamma_saddle_theta_u_2.append(Gamma_trans_theta_u_temp)
    Gamma_saddle_theta_v_2.append(Gamma_trans_theta_v_temp)    
# Define a temporary time vector since we need to solve this ODE for a longer time
t_trans_2 = linspace(0, 3.5, 500)              # Time
# Solve to get the transformed angle solution
X2_angle, infodict = integrate.odeint(dX_dt, array([Gamma_saddle_theta_u_2[0][-1], Gamma_saddle_theta_v_2[0][-1]]), t_trans_2, args = (a_saddle,b_saddle,c_saddle,d_saddle),full_output=True)
# Extract the solution
u_saddle_angle_2, v_saddle_angle_2 = X2_angle.T
#=================================================================================
#=================================================================================
# Plot the linear system
#=================================================================================
#=================================================================================
#Define the first figure
f1, axs = p.subplots(2, 2, constrained_layout=True, figsize=(20, 8))
# Plot 1: Stable node with radial symmetry
axs[0][0].plot(u_1,v_1,color=(0/256,68/256,27/256),label="$(u,v)$",linewidth=4.0)
axs[0][0].plot(u_r_1,v_r_1,color=(0/256,109/256,44/256),label="$(\\hat{u},\\hat{v})$",linewidth=4.0)
axs[0][0].plot(u_r_2,v_r_2,color=(102/256,194/256,164/256),label="$(\\hat{\\hat{u}},\\hat{\\hat{v}})$",linewidth=4.0)
for index in range(len(Gamma_trans_r_u_1)):
    if index == 0:
        axs[0][0].plot(Gamma_trans_r_u_1[index], Gamma_trans_r_v_1[index],'--',color=(0,0,0),label="$\\Gamma_{\\epsilon}^{r}$",linewidth=0.5)
    else:
        axs[0][0].plot(Gamma_trans_r_u_1[index], Gamma_trans_r_v_1[index],'--',color=(0,0,0),linewidth=0.5)
for index in range(len(Gamma_trans_r_u_2)):
    axs[0][0].plot(Gamma_trans_r_u_2[index], Gamma_trans_r_v_2[index],'--',color=(0,0,0),linewidth=0.5)
axs[0][0].grid()
axs[0][0].legend(loc='best',prop={"size":20})
axs[0][0].set_xlabel(xlabel="$u(\\tau)$",fontsize=25)
axs[0][0].set_ylabel(ylabel="$v(\\tau)$",fontsize=25)
axs[0][0].tick_params(axis='both', which='major', labelsize=20)
axs[0][0].tick_params(axis='both', which='minor', labelsize=20)
# Plot 2: Stable node with angular symmetry
axs[0][1].plot(u_1,v_1,color=(103/256,0/256,31/256),label="$(u,v)$",linewidth=4.0)
axs[0][1].plot(u_angle_1,v_angle_1,color=(206/256,18/256,86/256),label="$(\\hat{u},\\hat{v})$",linewidth=4.0)
axs[0][1].plot(u_angle_2,v_angle_2,color=(223/256,101/256,176/256),label="$(\\hat{\\hat{u}},\\hat{\\hat{v}})$",linewidth=4.0)
for index in range(len(Gamma_trans_theta_u_1)):
    if index == 0:
        axs[0][1].plot(Gamma_trans_theta_u_1[index], Gamma_trans_theta_v_1[index],'--',color=(0,0,0),label="$\\Gamma_{\\epsilon}^{\\theta}$",linewidth=0.5)
    else:
        axs[0][1].plot(Gamma_trans_theta_u_1[index], Gamma_trans_theta_v_1[index],'--',color=(0,0,0),linewidth=0.5)
for index in range(len(Gamma_trans_theta_u_2)):
    axs[0][1].plot(Gamma_trans_theta_u_2[index], Gamma_trans_theta_v_2[index],'--',color=(0,0,0),linewidth=0.5)
axs[0][1].grid()
axs[0][1].legend(loc='best',prop={"size":20})
axs[0][1].set_xlabel(xlabel="$u(\\tau)$",fontsize=25)
axs[0][1].set_ylabel(ylabel="$v(\\tau)$",fontsize=25)
axs[0][1].tick_params(axis='both', which='major', labelsize=20)
axs[0][1].tick_params(axis='both', which='minor', labelsize=20)
# Plot 3: Saddle point with radial symmetry
axs[1][0].plot(u_2,v_2,color=(2/256,56/256,88/256),label="$(u,v)$",linewidth=4.0)
axs[1][0].plot(u_saddle_r_1,v_saddle_r_1,color=(5/256,112/256,176/256),label="$(\\hat{u},\\hat{v})$",linewidth=4.0)
axs[1][0].plot(u_saddle_r_2,v_saddle_r_2,color=(116/256,169/256,207/256),label="$(\\hat{\\hat{u}},\\hat{\\hat{v}})$",linewidth=4.0)
for index in range(len(Gamma_saddle_r_u_1)):
    if index == 0:
        axs[1][0].plot(Gamma_saddle_r_u_1[index], Gamma_saddle_r_v_1[index],'--',color=(0,0,0),label="$\\Gamma_{\\epsilon}^{r}$",linewidth=0.5)
    else:
        axs[1][0].plot(Gamma_saddle_r_u_1[index], Gamma_saddle_r_v_1[index],'--',color=(0,0,0),linewidth=0.5)
for index in range(len(Gamma_saddle_r_u_2)):
    axs[1][0].plot(Gamma_saddle_r_u_2[index], Gamma_saddle_r_v_2[index],'--',color=(0,0,0),linewidth=0.5)
axs[1][0].grid()
axs[1][0].legend(loc='best',prop={"size":20})
axs[1][0].set_xlabel(xlabel="$u(\\tau)$",fontsize=25)
axs[1][0].set_ylabel(ylabel="$v(\\tau)$",fontsize=25)
axs[1][0].tick_params(axis='both', which='major', labelsize=20)
axs[1][0].tick_params(axis='both', which='minor', labelsize=20)
# Plot 4: Saddle point with angular symmetry
axs[1][1].plot(u_2,v_2,color=(102/256,37/256,6/256),label="$(u,v)$",linewidth=4.0)
axs[1][1].plot(u_saddle_angle_1,v_saddle_angle_1,color=(204/256,76/256,2/256),label="$(\\hat{u},\\hat{v})$",linewidth=4.0)
axs[1][1].plot(u_saddle_angle_2,v_saddle_angle_2,color=(254/256,153/256,41/256),label="$(\\hat{\\hat{u}},\\hat{\\hat{v}})$",linewidth=4.0)
for index in range(len(Gamma_saddle_theta_u_1)):
    if index == 0:
        axs[1][1].plot(Gamma_saddle_theta_u_1[index], Gamma_saddle_theta_v_1[index],'--',color=(0,0,0),label="$\\Gamma_{\\epsilon}^{\\theta}$",linewidth=0.5)
    else:
        axs[1][1].plot(Gamma_saddle_theta_u_1[index], Gamma_saddle_theta_v_1[index],'--',color=(0,0,0),linewidth=0.5)
for index in range(len(Gamma_saddle_theta_u_2)):
    axs[1][1].plot(Gamma_saddle_theta_u_2[index], Gamma_saddle_theta_v_2[index],'--',color=(0,0,0),linewidth=0.5)
axs[1][1].grid()
axs[1][1].legend(loc='best',prop={"size":20})
axs[1][1].set_xlabel(xlabel="$u(\\tau)$",fontsize=25)
axs[1][1].set_ylabel(ylabel="$v(\\tau)$",fontsize=25)
axs[1][1].tick_params(axis='both', which='major', labelsize=20)
axs[1][1].tick_params(axis='both', which='minor', labelsize=20)
# Show the plot in the end
#plt.show()


# Plot in LaTeX as well...
# STABLE POINT RADIAL SYMMETRY
plot_LaTeX_2D(u_1,v_1,"../Figures/linear_symmetries/Input/stable_r.tex","color=stab_r_1,line width=1.5pt,","$(u,v)$")
plot_LaTeX_2D(u_r_1,v_r_1,"../Figures/linear_symmetries/Input/stable_r.tex","color=stab_r_2,line width=1.5pt,","$(\\hat{u},\\hat{v})$")
plot_LaTeX_2D(u_r_2,v_r_2,"../Figures/linear_symmetries/Input/stable_r.tex","color=stab_r_3,line width=1.5pt,","$(\\hat{\\hat{u}},\\hat{\\hat{v}})$")
for index in range(0,round(len(Gamma_trans_r_u_1)/2.5)):
    if index == 0:
        plot_LaTeX_2D(Gamma_trans_r_u_1[index][0:-2], Gamma_trans_r_v_1[index][0:-2],"../Figures/linear_symmetries/Input/stable_r.tex","color=black,->,>=latex,densely dashed,line width=1.0pt","$\\Gamma_{\\epsilon}^{r}$")
    else:
        plot_LaTeX_2D(Gamma_trans_r_u_1[index][0:-2], Gamma_trans_r_v_1[index][0:-2],"../Figures/linear_symmetries/Input/stable_r.tex","color=black,->,>=latex,densely dashed,line width=1.0pt",[])        
for index in range(0,round(len(Gamma_trans_r_u_2)/2.5)):
    plot_LaTeX_2D(Gamma_trans_r_u_2[index][0:-2], Gamma_trans_r_v_2[index][0:-2],"../Figures/linear_symmetries/Input/stable_r.tex","color=black,->,>=latex,densely dashed,line width=1.0pt",[])            
# STABLE POINT ANGULAR SYMMETRY
plot_LaTeX_2D(u_1,v_1,"../Figures/linear_symmetries/Input/stable_theta.tex","color=stab_theta_1,line width=1.5pt,","$(u,v)$")
plot_LaTeX_2D(u_angle_1,v_angle_1,"../Figures/linear_symmetries/Input/stable_theta.tex","color=stab_theta_2,line width=1.5pt,","$(\\hat{u},\\hat{v})$")
plot_LaTeX_2D(u_angle_2,v_angle_2,"../Figures/linear_symmetries/Input/stable_theta.tex","color=stab_theta_3,line width=1.5pt,","$(\\hat{\\hat{u}},\\hat{\\hat{v}})$")
for index in range(0,round(len(Gamma_trans_theta_u_1)/2.5)):
    if index == 0:
        plot_LaTeX_2D(Gamma_trans_theta_u_1[index][0:-1], Gamma_trans_theta_v_1[index][0:-1],"../Figures/linear_symmetries/Input/stable_theta.tex","color=black,->,>=latex,densely dashed,line width=1.0pt","$\\Gamma_{\\epsilon}^{\\theta}$")
    else:
        plot_LaTeX_2D(Gamma_trans_theta_u_1[index][0:-1], Gamma_trans_theta_v_1[index][0:-1],"../Figures/linear_symmetries/Input/stable_theta.tex","color=black,->,>=latex,densely dashed,line width=1.0pt",[])        
for index in range(0,round(len(Gamma_trans_theta_u_2)/2.5)):
    plot_LaTeX_2D(Gamma_trans_theta_u_2[index][0:-1], Gamma_trans_theta_v_2[index][0:-1],"../Figures/linear_symmetries/Input/stable_theta.tex","color=black,->,>=latex,densely dashed,line width=1.0pt",[])            
# SADDLE POINT RADIAL SYMMETRY
plot_LaTeX_2D(u_2,v_2,"../Figures/linear_symmetries/Input/saddle_r.tex","color=saddle_r_1,line width=1.5pt,","$(u,v)$")
plot_LaTeX_2D(u_saddle_r_1,v_saddle_r_1,"../Figures/linear_symmetries/Input/saddle_r.tex","color=saddle_r_2,line width=1.5pt,","$(\\hat{u},\\hat{v})$")
plot_LaTeX_2D(u_saddle_r_2,v_saddle_r_2,"../Figures/linear_symmetries/Input/saddle_r.tex","color=saddle_r_3,line width=1.5pt,","$(\\hat{\\hat{u}},\\hat{\\hat{v}})$")
for index in range(0,round(len(Gamma_saddle_r_u_1)/2.5)):
    if index == 0:
        plot_LaTeX_2D(Gamma_saddle_r_u_1[index][0:-2], Gamma_saddle_r_v_1[index][0:-2],"../Figures/linear_symmetries/Input/saddle_r.tex","color=black,->,>=latex,densely dashed,line width=1.0pt","$\\Gamma_{\\epsilon}^{r}$")
    else:
        plot_LaTeX_2D(Gamma_saddle_r_u_1[index][0:-2], Gamma_saddle_r_v_1[index][0:-2],"../Figures/linear_symmetries/Input/saddle_r.tex","color=black,->,>=latex,densely dashed,line width=1.0pt",[])        
for index in range(0,round(len(Gamma_saddle_r_u_2)/2.5)):
    plot_LaTeX_2D(Gamma_saddle_r_u_2[index][0:-2], Gamma_saddle_r_v_2[index][0:-2],"../Figures/linear_symmetries/Input/saddle_r.tex","color=black,->,>=latex,densely dashed,line width=1.0pt",[])
# SADDLE POINT ANGULAR SYMMETRY
plot_LaTeX_2D(u_2,v_2,"../Figures/linear_symmetries/Input/saddle_theta.tex","color=saddle_theta_1,line width=1.5pt,","$(u,v)$")
plot_LaTeX_2D(u_saddle_angle_1,v_saddle_angle_1,"../Figures/linear_symmetries/Input/saddle_theta.tex","color=saddle_theta_2,line width=1.5pt,","$(\\hat{u},\\hat{v})$")
plot_LaTeX_2D(u_saddle_angle_2,v_saddle_angle_2,"../Figures/linear_symmetries/Input/saddle_theta.tex","color=saddle_theta_3,line width=1.5pt,","$(\\hat{\\hat{u}},\\hat{\\hat{v}})$")
for index in range(1,round(len(Gamma_saddle_theta_u_1)*0.7)):
    if index == 1:
        plot_LaTeX_2D(Gamma_saddle_theta_u_1[index][0:-2], Gamma_saddle_theta_v_1[index][0:-2],"../Figures/linear_symmetries/Input/saddle_theta.tex","color=black,->,>=latex,densely dashed,line width=1.0pt","$\\Gamma_{\\epsilon}^{\\theta}$")
    else:
        plot_LaTeX_2D(Gamma_saddle_theta_u_1[index][0:-2], Gamma_saddle_theta_v_1[index][0:-2],"../Figures/linear_symmetries/Input/saddle_theta.tex","color=black,->,>=latex,densely dashed,line width=1.0pt",[])        
for index in range(1,round(0.55*len(Gamma_saddle_theta_u_2))):
    plot_LaTeX_2D(Gamma_saddle_theta_u_2[index][0:-2], Gamma_saddle_theta_v_2[index][0:-2],"../Figures/linear_symmetries/Input/saddle_theta.tex","color=black,->,>=latex,densely dashed,line width=1.0pt",[])


#=================================================================================
#=================================================================================
# Transforming the biological oscillator
#=================================================================================
#=================================================================================
# ---------------------------------------------------------------------------
# Overall properties
# ---------------------------------------------------------------------------
# The indices we wish to transform
lin_indices = arange(0,60,6)
# Plot the action of the symmetry
epsilon_vec = arange(0,epsilon_r-0.005,epsilon_r/50)
# The frequency rate
omega = 1
# Define the time vector
t = linspace(0, 10, 500)              # Time
X0 = array([2, 2])                  # ICs
X1_os, infodict = integrate.odeint(dX_dt_os, X0, t, args = ((omega),),full_output=True)
# Extract the original solution with the defined parameters and initial conditions
u_3, v_3 = X1_os.T
# Allocate our empty lists
Gamma_os_r_u_1 = []
Gamma_os_r_v_1 = []
# Loop over our indices and find the transformations
for lin_index in lin_indices:
    # Get our transformations
    X_r, infodict = integrate.odeint(Gamma_r_os, array([u_3[lin_index],v_3[lin_index]]), epsilon_vec, args = ((omega),),full_output=True)
    infodict['message']                     # >>> 'Integration successful.'
    # Extract transformations
    Gamma_os_r_u_temp, Gamma_os_r_v_temp = X_r.T
    # Save our transformations
    Gamma_os_r_u_1.append(Gamma_os_r_u_temp)
    Gamma_os_r_v_1.append(Gamma_os_r_v_temp)
# Solve to get the transformed angle solution
X1_r, infodict = integrate.odeint(dX_dt_os, array([Gamma_os_r_u_1[0][-1], Gamma_os_r_v_1[0][-1]]), t, args = ((omega,)),full_output=True)    
# Extract the solution
u_os_r_2, v_os_r_2 = X1_r.T
# Plot the action of the symmetry
epsilon_vec = arange(0,epsilon_theta-0.005,epsilon_theta/50)
# Allocate our empty lists
Gamma_os_theta_u_1 = []
Gamma_os_theta_v_1 = []
# Loop over our indices and find the transformations
for lin_index in lin_indices:
    # Get our transformations
    X_theta, infodict = integrate.odeint(Gamma_theta_os, array([u_3[lin_index],v_3[lin_index]]), epsilon_vec, args = ((omega),),full_output=True)
    infodict['message']                     # >>> 'Integration successful.'
    # Extract transformations
    Gamma_os_theta_u_temp, Gamma_os_theta_v_temp = X_theta.T
    # Save our transformations
    Gamma_os_theta_u_1.append(Gamma_os_theta_u_temp)
    Gamma_os_theta_v_1.append(Gamma_os_theta_v_temp)
# Solve to get the transformed angle solution
X1_theta, infodict = integrate.odeint(dX_dt_os, array([Gamma_os_theta_u_1[0][-1], Gamma_os_theta_v_1[0][-1]]), t, args = ((omega,)),full_output=True)    
# Extract the solution
u_os_theta_2, v_os_theta_2 = X1_theta.T
#=================================================================================
#=================================================================================
# Plot the oscillatory system
#=================================================================================
#=================================================================================
#Define the first figure
f2, axs = p.subplots(1, 2, constrained_layout=True, figsize=(20, 8))
# Plot 1: Radial symmetry on the biological oscillator
axs[0].plot(u_3,v_3,color=(0/256,68/256,27/256),label="$(u,v)$",linewidth=4.0)
axs[0].plot(u_os_r_2,v_os_r_2,color=(0/256,109/256,44/256),label="$(\\hat{u},\\hat{v})$",linewidth=4.0)
for index in range(len(Gamma_os_r_u_1)):
    if index == 0:
        axs[0].plot(Gamma_os_r_u_1[index], Gamma_os_r_v_1[index],'--',color=(0,0,0),label="$\\Gamma_{\\epsilon}^{r}$",linewidth=0.5)
    else:
        axs[0].plot(Gamma_os_r_u_1[index], Gamma_os_r_v_1[index],'--',color=(0,0,0),linewidth=0.5)
axs[0].grid()
axs[0].legend(loc='best',prop={"size":20})
axs[0].set_xlabel(xlabel="$u(\\tau)$",fontsize=25)
axs[0].set_ylabel(ylabel="$v(\\tau)$",fontsize=25)
axs[0].tick_params(axis='both', which='major', labelsize=20)
axs[0].tick_params(axis='both', which='minor', labelsize=20)
# Plot 2: Angular symmetry on the biological oscillator
axs[1].plot(u_3,v_3,color=(103/256,0/256,31/256),label="$(u,v)$",linewidth=4.0)
axs[1].plot(u_os_theta_2,v_os_theta_2,color=(206/256,18/256,86/256),label="$(\\hat{u},\\hat{v})$",linewidth=4.0)
for index in range(len(Gamma_os_theta_u_1)):
    if index == 0:
        axs[1].plot(Gamma_os_theta_u_1[index], Gamma_os_theta_v_1[index],'--',color=(0,0,0),label="$\\Gamma_{\\epsilon}^{\\theta}$",linewidth=0.5)
    else:
        axs[1].plot(Gamma_os_theta_u_1[index], Gamma_os_theta_v_1[index],'--',color=(0,0,0),linewidth=0.5)
axs[1].grid()
axs[1].legend(loc='best',prop={"size":20})
axs[1].set_xlabel(xlabel="$u(\\tau)$",fontsize=25)
axs[1].set_ylabel(ylabel="$v(\\tau)$",fontsize=25)
axs[1].tick_params(axis='both', which='major', labelsize=20)
axs[1].tick_params(axis='both', which='minor', labelsize=20)
# Show the plot in the end
plt.show()

# Plot in LaTeX as well...
# RADIAL SYMMETRY
for index in range(0,len(Gamma_os_r_u_1)):
    if index == 0:
        plot_LaTeX_2D(Gamma_os_r_u_1[index][0:-1], Gamma_os_r_v_1[index][0:-1],"../Figures/oscillator_symmetries/Input/r.tex","color=black,->,>=latex,densely dashed,line width=1.0pt","$\\Gamma_{\\epsilon}^{r}$")
    else:
        plot_LaTeX_2D(Gamma_os_r_u_1[index][0:-1], Gamma_os_r_v_1[index][0:-1],"../Figures/oscillator_symmetries/Input/r.tex","color=black,->,>=latex,densely dashed,line width=1.0pt",[])
plot_LaTeX_2D(u_3,v_3,"../Figures/oscillator_symmetries/Input/r.tex","color=r_1,line width=1.5pt,","$(u,v)$")
plot_LaTeX_2D(u_os_r_2,v_os_r_2,"../Figures/oscillator_symmetries/Input/r.tex","color=r_2,line width=1.5pt,","$(\\hat{u},\\hat{v})$")

# ANGULAR SYMMETRY
for index in range(0,len(Gamma_os_theta_u_1)):
    if index == 0:
        plot_LaTeX_2D(Gamma_os_theta_u_1[index][0:-1], Gamma_os_theta_v_1[index][0:-1],"../Figures/oscillator_symmetries/Input/theta.tex","color=black,->,>=latex,densely dashed,line width=1.0pt","$\\Gamma_{\\epsilon}^{\\theta}$")
    else:
        plot_LaTeX_2D(Gamma_os_theta_u_1[index][0:-1], Gamma_os_theta_v_1[index][0:-1],"../Figures/oscillator_symmetries/Input/theta.tex","color=black,->,>=latex,densely dashed,line width=1.0pt",[])
plot_LaTeX_2D(u_3,v_3,"../Figures/oscillator_symmetries/Input/theta.tex","color=theta_1,line width=1.5pt,","$(u,v)$")
plot_LaTeX_2D(u_os_theta_2,v_os_theta_2,"../Figures/oscillator_symmetries/Input/theta.tex","color=theta_2,line width=1.5pt,","$(\\hat{u},\\hat{v})$")

