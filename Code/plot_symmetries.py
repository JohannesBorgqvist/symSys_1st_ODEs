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
def dX_dt(X, t=0,a=1):
    """ Return the growth rate of fox and rabbit populations. """
    return array([ X[0]*(1-X[1]) ,
                   a*X[1]*(X[0]-1)])
# Function 2: u_transf
def u_transf(u, epsilon):
    func = lambda u_hat :  u - log(u) + epsilon - (u_hat - log(u_hat))
    if abs(u-1)<0.000000000001:
        u_hat_solution = u
    else:
        u_hat_solution = fsolve(func, u)[0]
    return u_hat_solution
# Function 3: v_transf
def v_transf(v, epsilon, alpha):
    func = lambda v_hat :  log(v**(1/alpha)) - (v/alpha) + epsilon - (log(v_hat**(1/alpha)) - (v_hat/alpha))
    if abs(v-1)<0.00000000001:
        v_hat_solution = v
    else:
        v_hat_solution = fsolve(func, v)[0]
    return v_hat_solution
# Function 4: S_transf
def S_transf(S, epsilon, a, r):
    #func = lambda S_hat :  S - (1/r)*log(a-r*S)+epsilon-S_hat+(1/r)*log(a-r*S_hat)
    func = lambda S_hat :  S - (a/r)*log(S)+epsilon-(S_hat - (a/r)*log(S_hat))
    S_hat_solution = fsolve(func, S)[0]
    return S_hat_solution
# Function 5: IC
def transform_IC_LV(u0, v0, alpha, H):
    func = lambda v:  v+(alpha*u0)-log((u0**alpha)*v) - H
    v0 = fsolve(func, v0)[0]
    return v0
# Function 6: plot_LaTeX_2D
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
def plot_LaTeX_3D(data,file_str,plot_str,legend_str,surfaceNotCurve):
    # Open a file with the append option
    # so that we can write to the same
    # file multiple times
    f = open(file_str, "a")
    # Create a temporary string which
    # is the one that does the plotting.
    # Here we incorporate the input plot_str
    # which contains the color, and the markers
    # of the plot at hand
    if surfaceNotCurve:
        if len(legend_str)==0:
            temp_str = "\\addplot3[forget plot," + plot_str+ "]\n"
        else:
            temp_str = "\\addplot3[" + plot_str+ "]\n"
    else:
        if len(legend_str)==0:
            temp_str = "\\addplot3+[forget plot," + plot_str+ "]\n"
        else:
            temp_str = "\\addplot3+[" + plot_str+ "]\n"        
    # Add the coordinates
    temp_str += "coordinates {%\n"
    # Loop over the input files and add
    # them to the file
    for index in range(len(data)):
        temp_str += "("+str(data[index][0]) + "," + str(data[index][1]) + "," + str(data[index][2]) + ")"
        if index>0:
            if index < len(data)-1 and data[index][1] < data[index+1][1]:
                temp_str += "\n"
            elif index == len(data)-1:
                temp_str += "\n"
            else:
                temp_str += "  "
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
# Solving the ODEs for the LV model
#=================================================================================
#=================================================================================
t = linspace(0, 10, 500)              # time
X0 = array([1, 0.10])                     # initials conditions: 10 rabbits and 5 foxes
#X0_2 = array([20, 10])                     # initials conditions: 20 rabbits and 10 foxes
X1, infodict = integrate.odeint(dX_dt, X0, t, args = (1,),full_output=True)
infodict['message']                     # >>> 'Integration successful.'

#!python
#X2, infodict = integrate.odeint(dX_dt, X0, t, args = (2,),full_output=True)
#infodict['message']                     # >>> 'Integration successful.'

#!python
u, v = X1.T
epsilon=0.5
alpha = 1
# Testing to find ICs
# Calculate internal energy
H= v[0] + alpha*u[0] - log((u[0]**alpha)*v[0])
# Define the u-coordinate in the new initial condition
u0 = X0[0]
v0 = X0[0]
# Find the other v-coordinate in the transformed initial condition
v0_new = transform_IC_LV(u0, v0, alpha, H-alpha*epsilon)
v0_new_2 = transform_IC_LV(u0, v0, alpha, H+alpha*epsilon)
# Print our new initial conditions
print("The original IC:\t(u0,v0)\t=\t(%0.3f,%0.3f)"%(u0,v0))
print("The original energy:\tH\t=\t%0.3f"%(H))
print("The transformed IC:\t(u0,v0)\t=\t(%0.3f,%0.3f)"%(u0,v0_new))
print("The transformed energy:\tH\t=\t%0.3f"%(H-alpha*epsilon))
# Try to solve the LV-model again with these new ICs
X2, infodict = integrate.odeint(dX_dt, array([u0,v0_new]), t, args = (1,),full_output=True)
u_transformed, v_transformed = X2.T
X3, infodict = integrate.odeint(dX_dt, array([u0,v0_new_2]), t, args = (1,),full_output=True)
u_transformed_2, v_transformed_2 = X3.T
# Transformed solutions
t_trans = asarray([t_temp + epsilon for t_temp in list(t)])
#u_transformed = asarray([u_transf(u_temp, epsilon) for u_temp in list(u)])
#v_transformed = asarray([v_transf(v_temp, epsilon,alpha) for v_temp in list(v)])
# Transformations time translation
t_sym = []
epsilon_vec = arange(0,epsilon-0.005,epsilon/50)
t_indices = arange(60,78,3)
for t_index in list(t_indices):
    trans_vec = [t[t_index]+epsilon_temp for epsilon_temp in list(epsilon_vec)]
    t_sym.append(trans_vec)
# Transformations u 
u_sym = []
u_indices = arange(70,125,3)
for u_index in list(u_indices):
    trans_vec = [u_transf(u[u_index],epsilon_temp) for epsilon_temp in list(epsilon_vec)]
    u_sym.append(trans_vec)
# Transformations v
v_sym = []
v_indices = arange(107,150,2)
for v_index in list(v_indices):
    trans_vec = [v_transf(v[v_index],epsilon_temp,alpha) for epsilon_temp in list(epsilon_vec)]
    v_sym.append(trans_vec)
# Update so that we remove all non nan number
#t = asarray([t[u_index] for u_index,u_temp in enumerate(u) if not math.isnan(u_temp)])
#u = asarray([u_temp for u_index,u_temp in enumerate(u) if not math.isnan(u_temp)])
#v = asarray([v[u_index] for u_index,u_temp in enumerate(u) if not math.isnan(u_temp)])
#t = asarray([t[v_index] for v_index,v_temp in enumerate(v) if not math.isnan(v_temp)])
#u = asarray([u[v_index] for v_index,v_temp in enumerate(v) if not math.isnan(v_temp)])
#v = asarray([v_temp for v_index,v_temp in enumerate(v) if not math.isnan(v_temp)])
#=================================================================================
#=================================================================================
# Plotting the solutions
#=================================================================================
#=================================================================================
#Define the first figure
f1, axs = p.subplots(1, 3, constrained_layout=True, figsize=(20, 8))
axs[0].plot(asarray(t_sym[0]),asarray([u[t_indices[0]]*((index+1)/(index+1)) for index in range(len(epsilon_vec))]), '--', label="$\\left.\\Gamma_{\\epsilon}^{\\tau}\\right|_{\\epsilon=" +str(epsilon) + "}$" ,color=(0,0,0),linewidth=2.0)
for index,t_index in enumerate(list(t_indices)):
    axs[0].plot(asarray(t_sym[index]),asarray([u[t_index]*((index+1)/(index+1)) for index in range(len(epsilon_vec))]), '--' ,color=(0,0,0),linewidth=2.0)
    axs[0].plot(asarray(t_sym[index]),asarray([v[t_index]*((index+1)/(index+1)) for index in range(len(epsilon_vec))]), '--' ,color=(0,0,0),linewidth=2.0)    
axs[0].plot(t, u, '-', label="$u(\\tau)$" ,color=(0/256,68/256,27/256),linewidth=3.0)
axs[0].plot(t, v  , '-', label='$v(\\tau)$',color=(153/256,216/256,201/256),linewidth=3.0)
axs[0].plot(t_trans, u, '-', label='$u(\\hat{\\tau})$',color=(77/256,0/256,75/256),linewidth=3.0)
axs[0].plot(t_trans, v  , '-', label='$v(\\hat{\\tau})$',color=(140/256,150/256,198/256),linewidth=3.0)
axs[0].grid()
axs[0].legend(loc='best',prop={"size":20})
axs[1].plot(u, v,'-', label='Original population, $(u,v)$',color=(0/256,68/256,27/256),linewidth=3.0)
axs[1].plot(u_transformed_2,v_transformed_2, '-', label='Transformed population, $(\\hat{u},v)$',color=(77/256,0/256,75/256),linewidth=3.0)
axs[1].plot(asarray(u_sym[0]),asarray([v[u_indices[0]]*index/index for index in range(len(epsilon_vec))]), '--', label="$\\left.\\Gamma_{\\epsilon}^{u}\\right|_{\\epsilon=" +str(epsilon) + "}$" ,color=(0,0,0),linewidth=2.0)
for index,u_index in enumerate(list(u_indices)):
    axs[1].plot(asarray(u_sym[index]),asarray([v[u_index]*((index+1)/(index+1)) for index in range(len(epsilon_vec))]), '--' ,color=(0,0,0),linewidth=2.0)    
axs[1].grid()
axs[1].legend(loc='best',prop={"size":20})
axs[2].plot(u, v,'-', label='Original population, $(u,v)$',color=(0/256,68/256,27/256),linewidth=3.0)
axs[2].plot(u_transformed,v_transformed, '-', label='Transformed population, $(u,\\hat{v})$',color=(77/256,0/256,75/256),linewidth=3.0)
axs[2].plot(asarray([u[v_indices[0]]*((index+1)/(index+1)) for index in range(len(epsilon_vec))]), asarray(v_sym[0]), '--', label="$\\left.\\Gamma_{\\epsilon}^{v}\\right|_{\\epsilon=" +str(epsilon) + "}$" ,color=(0,0,0),linewidth=2.0)
for index,v_index in enumerate(list(v_indices)):
    axs[2].plot(asarray([u[v_index]*((index+1)/(index+1)) for index in range(len(epsilon_vec))]),asarray(v_sym[index]), '--' ,color=(0,0,0),linewidth=2.0)    
axs[2].grid()
axs[2].legend(loc='best',prop={"size":20})
# Set fontsize of labels
axs[0].set_xlabel(xlabel='Time, $\\tau$',fontsize=25)
axs[0].set_ylabel(ylabel='Population size',fontsize=25)
axs[1].set_xlabel(xlabel='Rabbits, $u(\\tau)$',fontsize=25)
axs[1].set_ylabel(ylabel='Foxes, $v(\\tau)$',fontsize=25)
axs[2].set_xlabel(xlabel='Rabbits, $u(\\tau)$',fontsize=25)
axs[2].set_ylabel(ylabel='Foxes, $v(\\tau)$',fontsize=25)
# Change the size of the ticks
axs[0].tick_params(axis='both', which='major', labelsize=20)
axs[0].tick_params(axis='both', which='minor', labelsize=20)
axs[1].tick_params(axis='both', which='major', labelsize=20)
axs[1].tick_params(axis='both', which='minor', labelsize=20)
axs[2].tick_params(axis='both', which='major', labelsize=20)
axs[2].tick_params(axis='both', which='minor', labelsize=20)
# Title and saving the figure
f1.suptitle('Symmetries of the Lotka-Volterra model',fontsize=30,weight='bold')
#f1.savefig('../Figures/LV_symmetries.png')
#p.show()


#=================================================================================
#=================================================================================
# Plotting the solutions in LaTeX
#=================================================================================
#=================================================================================
# PLOT THE ACTION OF THE SYMMETRY
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Translations in tau
for index,t_index in enumerate(list(t_indices)):
    if index == 0:
        plot_LaTeX_2D(asarray(t_sym[index]),asarray([u[t_index]*((index+1)/(index+1)) for index in range(len(epsilon_vec))]),"../Figures/LV_symmetries/Input/tau_trans.tex","color=black,->,>=latex,densely dashed","$\\Gamma^{\mathrm{LV},\\tau}_{\\epsilon}$")
        plot_LaTeX_2D(asarray(t_sym[index]),asarray([v[t_index]*((index+1)/(index+1)) for index in range(len(epsilon_vec))]),"../Figures/LV_symmetries/Input/tau_trans.tex","color=black,->,>=latex,densely dashed",[])         
    else:
        plot_LaTeX_2D(asarray(t_sym[index]),asarray([u[t_index]*((index+1)/(index+1)) for index in range(len(epsilon_vec))]),"../Figures/LV_symmetries/Input/tau_trans.tex","color=black,->,>=latex,densely dashed",[])
        plot_LaTeX_2D(asarray(t_sym[index]),asarray([v[t_index]*((index+1)/(index+1)) for index in range(len(epsilon_vec))]),"../Figures/LV_symmetries/Input/tau_trans.tex","color=black,->,>=latex,densely dashed",[]) 
plot_LaTeX_2D(t,u,"../Figures/LV_symmetries/Input/tau_trans.tex","color=clr_1,line width=1.5pt,","$u(\\tau)$")
plot_LaTeX_2D(t,v,"../Figures/LV_symmetries/Input/tau_trans.tex","color=clr_3,line width=1.5pt,","$v(\\tau)$")
plot_LaTeX_2D(t_trans,u,"../Figures/LV_symmetries/Input/tau_trans.tex","color=clr_2,line width=1.5pt,","$\\hat{u}(\\tau)$")
plot_LaTeX_2D(t_trans,v,"../Figures/LV_symmetries/Input/tau_trans.tex","color=clr_4,line width=1.5pt,","$\\hat{v}(\\tau)$")
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Scaled translations in u
plot_LaTeX_2D(asarray([u_sym_temp for u_sym_temp in u_sym[0]]),asarray([v[u_indices[0]]*((index+1)/(index+1)) for index in range(len(epsilon_vec))]),"../Figures/LV_symmetries/Input/u_trans.tex","color=black,->,>=latex,densely dashed","$\\Gamma^{\mathrm{LV},u}_{\\epsilon}$")
for index,u_index in enumerate(list(u_indices)):
    plot_LaTeX_2D(asarray([u_sym_temp for u_sym_temp in u_sym[index]]),asarray([v[u_index]*((index+1)/(index+1)) for index in range(len(epsilon_vec))]),"../Figures/LV_symmetries/Input/u_trans.tex","color=black,->,>=latex,densely dashed",[])
plot_LaTeX_2D(u,v,"../Figures/LV_symmetries/Input/u_trans.tex","color=clr_1,line width=1.5pt,","$(u,v)$")
plot_LaTeX_2D(u_transformed_2,v_transformed_2,"../Figures/LV_symmetries/Input/u_trans.tex","color=clr_2,line width=1.5pt,","$(\\hat{u},\\hat{v})$")
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Scaled translations in v
#plot_LaTeX_2D(asarray([u[v_indices[0]]*index/index for index in range(len(epsilon_vec))]),asarray([v_sym_temp[0] for v_sym_temp in v_sym[v_indices[0]]]),"../Figures/LV_symmetries/Input/v_trans.tex","color=black,->,>=latex,densely dashed","$\\Gamma^{\mathrm{LV},v}_{\\epsilon}$")
for index,v_index in enumerate(list(v_indices)):
    if index==0:
        plot_LaTeX_2D(asarray([u[v_index]*index/index for index in range(len(epsilon_vec))]),asarray([v_sym_temp for v_index, v_sym_temp in enumerate(v_sym[index])]),"../Figures/LV_symmetries/Input/v_trans.tex","color=black,->,>=latex,densely dashed","$\\Gamma^{\mathrm{LV},v}_{\\epsilon}$")
    else:
        plot_LaTeX_2D(asarray([u[v_index]*index/index for index in range(len(epsilon_vec))]),asarray([v_sym_temp for v_index, v_sym_temp in enumerate(v_sym[index])]),"../Figures/LV_symmetries/Input/v_trans.tex","color=black,->,>=latex,densely dashed",[])
        
plot_LaTeX_2D(u,v,"../Figures/LV_symmetries/Input/v_trans.tex","color=clr_3,line width=1.5pt,","$(u,v)$")
plot_LaTeX_2D(u_transformed,v_transformed,"../Figures/LV_symmetries/Input/v_trans.tex","color=clr_4,line width=1.5pt,","$(\\hat{u},\\hat{v})$")
#plot_LaTeX_2D(asarray([u[index] for index, v_temp in enumerate(v_transformed) if index>280 and index<456]),asarray([v_temp for index, v_temp in enumerate(v_transformed) if index >280 and index<456]),"../Figures/LV_symmetries/Input/v_trans.tex","color=clr_4,line width=1.5pt,","Transformed population")
#plot_LaTeX_2D(asarray([u[index] for index, v_temp in enumerate(v_transformed) if index>105 and index<175]),asarray([v_temp for index, v_temp in enumerate(v_transformed) if index >105 and index<175]),"../Figures/LV_symmetries/Input/v_trans.tex","color=clr_4,line width=1.5pt,",[])



#=================================================================================
#=================================================================================
# Plotting the SIR model
#=================================================================================
#=================================================================================
# DEFINE PARAMETER VALUES
r = 0.00218 # In units per days
p = 202 # In units days (defined as p=a/r)
a = r*p # In units per individuals per days
#S0_ori = 300 # Initial condition for S in units individuals
N = 763 # Total population size in units individuals
S0 = N-1 # Initial condition for S in units individuals
I0 = N-S0 # Initial condition for I in units individuals
# CALCULATE PROPERTIES OF THE SIR MODEL
# Define our constant for the solution trajectory
C = I0+S0-p*log(S0)
# Define I_max
I_max = N-p+p*log(p/S0)
# Define a large epsilon
epsilon = 100
# Allocate memory for the S-vector
S = linspace(0,S0,500)
I = asarray([C-S_temp +p*log(S_temp) for S_temp in S])
I_trans = asarray([(C+epsilon)-S_temp +p*log(S_temp) for S_temp in S])
I_trans_2 = asarray([(C+2*epsilon)-S_temp +p*log(S_temp) for S_temp in S])
# Calculate the transformed solution
S_transformed = asarray([S_transf(S_temp, epsilon,a,r) for S_temp in list(S)])
# Modify S and I so that we only plot up until the threshold N
S_new = asarray([S_temp for index,S_temp in enumerate(S) if (S_temp+I[index])<=N and I[index]>0])
I_new = asarray([I[index] for index,S_temp in enumerate(S) if (S_temp+I[index])<=N and I[index]>0])
# The same goes for the transformed solution
S_trans_new = asarray([S_temp for index,S_temp in enumerate(S) if (S_temp+I_trans[index])<=N and I_trans[index]>0])
I_trans_new = asarray([I_trans[index] for index,S_temp in enumerate(S) if (S_temp+I_trans[index])<=N and I_trans[index]>0])
S_trans_new_2 = asarray([S_temp for index,S_temp in enumerate(S) if (S_temp+I_trans_2[index])<=N and I_trans_2[index]>0])
I_trans_new_2 = asarray([I_trans_2[index] for index,S_temp in enumerate(S) if (S_temp+I_trans_2[index])<=N and I_trans_2[index]>0])
#I_trans = asarray([I[index] for index,S_temp in enumerate(S_new) if (S_temp+I[index])<=N and I[index]>0])
# Plot the action of the symmetry
epsilon_vec = arange(0,epsilon-0.005,epsilon/50)
# Transformations S
S_sym = []
S_indices = arange(35,65,5)
for S_index in list(S_indices):
    trans_vec = [S_transf(S_new[S_index],epsilon_temp,a,r) for epsilon_temp in list(epsilon_vec)]
    S_sym.append(list(trans_vec))
# Add the other ones as well
S_sym_2 = []
S_indices_2 = arange(21,30,2)
for S_index in list(S_indices_2):
    trans_vec = [S_transf(S_trans_new[S_index],epsilon_temp,a,r) for epsilon_temp in list(epsilon_vec)]
    S_sym_2.append(list(trans_vec))    
# Plot the solutions
f2, ax_2 = plt.subplots(1,1,constrained_layout=True, figsize=(20, 8))
#ax_2.plot(S,I, '--', label="Original solution, $(S,I)$" ,color=(0/256,68/256,27/256),linewidth=3.0)
ax_2.plot(S_new,I_new, '-', label="Original solution, $(S,I)$" ,color=(0/256,68/256,27/256),linewidth=3.0)
ax_2.plot(S_trans_new,I_trans_new, '-', label="Transformed solution, $(\hat{S},I)$" ,color=(35/256,139/256,69/256),linewidth=3.0)
ax_2.plot(S_trans_new_2,I_trans_new_2, '-', label="Transformed solution, $(\hat{\hat{S}},I)$" ,color=(102/256,194/256,164/256),linewidth=3.0)
for index,S_index in enumerate(list(S_indices)):
    if index == 0:
        ax_2.plot(asarray(S_sym[index]),asarray([I_new[S_index]*((temp_index+1)/(temp_index+1)) for temp_index in range(len(epsilon_vec))]), '--' ,color=(0,0,0),label="$\\Gamma^{\\mathrm{SIR},S}_{\\epsilon}$",linewidth=2.0)
    else:
        ax_2.plot(asarray(S_sym[index]),asarray([I_new[S_index]*((temp_index+1)/(temp_index+1)) for temp_index in range(len(epsilon_vec))]), '--' ,color=(0,0,0),linewidth=2.0)
#for index,S_val in enumerate(S_sym_2):
for index in range(len(S_sym_2)):
    ax_2.plot(asarray(S_sym_2[index]),asarray([I_trans_new[S_indices_2[index]]*((temp_index+1)/(temp_index+1)) for temp_index in range(len(epsilon_vec))]), '--' ,color=(0,0,0),linewidth=2.0)    
ax_2.plot(asarray([N, 0]),asarray([0, N]),label="Total population size, $N=I+S$",color=(0,0,0),linewidth=3.0)
ax_2.plot(asarray([0, p]),asarray([I_max, I_max]), '--' ,color=(0,0,0),linewidth=3.0)
ax_2.plot(asarray([0, p]),asarray([I_max+epsilon, I_max+epsilon]), '--' ,color=(0,0,0),linewidth=3.0)
ax_2.plot(asarray([0, p]),asarray([I_max+2*epsilon, I_max+2*epsilon]), '--' ,color=(0,0,0),linewidth=3.0)
ax_2.plot(asarray([p, p]),asarray([0, I_max+2*epsilon]), '--' ,color=(0,0,0),linewidth=3.0)
ax_2.tick_params(axis='both', which='major', labelsize=20)
ax_2.tick_params(axis='both', which='minor', labelsize=20)
ax_2.grid()
ax_2.legend(loc='best',prop={"size":20})
ax_2.set_xlabel(xlabel='Susceptibles, $S(t)$',fontsize=25)
ax_2.set_ylabel(ylabel='Infected, $I(t)$',fontsize=25)
plt.show()


#=================================================================================
#=================================================================================
# Plotting the solutions in LaTeX
#=================================================================================
#=================================================================================
for index,S_index in enumerate(list(S_indices)):
    if index == 0:
        plot_LaTeX_2D(asarray(S_sym[index]),asarray([((I_new[S_index]*(temp_index+1))/(temp_index+1)) for temp_index in range(len(epsilon_vec))]),"../Figures/SIR_symmetry/Input/S_trans.tex","color=black,->,>=latex,densely dashed,line width=1.0pt","$\\Gamma^{\mathrm{SIR},S}_{\\epsilon}$")
    elif index<5:
        plot_LaTeX_2D(asarray(S_sym[index]),asarray([((I_new[S_index]*(temp_index+1))/(temp_index+1)) for temp_index in range(len(epsilon_vec))]),"../Figures/SIR_symmetry/Input/S_trans.tex","color=black,->,>=latex,densely dashed",[])
for index in range(len(S_sym_2)):
    plot_LaTeX_2D(asarray(S_sym_2[index]),asarray([I_trans_new[S_indices_2[index]]*((temp_index+1)/(temp_index+1)) for temp_index in range(len(epsilon_vec))]),"../Figures/SIR_symmetry/Input/S_trans.tex","color=black,->,>=latex,densely dashed,line width=1.0pt",[])
plot_LaTeX_2D([N*x for x in linspace(0,1,50)],[N*(1-x) for x in linspace(0,1,B50)],"../Figures/SIR_symmetry/Input/S_trans.tex","color=black,mark=diamond*,only marks,line width=0.75pt,","$N=I+S$")
plot_LaTeX_2D([0, p],[I_max,I_max],"../Figures/SIR_symmetry/Input/S_trans.tex","color=black,densely dashed,line width=1.0pt,",[])
plot_LaTeX_2D([0, p],[I_max+epsilon,I_max+epsilon],"../Figures/SIR_symmetry/Input/S_trans.tex","color=black,densely dashed,line width=1.0pt,",[])
plot_LaTeX_2D([0, p],[I_max+2*epsilon,I_max+2*epsilon],"../Figures/SIR_symmetry/Input/S_trans.tex","color=black,densely dashed,line width=1.0pt,",[])
plot_LaTeX_2D([p, p],[0,I_max+2*epsilon],"../Figures/SIR_symmetry/Input/S_trans.tex","color=black,densely dashed,line width=1.0pt,",[])
# We add the solutions last
plot_LaTeX_2D(S_new,I_new,"../Figures/SIR_symmetry/Input/S_trans.tex","color=clr_1,line width=2.0pt,","$I(S)$")
plot_LaTeX_2D(S_trans_new,I_trans_new,"../Figures/SIR_symmetry/Input/S_trans.tex","color=clr_2,line width=2.0pt,","$\\hat{I}(S;\\epsilon)$")
plot_LaTeX_2D(S_trans_new_2,I_trans_new_2,"../Figures/SIR_symmetry/Input/S_trans.tex","color=clr_3,line width=2.0pt,","$\\hat{\\hat{I}}(S;\\epsilon)=\\hat{I}(S;2\\epsilon)$")

print("\n\n\tHere are some nice properties:")
print("\t\tN\t=\t%d"%(N))
print("\t\tp\t=\t%d"%(p))
print("\t\tI_max\t=\t%d"%(I_max))
