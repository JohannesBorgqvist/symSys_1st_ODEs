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
    u_hat_solution = fsolve(func, u)
    return u_hat_solution
# Function 3: v_transf
def v_transf(v, epsilon, alpha):
    func = lambda v_hat :  log(v**(1/alpha)) - (v/alpha) + epsilon - (log(v_hat**(1/alpha)) - (v_hat/alpha))
    v_hat_solution = fsolve(func, v)
    return v_hat_solution
# Function 2: plot_LaTeX_2D
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
# Solving the ODEs
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
# Transformed solutions
t_trans = asarray([t_temp + epsilon for t_temp in list(t)])
u_transformed = asarray([u_transf(u_temp, epsilon) for u_temp in list(u)])
v_transformed = asarray([v_transf(v_temp, epsilon,alpha) for v_temp in list(v)])
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
t = asarray([t[u_index] for u_index,u_temp in enumerate(u) if not math.isnan(u_temp)])
u = asarray([u_temp for u_index,u_temp in enumerate(u) if not math.isnan(u_temp)])
v = asarray([v[u_index] for u_index,u_temp in enumerate(u) if not math.isnan(u_temp)])
t = asarray([t[v_index] for v_index,v_temp in enumerate(v) if not math.isnan(v_temp)])
u = asarray([u[v_index] for v_index,v_temp in enumerate(v) if not math.isnan(v_temp)])
v = asarray([v_temp for v_index,v_temp in enumerate(v) if not math.isnan(v_temp)])
#=================================================================================
#=================================================================================
# Plotting the solutions
#=================================================================================
#=================================================================================
#Define the first figure
f1, axs = p.subplots(1, 3, constrained_layout=True, figsize=(20, 8))
axs[0].plot(asarray(t_sym[0]),asarray([u[0]*index/index for index in range(len(epsilon_vec))]), '--', label="$\\left.\\Gamma_{\\epsilon}^{\\tau}\\right|_{\\epsilon=" +str(epsilon) + "}$" ,color=(0,0,0),linewidth=2.0)
for index,t_index in enumerate(list(t_indices)):
    axs[0].plot(asarray(t_sym[index]),asarray([u[t_index]*index/index for index in range(len(epsilon_vec))]), '--' ,color=(0,0,0),linewidth=2.0)
    axs[0].plot(asarray(t_sym[index]),asarray([v[t_index]*index/index for index in range(len(epsilon_vec))]), '--' ,color=(0,0,0),linewidth=2.0)    
axs[0].plot(t, u, '-', label="$u(\\tau)$" ,color=(0/256,68/256,27/256),linewidth=3.0)
axs[0].plot(t, v  , '-', label='$v(\\tau)$',color=(153/256,216/256,201/256),linewidth=3.0)
axs[0].plot(t_trans, u, '-', label='$u(\\hat{\\tau})$',color=(77/256,0/256,75/256),linewidth=3.0)
axs[0].plot(t_trans, v  , '-', label='$v(\\hat{\\tau})$',color=(140/256,150/256,198/256),linewidth=3.0)
axs[0].grid()
axs[0].legend(loc='best',prop={"size":20})
axs[1].plot(u, v,'-', label='Original population, $(u,v)$',color=(0/256,68/256,27/256),linewidth=3.0)
axs[1].plot(u_transformed,v, '-', label='Transformed population, $(\\hat{u},v)$',color=(77/256,0/256,75/256),linewidth=3.0)
axs[1].plot(asarray(u_sym[0]),asarray([v[0]*index/index for index in range(len(epsilon_vec))]), '--', label="$\\left.\\Gamma_{\\epsilon}^{u}\\right|_{\\epsilon=" +str(epsilon) + "}$" ,color=(0,0,0),linewidth=2.0)
for index,u_index in enumerate(list(u_indices)):
    axs[1].plot(asarray(u_sym[index]),asarray([v[u_index]*index/index for index in range(len(epsilon_vec))]), '--' ,color=(0,0,0),linewidth=2.0)    
axs[1].grid()
axs[1].legend(loc='best',prop={"size":20})
axs[2].plot(u, v,'-', label='Original population, $(u,v)$',color=(0/256,68/256,27/256),linewidth=3.0)
axs[2].plot(u,v_transformed, '-', label='Transformed population, $(u,\\hat{v})$',color=(77/256,0/256,75/256),linewidth=3.0)
axs[2].plot(asarray([u[0]*index/index for index in range(len(epsilon_vec))]), asarray(v_sym[0]), '--', label="$\\left.\\Gamma_{\\epsilon}^{v}\\right|_{\\epsilon=" +str(epsilon) + "}$" ,color=(0,0,0),linewidth=2.0)
for index,v_index in enumerate(list(v_indices)):
    axs[2].plot(asarray([u[v_index]*index/index for index in range(len(epsilon_vec))]),asarray(v_sym[index]), '--' ,color=(0,0,0),linewidth=2.0)    
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
f1.savefig('../Figures/LV_symmetries.png')
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
        plot_LaTeX_2D(asarray(t_sym[index]),asarray([u[t_index]*index/index for index in range(len(epsilon_vec))]),"../Figures/LV_symmetries/Input/tau_trans.tex","color=black,->,>=latex,densely dashed","$\\Gamma^{\mathrm{LV},\\tau}_{\\epsilon}$")
        plot_LaTeX_2D(asarray(t_sym[index]),asarray([v[t_index]*index/index for index in range(len(epsilon_vec))]),"../Figures/LV_symmetries/Input/tau_trans.tex","color=black,->,>=latex,densely dashed",[])         
    else:
        plot_LaTeX_2D(asarray(t_sym[index]),asarray([u[t_index]*index/index for index in range(len(epsilon_vec))]),"../Figures/LV_symmetries/Input/tau_trans.tex","color=black,->,>=latex,densely dashed",[])
        plot_LaTeX_2D(asarray(t_sym[index]),asarray([v[t_index]*index/index for index in range(len(epsilon_vec))]),"../Figures/LV_symmetries/Input/tau_trans.tex","color=black,->,>=latex,densely dashed",[]) 
plot_LaTeX_2D(t,u,"../Figures/LV_symmetries/Input/tau_trans.tex","color=clr_1,line width=1.5pt,","$u(\\tau)$")
plot_LaTeX_2D(t,v,"../Figures/LV_symmetries/Input/tau_trans.tex","color=clr_3,line width=1.5pt,","$v(\\tau)$")
plot_LaTeX_2D(t_trans,u,"../Figures/LV_symmetries/Input/tau_trans.tex","color=clr_2,line width=1.5pt,","$\\hat{u}(\\tau)$")
plot_LaTeX_2D(t_trans,v,"../Figures/LV_symmetries/Input/tau_trans.tex","color=clr_4,line width=1.5pt,","$\\hat{v}(\\tau)$")
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Scaled translations in u
plot_LaTeX_2D(asarray([u_sym_temp[0] for u_sym_temp in u_sym[0]]),asarray([v[0]*index/index for index in range(len(epsilon_vec))]),"../Figures/LV_symmetries/Input/u_trans.tex","color=black,->,>=latex,densely dashed","$\\Gamma^{\mathrm{LV},u}_{\\epsilon}$")
for index,u_index in enumerate(list(u_indices)):
    plot_LaTeX_2D(asarray([u_sym_temp[0] for u_sym_temp in u_sym[index]]),asarray([v[u_index]*index/index for index in range(len(epsilon_vec))]),"../Figures/LV_symmetries/Input/u_trans.tex","color=black,->,>=latex,densely dashed",[])
plot_LaTeX_2D(u,v,"../Figures/LV_symmetries/Input/u_trans.tex","color=clr_1,line width=1.5pt,","Original population")
plot_LaTeX_2D(asarray([u_temp[0] for u_temp in u_transformed]),v,"../Figures/LV_symmetries/Input/u_trans.tex","color=clr_2,line width=1.5pt,","Transformed population")
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Scaled translations in u
#plot_LaTeX_2D(asarray([u[v_indices[0]]*index/index for index in range(len(epsilon_vec))]),asarray([v_sym_temp[0] for v_sym_temp in v_sym[v_indices[0]]]),"../Figures/LV_symmetries/Input/v_trans.tex","color=black,->,>=latex,densely dashed","$\\Gamma^{\mathrm{LV},v}_{\\epsilon}$")
for index,v_index in enumerate(list(v_indices)):
    if index==0:
        plot_LaTeX_2D(asarray([u[v_index]*index/index for index in range(len(epsilon_vec))]),asarray([v_sym_temp[0] for v_index, v_sym_temp in enumerate(v_sym[index])]),"../Figures/LV_symmetries/Input/v_trans.tex","color=black,->,>=latex,densely dashed","$\\Gamma^{\mathrm{LV},v}_{\\epsilon}$")
    else:
        plot_LaTeX_2D(asarray([u[v_index]*index/index for index in range(len(epsilon_vec))]),asarray([v_sym_temp[0] for v_index, v_sym_temp in enumerate(v_sym[index])]),"../Figures/LV_symmetries/Input/v_trans.tex","color=black,->,>=latex,densely dashed",[])
        
plot_LaTeX_2D(u,v,"../Figures/LV_symmetries/Input/v_trans.tex","color=clr_3,line width=1.5pt,","Original population")
plot_LaTeX_2D(asarray([u[index] for index, v_temp in enumerate(v_transformed) if index>280 and index<456]),asarray([v_temp[0] for index, v_temp in enumerate(v_transformed) if index >280 and index<456]),"../Figures/LV_symmetries/Input/v_trans.tex","color=clr_4,line width=1.5pt,","Transformed population")
plot_LaTeX_2D(asarray([u[index] for index, v_temp in enumerate(v_transformed) if index>105 and index<175]),asarray([v_temp[0] for index, v_temp in enumerate(v_transformed) if index >105 and index<175]),"../Figures/LV_symmetries/Input/v_trans.tex","color=clr_4,line width=1.5pt,",[])
