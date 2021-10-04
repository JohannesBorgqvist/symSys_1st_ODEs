from sympy import *


# Create an arbitrary coefficient
k = Symbol('k')
# Create the variable which we integrate with
t = Symbol('t')
# Create the dummy variable with which we integrate with
s = Symbol('s')
# Create an arbitrary function
f = symbols('f',cls=Function)
# Define an expression for the integrand
integrand = Derivative(f(t),t)*exp(k*t)
# Create an integral
integral = Integral(integrand,t)
# Print the integral
print("\n\tInitial Integral:\n\t\t%s\n"%(latex(integral)))
# Evaluate integral
integral = integral.doit()
# Print the evaluated integral
print("\n\tEvaluated Integral:\n\t\t%s\n"%(latex(integral)))
