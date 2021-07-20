using SymPoly
using SymbolicUtils
using DynamicPolynomials

@syms y

@polyvar x

eq = y^2 - 4y + 3

typeof(eq)

p = poly(eq, y => x)

typeof(p)

sym(p, x => y)
