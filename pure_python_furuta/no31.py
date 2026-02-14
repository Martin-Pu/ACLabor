import sympy as sp
x1, x2, x3, x4, x2dot, x4dot = sp.symbols('x1 x2 x3 x4 x2dot x4dot')
m2, L1, l2, J0, J2, b1, tau, b2, g = sp.symbols('m2 L1 l2 J0 J2 b1 tau b2 g')
eqn1 = sp.Eq(m2*L1*l2*(sp.cos(x3)*x4dot - sp.sin(x3)*x4**2) + J0*x2dot + \
J2*(sp.sin(x3)**2*x2dot + sp.sin(2*x3)*x2*x4) + b1 , tau)
eqn2 = sp.Eq(m2*L1*l2*sp.cos(x3)*x2dot - J2*(0.5*sp.sin(2*x3)*x2**2 - x4dot) \
+ b2 + m2*l2*g*sp.sin(x3), 0)

# solve eqn2 for x2dot
x2dot_expr = sp.solve(eqn2, x2dot)[0]

# substitute into eqn1
eqn1_sub = eqn1.subs(x2dot, x2dot_expr)

# solve for x4dot
x4dot_solutions = sp.solve(eqn1_sub, x4dot)

# back-substitute to get x2dot solutions
solutions = [
    {
        x4dot: sol,
        x2dot: x2dot_expr.subs(x4dot, sol)
    }
    for sol in x4dot_solutions
]

print(solutions)
