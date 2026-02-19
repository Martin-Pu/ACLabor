import numpy as np
import sympy as sp

# Symbole
th1, th1d, th2, th2d, tau = sp.symbols("th1 th1d th2 th2d tau", real=True)

def make_non_linear_furuta_accel_expr(
    *,
    g, L1, l2, J0_hat, J2_hat, m2,
    mu_V1, mu_H1, mu_V2, mu_H2,
    epsilon,
    simplify_expr: bool = True
):
    """
    Gives the sympy expressions for theta1_ddot and theta2_ddot
    """

    # Reibung: mu_V*v + mu_H*tanh(v/epsilon)
    b1 = mu_V1 * th1d + mu_H1 * sp.tanh(th1d / epsilon)
    b2 = mu_V2 * th2d + mu_H2 * sp.tanh(th2d / epsilon)
    
    # Massenmatrix
    M11 = J0_hat + J2_hat * sp.sin(th2) ** 2
    M12 = m2 * L1 * l2 * sp.cos(th2)
    M = sp.Matrix([[M11, M12],
                   [M12, J2_hat]])

    # Rechte Seite
    r1 = (
        tau
        + m2 * L1 * l2 * sp.sin(th2) * th2d**2
        - J2_hat * sp.sin(2 * th2) * th1d * th2d
        - b1
    )
    r2 = (
        sp.Rational(1, 2) * J2_hat * sp.sin(2 * th2) * th1d**2
        - m2 * g * l2 * sp.sin(th2)
        - b2
    )
    r = sp.Matrix([r1, r2])

    # Lösung ohne inv(): stabiler
    th1dd, th2dd = M.LUsolve(r)

    if simplify_expr:
        th1dd = sp.simplify(th1dd)
        th2dd = sp.simplify(th2dd)

    return (th1dd, th2dd)

def make_furuta_accel_lambdas(accel_exprs):
    """
    Lambdifys sympy expressions with the symbols th1, th1d, th2, th2d and tau
    """
    return sp.lambdify((th1, th1d, th2, th2d, tau), accel_exprs, modules="numpy")
    