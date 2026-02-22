import numpy as np
import sympy as sp
import furuta_params

# Symbole
th1, th1d, th2, th2d, tau = sp.symbols("th1 th1d th2 th2d tau", real=True)

def make_non_linear_furuta_accel_expr(
    params,
    simplify_expr: bool = True
):
    """
    Gives the sympy expressions for theta1_ddot and theta2_ddot
    """

    if not isinstance(params, furuta_params.FurutaParams):
        raise TypeError(f"params must be FurutaParams, got {type(params).__name__}")

    g = params.g

    # Geometrie
    L1 = params.L1
    L2 = params.L2
    lP = params.lP
    l2 = params.l2
    
    # Massen
    mS = params.mS
    mP = params.mP
    m2 = params.m2
    
    # Trägheitsmomente
    J1_hat = params.J1_hat
    J0_hat = params.J0_hat
    J2_hat = params.J2_hat
    
    # Reibungsparameter
    mu_V1 = params.mu_V1
    mu_H1 = params.mu_H1
    mu_V2 = params.mu_V2
    mu_H2 = params.mu_H2
    
    # Sonstiges
    epsilon = params.epsilon
    tau_max = params.tau_max

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

def get_non_linear_furuta_system(params):
    """
    Lambdifys sympy expressions with the symbols th1, th1d, th2, th2d and tau
    """

    # non linear expressions
    th1dd_expr, th2dd_expr = make_non_linear_furuta_accel_expr(
        params,
        simplify_expr=False,   # wenn's beim Start zu lange dauert -> False
    )
    
    th1dd_lambda = sp.lambdify((th1, th1d, th2, th2d, tau), th1dd_expr, modules="numpy")
    th2dd_lambda = sp.lambdify((th1, th1d, th2, th2d, tau), th2dd_expr, modules="numpy")

    def non_linear_sys(t, x, u):
        x1, x2, x3, x4 = x
        x1d = x2
        x2d = th1dd_lambda(x1, x2, x3, x4, u)
        x3d = x4
        x4d = th2dd_lambda(x1, x2, x3, x4, u)
        return float(x1d), float(x2d), float(x3d), float(x4d)

    return non_linear_sys

def get_linear_furuta(params, operating_point):
    """
    Returns A Matrix and b Vektor for a specified operating point
    ordered like: (th1_op, th1d_op, th2_op, th2d_op, tau_op)
    """
    th1_op, th1d_op, th2_op, th2d_op, tau_op = operating_point
    
    th1dd_expr, th2dd_expr = make_non_linear_furuta_accel_expr(
        params,
        simplify_expr=False,   # wenn's beim Start zu lange dauert -> False
    )
    
    f1 = th1d
    f2 = th1dd_expr
    f3 = th2d
    f4 = th2dd_expr
    
    f = sp.Matrix([f1, f2, f3, f4])
    x = sp.Matrix([th1, th1d, th2, th2d])
    u = sp.Matrix([tau])
    
    A = f.jacobian(x)
    b = f.jacobian(u)
    
    A_lin = A.subs({th1:th1_op, th1d:th1d_op, th2:th2_op, th2d:th2d_op, tau:tau_op})
    b_lin = b.subs({th1:th1_op, th1d:th1d_op, th2:th2_op, th2d:th2d_op, tau:tau_op})

    A_np = np.array(A_lin.tolist(), dtype=float)
    b_np = np.array(b_lin.tolist(), dtype=float)

    return A_np, b_np


def get_RNF_furuta(params, operating_point):
    """
    Returns A_RNF and Q Matrix and  for a specified operating point
    """
    A, b = get_linear_furuta(params, operating_point)

    n = A.shape[0]
    
    # Steuerbarkeitsmatrix S = [B, AB, A^2B, ...]
    S_blocks = [np.linalg.matrix_power(A, k) @ b for k in range(n)] #liste von spalten vektoren
    S = np.hstack(S_blocks) #spaletenvektoren werden zu matrix zusammen gefügt
    S_inv = np.linalg.inv(S)

    q1T = S_inv[-1, :].reshape(1, -1) #reshape to explicitly make it a row vector

    Q_blocks = [q1T @ np.linalg.matrix_power(A, k)for k in range(n)] #liste von zeilen vektoren
    Q = np.vstack(Q_blocks)
    Q_inv = np.linalg.inv(Q)

    tol = 1e-13   # wählbar, z.B. 1e-12, 1e-9, je nach Skalierung

    A_RNF = Q @ A @ Q_inv
    A_RNF[np.abs(A_RNF) < tol] = 0.0

    return A_RNF, Q

def get_a_coefficients_furuta(params, operating_point):
    """
    Returns the a coefficients of the system for a specified operating point
    """
    A_RNF, Q = get_RNF_furuta(params, operating_point)
    a0, a1, a2, a3 = -A_RNF[-1]
    return a0, a1, a2, a3

def get_b_coefficients_furuta(cT, params, operating_point):
    """
    Returns the b coefficients of the system for a specified operating point
    and a specified output vector cT
    """
    A_RNF, Q = get_RNF_furuta(params, operating_point)
    Q_inv = np.linalg.inv(Q)
    
    cT = cT[-1, :].reshape(1, -1) #reshape to explicitly make it a row vector

    cT_RNF = (cT @ Q_inv)
    b0, b1, b2, b3 = cT_RNF[0]

    return b0, b1, b2, b3



    
    