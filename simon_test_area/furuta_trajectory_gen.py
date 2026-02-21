import furuta_systems
import furuta_params
import sympy as sp
import numpy as np

def generate_trajectorys(eta_t0, eta_t1, t0, t1, params, operating_point):
    """
    generates trajectorys for u and x that are only valid from t0 to t1
    """
    a0, a1, a2, a3 = furuta_systems.get_a_coefficients_furuta(params, operating_point)
    _, Q = furuta_systems.get_RNF_furuta(params, operating_point)
    Q_inv = np.linalg.inv(Q)

    # Symbol definieren
    tau, t_sym = sp.symbols('tau t')
    
    # phi_4 definieren
    phi = 126*tau**5 - 420*tau**6 + 540*tau**7 - 315*tau**8 + 70*tau**9
    eta = eta_t0 + (eta_t1 - eta_t0) * phi
    eta = eta.subs(tau, (t_sym - t0)/(t1 - t0))
    
    eta_d1 = sp.diff(eta, t_sym, 1)
    eta_d2 = sp.diff(eta_d1, t_sym, 1)
    eta_d3 = sp.diff(eta_d2, t_sym, 1)
    eta_d4 = sp.diff(eta_d3, t_sym, 1)
    
    eta_ref = sp.lambdify(t_sym, eta, 'numpy')
    eta_d1_ref = sp.lambdify(t_sym, eta_d1, 'numpy')
    eta_d2_ref = sp.lambdify(t_sym, eta_d2, 'numpy')
    eta_d3_ref = sp.lambdify(t_sym, eta_d3, 'numpy')
    eta_d4_ref = sp.lambdify(t_sym, eta_d4, 'numpy')
    
    
    u_ref = lambda t: eta_d4_ref(t) + a3 * eta_d3_ref(t) + a2 * eta_d2_ref(t) + a1 * eta_d1_ref(t) + a0 * eta_ref(t)
    
    x_ref = lambda t: Q_inv @ np.array([eta_ref(t), eta_d1_ref(t), eta_d2_ref(t), eta_d3_ref(t)])

    return u_ref, x_ref