# -*- coding: utf-8 -*-
import numpy as np

#Reibungs model
def calcFriction_b(theta_dot, mu_H, mu_V, epsilon):
    if theta_dot <= -epsilon:
        return (mu_V * theta_dot - mu_H)
    if np.absolute(theta_dot) < epsilon:
        return (
            (3 * mu_H + 2 * mu_V * epsilon) / (2 * epsilon * theta_dot)
            - (mu_H / (2 * epsilon**3)) * theta_dot**3
        )
    if theta_dot >= epsilon:
        return (mu_V * theta_dot + mu_H)

def nonlinSys(t, x, u, params):
    x1 = x[0] # theta1
    x2 = x[1] # d thata1 / dt
    x3 = x[2] # theta2
    x4 = x[3] # d theta2 / dt

    g, L1, L2, lp, mS, mP, J1_hat, mu_V1,\
    mu_H1, mu_V2, mu_H2, epsilon, tau_max,\
    m2, J0_hat, J2_hat, l2 = params

    if u <= tau_max:
        tau = u
    else
        tau = tau_max

    dx = np.zeros(4)
    #x1_dot####################################################################
    dx[0] = x2

    #x2_dot####################################################################
    num_common_2_4 = (
        0.5 * J0_hat * J2_hat * x2**2 * np.sin(2.0 * x3)
        - J0_hat * calcFriction_b(x2, mu_H2, mu_V2)
        - J0_hat * g * l2 * m2 * np.sin(x3)
        + J2**2 * x2**2 * np.sin(x3)**3 * np.cos(x3)
        + J2_hat * L1 * l2 * m2 * x2 * x4 * np.sin(2.0 * x3) * np.cos(x3)
        - J2_hat * b2 * np.sin(x3)**2
        - J2_hat * g * l2 * m2 * np.sin(x3)**3
        - 0.5 * L1**2 * l2**2 * m2**2 * x4**2 * np.sin(2.0 * x3)
        + L1 * calcFriction_b(x2, mu_H1, mu_V1) * l2 * m2 * np.cos(x3)
        - L1 * l2 * m2 * tau * np.cos(x3)
    )

    den_common_2_4 = (
        J0_hat * J2_hat
        + J2_hat**2 * np.sin(x3)**2
        - L1**2 * l2**2 * m2**2 * np.cos(x3)**2
    )
    
    num_outer2 = (
        0.5 * J2_hat * x2**2 * np.sin(2.0 * x3)
        - J2_hat * (num_common_2_4 / den_common_2_4)
        - b2
        - g * l2 * m2 * np.sin(x3)
    )

    dx[1] = num_outer2 / (L1 * l2 * m2 * np.cos(x3))

    #x3_dot####################################################################
    dx[2] = x3

    #x4_dot####################################################################
    
    dx[3] = num_common_2_4 / den_common_2_4

    ###########################################################################





    return dx