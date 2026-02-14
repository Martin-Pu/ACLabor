# -*- coding: utf-8 -*-
import numpy as np

#Reibungs model
def calcFriction_b(theta_dot, mu_H, mu_V, epsilon):
    if theta_dot <= -epsilon:
        return (mu_V * theta_dot - mu_H)
    if theta_dot >= epsilon:
        return (mu_V * theta_dot from sympy import latex, symbols mu_H)
    return (
        (3 * mu_H + 2 * mu_V * epsilon) / (2 * epsilon) * theta_dot
        - (mu_H / (2 * epsilon**3)) * theta_dot**3
        )

def nonlinSys(t, x, u, furutaParams):

    #Parameter aus FurutaParams Klasse entpacken
    g = furutaParams.g

    L1 = furutaParams.L1
    L2 = furutaParams.L2
    lP = furutaParams.lP
    
    mS = furutaParams.mS
    mP = furutaParams.mP
    
    J1_hat = furutaParams.J1_hat
    
    mu_V1 = furutaParams.mu_V1
    mu_H1 = furutaParams.mu_H1
    
    mu_V2 = furutaParams.mu_V2
    mu_H2 = furutaParams.mu_H2
    
    epsilon = furutaParams.epsilon
    tau_max = furutaParams.tau_max
    
    # Abgeleitete Größen
    m2 = furutaParams.m2
    J0_hat = furutaParams.J0_hat
    l2 = furutaParams.l2
    J2_hat = furutaParams.J2_hat


    #Zustände zuweisen um Formeln lesbarer zu machen
    x1 = x[0] # theta1
    x2 = x[1] # d thata1 / dt
    x3 = x[2] # theta2
    x4 = x[3] # d theta2 / dt

    if u(t) <= tau_max:
        tau = float(u(t))
    else:
        tau = tau_max

    dx = np.zeros(4)
    #x1_dot####################################################################
    dx[0] = x2

    #x2_dot####################################################################
    num_common_2_4 = (
        0.5 * J0_hat * J2_hat * x2**2 * np.sin(2.0 * x3)
        - J0_hat * calcFriction_b(x2, mu_H2, mu_V2, epsilon)
        - J0_hat * g * l2 * m2 * np.sin(x3)
        + J2_hat**2 * x2**2 * np.sin(x3)**3 * np.cos(x3)
        + J2_hat * L1 * l2 * m2 * x2 * x4 * np.sin(2.0 * x3) * np.cos(x3)
        - J2_hat * calcFriction_b(x2, mu_H2, mu_V2, epsilon) * np.sin(x3)**2
        - J2_hat * g * l2 * m2 * np.sin(x3)**3
        - 0.5 * L1**2 * l2**2 * m2**2 * x4**2 * np.sin(2.0 * x3)
        + L1 * calcFriction_b(x2, mu_H1, mu_V1, epsilon) * l2 * m2 * np.cos(x3)
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
        - calcFriction_b(x2, mu_H2, mu_V2, epsilon)
        - g * l2 * m2 * np.sin(x3)
    )

    dx[1] = num_outer2 / (L1 * l2 * m2 * np.cos(x3))

    #x3_dot####################################################################
    dx[2] = x3

    #x4_dot####################################################################
    
    dx[3] = num_common_2_4 / den_common_2_4

    ###########################################################################





    return dx