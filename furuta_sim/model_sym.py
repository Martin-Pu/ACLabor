# -*- coding: utf-8 -*-
import numpy as np
import params as st
import sympy as sp

class FurutaPendel:
    def __init__(self):
        self.params = st.FurutaParams()
        self.x2dotFunc = None # will be set by lambdify
        self.x4dotFunc = None # will be set by lambdify
        self.updateModelFunctions()

    """updates the functions for x2dot and x4dot"""
    def updateModelFunctions(self):
        x1, x2, x3, x4, x2dot, x4dot = sp.symbols('x1 x2 x3 x4 x2dot x4dot')
        m2, L1, l2, J0_hat, J2_hat, b1, tau, b2, g = sp.symbols('m2 L1 l2 J0 J2 b1 tau b2 g')
        eqn1 = sp.Eq(m2*L1*l2*(sp.cos(x3)*x4dot - sp.sin(x3)*x4**2) + J0_hat*x2dot + J2_hat*(sp.sin(x3)**2*x2dot + sp.sin(2*x3)*x2*x4) + b1 , tau)
        eqn2 = sp.Eq(m2*L1*l2*sp.cos(x3)*x2dot - J2_hat*(0.5*sp.sin(2*x3)*x2**2 - x4dot) + b2 + m2*l2*g*sp.sin(x3), 0)

        # solve eqn2 for x2dot
        x2dot_expr = sp.solve(eqn2, x2dot)[0] 
        # substitute into eqn1
        eqn1_sub = eqn1.subs(x2dot, x2dot_expr)       
        # solve for x4dot
        x4dot_solution = sp.solve(eqn1_sub, x4dot)[0]
        # back-substitude into x2dot
        x2dot_solution = x2dot_expr.subs(x4dot, x4dot_solution)

        #set constant values for lambda function and lambdify
        x4dot_solution_const = x4dot_solution.subs({
            m2: self.params.g,
            L1: self.params.L1,
            l2: self.params.l2,
            J0_hat: self.params.J0_hat,
            J2_hat: self.params.J2_hat,
            g: self.params.g
        })
        self.x4dotFunc = sp.lambdify((x1, x2, x3, x4, b1, b2, tau), x4dot_solution_const)

        #set constant values for lambda function and lambdify
        x2dot_solution_const = x2dot_solution.subs({
            m2: self.params.g,
            L1: self.params.L1,
            l2: self.params.l2,
            J0_hat: self.params.J0_hat,
            J2_hat: self.params.J2_hat,
            g: self.params.g
        })
        self.x2dotFunc = sp.lambdify((x1, x2, x3, x4, b1, b2, tau), x2dot_solution_const)

    """friction model of the pendulum"""
    def calcFriction_b(self, theta_dot, mu_H, mu_V, epsilon):
        if theta_dot <= -epsilon:
            return (mu_V * theta_dot - mu_H)
        if theta_dot >= epsilon:
            return (mu_V * theta_dot + mu_H)
        return (
            (3 * mu_H + 2 * mu_V * epsilon) / (2 * epsilon) * theta_dot
            - (mu_H / (2 * epsilon**3)) * theta_dot**3
            )
    
    def calc_dx(self, t, x, u):
    
        #Zustände zuweisen um Formeln lesbarer zu machen
        x1 = x[0] # theta1
        x2 = x[1] # d thata1 / dt
        x3 = x[2] # theta2
        x4 = x[3] # d theta2 / dt
    
        #stellgrößen beschränkung
        if u(t) <= self.params.tau_max:
            tau = u(t)
        else:
            tau = tau_max
    
        #calc friction
        b1 = self.calcFriction_b(x2, self.params.mu_H1, self.params.mu_V1, self.params.epsilon)
        b2 = self.calcFriction_b(x4, self.params.mu_H2, self.params.mu_V2, self.params.epsilon)
    
        dx = np.zeros(4)
        dx[0] = x2
        dx[1] = self.x2dotFunc(x1, x2, x3, x4, b1, b2, tau)
        dx[2] = x4
        dx[3] = self.x4dotFunc(x1, x2, x3, x4, b1, b2, tau)
    
        return dx