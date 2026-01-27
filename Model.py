# -*- coding: utf-8 -*-
#BeispielCode aus Brückenkran
import numpy as np

def nonlinSys(t, x, u, params):
    D = x[0]
    dD = x[1]
    th = x[2]
    dth = x[3]

    cs = np.cos(th)
    sn = np.sin(th)

    g, M, m , l = params

    dx = np.zeros(4)
    h = l * dth ** 2 + g * cs
    dx[0] = dD
    dx[1] = (u(t) + m * sn * h) / (M + m * sn ** 2)
    dx[2] = dth
    dx[3] = -(dx[1] * cs + g * sn) / l

    return dx

def linSys(t, x, u, params):
    D = x[0]
    dD = x[1]
    th = x[2]
    dth = x[3]

    g, M, m , l = params

    dx = np.zeros(4)
    dx[0] = dD;
    dx[1] = (u(t) + m * g * th) / M;
    dx[2] = dth;
    dx[3] = -(dx[1] + g * th) / l;

    return dx