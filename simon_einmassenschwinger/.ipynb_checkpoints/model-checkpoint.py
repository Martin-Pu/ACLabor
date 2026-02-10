# -*- coding: utf-8 -*-
import numpy as np

def linSys(t, x, u, params):
    x1 = x[0]
    x2 = x[1]
    d, k, m, _ = params

    dx = np.zeros(2)
    dx[0] = x2
    dx[1] = 1/m * ( -k*x1 -d*x2 +u(t) )

    return dx
