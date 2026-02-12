# -*- coding: utf-8 -*-
import params as st
import sympy as sp
import numpy as np

def stepOnOff(t_on = 0.4, t_off = 0.5, u_on = 10):
    uIn = lambda t, t_on=t_on, t_off=t_off, u=u_on: u if (t_on <= t < t_off) else 0.0
    return uIn

#validation with csv data
def u_validation(t):
    if 0 < t and t < 1:
        return 0.62*t
    if 1 <= t and t<=1.4:
        return -4.12*t + 4.74
    if 1.4 < t and t<2:
        return 1.72*t - 3.43
    return 0
        


