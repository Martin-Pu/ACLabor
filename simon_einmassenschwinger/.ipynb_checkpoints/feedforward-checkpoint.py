# -*- coding: utf-8 -*-
import params as st
import sympy as sp
import numpy as np

def stepOnOff(t_on = 100, t_off = 2000, u_on = 10):
    uIn = lambda t, t_on=t_on, t_off=t_off, u=u_on: u if (t_on <= t < t_off) else 0.0
    return uIn


