# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 08:25:34 2020

@author: Sookyung Lee
MA 514 Project 3 Problem 1
"""
import math
import numpy

def int_taylor_series(t):
    func = 0
    index = 0
    term = 1
    while abs(term) > 10**(-8):
        term = (-1)**index*t**((index*2)+1)/(math.factorial(index) * ((index*2)+1))
        func += term
        index += 1
        
    return func
    

def erf(t):
    func = 0
    index = 0
    term = 1
    while abs(term) > 10**(-8):
        term = (-1)**index*t**((index*2)+1)/(math.factorial(index) * ((index*2)+1))
        func += term
        index += 1
    func = 2/math.sqrt(math.pi) * func
    return func
#Run the numerical experiment
for t in numpy.arange (0, 1.1, 0.1):
    # print('t:', t)
    # print('Integral value:', int_taylor_series(t))
    # print('erf(x): ', erf(t))
    print(int_taylor_series(t) - erf(t))