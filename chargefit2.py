#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 10:11:46 2021

@author: yavar
"""


import matplotlib.pyplot as plt
import numpy as np


from CUBE_analyzer import cube1



from scipy.optimize import curve_fit



def funcradial(x,  n, alpha, beta,epsilon):
    return np.exp(-alpha * (x**2))*((beta*x+epsilon)**-n)

test = np.loadtxt("data200_half.txt")
x = test[:,0]
y = test[:,1]

yf = funcradial(x,  10.91153357, -0.66923204, 1.0, 0.6130577)


plt.plot(x, y , lw = 3)
plt.plot(x, yf, '--', lw=2)



popt, pcov = curve_fit(funcradial, x, y, p0=[10, -0.76923204, 1.2, 0.6130577], maxfev=100000)


fitted = funcradial(x, *popt)

diff = fitted - y 
print(np.sum(diff))
print(popt)

plt.plot(x, y , lw = 7, alpha=0.7)
plt.plot(x, yf, '--', lw=2)
plt.plot(x, fitted, '.', lw=3)
plt.show()