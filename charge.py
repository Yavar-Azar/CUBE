#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 15 23:54:07 2021

@author: yavar
"""

import matplotlib.pyplot as plt
import numpy as np


from CUBE_analyzer import cube1





test = cube1("chargedensity.cube")

data = test.grid3d

gridnumber = len(data[0,0,:])

halfnum = int((gridnumber-1)/2)

x = np.linspace(-5,5, gridnumber)[halfnum :]

xdata = data[halfnum :, halfnum, halfnum]

chargesum = test.diffvol*np.sum(data)

print("total charge is :  "+str(chargesum))






