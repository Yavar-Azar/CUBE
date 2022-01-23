#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 15 23:54:07 2021

@author: yavar
"""

import matplotlib.pyplot as plt
import numpy as np


from CUBE_analyzer import cube1



test1 = cube1("charge1.cube")
test2 = cube1("charge2.cube")




data1= test1.grid3d
data2 =test2.grid3d


res1 = np.sum(data1)*test1.diffvol
res2 = np.sum(data2)*test2.diffvol

print(res1, np.max(data1))
print(res2, np.max(data2))


