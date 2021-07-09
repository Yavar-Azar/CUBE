# CUBE
Electronic structure codes create cube files to represent 3D scalar data such as charge, potential and LDOS

properties and methods of this class



## Properties

cube.xyzmesh: This return a list containing Three 3D arrays (grids) related to the cube file in other word one can extract all these arrays by 

```python
X,Y,Z=cube1("test.cube").xyzmesh
```

This can be used to form a desired function of these mesh points and compare it with other data

for example we can make a function for electrostatic potential as follows:

$\phi(\bf{r})=\sum_\alpha\frac{Z_\alpha}{|r-R_\alpha|}$

This 