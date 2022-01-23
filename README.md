# CUBE
Electronic structure codes create cube files to represent 3D scalar data such as charge, potential and LDOS

properties and methods of this class





## Discussions

There are several challenges against using cube files as reliable outputs

- First of all the maximum value are strongly dependent on the origin of grid

- The summation of all values is also grid dependent (making denser grid can solve problem somehow)

  

  

 

## Properties

cube.xyz mesh: This return a list containing Three 3D arrays (grids) related to the cube file in other word one can extract all these arrays by 

```python
X,Y,Z=cube1("test.cube").xyzmesh
```

This can be used to form a desired function of these mesh points and compare it with other data

for example we can make a function for electrostatic potential as follows:

$\phi(\bf{r})=\sum_\alpha\frac{Z_\alpha}{|r-R_\alpha|}$

This 





## TECHNICALITIES

- gaussian cubegen

  you need an input file like this 

  ```tex
   -1   -6.512752   -6.512752   -6.512752
   -83    0.168040    0.000000    0.000000
   79    0.000000    0.168040    0.000000
   79    0.000000    0.000000    0.168040
  ```

  

to create your own cube file with cubegen by following command

```bash
cubegen 1 fdensity=scf Test.FChk result.cube -1 < value.txt
```





- NWCHEM input to generate inputs

  In the case of nwchem we can change 

  ```shell
  start n1  
  geometry  units au noautosym noautoz nocenter
    C  0.0 0.0 0.0
    O  2.1 0.0 0.0
  end  
  basis; * library cc-pvdz;end 
  dft
   iterations 50
   vectors  output co.movecs
   print  kinetic_energy  
   xc xpbe96 cpbe96  
  end  
  dplot  
    TITLE charge  
    vectors co.movecs  
  LimitXYZ units au  
  -5. 5.0 79    
  -5. 5.0 79   
  -5. 5.0 79 
   spin total 
   gaussian  
   output chargedensity2.cube  
  end  
  
  ```

  