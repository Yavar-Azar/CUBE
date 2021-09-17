#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 29 10:15:04 2019

@author: yavar001
"""

import atexit,os
import sys, string
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import re
import collections
from ase import Atoms, Atom
from ase.units import Bohr
from ase.io import read, write
from ase.build import make_supercell




from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster import KMeans
from sklearn.mixture import  GaussianMixture

from sklearn import datasets





cpk_colors={
    "1"	:	"#FFFFFF",  
"2"	:	"#D9FFFF",
"3"	:	"#CC80FF",
"4"	:        "#C2FF00",
"5"	:	"#FFB5B5",
"6"	:	"#909090",
"7"	:        "#3050F8",
"8"	:        "#FF0D0D",
"9"	:	"#90E050",
"10"    :  	"#B3E3F5",
"11"	:	"#AB5CF2",
"12"	:        "#8AFF00",
"13"	:	"#BFA6A6",
"14"	:	"#F0C8A0",
"15"	:        "#FF8000",
"16"	:	"#FFFF30",
"17"	:        "#1FF01F",
"18"	:	"#80D1E3",
"19"	:	"#8F40D4",
"20"	:        "#3DFF00",
"21"	:	"#E6E6E6",
"22"	:	"#BFC2C7",
"23"	:	"#A6A6AB",
"24"	:	"#8A99C7",
"25"	: 	"#9C7AC7",
"26"	:	"#E06633",
"27"    :  	"#F090A0",
"28"	:        "#50D050",
"29"	:	"#C88033",
"30"	:	"#7D80B0",
"31"	:	"#C28F8F",
"32"	:	"#668F8F",
"33"	:	"#BD80E3",
"34"	:        "#FFA100",
"35"	:        "#A62929",
"36"	:	"#5CB8D1",
"37"	:	"#702EB0",
"38"	:        "#00FF00",
"39"	:	"#94FFFF",
"40"	:	"#94E0E0",
"41"	:	"#73C2C9",
"42"	:	"#54B5B5",
"43"	:	"#3B9E9E",
"44"	:	"#248F8F",
"45"	:	"#0A7D8C",
"46"	:        "#006985",
"47"	:	"#C0C0C0",
"48"	:        "#FFD98F",
"49"	:	"#A67573",
"50"	:	"#668080",
"51"	:	"#9E63B5",
"52"	:        "#D47A00",
"53"	:        "#940094",
"54"	:	"#429EB0",
"55"	:        "#57178F",
"56"	:        "#00C900",
"57"	: 	"#70D4FF",
"58"	:	"#FFFFC7",
"59"	:	"#D9FFC7",
"60"	:	"#C7FFC7",
"61"	:	"#A3FFC7",
"62"	:	"#8FFFC7",
"63"	:	"#61FFC7",
"64"	:	"#45FFC7",
"65"	:	"#30FFC7",
"66"	:	"#1FFFC7",
"67"	:        "#00FF9C",
"68"	:        "#00E675",
"69"	:        "#00D452",
"70"	:        "#00BF38",
"71"	:        "#00AB24",
"72"	:	"#4DC2FF",
"73"	:	"#4DA6FF",
"74"	:	"#2194D6",
"75"	:	"#267DAB",
"76"	:	"#266696",
"77"	:        "#175487",
"78"	:	"#D0D0E0",
"79"	:	"#FFD123",
"80"	:	"#B8B8D0",
"81"	:        "#A6544D",
"82"	:        "#575961",
"83"	:	"#9E4FB5",
"84"	:        "#AB5C00",
"85"	:        "#754F45",
"86"	:	"#428296",
"87"	:        "#420066",
"88"	:        "#007D00",
"89"	:	"#70ABFA",
"90"	:        "#00BAFF",
"91"	:        "#00A1FF",
"92"	:        "#008FFF",
"93"	:        "#0080FF",
"94"	:        "#006BFF",
"95"	:        "#545CF2",
"96"	:	"#785CE3",
"97"	:	"#8A4FE3",
"98"	:	"#A136D4",
"99"	:	"#B31FD4",
"100"	:	"#B31FBA",
"101"	:	"#B30DA6",
"102"	:	"#BD0D87",
"103"	:	"#C70066",
"104"	:	"#CC0059",
"105"	:	"#D1004F",
"106"	:	"#D90045",
"107"	:	"#E00038",
"108"	:	"#E6002E",
"109"	:	"#EB0026"

    }


















class cube1(object):
      
    def __init__(self,cube):
        self.cube=cube
        
 #   @staticmethod
 #   def validate_log(self):
        if os.path.exists(self.cube):
            print("your file is in path")
            if os.stat(self.cube).st_size<5:
                raise ValueError("File Size NOT Normal")
            else:
                print("it seems normal file")
        else:
            print("!!!!!!!!!!file not found, Enter correct path..........")
        
    @staticmethod            
    def readcube(filename):
        
        
        #ENTER THE NAME OF CUBE_DIFF FILE
        N=np.zeros(3,dtype=int)
        v=np.zeros((3,3))
        cell=np.zeros((3,3))
        f = open(filename,'r')
 #put here the interval you want
        next(f)
        next(f)
        line = next(f)
          #Get the number of atoms if we want to store it
        natoms = int(line.split()[0])
        origx,origy,origz=line.split()[1],line.split()[2],line.split()[3]
        origx=float(origx)
        origy=float(origy)
        origz=float(origz)
          #This gets the nx,ny,nz info of the charge density
          #As well as the differential volume
        for i in range(0,3):
            line = next(f).split()
            N[i] = int(line[0])
            for j in range(1,4):
                v[i][j-1] = float(line[j])  
            cell[i]=v[i,:]*N[i]                
        dv=v[0][0]*v[1][1]*v[2][2]
        mesh=[v[0][0],v[1][1],v[2][2]]
        #  #As of now we dont care about the positions of the atoms,
        #  #But if you did you could read them here:
 
        
        coords=np.zeros((np.abs(natoms),4))
        for i in range(0,natoms):
            line=next(f).split()
            for k in range(0,4):
                coords[i,k]=float(line[k+1])#
                
        #coords=np.loadtxt("mynuc.txt",usecols=(0,2,3,4))
        sumposcharge=np.sum(coords[:,0])
                
        origin=[origx, origy, origz]
        #xposcenter=np.sum(coords[:,0]*coords[:,1]/sumposcharge)
        #yposcenter=np.sum(coords[:,0]*coords[:,2]/sumposcharge)
        #zposcenter=np.sum(coords[:,0]*coords[:,3]/sumposcharge)
        
        # 
        #
        #This reads the data into a 1D array of size nx*ny*nz
        grid=np.zeros((N[0]*N[1]*N[2],1))
        count = 0
        for i in f:
            for j in i.split():
                grid[count,0] = float(j)
                count+=1
        f.close()
        
        
        
    
        f2 = open(filename, 'r')
        heads=f2.readlines()[0:6+natoms]
        f2.close()
        
           
        
        grid3d=np.reshape(grid, (N[0], N[1], N[2]))
        
        
        ############  POSITION VECTOR FOR GRID POINTS
        gridvect=np.zeros((N[0]*N[1]*N[2], 3))
        for i in range(N[0]):
            for j in range(N[1]):
                for k in range(N[2]):
                    numb=(N[2]*N[1]*i+N[2]*j+k)
                    gridvect[numb]=[i*mesh[0],j*mesh[1],k*mesh[2]]
        
        
        
                    
            
        
        gdims=[N[0], N[1], N[2]]
        
        ######## !!!!!!!!!!!!!!
        #####   ONLY ORTHOROMBIC CELL
    
        
        
        stepx=cell[0,0]/N[0]
        stepy=cell[1,1]/N[1]
        stepz=cell[2,2]/N[2]
        
        gridpx=np.arange(0, N[0])*stepx+origx
        gridpy=np.arange(0, N[1])*stepy+origy
        gridpz=np.arange(0, N[2])*stepz+origz
        
        
        spacemesh=np.meshgrid(gridpx, gridpy, gridpz)
        
        xyz=[]
        xyz.append(natoms)
        xyz.append('  ')
        for i in range(natoms):
            xyz.append([Atom(int(coords[i,0])).symbol, coords[i,1], coords[i,2], coords[i,3]])        
        
        
        
        datasum = dv*np.sum(grid3d)
            
        data={'grid3d':grid3d}
        data.update({'1dgrid': grid})
        data.update({'Origin': origin})
        data.update({'Cell': cell})
        data.update({'nuclcharge': sumposcharge})
        data.update({'positions': coords})
        data.update({'natoms': natoms})
        data.update({'xyzlist':xyz})
        data.update({'gdims':gdims})
        data.update({'spacemesh': spacemesh})
        data.update({'diffvol':dv})
        data.update({'diffxyz':mesh})
        data.update({'gridvector':gridvect})
        data.update({'headlines':heads})
        data.update({'datasum':datasum})

        
        
        return data

    

        
    @property     
    def natoms(self):
        return self.readcube(self.cube)['natoms']

    @property     
    def diffvol(self):
        return self.readcube(self.cube)['diffvol']
        
        
    @property
    def origin(self):
        return self.readcube(self.cube)['Origin']
    
    @property
    def grid3d(self):
        return self.readcube(self.cube)['grid3d']
        
    @property
    def cell(self):
        return self.readcube(self.cube)['Cell']


    @property
    def headlines(self):
        return self.readcube(self.cube)['headlines']
    
    @property
    def position(self):
        return self.readcube(self.cube)['positions']   
    
    @property
    def xyz(self):
        return self.readcube(self.cube)['xyzlist'] 

    @property
    def gdims(self):
        return self.readcube(self.cube)['gdims'] 

    @property
    def spacemesh(self):
        return self.readcube(self.cube)['spacemesh'] 
    
    @property
    def diffxyz(self):
        return self.readcube(self.cube)['diffxyz']    

    @property
    def grid1d(self):
        return self.readcube(self.cube)['1dgrid']    
    @property
    def datasum(self):
        return self.readcube(self.cube)['datasum']   
   
    @property
    def gridvectors(self):
        return self.readcube(self.cube)['gridvector']    
           
    

    def writecube(headlines, data, N1, N2, N3):
        
        
        """
        
        Parameters
        ----------
        headlines : list
            first natoms+6 lines of cube file imported from other cube 
        data : ndarray
            1d array of data with (n1n2n3) shape.
        N1 : int
            number of mesh grid in x direction.
        N2 : int
            number of mesh grid in y direction
        N3 : int
            number of mesh grid in z direction.

        Returns
        -------
        writes a cube file

        """
        os.remove("res.cube")
        with open("res.cube", "a") as outfile:
            for line in headlines:
                outfile.write(line)
            rownumb=int(N3/6)
            b=N3%6
            if b > 0:
                rownumb=rownumb+1
            
        
            
            print("row numb is ", str(rownumb))    
        
            for j in range(N1*N2):
                blocks=data[N3*j:N3*(j+1)]
                for i in range(rownumb):
                    a=blocks[6*i:(i+1)*6]
                    np.savetxt(outfile, a.T, fmt='  %-.5E', delimiter='')

            outfile.close()
        
        return
    

    @staticmethod
    def sumcubes(file1, file2):
        
        
        """
        Parameters
        ----------
        file1 : cubefile
            first cube.
        file2 : cubefile
            second cube.

        Returns
        -------
        sum of two cubes

        """
        
        
        cub1=cube1(file1)
        cub2=cube1(file2)
        
        
        head1=cub1.headlines
        head2=cub2.headlines
        
        
        if head1 != head2:
            raise Exception("Sorry, try similar cube files (same atom number)")

        
        data1=cub1.grid1d
        data2=cub2.grid1d
        
        
        
        sumdata=data1+data2
        N1, N2, N3= cub1.gdims[0], cub1.gdims[1], cub1.gdims[2]
        
        
        cube1.writecube(head1, sumdata, N1, N2, N3 )
        
        
        return
    
    
    
    
    
    
    @staticmethod
    def subtractcubes(file1, file2):
        
        
        """
        Parameters
        ----------
        file1 : cubefile
            first cube file.
        file2 : cubefile
            second cube file.

        Returns
        -------
        sum of two cubes

        """
        
        
        cub1=cube1(file1)
        cub2=cube1(file2)
        
        
        head1=cub1.headlines
        head2=cub2.headlines
        
        
        if head1 != head2:
            raise Exception("Sorry, try similar cube files (same atom number)")

        
        data1=cub1.grid1d
        data2=cub2.grid1d
        
        
        
        diffdata=data1-data2
        N1, N2, N3= cub1.gdims[0], cub1.gdims[1], cub1.gdims[2]
        
        
        cube1.writecube(head1, diffdata, N1, N2, N3 )
        
        
        return
        



    @staticmethod
    def multiplycubes(file1, file2):
        
        
        """
        Parameters
        ----------
        file1 : cubefile
            first cube file.
        file2 : cubefile
            second cube file.

        Returns
        -------
        sum of two cubes

        """
        
        
        cub1=cube1(file1)
        cub2=cube1(file2)
        
        
        head1=cub1.headlines
        head2=cub2.headlines
        
        
        if head1 != head2:
            raise Exception("Sorry, try similar cube files (same atom number)")

        
        data1=cub1.grid1d
        data2=cub2.grid1d
        
        
        
        mdata=data1*data2
        N1, N2, N3= cub1.gdims[0], cub1.gdims[1], cub1.gdims[2]
        
        
        cube1.writecube(head1, mdata, N1, N2, N3 )
        
        
        return
        






    @staticmethod
    def squarecubes(file1):
        
        
        """
        Parameters
        ----------
        file1 : cubefile
            name of cube file.
        
        Returns
        -------
        sum of two cubes

        """
        
        
        cub1=cube1(file1)
        
            
        head1=cub1.headlines
        
        
               
        data1=cub1.grid1d
        
        
        
        
        sqdata=data1*data1
        N1, N2, N3= cub1.gdims[0], cub1.gdims[1], cub1.gdims[2]
        
        
        cube1.writecube(head1, sqdata, N1, N2, N3 )
        
        
        return
        





        
        
        
        
        
        
        
        

    def plotxy(self, slicenumb):
        plt.xticks([])
        plt.yticks([])
        data2d=self.grid3d[:,:,slicenumb]
        plt.imshow(data2d, interpolation="gaussian")
        plt.savefig("xy_slice#"+str(slicenumb)+".png", dpi=300)
        
    
    def plotxz(self, slicenumb):
        plt.xticks([])
        plt.yticks([])
        data2d=self.grid3d[:,slicenumb,:]
        plt.imshow(data2d, interpolation="gaussian")
        plt.savefig("xz_slice#"+str(slicenumb)+".png", dpi=300)
        
                    
    def plotyz(self, slicenumb):
        """
        plots this data at given slice
        """
        plt.xticks([])
        plt.yticks([])
        data2d=self.grid3d[slicenumb,:,:]
        plt.imshow(data2d, interpolation="gaussian")
        plt.savefig("yz_slice#"+str(slicenumb)+".png", dpi=300)
        
    @property
    def xyplanaravg(self):
        area0=np.cross(self.cell[0,:], self.cell[1,:])
        area=np.linalg.norm(area0)
        xyavg=np.zeros((self.gdims[2],1))
        z=np.zeros((self.gdims[2],1))
        for i in range(self.gdims[2]):
            z[i]=i*np.linalg.norm(self.cell[2,:])/self.gdims[2]
            xyavg[i]=(1.0/area)*np.sum(self.grid3d[:,:,i])*self.diffxyz[1]*self.diffxyz[2]
        xydata=np.hstack((z,xyavg))
        plt.plot(z, xyavg, lw=3.0)
        return xydata
        
    
    @property
    def xzplanaravg(self):
        area0=np.cross(self.cell[0,:], self.cell[2,:])
        area=np.linalg.norm(area0)
        xzavg=np.zeros((self.gdims[1],1))
        y=np.zeros((self.gdims[1],1))
        for i in range(self.gdims[1]):
            y[i]=i*np.linalg.norm(self.cell[1,:])/self.gdims[1]
            xzavg[i]=(1.0/area)*np.sum(self.grid3d[:,i,:])*self.diffxyz[0]*self.diffxyz[2]
        xzdata=np.hstack((y,xzavg))
        plt.plot(y, xzavg, lw=3.0)
        return xzdata

    @property
    def yzplanaravg(self):
        area0=np.cross(self.cell[1,:], self.cell[2,:])
        area=np.linalg.norm(area0)
        yzavg=np.zeros((self.gdims[0],1))
        x=np.zeros((self.gdims[0],1))
        for i in range(self.gdims[0]):
            x[i]=i*np.linalg.norm(self.cell[0,:])/self.gdims[0]
            yzavg[i]=(1.0/area)*np.sum(self.grid3d[i,:,:])*self.diffxyz[1]*self.diffxyz[2]
        yzdata=np.hstack((x,yzavg))
        plt.plot(x, yzavg, lw=3.0)
        return yzdata          
        
####################  CHECK CHARGE
##################new=test.xzplanaravg
        ###### 
        ##########  np.sum(new[:,1])*new[1,0]*test.cell[0,0]*test.cell[2,2]
    
        
#############################################################
#################    CHARACTERS    ##########################


    

####    return mesh of X,Y,Z from all points in grid
    @property    
    def xyzmesh(self):
        N1, N2, N3 = self.gdims[0], self.gdims[1], self.gdims[2]
        a1, a2, a3= self.cell[0][0], self.cell[1][1], self.cell[2][2]
        x = np.linspace(0, a1, N1+1)[:-1]
        y = np.linspace(0, a2, N2+1)[:-1]
        z = np.linspace(0, a3, N3+1)[:-1]
        xyzmesh=np.meshgrid(x,y,z, indexing='ij')
        return xyzmesh
                    
    
####    return distance of all points from a refrence
    def distfromref(self,xx,yy,zz):
        X,Y,Z=self.xyzmesh
        rr=(X-xx)*(X-xx)+(Y-yy)*(Y-yy)+(Z-zz)*(Z-zz)
        rs=np.sqrt(rr)
        return rs
    
    
    def dist2from(self, xx,yy,zz):
        N1, N2, N3 = self.gdims[0], self.gdims[1], self.gdims[2]
        xe, ye, ze = self.diffxyz[0], self.diffxyz[1], self.diffxyz[2]
        dd=np.zeros((N1, N2, N3))
        for i in  range(N1):
            for j in range(N2):
                for k in range(N3):
                    x=(i*xe-xx)
                    y=(j*ye-yy)
                    z=(k*ze-zz)
                    dd[i,j,k]=np.sqrt(x**2+y**2+z**2)
        return dd


    @property
    def conv2ase(self):
        temp=[x[0] for x in self.xyz[2:]]
        symb0=''
        for name in temp:
            symb0=symb0+name
        pos1=np.array([x[1:] for x in self.position])
        cell1=self.cell
        cell1=cell1*Bohr
        mystr=Atoms(symb0)
        pos1=pos1*Bohr
        mystr.cell=cell1
        mystr.positions=pos1
        mystr.set_positions(pos1)
        mystr.set_cell(cell1)
 #       mystr.set_positions(pos1*Bohr)
        return mystr
    
    
    # returns all vectors in the translational group of special basis
    # and their norm in 4*N matrix with criteria norm < 25 Bohr!!!
    def supervectors(self, xatom,yatom,zatom):
        rpos=[[xatom,yatom,zatom, np.linalg.norm([xatom,yatom,zatom])]]
        a=self.cell[0]
        b=self.cell[1]
        c=self.cell[2]
        for i in range(-2,3):
            for j in range(-2,3):
                for k in range(-2,3):
                    temp=i*a+j*b+k*c+np.array([xatom,yatom,zatom])
                    temp0=np.append(temp, np.linalg.norm(temp))
                    temp0=np.array([temp0])
                    rpos=np.append(rpos, temp0, axis=0)
        filter1=rpos[np.argsort(rpos[:,3])]
        filter2=filter1[1:,:]
        filter3=filter2[filter2[:,3] < 25]
        return filter3
    
    @property
    def supercell(self):
        posit=self.position
        posit=np.float32(posit)
        cell=self.cell
        cell=np.float32(cell)
        supcell=posit
        for i in range(-3,4):
            for j in range(-3,4):
                for k in range(-3,4):
                    vector=i*cell[0]+j*cell[1]+k*cell[2]
                    vectmodified=np.hstack((np.array([0]), vector))
                    new=posit+vectmodified
                    if ( np.linalg.norm(vector) > 0.0 and np.linalg.norm(vector) < 34):
                        supcell=np.append(supcell, new, axis=0)   
        supcell=np.float32(supcell)
        return supcell
                    
                    
                    
        
        
    
    @property
    def cube2posdens(self):
        
        '''
         Parameters
        ----------
        
        Returns
        -------
        4 column data including xyz and density
        
        
        '''
        
        
      
        
        
        ggg=self.grid1d
        vvv=self.gridvectors
        denspos=np.hstack((vvv,ggg))
        return denspos
            
    
    
    
    def filterdenspos(self,threshold=0.0000002):
        ggg=self.grid1d
        vvv=self.gridvectors
        denspos=np.hstack((ggg,vvv))
        denspos1=denspos[np.argsort(denspos[:,0])]
        aa=np.array([denspos1[0]])
        for i in range(len(ggg)-1):
            if denspos1[i+1,0]-denspos1[i,0] > threshold:
                aa=np.append(aa, np.array([denspos1[i+1]]), axis=0)
        return aa
        
    
    
    
        
    @staticmethod
    def Skomegafunc(rnorm,Rc,sigmak):
        '''
        Parameters
        ----------
        rnorm : rnormvector
            distance between atom_i and grid point.
        Rc : float
            cutoff_raddi in Bohr.
        sigmak : gaussian width

        Returns
        -------
        scalar fuction Skomrgi for a given atom and grid distance.
        
        see paper J. Phys. Chem. C 2019, 123, 15859−15866.
        '''
        fcut=0.5*(np.cos(rnorm*np.pi/Rc)+1)
        g=np.exp(-np.power(rnorm, 2.) / (2 * np.power(sigmak, 2.))) 
        Sfunc0=g*fcut*(1./(sigmak*np.sqrt(np.pi)))**3
        Sfunc=np.sum(Sfunc0)
        return Sfunc
        
        

    def Skomeg1point(self, xg,yg,zg, sigmaval, Rc=18.0):
        '''
        
        Parameters
        ----------
        xg : float
            x of grid point.
        yg : float
            y of grid point.
        zg : float
            z of grid point.
        sigmaval : float
            width of gaussian.
        Rc : float, optional
            Width of gaussian. The default is 18.0.

        Returns
        -------
        Sfork : vector (ndarray)
            vector of S for each grid point.

        '''
        
        
        supercell=self.supercell
        sortedpositions=supercell[np.argsort(supercell[:,0])]
        vecmod=np.array([0.0,xg,yg,zg])
        
        sortpossub=sortedpositions-vecmod
        
        sorted2=sortpossub[sortpossub[:,1]**2+sortpossub[:,2]**2+sortpossub[:,3]**2 < Rc**2]
        normsort=np.sum(sorted2[:,1:3]**2, axis=-1)**(1./2)
        
        
        species_dict=collections.Counter(sorted2[:,0])
        
        typesymb=list(species_dict.keys())
        typefreq=list(species_dict.values())
                
        #ntype return number of atom types
        ntypes=len(typesymb)
        
        #define an array to store data
  
        
        
        
        #cut separate blocks of sortedpositionmat
        j=0
        Sfork=np.zeros((1,ntypes))
        for i in range(ntypes):
            nni=typefreq[i]
            matblock=normsort[j:j+nni]
            Sfork[0,i]=cube1.Skomegafunc(matblock, Rc, sigmaval)
            j=j+nni
        
              
        return Sfork
  
               
        
    def Skomeg1points(self, gridvector, sigmak, Rc=18.0):
        
        '''
        
        Parameters
        ----------
 
        gridvector : array
            nx3 array of grid point coordinates.
        sigmak : float
            width of gaussian.
        Rc : float, optional
            Width of gaussian. The default is 18.0.

        Returns
        -------
        Sfork : vector (ndarray)
            vector of S for each grid point.

        '''
        
        
        #extract number of points in grid
        
        gridnumb=np.shape(gridvector)[0]
        
        
        
        #
        
        
#       extract and sort supercell positions by atomic number
#        
        supercell=self.supercell
        sortedpositions=supercell[np.argsort(supercell[:,0])]

#       extract atomic number and frequency of each species  
        species_dict=collections.Counter(sortedpositions[:,0])
        
        
        
        typesymb=list(species_dict.keys())
        typefreq=list(species_dict.values())
        
        
        typeind=np.cumsum(typefreq)
        typeindex=[0, *typeind]
                
        #ntype return number of atom types
        ntypes=len(typesymb)
        
        # return number of rows in supercell
        supnum=np.shape(sortedpositions)[0]
        
        # reserve a mat for repeating all grid points over positions
        repeatovergrid=np.zeros((gridnumb*supnum,3), dtype=np.float32)
        
        Sfork=np.zeros((gridnumb,ntypes),dtype=np.float32)

        
        for i in range(gridnumb):
            sortpossub=sortedpositions[:,1:4]-gridvector[i]
            repeatovergrid[i*supnum:(i+1)*supnum, :]=sortpossub
        
       
        
       
       #  size of repeatedgrid is around 1.3 GB !!!!!!!!!!!! for 50 50 50 
       
        
        # right distance and vector and tensor here
        rnorm=np.sum(repeatovergrid[:,0:3]**2, axis=-1)**(1./2)
        cutfunc=np.heaviside(Rc-rnorm,0.0)
        
        normcoef=(1./(sigmak*2*np.sqrt(np.pi)))**3
        
        
        
        Scalar0=np.exp(-np.power(rnorm, 2.) / (2 * np.power(sigmak, 2.)))*0.5*(np.cos(rnorm*np.pi/Rc)+1)*normcoef*cutfunc
        
        #define an array to store data
  
        
  
        # Now we have a #gridnum blocls of repeated supcell# dim matrix
        # access to each block
        for p in range(gridnumb):
            for j in range(ntypes):
                Sfork[p,j]=np.sum(Scalar0[p*supnum+typeindex[j]:p*supnum+typeindex[j+1]])
                
                
                
                
        return Sfork    
                    
                    

                    
#!!!!!!!  #############This should be edited
        
    def Skomeg1molecule(self, gridvector, sigmak, Rc=18.0):
        
        '''
        
        Parameters
        ----------
 
        gridvector : array
            nx3 array of grid point coordinates.
        sigmak : float
            width of gaussian.
        Rc : float, optional
            Width of gaussian. The default is 18.0.

        Returns
        -------
        Sfork : vector (ndarray)
            vector of S for each grid point.

        '''
        
        
        #extract number of points in grid
        
        gridnumb=np.shape(gridvector)[0]
        
        
        
        #
        
        
#       extract and sort supercell positions by atomic number
#        
        supercell=self.supercell
        sortedpositions=supercell[np.argsort(supercell[:,0])]

#       extract atomic number and frequency of each species  
        species_dict=collections.Counter(sortedpositions[:,0])
        
        
        
        typesymb=list(species_dict.keys())
        typefreq=list(species_dict.values())
        
        
        typeind=np.cumsum(typefreq)
        typeindex=[0, *typeind]
                
        #ntype return number of atom types
        ntypes=len(typesymb)
        
        # return number of rows in supercell
        supnum=np.shape(sortedpositions)[0]
        
        # reserve a mat for repeating all grid points over positions
        repeatovergrid=np.zeros((gridnumb*supnum,3), dtype=np.float32)
        
        Sfork=np.zeros((gridnumb,ntypes),dtype=np.float32)

        
        for i in range(gridnumb):
            sortpossub=sortedpositions[:,1:4]-gridvector[i]
            repeatovergrid[i*supnum:(i+1)*supnum, :]=sortpossub
        
       
        
       
       #  size of repeatedgrid is around 1.3 GB !!!!!!!!!!!! for 50 50 50 
       
        
        # right distance and vector and tensor here
        rnorm=np.sum(repeatovergrid[:,0:3]**2, axis=-1)**(1./2)
        cutfunc=np.heaviside(Rc-rnorm,0.0)
        
        normcoef=(1./(sigmak*2*np.sqrt(np.pi)))**3
        
        
        
        Scalar0=np.exp(-np.power(rnorm, 2.) / (2 * np.power(sigmak, 2.)))*0.5*(np.cos(rnorm*np.pi/Rc)+1)*normcoef*cutfunc
        
        #define an array to store data
  
        
  
        # Now we have a #gridnum blocls of repeated supcell# dim matrix
        # access to each block
        for p in range(gridnumb):
            for j in range(ntypes):
                Sfork[p,j]=np.sum(Scalar0[p*supnum+typeindex[j]:p*supnum+typeindex[j+1]])
                
                
                
                
        return Sfork    
                    
                    







                    
    @staticmethod
    def vec2grid(vector, N0, N1, N2):
        """
        Parameters
        ----------
        vector : ndarray
            1d array
        N0 : int
            DESCRIPTION.
        N1 : int
            DESCRIPTION.
        N2 : int
            DESCRIPTION.

        Returns
        -------
        grid : ndarray
            3d array of converted vector.

        """
        
        grid=np.zeros((N0,N1,N2))
        for i in range(N0):
            for j in range(N1):
                for k in range(N2):
                    grid[i,j,k]=vector[i*N1*N2+j*N2+k]    
        return grid
    
    
    
    @staticmethod
    def createcubegrid(N0, N1, N2, a, b, c):
        dx=a/N0
        dy=b/N1
        dz=c/N2
        grid=np.zeros((3, N0, N1, N2))
        for i in range(N0):
            for j in range(N1):
                for k in range(N2):
                    grid[:,i,j,k]=[i*dx, j*dy, k*dz]
                    
                    

        

    @staticmethod
    def createvector(N0,N1,N2, a, b, c):
        myvector0=np.zeros((N0*N1*N2,3))
        dx=a/N0
        dy=b/N1
        dz=c/N2
        for i in range(N0): 
            for j in range(N1): 
                for k in range(N2): 
                    myvector0[i*N1*N2+j*N2+k]=np.array([i*dx,j*dy,k*dz]) 
        return myvector0
        
    
    
    
##  When do you need float16 (4 decimal), float32 or float64 (default)
       
            
    
 ############  Draft codes
# How to return vector of meshgrid  ????????????????????

##  X,Y,Z=np.meshgrid(..............)
##a=np.array([X,Y,Z]) 
 
#  In [25]: a[:,3,3,2]                                                                     
# Out[25]: array([ 1.69230769,  3.61538462, -2.69230769])

# In [26]: vect=np.zeros((64000,3))                                                       

# In [27]: for i in range (40): 
#     ...:     for j in range(40): 
#     ...:         for k in range(40): 
#     ...:             vect[i*1600+j*40+k]=a[:,i,j,k] 
#     ...:                                                                                

# In [28]: vect[63000]                                                                    
# Out[28]: array([ 4.46153846, 11.        , -3.        ])

   




#-------------------------------------------
    
        
#        Ten Gaussians with width varying from 0.75 to 8 Å on a logarithmic scale
#        kvalues=np.logspace(np.log10(0.75/Bohr),np.log10(8/Bohr),10)
        
        
        
        
        
    
    
    
        
#--------make vector of [x,y,z, n          ]  with ngridx4
        
        
                    
                    


                
        
########  CHARGE partitioning with ML techniques  k_means ##############

    @staticmethod
    def gmm_partition(densityfile):
        
        basename=densityfile.split(".")[0]
        denscube=cube1(densityfile)
        
        
        natoms=denscube.natoms
        rhovector=denscube.cube2posdens
        
        
        threshold1=0.01
         

#Enter filter

        rhodata=rhovector[rhovector[:,3] > threshold1]  
        
        
        #   density data
        rho=rhodata[:,3]
        
        
        #grid xyz
        coords=rhodata[:,:3]
        
        
       
        
        rho1d=np.reshape(rho, (np.shape(rho)[0],))
        
        
        
        
        diffv=denscube.diffvol
                    
        sum1=np.sum(rho)*diffv
        
        sumatomicnumbers=sum(denscube.position[:,0])
        
        lerr=(sum1-sumatomicnumbers)/sumatomicnumbers
        renorm=sum1/sumatomicnumbers
        
        rrho=rho/renorm
        
        
        
        
        print("diffv is   ", diffv)
        print("charge leakage is  ", lerr*100 ," %" )
        if abs(lerr) > 0.1:
            print("there is big difference between cube and atomic charge")
            return  # 'exit' function and return to caller
        
        partialcharges=np.zeros((natoms,2))
        
        
        
        xyzvector=rhodata[:,:3]
#        weight1=((elfvector[:,3]))
        
 #       weightelf=1/weight1
 #       weight=weightelf/max(weightelf)
        
        atomcenters=denscube.position[:,1:]-denscube.origin
        atomicnumbers=denscube.position[:,0]
        #
        atomicweights=atomicnumbers/np.sum(atomicnumbers)
        
        
        
        
        #calculate weighted by elf
  #      func1=KMeans(n_clusters=natoms,  init=atomcenters,tol=0.00001,
    #                 verbose=1, max_iter=10000)
  #      func1=DBSCAN()
        func1=GaussianMixture(n_components=natoms, max_iter=100, n_init=1, 
                              weights_init=atomicweights, means_init=atomcenters,  verbose=1, verbose_interval=2)
  
  #      func1.fit(xyzvector, sample_weight=weight)
        labels=func1.fit_predict(xyzvector)
        

        
        labels01=np.reshape(labels, (np.shape(xyzvector)[0],1) )
        xyzlabeled=np.hstack((xyzvector,labels01 ))
            
        partialcharges[:,0]=denscube.position[:,0]
        for i in range(natoms):
            partialcharges[i,1]=np.sum(rrho[labels==i])
#            partialcharges[i,2]=np.sum(rho[labelsrho==i])
            
        
        partialcharges[:,1]=atomicnumbers-(partialcharges[:,1]*diffv)
        
        
        
        
        np.savetxt(basename+".gmm", partialcharges, fmt="%5d %8.4f")
        
###  uncomment for more details        
#        partialcharges[:,2]=partialcharges[:,2]*diffv
 
        
#        np.savetxt( "data.txt", xyzlabeled, fmt="%8.4f")
        
#        fig = plt.figure(figsize=(8, 6))
#        ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=48, azim=134)
        
        
#        ax.scatter(coords[:,0],coords[:,1], coords[:,2], 
#               c=labels.astype(float), alpha=0.7)
        
        
        
        
        return partialcharges
            
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        


        
    
    