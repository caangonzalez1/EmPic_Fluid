import numpy as np                                                  
import struct                                                       
import matplotlib.pyplot as plt                                     
from matplotlib import cm
from matplotlib.gridspec import GridSpec
import glob
import sys
 

#read particles from processor 0
files = sorted(glob.glob('bx_n*'))
nout = len(files)
bx_n = np.zeros([Nx+1,Ny+1,nout])
for i in range(nout):
      print(files[i])
      bz0 =  np.loadtxt(files[i])
      tmp = np.reshape(bz0,[Nx+1,Ny+1],order = 'C')
      bx_n[:,:,i] = tmp
files = sorted(glob.glob('by_n*'))
nout = len(files)
by_n = np.zeros([Nx+1,Ny+1,nout])
for i in range(nout):
      print(files[i])
      bz0 =  np.loadtxt(files[i])
      tmp = np.reshape(bz0,[Nx+1,Ny+1],order = 'C')
      by_n[:,:,i] = tmp

files = sorted(glob.glob('bz_n*'))
nout = len(files)
bz_n = np.zeros([Nx+1,Ny+1,nout])
for i in range(nout):
      print(files[i])
      bz0 =  np.loadtxt(files[i])
      tmp = np.reshape(bz0,[Nx+1,Ny+1],order = 'C')
      bz_n[:,:,i] = tmp


files = sorted(glob.glob('ex_n*'))
nout = len(files)
ex_n = np.zeros([Nx+1,Ny+1,nout])
for i in range(nout):
      print(files[i])
      bz0 =  np.loadtxt(files[i])
      tmp = np.reshape(bz0,[Nx+1,Ny+1],order = 'C')
      ex_n[:,:,i] = tmp
files = sorted(glob.glob('ey_n*'))
nout = len(files)
ey_n = np.zeros([Nx+1,Ny+1,nout])

for i in range(nout):
      print(files[i])
      bz0 =  np.loadtxt(files[i])
      tmp = np.reshape(bz0,[Nx+1,Ny+1],order = 'C')
      ey_n[:,:,i] = tmp
files = sorted(glob.glob('ez_n*'))
nout = len(files)
ez_n = np.zeros([Nx+1,Ny+1,nout])
for i in range(nout):
      print(files[i])
      bz0 =  np.loadtxt(files[i])
      tmp = np.reshape(bz0,[Nx+1,Ny+1],order = 'C')
      ez_n[:,:,i] = tmp
      
      
    
files = sorted(glob.glob('jx_n*'))
nout = len(files)
jx_n = np.zeros([Nx+1,Ny+1,nout])
for i in range(nout):
      print(files[i])
      bz0 =  np.loadtxt(files[i])
      tmp = np.reshape(bz0,[Nx+1,Ny+1],order = 'C')
      jx_n[:,:,i] = tmp
files = sorted(glob.glob('jy_n*'))
nout = len(files)
jy_n = np.zeros([Nx+1,Ny+1,nout])
for i in range(nout):
      print(files[i])
      bz0 =  np.loadtxt(files[i])
      tmp = np.reshape(bz0,[Nx+1,Ny+1],order = 'C')
      jy_n[:,:,i] = tmp
files = sorted(glob.glob('jz_n*'))
nout = len(files)
jz_n = np.zeros([Nx+1,Ny+1,nout])
for i in range(nout):
      print(files[i])
      bz0 =  np.loadtxt(files[i])
      tmp = np.reshape(bz0,[Nx+1,Ny+1],order = 'C')
      jz_n[:,:,i] = tmp
