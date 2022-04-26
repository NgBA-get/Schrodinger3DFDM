import numpy as np
'''
For 2D set Nx,Nz for grids and Lx, Lz for box size along the respective direction 
'''
# No of grids points along the different direction 
N1 = 64
Nx = N1
Ny = N1
Nz = N1

# Box size along the different directions 
Lx=1
Ly=1
Lz=1

dimension=3 #set dimension=2 for 2d and dimension=3 for 3d

# Non linearity for the GPE
Nlin=0

# Time at which program will calculte Norm, Energy and Chemical potential
tstep=np.array([0.001,1.0,4.0,6.0,10.0]) 

# Time upto which we have to run the code
tMax=5

dt=0.0005

# No of time steps required to reach tMax 
Ngt=int(tMax/dt)

x=np.arange(-Nx//2+1,Nx//2+1)*Lx/Nx
dx=x[1]-x[0]

if dimension==3:
    y=np.arange(-Ny//2+1,Ny//2+1)*Ly/Ny
    dy=y[1]-y[0]

z=np.arange(-Nz//2+1,Nz//2+1)*Lz/Nz
dz=z[1]-z[0] 


# Meshgrid 
if dimension ==2:
    x_mesh,z_mesh=np.meshgrid(x,z,indexing='ij')

elif dimension==3:
    x_mesh,y_mesh,z_mesh=np.meshgrid(x,y,z,indexing='ij')
