import numpy as np
import para
import copy


def laplacian3D(psi):
    temp=copy.deepcopy(psi) 
    temp[1:-1,:,:]= (psi[2:,:,:]-2*psi[1:-1,:,:]+psi[:-2,:,:])/para.dx**2
    temp[:,1:-1,:]+=(psi[:,2:,:]-2*psi[:,1:-1,:]+psi[:,:-2,:])/para.dy**2
    temp[:,:,1:-1]+=(psi[:,:,2:]-2*psi[:,:,1:-1]+psi[:,:,:-2])/para.dz**2
    return temp

def laplacian2D(psi):
    temp=copy.deepcopy(psi) 
    temp[1:-1,:]= (psi[2:,:]-2*psi[1:-1,:]+psi[:-2,:])/para.dx**2
    temp[:,1:-1]+=(psi[:,2:]-2*psi[:,1:-1]+psi[:,:-2])/para.dz**2
    return temp
'''
def derivative(psi,dx):     ## Richardson extrapolation
    temp=(np.roll(psi,2)-8*np.roll(psi,1)+8*np.roll(psi,-1)-np.roll(psi,-2))/(12*dx)
    temp[0]=0
    temp[1]=(psi[2]-psi[0])/(2*dx)
    temp[-1]=0
    return temp
'''