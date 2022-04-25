import numpy as np
import para
import copy

if para.dimension==2:
    from derivative import laplacian2D as laplacian 
    from fns import my_integral2D as my_integral
elif para.dimension==3:
    from derivative import laplacian3D as laplacian 
    from fns import my_integral3D as my_integral

class Wavefunction:
    def __init__(self):
        self.psi=[]
        self.V=[]
        self.temp=[] 
    
    def set_arrays(self):
        if para.dimension==2:
            self.psi=np.zeros((para.Nx,para.Nz),dtype=np.complex128)
            self.temp=np.zeros((para.Nx,para.Nz),dtype=np.complex128)
            self.V=np.zeros((para.Nx,para.Nz),dtype=np.float64)
        elif para.dimension==3:
            self.psi=np.zeros((para.Nx,para.Ny,para.Nz),dtype=np.complex128)
            self.temp=np.zeros((para.Nx,para.Ny,para.Nz),dtype=np.complex128)
            self.V=np.zeros((para.Nx,para.Ny,para.Nz),dtype=np.float64)
        
    def energy(self):
        z=np.conjugate(self.psi)*(-1/128*laplacian(self.psi)+(self.V+para.Nlin*np.abs(self.psi)**2/2)*self.psi)
        return my_integral(z.real)
    

    def norm(self):
        return my_integral(np.abs(self.psi)**2)
    
     
    def chempot(self):  # When Non-lineartity is present
        z=np.conjugate(self.psi)*(-1/128*laplacian(self.psi)+(self.V+para.Nlin*np.abs(self.psi)**2)*self.psi)
        return my_integral(z.real)
    

    def compute_RHS(self):
        self.psi=-1j*(-1/128*laplacian(self.psi)+(self.V+para.Nlin*np.abs(self.psi)**2)*self.psi)
    
    # Eular scheme for single step time evolution
    def sstep_eular(self):
        self.temp=copy.deepcopy(self.psi)
        self.compute_RHS()
        self.psi=self.temp+self.psi*para.dt

    # RK2 schemes for single step time evolution
    def sstep_RK2(self):
        self.temp=copy.deepcopy(self.psi)
        self.compute_RHS()
        self.psi=self.temp+self.psi*para.dt/2
        self.compute_RHS()
        self.psi=self.temp+para.dt*self.psi