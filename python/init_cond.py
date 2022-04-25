import para
import numpy as np
import copy


def initcond(G):
    if para.dimension==2:
        '''
        Initial wavefunction
        '''
        G.psi = ((64/np.pi)**(1/2))*np.exp(-64*(para.z_mesh**2/2+para.x_mesh**2/2))+0j    #2D initial condition
        G.psi=np.array(G.psi,dtype=np.complex128)
        G.temp=copy.deepcopy(G.psi)
        

        '''
        potential
        '''
        G.V=64*(para.x_mesh**2+para.z_mesh**2)/2
        del para.x_mesh,para.z_mesh


    elif para.dimension==3:
        '''
        Initial wavefunction
        '''
        G.psi = ((64/np.pi)**(3/4))*np.exp(-64*(para.z_mesh**2+para.y_mesh**2+para.x_mesh**2)/2)+0j  # 3D initial condition
        G.psi=np.array(G.psi,dtype=np.complex128)
        G.temp=copy.deepcopy(G.psi)


        '''
        potential
        '''
        G.V=64*(para.x_mesh**2+para.y_mesh**2+para.z_mesh**2)/2
        del para.x_mesh,para.y_mesh,para.z_mesh
