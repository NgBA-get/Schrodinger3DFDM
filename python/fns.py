import para
from scipy.integrate import simpson as sn

def my_integral3D(psi):
    return sn(sn(sn(psi,para.z),para.y),para.x)
def my_integral2D(psi):
    return sn(sn(psi,para.z),para.x)