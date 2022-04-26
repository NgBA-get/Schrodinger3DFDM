from sch import Wavefunction
import numpy as np
import para 
import init_cond
import time

G=Wavefunction()
G.set_arrays()
init_cond.initcond(G)
print('norm: ',G.norm())
t=para.dt
j=0

start = time.time()
for i in range(para.Ngt):
    
    G.sstep_RK2()

    if (para.tstep[j]-t)/para.dt<para.dt:
        print('\n----------------------')
        print('t: ',t)
        print('norm: ',G.norm())
        print('energy: ',G.energy())
        print('chemical pot: ',G.chempot())
        j=j+1
    t =t+para.dt
print(time.time()-start)
