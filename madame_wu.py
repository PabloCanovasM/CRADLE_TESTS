import numpy as np
import subprocess

mu_60Co = 3.799*3.152451e-8 #ev.T^-1
kb = 8.617e-5 #ev.K^-1
B = 1 #T

T = np.linspace(0.001,0.01,10)

mj_states = np.arange(-5,6)

p_states = np.exp(mu_60Co*B*mj_states[:,np.newaxis]/5/kb/T[np.newaxis,:])
#print(p_states.shape)
Z = p_states.sum(axis=0)
#print(Z.shape)
p_states /= Z[np.newaxis,:]
#for i in range(10):
#    print(p_states[:,i])


polz_list = np.sum(p_states*mj_states[:,np.newaxis]/5,axis=0)
Jzz = np.sum(p_states*mj_states[:,np.newaxis]**2,axis=0)
align_list = (3*Jzz-5*6)/5/9

polz_list = np.insert(polz_list,0,1)
align_list = np.insert(align_list,0,1)

print(polz_list)
print(align_list)

for i, polz, align in zip(range(polz_list.size),polz_list,align_list):
    print(i, polz, align)
    subprocess.run(["./CRADLE++", "nucleus", "-n", "60Co", "-Z", "27", "-A", "60", "general", "-l", "200000", "-o", f"build/60Co_polZWu{i}.txt", "-t", "10", "betaDecay", "-r", "0", "--PolarisationZ", "1", "--PolarisationMag", str(polz), "--Alignment", str(align)]) 
