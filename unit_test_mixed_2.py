import subprocess
import numpy as np


ccNames = ["cs","csp","ct","ctp","cv","cvp","ca","cap"]
ccTag = ["--CS","--CSP","--CT","--CTP","--CV","--CVP","--CA","--CAP"]
real_interest = [(0,1),(2,3),(4,5),(6,7),(0,3),(0,6),(0,7),(1,2),(1,6),(1,7),(2,4),(2,5),(2,7),(3,4),(3,5),(3,6),(4,7),(5,6)]
complex_interest = [(0,1),(2,3),(4,5),(6,7),(0,2),(0,7),(1,3),(1,6),(2,5),(2,7),(3,4),(3,6),(4,6),(5,7)]

Q = 5502.5 #kev
j_i = 1.5
Z = -20

for i in range(8):
    for j in range(i+1,8):
        interesting = False
        for pair in real_interest:
            if (pair[0] == i and pair[1] == j):
                interesting = True
        if interesting:
            couplingConstants = np.zeros(8)
            couplingConstants[i] = 1
            couplingConstants[j] = 1
            run_array = ["build/Test", "minmax_f", "-Q", str(Q), "-M", str(1), str(1), "-Z", str(Z), "-J", str(j_i), str(j_i), "-o", f"build/39Ca_{ccNames[i]}{ccNames[j]}_minmaxf_real.txt", "-C"]
            for val in couplingConstants:
                run_array.append(str(val))
                run_array.append(str(0))
            print(ccNames[i],ccNames[j],"real")
            subprocess.run(run_array)
            print(ccNames[i],ccNames[j],"real: completed")
        interesting = False
        for pair in complex_interest:
            if (pair[0] == i and pair[1] == j):
                interesting = True
        if interesting:
            couplingConstants = np.zeros(8).astype(np.cdouble)
            couplingConstants[i] = 0.8+0.6j
            couplingConstants[j] = 0.6-0.8j
            run_array = ["build/Test", "minmax_f", "-Q", str(Q), "-M", str(1), str(1), "-Z", str(Z), "-J", str(j_i), str(j_i), "-o", f"build/39Ca_{ccNames[i]}{ccNames[j]}_minmaxf_imag.txt", "-C"]
            for val in couplingConstants:
                run_array.append(str(val.real))
                run_array.append(str(val.imag))        
            print(ccNames[i],ccNames[j],"imag")
            subprocess.run(run_array)
            print(ccNames[i],ccNames[j],"imag: completed")
