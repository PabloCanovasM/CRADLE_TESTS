import subprocess
import numpy as np


ccNames = ["cs","csp","ct","ctp","cv","cvp","ca","cap"]
ccTag = ["--CS","--CSP","--CT","--CTP","--CV","--CVP","--CA","--CAP"]
real_interest = [(0,1),(2,3),(4,5),(6,7),(0,3),(0,6),(0,7),(1,2),(1,6),(1,7),(2,4),(2,5),(2,7),(3,4),(3,5),(3,6),(4,7),(5,6)]
complex_interest = [(0,1),(2,3),(4,5),(6,7),(0,2),(0,7),(1,3),(1,6),(2,5),(2,7),(3,4),(3,6),(4,6),(5,7)]

#account for the value of M_GT, makig it so that M_GT*C(')_T and M_GT*C(')_A have norm 1
correction = 1/0.6606
val_corr = [2,3,6,7]

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
            couplingConstants[val_corr] *= correction
            print(couplingConstants[[i,j]])
            run_array = ["./CRADLE++", "nucleus", "-n", "39Ca", "-Z", "20", "-A", "39", "general", "-l", "1000000", "-t", "20", "-v", "1", "-o", f"build/39Ca_{ccNames[i]}{ccNames[j]}_real.txt", "betaDecay","-E", "0", "-r", "0", "--PolarisationZ", "1", "--PolarisationMag", str(1), "--Alignment", str(1), "coupling"]
            for name, val in zip(ccTag,couplingConstants):
                run_array.append(name)
                run_array.append(str(val))
            print(ccNames[i],ccNames[j],"real")
            subprocess.run(run_array)
        interesting = False
        for pair in complex_interest:
            if (pair[0] == i and pair[1] == j):
                interesting = True
        if interesting:
            couplingConstants = np.zeros(8).astype(np.cdouble)
            couplingConstants[i] = 0.8+0.6j
            couplingConstants[j] = 0.6-0.8j
            couplingConstants[val_corr] *= correction
            print(couplingConstants[[i,j]])
            run_array = ["./CRADLE++", "nucleus", "-n", "39Ca", "-Z", "20", "-A", "39", "general", "-l", "1000000", "-t", "20", "-v", "1", "-o", f"build/39Ca_{ccNames[i]}{ccNames[j]}_imag.txt", "betaDecay","-E", "0", "-r", "0", "--PolarisationZ", "1", "--PolarisationMag", str(1), "--Alignment", str(1), "coupling"]
            for name, val in zip(ccTag,couplingConstants):
                run_array.append(name)
                if val.imag == 0:
                    run_array.append(str(val.real)) 
                else:
                    sign = "-" if val.imag < 0 else "+"
                    run_array.append(str(val.real)+sign+str(abs(val.imag))+"j")        
            print(ccNames[i],ccNames[j],"imag")
            subprocess.run(run_array)
        
