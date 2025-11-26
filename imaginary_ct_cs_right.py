import subprocess
import numpy as np


l_list = np.array([4e5]).astype(int) #sizes

cs_list = np.linspace(0.4,0.4,1)
ct_list = np.linspace(0.4,0.4,1)
ccTag = ["--CS","--CSP","--CT","--CTP","--CV","--CVP","--CA","--CAP"]

for l in l_list:
    for cs in cs_list:
        for ct in ct_list:
            couplingConstants = np.zeros(8).astype(np.cdouble)
            couplingConstants[[4,5]] = 1
            couplingConstants[[6,7]] = 1
            couplingConstants[:4] = [cs*1.j, -cs*1.j, ct, -ct]
            print(couplingConstants)
            couplingConstans_array = []
            for name, val in zip(ccTag,couplingConstants):
                couplingConstans_array.append(name)
                if val.imag == 0:
                    couplingConstans_array.append(str(val.real)) 
                else:
                    sign = "-" if val.imag < 0 else "+"
                    couplingConstans_array.append(str(val.real)+sign+str(abs(val.imag))+"j")     
            run_array_39Ca = ["./CRADLE++", "nucleus", "-n", "39Ca", "-Z", "20", "-A", "39", "general", "-l", str(l), "-t", "20", "-v", "1", "-o", f"build/39Ca_cs_{cs:0.2f}jct_{ct:0.2f}_R_2.txt", "betaDecay","-E", "0", "-r", "0", "--PolarisationZ", "1", "--PolarisationMag", str(1), "--Alignment", str(1), "coupling"]
            run_array_39Ca.extend(couplingConstans_array)
            subprocess.run(run_array_39Ca)
            run_array_23Mg = ["./CRADLE++", "nucleus", "-n", "23Mg", "-Z", "12", "-A", "23", "general", "-l", str(l), "-t", "20", "-v", "1", "-o", f"build/23Mg_cs_{cs:0.2f}jct_{ct:0.2f}_R.txt", "betaDecay","-E", "0", "-r", "0", "--PolarisationZ", "1", "--PolarisationMag", str(1), "--Alignment", str(1), "coupling"]
            run_array_23Mg.extend(couplingConstans_array)
            subprocess.run(run_array_23Mg)
