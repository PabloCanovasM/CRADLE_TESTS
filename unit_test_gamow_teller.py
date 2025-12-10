import subprocess
import numpy as np

ccTag = ["--CT","--CTP"]
fileNames = ["","_ctposRe","_ctposIm","_ctnegRe","_ctnegIm","_ctoppRe","_ctoppIm"]
#account for the value of M_GT, makig it so that M_GT*C(')_T and M_GT*C(')_A have norm 1

couplingConstants = np.zeros((2,7),np.cdouble)
couplingConstants[:,1] = np.ones(2)*1.2754/np.sqrt(2)
for i in range(2,7):
    couplingConstants[:,i] = couplingConstants[:,i-1]*1.0j
couplingConstants[1,5:] *= -1
print(couplingConstants)
for i, fileName in enumerate(fileNames):
    run_array = ["./CRADLE++", "nucleus", "-n", "60Co", "-Z", "27", "-A", "60", "general", "-l", "500000", "-t", "20", "-v", "1", "-o", f"build/60Co_polZpos{fileName}.txt", "betaDecay","-E", "0", "-r", "0", "--PolarisationZ", "1", "--PolarisationMag", str(1), "--Alignment", str(1), "coupling"]
    for name, val in zip(ccTag,couplingConstants[:,i]):
        run_array.append(name)
        if val.imag == 0:
            run_array.append(str(val.real)) 
        else:
            sign = "-" if val.imag < 0 else "+"
            run_array.append(str(val.real)+sign+str(abs(val.imag))+"j")   
    print("SM" if i == 0 else fileName[1:])
    subprocess.run(run_array)
    run_array = ["./CRADLE++", "nucleus", "-n", "18F", "-Z", "9", "-A", "18", "general", "-l", "500000", "-t", "20", "-v", "1", "-o", f"build/18F_polZpos{fileName}.txt", "betaDecay","-E", "0", "-r", "0", "--PolarisationZ", "1", "--PolarisationMag", str(1), "--Alignment", str(1), "coupling"]
    for name, val in zip(ccTag,couplingConstants[:,i]):
        run_array.append(name)
        if val.imag == 0:
            run_array.append(str(val.real)) 
        else:
            sign = "-" if val.imag < 0 else "+"
            run_array.append(str(val.real)+sign+str(abs(val.imag))+"j")   
    subprocess.run(run_array)
        
