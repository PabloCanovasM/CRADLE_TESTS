import subprocess
import numpy as np

ctctp_vals = np.array([[0.0,0.0],[1.0+0.0j, 1.0+0.0j],[-1.0+0.0j, -1.0+0.0j],
                       [0.0+1.0j, 0.0+1.0j],[0.0-1.0j, 0.0-1.0j],
                       [1.0+0.0j, -1.0+0.0j],[0.0+1.0j, 0.0-1.0j]])
ctctp_vals /= np.sqrt(2)

ccTag = ["--CS","--CSP","--CT","--CTP","--CV","--CVP","--CA","--CAP"]
names = ["","_ctposRe","_ctnegRe","_ctposIm","_ctnegIm","_ctoppRe","_ctoppIm"]

for name, pair in zip(names,ctctp_vals):
    couplingConstants = np.zeros(8).astype(np.cdouble)
    couplingConstants[[6,7]] = [1.0, 1.0]
    couplingConstants[2] = pair[0]
    couplingConstants[3] = pair[1]
    run_array = ["./CRADLE++", "nucleus", "-n", "60Co", "-Z", "27", "-A", "60", "general", "-l", "200000", "-t", "20", "-v", "1", "-o", f"build/60Co_polZpos{name}.txt", "betaDecay","-E", "0", "-r", "0", "--PolarisationZ", "1", "--PolarisationMag", str(1), "--Alignment", str(1), "coupling"]
    for tag, val in zip(ccTag,couplingConstants):
        run_array.append(tag)
        if val.imag == 0:
            run_array.append(str(val.real)) 
        else:
            sign = "-" if val.imag < 0 else "+"
            run_array.append(str(val.real)+sign+str(abs(val.imag))+"j")  
    subprocess.run(run_array)

