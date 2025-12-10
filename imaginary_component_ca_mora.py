import subprocess
import numpy as np

angles = np.around(np.arcsin(np.linspace(0,1,6))*180/np.pi)
ca_list = 1.0*np.exp(1j*angles*np.pi/180)
ccTag = ["--CS","--CSP","--CT","--CTP","--CV","--CVP","--CA","--CAP"]

for ang, ca in zip(angles,ca_list):
    print(ang)
    couplingConstants = np.zeros(8).astype(np.cdouble)
    couplingConstants[[4,5]] = 1
    couplingConstants[[6,7]] = ca
    print(couplingConstants)
    couplingConstans_array = []
    for name, val in zip(ccTag,couplingConstants):
        couplingConstans_array.append(name)
        if val.imag == 0:
            couplingConstans_array.append(str(val.real)) 
        else:
            sign = "-" if val.imag < 0 else "+"
            couplingConstans_array.append(str(val.real)+sign+str(abs(val.imag))+"j")     
    run_array_39Ca = ["./CRADLE++", "nucleus", "-n", "39Ca", "-Z", "20", "-A", "39", "general", "-l", str(int(1e6)), "-t", "20", "-v", "1", "-o", f"build/39Ca_MORA_ang_{ang}.txt", "betaDecay","-E", "0", "-r", "0", "--PolarisationZ", "1", "--PolarisationMag", str(1), "--Alignment", str(1), "coupling"]
    run_array_39Ca.extend(couplingConstans_array)
    subprocess.run(run_array_39Ca)
