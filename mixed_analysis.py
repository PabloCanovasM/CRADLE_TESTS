#cscsp

alpha = 0.0072973525664
Z = 19 #need Z of daughter nuclei, 39K
gamma = np.sqrt(1-Z**2*alpha**2)
lambda_jj = 0.4
Lambda_jj = -4/5
sqrt_jj = np.sqrt(3/5)

def cscsp_real(E):
    xi = 2
    a = -2
    b = 0
    c = 0
    A = 0
    B = 0
    D = 0
    return (a/xi, b/xi, c/xi, A/xi, B/xi, D/xi)

print("cs csp real")
SM_39Ca_polZpos = Ca39_data_analysis("build/39Ca_cscsp_real.txt")
SM_39Ca_polZpos.event_mask()
SM_39Ca_polZpos.kinematical_data()
SM_39Ca_polZpos.kinematical_data_ref()
SM_39Ca_polZpos.theoretical_dist(cscsp_real)
SM_39Ca_polZpos.plot_distribution_histograms(showTheory=True,bins=100,fileName="39Ca_cscsp_real",lwidth=1)

#cscsp
print("cs csp imaginary")

def cscsp_imag(E):
    xi = 2
    a = -2
    b = 0
    c = 0
    A = 0
    B = 0
    D = 0
    return (a/xi, b/xi, c/xi, A/xi, B/xi, D/xi)

SM_39Ca_polZpos = Ca39_data_analysis("build/39Ca_cscsp_imag.txt")
SM_39Ca_polZpos.event_mask()
SM_39Ca_polZpos.kinematical_data()
SM_39Ca_polZpos.kinematical_data_ref()
SM_39Ca_polZpos.theoretical_dist(cscsp_imag)
SM_39Ca_polZpos.plot_distribution_histograms(showTheory=True,bins=100,fileName="39Ca_cscsp_imag",lwidth=1)

#csct
print("cs ct imaginary")
    
def csct_imag(E):
    xi = 2
    a = -1+1/3
    b = 0
    c = -Lambda_jj
    A = 0
    B = 0
    D = 2*sqrt_jj
    return (a/xi, b/xi, c/xi, A/xi, B/xi, D/xi)

SM_39Ca_polZpos = Ca39_data_analysis("build/39Ca_csct_imag.txt")
SM_39Ca_polZpos.event_mask()
SM_39Ca_polZpos.kinematical_data()
SM_39Ca_polZpos.kinematical_data_ref()
SM_39Ca_polZpos.theoretical_dist(csct_imag)
SM_39Ca_polZpos.plot_distribution_histograms(showTheory=True,bins=100,fileName="39Ca_csct_imag",lwidth=1)

#csctp

print("cs ctp real")

def csctp_real(E):
    xi = 2
    a = -1+1/3
    b = 0
    c = -Lambda_jj
    A = sqrt_jj*2
    B = -sqrt_jj*2
    D = 0
    return (a/xi, b/xi, c/xi, A/xi, B/xi, D/xi)

SM_39Ca_polZpos = Ca39_data_analysis("build/39Ca_csctp_real.txt")
SM_39Ca_polZpos.event_mask()
SM_39Ca_polZpos.kinematical_data()
SM_39Ca_polZpos.kinematical_data_ref()
SM_39Ca_polZpos.theoretical_dist(csctp_real)
SM_39Ca_polZpos.plot_distribution_histograms(showTheory=True,bins=100,fileName="39Ca_csctp_real",lwidth=1)

#csca
print("cs ca real")

def csca_real(E):
    coulombCorr = alpha*Z/np.sqrt(E**2/mass_e**2-1)
    xi = 2
    a = -1-1/3
    b = 0
    c = Lambda_jj
    A = 0
    B = 0
    D = sqrt_jj*2*coulombCorr
    return (a/xi, b/xi, c/xi, A/xi, B/xi, D/xi)

SM_39Ca_polZpos = Ca39_data_analysis("build/39Ca_csca_real.txt")
SM_39Ca_polZpos.event_mask()
SM_39Ca_polZpos.kinematical_data()
SM_39Ca_polZpos.kinematical_data_ref()
SM_39Ca_polZpos.theoretical_dist(csca_real)
SM_39Ca_polZpos.plot_distribution_histograms(showTheory=True,bins=100,fileName="39Ca_csca_real",lwidth=1)

#cscap

print("cs cap real")

def cscap_real(E):
    coulombCorr = alpha*Z/np.sqrt(E**2/mass_e**2-1)
    xi = 2
    a = -1-1/3
    b = 0
    c = Lambda_jj
    A = 0
    B = sqrt_jj*2*gamma*mass_e/E
    D = 0
    return (a/xi, b/xi, c/xi, A/xi, B/xi, D/xi)

SM_39Ca_polZpos = Ca39_data_analysis("build/39Ca_cscap_real.txt")
SM_39Ca_polZpos.event_mask()
SM_39Ca_polZpos.kinematical_data()
SM_39Ca_polZpos.kinematical_data_ref()
SM_39Ca_polZpos.theoretical_dist(cscap_real)
SM_39Ca_polZpos.plot_distribution_histograms(showTheory=True,bins=100,fileName="39Ca_cscap_real",lwidth=1)

#cscap

print("cs cap imaginary")
Q = 6524.49 - 2*mass_e

def beta_E(E):
    return np.sqrt(1-mass_e**2/E**2)
    
def cscap_imag(E):
    coulombCorr = alpha*Z/np.sqrt(E**2/mass_e**2-1)
    xi = 2
    a = -1-1/3
    b = 0
    c = Lambda_jj
    A = -sqrt_jj*2*coulombCorr
    B = 0
    D = 0
    return (a/xi, b/xi, c/xi, A/xi, B/xi, D/xi)

print(cscap_imag(mass_e+Q/3),beta_E(mass_e+Q/3))

SM_39Ca_polZpos = Ca39_data_analysis("build/39Ca_cscap_imag.txt")
SM_39Ca_polZpos.event_mask()
SM_39Ca_polZpos.kinematical_data()
SM_39Ca_polZpos.kinematical_data_ref()
SM_39Ca_polZpos.theoretical_dist(cscap_imag)
SM_39Ca_polZpos.plot_distribution_histograms(showTheory=True,bins=100,fileName="39Ca_cscap_imag",lwidth=1)

#cspct

print("csp ct real")

def cspct_real(E):
    xi = 2
    a = -1+1/3
    b = 0
    c = -Lambda_jj
    A = sqrt_jj*2
    B = -sqrt_jj*2
    D = 0
    return (a/xi, b/xi, c/xi, A/xi, B/xi, D/xi)

SM_39Ca_polZpos = Ca39_data_analysis("build/39Ca_cspct_real.txt")
SM_39Ca_polZpos.event_mask()
SM_39Ca_polZpos.kinematical_data()
SM_39Ca_polZpos.kinematical_data_ref()
SM_39Ca_polZpos.theoretical_dist(cspct_real)
SM_39Ca_polZpos.plot_distribution_histograms(showTheory=True,bins=100,fileName="39Ca_cspct_real",lwidth=1)

#cspctp

print("cs ctp imaginary")

def cspctp_imag(E):
    xi = 2
    a = -1+1/3
    b = 0
    c = -Lambda_jj
    A = 0
    B = 0
    D = 2*sqrt_jj
    return (a/xi, b/xi, c/xi, A/xi, B/xi, D/xi)

SM_39Ca_polZpos = Ca39_data_analysis("build/39Ca_cspctp_imag.txt")
SM_39Ca_polZpos.event_mask()
SM_39Ca_polZpos.kinematical_data()
SM_39Ca_polZpos.kinematical_data_ref()
SM_39Ca_polZpos.theoretical_dist(cspctp_imag)
SM_39Ca_polZpos.plot_distribution_histograms(showTheory=True,bins=100,fileName="39Ca_cspctp_imag",lwidth=1)

#cspcap

print("csp cap real")

def cspcap_real(E):
    coulombCorr = alpha*Z/np.sqrt(E**2/mass_e**2-1)
    xi = 2
    a = -1-1/3
    b = 0
    c = Lambda_jj
    A = 0
    B = 0
    D = sqrt_jj*2*coulombCorr
    return (a/xi, b/xi, c/xi, A/xi, B/xi, D/xi)

SM_39Ca_polZpos = Ca39_data_analysis("build/39Ca_cspcap_real.txt")
SM_39Ca_polZpos.event_mask()
SM_39Ca_polZpos.kinematical_data()
SM_39Ca_polZpos.kinematical_data_ref()
SM_39Ca_polZpos.theoretical_dist(cspcap_real)
SM_39Ca_polZpos.plot_distribution_histograms(showTheory=True,bins=100,fileName="39Ca_cspcap_real",lwidth=1)

#cspca

print("csp ca real")

def cspca_real(E):
    coulombCorr = alpha*Z/np.sqrt(E**2/mass_e**2-1)
    xi = 2
    a = -1-1/3
    b = 0
    c = Lambda_jj
    A = 0
    B = sqrt_jj*2*gamma*mass_e/E
    D = 0
    return (a/xi, b/xi, c/xi, A/xi, B/xi, D/xi)

SM_39Ca_polZpos = Ca39_data_analysis("build/39Ca_cspca_real.txt")
SM_39Ca_polZpos.event_mask()
SM_39Ca_polZpos.kinematical_data()
SM_39Ca_polZpos.kinematical_data_ref()
SM_39Ca_polZpos.theoretical_dist(cspca_real)
SM_39Ca_polZpos.plot_distribution_histograms(showTheory=True,bins=100,fileName="39Ca_cspca_real",lwidth=1)

#cspca

print("csp ca imaginary")
    
def cspca_imag(E):
    coulombCorr = alpha*Z/np.sqrt(E**2/mass_e**2-1)
    xi = 2
    a = -1-1/3
    b = 0
    c = Lambda_jj
    A = -sqrt_jj*2*coulombCorr
    B = 0
    D = 0
    return (a/xi, b/xi, c/xi, A/xi, B/xi, D/xi)

print(cscap_imag(mass_e+Q/3),beta_E(mass_e+Q/3))

SM_39Ca_polZpos = Ca39_data_analysis("build/39Ca_cspca_imag.txt")
SM_39Ca_polZpos.event_mask()
SM_39Ca_polZpos.kinematical_data()
SM_39Ca_polZpos.kinematical_data_ref()
SM_39Ca_polZpos.theoretical_dist(cspca_imag)
SM_39Ca_polZpos.plot_distribution_histograms(showTheory=True,bins=100,fileName="39Ca_cspca_imag",lwidth=1)

#ctctp
print("ct ctp real")

def ctctp_real(E):
    xi = 2
    a = 2/3
    b = 0
    c = -2*Lambda_jj
    A = -2*lambda_jj
    B = -2*lambda_jj
    D = 0
    return (a/xi, b/xi, c/xi, A/xi, B/xi, D/xi)

SM_39Ca_polZpos = Ca39_data_analysis("build/39Ca_ctctp_real.txt")
SM_39Ca_polZpos.event_mask()
SM_39Ca_polZpos.kinematical_data()
SM_39Ca_polZpos.kinematical_data_ref()
SM_39Ca_polZpos.theoretical_dist(ctctp_real)
SM_39Ca_polZpos.plot_distribution_histograms(showTheory=True,bins=100,fileName="39Ca_ctctp_real",lwidth=1)

#ctctp

print("ct ctp imaginary")

def ctctp_imag(E):
    xi = 2
    a = 2/3
    b = 0
    c = -2*Lambda_jj
    A = 0
    B = 0
    D = 0
    return (a/xi, b/xi, c/xi, A/xi, B/xi, D/xi)

SM_39Ca_polZpos = Ca39_data_analysis("build/39Ca_ctctp_imag.txt")
SM_39Ca_polZpos.event_mask()
SM_39Ca_polZpos.kinematical_data()
SM_39Ca_polZpos.kinematical_data_ref()
SM_39Ca_polZpos.theoretical_dist(ctctp_imag)
SM_39Ca_polZpos.plot_distribution_histograms(showTheory=True,bins=100,fileName="39Ca_ctctp_imag",lwidth=1)

#ctcv

print("ct cv real")

def ctcv_real(E):
    coulombCorr = alpha*Z/np.sqrt(E**2/mass_e**2-1)
    xi = 2
    a = 1+1/3
    b = 0
    c = -Lambda_jj
    A = 0
    B = 0
    D = -sqrt_jj*2*coulombCorr
    return (a/xi, b/xi, c/xi, A/xi, B/xi, D/xi)

SM_39Ca_polZpos = Ca39_data_analysis("build/39Ca_ctcv_real.txt")
SM_39Ca_polZpos.event_mask()
SM_39Ca_polZpos.kinematical_data()
SM_39Ca_polZpos.kinematical_data_ref()
SM_39Ca_polZpos.theoretical_dist(ctcv_real)
SM_39Ca_polZpos.plot_distribution_histograms(showTheory=True,bins=100,fileName="39Ca_ctcv_real",lwidth=1)

#ctcvp

print("ct cvp real")

def ctcvp_real(E):
    coulombCorr = alpha*Z/np.sqrt(E**2/mass_e**2-1)
    xi = 2
    a = +1+1/3
    b = 0
    c = -Lambda_jj
    A = 0
    B = sqrt_jj*2*gamma*mass_e/E
    D = 0
    return (a/xi, b/xi, c/xi, A/xi, B/xi, D/xi)

SM_39Ca_polZpos = Ca39_data_analysis("build/39Ca_ctcvp_real.txt")
SM_39Ca_polZpos.event_mask()
SM_39Ca_polZpos.kinematical_data()
SM_39Ca_polZpos.kinematical_data_ref()
SM_39Ca_polZpos.theoretical_dist(ctcvp_real)
SM_39Ca_polZpos.plot_distribution_histograms(showTheory=True,bins=100,fileName="39Ca_ctcvp_real",lwidth=1)

#ctcvp

print("ct cvp imaginary")

def ctcvp_imag(E):
    coulombCorr = alpha*Z/np.sqrt(E**2/mass_e**2-1)
    xi = 2
    a = +1+1/3
    b = 0
    c = -Lambda_jj
    A = -sqrt_jj*2*coulombCorr
    B = 0
    D = 0
    return (a/xi, b/xi, c/xi, A/xi, B/xi, D/xi)

print(ctcvp_imag(mass_e+Q/3),beta_E(mass_e+Q/3))

SM_39Ca_polZpos = Ca39_data_analysis("build/39Ca_ctcvp_imag.txt")
SM_39Ca_polZpos.event_mask()
SM_39Ca_polZpos.kinematical_data()
SM_39Ca_polZpos.kinematical_data_ref()
SM_39Ca_polZpos.theoretical_dist(ctcvp_imag)
SM_39Ca_polZpos.plot_distribution_histograms(showTheory=True,bins=100,fileName="39Ca_ctcvp_imag",lwidth=1)

#ctpcvp

print("ctp cvp real")

def ctpcvp_real(E):
    coulombCorr = alpha*Z/np.sqrt(E**2/mass_e**2-1)
    xi = 2
    a = 1+1/3
    b = 0
    c = -Lambda_jj
    A = 0
    B = 0
    D = -sqrt_jj*2*coulombCorr
    return (a/xi, b/xi, c/xi, A/xi, B/xi, D/xi)

SM_39Ca_polZpos = Ca39_data_analysis("build/39Ca_ctpcvp_real.txt")
SM_39Ca_polZpos.event_mask()
SM_39Ca_polZpos.kinematical_data()
SM_39Ca_polZpos.kinematical_data_ref()
SM_39Ca_polZpos.theoretical_dist(ctpcvp_real)
SM_39Ca_polZpos.plot_distribution_histograms(showTheory=True,bins=100,fileName="39Ca_ctpcvp_real",lwidth=1)

#ctpcv

print("ctp cv real")

alpha = 0.0072973525664
Z = 19 #need Z of daughter nuclei, 39K
gamma = np.sqrt(1-Z**2*alpha**2)
lambda_jj = 0.4
Lambda_jj = -4/5
sqrt_jj = np.sqrt(3/5)

def ctpcv_real(E):
    coulombCorr = alpha*Z/np.sqrt(E**2/mass_e**2-1)
    xi = 2
    a = +1+1/3
    b = 0
    c = -Lambda_jj
    A = 0
    B = sqrt_jj*2*gamma*mass_e/E
    D = 0
    return (a/xi, b/xi, c/xi, A/xi, B/xi, D/xi)

SM_39Ca_polZpos = Ca39_data_analysis("build/39Ca_ctpcv_real.txt")
SM_39Ca_polZpos.event_mask()
SM_39Ca_polZpos.kinematical_data()
SM_39Ca_polZpos.kinematical_data_ref()
SM_39Ca_polZpos.theoretical_dist(ctpcv_real)
SM_39Ca_polZpos.plot_distribution_histograms(showTheory=True,bins=100,fileName="39Ca_ctpcv_real",lwidth=1)

#ctpcv

print("ctp cv imaginary")

def ctpcv_imag(E):
    coulombCorr = alpha*Z/np.sqrt(E**2/mass_e**2-1)
    xi = 2
    a = +1+1/3
    b = 0
    c = -Lambda_jj
    A = -sqrt_jj*2*coulombCorr
    B = 0
    D = 0
    return (a/xi, b/xi, c/xi, A/xi, B/xi, D/xi)

print(ctpcv_imag(mass_e+Q/3),beta_E(mass_e+Q/3))

SM_39Ca_polZpos = Ca39_data_analysis("build/39Ca_ctpcv_imag.txt")
SM_39Ca_polZpos.event_mask()
SM_39Ca_polZpos.kinematical_data()
SM_39Ca_polZpos.kinematical_data_ref()
SM_39Ca_polZpos.theoretical_dist(ctcvp_imag)
SM_39Ca_polZpos.plot_distribution_histograms(showTheory=True,bins=100,fileName="39Ca_ctpcv_imag",lwidth=1)

#ctcap

print("ct cap real")

alpha = 0.0072973525664
Z = 19 #need Z of daughter nuclei, 39K
gamma = np.sqrt(1-Z**2*alpha**2)
lambda_jj = 0.4
Lambda_jj = -4/5

def ctcap_real(E):
    xi = 2
    a = 0
    b = -2*gamma
    c = 0
    A = 0
    B = 2*lambda_jj*gamma*mass_e/E
    D = 0
    return (a/xi, b/xi, c/xi, A/xi, B/xi, D/xi)

SM_39Ca_polZpos = Ca39_data_analysis("build/39Ca_ctcap_real.txt")
SM_39Ca_polZpos.event_mask()
SM_39Ca_polZpos.kinematical_data()
SM_39Ca_polZpos.kinematical_data_ref()
SM_39Ca_polZpos.theoretical_dist(ctcap_real)
SM_39Ca_polZpos.plot_distribution_histograms(showTheory=True,bins=100,fileName="39Ca_ctcap_real",lwidth=1)

#ctcap

print("ct cap imaginary")

def ctcap_imag(E):
    coulombCorr = alpha*Z/np.sqrt(E**2/mass_e**2-1)
    xi = 2
    a = -2/3*coulombCorr
    b = 0
    c = 2*coulombCorr*Lambda_jj
    A = 2*coulombCorr*lambda_jj
    B = 0
    D = 0
    return (a/xi, b/xi, c/xi, A/xi, B/xi, D/xi)

SM_39Ca_polZpos = Ca39_data_analysis("build/39Ca_ctcap_imag.txt")
SM_39Ca_polZpos.event_mask()
SM_39Ca_polZpos.kinematical_data()
SM_39Ca_polZpos.kinematical_data_ref()
SM_39Ca_polZpos.theoretical_dist(ctcap_imag)
SM_39Ca_polZpos.plot_distribution_histograms(showTheory=True,bins=100,fileName="39Ca_ctcap_imag",lwidth=1)

#cvcvp

print("cv cvp real")

def cvcvp_real(E):
    xi = 2
    a = 2
    b = 0
    c = 0
    A = 0
    B = 0
    D = 0
    return (a/xi, b/xi, c/xi, A/xi, B/xi, D/xi)

SM_39Ca_polZpos = Ca39_data_analysis("build/39Ca_cvcvp_real.txt")
SM_39Ca_polZpos.event_mask()
SM_39Ca_polZpos.kinematical_data()
SM_39Ca_polZpos.kinematical_data_ref()
SM_39Ca_polZpos.theoretical_dist(cvcvp_real)
SM_39Ca_polZpos.plot_distribution_histograms(showTheory=True,bins=100,fileName="39Ca_cvcvp_real",lwidth=1)

#cvcvp

print("cv cvp imaginary")

def cvcvp_imag(E):
    xi = 2
    a = 2
    b = 0
    c = 0
    A = 0
    B = 0
    D = 0
    return (a/xi, b/xi, c/xi, A/xi, B/xi, D/xi)

SM_39Ca_polZpos = Ca39_data_analysis("build/39Ca_cvcvp_imag.txt")
SM_39Ca_polZpos.event_mask()
SM_39Ca_polZpos.kinematical_data()
SM_39Ca_polZpos.kinematical_data_ref()
SM_39Ca_polZpos.theoretical_dist(cvcvp_imag)
SM_39Ca_polZpos.plot_distribution_histograms(showTheory=True,bins=100,fileName="39Ca_cvcvp_imag",lwidth=1)

#cacap

print("ca cap real")

def cacap_real(E):
    xi = 2
    a = -2/3
    b = 0
    c = 2*Lambda_jj
    A = 2*lambda_jj
    B = -2*lambda_jj
    D = 0
    return (a/xi, b/xi, c/xi, A/xi, B/xi, D/xi)

SM_39Ca_polZpos = Ca39_data_analysis("build/39Ca_cacap_real.txt")
SM_39Ca_polZpos.event_mask()
SM_39Ca_polZpos.kinematical_data()
SM_39Ca_polZpos.kinematical_data_ref()
SM_39Ca_polZpos.theoretical_dist(cacap_real)
SM_39Ca_polZpos.plot_distribution_histograms(showTheory=True,bins=100,fileName="39Ca_cacap_real",lwidth=1)

#cacap

print("ca cap imaginary")

def cacap_imag(E):
    xi = 2
    a = -2/3
    b = 0
    c = 2*Lambda_jj
    A = 0
    B = 0
    D = 0
    return (a/xi, b/xi, c/xi, A/xi, B/xi, D/xi)

SM_39Ca_polZpos = Ca39_data_analysis("build/39Ca_cacap_imag.txt")
SM_39Ca_polZpos.event_mask()
SM_39Ca_polZpos.kinematical_data()
SM_39Ca_polZpos.kinematical_data_ref()
SM_39Ca_polZpos.theoretical_dist(cacap_imag)
SM_39Ca_polZpos.plot_distribution_histograms(showTheory=True,bins=100,fileName="39Ca_cacap_imag",lwidth=1)

#cvcap

print("cv cap real")

def cvcap_real(E):
    xi = 2
    a = 1-1/3
    b = 0
    c = Lambda_jj
    A = -2*sqrt_jj
    B = -2*sqrt_jj
    D = 0
    return (a/xi, b/xi, c/xi, A/xi, B/xi, D/xi)

SM_39Ca_polZpos = Ca39_data_analysis("build/39Ca_cvcap_real.txt")
SM_39Ca_polZpos.event_mask()
SM_39Ca_polZpos.kinematical_data()
SM_39Ca_polZpos.kinematical_data_ref()
SM_39Ca_polZpos.theoretical_dist(cvcap_real)
SM_39Ca_polZpos.plot_distribution_histograms(showTheory=True,bins=100,fileName="39Ca_cvcap_real",lwidth=1)


#cvcap

print("cv cap real")

def cvcap_imaginary(E):
    xi = 2
    a = 1-1/3
    b = 0
    c = Lambda_jj
    A = 0
    B = 0
    D = -2*sqrt_jj
    return (a/xi, b/xi, c/xi, A/xi, B/xi, D/xi)

SM_39Ca_polZpos = Ca39_data_analysis("build/39Ca_cvcap_imag.txt")
SM_39Ca_polZpos.event_mask()
SM_39Ca_polZpos.kinematical_data()
SM_39Ca_polZpos.kinematical_data_ref()
SM_39Ca_polZpos.theoretical_dist(cvcap_real)
SM_39Ca_polZpos.plot_distribution_histograms(showTheory=True,bins=100,fileName="39Ca_cvcap_imag",lwidth=1)
