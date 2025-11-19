
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma, loggamma, spence
from scipy.special import factorial, factorial2
   
def histParam_kin(p_names,names,data,filename):
    for name in p_names:
        var_names = ["px","py","pz"],
        data_mom = data[names == name,5:]
        for j, var_name in enumerate(var_names):
            plt.figure()
            plt.hist(data_mom[:,j], bins=100)
            plt.xlabel(f'{var_name} (keV)')
            plt.ylabel('Counts')
            #plt.yscale('log')
            plt.title(name)
            plt.savefig(f'plots/{filename[6:-4]}_{name}_{var_name}.pdf')
        data_E = data[names == name,4]
        data_v = data_mom/data_E[:,np.newaxis]
        var_names2 = ["vx","vy","vz"]
        for j, var_name in enumerate(var_names2):
            plt.figure()
            plt.hist(data_v[:,j], bins=100)
            plt.xlabel(f'{var_name} (c)')
            plt.ylabel('Counts')
            #plt.yscale('log')
            plt.title(name)
            plt.savefig(f'plots/{filename[6:-4]}_{name}_{var_name}.pdf')
        if name in ['e+','e-']:
            data_dir = data_v/np.tile(np.linalg.norm(data_v,axis=1),(3,1)).T
            var_names3 = ["x_dir","y_dir","z_dir"]
            for j, var_name in enumerate(var_names3):
                plt.figure()
                plt.hist(data_dir[:,j], bins=100)
                plt.xlabel(f'{var_name}')
                plt.ylabel('Counts')
                #plt.yscale('log')
                plt.title(name)
                plt.savefig(f'plots/{filename[6:-4]}_{name}_{var_name}.pdf')
    
    chargedl_data = data[np.logical_or(names == 'e+',names == 'e-'),5:8]
    neutrino_data = data[np.logical_or(names == 'enu',names == 'enubar'),5:8]
    chargedl_data = chargedl_data/np.tile(np.linalg.norm(chargedl_data,axis=1),(3,1)).T
    neutrino_data = neutrino_data/np.tile(np.linalg.norm(neutrino_data,axis=1),(3,1)).T
    angle = np.sum(chargedl_data*neutrino_data,axis=1)
    plt.figure()
    plt.hist(angle, bins=100)
    plt.xlabel('cos(theta_{e,enu})')
    plt.ylabel('Counts')
    #plt.yscale('log')
    plt.savefig(f'plots/{filename[6:-4]}_angular_correlation_enu_e.pdf')
        
    
def phase_space(W, W0, **kwargs):
    """Phase space,
    :param W: Electron energy in iunits of me c^2
    :param W0: Electron endpoint energy in units of me c^2
    """
    return np.sqrt(W**2-1)*(W-W0)**2*W
    
def fermi_function(W, Z, A, **kwargs):
    """Fermi Function
    
    :param Z: Proton number of the final nuclear state
    :param W: Electron energy in units of me c^2
    :param A: Mass number of the final (and initial) state
    """
    f = 1
    ALPHA = 0.0072973525643
        
    # R Nuclear radius in units of the electron Compton wavelength
    
    R = (1.15+1.8*A**(-2./3.)-1.2*A**(-4./3.))*A**(1./3.)*1e-15;
    R *= np.sqrt(5/3)/3.86159268e-13 
        
    if Z == 0:
        return f
    
    g = np.sqrt(1-(ALPHA*Z)**2)
    p = np.sqrt(W**2-1)
    y = ALPHA*Z*W/p
        
    #We use the Fermi function with 2(1+gamma), consistent with CRADLE
    f = (2*(g+1)
        *np.power(2*p*R, 2*(g-1))
        *np.exp(np.pi*y)
        *(np.abs(gamma(g+1.j*y))/(gamma(1+2*g)))**2
        )
    return f


mass_e = 510.9989461 #keV

class data_analysis:
    def __init__(self,fileName,betaType):
        data = np.genfromtxt(fileName, dtype=str)
        self.names = data[:,2]
        self.data = data[:, np.arange(data.shape[1]) != 2].astype(dtype=np.float32)
        self.events = data[:,0].astype(int)
        self.data[np.isnan(self.data)] = 0
        self.c_lepton = "e-" if (betaType == 1) else "e+"
        self.neutrino = "enubar" if (betaType == 1) else "enu"
    def __len__(self):
        return self.data.shape[0]
    def event_mask(self):
        pass
    def kinematical_data(self):
        #energy
        self.data_E_el = self.data[self.names == self.c_lepton ,4]
        
        #directions
        self.data_dir_el = self.data[self.names == self.c_lepton ,5:8]
        self.data_dir_el /= np.tile(np.linalg.norm(self.data_dir_el,axis=1),(3,1)).T
        self.data_dir_enu = self.data[self.names == self.neutrino ,5:8]
        self.data_dir_enu /= np.tile(np.linalg.norm(self.data_dir_enu,axis=1),(3,1)).T
        
        #angle between e and enu
        self.data_angle_elenu = np.sum(self.data_dir_el*self.data_dir_enu,axis=1)
     
    def kinematical_data_ref(self,polDir = np.array([0,0,1])):
        #angles in reference to polDir
        z_cross_elenu = np.sum(np.cross(self.data_dir_el,self.data_dir_enu)*polDir[np.newaxis,:],axis=1)
        self.cos_el = np.sum(self.data_dir_el*polDir[np.newaxis,:],axis=1)
        self.cos_enu = np.sum(self.data_dir_enu*polDir[np.newaxis,:],axis=1)
        sin_el = np.sqrt(1-self.cos_el**2)
        sin_enu = np.sqrt(1-self.cos_enu**2)
        data_cosphi = (self.data_angle_elenu-self.cos_el*self.cos_enu)
        data_sinphi = z_cross_elenu
        self.data_phi = np.arctan2(data_sinphi,data_cosphi)
    def plot_E_dist(self):
        print(self.data_E_el.size)
        plt.hist(self.data_E_el, bins=100)
        plt.xlabel('E')
        plt.ylabel('Counts')
        #plt.yscale('log')
        plt.title("e-")
        plt.show()
    def theoretical_dist(self, coeffs,numComp=False,bins=100):
        #note a, c, A, B and D need to be functions dependent on energy
        E_hist, E_bin = np.histogram(self.data_E_el, bins= 200)
        E_vals = (E_bin[1:]+E_bin[:-1])/2
        N = self.data_E_el.size

        z_el = np.linspace(-1,1,201)
        z_nu = z_el.copy()
        z_elenu = z_el.copy()
        phi = np.linspace(-np.pi,np.pi,501)

        self.z_el = z_el
        self.z_nu = z_nu
        self.z_elenu = z_elenu
        self.phi = phi
        
        self.z_el_theory = np.zeros_like(z_el)
        self.z_nu_theory = np.zeros_like(z_nu)
        self.z_elenu_theory = np.zeros_like(z_elenu)
        self.phi_theory = np.zeros_like(phi)

        for E, dist in zip(E_vals, E_hist):
            a, b, c, A, B, D = coeffs(E)
            beta = np.sqrt(1-mass_e**2/E**2)
            self.z_el_theory += dist/N*(1+b*mass_e/E+A*beta*z_el)/2/(1+b*mass_e/E) 
            self.z_elenu_theory += dist/N*(1+b*mass_e/E+a*beta*z_elenu)/2/(1+b*mass_e/E)
            self.z_nu_theory += dist/N*(1+b*mass_e/E+B*z_nu)/2/(1+b*mass_e/E)
            phi_theory_E = dist/N*(1 + b*mass_e/E + beta*((a+c/3)*np.cos(phi)+D*np.sin(phi))*(np.pi/4)**2)
            self.phi_theory += phi_theory_E/2/np.pi/(1+b*mass_e/E)
            
        if numComp:
            cos_el_hist, cos_el_bins = np.histogram(self.cos_el,bins=bins)
            cos_enu_hist, cos_enu_bins = np.histogram(self.cos_enu,bins=bins)
            cos_elenu_hist, cos_elenu_bins = np.histogram(self.data_angle_elenu,bins=bins)
            phi_hist, phi_bins = np.histogram(self.data_phi,bins=bins)
            
            self.cos_el_bins = (cos_el_bins[1:]+cos_el_bins[:-1])/2
            self.cos_enu_bins = (cos_enu_bins[1:]+cos_enu_bins[:-1])/2
            self.cos_elenu_bins = (cos_elenu_bins[1:]+cos_elenu_bins[:-1])/2
            self.phi_bins = (phi_bins[1:]+phi_bins[:-1])/2
            
            cos_el_theory = np.zeros_like(self.cos_el_bins)
            cos_enu_theory = np.zeros_like(self.cos_enu_bins)
            cos_elenu_theory = np.zeros_like(self.cos_elenu_bins)
            cphi_theory = np.zeros_like(self.phi_bins)
            
            for E, dist in zip(E_vals, E_hist):
                a, b, c, A, B, D = coeffs(E)
                beta = np.sqrt(1-mass_e**2/E**2)
                cos_el_theory += dist/bins*(1+b*mass_e/E+A*beta*self.cos_el_bins)/(1+b*mass_e/E) 
                cos_elenu_theory += dist/bins*(1+b*mass_e/E+a*beta*self.cos_elenu_bins)/(1+b*mass_e/E)
                cos_enu_theory += dist/bins*(1+b*mass_e/E+B*self.cos_enu_bins)/(1+b*mass_e/E)
                phi_theory_E = 1 + b*mass_e/E + beta*(np.pi/4)**2*((a+c/3)*np.cos(self.phi_bins)+
                                                                   D*np.sin(self.phi_bins))
                cphi_theory += dist/bins*phi_theory_E/(1+b*mass_e/E)
                           
            self.res_cos_el = (cos_el_hist-cos_el_theory)/np.sqrt(cos_el_hist)
            self.res_cos_enu = (cos_enu_hist-cos_enu_theory)/np.sqrt(cos_enu_hist)
            self.res_cos_elenu = (cos_elenu_hist-cos_elenu_theory)/np.sqrt(cos_elenu_hist)
            self.res_phi = (phi_hist-cphi_theory)/np.sqrt(phi_hist)
                                 
            self.chi2_cos_el = np.sum(self.res_cos_el**2)
            self.chi2_cos_enu = np.sum(self.res_cos_enu**2)
            self.chi2_cos_elenu = np.sum(self.res_cos_elenu**2)
            self.chi2_phi = np.sum(self.res_phi**2)
            print(self.chi2_cos_el/bins,self.chi2_cos_enu/bins,self.chi2_cos_elenu/bins,self.chi2_phi/bins)
    def theoretical_dist2(self, coeffs, numComp=False,bins=100):
        #note a, c, A, B and D need to be functions dependent on energy
        E_hist, E_bin = np.histogram(self.data_E_el, bins= 200)
        E_vals = (E_bin[1:]+E_bin[:-1])/2
        N = self.data_E_el.size

        z_el = np.linspace(-1,1,201)
        z_nu = z_el.copy()
        z_elenu = z_el.copy()
        phi = np.linspace(-np.pi,np.pi,501)

        self.z_el = z_el
        self.z_nu = z_nu
        self.z_elenu = z_elenu
        self.phi = phi
        
        self.z_el_theory2 = np.zeros_like(z_el)
        self.z_nu_theory2 = np.zeros_like(z_nu)
        self.z_elenu_theory2 = np.zeros_like(z_elenu)
        self.phi_theory2 = np.zeros_like(phi)

        for E, dist in zip(E_vals, E_hist):
            a, b, c, A, B, D = coeffs(E)
            beta = np.sqrt(1-mass_e**2/E**2)
            self.z_el_theory2 += dist/N*(1+b*mass_e/E+A*beta*z_el)/2/(1+b*mass_e/E) 
            self.z_elenu_theory2 += dist/N*(1+b*mass_e/E+a*beta*z_elenu)/2/(1+b*mass_e/E)
            self.z_nu_theory2 += dist/N*(1+b*mass_e/E+B*z_nu)/2/(1+b*mass_e/E)
            phi_theory_E = dist/N*(1 + b*mass_e/E + beta*((a+c/3)*np.cos(phi)+D*np.sin(phi))*(np.pi/4)**2)
            self.phi_theory2 += phi_theory_E/2/np.pi/(1+b*mass_e/E)
        
        if numComp:
            cos_el_hist, cos_el_bins = np.histogram(self.cos_el,bins=bins)
            cos_enu_hist, cos_enu_bins = np.histogram(self.cos_enu,bins=bins)
            cos_elenu_hist, cos_elenu_bins = np.histogram(self.data_angle_elenu,bins=bins)
            phi_hist, phi_bins = np.histogram(self.data_phi,bins=bins)
            
            self.cos_el_bins = (cos_el_bins[1:]+cos_el_bins[:-1])/2
            self.cos_enu_bins = (cos_enu_bins[1:]+cos_enu_bins[:-1])/2
            self.cos_elenu_bins = (cos_elenu_bins[1:]+cos_elenu_bins[:-1])/2
            self.phi_bins = (phi_bins[1:]+phi_bins[:-1])/2
            
            cos_el_theory = np.zeros_like(self.cos_el_bins)
            cos_enu_theory = np.zeros_like(self.cos_enu_bins)
            cos_elenu_theory = np.zeros_like(self.cos_elenu_bins)
            cphi_theory = np.zeros_like(self.phi_bins)
            
            for E, dist in zip(E_vals, E_hist):
                a, b, c, A, B, D = coeffs(E)
                beta = np.sqrt(1-mass_e**2/E**2)
                cos_el_theory += dist/bins*(1+b*mass_e/E+A*beta*self.cos_el_bins)/(1+b*mass_e/E) 
                cos_elenu_theory += dist/bins*(1+b*mass_e/E+a*beta*self.cos_elenu_bins)/(1+b*mass_e/E)
                cos_enu_theory += dist/bins*(1+b*mass_e/E+B*self.cos_enu_bins)/(1+b*mass_e/E)
                phi_theory_E = 1 + b*mass_e/E + beta*(np.pi/4)**2*((a+c/3)*np.cos(self.phi_bins)+
                                                                   D*np.sin(self.phi_bins))
                cphi_theory += dist/bins*phi_theory_E/(1+b*mass_e/E)
                           
            self.res_cos_el_2 = (cos_el_hist-cos_el_theory)/np.sqrt(cos_el_hist)
            self.res_cos_enu_2 = (cos_enu_hist-cos_enu_theory)/np.sqrt(cos_enu_hist)
            self.res_cos_elenu_2 = (cos_elenu_hist-cos_elenu_theory)/np.sqrt(cos_elenu_hist)
            self.res_phi_2 = (phi_hist-cphi_theory)/np.sqrt(phi_hist)
                                 
            self.chi2_cos_el_2 = np.sum(self.res_cos_el_2**2)
            self.chi2_cos_enu_2 = np.sum(self.res_cos_enu_2**2)
            self.chi2_cos_elenu_2 = np.sum(self.res_cos_elenu_2**2)
            self.chi2_phi_2 = np.sum(self.res_phi_2**2)
            print(self.chi2_cos_el_2/bins,self.chi2_cos_enu_2/bins,self.chi2_cos_elenu_2/bins,self.chi2_phi_2/bins)

    def plot_kin_histograms(self,fileName):
        p_names = (set(self.names)&{'e+','e-','enubar','enu'})
        histParam_kin(p_names,self.names,self.data,fileName)
        plt.show()
    def plot_distribution_histograms(self,showTheory=True,compareTheory=False,bins=100,fileName=None,lwidth=2):
        fig, axs = plt.subplots(2,2,figsize=(10,7))
        for ax in axs.flatten():
            ax.tick_params('both',labelsize = 12)
        axs[0,0].hist(self.cos_el, bins=bins,density=True) #z component
        if showTheory:
            axs[0,0].plot(self.z_el,self.z_el_theory,linewidth = lwidth)
            if compareTheory:
                axs[0,0].plot(self.z_el,self.z_el_theory2,linewidth = lwidth)
        axs[0,0].set_xlabel("$\\cos \\theta_{e,j}$",size = 16)
        axs[0,0].set_ylabel('pdf',size = 16)
        #plt.yscale('log')
        axs[1,0].hist(self.cos_enu, bins=bins,density=True) #z component
        if showTheory: 
            axs[1,0].plot(self.z_nu,self.z_nu_theory,linewidth = lwidth)
            if compareTheory:
                axs[1,0].plot(self.z_nu,self.z_nu_theory2,linewidth = lwidth)
        axs[1,0].set_xlabel("$\\cos \\theta_{\\nu,j}$",size = 16)
        axs[1,0].set_ylabel('pdf',size = 16)
        #plt.yscale('log')
        axs[0,1].hist(self.data_angle_elenu, bins=bins,density=True)
        if showTheory: 
            axs[0,1].plot(self.z_elenu,self.z_elenu_theory,linewidth = lwidth)
            if compareTheory:
                axs[0,1].plot(self.z_elenu,self.z_elenu_theory2,linewidth = lwidth)
        axs[0,1].set_xlabel("$\\cos \\theta_{e,\\nu}$",size = 16)
        axs[0,1].set_ylabel('pdf',size = 16)
        #plt.yscale('log')
        axs[1,1].hist(self.data_phi, bins=bins,density=True)
        if showTheory: 
            axs[1,1].plot(self.phi,self.phi_theory,linewidth = lwidth)
            if compareTheory:
                axs[1,1].plot(self.phi,self.phi_theory2,linewidth = lwidth)
        axs[1,1].set_xlabel("$\\phi$",size = 16)
        axs[1,1].set_ylabel('pdf',size = 16)
        #plt.yscale('log')
        fig.tight_layout()
        if fileName != None:
            fig.savefig(f"plots/{fileName}.png")
        plt.show()
    def plot_residuals(self,fileName=None):
        fig, axs = plt.subplots(2,2,figsize=(10,7))
        for ax in axs.flatten():
            ax.tick_params('both',labelsize = 14)
        axs[0,0].plot(self.cos_el_bins,self.res_cos_el,"x") #z component
        axs[0,0].set_xlabel("$\\cos \\theta_{e,j}$",size = 16)
        axs[0,0].set_ylabel('Residual',size = 16)
        #plt.yscale('log')
        axs[1,0].plot(self.cos_enu_bins,self.res_cos_enu,"x") #z component
        axs[1,0].set_xlabel("$\\cos \\theta_{\\nu,j}$",size = 16)
        axs[1,0].set_ylabel('Residual',size = 16)
        #plt.yscale('log')
        axs[0,1].plot(self.cos_elenu_bins,self.res_cos_elenu,"x") #z component
        axs[0,1].set_xlabel("$\\cos \\theta_{e,\\nu}$",size = 16)
        axs[0,1].set_ylabel('Residual',size = 16)
        #plt.yscale('log')
        axs[1,1].plot(self.phi_bins,self.res_phi,"x") #z component
        axs[1,1].set_xlabel("$\\phi$",size = 16)
        axs[1,1].set_ylabel('Residual',size = 16)
        #plt.yscale('log')
        fig.tight_layout()
        if fileName != None:
            fig.savefig(f"plots/{fileName}.png")
        plt.show()
    def plot_2dhist(self,bins=20):
        plt.hist2d(self.cos_el,self.cos_enu,bins=bins)
        plt.colorbar()
        plt.xlabel("$\\cos \\theta_{e,j}$")
        plt.ylabel("$\\cos \\theta_{\\nu,j}$")
        plt.show()
        plt.hist2d(self.cos_el,self.data_phi,bins=bins)
        plt.colorbar()
        plt.xlabel("$\\cos \\theta_{e,j}$")
        plt.ylabel("$\\phi$")
        plt.show()
        plt.hist2d(self.cos_enu,self.data_phi,bins=bins)
        plt.colorbar()
        plt.xlabel("$\\cos \\theta_{\\nu,j}$")
        plt.ylabel("$\\phi$")
        plt.show()
    def plot_residuals2(self,fileName=None):
        fig, axs = plt.subplots(2,2,figsize=(10,7))
        for ax in axs.flatten():
            ax.tick_params('both',labelsize = 14)
        axs[0,0].plot(self.cos_el_bins,self.res_cos_el_2,"x") #z component
        axs[0,0].set_xlabel("$\\cos \\theta_{e,j}$",size = 16)
        axs[0,0].set_ylabel('Residual',size = 16)
        #plt.yscale('log')
        axs[1,0].plot(self.cos_enu_bins,self.res_cos_enu_2,"x") #z component
        axs[1,0].set_xlabel("$\\cos \\theta_{\\nu,j}$",size = 16)
        axs[1,0].set_ylabel('Residual',size = 16)
        #plt.yscale('log')
        axs[0,1].plot(self.cos_elenu_bins,self.res_cos_elenu_2,"x") #z component
        axs[0,1].set_xlabel("$\\cos \\theta_{e,\\nu}$",size = 16)
        axs[0,1].set_ylabel('Residual',size = 16)
        #plt.yscale('log')
        axs[1,1].plot(self.phi_bins,self.res_phi_2,"x") #z component
        axs[1,1].set_xlabel("$\\phi$",size = 16)
        axs[1,1].set_ylabel('Residual',size = 16)
        #plt.yscale('log')
        fig.tight_layout()
        if fileName != None:
            fig.savefig(f"plots/{fileName}.png")
        plt.show()        