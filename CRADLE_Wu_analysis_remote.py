import data_analysis as DA
import numpy as np
import matplotlib.pyplot as plt

mu_60Co = 3.799*3.152451e-8 #ev.T^-1
kb = 8.617e-5 #ev.K^-1
B = 1 #T
   
T = np.linspace(0.001,0.01,10)
   
mj_states = np.arange(-5,6)
    
p_states = np.exp(mu_60Co*B*mj_states[:,np.newaxis]/5/kb/T[np.newaxis,:])
print(p_states.shape)
Z = p_states.sum(axis=0)
print(Z.shape)
p_states /= Z[np.newaxis,:]
for i in range(10):
    print(p_states[:,i])
    
polz_list = np.sum(p_states*mj_states[:,np.newaxis]/5,axis=0)
Jzz = np.sum(p_states*mj_states[:,np.newaxis]**2,axis=0)
align_list = (3*Jzz-5*6)/5/9
    
print(polz_list)
print(align_list)
    
p_states_ext = np.exp(mu_60Co*B*mj_states/kb/1e-5)
Z_ext = p_states_ext.sum()
p_states_ext /= Z_ext
polz_ext = np.sum(p_states_ext*mj_states)/5
Jzz_ext = np.sum(p_states_ext*mj_states**2)
align_ext = (3*Jzz_ext-5*6)/5/9
print(p_states_ext,polz_ext,Jzz_ext,align_ext)


class Co60_data_analysis(DA.data_analysis):
    def __init__(self, fileName):
        super().__init__(fileName, 1)
    def event_mask(self):
        part_events = np.bincount(self.events)
        events_mask = np.arange(part_events.size)[part_events == 5] #60Ni + e- + enu + 2 gamma
        self.data = self.data[np.isin(self.events, events_mask),:]
        self.names = self.names[np.isin(self.events, events_mask)]
        
#Standard Model Gamov-Teller

def SM_coeffs(E):
    return (-1/3, 0, 1, -1, 1, 0)

cos_lim = np.cos(15*np.pi/180)
   
print("T = 0")
SM_60Co_polZpos = Co60_data_analysis("build/60Co_polZWu0.txt")
print("Data read")
SM_60Co_polZpos.event_mask()
print("Mask applied")
SM_60Co_polZpos.kinematical_data()
SM_60Co_polZpos.kinematical_data_ref()
print("Kinematical data interpreted")
SM_60Co_polZpos.theoretical_dist(SM_coeffs)
print("Generating plots")
SM_60Co_polZpos.plot_distribution_histograms(showTheory=True,bins=100)
el_zdist = getattr(SM_60Co_polZpos, "cos_el")
event_count = [getattr(SM_60Co_polZpos,"data_E_el").size]
count_par = [np.sum(el_zdist > cos_lim)]
count_antipar = [np.sum(el_zdist < -cos_lim)]


for i in range(10):
    def SM_coeffs(E):
        return (-1/3, 0, align_list[i], -1*polz_list[i], polz_list[i], 0)
    print(i+1,T[i],np.exp(-mu_60Co*B/kb/T[i]),polz_list[i],align_list[i])
    
    SM_60Co_polZpos = Co60_data_analysis(f"build/60Co_polZWu{i+1}.txt")
    SM_60Co_polZpos.event_mask()
    SM_60Co_polZpos.kinematical_data()
    SM_60Co_polZpos.kinematical_data_ref()
    #SM_60Co_polZpos.theoretical_dist(SM_coeffs)
    #SM_60Co_polZpos.plot_distribution_histograms(showTheory=True,bins=100)
    el_zdist = getattr(SM_60Co_polZpos, "cos_el")
    event_count.append(getattr(SM_60Co_polZpos,"data_E_el").size)
    count_par.append(np.sum(el_zdist > cos_lim))
    count_antipar.append(np.sum(el_zdist < -cos_lim))

T_theory = np.linspace(0.00001,0.01,1000)
p_states_theory = np.exp(mu_60Co*B*mj_states[:,np.newaxis]/5/kb/T_theory[np.newaxis,:])
Z_theory = p_states_theory.sum(axis=0)
p_states_theory /= Z_theory[np.newaxis,:]
polz_list_theory = np.sum(p_states_theory*mj_states[:,np.newaxis]/5,axis=0)
print(polz_list_theory[0])

#calculation of <beta>
W_0 = 317.05/DA.mass_e + 1
E_range = np.linspace(DA.mass_e+0.1,317.05+DA.mass_e,800)
dist_E = []
av_beta = 0
for E in E_range:
    W = E/DA.mass_e
    dist = DA.phase_space(W,W_0)*DA.fermi_function(W,28,60)
    dist_E.append(dist)
    av_beta += dist*np.sqrt(1-1/W**2)
dist_E = np.array(dist_E)
av_beta /= dist_E.sum()

polz_list_theory *= av_beta
count_par_theory = (2-polz_list_theory*(1+cos_lim))/2
count_antipar_theory = (2+polz_list_theory*(1+cos_lim))/2

T_plot = np.linspace(0,10,11)
fig, ax = plt.subplots()
error_count_par = 2*np.sqrt(np.array(count_par))/np.array(event_count)/(1-cos_lim)
count_par = 2*np.array(count_par)/np.array(event_count)/(1-cos_lim)
error_count_antipar = 2*np.sqrt(np.array(count_antipar))/np.array(event_count)/(1-cos_lim)
count_antipar = 2*np.array(count_antipar)/np.array(event_count)/(1-cos_lim)

ax.tick_params("both",labelsize=14)
ax.set_xlabel("T (mK)",size=16)
ax.set_ylabel("$\\frac{Counts}{Counts (T>>0)}$",size=16)
ax.errorbar(T_plot,count_par,yerr=error_count_par, linewidth = 0, elinewidth = 1, 
            marker = "x",markersize = 4,label = "Parallel to J")
ax.plot(T_theory*1e3,count_par_theory,label = "Parallel to J Theory")
ax.errorbar(T_plot,count_antipar,yerr=error_count_antipar, linewidth = 0, elinewidth = 1, 
            marker = "x",markersize = 4 ,label = "Antiparallel to J")
ax.plot(T_theory*1e3,count_antipar_theory,label = "Antiparallel to J Theory")
ax.legend()
fig.tight_layout()
fig.savefig("plots/Wu_experiment.png")

