import ROOT
import vector
import matplotlib.pyplot as plt
import uproot
import numpy as np
import random
from scipy.optimize import curve_fit
from scipy.special import wofz
import awkward as ak
import math


# file : https://opendata.cern.ch/record/12361
# load root file (simulation data of Higgs to four leptons) using uproot

file = uproot.open('SMHiggsToZZTo4L.root')

# bring key names and value from 'Events'  >> make numpy data
tree = file['Events']
tree.keys()
branches = tree.arrays()

# pt_cut
E_pt_cut = 15
M_pt_cut = 24


# fourmuons p4 : Higgs to four muons

four_muons_mask = (branches['nMuon'] == 4 )
muon_p4 = vector.zip({'pt': branches['Muon_pt'], 'eta': branches['Muon_eta'], 'phi': branches['Muon_phi'], 'mass': branches['Muon_mass']})   # make muon_p4 vector with specific data branches (pt, eta, phi,mass)
###


four_muons_p4_0 = muon_p4[four_muons_mask]   # make new data applying 4-muons filter

# muon charge cut_ 4 muons


four_muons_charges = branches['Muon_charge'][four_muons_mask]    # data with charges of four muons

opposite_sign_muons_1_mask = (four_muons_charges[:,0] != four_muons_charges[:,1]) & (four_muons_charges[:,2] != four_muons_charges[:,3]) # conditions : muon charge with diffrent sign

opposite_sign_muons_2_mask = (four_muons_charges[:,0] != four_muons_charges[:,2]) & (four_muons_charges[:,3] != four_muons_charges[:,1])
opposite_sign_muons_3_mask = (four_muons_charges[:,0] != four_muons_charges[:,3]) & (four_muons_charges[:,1] != four_muons_charges[:,2])


fourmuons_1_p4 = four_muons_p4_0[opposite_sign_muons_1_mask]     # make new data applying muons charge filter (+,-) or (-,+)
fourmuons_2_p4 = four_muons_p4_0[opposite_sign_muons_2_mask]

fourmuons_3_p4 = four_muons_p4_0[opposite_sign_muons_3_mask]


fourmuons_p4_1 = np.concatenate((fourmuons_1_p4, fourmuons_2_p4, fourmuons_3_p4), axis=0)     # combine two data > charge 1 mask , charge 2 mask ,charge 3 mask

# select only events where all muons have 20 < pt
mask1_4M_pt = ak.all(fourmuons_p4_1.pt > M_pt_cut, axis=1)

four_muons_p4_2 = fourmuons_p4_1[mask1_4M_pt]

mask2_4M_eta = ak.all( abs(four_muons_p4_2.eta) <3, axis=1)

four_muons_p4 = four_muons_p4_2[mask2_4M_eta]


first_muon_p4 = four_muons_p4[:,0]        # To sum fourvector of 4-leptons, devide it to first, second, third, fourth muons fourvectors
second_muon_p4 = four_muons_p4[:,1]
third_muon_p4 = four_muons_p4[:,2]
fourth_muon_p4 = four_muons_p4[:,3]

sum_muons_p4 = first_muon_p4 + second_muon_p4 + third_muon_p4 +fourth_muon_p4    # sum fourvector of 4-muons


# fourelectorns : Higgs to four electrons

four_electrons_mask = branches['nElectron'] == 4   # make a mask to filter four electrons events

electrons_p4 = vector.zip({'pt': branches['Electron_pt'], 'eta': branches['Electron_eta'], 'phi': branches['Electron_phi'], 'mass': branches['Electron_mass']})    # make electron_p4 vector with specific data branches (pt, eta, phi,mass)

four_electrons_p4_0 = electrons_p4[four_electrons_mask]

# electron charge cut _ 4 electrons

four_electrons_charges = branches['Electron_charge'][four_electrons_mask] # make a mask to filter Electorn charge (+, -)

# conditions : muon charge with diffrent sign
opposite_sign_electrons_1_mask = (four_electrons_charges[:,0] != four_electrons_charges[:,1]) & (four_electrons_charges[:,2] != four_electrons_charges[:,3])
    
opposite_sign_electrons_2_mask = (four_electrons_charges[:,0] != four_electrons_charges[:,2]) & (four_electrons_charges[:,3] != four_electrons_charges[:,1])

opposite_sign_electrons_3_mask = (four_electrons_charges[:,0] != four_electrons_charges[:,3]) & (four_electrons_charges[:,1] != four_electrons_charges[:,2])

fourelectrons_1_p4 = four_electrons_p4_0[opposite_sign_electrons_1_mask]    # make new data applying electrons charge filter (+,-) or (-,+)
fourelectrons_2_p4 = four_electrons_p4_0[opposite_sign_electrons_2_mask]
fourelectrons_3_p4 = four_electrons_p4_0[opposite_sign_electrons_3_mask]
        

four_electrons_p4_1 = np.concatenate((fourelectrons_1_p4, fourelectrons_2_p4, fourelectrons_3_p4), axis=0)   # combine two data > charge 1 mask , charge 2 mask

mask1_4E_pt = ak.all(four_electrons_p4_1.pt > E_pt_cut, axis=1)

four_electrons_p4_2 = four_electrons_p4_1[mask1_4E_pt]

mask2_4E_eta = ak.all( abs(four_electrons_p4_2.eta) <3, axis=1)

four_electrons_p4 = four_electrons_p4_2[mask2_4E_eta]


first_electron_p4 = four_electrons_p4[:, 0]    # To sum fourvector of 4-leptons, devide it to first, second, third, fourth electrons fourvectors
second_electron_p4 = four_electrons_p4[:, 1]
third_electron_p4 = four_electrons_p4[:, 2]
fourth_electron_p4 = four_electrons_p4[:, 3]


sum_electrons_p4 = first_electron_p4 + second_electron_p4 + third_electron_p4 +fourth_electron_p4



# twoE twoM

twoEtwoM_mask = ( branches['nMuon'] == 2) & (branches['nElectron'] == 2)   # make a mask to filter two electrons and two electorns events
 # twoE_mask = ( branches['nElectron'] == 2 )

# make muon_p4 vector with specific data branches (pt, eta, phi,mass)

lep1_p4 = vector.zip({'pt': branches['Muon_pt'], 'eta': branches['Muon_eta'], 'phi': branches['Muon_phi'], 'mass': branches['Muon_mass']})

# make electron_p4 vector with specific data branches (pt, eta, phi,mass)

lep2_p4 = vector.zip({'pt': branches['Electron_pt'], 'eta': branches['Electron_eta'], 'phi': branches['Electron_phi'], 'mass': branches['Electron_mass']})



twoM_p4 = lep1_p4[twoEtwoM_mask]     # make data applying 2-muons filter
twoE_p4 = lep2_p4[twoEtwoM_mask]     # make data applying 2-electrons filter


two_electrons_charges = branches['Electron_charge'][twoEtwoM_mask]
two_muons_charges = branches['Muon_charge'][twoEtwoM_mask]

twoEtwoM_step2_charge = (two_electrons_charges[:,0] != two_electrons_charges[:,1]) & (two_muons_charges[:,0] != two_muons_charges[:,1])

twoelectrons_p4_0 = twoE_p4[twoEtwoM_step2_charge]
twomuons_p4_0 = twoM_p4[twoEtwoM_step2_charge]


mask1_2E_pt = ak.all(twoelectrons_p4_0.pt > E_pt_cut, axis=1)
mask1_2M_pt = ak.all(twomuons_p4_0.pt > M_pt_cut, axis=1)

two_muons_p4_1 = twoelectrons_p4_0[mask1_2M_pt & mask1_2E_pt]
two_electrons_p4_1 = twomuons_p4_0[mask1_2M_pt & mask1_2E_pt]

mask2_2E_eta = ak.all( abs(two_electrons_p4_1.eta) <3, axis=1)
mask2_2M_eta = ak.all( abs(two_muons_p4_1.eta) <3, axis=1)


twomuons_p4 = two_muons_p4_1[mask2_2M_eta & mask2_2E_eta]
twoelectrons_p4 = two_electrons_p4_1[mask2_2E_eta & mask2_2M_eta]

 # charge cut _ 2 electrons and 2 muons
 


firstM_p4 = twomuons_p4[:, 0]
secondM_p4 = twomuons_p4[:, 1]
firstE_p4 = twoelectrons_p4[:, 0]
secondE_p4 = twoelectrons_p4[:, 1]
        
sum_2electrons_p4 = firstE_p4 + secondE_p4    ## >>>> 둘다 20 이내의 mass 인데 합치면

sum_2muons_p4 = firstM_p4 + secondM_p4

#sum_2E_p4 = sum_2electrons_p4[0:28511]      # To sum 2 data, set size of data equally


sum_2E2M_p4 = sum_2electrons_p4 + sum_2muons_p4        # sum fourvector of 2-electrons, 2-muons

# combine all data (twoMuons twoElectrons, fourElectrons, fourMuons)

sum_all_p4 = np.concatenate((sum_2E2M_p4,sum_muons_p4,sum_electrons_p4), axis =0)


# Gaussian fitting


# Define gaussian function
def gaussian(x, amp, cen, wid):
    return amp * np.exp(-(x - cen) ** 2 / wid)

# make histogram with Higgs mass
mass_histogram1 = plt.hist(sum_all_p4.mass, bins=200, range=(75,170))

def Voigt(x, x0, y0, a, sigma, gamma):
    #sigma = alpha / np.sqrt(2 * np.log(2))

    return y0 + a * np.real(wofz((x - x0 + 1j*gamma)/sigma/np.sqrt(2))) / sigma /np.sqrt(2*np.pi)


# save data in hist_N and hist_x
hist1_N = mass_histogram1[0]
hist1_x = []
for i in range(len(mass_histogram1[1])-1):
    mass = mass_histogram1[1]
    hist1_x.append((mass[i] + mass[i+1]) /2.)
    
hist_xx=np.linspace(np.min(hist1_x),np.max(hist1_x),200)

popt1, pcov1 = curve_fit(Voigt, hist1_x, hist1_N,p0=[125, np.max(hist1_N), -(np.max(hist1_N)-np.min(hist1_N)), 2.4,2.4] ,maxfev=5000)
mean_f1 = popt1[0]
sigma_f1 = popt1[3]
gamma_f1 = popt1[4]

covariance_mean_1 = pcov1[0,0]


fv1 = 0.5346 * 2* 2*gamma_f1 + (0.2166*4*gamma_f1**2 + 4*sigma_f1**2*2*np.log(2))**(1/2)


plt.plot(hist_xx, Voigt(hist1_x,*popt1), 'r-', label='Voigtian fit')
plt.grid()   # grid on
plt.legend()
plt.text(80, 520,'mean :'+ str(round(mean_f1,2)) +'±' +str(round(math.sqrt(covariance_mean_1),4)) )        # describe mean on plot
plt.text(80, 470,'sigma :'+ str(round(fv1,2)) )      # describe sigma on plot
plt.xlabel('invariant mass [GeV]')
plt.ylabel('Number of events')
plt.show

plt.savefig('./plot_voigtian_fitting_all.png')   # save plot/

plt.clf()
# 4_muons


mass_histogram2 = plt.hist(sum_muons_p4.mass, bins=200, range=(75,170))
hist2_N = mass_histogram2[0]
hist2_x = []
for i in range(len(mass_histogram2[1])-1):
    mass = mass_histogram2[1]
    hist2_x.append((mass[i] + mass[i+1]) /2.)
    
popt2, pcov2 = curve_fit(Voigt, hist2_x, hist2_N,p0=[125, np.max(hist2_N), -(np.max(hist2_N)-np.min(hist2_N)), 2.4,2.4] ,maxfev=5000)


mean_f2 = popt2[0]
sigma_f2 = popt2[3]
gamma_f2 = popt2[4]

covariance_mean_2 = pcov2[0,0]

fv2 = 0.5346 * 2* 2*gamma_f2 + (0.2166*4*gamma_f2**2 + 4*sigma_f2**2*2*np.log(2))**(1/2)

plt.plot(hist_xx, Voigt(hist2_x,*popt2), 'r-', label='Voigtian fit_4M')
plt.grid()   # grid on
plt.legend()
plt.text(75, 325,'mean_4M :   '+ str(round(mean_f2,2)) +'±' +str(round(math.sqrt(covariance_mean_2),4)) )        # describe mean on plot
plt.text(75,310,'sigma_4M :   '+ str(round(fv2,2)) )      # describe sigma on plot
plt.xlabel('invariant mass [GeV]')
plt.ylabel('Number of events')


#plt.savefig('./plot_voigtian_fitting_4muons.png')   # save plot/

# 4_electrons


mass_histogram3 = plt.hist(sum_electrons_p4.mass, bins=200, range=(75,170))
hist3_N = mass_histogram3[0]
hist3_x = []
for i in range(len(mass_histogram3[1])-1):
    mass = mass_histogram3[1]
    hist3_x.append((mass[i] + mass[i+1]) /2.)
    
popt3, pcov3 = curve_fit(Voigt, hist3_x, hist3_N,p0=[125, np.max(hist3_N), -(np.max(hist3_N)-np.min(hist3_N)), 2.4,2.4] ,maxfev=5000)

mean_f3 = popt3[0]
sigma_f3 = popt3[3]
gamma_f3 = popt3[4]

covariance_mean_3 = pcov3[0,0]

fv3 = 0.5346 * 2* 2*gamma_f3 + (0.2166*4*gamma_f3**2 + 4*sigma_f3**2*2*np.log(2))**(1/2)

plt.plot(hist_xx, Voigt(hist3_x,*popt3), 'b-', label='Voigtian fit_4E')
plt.grid()   # grid on
plt.legend()
plt.text(75, 295,'mean_4E :   '+ str(round(mean_f3,2)) +'±' +str(round(math.sqrt(covariance_mean_3),4)) )        # describe mean on plot
plt.text(75,280 ,'sigma_4E :   '+ str(round(fv3,2)) )      # describe sigma on plot
plt.xlabel('invariant mass [GeV]')
plt.ylabel('Number of events')


#plt.savefig('./plot_voigtian_fitting_4electrons.png')   # save plot/


# 2_electrons and 2_muons

mass_histogram4 = plt.hist(sum_2E2M_p4.mass, bins=200, range=(75,170))
hist4_N = mass_histogram4[0]
hist4_x = []
for i in range(len(mass_histogram4[1])-1):
    mass = mass_histogram4[1]
    hist4_x.append((mass[i] + mass[i+1]) /2.)
    
popt4, pcov4 = curve_fit(Voigt, hist4_x, hist4_N,p0=[125, np.max(hist4_N), -(np.max(hist4_N)-np.min(hist4_N)), 2.4,2.4] ,maxfev=5000)

mean_f4 = popt4[0]
sigma_f4 = popt4[3]
gamma_f4 = popt4[4]

covariance_mean_4 = pcov4[0,0]

fv4 = 0.5346 * 2* 2*gamma_f4 + (0.2166*4*gamma_f4**2 + 4*sigma_f4**2*2*np.log(2))**(1/2)

plt.plot(hist_xx, Voigt(hist4_x,*popt4), 'k-', label='Voigtian fit_2M2E')
plt.grid()   # grid on
plt.legend()
plt.text(75, 265,'mean_2E2M :'+ str(round(mean_f4,2)) +'±' +str(round(math.sqrt(covariance_mean_4),4)) )        # describe mean on plot
plt.text(75, 250,'sigma_2E2M :'+ str(round(fv4,2)) )      # describe sigma on plot
plt.xlabel('invariant mass [GeV]')
plt.ylabel('Number of events')


plt.savefig('./plot_voigtian_fitting_final.png')   # save plot/




