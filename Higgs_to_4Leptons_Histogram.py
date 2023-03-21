import ROOT
import vector
import matplotlib.pyplot as plt
import uproot
import numpy as np
import random
from scipy.optimize import curve_fit


# file : https://opendata.cern.ch/record/12361
# load root file (simulation data of Higgs to four leptons) using uproot

file = uproot.open('SMHiggsToZZTo4L.root')

# bring key names and value from 'Events'  >> make numpy data
tree = file['Events']
tree.keys()
branches = tree.arrays()

# fourmuons p4 : Higgs to four muons

four_muons_mask = branches['nMuon'] == 4   # make a mask to filter four muons events


muon_p4 = vector.zip({'pt': branches['Muon_pt'], 'eta': branches['Muon_eta'], 'phi': branches['Muon_phi'], 'mass': branches['Muon_mass']})   # make muon_p4 vector with specific data branches (pt, eta, phi,mass)

four_muons_p4 = muon_p4[four_muons_mask]    # make new data applying 4-muons filter

first_muon_p4 = four_muons_p4[:, 0]        # To sum fourvector of 4-leptons, devide it to first, second, third, fourth muons fourvectors
second_muon_p4 = four_muons_p4[:, 1]
third_muon_p4 = four_muons_p4[:, 2]
fourth_muon_p4 = four_muons_p4[:, 3]

sum_muons_p4 = first_muon_p4 + second_muon_p4 + third_muon_p4 +fourth_muon_p4    # sum fourvector of 4-muons


# muon charge cut_ 4 muons

four_muons_charges = branches['Muon_charge'][four_muons_mask]    # data with charges of four muons

opposite_sign_muons_1_mask = (four_muons_charges[:,0] != four_muons_charges[:,1]) & (four_muons_charges[:,2] != four_muons_charges[:,3])  # conditions : muon charge with diffrent sign
    
opposite_sign_muons_2_mask = (four_muons_charges[:,0] != four_muons_charges[:,2]) & (four_muons_charges[:,3] != four_muons_charges[:,1])

opposite_sign_muons_3_mask = (four_muons_charges[:,0] != four_muons_charges[:,3]) & (four_muons_charges[:,1] != four_muons_charges[:,2])

fourmuons_1_p4 = sum_muons_p4[opposite_sign_muons_1_mask]     # make new data applying muons charge filter (+,-) or (-,+)
fourmuons_2_p4 = sum_muons_p4[opposite_sign_muons_2_mask]

fourmuons_3_p4 = sum_muons_p4[opposite_sign_muons_3_mask]

fourmuons_p4 = np.concatenate((fourmuons_1_p4, fourmuons_2_p4, fourmuons_3_p4), axis=0)     # combine two data > charge 1 mask , charge 2 mask ,charge 3 mask



# fourelectorns : Higgs to four electrons

four_electrons_mask = branches['nElectron'] == 4   # make a mask to filter four electrons events
four_electrons_charges = branches['Electron_charge'][four_electrons_mask] # make a mask to filter Electorn charge (+, -)

electron_p4 = vector.zip({'pt': branches['Electron_pt'], 'eta': branches['Electron_eta'], 'phi': branches['Electron_phi'], 'mass': branches['Electron_mass']})    # make electron_p4 vector with specific data branches (pt, eta, phi,mass)


four_electrons_p4 = electron_p4[four_electrons_mask]   # make new data applying 4-electrons filter


first_electron_p4 = four_electrons_p4[:, 0]    # To sum fourvector of 4-leptons, devide it to first, second, third, fourth electrons fourvectors
second_electron_p4 = four_electrons_p4[:, 1]
third_electron_p4 = four_electrons_p4[:, 2]
fourth_electron_p4 = four_electrons_p4[:, 3]


sum_electrons_p4 = first_electron_p4 + second_electron_p4 + third_electron_p4 +fourth_electron_p4    # sum fourvector of 4-electrons



# electron charge cut _ 4 electrons


# conditions : muon charge with diffrent sign
opposite_sign_electrons_1_mask = (four_electrons_charges[:,0] != four_electrons_charges[:,1]) & (four_electrons_charges[:,2] != four_electrons_charges[:,3])
    
opposite_sign_electrons_2_mask = (four_electrons_charges[:,0] != four_electrons_charges[:,2]) & (four_electrons_charges[:,3] != four_electrons_charges[:,1])

opposite_sign_electrons_3_mask = (four_electrons_charges[:,0] != four_electrons_charges[:,3]) & (four_electrons_charges[:,1] != four_electrons_charges[:,2])

fourelectrons_1_p4 = sum_electrons_p4[opposite_sign_electrons_1_mask]    # make new data applying electrons charge filter (+,-) or (-,+)
fourelectrons_2_p4 = sum_electrons_p4[opposite_sign_electrons_2_mask]
fourelectrons_3_p4 = sum_electrons_p4[opposite_sign_electrons_3_mask]
        

fourelectrons_p4 = np.concatenate((fourelectrons_1_p4, fourelectrons_2_p4, fourelectrons_3_p4), axis=0)   # combine two data > charge 1 mask , charge 2 mask




# twoE twoM

twoEtwoM_mask = ( branches['nMuon'] == 2) & (branches['nElectron'] == 2)   # make a mask to filter two electrons and two electorns events
 # twoE_mask = ( branches['nElectron'] == 2 )
two_electrons_charges = branches['Electron_charge'][twoEtwoM_mask]
two_muons_charges = branches['Muon_charge'][twoEtwoM_mask]

# make muon_p4 vector with specific data branches (pt, eta, phi,mass)

lep1_p4 = vector.zip({'pt': branches['Muon_pt'], 'eta': branches['Muon_eta'], 'phi': branches['Muon_phi'], 'mass': branches['Muon_mass']})

# make electron_p4 vector with specific data branches (pt, eta, phi,mass)

lep2_p4 = vector.zip({'pt': branches['Electron_pt'], 'eta': branches['Electron_eta'], 'phi': branches['Electron_phi'], 'mass': branches['Electron_mass']})



twoM_p4 = lep1_p4[twoEtwoM_mask]     # make data applying 2-muons filter
twoE_p4 = lep2_p4[twoEtwoM_mask]     # make data applying 2-electrons filter


 # charge cut _ 2 electrons and 2 muons
 
two_electrons_charges = branches['Electron_charge'][twoEtwoM_mask]
two_muons_charges = branches['Muon_charge'][twoEtwoM_mask]

twoEtwoM_step2_charge = (two_electrons_charges[:,0] != two_electrons_charges[:,1]) & (two_muons_charges[:,0] != two_muons_charges[:,1])

twoelectrons_p4 = twoE_p4[twoEtwoM_step2_charge]
twomuons_p4 = twoM_p4[twoEtwoM_step2_charge]



firstM_p4 = twomuons_p4[:, 0]
secondM_p4 = twomuons_p4[:, 1]
firstE_p4 = twoelectrons_p4[:, 0]
secondE_p4 = twoelectrons_p4[:, 1]
        
sum_2electrons_p4 = firstE_p4 + secondE_p4
sum_2muons_p4 = firstM_p4 + secondM_p4

#sum_2E_p4 = sum_2electrons_p4[0:28511]      # To sum 2 data, set size of data equally
#sum_2M_p4 = sum_2muons_p4[0:28511]


sum_2E2M_p4 = sum_2electrons_p4 + sum_2muons_p4
# sum_2electrons_p4 + sum_2muons_p4        # sum fourvector of 2-electrons, 2-muons

# combine all data (twoMuons twoElectrons, fourElectrons, fourMuons)

sum_all_p4 = np.concatenate((sum_2E2M_p4,fourmuons_p4,fourelectrons_p4), axis =0)



# Gaussian fitting


# Define gaussian function
def gaussian(x, amp, cen, wid):
    return amp * np.exp(-(x - cen) ** 2 / wid)

# make histogram with Higgs mass
mass_histogram = plt.hist(sum_all_p4.mass, bins=300, range=(75,170))

# save data in hist_N and hist_x
hist_N = mass_histogram[0]
hist_x = []
for i in range(len(mass_histogram[1])-1):
    mass = mass_histogram[1]
    hist_x.append((mass[i] + mass[i+1]) /2.)


# Gaussian fit

# Calculate roughly mean and sigma

mean = sum(hist_N*hist_x)/sum(hist_N)
sigma = sum(hist_N*(hist_x-mean)**2)/sum(hist_N)

# curve_fit >> gaussian fitting
popt, pcov = curve_fit(gaussian, hist_x, hist_N,p0=[max(hist_N),mean,sigma],maxfev=5000)
# A = popt[0] , mu = popt[1] , sigma = popt[2]


hist_xx=np.linspace(np.min(hist_x),np.max(hist_x),300)


# plotting

#plt.hist(sum_all_p4.mass, bins=200, range=(75,170), label = 'Higgs Data')    # make histogram with mass

plt.plot(hist_xx, gaussian(hist_x,*popt), 'r-', label='Gaussian fit')
plt.grid()   # grid on
plt.legend()
plt.text(149, 4750,'mean :'+ str(round(popt[1],2)) )        # describe mean on plot
plt.text(149, 4550,'sigma :'+ str(round(popt[2],2)) )      # describe sigma on plot
plt.xlabel('invariant mass [GeV]')
plt.ylabel('Number of events')
plt.show

plt.savefig('./plot_gaussian_fitting.png')   # save plot/


