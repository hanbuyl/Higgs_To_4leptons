import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = False

df = ROOT.ROOT.RDataFrame("Events", "SMHiggsToZZTo4L.root")
samples_mu = ["Muon_eta","Muon_phi","Muon_pt","Muon_mass"]
samples_el = ["Electron_eta","Electron_phi","Electron_pt","Electron_mass"]
df_1 = df.Filter("Muon_charge[0] != Muon_charge[1]", "Muons with opposite charge")
df_2 = df.Filter("Electron_charge[0] != Electron_charge[1]", "Electrons with opposite charge")


'''fourMuon = df.Filter("nMuon == 4" )
df_os_mu = fourMuon.Filter("Muon_charge[0] != Muon_charge[1]", "Muons with opposite charge")

fourElectron = df.Filter("nElectron == 4")
df_os_el = fourMuon.Filter("Electron_charge[0] != Electron_charge[1]", "Electrons with opposite charge")


twoMu = df.Filter("nMuon == 2")
twoMutwoEl = twoMu.filter("nElectron == 2" )
df_os_mu_el = twoMutwoEl.Filter("Muon_charge[0] != Muon_charge[1]", "Muons with opposite charge" and "Electron_charge[0] != Electron_charge[1]", "Electrons with opposite charge" )'''


ROOT.gInterpreter.Declare("""
using cRVecF = const ROOT::RVecF &;
bool GoodElectronsAndMuons(const ROOT::RVecI & type, cRVecF pt, cRVecF eta, cRVecF phi, cRVecF e, cRVecF trackd0pv, cRVecF tracksigd0pv, cRVecF z0)
{
    for (size_t i = 0; i < type.size(); i++) {
        ROOT::Math::PtEtaPhiEVector p(pt[i] / 1000.0, eta[i], phi[i], e[i] / 1000.0);
        if (type[i] == 11) {
            if (pt[i] < 7000 || abs(eta[i]) > 2.47 || abs(trackd0pv[i] / tracksigd0pv[i]) > 5 || abs(z0[i] * sin(p.Theta())) > 0.5) return false;
        } else {
            if (abs(trackd0pv[i] / tracksigd0pv[i]) > 5 || abs(z0[i] * sin(p.Theta())) > 0.5) return false;
        }
    }
    return true;
}
""")

df_good_4mu = df_1.Filter("df_good_4mu","abs(Muon_eta) < 2.5 && Muon_pt > 5000 && nMuon == 4 ")
df_good_4el = df_2.Filter("df_good_4el","abs(Electron_eta) < 2.5 && Electron_pt > 5000 && nElectron == 4 ")
df_good_2el_2mu = df.Filter("abs(Muon_eta) < 2.5 && Muon_pt > 5000 && nMuon == 2 && abs(Electron_eta) < 2.5 && Electron_pt > 5000 && nElectron == 4 && Muon_charge[0] != Muon_charge[1] && nElectron == 2 && Electron_charge[0] != Electron_charge[1]")

'''
for s in samples:

    # Select events with exactly four good leptons conserving charge and lepton numbers
    # Note that all collections are RVecs and good_lep is the mask for the good leptons.
    # The lepton types are PDG numbers and set to 11 or 13 for an electron or muon
    # irrespective of the charge.
    df[s] = df[s].Define("good_lep", "abs(lep_eta) < 2.5 && lep_pt > 5000 && lep_ptcone30 / lep_pt < 0.3 && lep_etcone20 / lep_pt < 0.3")\
                 .Filter("Sum(good_lep) == 4")\
                 .Filter("Sum(lep_charge[good_lep]) == 0")\
                 .Define("goodlep_sumtypes", "Sum(lep_type[good_lep])")\
                 .Filter("goodlep_sumtypes == 44 || goodlep_sumtypes == 52 || goodlep_sumtypes == 48")

    # Apply additional cuts depending on lepton flavour
    df[s] = df[s].Filter("GoodElectronsAndMuons(lep_type[good_lep], lep_pt[good_lep], lep_eta[good_lep], lep_phi[good_lep], lep_E[good_lep], lep_trackd0pvunbiased[good_lep], lep_tracksigd0pvunbiased[good_lep], lep_z0[good_lep])")

    # Create new columns with the kinematics of good leptons
    df[s] = df[s].Define("goodlep_pt", "lep_pt[good_lep]")\
                 .Define("goodlep_eta", "lep_eta[good_lep]")\
                 .Define("goodlep_phi", "lep_phi[good_lep]")\
                 .Define("goodlep_E", "lep_E[good_lep]")

    # Select leptons with high transverse momentum
    df[s] = df[s].Filter("goodlep_pt[0] > 25000 && goodlep_pt[1] > 15000 && goodlep_pt[2] > 10000")
'''


# sum of four leptons
ROOT.gInterpreter.Declare("""
float ComputeInvariantMass(cRVecF pt, cRVecF eta, cRVecF phi, cRVecF e)
{
    ROOT::Math::PtEtaPhiEVector p1(pt[0], eta[0], phi[0], e[0]);
    ROOT::Math::PtEtaPhiEVector p2(pt[1], eta[1], phi[1], e[1]);
    ROOT::Math::PtEtaPhiEVector p3(pt[2], eta[2], phi[2], e[2]);
    ROOT::Math::PtEtaPhiEVector p4(pt[3], eta[3], phi[3], e[3]);
    return (p1 + p2 + p3 + p4).M() / 1000;
}
""")

def ComputeInvariantMass(pt, eta, phi, e):
    ROOT.Math.Pt

# make histograms

histos = {}
for s in samples_mu:
    df_good_4mu[s] = df_good_4mu[s].Define("m4l", "ComputeInvariantMass(Muon_pt, Muon_eta, Muon_phi, Muon_mass)")
    histos[s] = df_good_4mu[s].Histo1D(ROOT.RDF.TH1DModel(s, "m4l", 24, 80, 170), "m4l", "weight")



'''
#mass
df_mass_mu = df_os_mu.Define("Dimuon_mass", "InvariantMass(Muon_pt, Muon_eta, Muon_phi, Muon_mass)")
df_mass_el = df_os_el.Define("DiElectron_mass", "InvariantMass(Electron_pt, Electron_eta, Electron_phi, Electron_mass)")
df_mass_el = df_os_mu_el.Define("DiElectron_mass", "InvariantMass(Electron_pt, Electron_eta, Electron_phi, Electron_mass)")
df_mass_mu = df_os_mu_el.Define("DiElectron_mass", "InvariantMass(Muon_pt, Muon_eta, Muon_phi, Muon_mass)")

'''


#h = df_mass_mu.Histo1D(("Dimuon_mass", "Dimuon mass;m_{#mu#mu} (GeV);N_{Events}", 30000, 0.25, 300), "Dimuon_mass")
#h1 = df_mass_el.Histo1D(("DiElectron_mass", "DiElectron mass;m_{#el#el} (GeV);N_{Events}", 30000, 0.25, 300), "DiElectron_mass")







c = ROOT.TCanvas("c", "", 600, 600)
pad = ROOT.TPad("upper_pad", "", 0, 0, 1, 1)
pad.SetTickx(False)
pad.SetTicky(False)
pad.Draw()
pad.cd()



















#plots
h.SetTitle("")
h.GetXaxis().SetTitleSize(0.04)
h.GetYaxis().SetTitleSize(0.04)
h.Draw()

label = ROOT.TLatex(); label.SetNDC(True)


label.DrawLatex(0.755, 0.680, "Z")
label.SetTextSize(0.040); label.DrawLatex(0.100, 0.920, "#bf{CMS Open Data}")
label.SetTextSize(0.030); label.DrawLatex(0.630, 0.920, "#sqrt{s} = 8 TeV, L_{int} = 11.6 fb^{-1}")

c.SaveAs("dimuon_spectrum.pdf")

