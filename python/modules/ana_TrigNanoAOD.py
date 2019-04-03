#!/usr/bin/env python
import os, sys
import math
import array
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class TrigEffAnalysis(Module):
    def __init__(self, era, dataset):
	self.Region = "signal"
	self.Region = "QCD"
	#self.baseline = "loose_baseline"
	self.baseline = "highdm"
	#self.baseline = "lowdm"
	self.maxEvents = -1
	self.Year = era
	#print "type(dataset) is ", type(dataset), "dataset is", dataset
	self.Dataset = dataset
	self.nEvents = 0
        self.writeHistFile=True

    def beginJob(self,histFile=None,histDirName=None):
        Module.beginJob(self,histFile,histDirName)
        met_bin = array.array('f',[100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,350,400,500,800])
        met_bin_len = len(met_bin) - 1 
        mu_bin = array.array('f',[0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,200,250,300,350,400,500])
        mu_bin_len = len(mu_bin) - 1 
        mu_eta_bin = array.array('f',[-2.4, -2.1, -1.6, -1.2, -0.9, -0.3, -0.2, 0.2, 0.3, 0.9, 1.2, 1.6, 2.1, 2.4])
        mu_eta_bin_len = len(mu_eta_bin) - 1 
        eg_eta_bin = array.array('f',[-2.5, -2, -1.566, -1.444, -1, -0.5, 0, 0.5, 1, 1.444, 1.566, 2, 2.5])
        eg_eta_bin_len = len(eg_eta_bin) - 1 
        ele_bin = array.array('f',[0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,250,300,350,400,500])
        ele_bin_len = len(ele_bin) - 1 
        photon_bin = array.array('f',[100,120,140,160,180,200,220,240,260,280,300,320,340,360,380,400,500,600,800])
        photon_bin_len = len(photon_bin) - 1 
        z_bin = array.array('f',[0,20,40,60,80,100,120,140,160,180,200,300,500,800])
        z_bin_len = len(z_bin) - 1 

        self.h_pass_baseline  = ROOT.TH1F("h_pass_baseline" , "; passed baseline" , 2 , 0. , 2. )
        self.h_pass_reftrig  = ROOT.TH1F("h_pass_reftrig" , "; passed ref trigger" , 2 , 0. , 2. )

        self.h_met_all      = ROOT.TH1F("h_met_all" , "; E_{T}^{miss} [GeV]" , met_bin_len, met_bin)
        self.h_met_passtrig = ROOT.TH1F("h_met_passtrig" , "; E_{T}^{miss} [GeV]" , met_bin_len, met_bin)
        self.h_mu_all      = ROOT.TH1F("h_mu_all" , "; Muon p_{T} [GeV]" , mu_bin_len, mu_bin)
        self.h_mu_passtrig = ROOT.TH1F("h_mu_passtrig" , "; Muon p_{T} [GeV]" , mu_bin_len, mu_bin)
        self.h_mu_all_eta      = ROOT.TH1F("h_mu_all_eta" , "; Muon Eta, p_{T} > 50 [GeV]" , mu_eta_bin_len, mu_eta_bin)
        self.h_mu_passtrig_eta = ROOT.TH1F("h_mu_passtrig_eta" , "; Muon Eta, p_{T} > 50 [GeV]" , mu_eta_bin_len, mu_eta_bin)
        self.h_ele_all      = ROOT.TH1F("h_ele_all" , "; Electron p_{T} [GeV]" , ele_bin_len, ele_bin)
        self.h_ele_passtrig = ROOT.TH1F("h_ele_passtrig" , "; Electron p_{T} [GeV]" , ele_bin_len, ele_bin)
        self.h_ele_all_eta      = ROOT.TH1F("h_ele_all_eta" , "; Electron Eta, p_{T} > 40 [GeV]" , eg_eta_bin_len, eg_eta_bin)
        self.h_ele_passtrig_eta = ROOT.TH1F("h_ele_passtrig_eta" , "; Electron Eta, p_{T} > 40 [GeV]" , eg_eta_bin_len, eg_eta_bin)
        self.h_photon_all      = ROOT.TH1F("h_photon_all" , "; Photon p_{T} [GeV]" , photon_bin_len, photon_bin)
        self.h_photon_passtrig = ROOT.TH1F("h_photon_passtrig" , "; Photon p_{T} [GeV]" , photon_bin_len, photon_bin)
        self.h_photon_all_eta      = ROOT.TH1F("h_photon_all_eta" , "; Photon Eta, p_{T} > 200 [GeV]" , eg_eta_bin_len, eg_eta_bin)
        self.h_photon_passtrig_eta = ROOT.TH1F("h_photon_passtrig_eta" , "; Photon Eta, p_{T} > 200 [GeV]" , eg_eta_bin_len, eg_eta_bin)

        self.h_met_all_mid      = ROOT.TH1F("h_met_all_mid" , "; E_{T}^{miss} [GeV]" , met_bin_len, met_bin)
        self.h_met_passtrig_mid  = ROOT.TH1F("h_met_passtrig_mid" , "; E_{T}^{miss} [GeV]" , met_bin_len, met_bin)
        self.h_mu_all_mid       = ROOT.TH1F("h_mu_all_mid" , "; Muon p_{T} [GeV]" , mu_bin_len, mu_bin)
        self.h_mu_passtrig_mid  = ROOT.TH1F("h_mu_passtrig_mid" , "; Muon p_{T} [GeV]" , mu_bin_len, mu_bin)
        self.h_mu_all_eta_mid       = ROOT.TH1F("h_mu_all_eta_mid" , "; Muon Eta, p_{T} > 50 [GeV]" , mu_eta_bin_len, mu_eta_bin)
        self.h_mu_passtrig_eta_mid  = ROOT.TH1F("h_mu_passtrig_eta_mid" , "; Muon Eta, p_{T} > 50 [GeV]" , mu_eta_bin_len, mu_eta_bin)
        self.h_ele_all_mid       = ROOT.TH1F("h_ele_all_mid" , "; Electron p_{T} [GeV]" , ele_bin_len, ele_bin)
        self.h_ele_passtrig_mid  = ROOT.TH1F("h_ele_passtrig_mid" , "; Electron p_{T} [GeV]" , ele_bin_len, ele_bin)
        self.h_ele_all_eta_mid       = ROOT.TH1F("h_ele_all_eta_mid" , "; Electron Eta, p_{T} > 40 [GeV]" , eg_eta_bin_len, eg_eta_bin)
        self.h_ele_passtrig_eta_mid = ROOT.TH1F("h_ele_passtrig_eta_mid" , "; Electron Eta, p_{T} > 40 [GeV]" , eg_eta_bin_len, eg_eta_bin)
        self.h_photon_all_mid      = ROOT.TH1F("h_photon_all_mid" , "; Photon p_{T} [GeV]" , photon_bin_len, photon_bin)
        self.h_photon_passtrig_mid = ROOT.TH1F("h_photon_passtrig_mid" , "; Photon p_{T} [GeV]" , photon_bin_len, photon_bin)
        self.h_photon_all_eta_mid      = ROOT.TH1F("h_photon_all_eta_mid" , "; Photon Eta, p_{T} > 200 [GeV]" , eg_eta_bin_len, eg_eta_bin)
        self.h_photon_passtrig_eta_mid = ROOT.TH1F("h_photon_passtrig_eta_mid" , "; Photon Eta, p_{T} > 200 [GeV]" , eg_eta_bin_len, eg_eta_bin)

        self.h_ele_all_tight       = ROOT.TH1F("h_ele_all_tight" , "; Electron p_{T} [GeV]" , ele_bin_len, ele_bin)
        self.h_ele_passtrig_tight  = ROOT.TH1F("h_ele_passtrig_tight" , "; Electron p_{T} [GeV]" , ele_bin_len, ele_bin)
        self.h_ele_all_eta_tight       = ROOT.TH1F("h_ele_all_eta_tight" , "; Electron Eta, p_{T} > 40 [GeV]" , eg_eta_bin_len, eg_eta_bin)
        self.h_ele_passtrig_eta_tight = ROOT.TH1F("h_ele_passtrig_eta_tight" , "; Electron Eta, p_{T} > 40 [GeV]" , eg_eta_bin_len, eg_eta_bin)

        self.h_zee_all_mid       = ROOT.TH1F("h_zee_all_mid" , "; Z(ee) p_{T} [GeV]" , z_bin_len, z_bin)
        self.h_zee_passtrig_mid       = ROOT.TH1F("h_zee_passtrig_mid" , "; Z(ee) p_{T} [GeV]" , z_bin_len, z_bin)
        self.h_zmumu_all_mid       = ROOT.TH1F("h_zmumu_all_mid" , "; Z(mumu) p_{T} [GeV]" , z_bin_len, z_bin)
        self.h_zmumu_passtrig_mid       = ROOT.TH1F("h_zmumu_passtrig_mid" , "; Z(mumu) p_{T} [GeV]" , z_bin_len, z_bin)

        self.addObject(self.h_pass_baseline )
        self.addObject(self.h_pass_reftrig )

        self.addObject(self.h_met_all )
        self.addObject(self.h_met_passtrig )
        self.addObject(self.h_mu_all )
        self.addObject(self.h_mu_passtrig )
        self.addObject(self.h_mu_all_eta )
        self.addObject(self.h_mu_passtrig_eta )
        self.addObject(self.h_ele_all )
        self.addObject(self.h_ele_passtrig )
        self.addObject(self.h_ele_all_eta )
        self.addObject(self.h_ele_passtrig_eta )
        self.addObject(self.h_photon_all )
        self.addObject(self.h_photon_passtrig )
        self.addObject(self.h_photon_all_eta )
        self.addObject(self.h_photon_passtrig_eta )

        self.addObject(self.h_met_all_mid )
        self.addObject(self.h_met_passtrig_mid )
        self.addObject(self.h_mu_all_mid )
        self.addObject(self.h_mu_passtrig_mid )
        self.addObject(self.h_mu_all_eta_mid )
        self.addObject(self.h_mu_passtrig_eta_mid )
        self.addObject(self.h_ele_all_mid )
        self.addObject(self.h_ele_passtrig_mid )
        self.addObject(self.h_ele_all_eta_mid )
        self.addObject(self.h_ele_passtrig_eta_mid )
        self.addObject(self.h_photon_all_mid )
        self.addObject(self.h_photon_passtrig_mid )
        self.addObject(self.h_photon_all_eta_mid )
        self.addObject(self.h_photon_passtrig_eta_mid )

        self.addObject(self.h_ele_all_tight )
        self.addObject(self.h_ele_passtrig_tight )
        self.addObject(self.h_ele_all_eta_tight )
        self.addObject(self.h_ele_passtrig_eta_tight )

        self.addObject(self.h_zee_all_mid )
        self.addObject(self.h_zee_passtrig_mid )
        self.addObject(self.h_zmumu_all_mid )
        self.addObject(self.h_zmumu_passtrig_mid )

    def mygetattr(self, my_obj, my_branch, default_bool):
	try: getattr(my_obj, my_branch)
	except RuntimeError:
	    #print my_branch, "not found"
	    return default_bool
	else: return getattr(my_obj, my_branch)

    #follow root dPhi calculation: https://root.cern.ch/root/html/src/ROOT__Math__VectorUtil.h.html#60
    def mydPhi(self, phi1, phi2):
	dPhi = phi1 - phi2
	if(dPhi > math.pi): dPhi = dPhi - 2 * math.pi
	elif(dPhi <= -math.pi): dPhi = dPhi + 2 * math.pi
	return abs(dPhi)

    def SelIsotrack(self, isk, met):
        iso = isk.pfRelIso03_chg
        if abs(isk.pdgId) == 11 or abs(isk.pdgId) == 13:
            if isk.pt < 5 or iso > 0.2:
                return False
        if abs(isk.pdgId) == 211:
            if isk.pt < 10 or iso > 0.1:
                return False
        mtW = math.sqrt( 2 * met.pt * isk.pt * (1 - math.cos(met.phi-isk.phi)))
        if mtW  > 100:
            return False
        return True

    def SelPhotons(self, photon):
        #if photon.pt < 200:
        #    return False
        abeta = math.fabs(photon.eta)
        if (abeta > 1.442 and abeta < 1.566) or (abeta > 2.5):
            return False
        ## cut-base ID, 2^0 loose ID
        cutbase =  photon.cutBasedBitmap  if self.Year != "2016" else photon.cutBased
        if not cutbase & 0b1:
            return False
        return True

    def SelJets(self, jet):
        if jet.pt < 20 or math.fabs(jet.eta) > 2.4 :
            return False
        return True

    def analyze(self, event):
    	self.nEvents += 1
	if (self.maxEvents != -1 and self.nEvents > self.maxEvents):
	    return False

	deepCSV_cut = 0.4941
	if(self.Year == "2016"): deepCSV_cut = 0.6324

        met       = Object(event, "MET")
        hlt       = Object(event, "HLT")
        flags     = Object(event, "Flag")
	electrons = Collection(event, "Electron")
	muons	  = Collection(event, "Muon")
	isks 	  = Collection(event, "IsoTrack")
	jets 	  = Collection(event, "Jet")
	photons   = Collection(event, "Photon")
	stop0l    = Object(event, "Stop0l")

	ele_veto = []
	ele_mid = []
	ele_tight = []
	for ele in electrons:
		if (ele.pt > 5 and abs(ele.eta) < 2.5 and ele.cutBasedNoIso >= 1 and ele.miniPFRelIso_all < 0.1):
			ele_veto.append(ele)
			if (ele.cutBasedNoIso >= 3):
				ele_mid.append(ele)
			if (ele.cutBasedNoIso >= 4):
				ele_tight.append(ele)
	n_ele = len(ele_veto)
	n_ele_mid = len(ele_mid)
	n_ele_tight = len(ele_tight)

	zee_mid = []
	if (n_ele_mid == 2 and ele_mid[0].pt > 40 and ele_mid[1].pt > 20 and (ele_mid[0].charge + ele_mid[1].charge) == 0):
		zee_cand = ROOT.TLorentzVector()
		zee_cand = ele_mid[0].p4() + ele_mid[1].p4()
		#Z mass = 91 GeV
		#print "zee_cand mass = ", zee_cand.M()
		if (zee_cand.M() > 81 and zee_cand.M() < 101):
			zee_mid.append(zee_cand)
		
	mu_loose = []
	mu_mid = []
	for mu in muons:
		if (mu.pt > 5 and abs(mu.eta) < 2.4 and mu.miniPFRelIso_all < 0.2):
			mu_loose.append(mu)
			if (mu.mediumId):
				mu_mid.append(mu)
	n_mu = len(mu_loose)
	n_mu_mid = len(mu_mid)

	zmumu_mid = []
	if (n_mu_mid == 2 and mu_mid[0].pt > 50 and mu_mid[1].pt > 20 and (mu_mid[0].charge + mu_mid[1].charge) == 0):
		zmumu_cand = ROOT.TLorentzVector()
		zmumu_cand = mu_mid[0].p4() + mu_mid[1].p4()
		#Z mass = 91 GeV
		#print "zmumu_cand mass = ", zmumu_cand.M()
		if (zmumu_cand.M() > 81 and zmumu_cand.M() < 101):
			zmumu_mid.append(zmumu_cand)

	#n_tau = 0
	#for isk in isks:
	#	if (self.SelIsotrack(isk, met)): n_tau += 1

	photon_loose = []
	for photon in photons:
		if (self.SelPhotons(photon)):
			photon_loose.append(photon)
	n_photon = len(photon_loose)

        n_jets = 0 
        ht = 0
	bjets = []
        for jet in jets:
                ht = ht + jet.pt
                if self.SelJets(jet):
                        n_jets += 1
                        if (jet.btagDeepB > deepCSV_cut):
				bjets.append(jet)
        n_bjets = len(bjets)

	Mtb = float('inf')
	btagidx = sorted(range(len(bjets)), key=lambda k: bjets[k].btagDeepB , reverse=True)
        for i in range(min(len(btagidx), 2)):
            bj = bjets[btagidx[i]]
            Mtb = min(Mtb, math.sqrt( 2 * met.pt * bj.pt * (1 - math.cos(met.phi-bj.phi))))

        if Mtb == float('inf'):
            Mtb = 0

	#print(self.nEvents)
	#print(electrons[0].p4().Pt())
	#for ele in electrons:
	#	print (ele.p4().Pt())
	#print(flag.HBHENoiseFilter)
	#print (self.mygetattr(hlt, 'PFHT180', False))
	#print ("n_ele = ", n_ele, "n_mu = ", n_mu, "n_tau = ", n_tau)

	pass_filter = (
		self.mygetattr(flags, 'goodVertices', True)
		and self.mygetattr(flags, 'HBHENoiseFilter', True)
		and self.mygetattr(flags, 'HBHENoiseIsoFilter', True)
		and self.mygetattr(flags, 'EcalDeadCellTriggerPrimitiveFilter', True)
		and self.mygetattr(flags, 'BadPFMuonFilter', True)
		and self.mygetattr(flags, 'BadChargedCandidateFilter', True)
		and self.mygetattr(flags, 'ecalBadCalibFilter', True)
		and self.mygetattr(flags, 'globalSuperTightHalo2016Filter', True)
		and self.mygetattr(flags, 'eeBadScFilter', True)
		)

	pass_dPhi = False
	pass_dPhi_highdm = False
	if (n_jets == 2):
		if(self.Region == "signal"):
			if (self.mydPhi(jets[0].phi, met.phi) > 0.5 and self.mydPhi(jets[1].phi, met.phi) > 0.15):
				pass_dPhi = True
		if(self.Region == "QCD"):
			if (self.mydPhi(jets[0].phi, met.phi) < 0.1 or self.mydPhi(jets[1].phi, met.phi) < 0.1):
				pass_dPhi = True
	if (n_jets > 2):
		if(self.Region == "signal"):
			if (self.mydPhi(jets[0].phi, met.phi) > 0.5 and self.mydPhi(jets[1].phi, met.phi) > 0.15 and self.mydPhi(jets[2].phi, met.phi) > 0.15):
				pass_dPhi = True
		if(self.Region == "QCD"):
			if (self.mydPhi(jets[0].phi, met.phi) < 0.1 or self.mydPhi(jets[1].phi, met.phi) < 0.1 or self.mydPhi(jets[2].phi, met.phi) < 0.1):
				pass_dPhi = True

	if (n_jets >= 5):
		if(self.Region == "signal"):
			if (self.mydPhi(jets[0].phi, met.phi) > 0.5 and self.mydPhi(jets[1].phi, met.phi) > 0.5 and self.mydPhi(jets[2].phi, met.phi) > 0.5 and self.mydPhi(jets[3].phi, met.phi) > 0.5):
				pass_dPhi_highdm = True
		if(self.Region == "QCD"):
			if (self.mydPhi(jets[0].phi, met.phi) < 0.1 or self.mydPhi(jets[1].phi, met.phi) < 0.1 or self.mydPhi(jets[2].phi, met.phi) < 0.1):
				pass_dPhi_highdm = True

	pass_loose_baseline = (pass_filter and pass_dPhi and ht > 300) 
	pass_highdm = (pass_loose_baseline and pass_dPhi_highdm and n_bjets >= 1)
	pass_lowdm = (pass_loose_baseline and Mtb < 175 and stop0l.nTop == 0 and stop0l.nW == 0
	and stop0l.nResolved == 0 and stop0l.ISRJetPt > 300 and met.pt / math.sqrt(ht) > 10)

	pass_loosejet = (pass_loose_baseline and pass_dPhi_highdm)
	pass_looseb = (pass_loose_baseline and n_bjets >= 1)

	pass_baseline = False
	if(self.baseline == "loose_baseline"): pass_baseline = pass_loose_baseline
	if(self.baseline == "highdm"): pass_baseline = pass_highdm
	if(self.baseline == "lowdm"): pass_baseline = pass_lowdm

	#pass_baseline = pass_loosejet
	#pass_baseline = pass_looseb

        self.h_pass_baseline.Fill(pass_baseline)
	if not (pass_baseline): return False

        refAccept = False
	if(self.Year == "2016"):
		if(self.Dataset == "SingleElectron"):
        		refAccept = hlt.Ele27_WPTight_Gsf
		if(self.Dataset == "SingleMuon"):
        		refAccept = hlt.Mu50
		if(self.Dataset == "SinglePhoton"):
        		refAccept = hlt.Photon175
		if(self.Dataset == "JetHT"):
			refAccept = (
			self.mygetattr(hlt, 'PFHT125', False)
			or self.mygetattr(hlt, 'PFHT200', False)
			or self.mygetattr(hlt, 'PFHT250', False)
			or self.mygetattr(hlt, 'PFHT300', False)
			or self.mygetattr(hlt, 'PFHT350', False)
			or self.mygetattr(hlt, 'PFHT400', False)
			or self.mygetattr(hlt, 'PFHT475', False)
			or self.mygetattr(hlt, 'PFHT600', False)
			or self.mygetattr(hlt, 'PFHT650', False)
			or self.mygetattr(hlt, 'PFHT800', False)
			or self.mygetattr(hlt, 'PFHT900', False)
			or self.mygetattr(hlt, 'CaloJet500_NoJetID', False)
			)

		if(self.Dataset == "MET"):
			refAccept = (
			self.mygetattr(hlt, 'PFMET100_PFMHT100_IDTight', False)
			or self.mygetattr(hlt, 'PFMET110_PFMHT110_IDTight', False)
			or self.mygetattr(hlt, 'PFMET120_PFMHT120_IDTight', False)
			or self.mygetattr(hlt, 'PFMETNoMu100_PFMHTNoMu100_IDTight', False)
			or self.mygetattr(hlt, 'PFMETNoMu110_PFMHTNoMu110_IDTight', False)
			or self.mygetattr(hlt, 'PFMETNoMu120_PFMHTNoMu120_IDTight', False)
			)

	if(self.Year == "2017" or self.Year == "2018"):
		if(self.Dataset == "SingleElectron"):
        		refAccept =(
			self.mygetattr(hlt, 'Ele32_WPTight_Gsf', False)
			or self.mygetattr(hlt, 'Ele35_WPTight_Gsf', False)
			)
		if(self.Dataset == "SingleMuon"):
        		refAccept = hlt.Mu50
		if(self.Dataset == "SinglePhoton"):
        		refAccept = hlt.Photon200
		if(self.Dataset == "JetHT"):
			refAccept =(
			self.mygetattr(hlt, 'PFHT1050', False)
			or self.mygetattr(hlt, 'PFHT180', False)
			or self.mygetattr(hlt, 'PFHT250', False)
			or self.mygetattr(hlt, 'PFHT350', False)
			or self.mygetattr(hlt, 'PFHT370', False)
			or self.mygetattr(hlt, 'PFHT430', False)
			or self.mygetattr(hlt, 'PFHT510', False)
			or self.mygetattr(hlt, 'PFHT590', False)
			or self.mygetattr(hlt, 'PFHT680', False)
			or self.mygetattr(hlt, 'PFHT780', False)
			or self.mygetattr(hlt, 'PFHT890', False)
			)

		if(self.Dataset == "MET"):
			refAccept = (
			self.mygetattr(hlt, 'PFMET110_PFMHT110_IDTight', False)
			or self.mygetattr(hlt, 'PFMET120_PFMHT120_IDTight', False)
			or self.mygetattr(hlt, 'PFMET130_PFMHT130_IDTight', False)
			or self.mygetattr(hlt, 'PFMET140_PFMHT140_IDTight', False)
			or self.mygetattr(hlt, 'PFMETNoMu110_PFMHTNoMu110_IDTight', False)
			or self.mygetattr(hlt, 'PFMETNoMu120_PFMHTNoMu120_IDTight', False)
			or self.mygetattr(hlt, 'PFMETNoMu130_PFMHTNoMu130_IDTight', False)
			or self.mygetattr(hlt, 'PFMETNoMu140_PFMHTNoMu140_IDTight', False)
			)
		    
        self.h_pass_reftrig.Fill(refAccept)
        if not refAccept:
            return False

        sigAccept_met = (
	self.mygetattr(hlt, 'PFMET100_PFMHT100_IDTight', False)
	or self.mygetattr(hlt, 'PFMET110_PFMHT110_IDTight', False)
	or self.mygetattr(hlt, 'PFMET120_PFMHT120_IDTight', False)
	or self.mygetattr(hlt, 'PFMET130_PFMHT130_IDTight', False)
	or self.mygetattr(hlt, 'PFMET140_PFMHT140_IDTight', False)
	or self.mygetattr(hlt, 'PFMETNoMu100_PFMHTNoMu100_IDTight', False)
	or self.mygetattr(hlt, 'PFMETNoMu110_PFMHTNoMu110_IDTight', False)
	or self.mygetattr(hlt, 'PFMETNoMu120_PFMHTNoMu120_IDTight', False)
	or self.mygetattr(hlt, 'PFMETNoMu130_PFMHTNoMu130_IDTight', False)
	or self.mygetattr(hlt, 'PFMETNoMu140_PFMHTNoMu140_IDTight', False)
	or self.mygetattr(hlt, 'PFMET100_PFMHT100_IDTight_PFHT60', False)
	or self.mygetattr(hlt, 'PFMET110_PFMHT110_IDTight_PFHT60', False)
	or self.mygetattr(hlt, 'PFMET120_PFMHT120_IDTight_PFHT60', False)
	or self.mygetattr(hlt, 'PFMET130_PFMHT130_IDTight_PFHT60', False)
	or self.mygetattr(hlt, 'PFMET140_PFMHT140_IDTight_PFHT60', False)
	or self.mygetattr(hlt, 'PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60', False)
	or self.mygetattr(hlt, 'PFMETNoMu110_PFMHTNoMu110_IDTight_PFHT60', False)
	or self.mygetattr(hlt, 'PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60', False)
	or self.mygetattr(hlt, 'PFMETNoMu130_PFMHTNoMu130_IDTight_PFHT60', False)
	or self.mygetattr(hlt, 'PFMETNoMu140_PFMHTNoMu140_IDTight_PFHT60', False)
	#or self.mygetattr(hlt, 'PFMET120_PFMHT120_IDTight_HFCleaned', False)
	#or self.mygetattr(hlt, 'PFMET120_PFMHT120_IDTight_PFHT60_HFCleaned', False)
	#or self.mygetattr(hlt, 'PFMETNoMu120_PFMHTNoMu120_IDTight_HFCleaned', False)
	)

        sigAccept_mu = (
	self.mygetattr(hlt, 'IsoMu20', False)
	or self.mygetattr(hlt, 'IsoMu22', False)
	or self.mygetattr(hlt, 'IsoMu24', False)
	or self.mygetattr(hlt, 'IsoMu27', False)
	or self.mygetattr(hlt, 'IsoMu22_eta2p1', False)
	or self.mygetattr(hlt, 'IsoMu24_eta2p1', False)
	or self.mygetattr(hlt, 'IsoTkMu22', False)
	or self.mygetattr(hlt, 'IsoTkMu24', False)
	or self.mygetattr(hlt, 'Mu50', False)
	or self.mygetattr(hlt, 'Mu55', False)
	)

        sigAccept_ele = (
	self.mygetattr(hlt, 'Ele105_CaloIdVT_GsfTrkIdT', False)
	or self.mygetattr(hlt, 'Ele115_CaloIdVT_GsfTrkIdT', False)
	or self.mygetattr(hlt, 'Ele135_CaloIdVT_GsfTrkIdT', False)
	or self.mygetattr(hlt, 'Ele145_CaloIdVT_GsfTrkIdT', False)
	or self.mygetattr(hlt, 'Ele25_eta2p1_WPTight_Gsf', False)
	or self.mygetattr(hlt, 'Ele20_eta2p1_WPLoose_Gsf', False)
	or self.mygetattr(hlt, 'Ele27_eta2p1_WPLoose_Gsf', False)
	or self.mygetattr(hlt, 'Ele27_WPTight_Gsf', False)
	or self.mygetattr(hlt, 'Ele35_WPTight_Gsf', False)
	or self.mygetattr(hlt, 'Ele20_WPLoose_Gsf', False)
	or self.mygetattr(hlt, 'Ele45_WPLoose_Gsf', False)
	or self.mygetattr(hlt, 'Ele23_Ele12_CaloIdL_TrackIdL_IsoVL', False)
	or self.mygetattr(hlt, 'Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ', False)
	or self.mygetattr(hlt, 'DoubleEle33_CaloIdL_GsfTrkIdVL', False)
	or self.mygetattr(hlt, 'DoubleEle33_CaloIdL_GsfTrkIdVL_MW', False)
	or self.mygetattr(hlt, 'DoubleEle25_CaloIdL_MW', False)
	or self.mygetattr(hlt, 'DoubleEle33_CaloIdL_MW', False)
	)

        sigAccept_photon = (
	self.mygetattr(hlt, 'Photon175', False)
	or self.mygetattr(hlt, 'Photon200', False)
	)

	#veto / loose ID
	if (self.Dataset == "SingleElectron" and n_ele >= 1):
		if (n_mu == 0):
        		self.h_met_all.Fill(met.pt)
        		if (sigAccept_met):
				self.h_met_passtrig.Fill(met.pt)
			if (n_photon >=1):
        			self.h_photon_all.Fill(photon_loose[0].pt)
        			if(photon_loose[0].pt > 200): self.h_photon_all_eta.Fill(photon_loose[0].eta)
        			if (sigAccept_photon):
					self.h_photon_passtrig.Fill(photon_loose[0].pt)
					if(photon_loose[0].pt > 200): self.h_photon_passtrig_eta.Fill(photon_loose[0].eta)
		if (n_mu >= 1):
        		self.h_mu_all.Fill(mu_loose[0].pt)
        		if (mu_loose[0].pt > 50): self.h_mu_all_eta.Fill(mu_loose[0].eta)
        		if (sigAccept_mu):
				self.h_mu_passtrig.Fill(mu_loose[0].pt)
        			if (mu_loose[0].pt > 50): self.h_mu_passtrig_eta.Fill(mu_loose[0].eta)

	if (self.Dataset == "SingleMuon" and n_mu >= 1):
		if (n_ele == 0):
        		self.h_met_all.Fill(met.pt)
        		if (sigAccept_met):
				self.h_met_passtrig.Fill(met.pt)
			if (n_photon >=1):
        			self.h_photon_all.Fill(photon_loose[0].pt)
        			if(photon_loose[0].pt > 200): self.h_photon_all_eta.Fill(photon_loose[0].eta)
        			if (sigAccept_photon):
					self.h_photon_passtrig.Fill(photon_loose[0].pt)
					if(photon_loose[0].pt > 200): self.h_photon_passtrig_eta.Fill(photon_loose[0].eta)
		if (n_ele >= 1):
        		self.h_ele_all.Fill(ele_veto[0].pt)
        		if (ele_veto[0].pt > 40): self.h_ele_all_eta.Fill(ele_veto[0].eta)
        		if (sigAccept_ele):
				self.h_ele_passtrig.Fill(ele_veto[0].pt)
        			if (ele_veto[0].pt > 40): self.h_ele_passtrig_eta.Fill(ele_veto[0].eta)

	if (self.Dataset == "SinglePhoton"):
		if (n_ele == 0 and n_mu == 0):
        		self.h_met_all.Fill(met.pt)
        		if (sigAccept_met):
				self.h_met_passtrig.Fill(met.pt)

	if (self.Dataset == "JetHT"):
		if (n_mu == 0 and n_ele == 0):
        		self.h_met_all.Fill(met.pt)
        		if (sigAccept_met):
				self.h_met_passtrig.Fill(met.pt)
			if (n_photon >=1):
        			self.h_photon_all.Fill(photon_loose[0].pt)
        			if(photon_loose[0].pt > 200): self.h_photon_all_eta.Fill(photon_loose[0].eta)
        			if (sigAccept_photon):
					self.h_photon_passtrig.Fill(photon_loose[0].pt)
					if(photon_loose[0].pt > 200): self.h_photon_passtrig_eta.Fill(photon_loose[0].eta)
		if (n_mu >= 1 and n_ele == 0):
        		self.h_mu_all.Fill(mu_loose[0].pt)
        		if (mu_loose[0].pt > 50): self.h_mu_all_eta.Fill(mu_loose[0].eta)
        		if (sigAccept_mu):
				self.h_mu_passtrig.Fill(mu_loose[0].pt)
        			if (mu_loose[0].pt > 50): self.h_mu_passtrig_eta.Fill(mu_loose[0].eta)
		if (n_mu == 0 and n_ele >= 1):
        		self.h_ele_all.Fill(ele_veto[0].pt)
        		if (ele_veto[0].pt > 40): self.h_ele_all_eta.Fill(ele_veto[0].eta)
        		if (sigAccept_ele):
				self.h_ele_passtrig.Fill(ele_veto[0].pt)
        			if (ele_veto[0].pt > 40): self.h_ele_passtrig_eta.Fill(ele_veto[0].eta)

	if (self.Dataset == "MET"):
		if (n_mu == 0 and n_ele == 0):
			if (n_photon >=1):
        			self.h_photon_all.Fill(photon_loose[0].pt)
        			if(photon_loose[0].pt > 200): self.h_photon_all_eta.Fill(photon_loose[0].eta)
        			if (sigAccept_photon):
					self.h_photon_passtrig.Fill(photon_loose[0].pt)
					if(photon_loose[0].pt > 200): self.h_photon_passtrig_eta.Fill(photon_loose[0].eta)
		if (n_mu >= 1 and n_ele == 0):
        		self.h_mu_all.Fill(mu_loose[0].pt)
        		if (mu_loose[0].pt > 50): self.h_mu_all_eta.Fill(mu_loose[0].eta)
        		if (sigAccept_mu):
				self.h_mu_passtrig.Fill(mu_loose[0].pt)
        			if (mu_loose[0].pt > 50): self.h_mu_passtrig_eta.Fill(mu_loose[0].eta)
		if (n_mu == 0 and n_ele >= 1):
        		self.h_ele_all.Fill(ele_veto[0].pt)
        		if (ele_veto[0].pt > 40): self.h_ele_all_eta.Fill(ele_veto[0].eta)
        		if (sigAccept_ele):
				self.h_ele_passtrig.Fill(ele_veto[0].pt)
        			if (ele_veto[0].pt > 40): self.h_ele_passtrig_eta.Fill(ele_veto[0].eta)

	#mid / tight ID
	if (self.Dataset == "SingleElectron" and n_ele_mid >= 1):
		if (n_mu == 0):
        		self.h_met_all_mid.Fill(met.pt)
        		if (sigAccept_met):
				self.h_met_passtrig_mid.Fill(met.pt)
			if (n_photon >=1):
        			self.h_photon_all_mid.Fill(photon_loose[0].pt)
        			if(photon_loose[0].pt > 200): self.h_photon_all_eta_mid.Fill(photon_loose[0].eta)
        			if (sigAccept_photon):
					self.h_photon_passtrig_mid.Fill(photon_loose[0].pt)
					if(photon_loose[0].pt > 200): self.h_photon_passtrig_eta_mid.Fill(photon_loose[0].eta)
		if (n_mu_mid >= 1):
        		self.h_mu_all_mid.Fill(mu_mid[0].pt)
        		if (mu_mid[0].pt > 50): self.h_mu_all_eta_mid.Fill(mu_mid[0].eta)
        		if (sigAccept_mu):
				self.h_mu_passtrig_mid.Fill(mu_mid[0].pt)
        			if (mu_mid[0].pt > 50): self.h_mu_passtrig_eta_mid.Fill(mu_mid[0].eta)
		if (len(zmumu_mid) == 1):
			self.h_zmumu_all_mid.Fill(zmumu_mid[0].Pt())
        		if (sigAccept_mu):
				self.h_zmumu_passtrig_mid.Fill(zmumu_mid[0].Pt())

	if (self.Dataset == "SingleMuon" and n_mu_mid >= 1):
		if (n_ele == 0):
        		self.h_met_all_mid.Fill(met.pt)
        		if (sigAccept_met):
				self.h_met_passtrig_mid.Fill(met.pt)
			if (n_photon >=1):
        			self.h_photon_all_mid.Fill(photon_loose[0].pt)
        			if(photon_loose[0].pt > 200): self.h_photon_all_eta_mid.Fill(photon_loose[0].eta)
        			if (sigAccept_photon):
					self.h_photon_passtrig_mid.Fill(photon_loose[0].pt)
					if(photon_loose[0].pt > 200): self.h_photon_passtrig_eta_mid.Fill(photon_loose[0].eta)
		if (n_ele_mid >= 1):
        		self.h_ele_all_mid.Fill(ele_mid[0].pt)
        		if (ele_mid[0].pt > 40): self.h_ele_all_eta_mid.Fill(ele_mid[0].eta)
        		if (sigAccept_ele):
				self.h_ele_passtrig_mid.Fill(ele_mid[0].pt)
        			if (ele_mid[0].pt > 40): self.h_ele_passtrig_eta_mid.Fill(ele_mid[0].eta)
		if (n_ele_tight >= 1):
        		self.h_ele_all_tight.Fill(ele_tight[0].pt)
        		if (ele_tight[0].pt > 40): self.h_ele_all_eta_tight.Fill(ele_tight[0].eta)
        		if (sigAccept_ele):
				self.h_ele_passtrig_tight.Fill(ele_tight[0].pt)
        			if (ele_tight[0].pt > 40): self.h_ele_passtrig_eta_tight.Fill(ele_tight[0].eta)
		if (len(zee_mid) == 1):
			self.h_zee_all_mid.Fill(zee_mid[0].Pt())
        		if (sigAccept_ele):
				self.h_zee_passtrig_mid.Fill(zee_mid[0].Pt())

	if (self.Dataset == "SinglePhoton"):
		if (n_ele == 0 and n_mu == 0):
        		self.h_met_all_mid.Fill(met.pt)
        		if (sigAccept_met):
				self.h_met_passtrig_mid.Fill(met.pt)

	if (self.Dataset == "JetHT"):
		if (n_mu == 0 and n_ele == 0):
        		self.h_met_all_mid.Fill(met.pt)
        		if (sigAccept_met):
				self.h_met_passtrig_mid.Fill(met.pt)
			if (n_photon >=1):
        			self.h_photon_all_mid.Fill(photon_loose[0].pt)
        			if(photon_loose[0].pt > 200): self.h_photon_all_eta_mid.Fill(photon_loose[0].eta)
        			if (sigAccept_photon):
					self.h_photon_passtrig_mid.Fill(photon_loose[0].pt)
					if(photon_loose[0].pt > 200): self.h_photon_passtrig_eta_mid.Fill(photon_loose[0].eta)
		if (n_mu_mid >= 1 and n_ele == 0):
        		self.h_mu_all_mid.Fill(mu_mid[0].pt)
        		if (mu_mid[0].pt > 50): self.h_mu_all_eta_mid.Fill(mu_mid[0].eta)
        		if (sigAccept_mu):
				self.h_mu_passtrig_mid.Fill(mu_mid[0].pt)
        			if (mu_mid[0].pt > 50): self.h_mu_passtrig_eta_mid.Fill(mu_mid[0].eta)
		if (n_ele == 0 and len(zmumu_mid) == 1):
			self.h_zmumu_all_mid.Fill(zmumu_mid[0].Pt())
        		if (sigAccept_mu):
				self.h_zmumu_passtrig_mid.Fill(zmumu_mid[0].Pt())
		if (n_mu == 0 and n_ele_mid >= 1):
        		self.h_ele_all_mid.Fill(ele_mid[0].pt)
        		if (ele_mid[0].pt > 40): self.h_ele_all_eta_mid.Fill(ele_mid[0].eta)
        		if (sigAccept_ele):
				self.h_ele_passtrig_mid.Fill(ele_mid[0].pt)
        			if (ele_mid[0].pt > 40): self.h_ele_passtrig_eta_mid.Fill(ele_mid[0].eta)
		if (n_mu == 0 and n_ele_tight >= 1):
        		self.h_ele_all_tight.Fill(ele_tight[0].pt)
        		if (ele_tight[0].pt > 40): self.h_ele_all_eta_tight.Fill(ele_tight[0].eta)
        		if (sigAccept_ele):
				self.h_ele_passtrig_tight.Fill(ele_tight[0].pt)
        			if (ele_tight[0].pt > 40): self.h_ele_passtrig_eta_tight.Fill(ele_tight[0].eta)
		if (n_mu == 0 and len(zee_mid) == 1):
			self.h_zee_all_mid.Fill(zee_mid[0].Pt())
        		if (sigAccept_ele):
				self.h_zee_passtrig_mid.Fill(zee_mid[0].Pt())

	if (self.Dataset == "MET"):
		if (n_mu == 0 and n_ele == 0):
			if (n_photon >=1):
        			self.h_photon_all_mid.Fill(photon_loose[0].pt)
        			if(photon_loose[0].pt > 200): self.h_photon_all_eta_mid.Fill(photon_loose[0].eta)
        			if (sigAccept_photon):
					self.h_photon_passtrig_mid.Fill(photon_loose[0].pt)
					if(photon_loose[0].pt > 200): self.h_photon_passtrig_eta_mid.Fill(photon_loose[0].eta)
		if (n_mu_mid >= 1 and n_ele == 0):
        		self.h_mu_all_mid.Fill(mu_mid[0].pt)
        		if (mu_mid[0].pt > 50): self.h_mu_all_eta_mid.Fill(mu_mid[0].eta)
        		if (sigAccept_mu):
				self.h_mu_passtrig_mid.Fill(mu_mid[0].pt)
        			if (mu_mid[0].pt > 50): self.h_mu_passtrig_eta_mid.Fill(mu_mid[0].eta)
		if (n_ele == 0 and len(zmumu_mid) == 1):
			self.h_zmumu_all_mid.Fill(zmumu_mid[0].Pt())
        		if (sigAccept_mu):
				self.h_zmumu_passtrig_mid.Fill(zmumu_mid[0].Pt())
		if (n_mu == 0 and n_ele_mid >= 1):
        		self.h_ele_all_mid.Fill(ele_mid[0].pt)
        		if (ele_mid[0].pt > 40): self.h_ele_all_eta_mid.Fill(ele_mid[0].eta)
        		if (sigAccept_ele):
				self.h_ele_passtrig_mid.Fill(ele_mid[0].pt)
        			if (ele_mid[0].pt > 40): self.h_ele_passtrig_eta_mid.Fill(ele_mid[0].eta)
		if (n_mu == 0 and n_ele_tight >= 1):
        		self.h_ele_all_tight.Fill(ele_tight[0].pt)
        		if (ele_tight[0].pt > 40): self.h_ele_all_eta_tight.Fill(ele_tight[0].eta)
        		if (sigAccept_ele):
				self.h_ele_passtrig_tight.Fill(ele_tight[0].pt)
        			if (ele_tight[0].pt > 40): self.h_ele_passtrig_eta_tight.Fill(ele_tight[0].eta)
		if (n_mu == 0 and len(zee_mid) == 1):
			self.h_zee_all_mid.Fill(zee_mid[0].Pt())
        		if (sigAccept_ele):
				self.h_zee_passtrig_mid.Fill(zee_mid[0].Pt())

        return True

#preselection=""
#files=["root://cmseos.fnal.gov//eos/uscms/store/user/lpcsusyhad/Stop_production/Fall17_94X_Mar_2018_NanAOD_Data/SingleElectron/2017_Data_Run2017B-31Mar2018-v1/190110_061741/0000/prod2017DATA_NANO_6.root"]
#files=["root://cmseos.fnal.gov//store/group/lpcsusyhad/hui/Run2017F_SingleElectron_Nano14Dec2018_test.root"]
#p=PostProcessor(".",files,cut=preselection,branchsel=None,modules=[TrigMETAnalysis()],noOut=True,histFileName="histos_METTrigNanoAOD.root",histDirName="metTrigAnalyzerMiniAOD")
#p.run()
