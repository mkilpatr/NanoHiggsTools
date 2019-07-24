import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import math
import numpy as np

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest

#2016 MC: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation80XReReco#Data_MC_Scale_Factors_period_dep
#2017 MC: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation94X

DeepCSVLooseWP = {
    "2016" : 0.2217,
    "2017" : 0.1522,
    "2018" : 0.1241
}

DeepCSVMediumWP ={
    "2016" : 0.6324,
    "2017" : 0.4941,
    "2018" : 0.4184
}

CSVv2MediumWP = {
    "2016" : 0.8484,
    "2017" : 0.8838,
    "2018" : 0.8838  # Not recommended, use 2017 as temp
}


class LLObjectsProducer(Module):
    def __init__(self, era, isData = False):
        self.era = era
	self.isData = isData
        self.metBranchName = "MET"
        # EE noise mitigation in PF MET
        # https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/1865.html
        if self.era == "2017":
            self.metBranchName = "METFixEE2017"

    def beginJob(self):
        pass
    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
	self.out.branch("Stop0l_nJets_HighPt",  		"I")
	self.out.branch("Stop0l_nJets30",  			"I")
	self.out.branch("Stop0l_nbtags_Loose",  		"I")
	self.out.branch("Stop0l_MtLepMET", 			"F")
	self.out.branch("nLeptonVeto",    			"I")
	self.out.branch("Stop0l_nIsoTracksLep", 		"I")
	self.out.branch("Stop0l_nIsoTracksHad", 		"I")
	self.out.branch("Stop0l_nIsoTracksHad_ptgeq20", 	"I")
	self.out.branch("Stop0l_nVetoElec", 			"I")
	self.out.branch("Stop0l_nVetoMuon", 			"I")
	self.out.branch("Stop0l_nVetoElecMuon", 		"I")
	self.out.branch("Stop0l_noMuonJet",			"O")
	self.out.branch("Pass_dPhiMETMedDM", 			"O")
	self.out.branch("Pass_dPhiQCD",				"O")
	self.out.branch("Stop0l_dPhiISRMET",			"F")
	self.out.branch("Stop0l_TauPOG",			"I")
	self.out.branch("Pass_HT30", 				"O")
	self.out.branch("Pass_dPhiMET30", 			"O")
	self.out.branch("Pass_dPhiMETLowDM30", 			"O")
	self.out.branch("Pass_dPhiMETHighDM30", 		"O")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def SelMtlepMET(self, ele, muon, isks, met):
	mt = 0.0
	for l in ele:
		if l.Stop0l: mt += math.sqrt( 2 * l.pt * met.pt * (1 - np.cos(deltaPhi(l.phi,met.phi))))
	for l in muon:
		if l.Stop0l: mt += math.sqrt( 2 * l.pt * met.pt * (1 - np.cos(deltaPhi(l.phi,met.phi))))
	for l in isks:
		if l.Stop0l: mt += math.sqrt( 2 * l.pt * met.pt * (1 - np.cos(deltaPhi(l.phi,met.phi))))
	return mt

    def SelJets(self, jet):
        if jet.pt < 20 or math.fabs(jet.eta) > 2.4 :
            return False
        return True

    def SelJets30(self, jet):
        if jet.pt < 30 or math.fabs(jet.eta) > 2.4 :
            return False
        return True

    def SelJetsHighPt(self, jet):
        if jet.pt < 75 or math.fabs(jet.eta) > 2.4 :
            return False
        return True

    def SelBtagJets(self, jet):
        global DeepCSVLooseWP
        if jet.btagDeepB >= DeepCSVLooseWP[self.era]:
            return True
        return False

    def GetJetSortedIdxVal(self, jets):
        ptlist = []
	etalist = []
        dphiMET = []
        for j in jets:
            if math.fabs(j.eta) > 4.7 or j.pt < 20:
                pass
            else:
                ptlist.append(j.pt)
		etalist.append(math.fabs(j.eta))
                dphiMET.append(j.dPhiMET)
	output = []
	for j in np.argsort(ptlist)[::-1]:
		if j != len(ptlist) - 1:
			if ptlist[j] == ptlist[j + 1]:
				#print "j: ", j
				#print "ptList: ", ptlist
				#print "etalist: ", etalist
				if etalist[j] <= etalist[j + 1]:
					output.append(dphiMET[j])
				else:
					output.append(dphiMET[j + 1])
		else:
			output.append(dphiMET[j])
	
	return output


    def GetJetSortedIdx(self, jets, jetpt = 20):
        ptlist = []
        dphiMET = []
        for j in jets:
            if math.fabs(j.eta) > 4.7 or j.pt < jetpt:
                pass
            else:
                ptlist.append(j.pt)
                dphiMET.append(j.dPhiMET)
        return [dphiMET[j] for j in np.argsort(ptlist)[::-1]]

    def PassdPhiVal(self, sortedPhi, dPhiCutsLow, dPhiCutsHigh):
	return all( (a < b and b < c) for a, b, c in zip(dPhiCutsLow, sortedPhi, dPhiCutsHigh))


    def PassdPhi(self, sortedPhi, dPhiCuts, invertdPhi =False):
        if invertdPhi:
            return any( a < b for a, b in zip(sortedPhi, dPhiCuts))
        else:
            return all( a > b for a, b in zip(sortedPhi, dPhiCuts))

    def SelTauPOG(self, tau):
	if tau.pt < 20 or abs(tau.eta) > 2.4 or not tau.idDecayMode or not (tau.idMVAoldDM2017v2 & 8):
		return False
	return True

    def SelGenTau(self, gentau):
	if gentau.pt < 10 or abs(gentau.eta) > 2.4:
		return False
	return True

    def CalHT(self, jets, jetpt):
	HT = sum([j.pt for i, j in enumerate(jets) if (self.Jet_Stop0l[i] and j.pt > jetpt)])
	return HT

    def SelNoMuon(self, jets, met):
	noMuonJet = True
	for j in jets:
		if j.pt > 200 and j.muEF > 0.5 and abs(deltaPhi(j.phi, met.phi)) > (math.pi - 0.4):
			noMuonJet = False
	return noMuonJet

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        ## Getting objects
        electrons = Collection(event, "Electron")
        muons     = Collection(event, "Muon")
        isotracks = Collection(event, "IsoTrack")
	muons     = Collection(event, "Muon")
	electrons = Collection(event, "Electron")
	jets      = Collection(event, "Jet")
	met       = Object(event, self.metBranchName)
	#tau	  = Collection(event, "Tau")
	if not self.isData:
		gentau	  = Collection(event, "GenVisTau")
	stop0l    = Object(event, "Stop0l")
	fatjets   = Collection(event, "FatJet")

        ## Selecting objects
	self.Jet_Stop0lHighPt= map(self.SelJetsHighPt, jets)
	self.Jet_Stop0l      = map(self.SelJets, jets)
	self.Jet_Stop0l30    = map(self.SelJets30, jets)
	local_BJet_Stop0l    = map(self.SelBtagJets, jets)
        self.BJet_Stop0l     = [a and b for a, b in zip(self.Jet_Stop0l, local_BJet_Stop0l )]
	mt 		     = self.SelMtlepMET(electrons, muons, isotracks, met)
	countIskLep 	     = sum([(i.Stop0l and (abs(i.pdgId) == 11 or abs(i.pdgId) == 13)) for i in isotracks])
	countIskHad 	     = sum([(i.Stop0l and abs(i.pdgId) == 211) for i in isotracks])
	countIskHad_ptgeq20  = sum([(i.Stop0l and abs(i.pdgId) == 211 and i.pt > 20) for i in isotracks])
	countEle	     = sum([e.Stop0l for e in electrons])
	countMuon	     = sum([m.Stop0l for m in muons])
	noMuonJet	     = self.SelNoMuon(jets, met)
	sortedPhiVal         = self.GetJetSortedIdxVal(jets)
	PassdPhiMedDM        = self.PassdPhiVal(sortedPhiVal, [0.15, 0.15, 0.15], [0.5, 4., 4.]) #Variable for LowDM Validation bins
	sortedPhi	     = self.GetJetSortedIdx(jets)
	PassdPhiQCD          = self.PassdPhi(sortedPhi, [0.1, 0.1, 0.1], invertdPhi =True)
	dphiISRMet	     = abs(deltaPhi(fatjets[stop0l.ISRJetIdx].phi, met.phi)) if stop0l.ISRJetIdx >= 0 else -1
	#self.Tau_Stop0l     = map(self.SelTauPOG, tau)
	#countTauPOG	     = sum(self.Tau_Stop0l)
	HT 		     = self.CalHT(jets, 30)
	PassHT30	     = HT >= 300
	sortedPhi30 	     = self.GetJetSortedIdx(jets, 30)
	PassdPhiLowDM30      = self.PassdPhi(sortedPhi30, [0.5, 0.15, 0.15])
	PassdPhiHighDM30     = self.PassdPhi(sortedPhi30, [0.5, 0.5, 0.5, 0.5])

        ### Store output
	self.out.fillBranch("Stop0l_nJets_HighPt",    	sum(self.Jet_Stop0lHighPt))
	self.out.fillBranch("Stop0l_nJets30",    	sum(self.Jet_Stop0l30))
	self.out.fillBranch("Stop0l_nbtags_Loose",   	sum(self.BJet_Stop0l))
	self.out.fillBranch("Stop0l_MtLepMET",  	mt)
	self.out.fillBranch("nLeptonVeto",    		countMuon + countEle + countIskLep)
	self.out.fillBranch("Stop0l_nIsoTracksLep",	countIskLep)
	self.out.fillBranch("Stop0l_nIsoTracksHad",	countIskHad)
	self.out.fillBranch("Stop0l_nIsoTracksHad_ptgeq20",	countIskHad_ptgeq20)
	self.out.fillBranch("Stop0l_nVetoElec", 	countEle)
	self.out.fillBranch("Stop0l_nVetoMuon", 	countMuon)
	self.out.fillBranch("Stop0l_nVetoElecMuon", 	countEle + countMuon)
	self.out.fillBranch("Stop0l_noMuonJet",		noMuonJet)
	self.out.fillBranch("Pass_dPhiMETMedDM", 	PassdPhiMedDM)
	self.out.fillBranch("Pass_dPhiQCD",		PassdPhiQCD)
	self.out.fillBranch("Stop0l_dPhiISRMET",	dphiISRMet)
	#self.out.fillBranch("Stop0l_TauPOG",		countTauPOG)
	self.out.fillBranch("Pass_HT30",		PassHT30)
	self.out.fillBranch("Pass_dPhiMET30", 		PassdPhiLowDM30)
	self.out.fillBranch("Pass_dPhiMETLowDM30", 	PassdPhiLowDM30)
	self.out.fillBranch("Pass_dPhiMETHighDM30", 	PassdPhiHighDM30)
	return True


 # define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
