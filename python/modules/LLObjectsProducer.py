import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import math
import numpy as np
from functools import reduce
import operator

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
    def __init__(self, era, isData = False, applyUncert=None):
        self.era = era
	self.isData = isData
        self.metBranchName = "MET"
        # EE noise mitigation in PF MET
        # https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/1865.html
        if self.era == "2017":
            self.metBranchName = "METFixEE2017"

	self.applyUncert = applyUncert
	self.suffix = ""
	
	if self.applyUncert == "JESUp":
	    self.suffix = "_JESUp"
	elif self.applyUncert == "JESDown":
	    self.suffix = "_JESDown"

    def beginJob(self):
        pass
    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
	self.out.branch("Stop0l_MtLepMET", 			"F")
	self.out.branch("Stop0l_nVetoElecMuon", 		"I")
	self.out.branch("Stop0l_noMuonJet",			"O")
	self.out.branch("Pass_dPhiQCD",				"O")
	self.out.branch("Pass_dPhiQCDSF",			"O")
	self.out.branch("Stop0l_dPhiISRMET",			"F")
	self.out.branch("Pass_HT30",				"O")
	self.out.branch("Pass_NJets30", 			"O")
	self.out.branch("Pass_dPhiMET30", 			"O")
        self.out.branch("Pass_dPhiMETLowDM30", 			"O")
        self.out.branch("Pass_dPhiMETMedDM30", 			"O")
        self.out.branch("Pass_dPhiMETHighDM30", 		"O")
	if not self.isData:
		self.out.branch("ElectronMedSF",		"F")
		self.out.branch("ElectronMedSFErr",		"F")
		self.out.branch("ElectronVetoSF",		"F")
		self.out.branch("ElectronVetoSFErr",		"F")
		self.out.branch("MuonLooseSF",			"F")
		self.out.branch("MuonLooseSFErr",		"F")
		self.out.branch("MuonMedSF",			"F")
		self.out.branch("MuonMedSFErr",			"F")
		self.out.branch("TauSF",			"F")
		self.out.branch("TauSF_Up",			"F")
		self.out.branch("TauSF_Down",			"F")
		self.out.branch("WtagSF",			"F")
		self.out.branch("WtagSFErr",			"F")
		self.out.branch("TopSF",			"F")
		self.out.branch("TopSFErr",			"F")
		self.out.branch("restopSF",			"F")
		self.out.branch("restopSF_Up",			"F")
		self.out.branch("restopSF_Down",		"F")
		self.out.branch("SoftBSF",			"F")
		self.out.branch("SoftBSFErr",			"F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def GetJetSortedIdx(self, jets, jetpt = 20, jeteta = 4.7):
        ptlist = []
	etalist = []
        dphiMET = []
        for j in jets:
            if math.fabs(j.eta) > jeteta or j.pt < jetpt:
                pass
            else:
		ptlist.append(-j.pt)
		etalist.append(math.fabs(j.eta))
                dphiMET.append(j.dPhiMET)

	sortIdx = np.lexsort((etalist, ptlist))

	return sortIdx, [dphiMET[j] for j in sortIdx]

    def PassdPhi(self, sortedPhi, dPhiCuts, invertdPhi =False):
        if invertdPhi:
            return any( a < b for a, b in zip(sortedPhi, dPhiCuts))
        else:
            return all( a > b for a, b in zip(sortedPhi, dPhiCuts))

    def PassdPhiVal(self, sortedPhi, dPhiCutsLow, dPhiCutsHigh):
            return all( (a < b and b < c) for a, b, c in zip(dPhiCutsLow, sortedPhi, dPhiCutsHigh))

    def SelNoMuon(self, jets, met):
	noMuonJet = True
	for j in jets:
		if j.pt > 200 and j.muEF > 0.5 and abs(deltaPhi(j.phi, met.phi)) > (math.pi - 0.4):
			noMuonJet = False
	return noMuonJet

    def CalHT(self, jets, jetpt = 20.):
        HT = sum([j.pt for j in jets if (j.Stop0l and j.pt >= jetpt)])
        return HT

    def PassNjets(self, jets, jetpt = 20.):
        countJets = sum([j.Stop0l for j in jets if j.pt >= jetpt])
        return countJets >= 2

    def ScaleFactorErrElectron(self, obj, kind="Medium"):
	sf = 1
	sfErr = 0
	for s in obj:
		if not s.Stop0l: continue
		if kind == "Medium":
			sf *= s.MediumSF
			sfErr += ((s.MediumSFErr)**2)
		elif kind == "Veto":
			sf *= s.VetoSF
			sfErr += ((s.VetoSFErr)**2)

	return sf, math.sqrt(sfErr)

    def ScaleFactorErrMuon(self, obj, kind="Loose"):
	sf = 1
	sfErr = 0
	for s in obj:
		if not s.Stop0l: continue
		if kind == "Loose":
			sf *= s.LooseSF
			sfErr += ((s.LooseSFErr)**2)
		elif kind == "Medium":
			sf *= s.MediumSF
			sfErr += ((s.MediumSFErr)**2)

	return sf, math.sqrt(sfErr)

    def ScaleFactorErrTau(self, obj):
	sf = 1
	sfUp = 0
	sfDown = 0
	for s in obj:
		if not s.Stop0l: continue
		sf *= s.MediumSF
		sfUp += ((s.MediumSF_Up)**2)
		sfDown += ((s.MediumSF_Down)**2)

	return sf, math.sqrt(sfUp), math.sqrt(sfDown)

    def ScaleFactorErrFatjet(self, obj, objType):
	sf = 1
	sfErr = 0
	for s in obj:
		if not (s.Stop0l == objType): continue
		sf *= s.SF
		sfErr += ((s.SFerr)**2)

	return sf, math.sqrt(sfErr)

    def ScaleFactorErrResTop(self, obj, objType):
	sf = 1
	sfUp = 0
	sfDown = 0
	for s, rT in zip(obj, objType):
		if not rT.Stop0l: continue
		sf *= s.sf
		sfUp += ((s.syst_Btag_Up)**2 + (s.syst_Pileup_Up)**2 + (s.syst_CSPur_Up)**2 + (s.syst_Stat_Up)**2)
		sfDown += ((s.syst_Btag_Down)**2 + (s.syst_Pileup_Down)**2 + (s.syst_CSPur_Down)**2 + (s.syst_Stat_Down)**2)
		if not s.genMatch:
			sfUp += (s.syst_Closure_Up)**2
			sfDown += (s.syst_Closure_Down)**2

	return sf, math.sqrt(sfUp), math.sqrt(sfDown)

    def ScaleFactorErrSoftB(self, obj):
	sf = 1
	sfErr = 0
	for s in obj:
		if not s.Stop0l: continue
		sf *= s.SF
		sfErr += ((s.SFerr)**2)

	return sf, math.sqrt(sfErr)

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        ## Getting objects
        electrons = Collection(event, "Electron")
        muons     = Collection(event, "Muon")
        isotracks = Collection(event, "IsoTrack")
	taus	  = Collection(event, "Tau")
	stop0l    = Object(event, "Stop0l")
	fatjets   = Collection(event, "FatJet")
	SB	  = Collection(event, "SB")

	if self.applyUncert == "JESUp":
	    restop    = Collection(event, "ResolvedTopCandidate_JESUp")
	    res	      = Collection(event, "ResolvedTop_JESUp", lenVar="nResolvedTopCandidate_JESUp")
	    jets      = CollectionRemapped(event, "Jet", replaceMap={"pt":"pt_jesTotalUp", "mass":"mass_jesTotalUp"})
	    met       = ObjectRemapped(event, self.metBranchName, replaceMap={"pt":"pt_jesTotalUp", "phi":"phi_jesTotalUp"})
	elif self.applyUncert == "JESDown":
	    restop    = Collection(event, "ResolvedTopCandidate_JESDown")
	    res	      = Collection(event, "ResolvedTop_JESDown", lenVar="nResolvedTopCandidate_JESDown")
	    jets      = CollectionRemapped(event, "Jet", replaceMap={"pt":"pt_jesTotalDown", "mass":"mass_jesTotalDown"})
	    met       = ObjectRemapped(event, self.metBranchName, replaceMap={"pt":"pt_jesTotalDown", "phi":"phi_jesTotalDown"})
	else:
	    restop    = Collection(event, "ResolvedTopCandidate")
	    res	      = Collection(event, "ResolvedTop", lenVar="nResolvedTopCandidate")
	    jets      = Collection(event, "Jet")
	    met       = Object(event, self.metBranchName)


        ## Selecting objects
	mt		     = sum([ e.MtW for e in electrons if e.Stop0l ] + [ m.MtW for m in muons if m.Stop0l ])
	countEle	     = sum([e.Stop0l for e in electrons])
	countMuon	     = sum([m.Stop0l for m in muons])
	noMuonJet	     = self.SelNoMuon(jets, met)
	sortedIdx, sortedPhi = self.GetJetSortedIdx(jets)
	PassdPhiQCD          = self.PassdPhi(sortedPhi, [0.1, 0.1, 0.1], invertdPhi =True)
	PassdPhiQCDSF        = self.PassdPhi(sortedPhi, [0.1, 0.1], invertdPhi =True)
	dphiISRMet	     = abs(deltaPhi(fatjets[stop0l.ISRJetIdx].phi, met.phi)) if stop0l.ISRJetIdx >= 0 else -1
	
	HT30		     = self.CalHT(jets, 30.) >= 300
	PassNjets30	     = self.PassNjets(jets, 30.)
	sortedIdx, sortedPhi = self.GetJetSortedIdx(jets, 30.)
        PassdPhiLowDM30      = self.PassdPhi(sortedPhi, [0.5, 0.15, 0.15])
	PassdPhiMedDM30      = self.PassdPhiVal(sortedPhi, [0.15, 0.15, 0.15], [0.5, 4., 4.]) #Variable for LowDM Validation bins
        PassdPhiHighDM30     = self.PassdPhi(sortedPhi, [0.5, 0.5, 0.5, 0.5])
        PassdPhiQCD30        = self.PassdPhi(sortedPhi, [0.1, 0.1, 0.1], invertdPhi =True)

	if not self.isData:
		electronMedSF, electronMedSFErr = self.ScaleFactorErrElectron(electrons, "Medium")
		electronVetoSF, electronVetoSFErr = self.ScaleFactorErrElectron(electrons, "Veto")
		muonLooseSF, muonLooseSFErr    	= self.ScaleFactorErrMuon(muons, "Loose")
		muonMedSF, muonMedSFErr    	= self.ScaleFactorErrMuon(muons, "Medium")
		tauSF, tauSFUp, tauSFDown 	= self.ScaleFactorErrTau(taus)
		## type top = 1, W = 2, else 0
		WtagSF, WtagSFErr	     	= self.ScaleFactorErrFatjet(fatjets, 2)
		topSF, topSFErr	     		= self.ScaleFactorErrFatjet(fatjets, 1)
		resSF, resSFUp, resSFDown 	= self.ScaleFactorErrResTop(restop, res)
		softBSF, softBSFErr    		= self.ScaleFactorErrSoftB(SB)

        ### Store output
	self.out.fillBranch("Stop0l_MtLepMET",  	mt)
	self.out.fillBranch("Stop0l_nVetoElecMuon", 	countEle + countMuon)
	self.out.fillBranch("Stop0l_noMuonJet",		noMuonJet)
	self.out.fillBranch("Pass_dPhiQCD",		PassdPhiQCD)
	self.out.fillBranch("Pass_dPhiQCDSF",		PassdPhiQCDSF)
	self.out.fillBranch("Stop0l_dPhiISRMET",	dphiISRMet)
	self.out.fillBranch("Pass_HT30",		HT30)
	self.out.fillBranch("Pass_NJets30", 		PassNjets30)
	self.out.fillBranch("Pass_dPhiMET30", 		PassdPhiLowDM30)
        self.out.fillBranch("Pass_dPhiMETLowDM30", 	PassdPhiLowDM30)
        self.out.fillBranch("Pass_dPhiMETMedDM30", 	PassdPhiMedDM30)
        self.out.fillBranch("Pass_dPhiMETHighDM30", 	PassdPhiHighDM30)
	
	if not self.isData:
		self.out.fillBranch("ElectronMedSF",	electronMedSF)
		self.out.fillBranch("ElectronMedSFErr",	electronMedSFErr)
		self.out.fillBranch("ElectronVetoSF",	electronVetoSF)
		self.out.fillBranch("ElectronVetoSFErr",electronVetoSFErr)
		self.out.fillBranch("MuonLooseSF",	muonLooseSF)
		self.out.fillBranch("MuonLooseSFErr",	muonLooseSFErr)
		self.out.fillBranch("MuonMedSF",	muonMedSF)
		self.out.fillBranch("MuonMedSFErr",	muonMedSFErr)
		self.out.fillBranch("TauSF",		tauSF)
		self.out.fillBranch("TauSF_Up",		tauSFUp)
		self.out.fillBranch("TauSF_Down",	tauSFDown)
		self.out.fillBranch("WtagSF",		WtagSF)
		self.out.fillBranch("WtagSFErr",	WtagSFErr)
		self.out.fillBranch("TopSF",		topSF)
		self.out.fillBranch("TopSFErr",		topSFErr)
		self.out.fillBranch("restopSF",		resSF)
		self.out.fillBranch("restopSF_Up",	resSFUp)
		self.out.fillBranch("restopSF_Down",	resSFDown)
		self.out.fillBranch("SoftBSF",		softBSF)
		self.out.fillBranch("SoftBSFErr",	softBSFErr)
	return True


 # define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
