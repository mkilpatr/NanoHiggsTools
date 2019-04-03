import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import math
import numpy as np

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest

#2016 MC: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation80XReReco#Data_MC_Scale_Factors_period_dep
#2017 MC: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation94X

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
    def __init__(self, era):
        self.era = era
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
	self.out.branch("nStop0l_MtLepMET", "I")
	self.out.branch("Stop0l_MtLepMET", "F",  lenVar="nStop0l_MtLepMET")
	self.out.branch("Jet_btagStop0l_pt1", "F")
	self.out.branch("Jet_btagStop0l_pt2", "F")
	self.out.branch("nLeptonVeto",    "I")
	self.out.branch("Pass_Trigger",   "O")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def PassLeptonVeto(self, eles, muons, isks, mva = False):
        countEle = sum([e.Stop0l for e in eles])
        countMu  = sum([m.Stop0l for m in muons])
        countIsk = sum([i.Stop0l for i in isks])
        if mva: return countEle + countMu
        else: return countEle + countMu + countIsk

    def SelMtlepMET(self, ele, muon, isks, met):
	mt = []
	for l in ele:
		if l.Stop0l: mt.append(math.sqrt( 2 * l.pt * met.pt * (1 - np.cos(deltaPhi(l.phi,met.phi)))))
	for l in muon:
		if l.Stop0l: mt.append(math.sqrt( 2 * l.pt * met.pt * (1 - np.cos(deltaPhi(l.phi,met.phi)))))
	for l in isks:
		if l.Stop0l: mt.append(math.sqrt( 2 * l.pt * met.pt * (1 - np.cos(deltaPhi(l.phi,met.phi)))))
	#print mt
	return mt

    def SelJets(self, jet):
        if jet.pt < 20 or math.fabs(jet.eta) > 2.4 :
            return False
        return True

    def SelBtagJets(self, jet):
        global DeepCSVMediumWP
        if jet.btagDeepB < DeepCSVMediumWP[self.era]:
            return False
        return True

    def CalMTbPTb(self, jets, met):
        Bjetpt = []

        # Getting bjet, ordered by pt
        bjets = [ j for i,j in enumerate(jets) if self.BJet_Stop0l[i]]
        # Getting btag index, ordered by b discriminator value
        btagidx = sorted(range(len(bjets)), key=lambda k: bjets[k].btagDeepB , reverse=True)

        for i in range(min(len(btagidx), 2)):
            bj = bjets[btagidx[i]]
            Bjetpt.append(bj.pt)
            if len(btagidx) == 1: Bjetpt.append(0.0)

        if len(btagidx) == 0:
                Bjetpt.append(0.0)
                Bjetpt.append(0.0)

        return Bjetpt

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        ## Getting objects
        electrons = Collection(event, "Electron")
        muons     = Collection(event, "Muon")
        isotracks = Collection(event, "IsoTrack")
	jets      = Collection(event, "Jet")
	met       = Object(event, self.metBranchName)
	hlt       = Object(event, "HLT")

        ## Selecting objects
	self.Jet_Stop0l      = map(self.SelJets, jets)
	local_BJet_Stop0l    = map(self.SelBtagJets, jets)
	self.BJet_Stop0l     = [a and b for a, b in zip(self.Jet_Stop0l, local_BJet_Stop0l )]
	bJetPt 		     = self.CalMTbPTb(jets, met)
	mt 		     = self.SelMtlepMET(electrons, muons, isotracks, met)
	PassLeptonVeto       = self.PassLeptonVeto(electrons, muons, isotracks)

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
	
        ### Store output
	self.out.fillBranch("nStop0l_MtLepMET", len(mt))
	self.out.fillBranch("Stop0l_MtLepMET",  mt)
	self.out.fillBranch("Jet_btagStop0l_pt1", bJetPt[0])
	self.out.fillBranch("Jet_btagStop0l_pt2", bJetPt[1])
	self.out.fillBranch("nLeptonVeto",    PassLeptonVeto)
	self.out.fillBranch("Pass_Trigger", sigAccept_met)
	return True


 # define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
