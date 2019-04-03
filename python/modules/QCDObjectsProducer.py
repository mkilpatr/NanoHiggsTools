#import ROOT
#ROOT.PyConfig.IgnoreCommandLineOptions = True
#import math
#import numpy as np
#
#from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
#from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
#from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest
#
##2016 MC: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation80XReReco#Data_MC_Scale_Factors_period_dep
##2017 MC: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation94X
#
#class QCDObjectsProducer(Module):
#    def __init__(self, fileName):
#        self.metBranchName = "MET"
#	self.process       = fileName
#
#    def beginJob(self):
#        pass
#    def endJob(self):
#        pass
#
#    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
#        self.out = wrappedOutputTree
#	self.out.branch("trueResp"             , "F")
#	self.out.branch("trueRespFlv"          , "I")
#	self.out.branch("trueRespGenPT"        , "F")
#	self.out.branch("pseudoResp"           , "F")
#	self.out.branch("pseudoRespCSV"        , "F")
#	self.out.branch("pseudoRespPseudoGenPT", "F")
#	self.out.branch("pseudoRespPassFilter" , "O")
#
#    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
#        pass
#
#    def addFourVec(self, obj1, obj2):
#	tot = ROOT.TLorentzVector()
#	v1 = ROOT.TLorentzVector()
#	v2 = ROOT.TLorentzVector()
#	v1.SetPtEtaPhiM(obj1.pt, 0, obj1.phi, 0)
#	v2.SetPtEtaPhiM(obj2.pt, obj2.eta, obj2.phi, obj2.mass)
#	tot = (v1 + v2)
#	return tot
#
#    def analyze(self, event):
#        """process event, return True (go to next module) or False (fail, go to next event)"""
#        ## Getting objects
#	jets	  = Collection(event, "Jet")
#	genjets   = Collection(event, "GenJet")
#        met       = Object(event, self.metBranchName)
#
#	jetNearMETInd, MMJetDPhi = -1
#	for iJ in xrange(len(jets)):
#		if iJ > 2 : continue
#		dPhi = deltaPhi(jets[iJ].eta, jets[iJ].phi, met.eta, met.phi)
#		if(MMJetDPhi < 0 or dPhi < MMJetDPhi):
#			MMJetDPhi = dPhi
#			jetNearMETInd = iJ
#
#	if jetNearMETInd < 0 : return True
#
#	pJ = jets[jetNearMETInd]
#	passFilter = len(jets) > jetNearMETInd
#	pseudoGenPT = self.addFourVec(met, pJ).Pt()
#	MMPseudoResp = pJ.pt/pseudoGenPT if pseudoGenPT > 0 else 999
#
#	# True response info
#	trueRespInd = jetAndMETCorrections.getQCDRespTailCorrector()->mmInd if QCD else -1
#	trueRespFlv = 99
#	trueRespGenPT = -1.0
#	if trueRespInd >= 0:
#		for iG in xrange(len(genjets)):
#			gjet = genjets[iG]
#			if iG != trueRespInd: continue
#			trueRespGenPT = gjet.pt
#			trueRespFlv = gjet.partonFlavor
#			break
#	
#	self.out.fillBranch("trueResp"     	   ,jetAndMETCorrections.getQCDRespTailCorrector()->mmResp if trueRespInd < 0 else -1
#	self.out.fillBranch("trueRespFlv"  	   ,trueRespFlv);
#	self.out.fillBranch("trueRespGenPT"	   ,trueRespGenPT);
#
#
#        ### Store output
#	self.out.fillBranch("pseudoResp"           , MMPseudoResp)
#	self.out.fillBranch("pseudoRespCSV"        , pJ.btagDeepB)
#	self.out.fillBranch("pseudoRespPseudoGenPT", pseudoGenPT)
#	self.out.fillBranch("pseudoRespPassFilter" , passFilter)
#	return True
#
#
# # define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
