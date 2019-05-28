import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import math
import numpy as np

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest

#2016 MC: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation80XReReco#Data_MC_Scale_Factors_period_dep
#2017 MC: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation94X

class QCDObjectsProducer(Module):
    def __init__(self, isQCD = False, isData = False, isQCDOrig = False):
        self.metBranchName = "MET"
	self.isQCD       = isQCD
	self.isData	 = isData
	self.isQCDOrig	 = isQCDOrig
	self.nBootstraps = 50

    def beginJob(self):
        pass
    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
	self.out.branch("trueResp"             , "F")
	self.out.branch("trueRespFlv"          , "I")
	self.out.branch("trueRespGenPT"        , "F")
	self.out.branch("pseudoResp"           , "F")
	self.out.branch("pseudoRespCSV"        , "F")
	self.out.branch("pseudoRespPseudoGenPT", "F")
	self.out.branch("pseudoRespPassFilter" , "O")
	if self.isQCDOrig:
		self.out.branch("nBootstrapWeight",        "I")
		self.out.branch("bootstrapWeight",         "I", lenVar="nBootstrapWeight")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def addFourVec(self, obj1, obj2):
	tot = ROOT.TLorentzVector()
	v1 = ROOT.TLorentzVector()
	v2 = ROOT.TLorentzVector()
	v1.SetPtEtaPhiM(obj1.pt, 0, obj1.phi, 0)
	v2.SetPtEtaPhiM(obj2.pt, obj2.eta, obj2.phi, obj2.mass)
	tot = (v1 + v2)
	return tot

    def getQCDRespTailCorrector(self, jets, genJets, met):
    
      	MM = -1
      	ind = -1
      	flv = -1
      	resp = -1
	mmout = []
	for iG in xrange(len(genJets)):
		if iG > 2: break
        	if genJets[iG].pt == 0: break
        	fpt = -1
		for rJ in xrange(len(jets)):
			if jets[rJ].genJetIdx == iG:
				fpt = jets[rJ].pt
				break
        	if fpt < 0: fpt = 9.5 #for the cases with no reco jet due to being below thresh
        	if(MM < 0 or abs(fpt - genJets[iG].pt) > MM):
			ind = iG
			resp =  fpt/genJets[iG].pt
			MM = abs(fpt - genJets[iG].pt)
			flv = genJets[iG].partonFlavour
    
		if ind >= 0:
			mmResp = resp
			mmInd = ind
			mmFlv = flv
      		else:
			mmResp = -1
			mmInd = -1
			mmFlv = -1

	mmout = [mmInd, mmResp, mmFlv]
	return mmout     		

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        ## Getting objects
	jets	  = Collection(event, "Jet")
	if self.isData == False: 
		genjets = Collection(event, "GenJet")
        met       = Object(event, self.metBranchName)

	jetNearMETInd, MMJetDPhi = -1, -1
	for iJ in xrange(len(jets)):
		if iJ > 2 : continue
		dPhi = abs(deltaPhi(jets[iJ].phi, met.phi))
		if(MMJetDPhi < 0 or dPhi < MMJetDPhi):
			MMJetDPhi = dPhi
			jetNearMETInd = iJ

	if jetNearMETInd < 0 : return True

	pJ = jets[jetNearMETInd]
	passFilter = len(jets) > jetNearMETInd
	pseudoGenPT = self.addFourVec(met, pJ).Pt()
	MMPseudoResp = pJ.pt/pseudoGenPT if pseudoGenPT > 0 else 999

	# True response info
	#print "isQCD: ", self.isQCD
	mmOut = []
	if self.isQCD == True:
		mmOut = self.getQCDRespTailCorrector(jets, genjets, met) 
	else:
		mmOut = [-1, -1.0, -1]
	trueRespInd, trueResp = mmOut[0], mmOut[1]
	#print "trueResp: ", trueRespInd, trueResp
	trueRespFlv = 99
	trueRespGenPT = -1.0
	if trueRespInd >= 0:
		for iG in xrange(len(genjets)):
			gjet = genjets[iG]
			if iG != trueRespInd: continue
			trueRespGenPT = gjet.pt
			trueRespFlv = gjet.partonFlavour
			break
	
	if self.isQCDOrig:
		b = []
		for iB in xrange(self.nBootstraps):
		        b.append(1)
		self.out.fillBranch("nBootstrapWeight",        self.nBootstraps)
		self.out.fillBranch("bootstrapWeight",         b)
	
        ### Store output
	self.out.fillBranch("pseudoResp"           , MMPseudoResp)
	self.out.fillBranch("pseudoRespCSV"        , pJ.btagDeepB)
	self.out.fillBranch("pseudoRespPseudoGenPT", pseudoGenPT)
	self.out.fillBranch("pseudoRespPassFilter" , passFilter)
	self.out.fillBranch("trueResp"     	   , trueResp if trueRespInd >= 0 else -1)
	self.out.fillBranch("trueRespFlv"  	   , trueRespFlv)
	self.out.fillBranch("trueRespGenPT"	   , trueRespGenPT)
	return True


 # define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
