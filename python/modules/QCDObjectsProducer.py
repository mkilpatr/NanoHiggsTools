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
    def __init__(self, isQCD = False, isData = False):
        self.metBranchName = "MET"
	self.isQCD       = isQCD
	self.isData	 = isData

    def beginJob(self):
        pass
    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
	self.out.branch("trueResp"             , "F")
	self.out.branch("trueRespFlv"          , "I")
	self.out.branch("trueRespFlvReco"      , "I")
	self.out.branch("trueRespGenPT"        , "F")
	self.out.branch("trueRespRecoPT"       , "F")
	self.out.branch("trueRespMM" 	       , "F")
	self.out.branch("trueRespGenIdx"       , "I")
	self.out.branch("trueRespRecoIdx"      , "I")
	self.out.branch("pseudoResp"           , "F")
	self.out.branch("pseudoRespCSV"        , "F")
	self.out.branch("pseudoRespPseudoGenPT", "F")
	self.out.branch("pseudoRespPassFilter" , "O")

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
    
      	MM = 0
      	ind = -1
      	flv = -1
      	recoflv = -1
      	resp = -1
	mmout = []
	recoInd = -1
	recoInd_ = -1
	recopt = -1
	for iG in range(0,len(genJets)):
		#if iG == 4: break
        	if genJets[iG].pt == 0: break
        	fpt = -1

		for rJ in xrange(len(jets)):
			if not jets[rJ].Stop0l: continue
			if jets[rJ].genJetIdx != iG: continue
			fpt = jets[rJ].pt
			#fpt = jets[rJ].pt*(1-jets[rJ].rawFactor)
			recoInd = rJ
			break

        	if fpt < 0: fpt = 9.5 #for the cases with no reco jet due to being below thresh

		#print("MM: %f, diff: %f, fpt: %f, genpt: %f" % (MM, abs(fpt - genJets[iG].pt), fpt, genJets[iG].pt))
        	if abs(fpt - genJets[iG].pt) > MM:
			ind = iG #Could have some error here. Need to find out what the function index() does in UCSB code.
			recoInd_ = recoInd
			resp =  float(fpt)/float(genJets[iG].pt)
			MM = abs(float(fpt - genJets[iG].pt))
			flv = abs(genJets[iG].hadronFlavour)
			if recoInd != -1: recoflv = abs(jets[recoInd].hadronFlavour)
			recopt = fpt
   
	#if flv == 5 and resp > 1: print("Jet index: %i, resp: %f, MM: %f, flv: %i, fpt: %f, genpt: %f, recopt: %f, MM_manual: %f" % (ind, resp, MM, flv, fpt, genJets[ind].pt, jets[recoInd_].pt))
	if ind >= 0:
		mmResp = resp
		mmInd = ind
		mmFlv = flv
		mmRecoFlv = recoflv
		mmRecoInd = recoInd_
		mmrecopt = recopt
	else:
		mmResp = -1
		mmInd = -1
		mmFlv = -1
		mmRecoFlv = -1
		mmRecoInd = -1
		mmrecopt = -1

	mmout = [mmInd, mmResp, mmFlv, MM, mmRecoInd, mmrecopt, mmRecoFlv]
	return mmout     		

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        ## Getting objects
	jets	  = Collection(event, "Jet")
	if self.isData == False: 
		genjets = Collection(event, "GenJet")
        met       = Object(event, self.metBranchName)

	mmOut = []
      	MM = []
      	ind = []
      	flv = []
      	resp = []
	recoInd_ = []
	if self.isQCD == True:
		mmOut = self.getQCDRespTailCorrector(jets, genjets, met) 
	else:
		mmOut = [-1, -1.0, -1, -1, -1, -1, -1]
	if mmOut[0] >= 0: trueRespGenPT = genjets[mmOut[0]].pt
	else: trueRespGenPT = 0

	self.out.fillBranch("trueResp"     	   , mmOut[1] if (mmOut[0] >= 0) else -1)
	self.out.fillBranch("trueRespFlv"  	   , mmOut[2])
	self.out.fillBranch("trueRespFlvReco"  	   , mmOut[6])
	self.out.fillBranch("trueRespGenPT"	   , trueRespGenPT)
	self.out.fillBranch("trueRespRecoPT"	   , mmOut[5])
	self.out.fillBranch("trueRespMM"	   , mmOut[3])
	self.out.fillBranch("trueRespGenIdx"	   , mmOut[0])
	self.out.fillBranch("trueRespRecoIdx"	   , mmOut[4])

	jetNearMETInd, MMJetDPhi = -1, -1
	for iJ in range(0,len(jets)):
		#if iJ == 4 : break
		if not jets[iJ].Stop0l: continue
		dPhi = abs(deltaPhi(jets[iJ].phi, met.phi))
		if(MMJetDPhi < 0 or dPhi < MMJetDPhi):
			MMJetDPhi = dPhi
			jetNearMETInd = iJ

	if jetNearMETInd < 0 : return True

	pJ = jets[jetNearMETInd]
	passFilter = len(jets) > jetNearMETInd #Could have error here
	pseudoGenPT = self.addFourVec(met, pJ).Pt()
	MMPseudoResp = pJ.pt/pseudoGenPT if pseudoGenPT > 0 else 999

        ### Store output
	self.out.fillBranch("pseudoResp"           , MMPseudoResp)
	self.out.fillBranch("pseudoRespCSV"        , pJ.btagDeepB)
	self.out.fillBranch("pseudoRespPseudoGenPT", pseudoGenPT)
	self.out.fillBranch("pseudoRespPassFilter" , passFilter)
	return True


 # define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
