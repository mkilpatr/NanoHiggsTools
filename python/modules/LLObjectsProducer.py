import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import math
import numpy as np

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest

#2016 MC: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation80XReReco#Data_MC_Scale_Factors_period_dep
#2017 MC: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation94X

class LLObjectsProducer(Module):
    def __init__(self):
        self.metBranchName = "MET"

    def beginJob(self):
        pass
    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
	self.out.branch("Stop0l_MtLepMET", "F",  lenVar="Pass_LeptonVeto")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def SelMtlepMET(self, ele, muon, isks, met):
	mt = []
	for l in ele:
		if l.Stop0l: mt.append(math.sqrt( 2 * met.pt * l.pt * (1 - math.cos(met.phi-l.phi))))
	for l in muon:
		if l.Stop0l: mt.append(math.sqrt( 2 * met.pt * l.pt * (1 - math.cos(met.phi-l.phi))))
	for l in isks:
		if l.Stop0l: mt.append(math.sqrt( 2 * met.pt * l.pt * (1 - math.cos(met.phi-l.phi))))
	return mt

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        ## Getting objects
        electrons = Collection(event, "Electron")
        muons     = Collection(event, "Muon")
        isotracks = Collection(event, "IsoTrack")
	met       = Object(event, self.metBranchName)

        ## Selecting objects
	mt = self.SelMtlepMET(electrons, muons, isotracks, met)
	
        ### Store output
	self.out.fillBranch("Stop0l_MtLepMET",  mt)
	return True


 # define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
