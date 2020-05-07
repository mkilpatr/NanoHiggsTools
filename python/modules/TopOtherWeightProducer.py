import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import math
import numpy as np
from functools import reduce
import operator
import os

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoSUSYTools.modules.datamodelRemap import ObjectRemapped, CollectionRemapped
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest

#2016 MC: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation80XReReco#Data_MC_Scale_Factors_period_dep
#2017 MC: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation94X

DeepCSVMediumWP ={
    "2016" : 0.6321,
    "2017" : 0.4941,
    "2018" : 0.4184
}

DeepCSVLooseWP ={
    "2016" : 0.2217,
    "2017" : 0.1522,
    "2018" : 0.1241
}

CSVv2MediumWP = {
    "2016" : 0.8484,
    "2017" : 0.8838,
    "2018" : 0.8838  # Not recommended, use 2017 as temp
}

class TopOtherWeightProducer(Module):
    def __init__(self, era, Process, isData = False):
        self.era = era
        self.sampleName = Process
        self.isData = isData
        self.metBranchName = "MET"
	self.xsRoot = os.environ["CMSSW_BASE"] + "/src/PhysicsTools/NanoSUSYTools/data/toppt/LostLepton_HM_topAK8_weight.root"
	self.xsRoot_ptlepbmet = os.environ["CMSSW_BASE"] + "/src/PhysicsTools/NanoSUSYTools/data/toppt/LostLepton_HM_ptlepbmet_weight.root"
        self.hist_AK8_2016 = "FatJet_pt[0]_singlelep-2016__llcr_hm_2016__num__"
        self.hist_AK8_2017 = "FatJet_pt[0]_singlelep-2017__llcr_hm_2017__num__"
        self.hist_AK8_2018 = "FatJet_pt[0]_singlelep-2018__llcr_hm_2018__num__"
        self.hist_Ptb_2016 = "Stop0l_PtLepMetB_singlelep-2016__llcr_hm_2016__over__bkgtotal"
        self.hist_Ptb_2017 = "Stop0l_PtLepMetB_singlelep-2017__llcr_hm_2017__over__bkgtotal"
        self.hist_Ptb_2018 = "Stop0l_PtLepMetB_singlelep-2018__llcr_hm_2018__over__bkgtotal"

        ## Gen pdgid for toppt weight
        self.p_tauminus = 15
        self.p_Z0       = 23
        self.p_Wplus    = 24
        self.p_gamma    = 22
        self.pfhfhad = 1
        self.pfem = 2
        self.pfelectron = 11
        self.p_nu_e = 12
        self.pfmuon = 13
        self.p_nu_mu = 14
        self.p_nu_tau = 16
        self.pfphoton = 22
        self.pfh0 = 130
        self.pfhplus = 211
        
    def beginJob(self):
        pass
    def endJob(self):
        pass

    def loadhisto(self,filename,hname):
        file =ROOT.TFile.Open(filename)
        hist_ = file.Get(hname)
        hist_.SetDirectory(0)
        return hist_

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        if "TTbar" in self.sampleName:
	    if self.era == '2016':   
                self.topAK8LLHMcr=self.loadhisto(self.xsRoot,self.hist_AK8_2016)
                self.topPtLepBMetLLHMcr=self.loadhisto(self.xsRoot_ptlepbmet,self.hist_Ptb_2016)
	    elif self.era == '2017': 
                self.topAK8LLHMcr=self.loadhisto(self.xsRoot,self.hist_AK8_2017)
                self.topPtLepBMetLLHMcr=self.loadhisto(self.xsRoot_ptlepbmet,self.hist_Ptb_2017)
	    elif self.era == '2018': 
                self.topAK8LLHMcr=self.loadhisto(self.xsRoot,self.hist_AK8_2018)
                self.topPtLepBMetLLHMcr=self.loadhisto(self.xsRoot_ptlepbmet,self.hist_Ptb_2018)
        self.out.branch("Stop0l_nMatchTopPt"                         ,  "I")
        self.out.branch("Stop0l_MatchTopPt"                          ,  "F", lenVar="Stop0l_nMatchTopPt") 
        self.out.branch("Stop0l_nMatchWPt"                           ,  "I")
        self.out.branch("Stop0l_MatchWPt"                            ,  "F", lenVar="Stop0l_nMatchTopPt") 
        self.out.branch("Stop0l_PtLepMetB"                           ,  "F")
        if not self.isData:
            self.out.branch('Stop0l_ntopPTLep', 			"I")
            self.out.branch('Stop0l_ntopPTHad', 			"I")
            self.out.branch('Stop0l_topPTLep', 			        "F", lenVar="Stop0l_ntopPTLep")
            self.out.branch('Stop0l_topPTHad', 			        "F", lenVar="Stop0l_ntopPTHad")
            self.out.branch('Stop0l_topAK8Weight', 			"F")
            self.out.branch('Stop0l_topPtLepBMetWeight', 		"F")

    def isA(self, particleID, p):
        return abs(p) == particleID

    def topPTLep(self, genpars):
        toppt_leptonic = []
        toppt_hadronic = []
        
        #for genpar in genpars :     #loop on genpars       
        #print len(genpars)
        for i in range(len(genpars)):
            genpar = genpars[i]
            if genpar.statusFlags & 8192 == 0: continue
            ## look at top daughters
            if abs(genpar.pdgId)==6:
                #print genpar.pt,genpar.pdgId
                for j in range(len(genpars)):
                    genpard = genpars[j]
                    if genpard.genPartIdxMother == i:
                        #print genpard.pt,genpard.pdgId
        
                        ## look at W daughters
                        if self.isA(self.p_Wplus, genpard.pdgId):
                            for k in range(len(genpars)):
                                genpardd = genpars[k]
                                if genpardd.genPartIdxMother == j:
                                    #print genpardd.pt,genpardd.pdgId
                                    if abs(genpardd.pdgId)>=11 and abs(genpardd.pdgId)<=16:
                                        toppt_leptonic.append(genpar.pt)
                                        break
                                    if abs(genpardd.pdgId)>=1 and abs(genpardd.pdgId)<=4:
                                        toppt_hadronic.append(genpar.pt)
                                        break

        toppt_leptonic.sort(reverse=True)
        toppt_hadronic.sort(reverse=True)
        #print("Leptonic: {0}".format(toppt_leptonic))
        #print("Hadronic: {0}".format(toppt_hadronic))

        return toppt_leptonic, toppt_hadronic 

    def topPTAK8Match(self, fatjet):
        topptWeight = 1.
        if len(fatjet) == 0: return topptWeight

        if fatjet[0].pt >= 150:
            topptWeight = self.topAK8LLHMcr.GetBinContent(self.topAK8LLHMcr.GetNbinsX()) if fatjet[0].pt >= 800 else self.topAK8LLHMcr.GetBinContent(self.topAK8LLHMcr.GetXaxis().FindBin(fatjet[0].pt))

        return topptWeight

    def SelBtagJets(self, jet):
        global DeepCSVMediumWP
        if jet.pt < 20 or math.fabs(jet.eta) > 2.4:
            return False
        if jet.btagDeepB >= DeepCSVMediumWP[self.era]:
            return True
        return False

    def ptlepmetb(self, electron, muon, tau, met, jets):
        topptWeight = 1.
        bjets = [ j for i,j in enumerate(jets) if self.BJet_Stop0l[i]]
        btagidx = sorted(range(len(bjets)), key=lambda k: bjets[k].btagDeepB , reverse=True)

        lep = electron if (len(electron) != 0 and len(muon) == 0) else muon
        if len(lep) == 0: lep_4vec = ROOT.TLorentzVector(0, 0, 0, 0)
        else:             lep_4vec = ROOT.TLorentzVector(lep[0].pt, lep[0].eta, lep[0].phi, lep[0].mass)

        Mtb = float('inf')

        met_4vec = ROOT.TLorentzVector()
        met_4vec.SetPtEtaPhiM(met.pt, 0, met.phi, 0)
        tot = ROOT.TLorentzVector()
        for b in range(min(len(btagidx), 2)):
            bj = bjets[btagidx[b]]
            bj_4vec = ROOT.TLorentzVector(bj.pt, bj.eta, bj.phi, bj.mass)
            if Mtb <= math.sqrt(2 * met.pt * bj.pt * (1 - math.cos(ROOT.TVector2.Phi_mpi_pi(met.phi-bj.phi)))): continue
            tot = bj_4vec + lep_4vec

        ptlepbmet = (met_4vec + tot)

        return ptlepbmet

    def topPTptlepbmetMatch(self, ptlepbmet):
        topptWeight = 1.
        if ptlepbmet == 0: return topptWeight

        topptWeight = self.topPtLepBMetLLHMcr.GetBinContent(self.topPtLepBMetLLHMcr.GetNbinsX()) if ptlepbmet >= 800 else self.topPtLepBMetLLHMcr.GetBinContent(self.topPtLepBMetLLHMcr.GetXaxis().FindBin(ptlepbmet))

        return topptWeight

    def FatJetMatchedPt(self, obj, objType):
	pt = []
	for s in obj:
		if not (s.Stop0l == objType): continue
		pt.append(s.pt)

	return pt

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        ## Getting objects
        if not self.isData: genpart   = Collection(event, "GenPart")
        electrons = Collection(event, "Electron")
        muons     = Collection(event, "Muon")
        taus      = Collection(event, "Tau")
        stop0l    = Object(event,     "Stop0l")
        fatjets   = Collection(event, "FatJet")
        jets      = Collection(event, "Jet")
        met       = Object(event, self.metBranchName)
        
        ## Selecting objects
        ## Adding matched Fatjet pT
        toppt                = self.FatJetMatchedPt(fatjets, 1)
        wpt                  = self.FatJetMatchedPt(fatjets, 2)
        self.BJet_Stop0l     = map(self.SelBtagJets, jets)
        ptlepmetb            = self.ptlepmetb(electrons, muons, taus, met, jets)
        
        self.out.fillBranch("Stop0l_nMatchTopPt",                          len(toppt))
        self.out.fillBranch("Stop0l_MatchTopPt",                           toppt)
        self.out.fillBranch("Stop0l_nMatchWPt",                            len(wpt))
        self.out.fillBranch("Stop0l_MatchWPt",                             wpt)
        self.out.fillBranch("Stop0l_PtLepMetB",                            ptlepmetb.Pt())
 
        if not self.isData:
            toppt_lepmatch, toppt_hadmatch = self.topPTLep(genpart)
            if "TTbar" in self.sampleName:
                topptAK8  = self.topPTAK8Match(fatjets) 
                topptPtlepbmet  = self.topPTptlepbmetMatch(ptlepmetb.Pt()) 
            else:
                topptAK8 = 1.
                topptPtlepbmet = 1.
                toppt_wgt = 1.
                toppt_only = 1.
                toppt_mgpow = 1.
                toppt_up = 1.
                toppt_dn = 1.
            #print("toppt_wgt: {0}".format(toppt_lep))
            self.out.fillBranch('Stop0l_ntopPTLep', 			len(toppt_lepmatch))
            self.out.fillBranch('Stop0l_ntopPTHad', 			len(toppt_hadmatch))
            self.out.fillBranch('Stop0l_topPTLep', 			toppt_lepmatch)
            self.out.fillBranch('Stop0l_topPTHad', 			toppt_hadmatch)
            self.out.fillBranch('Stop0l_topAK8Weight', 			topptAK8)
            self.out.fillBranch('Stop0l_topPtLepBMetWeight', 		topptPtlepbmet)
        return True


 # define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
