import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import math
import numpy as np

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaR

from PhysicsTools.NanoHiggsTools.modules.datamodelRemap import ObjectRemapped, CollectionRemapped

#2016 MC: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation2016Legacy
#2017 MC: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation94X
#2018 MC: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation102X
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

class Stop0lObjectsProducer(Module):
    def __init__(self, era, applyUncert = None):
        self.era = era
        self.metBranchName = "MET"

        self.applyUncert = applyUncert

        self.suffix = ""

        if self.applyUncert == "JESUp":
            self.suffix = "_JESUp"
        elif self.applyUncert == "METUnClustUp":
            self.suffix = "_METUnClustUp"
        elif self.applyUncert == "JESDown":
            self.suffix = "_JESDown"
        elif self.applyUncert == "METUnClustDown":
            self.suffix = "_METUnClustDown"

    def beginJob(self):
        pass
    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        if self.applyUncert == None:
            self.out.branch("Electron_Stop0l" + self.suffix, "O", lenVar="nElectron", title="cutBased Veto ID with miniISO < 0.1, pT > 5")
            self.out.branch("Muon_Stop0l"     + self.suffix, "O", lenVar="nMuon")
            self.out.branch("Tau_Stop0l"      + self.suffix, "O", lenVar="nTau")
            self.out.branch("Stop0l_nSoftb"   + self.suffix, "I")
        
        if self.applyUncert == None or "JES" in self.applyUncert or "METUnClust" in self.applyUncert:
            self.out.branch("Electron_MtW"    + self.suffix, "F", lenVar="nElectron")
            self.out.branch("Muon_MtW"        + self.suffix, "F", lenVar="nMuon")
            self.out.branch("IsoTrack_MtW"    + self.suffix, "F", lenVar="nIsoTrack")
            self.out.branch("Tau_MtW"	      + self.suffix, "F", lenVar="nTau")
            self.out.branch("Jet_dPhiMET"     + self.suffix, "F", lenVar="nJet")
            self.out.branch("IsoTrack_Stop0l" + self.suffix, "O", lenVar="nIsoTrack")
            self.out.branch("Stop0l_Mtb"      + self.suffix, "F")
            self.out.branch("Stop0l_METSig"   + self.suffix, "F")

        if self.applyUncert == None or "JES" in self.applyUncert:
            self.out.branch("Jet_Stop0l" + self.suffix,      "O", lenVar="nJet")
            self.out.branch("Jet_btagStop0l" + self.suffix,  "O", lenVar="nJet")
            self.out.branch("Stop0l_HT" + self.suffix,       "F")
            self.out.branch("Stop0l_Ptb" + self.suffix,      "F")
            self.out.branch("Stop0l_nJets" + self.suffix,    "I")
            self.out.branch("Stop0l_nbtags" + self.suffix,   "I")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass


    def SelEle(self, ele):
        if math.fabs(ele.eta) > 2.5 or ele.pt < 13:
            return False
        ## Veto ID electron
        if ele.cutBasedNoIso < 1:
            return False
        ## MiniIso < 0.1
        if ele.miniPFRelIso_all > 0.1:
            return False
        return True

    def SelMuon(self, mu):
        ## NanoAOD store loose ID Muon by default
        if math.fabs(mu.eta) > 2.4 or mu.pt < 10:
            return False
        ## MiniIso < 0.1
        if mu.miniPFRelIso_all > 0.2:
            return False
        return True

    def SelIsotrack(self, isk, met):
        iso = isk.pfRelIso03_chg
        if abs(isk.pdgId) == 11 or abs(isk.pdgId) == 13:
            if isk.pt < 5 or iso > 0.2:
                return False
        if abs(isk.pdgId) == 211:
            if isk.pt < 10 or iso > 0.1:
                return False
        mtW = self.CalMtW(isk, met)
        if mtW  > 100:
            return False
        return True

    def SelTauPOG(self, tau, met):
        if tau.pt < 30 or abs(tau.eta) > 2.3 or not tau.idDecayMode or not (tau.idMVAoldDM2017v2 & 8):
                return False
        return True

    def CalMtW(self, lep, met):
        return math.sqrt( 2 * met.pt * lep.pt * (1 - math.cos(ROOT.TVector2.Phi_mpi_pi(met.phi-lep.phi))))

    def SelBtagJets(self, jet):
        global DeepCSVMediumWP
        if jet.pt < 20 or math.fabs(jet.eta) > 2.4:
            return False
        if jet.btagDeepB >= DeepCSVMediumWP[self.era]:
            return True
        return False

    def SelSoftb(self, isv, jets):
        ## Select soft bs
        ## SV is not associate with selected jets
        for j in jets:
            if j.pt >= 20 and math.fabs(j.eta) <= 2.4 and deltaR(j.eta, j.phi, isv.eta, isv.phi) <= 0.4 :
                return False
        if isv.ntracks < 3 or math.fabs(isv.dxy)>3. or isv.dlenSig <4:
            return False
        if isv.DdotP < 0.98 :
            return False
        return True


    def SelJets(self, jet):
        if jet.pt < 30 or math.fabs(jet.eta) > 2.4:
            return False
        return True

    def CalHT(self, jets):
        HT = sum([j.pt for i, j in enumerate(jets) if self.Jet_Stop0l[i]])
        return HT

    def CalMTbPTb(self, jets, met):
        Mtb = float('inf')  #Max value
        Ptb = 0

        # Getting bjet, ordered by pt
        bjets = [ j for i,j in enumerate(jets) if self.BJet_Stop0l[i]]
        # Getting btag index, ordered by b discriminator value
        btagidx = sorted(range(len(bjets)), key=lambda k: bjets[k].btagDeepB , reverse=True)

        for i in range(min(len(bjets), 2)):
            Ptb += bjets[i].pt

        for i in range(min(len(btagidx), 2)):
            bj = bjets[btagidx[i]]
            Mtb = min(Mtb, math.sqrt( 2 * met.pt * bj.pt * (1 - math.cos(ROOT.TVector2.Phi_mpi_pi(met.phi-bj.phi)))))

        if Mtb == float('inf'):
            Mtb = 0
        return Mtb, Ptb

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        ## Getting objects
        electrons = Collection(event, "Electron")
        muons     = Collection(event, "Muon")
        isotracks = Collection(event, "IsoTrack")
        taus	  = Collection(event, "Tau")

        if self.applyUncert == "JESUp":
            jets      = CollectionRemapped(event, "Jet", replaceMap={"pt":"pt_jesTotalUp", "mass":"mass_jesTotalUp"})
            met       = ObjectRemapped(event, self.metBranchName, replaceMap={"pt":"pt_jesTotalUp", "phi":"phi_jesTotalUp"})
        elif self.applyUncert == "JESDown":
            jets      = CollectionRemapped(event, "Jet", replaceMap={"pt":"pt_jesTotalDown", "mass":"mass_jesTotalDown"})
            met       = ObjectRemapped(event, self.metBranchName, replaceMap={"pt":"pt_jesTotalDown", "phi":"phi_jesTotalDown"})
        elif self.applyUncert == "METUnClustUp":
            jets      = Collection(event, "Jet")
            met       = ObjectRemapped(event, self.metBranchName, replaceMap={"pt":"pt_unclustEnUp", "phi":"phi_unclustEnUp"})
        elif self.applyUncert == "METUnClustDown":
            jets      = Collection(event, "Jet")
            met       = ObjectRemapped(event, self.metBranchName, replaceMap={"pt":"pt_unclustEnDown", "phi":"phi_unclustEnDown"})
        else:
            jets      = Collection(event, "Jet")
            met       = Object(event, self.metBranchName)

        flags   = Object(event,     "Flag")
       
        ## Selecting objects
        self.Electron_Stop0l = map(self.SelEle, electrons)
        self.Muon_Stop0l     = map(self.SelMuon, muons)
        self.Tau_Stop0l      = map(lambda x : self.SelTauPOG(x, met), taus)
        self.Electron_MtW    = map(lambda x : self.CalMtW(x, met), electrons)
        self.Muon_MtW        = map(lambda x : self.CalMtW(x, met), muons)
        self.IsoTrack_MtW    = map(lambda x : self.CalMtW(x, met), isotracks)
        self.IsoTrack_Stop0l = map(lambda x : self.SelIsotrack(x, met), isotracks)
        self.Tau_MtW	     = map(lambda x : self.CalMtW(x, met), taus)
        self.Jet_Stop0l      = map(self.SelJets, jets)
        self.BJet_Stop0l     = map(self.SelBtagJets, jets)

        ## Jet variables
        jet_phi = np.asarray([jet.phi for jet in jets])
        Jet_dPhi = jet_phi - met.phi
        np.subtract(Jet_dPhi, 2*math.pi, out = Jet_dPhi, where= (Jet_dPhi >=math.pi))
        np.add(Jet_dPhi, 2*math.pi,  out =Jet_dPhi , where=(Jet_dPhi < -1*math.pi))
        np.fabs(Jet_dPhi, out=Jet_dPhi)
        ## TODO: Need to improve speed
        HT = self.CalHT(jets)
        Mtb, Ptb = self.CalMTbPTb(jets, met)

        ### Store output
        if self.applyUncert == None:
            self.out.fillBranch("Electron_Stop0l" + self.suffix, self.Electron_Stop0l)
            self.out.fillBranch("Muon_Stop0l" + self.suffix,     self.Muon_Stop0l)
            self.out.fillBranch("Tau_Stop0l" + self.suffix,	 self.Tau_Stop0l)

        if self.applyUncert == None or "JES" in self.applyUncert or "METUnClust" in self.applyUncert:
            self.out.fillBranch("IsoTrack_Stop0l" + self.suffix, self.IsoTrack_Stop0l)
            self.out.fillBranch("Electron_MtW" + self.suffix   , self.Electron_MtW)
            self.out.fillBranch("Muon_MtW" + self.suffix       , self.Muon_MtW)
            self.out.fillBranch("IsoTrack_MtW" + self.suffix   , self.IsoTrack_MtW)
            self.out.fillBranch("Tau_MtW" + self.suffix	       , self.Tau_MtW)
            self.out.fillBranch("Jet_dPhiMET" + self.suffix    , Jet_dPhi)
            self.out.fillBranch("Stop0l_Mtb" + self.suffix     , Mtb)
            self.out.fillBranch("Stop0l_METSig" + self.suffix  , met.pt / math.sqrt(HT) if HT > 0 else 0)

        if self.applyUncert == None or "JES" in self.applyUncert:
            self.out.fillBranch("Jet_btagStop0l" + self.suffix,  self.BJet_Stop0l)
            self.out.fillBranch("Jet_Stop0l" + self.suffix,      self.Jet_Stop0l)
            self.out.fillBranch("Stop0l_HT" + self.suffix,       HT)
            self.out.fillBranch("Stop0l_Ptb" + self.suffix,      Ptb)
            self.out.fillBranch("Stop0l_nJets" + self.suffix,    sum(self.Jet_Stop0l))
            self.out.fillBranch("Stop0l_nbtags" + self.suffix,   sum(self.BJet_Stop0l))

        return True


 # define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
