import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import math
import numpy as np

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from PhysicsTools.NanoSUSYTools.modules.datamodelRemap import ObjectRemapped, CollectionRemapped

class Stop0lBaselineProducer(Module):
    def __init__(self, era, isData = False, isFastSim=False, applyUncert=None):
        self.era = era
        self.isFastSim = isFastSim
        self.isData = isData

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

        self.branchMap = {}

    def beginJob(self):
        pass
    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        #self.out.branch("Pass_LeptonTauVeto", "O")
        self.out.branch("Pass_JetID"         + self.suffix, "O")
        self.out.branch("Pass_CaloMETRatio"  + self.suffix, "O", title="ICHEP16 Filter: pfMET/CaloMET < 5")
        self.out.branch("Pass_EventFilter"   + self.suffix, "O")
        self.out.branch("Pass_ElecVeto"      + self.suffix, "O")
        self.out.branch("Pass_MuonVeto"      + self.suffix, "O")
        self.out.branch("Pass_IsoTrkVeto"    + self.suffix, "O")
        self.out.branch("Pass_TauVeto"       + self.suffix, "O")
        self.out.branch("Pass_LeptonVeto"    + self.suffix, "O")
        self.out.branch("Pass_NJets30"       + self.suffix, "O")
        self.out.branch("Pass_MET"           + self.suffix, "O")
        self.out.branch("Pass_HT"            + self.suffix, "O")
        self.out.branch("Pass_Baseline"      + self.suffix, "O")
        # self.out.branch("Pass_HEMVeto20"     + self.suffix, "O", title="HEM Veto 2018: eta[-3, -1.4], phi[-1.57, -0.87], pt > 20")
        # self.out.branch("Pass_HEMVeto30"     + self.suffix, "O", title="HEM Veto 2018: eta[-3, -1.4], phi[-1.57, -0.87], pt > 30")
        self.out.branch("Pass_exHEMVeto20"   + self.suffix, "O", title="HEM Veto 2018: eta[-3.2, -1.2], phi[-1.77, -0.67], pt > 20")
        self.out.branch("Pass_exHEMVeto30"   + self.suffix, "O", title="HEM Veto 2018: eta[-3.2, -1.2], phi[-1.77, -0.67], pt > 30")
        self.out.branch("Jet_sortedIdx"      + self.suffix, "I", lenVar="Jet_nsortedIdx" + self.suffix)

        # Construct Stop0l map
        lob = wrappedOutputTree._branches.keys()
        for bn in lob:
            if self.suffix and "Stop0l" == bn[0:6] and self.suffix in bn:
                self.branchMap[bn[len("Stop0l_"):-len(self.suffix)]] = bn[len("Stop0l_"):]

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def calculateNLeptons(self, eles, muons, isks, taus):
        countEle = sum([e.Stop0l for e in eles])
        countMu  = sum([m.Stop0l for m in muons])
        countIsk = sum([i.Stop0l for i in isks])
	countTau = sum([t.Stop0l for t in taus])
        return countEle, countMu, countIsk, countTau

    def PassEventFilter(self, flags):
        # https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#2016_data
        passEventFilter = None

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 2016 ~~~~~
        if self.era == "2016":
            ## Common filters
            passEventFilter = flags.goodVertices and flags.HBHENoiseFilter and \
                    flags.HBHENoiseIsoFilter and flags.EcalDeadCellTriggerPrimitiveFilter \
                    and flags.BadPFMuonFilter 
            ## Only data
            if self.isData:
                passEventFilter = passEventFilter and flags.globalSuperTightHalo2016Filter and flags.eeBadScFilter
            elif not self.isFastSim:
                passEventFilter = passEventFilter and flags.globalSuperTightHalo2016Filter

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 2017 ~~~~~
        if self.era == "2017" or self.era == "2018":
            ## Common filters
            ## Missing the latest ecalBadCalibReducedMINIAODFilter, not in MiniAOD
            ## But still using the old ecalBadCalibFilter from MiniAOD
            passEventFilter = flags.goodVertices and flags.HBHENoiseFilter and \
                    flags.HBHENoiseIsoFilter and flags.EcalDeadCellTriggerPrimitiveFilter \
                    and flags.BadPFMuonFilter and flags.ecalBadCalibFilter  
            ## Only data
            if self.isData:
                passEventFilter = passEventFilter and flags.globalSuperTightHalo2016Filter and flags.eeBadScFilter
            elif not self.isFastSim:
                passEventFilter = passEventFilter and flags.globalSuperTightHalo2016Filter

        return passEventFilter

    def PassJetID(self, jets):
        # In case of fastsim, it has been observed with lower efficiency
        # https://hypernews.cern.ch/HyperNews/CMS/get/jet-algorithms/379.html
        # The conclusion is to ignore it but cover with systematics.
        if self.isFastSim:
            return True

        # https://twiki.cern.ch/twiki/bin/view/CMS/JetID#Recommendations_for_13_TeV_2017
        # For 2016, loose and tight ID is the same : https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016
        # For 2017, only tight ID available: https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2017
        # Select jet pt > 30GeV, which is used in jet ID study:
        # https://indico.cern.ch/event/592313/contributions/2398387/attachments/1384663/2106630/JetID_JMARmeeting_7_12_2016.pdf
        jetIDs = [j.jetId & 0b010 for j in jets if j.pt > 30]
        return (0 not in jetIDs)


    def PassNjets(self, jets):
        countJets = sum([j.Stop0l for j in jets])
        return countJets >= 2

    def GetJetSortedIdx(self, jets):
        ptlist = []
        etalist = []
        dphiMET = []
        for j in jets:
            if math.fabs(j.eta) > 4.7 or j.pt < 30:
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

    def PassHEMVeto(self, jets, etalow, etahigh, philow, phihigh, ptcut):
        # Calculating HEM veto for 2017 and 2018.
        # Including 2017 in case we need to use 2017 MC for 2018 Data
        if self.era == "2016":
            return True
        for j in jets:
            if (j.eta >= etalow and j.eta <= etahigh) and (j.phi >= philow and j.phi <= phihigh) and j.pt > ptcut:
                return False
        return True

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        if self.applyUncert == "JESUp":
            isotracks = CollectionRemapped(event, "IsoTrack", replaceMap={"Stop0l":"Stop0l_JESUp"})
            jets      = CollectionRemapped(event, "Jet", replaceMap={"pt":"pt_jesTotalUp", "mass":"mass_jesTotalUp", "dPhiMET":"dPhiMET_JESUp", "Stop0l":"Stop0l_JESUp"})
            met       = ObjectRemapped(event,     "MET", replaceMap={"pt":"pt_jesTotalUp", "phi":"phi_jesTotalUp"})
            stop0l    = ObjectRemapped(event,     "Stop0l", replaceMap=self.branchMap)
        elif self.applyUncert == "JESDown":
            isotracks = CollectionRemapped(event, "IsoTrack", replaceMap={"Stop0l":"Stop0l_JESDown"})
            jets      = CollectionRemapped(event, "Jet", replaceMap={"pt":"pt_jesTotalDown", "mass":"mass_jesTotalDown", "dPhiMET":"dPhiMET_JESDown", "Stop0l":"Stop0l_JESDown"})
            met       = ObjectRemapped(event,     "MET", replaceMap={"pt":"pt_jesTotalDown", "phi":"phi_jesTotalDown"})
            stop0l    = ObjectRemapped(event,     "Stop0l", replaceMap=self.branchMap)
        elif self.applyUncert == "METUnClustUp":
            isotracks = CollectionRemapped(event, "IsoTrack", replaceMap={"Stop0l":"Stop0l_METUnClustUp"})
            jets      = CollectionRemapped(event, "Jet", replaceMap={"dPhiMET":"dPhiMET_METUnClustUp"})
            met       = ObjectRemapped(event,     "MET", replaceMap={"pt":"pt_unclustEnUp", "phi":"phi_unclustEnUp"})
            stop0l    = ObjectRemapped(event,     "Stop0l", replaceMap=self.branchMap)
        elif self.applyUncert == "METUnClustDown":
            isotracks = CollectionRemapped(event, "IsoTrack", replaceMap={"Stop0l":"Stop0l_METUnClustDown"})
            jets      = CollectionRemapped(event, "Jet", replaceMap={"dPhiMET":"dPhiMET_METUnClustDown"})
            met       = ObjectRemapped(event,     "MET", replaceMap={"pt":"pt_unclustEnDown", "phi":"phi_unclustEnDown"})
            stop0l    = ObjectRemapped(event,     "Stop0l", replaceMap=self.branchMap)
        else:
            isotracks = Collection(event, "IsoTrack")
            jets      = Collection(event, "Jet")
            met       = Object(event,     "MET")
            stop0l    = Object(event,     "Stop0l")

        caloMET   = Object(event, "CaloMET")
        flags     = Object(event,     "Flag")
        electrons = Collection(event, "Electron")
        muons     = Collection(event, "Muon")
        isotracks = Collection(event, "IsoTrack")
        taus      = Collection(event, "Tau")

        ## Baseline Selection
        PassJetID       = self.PassJetID(jets)
        ## This was an old recommendation in ICHEP16, store this optional bit in case we need it
        ## https://twiki.cern.ch/twiki/bin/viewauth/CMS/SUSRecommendationsICHEP16 
        PassCaloMETRatio= (met.pt / caloMET.pt ) < 5 if caloMET.pt > 0 else True
        PassEventFilter = self.PassEventFilter(flags)

        countEle, countMu, countIsk, countTauPOG = self.calculateNLeptons(electrons, muons, isotracks, taus)
        PassElecVeto   = countEle == 0
        PassMuonVeto   = countMu == 0
        PassIsoTrkVeto = countIsk == 0
        PassTauVeto    = countTauPOG == 0
        PassLeptonVeto  = PassElecVeto and PassMuonVeto and PassIsoTrkVeto and PassTauVeto

        totlep = countEle + countMu
        PassNjets       = self.PassNjets(jets)
        PassMET         = met.pt >= 250
        PassHT          = stop0l.HT >= 300
        ## In case JEC changed jet pt order, resort jets
        sortedIdx, sortedPhi = self.GetJetSortedIdx(jets)

        PassBaseline    = PassEventFilter and PassJetID

        PassHEMVeto20   = self.PassHEMVeto(jets, -3, -1.4, -1.57, -0.87, 20)
        PassHEMVeto30   = self.PassHEMVeto(jets, -3, -1.4, -1.57, -0.87, 30)
        PassexHEMVeto20 = self.PassHEMVeto(jets, -3.2, -1.2, -1.77, -0.67, 20)
        PassexHEMVeto30 = self.PassHEMVeto(jets, -3.2, -1.2, -1.77, -0.67, 30)
        ### Store output
        #self.out.fillBranch("Pass_LeptonTauVeto", PassLeptonTauVeto)
        self.out.fillBranch("Pass_JetID"         + self.suffix, PassJetID)
        self.out.fillBranch("Pass_CaloMETRatio"  + self.suffix, PassCaloMETRatio)
        self.out.fillBranch("Pass_EventFilter"   + self.suffix, PassEventFilter)
        self.out.fillBranch("Pass_ElecVeto"      + self.suffix, PassElecVeto)
        self.out.fillBranch("Pass_MuonVeto"      + self.suffix, PassMuonVeto)
        self.out.fillBranch("Pass_IsoTrkVeto"    + self.suffix, PassIsoTrkVeto)
        self.out.fillBranch("Pass_TauVeto"       + self.suffix, PassTauVeto)
        self.out.fillBranch("Pass_LeptonVeto"    + self.suffix, PassLeptonVeto)
        self.out.fillBranch("Pass_NJets30"       + self.suffix, PassNjets)
        self.out.fillBranch("Pass_MET"           + self.suffix, PassMET)
        self.out.fillBranch("Pass_HT"            + self.suffix, PassHT)
        self.out.fillBranch("Pass_Baseline"      + self.suffix, PassBaseline)
        # self.out.fillBranch("Pass_HEMVeto20"     + self.suffix, PassHEMVeto20)
        # self.out.fillBranch("Pass_HEMVeto30"     + self.suffix, PassHEMVeto30)
        self.out.fillBranch("Pass_exHEMVeto20"   + self.suffix, PassexHEMVeto20)
        self.out.fillBranch("Pass_exHEMVeto30"   + self.suffix, PassexHEMVeto30)
        self.out.fillBranch("Jet_nsortedIdx"     + self.suffix, len(sortedIdx))
        self.out.fillBranch("Jet_sortedIdx"      + self.suffix, sortedIdx)
        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
# Stop0lBaseline = lambda : Stop0lBaselineProducer("2016", False)
