import os, sys
from os import path
import io
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import math
import numba
import numpy as np
import gzip
import json
from array import array
from importlib import import_module

from PhysicsTools.NanoSUSYTools.modules.HelperFunctions import *
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest
from PhysicsTools.NanoAODTools.postprocessing.framework.treeReaderArrayTools import *
from rootpy.tree import Tree, TreeModel, IntCol, FloatArrayCol

@numba.jit(nopython=True)
def recursiveMotherSearch(startIdx, targetIdx, GenPartCut_genPartIdxMother):
    if startIdx < 0:
        return False

    mom = GenPartCut_genPartIdxMother[startIdx]

    if mom < 0:
        return False
    elif startIdx == targetIdx:
        return True
    else:
        return recursiveMotherSearch(mom, targetIdx, GenPartCut_genPartIdxMother)

@numba.jit(nopython=True)
def genParticleAssociation(GenPart_genPartIdxMother, GenPart_pdgId, GenPart_statusFlags):
    GenPart_momPdgId = GenPart_pdgId[GenPart_genPartIdxMother]
    #set any particles with an index of -1 to id 0                                                                                                                      
    GenPart_momPdgId[GenPart_genPartIdxMother < 0] = 0

    genTauDaughters = []
    for iGP, pdgId in enumerate(GenPart_pdgId):
        if (GenPart_statusFlags[iGP] & 0x2100) != 0x2100:
            continue
        if abs(pdgId) == 15:
            gtd = []
            for iGP2, pdgId2 in enumerate(GenPart_pdgId):
                if (abs(pdgId2) == 11 or abs(pdgId2) == 13 or abs(pdgId2) == 211) and (GenPart_statusFlags[iGP2] & 0x2008) == 0x2008:
                    if recursiveMotherSearch(iGP2, iGP, GenPart_genPartIdxMother):
                        gtd.append(iGP2)
            genTauDaughters.extend(gtd)

    return genTauDaughters

def deltaRMatch(HiggsCandPdgId, HiggsCandEta, HiggsCandPhi, genTauDaughters_pdgId, genTauDaughters_eta, genTauDaughters_phi):

    matches = np.zeros(len(HiggsCandEta), dtype=int)
    matchesPdgId = np.zeros(len(HiggsCandEta), dtype=int)

    if len(genTauDaughters_eta):

        higgsEtaVals = np.array(np.meshgrid(HiggsCandEta, genTauDaughters_eta)).T.reshape(-1,2)
        higgsPhiVals = np.array(np.meshgrid(HiggsCandPhi, genTauDaughters_phi)).T.reshape(-1,2)
        higgsPdgIdVals = np.array(np.meshgrid(HiggsCandPdgId, genTauDaughters_pdgId)).T.reshape(-1,2)

        ## Using ufunc for vector operation
        deta = np.power(higgsEtaVals[:,0] - higgsEtaVals[:,1], 2)
        dPhi = higgsPhiVals[:,0] - higgsPhiVals[:,1]
        dR = np.sqrt((( abs(abs(dPhi)-np.pi)-np.pi )**2+(deta)**2)).reshape([-1,len(genTauDaughters_eta)/len(genTauDaughters_eta), len(genTauDaughters_eta)])
        ## Check PDG agreement
        matchPdgID = (higgsPdgIdVals[:,0] == abs(higgsPdgIdVals[:,1])).reshape([-1,len(genTauDaughters_eta)/len(genTauDaughters_eta), len(genTauDaughters_eta)])
        
        matches[dR.max(axis=2).min(axis=1) < 0.8] = 1
        matchesPdgId[matchPdgID.max(axis=2).min(axis=1)] = 1

    return matches, matchesPdgId

TauDecay = {
    -999. : 'kUnknown',
    -111. : 'kUnknown',
    1 : 'kTauToHadDecay',
    2 : 'kTauToElecDecay',
    3 : 'kTauToMuDecay',
}

tauDecayPdgId = {
    'kUnknown' : -1,
    'kTauToHadDecay' : 211,
    'kTauToElecDecay' : 11,
    'kTauToMuDecay' : 13,
}

class HiggsJSONProducer(Module):
    def __init__(self, processName):
        self.metBranchName = "MET"
        self.filename=processName + '.json.gz'

    def beginJob(self):
        pass
    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.fout = gzip.open(self.filename, 'w')

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.fout.close()

    def makeTLorentzVector(self, svfit, idx, isTau=""):
        v = ROOT.TLorentzVector()
        if isTau=="tau1":   v.SetPtEtaPhiM(svfit[idx].tau1Pt, svfit[idx].tau1Eta, svfit[idx].tau1Phi, svfit[idx].tau1Mass)
        elif isTau=="tau2": v.SetPtEtaPhiM(svfit[idx].tau2Pt, svfit[idx].tau2Eta, svfit[idx].tau2Phi, svfit[idx].tau2Mass)
        else:               v.SetPtEtaPhiM(svfit[idx].Pt, svfit[idx].Eta, svfit[idx].Phi, svfit[idx].Mass)
        return v

    def higg2json(self, svfit, idx):
        h = self.makeTLorentzVector(svfit, idx)
        t1 = self.makeTLorentzVector(svfit, idx, 'tau1')
        t2 = self.makeTLorentzVector(svfit, idx, 'tau2')

        
        j = [{'E':float(h.E()), 'px':float(h.Px()), 'py':float(h.Py()), 'pz':float(h.Pz())}]
        j.append({'E':float(t1.E()), 'px':float(t1.Px()), 'py':float(t1.Py()), 'pz':float(t1.Pz())})
        j.append({'E':float(t2.E()), 'px':float(t2.Px()), 'py':float(t2.Py()), 'pz':float(t2.Pz())})
        return j

    class TTreeReaderArrayWrapper:
        def __init__(self, ttarray):
            self.ttarray = ttarray

        def __iter__(self):
            for i in xrange(len(self.ttarray)):
                yield self.ttarray[i]
            return

    def HiggsGenMatch(self, event, HiggsCandEta, HiggsCandPhi, HiggsCandPdgId):
        
        GenPart_eta              = np.fromiter(self.TTreeReaderArrayWrapper(event.GenPart_eta),              dtype=float)
        GenPart_phi              = np.fromiter(self.TTreeReaderArrayWrapper(event.GenPart_phi),              dtype=float)

        GenPart_genPartIdxMother = np.fromiter(self.TTreeReaderArrayWrapper(event.GenPart_genPartIdxMother), dtype=int)
        GenPart_pdgId            = np.fromiter(self.TTreeReaderArrayWrapper(event.GenPart_pdgId),            dtype=int)
        GenPart_statusFlags      = np.fromiter(self.TTreeReaderArrayWrapper(event.GenPart_statusFlags),      dtype=int)
    
        genTauDaughters_list = genParticleAssociation(GenPart_genPartIdxMother, GenPart_pdgId, GenPart_statusFlags)
    
        genTauDaughters = np.array(genTauDaughters_list)

        if(len(genTauDaughters)):
            genTauDaughters_pdgid = GenPart_pdgId[genTauDaughters]
            genTauDaughters_eta = GenPart_eta[genTauDaughters]
            genTauDaughters_phi = GenPart_phi[genTauDaughters]
        else:
            genTauDaughters_pdgid = np.array([])
            genTauDaughters_eta = np.array([])
            genTauDaughters_phi = np.array([])
   
        matches, matchPdgID = deltaRMatch(HiggsCandPdgId, HiggsCandEta, HiggsCandPhi, genTauDaughters_pdgid, genTauDaughters_eta, genTauDaughters_phi)

        return matches, matchPdgID, genTauDaughters 

    def recursiveFindHiggs(self, startIdx, targetPdgId, genpart):
        if startIdx < 0:
            return -1
    
        mom = genpart[startIdx].genPartIdxMother
    
        if mom < 0:
            return -1
        elif genpart[startIdx].pdgId == targetPdgId:
            return mom
        else:
            return self.recursiveFindHiggs(mom, targetPdgId, genpart)

    def analyze(self, event):
        ## Getting objects
        svfit     = Collection(event, "SVFit")
        svfitmet  = Collection(event, "SVFitMET")
        genpart   = Collection(event, "GenPart")

        for iG in xrange(len(genpart)):
            g = genpart[iG]
            if isA(g.pdgId, 25):# and hasBit(g.statusFlags, 0) and hasBit(g.statusFlags, 13):
                print("Is a Higgs and index: {0}".format(iG))

        tau1Pt = np.fromiter(self.TTreeReaderArrayWrapper(event.SVFit_tau1Pt), dtype=float)
        tau1Eta = np.fromiter(self.TTreeReaderArrayWrapper(event.SVFit_tau1Eta), dtype=float)
        tau1Phi = np.fromiter(self.TTreeReaderArrayWrapper(event.SVFit_tau1Phi), dtype=float)
        tau1pdgId = np.fromiter(self.TTreeReaderArrayWrapper(event.SVFit_tau1pdgId), dtype=float)
        tau1pdgIdTranslate = [tauDecayPdgId[TauDecay[t]] for t in tau1pdgId]
        tau1GenMatch, tau1GenMatchPdgId, whichTau = self.HiggsGenMatch(event, tau1Eta, tau1Phi, tau1pdgIdTranslate)
        print("tau1: {0} {1} {2} {3}".format(tau1GenMatch, whichTau, tau1pdgIdTranslate, tau1GenMatchPdgId))
        for gD in whichTau:
            higgsIdx = self.recursiveFindHiggs(gD, 25, genpart)
            print("Matched Higgs Index: {0}".format(higgsIdx))
        for t1 in xrange(len(tau1GenMatch)):
            if tau1GenMatch[t1]:
                print("Decay Channel: {0}, Higgs mass: {1}".format(svfit[t1].channel, svfit[t1].Mass))
                print("matched pdgID: {0}".format(genpart[t1].pdgId))

        for idx in xrange(len(svfit)):
            if svfit[svfit[idx].Index].Pt > 0.:
                j = self.higg2json(svfit, svfit[idx].Index)
                self.fout.write((json.dumps(j, sort_keys=False)+'\n').encode('utf-8'))

        return True
        
        
 # define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
