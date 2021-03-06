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

from PhysicsTools.NanoHiggsTools.modules.HelperFunctions import *
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

def deltaRMatch(CandPdgId, CandEta, CandPhi, genTauDaughters_pdgId, genTauDaughters_eta, genTauDaughters_phi):

    matches = np.zeros(len(CandEta), dtype=int)
    candmatches = np.zeros(len(genTauDaughters_eta), dtype=int)
    matchesPdgId = np.zeros(len(CandEta), dtype=int)

    if len(genTauDaughters_eta):

        EtaVals = np.array(np.meshgrid(CandEta, genTauDaughters_eta)).T.reshape(-1,2)
        PhiVals = np.array(np.meshgrid(CandPhi, genTauDaughters_phi)).T.reshape(-1,2)
        PdgIdVals = np.array(np.meshgrid(CandPdgId, genTauDaughters_pdgId)).T.reshape(-1,2)

        #print("candEta: {0}, genTauEta: {1}".format(CandEta, genTauDaughters_eta))

        ## Using ufunc for vector operation
        deta = np.power(EtaVals[:,0] - EtaVals[:,1], 2)
        dPhi = PhiVals[:,0] - PhiVals[:,1]
        dR = np.sqrt((( abs(abs(dPhi)-np.pi)-np.pi )**2+(deta)**2)).reshape([-1, len(genTauDaughters_eta), 1])
        dRMatch = dR.T.reshape([-1, len(CandEta), 1])
        ## Check PDG agreement
        matchPdgID = (PdgIdVals[:,0] == abs(PdgIdVals[:,1])).reshape([-1, len(genTauDaughters_eta), 1])

        #print("dR: {0}, dRMatches: {1}".format(dR, dRMatch))

        matches[dR.max(axis=2).min(axis=1) < 0.8] = 1
        candmatches[dRMatch.max(axis=2).min(axis=1) < 0.8] = 1
        matchesPdgId[matchPdgID.max(axis=2).min(axis=1)] = 1
        #print("matches: {0}, candmatches: {1}".format(matches, candmatches))

    return matches, candmatches, matchesPdgId

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
    def __init__(self, processName, match, debug):
        self.metBranchName = "MET"
        self.processName = processName
        self.match = match
        self.debug = debug

    def beginJob(self):
        pass
    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.fout = gzip.open('genTau.json.gz', 'a')
        self.fgenout = gzip.open("genHiggs.json.gz", 'a')

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.fout.close()
        self.fgenout.close()

    def makeTLorentzVector(self, svfit, isTau=""):
        v = ROOT.TLorentzVector()
        if isTau=="tau1":   v.SetPtEtaPhiM(svfit.tau1Pt, svfit.tau1Eta, svfit.tau1Phi, svfit.tau1Mass)
        elif isTau=="tau2": v.SetPtEtaPhiM(svfit.tau2Pt, svfit.tau2Eta, svfit.tau2Phi, svfit.tau2Mass)
        elif isTau=="gen":  v.SetPtEtaPhiM(svfit.pt, svfit.eta, svfit.phi, svfit.mass)
        else:               v.SetPtEtaPhiM(svfit.Pt, svfit.Eta, svfit.Phi, svfit.Mass)
        return v

    def checkJSON(self, jin, jcheck):
        for j1, j in enumerate(jin):
            #if self.debug: print("{0} and {1}".format(j, jcheck[0]))
            if j['E'] == jcheck[0]['E'] and j['px'] == jcheck[0]['px'] and j['py'] == jcheck[0]['py'] and j['pz'] == jcheck[0]['pz']:
                return False
        return True

    def checkJSONNaN(self, jin):
        for j1, j in enumerate(jin):
            if math.isnan(j['E']) or math.isnan(j['px']) or math.isnan(j['py']) or math.isnan(j['pz']):
                return True
        return False

    def tau2json(self, svfit, tauGenMatch, tauGenPartMatch, svfitmet, tau = None, type = None):
        j = []
        isFirst = True
        for idx, sv in enumerate(svfit):
            #print("Pass Baseline: {0} Pass Lepton: {1}".format(sv.PassBaseline, sv.PassLepton))
            if not sv.PassBaseline or not sv.PassLepton: continue
            if (tauGenMatch[idx]):
                if isFirst: h = self.makeTLorentzVector(tau, "gen")
                t = self.makeTLorentzVector(sv, type)
                
                if isFirst:
                    jTot = [{'E':float(h.E()), 'px':float(h.Px()), 'py':float(h.Py()), 'pz':float(h.Pz())}]
                    j.append(jTot[0])
                if type == 'tau1': 
                    nu = [{'E':float(sv.tau1nuE), 'px':float(sv.tau1nuPx), 'py':float(sv.tau1nuPy), 'pz':float(sv.tau1nuPz)}]
                    j.append(nu[0])
                elif type == 'tau2':
                    nu = [{'E':float(sv.tau2nuE), 'px':float(sv.tau2nuPx), 'py':float(sv.tau2nuPy), 'pz':float(sv.tau2nuPz)}]
                    j.append(nu[0])
                cand = [{'E':float(t.E()), 'px':float(t.Px()), 'py':float(t.Py()), 'pz':float(t.Pz())}]

                if self.debug: 
                    print("DM1: {0}, DM2: {1}".format(sv.tau1DM, sv.tau2DM))
                    print("cand: {0}".format(cand))
                fill = self.checkJSON(j, cand)
                #if (type == 'tau1' and sv.tau1DM == -1.) or (type == 'tau2' and sv.tau2DM == -1.):
                #    fill = False
                
                if fill: j.append(cand[0])
                isFirst = False
        #if len(j) < 3: j = []
        if self.checkJSONNaN(j): j = []
        if self.debug: print("Output json: {0}".format(j))

        return j

    def higg2json(self, svfit, tau1GenMatch, tau2GenMatch, tau1GenMatchPdgId, tau2GenMatchPdgId, svfitmet, genHiggs = None):
        j = []
        isFirst = True
        for idx, sv in enumerate(svfit):
            if not sv.PassBaseline or not sv.PassLepton: continue
            if (tau1GenMatch[idx] and tau2GenMatch[idx]) or (tau1GenMatchPdgId[idx] and tau2GenMatchPdgId[idx]) or genHiggs == None:
                if isFirst: h = self.makeTLorentzVector(sv if genHiggs == None else genHiggs, "gen")
                t1 = self.makeTLorentzVector(sv, 'tau1')
                t2 = self.makeTLorentzVector(sv, 'tau2')
                
                if isFirst:
                    jTot = [{'E':float(h.E()), 'px':float(h.Px()), 'py':float(h.Py()), 'pz':float(h.Pz())}]
                    j.append(jTot[0])
                    nu1 = [{'E':float(sv.tau1nuE), 'px':float(sv.tau1nuPx), 'py':float(sv.tau1nuPy), 'pz':float(sv.tau1nuPz)}]
                    nu2 = [{'E':float(sv.tau2nuE), 'px':float(sv.tau2nuPx), 'py':float(sv.tau2nuPy), 'pz':float(sv.tau2nuPz)}]
                    j.append(nu1[0])
                    j.append(nu2[0])
                cand1 = [{'E':float(t1.E()), 'px':float(t1.Px()), 'py':float(t1.Py()), 'pz':float(t1.Pz())}]
                cand2 = [{'E':float(t2.E()), 'px':float(t2.Px()), 'py':float(t2.Py()), 'pz':float(t2.Pz())}]

                if self.debug: 
                    print("DM1: {0}, DM2: {1}".format(sv.tau1DM, sv.tau2DM))
                    print("cand1: {0}, cand2: {1}".format(cand1, cand2))
                fill1 = self.checkJSON(j, cand1)
                fill2 = self.checkJSON(j, cand2)
                #if sv.tau1DM == -1.: 
                #    fill1 = False
                #if sv.tau2DM == -1.:
                #    fill2 = False
                
                if fill1: j.append(cand1[0])
                if fill2: j.append(cand2[0])
                isFirst = False
        if self.checkJSONNaN(j): j = []
        if self.debug: print("Output json: {0}".format(j))

        return j

    class TTreeReaderArrayWrapper:
        def __init__(self, ttarray):
            self.ttarray = ttarray

        def __iter__(self):
            for i in xrange(len(self.ttarray)):
                yield self.ttarray[i]
            return

    def HiggsGenMatch(self, event, CandEta, CandPhi, CandPdgId, match = 'GenPart'):
        genTauDaughters_list = []
        genTaupdgId_list = []
    
        eta              = np.fromiter(self.TTreeReaderArrayWrapper(event.GenPart_eta),              dtype=float)
        phi              = np.fromiter(self.TTreeReaderArrayWrapper(event.GenPart_phi),              dtype=float)
        genPartIdxMother = np.fromiter(self.TTreeReaderArrayWrapper(event.GenPart_genPartIdxMother), dtype=int)
        pdgId            = np.fromiter(self.TTreeReaderArrayWrapper(event.GenPart_pdgId),            dtype=int)
        statusFlags      = np.fromiter(self.TTreeReaderArrayWrapper(event.GenPart_statusFlags),      dtype=int)
    
        genTauDaughters_list = genParticleAssociation(genPartIdxMother, pdgId, statusFlags)
    
        genTauDaughters = np.array(genTauDaughters_list)
        genTaupdgId = np.array(genTaupdgId_list)
    
        if(len(genTauDaughters)):
            genTauDaughters_pdgid = pdgId[genTauDaughters]
            genTauDaughters_eta = eta[genTauDaughters]
            genTauDaughters_phi = phi[genTauDaughters]
        elif len(genTauDaughters) == 0 and "qcd" in self.processName.lower():
            genTauDaughters = np.linspace(0, len(eta), num = len(eta) + 1)
            genTauDaughters_pdgid = pdgId #np.array([])
            genTauDaughters_eta = eta #np.array([])
            genTauDaughters_phi = phi #np.array([])
        else:
            genTauDaughters_pdgid = np.array([])
            genTauDaughters_eta = np.array([])
            genTauDaughters_phi = np.array([])
    
        matches, candmatches, matchPdgID = deltaRMatch(CandPdgId, CandEta, CandPhi, genTauDaughters_pdgid, genTauDaughters_eta, genTauDaughters_phi)
    
        return matches, candmatches, matchPdgID, genTauDaughters

    def recursiveFindHiggs(self, startIdx, startPdgId, genpart):
        if startIdx < 0:
            return -1
    
        mom = genpart[startIdx].genPartIdxMother
    
        if mom < 0:
            return startIdx
        elif genpart[mom].pdgId != startPdgId:
            return mom
        else:
            return self.recursiveFindHiggs(mom, startPdgId, genpart)

    def findwhichTau(self, tau, genpart):
        idx = -1
        PT = -999
        for i, t in enumerate(tau):
            pdgid = genpart[t].pdgId
            pt = genpart[t].pt

            if pt >= PT:
                PT = pt
                idx = i

        return i

    def firstRealValue(self, index, value):
        if index < 0 or index + 1 > len(value): return -1

        if value[index] < 0:
            if index + 1 == len(value): return -1
            return self.firstRealValue(index + 1, value)
        else:
            return value[index]

    def findMothers(self, tau1, tau2, higgs1, higgs2, wt1, wt2, genpart):
        t1 = t2 = h = -1
        if len(wt1) == 0 and len(wt2) == 0: return t1, t2, h        
        if self.debug: print("tau1: {0}, wt1: {4}, tau2: {1}, wt2: {5}, h1: {2}, h2: {3}".format(tau1, tau2, higgs1, higgs2, wt1, wt2))

        
        if len(tau1) > 0: t1 = tau1[0]
        elif len(wt1) > 0: 
            m1 = genpart[int(wt1[0])].genPartIdxMother
            m2 = genpart[m1].genPartIdxMother
            t1 = m1
            higgs1.append(m2)

        if len(tau2) > 0: t2 = tau2[0]
        elif len(wt2) > 0: 
            m1 = genpart[int(wt2[0])].genPartIdxMother
            m2 = genpart[m1].genPartIdxMother
            t2 = m1
            higgs2.append(m2)

        #if len(wt1) == 0:
        #    h1 = genpart[genpart[wt2[0]].genPartIdxMother].genPartIdxMother

        if len(higgs1) > 0 or len(higgs2) > 0: 
            if len(higgs1) != 0: h = self.firstRealValue(0, higgs1)
            if higgs1 == -1 and len(higgs2) != 0: h = self.firstRealValue(0, higgs2)

        if self.debug: print("Index of (tau1, tau2, higgs) = ({0}, {1}, {2})".format(t1, t2, h))

        return t1, t2, h

    def analyze(self, event):
        ## Getting objects
        svfit      = Collection(event, "SVFit")
        svfitmet   = Collection(event, "SVFitMET")
        genpart    = Collection(event, "GenPart")

        tau1Eta = np.fromiter(self.TTreeReaderArrayWrapper(event.SVFit_tau1Eta), dtype=float)
        tau1Phi = np.fromiter(self.TTreeReaderArrayWrapper(event.SVFit_tau1Phi), dtype=float)
        tau1pdgId = np.fromiter(self.TTreeReaderArrayWrapper(event.SVFit_tau1pdgId), dtype=float)
        tau1pdgIdTranslate = [tauDecayPdgId[TauDecay[t]] for t in tau1pdgId]
        tau1GenMatch, tau1GenPartMatch, tau1GenMatchPdgId, whichTau1 = self.HiggsGenMatch(event, tau1Eta, tau1Phi, tau1pdgIdTranslate, self.match)

        tau2Eta = np.fromiter(self.TTreeReaderArrayWrapper(event.SVFit_tau2Eta), dtype=float)
        tau2Phi = np.fromiter(self.TTreeReaderArrayWrapper(event.SVFit_tau2Phi), dtype=float)
        tau2pdgId = np.fromiter(self.TTreeReaderArrayWrapper(event.SVFit_tau2pdgId), dtype=float)
        tau2pdgIdTranslate = [tauDecayPdgId[TauDecay[t]] for t in tau2pdgId]
        tau2GenMatch, tau2GenPartMatch, tau2GenMatchPdgId, whichTau2 = self.HiggsGenMatch(event, tau2Eta, tau2Phi, tau2pdgIdTranslate, self.match)

        if self.debug:
            print("tau1: {0} {1} {2} {3}".format(tau1GenMatch, tau1pdgIdTranslate, whichTau1, tau1GenPartMatch))
            print("tau2: {0} {1} {2} {3}".format(tau2GenMatch, tau2pdgIdTranslate, whichTau2, tau2GenPartMatch))

        tauIdx1 = tauIdx2 = matchtau1 = matchtau2 = higgs1Idx = higgs2Idx = []
        for idx, (gD1, gD2, t1, t2) in enumerate(zip(whichTau1, whichTau2, tau1GenPartMatch, tau2GenPartMatch)):
            idx1 = self.recursiveFindHiggs(int(gD1), genpart[int(gD1)].pdgId, genpart) if bool(t1) else -1
            idx2 = self.recursiveFindHiggs(int(gD2), genpart[int(gD2)].pdgId, genpart) if bool(t2) else -1
            if idx1 >= 0:
                tauIdx1.append(idx1)
                higgs1Idx.append(genpart[idx1].genPartIdxMother)
                if self.debug: print("Matched Higgs Index: {0} --> tau1: {1} pdgId: {2}".format(genpart[idx1].genPartIdxMother, idx1, genpart[idx1].pdgId))
            if idx2 >= 0:
                tauIdx2.append(idx2)
                higgs2Idx.append(genpart[idx2].genPartIdxMother)
                if self.debug: print("Matched Higgs Index: {0} --> tau2: {1} pdgId: {2}".format(genpart[idx2].genPartIdxMother, idx2, genpart[idx2].pdgId))
        if self.debug: print("Higgs Vec: {0} = {1}, tau1 Vec: {2}, tau2 Vec: {3}".format(higgs1Idx, higgs2Idx, tauIdx1, tauIdx2))

        j = []
        t1, t2, h = self.findMothers(tauIdx1, tauIdx2, higgs1Idx, higgs2Idx, whichTau1, whichTau2, genpart)
        if t1 >= 0: 
            j = self.tau2json(svfit, tau1GenMatch, tau1GenPartMatch, svfitmet, genpart[t1], 'tau1')
            if len(j) > 2: self.fout.write((json.dumps(j, sort_keys=False)+'\n').encode('utf-8'))
        if t2 >= 0: 
            j = self.tau2json(svfit, tau2GenMatch, tau2GenPartMatch, svfitmet, genpart[2], 'tau2')
            if len(j) > 2: self.fout.write((json.dumps(j, sort_keys=False)+'\n').encode('utf-8'))
        if h >= 0:
            j = self.higg2json(svfit, tau1GenMatch, tau2GenMatch, tau1GenMatchPdgId, tau2GenMatchPdgId, svfitmet, genpart[h])
            if len(j) > 3: self.fgenout.write((json.dumps(j, sort_keys=False)+'\n').encode('utf-8'))

        return True
        
        
 # define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
