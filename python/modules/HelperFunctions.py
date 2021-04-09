import math
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

def hasBit(value,bit):
    """Check if i'th bit is set to 1, i.e. binary of 2^(i-1),
    from the right to the left, starting from position i=0."""
    # https://cms-nanoaod-integration.web.cern.ch/integration/master-102X/mc102X_doc.html#GenPart
    # Gen status flags, stored bitwise, are:
    #    0: isPrompt,                          8: fromHardProcess,
    #    1: isDecayedLeptonHadron,             9: isHardProcessTauDecayProduct,
    #    2: isTauDecayProduct,                10: isDirectHardProcessTauDecayProduct,
    #    3: isPromptTauDecayProduct,          11: fromHardProcessBeforeFSR,
    #    4: isDirectTauDecayProduct,          12: isFirstCopy,
    #    5: isDirectPromptTauDecayProduct,    13: isLastCopy,
    #    6: isDirectHadronDecayProduct,       14: isLastCopyBeforeFSR
    #    7: isHardProcess,
    ###return bin(value)[-bit-1]=='1'
    ###return format(value,'b').zfill(bit+1)[-bit-1]=='1'
    return (value & (1 << bit))>0

def isA(particleID, p):
    return abs(p) == particleID

def addFourVec(obj):
    tot = ROOT.TLorentzVector()
    v1 = ROOT.TLorentzVector()
    v2 = ROOT.TLorentzVector()

    if(len(obj) > 0): v1.SetPtEtaPhiM(obj[0].pt, obj[0].eta, obj[0].phi, obj[0].mass)
    else: v1.SetPtEtaPhiM(0, 0, 0, 0)
    if(len(obj) > 1): v2.SetPtEtaPhiM(obj[1].pt, obj[1].eta, obj[1].phi, obj[1].mass)
    else: v2.SetPtEtaPhiM(0, 0, 0, 0)
    tot = (v1 + v2)
    return tot
