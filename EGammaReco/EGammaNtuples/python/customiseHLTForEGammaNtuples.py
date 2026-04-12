from RecoEgamma.EgammaHLTProducers.hltEgammaHLTExtra_cfi import hltEgammaHLTExtra
from EGammaReco.EGammaNtuples.EGammaNtupler_cfi import EGammaNtuplesSequence, EGammaNtuples
import FWCore.ParameterSet.Config as cms

def customiseHLTForPhotonNtuples(process):

    process.EGammaNtuples = EGammaNtuples.clone(
        pType = "pho",
        ptThreshold = 20.0
    )

    from Configuration.ProcessModifiers.ticl_v5_cff import ticl_v5
    ticl_v5.toModify(process.EGammaNtuples, scHGCalL1Seeded = cms.InputTag('hltTiclEGammaSuperClusterProducerUnseeded'))

    process.FEVTDEBUGHLToutput_step = cms.EndPath(
        process.FEVTDEBUGHLToutput
        * process.EGammaNtuples)
    return process

def customiseHLTForElectronNtuples(process):

    process.EGammaNtuples = EGammaNtuples.clone(
        pType = "ele",
        ptThreshold =  10.0
    )

    from Configuration.ProcessModifiers.ticl_v5_cff import ticl_v5
    ticl_v5.toModify(process.EGammaNtuples, scHGCalL1Seeded = cms.InputTag('hltTiclEGammaSuperClusterProducerUnseeded'))

    process.FEVTDEBUGHLToutput_step = cms.EndPath(
        process.FEVTDEBUGHLToutput 
        * process.EGammaNtuples)
    return process
