import FWCore.ParameterSet.Config as cms
from RecoEgamma.EgammaHLTProducers.hltEgammaHLTExtra_cfi import hltEgammaHLTExtra


# EGammaNtuples module
EGammaNtuples = cms.EDAnalyzer(
    "EGammaNtuples",
    genParticles = cms.InputTag("genParticles"),
    scBarrelL1Seeded = cms.InputTag("hltParticleFlowSuperClusterECALUnseeded", "particleFlowSuperClusterECALBarrel"),
    scHGCalL1Seeded = cms.InputTag("hltParticleFlowSuperClusterHGCalFromTICLUnseeded", ""),
    ebRecHits = cms.InputTag("hltEgammaHLTExtra", "EcalRecHitsEB"),
    HGCeeRecHits = cms.InputTag( "hltHGCalRecHit", "HGCEERecHits", "HLTX"),
    HGChebRecHits = cms.InputTag( "hltHGCalRecHit", "HGCHEBRecHits", "HLTX"),
    HGChefRecHits = cms.InputTag( "hltHGCalRecHit", "HGCHEFRecHits", "HLTX"),
    sigmaIEtaIEta = cms.InputTag("hltEgammaClusterShapeUnseeded", "sigmaIEtaIEta5x5"),
    sigmaIPhiIPhi = cms.InputTag("hltEgammaClusterShapeUnseeded", "sigmaIPhiIPhi5x5"),
    sigmaIEtaIEtaNoiseCleaned = cms.InputTag("hltEgammaClusterShapeUnseeded", "sigmaIEtaIEta5x5NoiseCleaned"),
    sigmaIPhiIPhiNoiseCleaned = cms.InputTag("hltEgammaClusterShapeUnseeded", "sigmaIPhiIPhi5x5NoiseCleaned"),
    nrHitsEB1GeV = cms.InputTag("hltEgammaHLTExtra", "countEcalRecHitsEcalRecHitsEBThres1GeV"),
    nrHitsEE1GeV = cms.InputTag("hltEgammaHLTExtra", "countEcalRecHitsHGCalRecHitsThres1GeV"),
    eGammaCandidates = cms.InputTag("hltEgammaCandidatesUnseeded"),
    vertices = cms.InputTag("offlinePrimaryVertices"),
    pfHGCALRecHits = cms.InputTag("particleFlowRecHitHGC"),
    trackster = cms.InputTag("hltTiclCandidate", "", "HLTX"),
    layerClusters = cms.InputTag("hltMergeLayerClusters", "", "HLTX"),
    timelayerClusters = cms.InputTag("hltMergeLayerClusters","timeLayerCluster", "HLTX"),
    pType = cms.string("ele"),
    ptThreshold = cms.double(10.0),
    genMatchDR = cms.double(0.1),
    maxDR      = cms.double(0.4),
)

# Sequence bundling both
EGammaNtuplesSequence = cms.Sequence(
    hltEgammaHLTExtra *
    EGammaNtuples
)

