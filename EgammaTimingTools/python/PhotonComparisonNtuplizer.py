import os

import FWCore.ParameterSet.Config as cms
from Configuration.AlCa.GlobalTag import GlobalTag
from Configuration.StandardSequences.Eras import eras

process = cms.Process("PhotonComparisonNtuplizer", eras.Phase2C17I13M9)

process.load("Configuration.Geometry.GeometryExtended2026D110Reco_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag = GlobalTag(process.GlobalTag, "auto:phase2_realistic_T21", "")

import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing("analysis")
options.register("sourceFile", "", VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "File containing list of input files")
options.register("photonLabel", "photonsHGC", VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Offline photon collection")
options.register("offlineProcess", "reRECO", VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Offline process name")
options.register("onlineLabel", "hltEgammaHLTExtra:Unseeded:HLT", VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Online EGammaObject collection")
options.register("onlineCandidateLabel", "hltEgammaCandidatesUnseeded::HLTX", VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Online RecoEcalCandidate collection")
options.register("onlineTracksterLabel", "hltTiclCandidate::HLTX", VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Online HLT trackster collection")
options.register("outDir", "", VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Output directory")
options.register("outFileNumber", -1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "Output file number")
options.parseArguments()

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))

if len(options.sourceFile):
  source_file = options.sourceFile

file_names = []
if len(options.inputFiles):
  file_names = options.inputFiles
else:
  with open(source_file) as handle:
    file_names = handle.readlines()

file_names = [name.strip() for name in file_names if name.strip() and name[0] != "#"]
for index, file_name in enumerate(file_names):
  if "file:" not in file_name and "root:" not in file_name:
    file_names[index] = "file:%s" % (file_name)

process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(file_names))

out_file_suffix = ""
if options.outFileNumber >= 0:
  out_file_suffix = "_%d" % (options.outFileNumber)

out_file = "comparisonNtuple%s.root" % (out_file_suffix)
if len(options.outDir):
  os.system("mkdir -p %s" % (options.outDir))
  out_file = "%s/%s" % (options.outDir, out_file)

process.TFileService = cms.Service("TFileService", fileName=cms.string(out_file))

process.ntuplizer = cms.EDAnalyzer(
    "PhotonComparisonNtuplizer",
    photonProducer=cms.InputTag(options.photonLabel, "", options.offlineProcess),
    onlineProducer=cms.InputTag(*options.onlineLabel.split(":")),
    onlineCandidateProducer=cms.InputTag(*options.onlineCandidateLabel.split(":")),
    onlineTracksterSrc=cms.InputTag(*options.onlineTracksterLabel.split(":")),
    onlineLayerClusterSrc=cms.InputTag("hltMergeLayerClusters", "", "HLTX"),
    onlineTimeLayerClusterSrc=cms.InputTag("hltMergeLayerClusters", "timeLayerCluster", "HLTX"),
    onlineRecHitsEE_Src=cms.InputTag("hltHGCalRecHit", "HGCEERecHits", "HLTX"),
    onlineRecHitsFH_Src=cms.InputTag("hltHGCalRecHit", "HGCHEFRecHits", "HLTX"),
    onlineRecHitsBH_Src=cms.InputTag("hltHGCalRecHit", "HGCHEBRecHits", "HLTX"),
    trackProducer=cms.InputTag("generalTracks", "", options.offlineProcess),
    mtdt0=cms.InputTag("tofPID", "t0", options.offlineProcess),
    mtdSigmat0=cms.InputTag("tofPID", "sigmat0", options.offlineProcess),
    mtdTrkQualMVA=cms.InputTag("mtdTrackQualityMVA", "mtdQualMVA", options.offlineProcess),
    BeamspotProducer=cms.InputTag("offlineBeamSpot", "", options.offlineProcess),
    vertexProducer=cms.InputTag("offlineSlimmedPrimaryVertices4D", "", options.offlineProcess),
    genParticles=cms.InputTag("prunedGenParticles"),
    tracksterSrc=cms.InputTag("ticlTrackstersMerge", "", options.offlineProcess),
    LayerClusterSrc=cms.InputTag("hgcalMergeLayerClusters", "", options.offlineProcess),
    RecHitsEE_Src=cms.InputTag("HGCalRecHit", "HGCEERecHits", options.offlineProcess),
    RecHitsFH_Src=cms.InputTag("HGCalRecHit", "HGCHEFRecHits", options.offlineProcess),
    RecHitsBH_Src=cms.InputTag("HGCalRecHit", "HGCHEBRecHits", options.offlineProcess),
    ptMin=cms.double(1.0),
    offlinePtMin=cms.double(15.0),
    onlinePtMin=cms.double(15.0),
    intRadiusBarrel=cms.double(0.01),
    intRadiusEndcap=cms.double(0.01),
    stripBarrel=cms.double(0.01),
    stripEndcap=cms.double(0.01),
    extRadius=cms.double(0.3),
    gen_deltaR=cms.double(0.1),
    isMC=cms.bool(True),
)

process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.StandardSequences.L1Reco_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.RecoSim_cff")
process.load("CalibTracker.SiStripESProducers.SiStripLorentzAngleDepESProducer_cfi")

process.p = cms.Path(process.ntuplizer)
process.schedule = cms.Schedule(process.p)
