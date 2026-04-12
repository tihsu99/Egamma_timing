import os

import FWCore.ParameterSet.Config as cms
from Configuration.AlCa.GlobalTag import GlobalTag
from Configuration.StandardSequences.Eras import eras

process = cms.Process("ElectronComparisonNtuplizer", eras.Phase2C17I13M9)

process.load("Configuration.Geometry.GeometryExtended2026D110Reco_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag = GlobalTag(process.GlobalTag, "auto:phase2_realistic_T21", "")

import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing("analysis")
options.register("sourceFile", "", VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "File containing list of input files")
options.register("electronLabel", "ecalDrivenGsfElectronsHGC", VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Offline electron collection")
options.register("onlineLabel", "hltEgammaHLTExtra:Unseeded:HLT", VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Online EGammaObject collection")
options.register("outDir", "", VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Output directory")
options.register("outFileNumber", -1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "Output file number")
options.register("useAOD", 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "Unused compatibility flag")
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
    "ElectronComparisonNtuplizer",
    electronProducer=cms.InputTag(options.electronLabel),
    onlineProducer=cms.InputTag(*options.onlineLabel.split(":")),
    trackProducer=cms.InputTag("generalTracks"),
    mtdt0=cms.InputTag("tofPID:t0"),
    mtdSigmat0=cms.InputTag("tofPID:sigmat0"),
    mtdTrkQualMVA=cms.InputTag("mtdTrackQualityMVA:mtdQualMVA"),
    BeamspotProducer=cms.InputTag("offlineBeamSpot"),
    vertexProducer=cms.InputTag("offlineSlimmedPrimaryVertices4D"),
    genParticles=cms.InputTag("prunedGenParticles"),
    tracksterSrc=cms.InputTag("ticlTrackstersMerge"),
    LayerClusterSrc=cms.InputTag("hgcalMergeLayerClusters"),
    RecHitsEE_Src=cms.InputTag("HGCalRecHit", "HGCEERecHits"),
    RecHitsFH_Src=cms.InputTag("HGCalRecHit", "HGCHEFRecHits"),
    RecHitsBH_Src=cms.InputTag("HGCalRecHit", "HGCHEBRecHits"),
    ptMin=cms.double(1.0),
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
