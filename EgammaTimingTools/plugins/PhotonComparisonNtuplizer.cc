#include "ComparisonNtuplizerUtils.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaTrackSelector.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include <Math/VectorUtil.h>
#include <TTree.h>

#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

class PhotonComparisonNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit PhotonComparisonNtuplizer(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void setDzOption(const std::string& option) {
    if (!option.compare("dz")) {
      dzOption_ = egammaisolation::EgammaTrackSelector::dz;
    } else if (!option.compare("vz")) {
      dzOption_ = egammaisolation::EgammaTrackSelector::vz;
    } else if (!option.compare("bs")) {
      dzOption_ = egammaisolation::EgammaTrackSelector::bs;
    } else if (!option.compare("vtx")) {
      dzOption_ = egammaisolation::EgammaTrackSelector::vtx;
    } else {
      dzOption_ = egammaisolation::EgammaTrackSelector::dz;
    }
  }

  void resetOfflineVectors();
  void resetOnlineBranches();
  void resetOnlineVectors();
  void fillOnlineBranches(const reco::RecoEcalCandidate&,
                         int candidateIndex,
                         const trigger::EgammaObject*,
                         int objectIndex,
                         const egammaTiming::GenMatchResult&);
  void fillGenBranches(int truthClass,
                       int genIndex,
                       const reco::GenParticle&,
                       const std::pair<int, double>& offlineMatch,
                       const std::pair<int, double>& onlineMatch,
                        const std::vector<reco::Photon>&,
                        const std::vector<reco::RecoEcalCandidate>&,
                        const std::vector<trigger::EgammaObject>&);

  int nEvent_;
  int nRun_;
  int nLumi_;
  int vtxN_;
  int sigTrkIdx_;
  int genIdx_;
  double sigTrkTime_;
  double sigTrkTimeErr_;
  double sigTrkMtdMva_;
  double Pho_pt_;
  double Pho_eta_;
  double Pho_phi_;
  double PVTime_;
  double PVTimeErr_;
  double matchedToGenPho_;
  double gen_pt_;
  double gen_eta_;
  double gen_phi_;
  double dr_gen_offline_;

  std::vector<double> Track_pt_;
  std::vector<double> Track_eta_;
  std::vector<double> Track_phi_;
  std::vector<double> Track_Time_;
  std::vector<double> Track_TimeErr_;
  std::vector<double> Track_MtdMva_;
  std::vector<double> Track_dz_;
  std::vector<double> Track_dxy_;

  std::vector<double> Trackster_pt_;
  std::vector<double> Trackster_eta_;
  std::vector<double> Trackster_phi_;
  std::vector<double> Trackster_dr_;
  std::vector<double> Trackster_Time_;
  std::vector<double> Trackster_TimeErr_;
  std::vector<int> Trackster_TrackIdx_;
  std::vector<double> Trackster_Trackpt_;
  std::vector<double> Trackster_Tracketa_;
  std::vector<double> Trackster_Trackphi_;
  std::vector<bool> Trackster_isSeed_;
  std::vector<bool> Trackster_inCluster_;

  std::vector<double> RecHit_dr_;
  std::vector<double> RecHit_eta_;
  std::vector<double> RecHit_phi_;
  std::vector<int> RecHit_layer_;
  std::vector<double> RecHit_energy_;
  std::vector<double> RecHit_time_;
  std::vector<double> RecHit_time_error_;
  std::vector<bool> RecHit_isSeeded_;
  std::vector<bool> RecHit_inCluster_;

  int online_genIdx_;
  int online_objIdx_;
  double online_matchedToGen_;
  double online_dr_gen_;
  double online_gen_pt_;
  double online_gen_eta_;
  double online_gen_phi_;
  double online_pt_;
  double online_eta_;
  double online_phi_;
  double online_energy_;
  double online_ecalPFIsol_;
  double online_hcalPFIsol_;
  double online_hgcalPFIsol_;
  double online_trkIsol_;
  double online_r9_;
  double online_sigmaIEtaIEta_;
  double online_rVar_;
  double online_hForHoverE_;
  int online_candIdx_;

  std::vector<double> online_Trackster_pt_;
  std::vector<double> online_Trackster_eta_;
  std::vector<double> online_Trackster_phi_;
  std::vector<double> online_Trackster_dr_;
  std::vector<double> online_Trackster_Time_;
  std::vector<double> online_Trackster_TimeErr_;
  std::vector<bool> online_Trackster_isSeed_;
  std::vector<bool> online_Trackster_inCluster_;

  std::vector<double> online_LayerCluster_energy_;
  std::vector<double> online_LayerCluster_eta_;
  std::vector<double> online_LayerCluster_phi_;
  std::vector<double> online_LayerCluster_dr_;
  std::vector<double> online_LayerCluster_Time_;
  std::vector<double> online_LayerCluster_TimeErr_;
  std::vector<bool> online_LayerCluster_isSeed_;
  std::vector<bool> online_LayerCluster_inCluster_;

  std::vector<double> online_RecHit_dr_;
  std::vector<double> online_RecHit_eta_;
  std::vector<double> online_RecHit_phi_;
  std::vector<int> online_RecHit_layer_;
  std::vector<double> online_RecHit_energy_;
  std::vector<double> online_RecHit_time_;
  std::vector<double> online_RecHit_time_error_;
  std::vector<bool> online_RecHit_isSeeded_;
  std::vector<bool> online_RecHit_inCluster_;

  int genTree_idx_;
  int genTree_truth_;
  int genTree_offline_idx_;
  int genTree_online_idx_;
  double genTree_pt_;
  double genTree_eta_;
  double genTree_phi_;
  double genTree_offline_exists_;
  double genTree_offline_pt_;
  double genTree_offline_eta_;
  double genTree_offline_phi_;
  double genTree_offline_dr_;
  double genTree_online_exists_;
  double genTree_online_pt_;
  double genTree_online_eta_;
  double genTree_online_phi_;
  double genTree_online_dr_;
  double genTree_online_hgcalIso_;

  const edm::EDGetTokenT<std::vector<reco::Photon>> photonProducer_;
  const edm::EDGetTokenT<std::vector<trigger::EgammaObject>> onlineProducer_;
  const edm::EDGetTokenT<std::vector<reco::RecoEcalCandidate>> onlineCandidateProducer_;
  edm::EDGetTokenT<std::vector<ticl::Trackster>> onlineTracksterToken_;
  edm::EDGetTokenT<std::vector<reco::CaloCluster>> onlineLayerClusterToken_;
  const edm::EDGetTokenT<edm::ValueMap<std::pair<float, float>>> onlineLayerClusterTimeToken_;
  const edm::EDGetTokenT<HGCRecHitCollection> onlineHitsEE_;
  const edm::EDGetTokenT<HGCRecHitCollection> onlineHitsFH_;
  const edm::EDGetTokenT<HGCRecHitCollection> onlineHitsBH_;
  const edm::EDGetTokenT<reco::TrackCollection> trackProducer_;
  const edm::EDGetTokenT<reco::BeamSpot> beamspotProducer_;
  const edm::EDGetTokenT<std::vector<reco::Vertex>> vertexProducer_;
  const edm::EDGetTokenT<edm::View<reco::GenParticle>> genParticleProducer_;
  edm::EDGetTokenT<std::vector<ticl::Trackster>> tracksterToken_;
  edm::EDGetTokenT<std::vector<reco::CaloCluster>> layerClusterToken_;
  const edm::EDGetTokenT<HGCRecHitCollection> hitsEE_;
  const edm::EDGetTokenT<HGCRecHitCollection> hitsFH_;
  const edm::EDGetTokenT<HGCRecHitCollection> hitsBH_;
  const edm::EDGetTokenT<edm::ValueMap<float>> mtdt0Token_;
  const edm::EDGetTokenT<edm::ValueMap<float>> mtdSigmat0Token_;
  const edm::EDGetTokenT<edm::ValueMap<float>> mtdTrkQualMVAToken_;
  const edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometryToken_;

  const double ptMin_;
  const double offlinePtMin_;
  const double onlinePtMin_;
  const double intRadiusBarrel_;
  const double intRadiusEndcap_;
  const double stripBarrel_;
  const double stripEndcap_;
  const double extRadius_;
  const double gen_deltaR_;
  int dzOption_;
  const bool isMC_;

  hgcal::RecHitTools recHitTools_;

  TTree* offlineTree_;
  TTree* onlineTree_;
  TTree* genTree_;
};

PhotonComparisonNtuplizer::PhotonComparisonNtuplizer(const edm::ParameterSet& config)
    : photonProducer_{consumes(config.getParameter<edm::InputTag>("photonProducer"))},
      onlineProducer_{consumes(config.getParameter<edm::InputTag>("onlineProducer"))},
      onlineCandidateProducer_{consumes(config.getParameter<edm::InputTag>("onlineCandidateProducer"))},
      onlineTracksterToken_{consumes(config.getParameter<edm::InputTag>("onlineTracksterSrc"))},
      onlineLayerClusterToken_{consumes(config.getParameter<edm::InputTag>("onlineLayerClusterSrc"))},
      onlineLayerClusterTimeToken_{consumes(config.getParameter<edm::InputTag>("onlineTimeLayerClusterSrc"))},
      onlineHitsEE_{consumes(config.getParameter<edm::InputTag>("onlineRecHitsEE_Src"))},
      onlineHitsFH_{consumes(config.getParameter<edm::InputTag>("onlineRecHitsFH_Src"))},
      onlineHitsBH_{consumes(config.getParameter<edm::InputTag>("onlineRecHitsBH_Src"))},
      trackProducer_{consumes(config.getParameter<edm::InputTag>("trackProducer"))},
      beamspotProducer_{consumes(config.getParameter<edm::InputTag>("BeamspotProducer"))},
      vertexProducer_{consumes(config.getParameter<edm::InputTag>("vertexProducer"))},
      genParticleProducer_{consumes<edm::View<reco::GenParticle>>(config.getParameter<edm::InputTag>("genParticles"))},
      tracksterToken_{consumes(config.getParameter<edm::InputTag>("tracksterSrc"))},
      layerClusterToken_{consumes(config.getParameter<edm::InputTag>("LayerClusterSrc"))},
      hitsEE_{consumes(config.getParameter<edm::InputTag>("RecHitsEE_Src"))},
      hitsFH_{consumes(config.getParameter<edm::InputTag>("RecHitsFH_Src"))},
      hitsBH_{consumes(config.getParameter<edm::InputTag>("RecHitsBH_Src"))},
      mtdt0Token_{consumes(config.getParameter<edm::InputTag>("mtdt0"))},
      mtdSigmat0Token_{consumes(config.getParameter<edm::InputTag>("mtdSigmat0"))},
      mtdTrkQualMVAToken_{consumes(config.getParameter<edm::InputTag>("mtdTrkQualMVA"))},
      caloGeometryToken_(esConsumes()),
      ptMin_{config.getParameter<double>("ptMin")},
      offlinePtMin_{config.getParameter<double>("offlinePtMin")},
      onlinePtMin_{config.getParameter<double>("onlinePtMin")},
      intRadiusBarrel_{config.getParameter<double>("intRadiusBarrel")},
      intRadiusEndcap_{config.getParameter<double>("intRadiusEndcap")},
      stripBarrel_{config.getParameter<double>("stripBarrel")},
      stripEndcap_{config.getParameter<double>("stripEndcap")},
      extRadius_{config.getParameter<double>("extRadius")},
      gen_deltaR_{config.getParameter<double>("gen_deltaR")},
      isMC_{config.getParameter<bool>("isMC")} {
  setDzOption("vz");
  usesResource(TFileService::kSharedResource);

  edm::Service<TFileService> fs;
  offlineTree_ = fs->make<TTree>("offlineTree", "offlineTree");
  onlineTree_ = fs->make<TTree>("onlineTree", "onlineTree");
  genTree_ = fs->make<TTree>("genTree", "genTree");

  offlineTree_->Branch("nEvent", &nEvent_);
  offlineTree_->Branch("nRun", &nRun_);
  offlineTree_->Branch("nLumi", &nLumi_);
  offlineTree_->Branch("vtxN", &vtxN_);
  offlineTree_->Branch("Pho_pt", &Pho_pt_);
  offlineTree_->Branch("Pho_eta", &Pho_eta_);
  offlineTree_->Branch("Pho_phi", &Pho_phi_);
  offlineTree_->Branch("matchedToGenPho", &matchedToGenPho_);
  offlineTree_->Branch("genIdx", &genIdx_);
  offlineTree_->Branch("gen_pt", &gen_pt_);
  offlineTree_->Branch("gen_eta", &gen_eta_);
  offlineTree_->Branch("gen_phi", &gen_phi_);
  offlineTree_->Branch("dr_gen", &dr_gen_offline_);
  offlineTree_->Branch("PV_Time", &PVTime_);
  offlineTree_->Branch("PV_TimeErr", &PVTimeErr_);
  offlineTree_->Branch("sigTrkIdx", &sigTrkIdx_);
  offlineTree_->Branch("sigTrkTime", &sigTrkTime_);
  offlineTree_->Branch("sigTrkTimeErr", &sigTrkTimeErr_);
  offlineTree_->Branch("sigTrkMtdMva", &sigTrkMtdMva_);
  offlineTree_->Branch("Track_pt", &Track_pt_);
  offlineTree_->Branch("Track_eta", &Track_eta_);
  offlineTree_->Branch("Track_phi", &Track_phi_);
  offlineTree_->Branch("Track_Time", &Track_Time_);
  offlineTree_->Branch("Track_TimeErr", &Track_TimeErr_);
  offlineTree_->Branch("Track_MtdMva", &Track_MtdMva_);
  offlineTree_->Branch("Track_dz", &Track_dz_);
  offlineTree_->Branch("Track_dxy", &Track_dxy_);
  offlineTree_->Branch("Trackster_dr", &Trackster_dr_);
  offlineTree_->Branch("Trackster_isSeed", &Trackster_isSeed_);
  offlineTree_->Branch("Trackster_inCluster", &Trackster_inCluster_);
  offlineTree_->Branch("Trackster_pt", &Trackster_pt_);
  offlineTree_->Branch("Trackster_eta", &Trackster_eta_);
  offlineTree_->Branch("Trackster_phi", &Trackster_phi_);
  offlineTree_->Branch("Trackster_Time", &Trackster_Time_);
  offlineTree_->Branch("Trackster_TimeErr", &Trackster_TimeErr_);
  offlineTree_->Branch("Trackster_TrackIdx", &Trackster_TrackIdx_);
  offlineTree_->Branch("Trackster_Trackpt", &Trackster_Trackpt_);
  offlineTree_->Branch("Trackster_Tracketa", &Trackster_Tracketa_);
  offlineTree_->Branch("Trackster_Trackphi", &Trackster_Trackphi_);
  offlineTree_->Branch("HGCRecHit_dr", &RecHit_dr_);
  offlineTree_->Branch("HGCRecHit_eta", &RecHit_eta_);
  offlineTree_->Branch("HGCRecHit_phi", &RecHit_phi_);
  offlineTree_->Branch("HGCRecHit_layer", &RecHit_layer_);
  offlineTree_->Branch("HGCRecHit_energy", &RecHit_energy_);
  offlineTree_->Branch("HGCRecHit_Time", &RecHit_time_);
  offlineTree_->Branch("HGCRecHit_TimeError", &RecHit_time_error_);
  offlineTree_->Branch("HGCRecHit_isSeed", &RecHit_isSeeded_);
  offlineTree_->Branch("HGCRecHit_inCluster", &RecHit_inCluster_);

  onlineTree_->Branch("nEvent", &nEvent_);
  onlineTree_->Branch("nRun", &nRun_);
  onlineTree_->Branch("nLumi", &nLumi_);
  onlineTree_->Branch("Pho_pt", &online_pt_);
  onlineTree_->Branch("Pho_eta", &online_eta_);
  onlineTree_->Branch("Pho_phi", &online_phi_);
  onlineTree_->Branch("Pho_energy", &online_energy_);
  onlineTree_->Branch("matchedToGenPho", &online_matchedToGen_);
  onlineTree_->Branch("genIdx", &online_genIdx_);
  onlineTree_->Branch("gen_pt", &online_gen_pt_);
  onlineTree_->Branch("gen_eta", &online_gen_eta_);
  onlineTree_->Branch("gen_phi", &online_gen_phi_);
  onlineTree_->Branch("dr_gen", &online_dr_gen_);
  onlineTree_->Branch("objIdx", &online_objIdx_);
  onlineTree_->Branch("candIdx", &online_candIdx_);
  onlineTree_->Branch("ecalPFIsol_default", &online_ecalPFIsol_);
  onlineTree_->Branch("hcalPFIsol_default", &online_hcalPFIsol_);
  onlineTree_->Branch("hgcalPFIsol_default", &online_hgcalPFIsol_);
  onlineTree_->Branch("trkIsol_default", &online_trkIsol_);
  onlineTree_->Branch("r9_default", &online_r9_);
  onlineTree_->Branch("sigmaIEtaIEta_default", &online_sigmaIEtaIEta_);
  onlineTree_->Branch("rVar_default", &online_rVar_);
  onlineTree_->Branch("hForHoverE_default", &online_hForHoverE_);
  onlineTree_->Branch("Trackster_dr", &online_Trackster_dr_);
  onlineTree_->Branch("Trackster_isSeed", &online_Trackster_isSeed_);
  onlineTree_->Branch("Trackster_inCluster", &online_Trackster_inCluster_);
  onlineTree_->Branch("Trackster_pt", &online_Trackster_pt_);
  onlineTree_->Branch("Trackster_eta", &online_Trackster_eta_);
  onlineTree_->Branch("Trackster_phi", &online_Trackster_phi_);
  onlineTree_->Branch("Trackster_Time", &online_Trackster_Time_);
  onlineTree_->Branch("Trackster_TimeErr", &online_Trackster_TimeErr_);
  onlineTree_->Branch("LayerCluster_dr", &online_LayerCluster_dr_);
  onlineTree_->Branch("LayerCluster_energy", &online_LayerCluster_energy_);
  onlineTree_->Branch("LayerCluster_eta", &online_LayerCluster_eta_);
  onlineTree_->Branch("LayerCluster_phi", &online_LayerCluster_phi_);
  onlineTree_->Branch("LayerCluster_Time", &online_LayerCluster_Time_);
  onlineTree_->Branch("LayerCluster_TimeErr", &online_LayerCluster_TimeErr_);
  onlineTree_->Branch("LayerCluster_isSeed", &online_LayerCluster_isSeed_);
  onlineTree_->Branch("LayerCluster_inCluster", &online_LayerCluster_inCluster_);
  onlineTree_->Branch("HGCRecHit_dr", &online_RecHit_dr_);
  onlineTree_->Branch("HGCRecHit_eta", &online_RecHit_eta_);
  onlineTree_->Branch("HGCRecHit_phi", &online_RecHit_phi_);
  onlineTree_->Branch("HGCRecHit_layer", &online_RecHit_layer_);
  onlineTree_->Branch("HGCRecHit_energy", &online_RecHit_energy_);
  onlineTree_->Branch("HGCRecHit_Time", &online_RecHit_time_);
  onlineTree_->Branch("HGCRecHit_TimeError", &online_RecHit_time_error_);
  onlineTree_->Branch("HGCRecHit_isSeed", &online_RecHit_isSeeded_);
  onlineTree_->Branch("HGCRecHit_inCluster", &online_RecHit_inCluster_);

  genTree_->Branch("nEvent", &nEvent_);
  genTree_->Branch("nRun", &nRun_);
  genTree_->Branch("nLumi", &nLumi_);
  genTree_->Branch("genIdx", &genTree_idx_);
  genTree_->Branch("truthClass", &genTree_truth_);
  genTree_->Branch("gen_pt", &genTree_pt_);
  genTree_->Branch("gen_eta", &genTree_eta_);
  genTree_->Branch("gen_phi", &genTree_phi_);
  genTree_->Branch("offline_exists", &genTree_offline_exists_);
  genTree_->Branch("offline_idx", &genTree_offline_idx_);
  genTree_->Branch("offline_pt", &genTree_offline_pt_);
  genTree_->Branch("offline_eta", &genTree_offline_eta_);
  genTree_->Branch("offline_phi", &genTree_offline_phi_);
  genTree_->Branch("offline_dr", &genTree_offline_dr_);
  genTree_->Branch("online_exists", &genTree_online_exists_);
  genTree_->Branch("online_idx", &genTree_online_idx_);
  genTree_->Branch("online_pt", &genTree_online_pt_);
  genTree_->Branch("online_eta", &genTree_online_eta_);
  genTree_->Branch("online_phi", &genTree_online_phi_);
  genTree_->Branch("online_dr", &genTree_online_dr_);
  genTree_->Branch("online_hgcalIso", &genTree_online_hgcalIso_);
}

void PhotonComparisonNtuplizer::resetOfflineVectors() {
  Track_pt_.clear();
  Track_eta_.clear();
  Track_phi_.clear();
  Track_Time_.clear();
  Track_TimeErr_.clear();
  Track_MtdMva_.clear();
  Track_dz_.clear();
  Track_dxy_.clear();
  Trackster_pt_.clear();
  Trackster_eta_.clear();
  Trackster_phi_.clear();
  Trackster_dr_.clear();
  Trackster_Time_.clear();
  Trackster_TimeErr_.clear();
  Trackster_TrackIdx_.clear();
  Trackster_Trackpt_.clear();
  Trackster_Tracketa_.clear();
  Trackster_Trackphi_.clear();
  Trackster_isSeed_.clear();
  Trackster_inCluster_.clear();
  RecHit_dr_.clear();
  RecHit_eta_.clear();
  RecHit_phi_.clear();
  RecHit_layer_.clear();
  RecHit_energy_.clear();
  RecHit_time_.clear();
  RecHit_time_error_.clear();
  RecHit_isSeeded_.clear();
  RecHit_inCluster_.clear();
}

void PhotonComparisonNtuplizer::resetOnlineBranches() {
  online_genIdx_ = -1;
  online_objIdx_ = -1;
  online_candIdx_ = -1;
  online_matchedToGen_ = 0.;
  online_dr_gen_ = 999.;
  online_gen_pt_ = -1.;
  online_gen_eta_ = -999.;
  online_gen_phi_ = -999.;
  online_pt_ = -1.;
  online_eta_ = -999.;
  online_phi_ = -999.;
  online_energy_ = -1.;
  online_ecalPFIsol_ = -999.;
  online_hcalPFIsol_ = -999.;
  online_hgcalPFIsol_ = -999.;
  online_trkIsol_ = -999.;
  online_r9_ = -999.;
  online_sigmaIEtaIEta_ = -999.;
  online_rVar_ = -999.;
  online_hForHoverE_ = -999.;
}

void PhotonComparisonNtuplizer::resetOnlineVectors() {
  online_Trackster_pt_.clear();
  online_Trackster_eta_.clear();
  online_Trackster_phi_.clear();
  online_Trackster_dr_.clear();
  online_Trackster_Time_.clear();
  online_Trackster_TimeErr_.clear();
  online_Trackster_isSeed_.clear();
  online_Trackster_inCluster_.clear();
  online_LayerCluster_energy_.clear();
  online_LayerCluster_eta_.clear();
  online_LayerCluster_phi_.clear();
  online_LayerCluster_dr_.clear();
  online_LayerCluster_Time_.clear();
  online_LayerCluster_TimeErr_.clear();
  online_LayerCluster_isSeed_.clear();
  online_LayerCluster_inCluster_.clear();
  online_RecHit_dr_.clear();
  online_RecHit_eta_.clear();
  online_RecHit_phi_.clear();
  online_RecHit_layer_.clear();
  online_RecHit_energy_.clear();
  online_RecHit_time_.clear();
  online_RecHit_time_error_.clear();
  online_RecHit_isSeeded_.clear();
  online_RecHit_inCluster_.clear();
}

void PhotonComparisonNtuplizer::fillOnlineBranches(const reco::RecoEcalCandidate& candidate,
                                                   int candidateIndex,
                                                   const trigger::EgammaObject* object,
                                                   int objectIndex,
                                                   const egammaTiming::GenMatchResult& genMatch) {
  resetOnlineBranches();
  resetOnlineVectors();

  online_candIdx_ = candidateIndex;
  online_objIdx_ = objectIndex;
  online_pt_ = candidate.et();
  online_eta_ = candidate.eta();
  online_phi_ = candidate.phi();
  online_energy_ = candidate.energy();
  if (object != nullptr) {
    online_ecalPFIsol_ = egammaTiming::getEgammaVar(*object, "hltEgammaEcalPFClusterIsoUnseeded");
    online_hcalPFIsol_ = egammaTiming::getEgammaVar(*object, "hltEgammaHcalPFClusterIsoUnseeded");
    online_hgcalPFIsol_ = egammaTiming::getEgammaVar(*object, "hltEgammaHGCalLayerClusterIsoUnseeded");
    online_trkIsol_ = egammaTiming::getEgammaVar(*object, "hltEgammaHollowTrackIsoUnseeded");
    online_r9_ = egammaTiming::getEgammaVar(*object, "hltEgammaR9Unseeded");
    online_sigmaIEtaIEta_ = egammaTiming::getEgammaVar(*object, "hltEgammaClusterShapeUnseeded_sigmaIEtaIEta5x5");
    online_rVar_ = egammaTiming::getEgammaVar(*object, "hltEgammaHGCALIDVarsUnseeded_rVar");
    online_hForHoverE_ = egammaTiming::getEgammaVar(*object, "hltEgammaHGCALIDVarsUnseeded_hForHOverE");
  }
  online_genIdx_ = genMatch.index;
  online_matchedToGen_ = genMatch.truthClass;
  online_dr_gen_ = genMatch.dr;
  online_gen_pt_ = genMatch.pt;
  online_gen_eta_ = genMatch.eta;
  online_gen_phi_ = genMatch.phi;
}

void PhotonComparisonNtuplizer::fillGenBranches(int truthClass,
                                                int genIndex,
                                                const reco::GenParticle& genParticle,
                                                const std::pair<int, double>& offlineMatch,
                                                const std::pair<int, double>& onlineMatch,
                                                const std::vector<reco::Photon>& offlinePhotons,
                                                const std::vector<reco::RecoEcalCandidate>& onlineCandidates,
                                                const std::vector<trigger::EgammaObject>& onlineObjects) {
  genTree_idx_ = genIndex;
  genTree_truth_ = truthClass;
  genTree_pt_ = genParticle.pt();
  genTree_eta_ = genParticle.eta();
  genTree_phi_ = genParticle.phi();

  genTree_offline_idx_ = offlineMatch.first;
  genTree_offline_exists_ = offlineMatch.first >= 0 ? 1. : 0.;
  genTree_offline_pt_ = offlineMatch.first >= 0 ? offlinePhotons[offlineMatch.first].pt() : -1.;
  genTree_offline_eta_ = offlineMatch.first >= 0 ? offlinePhotons[offlineMatch.first].eta() : -999.;
  genTree_offline_phi_ = offlineMatch.first >= 0 ? offlinePhotons[offlineMatch.first].phi() : -999.;
  genTree_offline_dr_ = offlineMatch.second;

  genTree_online_idx_ = onlineMatch.first;
  genTree_online_exists_ = onlineMatch.first >= 0 ? 1. : 0.;
  genTree_online_pt_ = onlineMatch.first >= 0 ? onlineCandidates[onlineMatch.first].et() : -1.;
  genTree_online_eta_ = onlineMatch.first >= 0 ? onlineCandidates[onlineMatch.first].eta() : -999.;
  genTree_online_phi_ = onlineMatch.first >= 0 ? onlineCandidates[onlineMatch.first].phi() : -999.;
  genTree_online_dr_ = onlineMatch.second;
  genTree_online_hgcalIso_ = -999.;
  if (onlineMatch.first >= 0) {
    const auto objectMatch =
        egammaTiming::findBestRecoMatch(onlineCandidates[onlineMatch.first].eta(), onlineCandidates[onlineMatch.first].phi(), onlineObjects, gen_deltaR_);
    if (objectMatch.first >= 0) {
      genTree_online_hgcalIso_ =
          egammaTiming::getEgammaVar(onlineObjects[objectMatch.first], "hltEgammaHGCalLayerClusterIsoUnseeded");
    }
  }

  genTree_->Fill();
}

void PhotonComparisonNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  nEvent_ = iEvent.id().event();
  nRun_ = iEvent.id().run();
  nLumi_ = iEvent.id().luminosityBlock();

  auto const& caloGeometry = iSetup.getData(caloGeometryToken_);
  recHitTools_.setGeometry(caloGeometry);

  const auto photonHandle = iEvent.getHandle(photonProducer_);
  const auto onlineHandle = iEvent.getHandle(onlineProducer_);
  const auto onlineCandidateHandle = iEvent.getHandle(onlineCandidateProducer_);
  const auto vertexHandle = iEvent.getHandle(vertexProducer_);
  const auto trackHandle = iEvent.getHandle(trackProducer_);
  const auto genParticles = iEvent.getHandle(genParticleProducer_);
  const auto hitsEEHandle = iEvent.getHandle(hitsEE_);
  const auto hitsFHHandle = iEvent.getHandle(hitsFH_);
  const auto hitsBHHandle = iEvent.getHandle(hitsBH_);
  const auto onlineHitsEEHandle = iEvent.getHandle(onlineHitsEE_);
  const auto onlineHitsFHHandle = iEvent.getHandle(onlineHitsFH_);
  const auto onlineHitsBHHandle = iEvent.getHandle(onlineHitsBH_);
  const auto& beamPoint = iEvent.get(beamspotProducer_).position();
  const auto& mtdt0 = iEvent.get(mtdt0Token_);
  const auto& mtdSigmat0 = iEvent.get(mtdSigmat0Token_);
  const auto& mtdTrkQualMVA = iEvent.get(mtdTrkQualMVAToken_);

  edm::Handle<std::vector<ticl::Trackster>> tracksterHandle;
  iEvent.getByToken(tracksterToken_, tracksterHandle);
  const std::vector<ticl::Trackster> tracksters = tracksterHandle.isValid() ? *tracksterHandle : std::vector<ticl::Trackster>();

  edm::Handle<std::vector<reco::CaloCluster>> layerClusterHandle;
  iEvent.getByToken(layerClusterToken_, layerClusterHandle);
  const std::vector<reco::CaloCluster>* layerClusters = layerClusterHandle.isValid() ? layerClusterHandle.product() : nullptr;

  std::vector<ticl::Trackster> onlineTracksters;
  {
    edm::Handle<std::vector<ticl::Trackster>> handle;
    iEvent.getByToken(onlineTracksterToken_, handle);
    if (handle.isValid()) {
      onlineTracksters = *handle;
    }
  }

  edm::Handle<std::vector<reco::CaloCluster>> onlineLayerClusterHandle;
  iEvent.getByToken(onlineLayerClusterToken_, onlineLayerClusterHandle);
  const std::vector<reco::CaloCluster>* onlineLayerClusters =
      onlineLayerClusterHandle.isValid() ? onlineLayerClusterHandle.product() : nullptr;

  edm::Handle<edm::ValueMap<std::pair<float, float>>> onlineLayerClusterTimeHandle;
  iEvent.getByToken(onlineLayerClusterTimeToken_, onlineLayerClusterTimeHandle);

  const std::vector<reco::Vertex> vertices = *vertexHandle;
  reco::Vertex primaryVertex = egammaTiming::getPrimaryVertex(vertices, vtxN_);
  PVTime_ = primaryVertex.t();
  PVTimeErr_ = primaryVertex.tError();
  sigTrkIdx_ = -1;
  sigTrkTime_ = -1.;
  sigTrkTimeErr_ = -1.;
  sigTrkMtdMva_ = -1.;

  std::vector<int> selectedGenIndices;
  std::vector<int> selectedTruthClasses;
  if (isMC_) {
    for (unsigned int index = 0; index < genParticles->size(); ++index) {
      const auto& particle = (*genParticles)[index];
      if (std::abs(particle.pdgId()) != 22 || particle.status() != 1) {
        continue;
      }
      selectedGenIndices.push_back(static_cast<int>(index));
      selectedTruthClasses.push_back(egammaTiming::classifyPhotonTruth(particle));
    }
  }

  const auto& offlinePhotons = *photonHandle;
  const auto& onlineObjects = *onlineHandle;
  const auto& onlineCandidates = *onlineCandidateHandle;

  for (unsigned int index = 0; index < selectedGenIndices.size(); ++index) {
    const int genIndex = selectedGenIndices[index];
    const auto& particle = (*genParticles)[genIndex];
    std::pair<int, double> offlineMatch{-1, 999.};
    for (unsigned int i = 0; i < offlinePhotons.size(); ++i) {
      if (offlinePhotons[i].pt() < offlinePtMin_) {
        continue;
      }
      const double dr = reco::deltaR(particle.eta(), particle.phi(), offlinePhotons[i].eta(), offlinePhotons[i].phi());
      if (dr < offlineMatch.second && dr < gen_deltaR_) {
        offlineMatch = {static_cast<int>(i), dr};
      }
    }
    std::pair<int, double> onlineMatch{-1, 999.};
    for (unsigned int i = 0; i < onlineCandidates.size(); ++i) {
      if (onlineCandidates[i].et() < onlinePtMin_) {
        continue;
      }
      const double dr = reco::deltaR(particle.eta(), particle.phi(), onlineCandidates[i].eta(), onlineCandidates[i].phi());
      if (dr < onlineMatch.second && dr < gen_deltaR_) {
        onlineMatch = {static_cast<int>(i), dr};
      }
    }
    fillGenBranches(selectedTruthClasses[index], genIndex, particle, offlineMatch, onlineMatch, offlinePhotons, onlineCandidates, onlineObjects);
  }

  for (unsigned int i = 0; i < offlinePhotons.size(); ++i) {
    const reco::Photon* photon = &(offlinePhotons.at(i));
    if (photon->pt() < offlinePtMin_) {
      continue;
    }
    const auto genMatch = isMC_ ? egammaTiming::findBestGenMatch(photon->eta(),
                                                                 photon->phi(),
                                                                 *genParticles,
                                                                 22,
                                                                 gen_deltaR_,
                                                                 egammaTiming::classifyPhotonTruth)
                                : egammaTiming::GenMatchResult();

    Pho_pt_ = photon->pt();
    Pho_eta_ = photon->eta();
    Pho_phi_ = photon->phi();
    matchedToGenPho_ = genMatch.truthClass;
    genIdx_ = genMatch.index;
    gen_pt_ = genMatch.pt;
    gen_eta_ = genMatch.eta;
    gen_phi_ = genMatch.phi;
    dr_gen_offline_ = genMatch.dr;

    const reco::CaloClusterPtr& photonSeed = photon->superCluster()->seed();
    const reco::CaloClusterPtrVector& photonClusters = photon->superCluster()->clusters();
    const auto& seedHitsAndFractions = photonSeed->hitsAndFractions();

    resetOfflineVectors();

    int trackIndex = -1;
    for (reco::TrackCollection::const_iterator itrTr = trackHandle->begin(); itrTr != trackHandle->end(); ++itrTr) {
      trackIndex++;
      const reco::TrackRef trackRef(trackHandle, trackIndex);
      const double thisPt = itrTr->pt();
      if (thisPt < ptMin_) {
        continue;
      }

      const double dr = ROOT::Math::VectorUtil::DeltaR(itrTr->momentum(), photon->momentum());
      const double deta = itrTr->eta() - Pho_eta_;
      const bool isBarrel = std::abs(Pho_eta_) < 1.479;
      const double intRadius = isBarrel ? intRadiusBarrel_ : intRadiusEndcap_;
      const double strip = isBarrel ? stripBarrel_ : stripEndcap_;
      if (!(dr < extRadius_ && dr >= intRadius && std::abs(deta) >= strip)) {
        continue;
      }

      double dz = 0.;
      switch (dzOption_) {
        case egammaisolation::EgammaTrackSelector::dz:
          dz = std::fabs(itrTr->dz() - photon->vertex().z());
          break;
        case egammaisolation::EgammaTrackSelector::vz:
          dz = std::fabs(itrTr->vz() - photon->vertex().z());
          break;
        case egammaisolation::EgammaTrackSelector::bs:
          dz = std::fabs(itrTr->dz(beamPoint) - photon->vertex().z());
          break;
        case egammaisolation::EgammaTrackSelector::vtx:
          dz = std::fabs(itrTr->dz(photon->vertex()));
          break;
        default:
          dz = std::fabs(itrTr->vz() - photon->vertex().z());
          break;
      }

      Track_pt_.push_back(thisPt);
      Track_eta_.push_back(itrTr->eta());
      Track_phi_.push_back(itrTr->phi());
      Track_Time_.push_back(mtdt0[trackRef]);
      Track_TimeErr_.push_back(mtdSigmat0[trackRef]);
      Track_MtdMva_.push_back(mtdTrkQualMVA[trackRef]);
      Track_dz_.push_back(dz);
      Track_dxy_.push_back(std::fabs(itrTr->dxy(beamPoint)));
    }

    if (layerClusters != nullptr) {
      for (const auto& trackster : tracksters) {
        if (trackster.vertices().empty()) {
          continue;
        }
        const auto momentum = trackster.eigenvectors(0);
        const double dr = ROOT::Math::VectorUtil::DeltaR(momentum, photon->momentum());
        if (dr >= extRadius_) {
          continue;
        }

        Trackster_pt_.push_back(trackster.raw_pt());
        Trackster_eta_.push_back(momentum.eta());
        Trackster_phi_.push_back(momentum.phi());
        Trackster_dr_.push_back(dr);
        Trackster_Time_.push_back(trackster.time());
        Trackster_TimeErr_.push_back(trackster.timeError());
        Trackster_TrackIdx_.push_back(trackster.trackIdx());
        if (trackster.trackIdx() >= 0) {
          const reco::TrackRef trkRef(trackHandle, trackster.trackIdx());
          Trackster_Trackpt_.push_back(trkRef->pt());
          Trackster_Tracketa_.push_back(trkRef->eta());
          Trackster_Trackphi_.push_back(trkRef->phi());
        } else {
          Trackster_Trackpt_.push_back(-1.);
          Trackster_Tracketa_.push_back(-999.);
          Trackster_Trackphi_.push_back(-999.);
        }

        bool isSeedTrackster = false;
        bool tracksterInCluster = false;
        for (const auto layerClusterId : trackster.vertices()) {
          if (isSeedTrackster && tracksterInCluster) {
            break;
          }
          const auto& hitAndFractions = (*layerClusters)[layerClusterId].hitsAndFractions();
          for (const auto& seedHit : seedHitsAndFractions) {
            for (const auto& hit : hitAndFractions) {
              if (hit.first == seedHit.first) {
                isSeedTrackster = true;
                break;
              }
            }
            if (isSeedTrackster) {
              break;
            }
          }

          for (auto clusterIt = photonClusters.begin(); clusterIt != photonClusters.end() && !tracksterInCluster; ++clusterIt) {
            const auto& clusterHits = (*clusterIt)->hitsAndFractions();
            for (const auto& clusterHit : clusterHits) {
              for (const auto& hit : hitAndFractions) {
                if (clusterHit.first == hit.first) {
                  tracksterInCluster = true;
                  break;
                }
              }
              if (tracksterInCluster) {
                break;
              }
            }
          }
        }
        Trackster_isSeed_.push_back(isSeedTrackster);
        Trackster_inCluster_.push_back(tracksterInCluster);
      }
    }

    auto fillRecHit = [this, &seedHitsAndFractions, &photonClusters](const HGCRecHit& hit) {
      const GlobalPoint position = recHitTools_.getPosition(hit.id());
      const float eta = recHitTools_.getEta(position, 0);
      const float phi = recHitTools_.getPhi(position);
      const float deltaR = reco::deltaR(Pho_eta_, Pho_phi_, eta, phi);
      if (deltaR >= extRadius_) {
        return;
      }

      bool isSeeded = false;
      bool inCluster = false;
      for (const auto& seedHit : seedHitsAndFractions) {
        if (seedHit.first == hit.id()) {
          isSeeded = true;
          break;
        }
      }

      for (auto clusterIt = photonClusters.begin(); clusterIt != photonClusters.end() && !inCluster; ++clusterIt) {
        const auto& clusterHits = (*clusterIt)->hitsAndFractions();
        for (const auto& clusterHit : clusterHits) {
          if (clusterHit.first == hit.id()) {
            inCluster = true;
            break;
          }
        }
      }

      RecHit_dr_.push_back(deltaR);
      RecHit_eta_.push_back(eta);
      RecHit_phi_.push_back(phi);
      RecHit_layer_.push_back(recHitTools_.getLayerWithOffset(hit.id()));
      RecHit_energy_.push_back(hit.energy());
      RecHit_time_.push_back(hit.time());
      RecHit_time_error_.push_back(hit.timeError());
      RecHit_isSeeded_.push_back(isSeeded);
      RecHit_inCluster_.push_back(inCluster);
    };

    for (const auto& hit : *hitsEEHandle) {
      fillRecHit(hit);
    }
    for (const auto& hit : *hitsFHHandle) {
      fillRecHit(hit);
    }
    for (const auto& hit : *hitsBHHandle) {
      fillRecHit(hit);
    }

    offlineTree_->Fill();
  }

  for (unsigned int i = 0; i < onlineCandidates.size(); ++i) {
    if (onlineCandidates[i].et() < onlinePtMin_) {
      continue;
    }
    const auto genMatch = isMC_ ? egammaTiming::findBestGenMatch(onlineCandidates[i].eta(),
                                                                 onlineCandidates[i].phi(),
                                                                 *genParticles,
                                                                 22,
                                                                 gen_deltaR_,
                                                                 egammaTiming::classifyPhotonTruth)
                                : egammaTiming::GenMatchResult();
    const auto onlineObjectMatch =
        egammaTiming::findBestRecoMatch(onlineCandidates[i].eta(), onlineCandidates[i].phi(), onlineObjects, gen_deltaR_);
    const trigger::EgammaObject* matchedObject = onlineObjectMatch.first >= 0 ? &onlineObjects[onlineObjectMatch.first] : nullptr;
    fillOnlineBranches(onlineCandidates[i], i, matchedObject, onlineObjectMatch.first, genMatch);

    if (onlineLayerClusters != nullptr) {
      const auto& candidate = onlineCandidates[i];
      const reco::SuperClusterRef scRef = candidate.superCluster();

      if (scRef.isNonnull() && std::abs(candidate.eta()) >= 1.479) {
        const reco::SuperCluster& hgcalSC = *scRef;
        const auto& seedHitsAndFractions = hgcalSC.seed()->hitsAndFractions();
        const auto& egClusters = hgcalSC.clusters();

        std::vector<bool> layerClusterIsSeed(onlineLayerClusters->size(), false);
        std::vector<bool> layerClusterInCluster(onlineLayerClusters->size(), false);

        for (size_t clusterIndex = 0; clusterIndex < onlineLayerClusters->size(); ++clusterIndex) {
          const auto& cluster = (*onlineLayerClusters)[clusterIndex];
          layerClusterIsSeed[clusterIndex] = egammaTiming::haveCommonHit(cluster.hitsAndFractions(), seedHitsAndFractions);

          for (const auto& egCluster : egClusters) {
            if (egammaTiming::haveCommonHit(cluster.hitsAndFractions(), egCluster->hitsAndFractions())) {
              layerClusterInCluster[clusterIndex] = true;
              break;
            }
          }

          const double clusterDr = reco::deltaR(online_eta_, online_phi_, cluster.eta(), cluster.phi());
          if (clusterDr >= extRadius_) {
            continue;
          }

          online_LayerCluster_dr_.push_back(clusterDr);
          online_LayerCluster_energy_.push_back(cluster.energy());
          online_LayerCluster_eta_.push_back(cluster.eta());
          online_LayerCluster_phi_.push_back(cluster.phi());
          if (onlineLayerClusterTimeHandle.isValid()) {
            edm::Ref<reco::CaloClusterCollection> clusterRef(onlineLayerClusterHandle, clusterIndex);
            const auto timePair = (*onlineLayerClusterTimeHandle)[clusterRef];
            online_LayerCluster_Time_.push_back(timePair.first);
            online_LayerCluster_TimeErr_.push_back(timePair.second);
          } else {
            online_LayerCluster_Time_.push_back(-99.);
            online_LayerCluster_TimeErr_.push_back(-1.);
          }
          online_LayerCluster_isSeed_.push_back(layerClusterIsSeed[clusterIndex]);
          online_LayerCluster_inCluster_.push_back(layerClusterInCluster[clusterIndex]);
        }

        for (const auto& trackster : onlineTracksters) {
          if (trackster.vertices().empty()) {
            continue;
          }
          const auto momentum = trackster.eigenvectors(0);
          const double tracksterDr = reco::deltaR(online_eta_, online_phi_, momentum.eta(), momentum.phi());
          if (tracksterDr >= extRadius_) {
            continue;
          }

          bool isSeedTrackster = false;
          bool tracksterInCluster = false;
          for (const auto layerClusterId : trackster.vertices()) {
            if (layerClusterId >= layerClusterIsSeed.size()) {
              continue;
            }
            isSeedTrackster = isSeedTrackster || layerClusterIsSeed[layerClusterId];
            tracksterInCluster = tracksterInCluster || layerClusterInCluster[layerClusterId];
          }

          online_Trackster_dr_.push_back(tracksterDr);
          online_Trackster_pt_.push_back(trackster.raw_pt());
          online_Trackster_eta_.push_back(momentum.eta());
          online_Trackster_phi_.push_back(momentum.phi());
          online_Trackster_Time_.push_back(trackster.time());
          online_Trackster_TimeErr_.push_back(trackster.timeError());
          online_Trackster_isSeed_.push_back(isSeedTrackster);
          online_Trackster_inCluster_.push_back(tracksterInCluster);
        }

        auto fillOnlineRecHit = [this, &seedHitsAndFractions, &egClusters](const HGCRecHit& hit) {
          const GlobalPoint position = recHitTools_.getPosition(hit.id());
          const float eta = recHitTools_.getEta(position, 0);
          const float phi = recHitTools_.getPhi(position);
          const float deltaR = reco::deltaR(online_eta_, online_phi_, eta, phi);
          if (deltaR >= extRadius_) {
            return;
          }

          bool isSeeded = false;
          bool inCluster = false;
          for (const auto& seedHit : seedHitsAndFractions) {
            if (seedHit.first == hit.id()) {
              isSeeded = true;
              break;
            }
          }

          for (auto clusterIt = egClusters.begin(); clusterIt != egClusters.end() && !inCluster; ++clusterIt) {
            for (const auto& clusterHit : (*clusterIt)->hitsAndFractions()) {
              if (clusterHit.first == hit.id()) {
                inCluster = true;
                break;
              }
            }
          }

          online_RecHit_dr_.push_back(deltaR);
          online_RecHit_eta_.push_back(eta);
          online_RecHit_phi_.push_back(phi);
          online_RecHit_layer_.push_back(recHitTools_.getLayerWithOffset(hit.id()));
          online_RecHit_energy_.push_back(hit.energy());
          online_RecHit_time_.push_back(hit.time());
          online_RecHit_time_error_.push_back(hit.timeError());
          online_RecHit_isSeeded_.push_back(isSeeded);
          online_RecHit_inCluster_.push_back(inCluster);
        };

        if (onlineHitsEEHandle.isValid()) {
          for (const auto& hit : *onlineHitsEEHandle) {
            fillOnlineRecHit(hit);
          }
        }
        if (onlineHitsFHHandle.isValid()) {
          for (const auto& hit : *onlineHitsFHHandle) {
            fillOnlineRecHit(hit);
          }
        }
        if (onlineHitsBHHandle.isValid()) {
          for (const auto& hit : *onlineHitsBHHandle) {
            fillOnlineRecHit(hit);
          }
        }
      }
    }
    onlineTree_->Fill();
  }
}

void PhotonComparisonNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("photonProducer", edm::InputTag("photonsHGC"));
  desc.add<edm::InputTag>("onlineProducer", edm::InputTag("hltEgammaHLTExtra", "Unseeded", "HLT"));
  desc.add<edm::InputTag>("onlineCandidateProducer", edm::InputTag("hltEgammaCandidatesUnseeded", "", "HLT"));
  desc.add<edm::InputTag>("onlineTracksterSrc", edm::InputTag("hltTiclCandidate", "", "HLT"));
  desc.add<edm::InputTag>("onlineLayerClusterSrc", edm::InputTag("hltMergeLayerClusters", "", "HLT"));
  desc.add<edm::InputTag>("onlineTimeLayerClusterSrc", edm::InputTag("hltMergeLayerClusters", "timeLayerCluster", "HLT"));
  desc.add<edm::InputTag>("onlineRecHitsEE_Src", edm::InputTag("hltHGCalRecHit", "HGCEERecHits", "HLT"));
  desc.add<edm::InputTag>("onlineRecHitsFH_Src", edm::InputTag("hltHGCalRecHit", "HGCHEFRecHits", "HLT"));
  desc.add<edm::InputTag>("onlineRecHitsBH_Src", edm::InputTag("hltHGCalRecHit", "HGCHEBRecHits", "HLT"));
  desc.add<edm::InputTag>("trackProducer", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("BeamspotProducer", edm::InputTag("offlineBeamSpot"));
  desc.add<edm::InputTag>("vertexProducer", edm::InputTag("offlineSlimmedPrimaryVertices4D"));
  desc.add<edm::InputTag>("genParticles", edm::InputTag("prunedGenParticles"));
  desc.add<edm::InputTag>("tracksterSrc", edm::InputTag("ticlTrackstersMerge"));
  desc.add<edm::InputTag>("LayerClusterSrc", edm::InputTag("hgcalMergeLayerClusters"));
  desc.add<edm::InputTag>("RecHitsEE_Src", edm::InputTag("HGCalRecHit", "HGCEERecHits"));
  desc.add<edm::InputTag>("RecHitsFH_Src", edm::InputTag("HGCalRecHit", "HGCHEFRecHits"));
  desc.add<edm::InputTag>("RecHitsBH_Src", edm::InputTag("HGCalRecHit", "HGCHEBRecHits"));
  desc.add<edm::InputTag>("mtdt0", edm::InputTag("tofPID:t0"));
  desc.add<edm::InputTag>("mtdSigmat0", edm::InputTag("tofPID:sigmat0"));
  desc.add<edm::InputTag>("mtdTrkQualMVA", edm::InputTag("mtdTrackQualityMVA:mtdQualMVA"));
  desc.add<double>("ptMin", 1.0);
  desc.add<double>("offlinePtMin", 15.0);
  desc.add<double>("onlinePtMin", 15.0);
  desc.add<double>("intRadiusBarrel", 0.01);
  desc.add<double>("intRadiusEndcap", 0.01);
  desc.add<double>("stripBarrel", 0.01);
  desc.add<double>("stripEndcap", 0.01);
  desc.add<double>("extRadius", 0.3);
  desc.add<double>("gen_deltaR", 0.1);
  desc.add<bool>("isMC", true);
  descriptions.add("photonComparisonNtuplizer", desc);
}

DEFINE_FWK_MODULE(PhotonComparisonNtuplizer);
