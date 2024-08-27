#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaTrackSelector.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/GeometryVector/interface/PV3DBase.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"

#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"


#include <Math/VectorUtil.h>
#include <TTree.h>
#include <TFile.h>

enum PhotonMatchType{
  FAKE_PHOTON,
  TRUE_PROMPT_PHOTON,
  TRUE_NON_PROMPT_PHOTON,
};


class TrackerNtuplizer_Photon : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit TrackerNtuplizer_Photon(const edm::ParameterSet&);
  int matchToTruth(reco::Photon const& photon, edm::View<reco::GenParticle> const& genParticles, double deltaR) const;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  void setDzOption(const std::string& s) {
      if (!s.compare("dz")) dzOption_ = egammaisolation::EgammaTrackSelector::dz;
      else if (!s.compare("vz")) dzOption_ = egammaisolation::EgammaTrackSelector::vz;
      else if (!s.compare("bs")) dzOption_ = egammaisolation::EgammaTrackSelector::bs;
      else if (!s.compare("vtx")) dzOption_ = egammaisolation::EgammaTrackSelector::vtx;
      else       dzOption_ = egammaisolation::EgammaTrackSelector::dz;
  }
private:
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void beginJob() override{};
  void endJob()   override{};

  // Global variables
  int nEvent_;
  int nRun_;
  int nLumi_;
  int getNpu_;
  int vtxN_;
  int sigTrkIdx_;
  double sigTrkTime_;
  double sigTrkTimeErr_;
  double sigTrkMtdMva_;

  double Pho_pt_;
  double Pho_eta_;
  double Pho_phi_;

  double PVTime_;
  double PVTimeErr_;

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
  std::vector<int>    Trackster_TrackIdx_;
  std::vector<double> Trackster_Trackpt_;
  std::vector<double> Trackster_Tracketa_;
  std::vector<double> Trackster_Trackphi_;
  std::vector<bool>   Trackster_isSeed_;

  std::vector<double> RecHit_dr_;
  std::vector<double> RecHit_eta_;
  std::vector<double> RecHit_phi_;
  std::vector<int>    RecHit_layer_;
  std::vector<double> RecHit_energy_;
  std::vector<double> RecHit_time_;
  std::vector<double> RecHit_time_error_;
  std::vector<bool>   RecHit_isSeeded_;
  // Tokens
  const edm::EDGetTokenT<std::vector<reco::Photon>> photonProducer_;
  const edm::EDGetTokenT<reco::TrackCollection> trackProducer_;
  const edm::EDGetTokenT<reco::BeamSpot> beamspotProducer_;
  const edm::EDGetTokenT<std::vector<reco::Vertex>> vertexProducer_;
  const edm::EDGetTokenT<edm::View<reco::GenParticle>> genParticleProducer_;

  edm::EDGetTokenT<std::vector<ticl::Trackster>> tracksterToken_;
  edm::Handle<std::vector<ticl::Trackster>>  trackstersH_;

  edm::EDGetTokenT<std::vector<reco::CaloCluster>> layer_cluster_Token_;
  edm::Handle<std::vector<reco::CaloCluster>>      layer_cluster_Handle_;


  const edm::EDGetTokenT<HGCRecHitCollection> hitsEE_;
  const edm::EDGetTokenT<HGCRecHitCollection> hitsFH_;
  const edm::EDGetTokenT<HGCRecHitCollection> hitsBH_;

  const edm::EDGetTokenT<edm::ValueMap<float>> mtdt0_H;
  const edm::EDGetTokenT<edm::ValueMap<float>> mtdSigmat0_H;
  const edm::EDGetTokenT<edm::ValueMap<float>> mtdTrkQualMVA_H;

  const edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometryToken_;

  // Cut Value
  const double ptMin_;
  const double intRadiusBarrel_;
  const double intRadiusEndcap_;
  const double stripBarrel_;
  const double stripEndcap_;
  const double extRadius_;
  const double gen_deltaR_;

  int dzOption_;

  const int dtRef_;
  const int dtType_;

// Other
  bool   isMC_;
  double matchedToGenPho_;
  hgcal::RecHitTools recHitTools_;

  TTree* tree_;
};

TrackerNtuplizer_Photon::TrackerNtuplizer_Photon(const edm::ParameterSet& config)
  : photonProducer_{consumes(config.getParameter<edm::InputTag>("photonProducer"))},
    trackProducer_{consumes(config.getParameter<edm::InputTag>("trackProducer"))},
    beamspotProducer_{consumes(config.getParameter<edm::InputTag>("BeamspotProducer"))},
    vertexProducer_{consumes(config.getParameter<edm::InputTag>("vertexProducer"))},
    genParticleProducer_{consumes<edm::View<reco::GenParticle>>(config.getParameter<edm::InputTag>("genParticles"))},
    tracksterToken_{consumes(config.getParameter<edm::InputTag>("tracksterSrc"))},
    layer_cluster_Token_{consumes(config.getParameter<edm::InputTag>("LayerClusterSrc"))},
    hitsEE_{consumes(config.getParameter<edm::InputTag>("RecHitsEE_Src"))},
    hitsFH_{consumes(config.getParameter<edm::InputTag>("RecHitsFH_Src"))},
    hitsBH_{consumes(config.getParameter<edm::InputTag>("RecHitsBH_Src"))},
    mtdt0_H{consumes(config.getParameter<edm::InputTag>("mtdt0"))},
    mtdSigmat0_H{consumes(config.getParameter<edm::InputTag>("mtdSigmat0"))},
    mtdTrkQualMVA_H{consumes(config.getParameter<edm::InputTag>("mtdTrkQualMVA"))},    
    caloGeometryToken_(esConsumes()),
    ptMin_{config.getParameter<double>("ptMin")},
    intRadiusBarrel_{config.getParameter<double>("intRadiusBarrel")},
    intRadiusEndcap_{config.getParameter<double>("intRadiusEndcap")},
    stripBarrel_{config.getParameter<double>("stripBarrel")},
    stripEndcap_{config.getParameter<double>("stripEndcap")},
    extRadius_{config.getParameter<double>("extRadius")},
    gen_deltaR_(config.getParameter<double>("gen_deltaR")),
    dtRef_{config.getParameter<int>("dtRef")},
    dtType_{config.getParameter<int>("dtType")},
    isMC_(config.getParameter<bool>("isMC"))
{
  setDzOption("vz");
  // Book Tree
  usesResource(TFileService::kSharedResource);
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree", "tree");
  tree_->Branch("nEvent", &nEvent_);
  tree_->Branch("nRun",   &nRun_);
  tree_->Branch("nLumi",  &nLumi_);
  tree_->Branch("vtxN",   &vtxN_);
  tree_->Branch("Pho_pt", &Pho_pt_);
  tree_->Branch("Pho_eta", &Pho_eta_);
  tree_->Branch("Pho_phi", &Pho_phi_);
  tree_->Branch("matchedToGenPho", &matchedToGenPho_);
  tree_->Branch("PV_Time", &PVTime_);
  tree_->Branch("PV_TimeErr", &PVTimeErr_);
  tree_->Branch("sigTrkIdx",  &sigTrkIdx_);
  tree_->Branch("sigTrkTime",  &sigTrkTime_);
  tree_->Branch("sigTrkTimeErr", &sigTrkTimeErr_);
  tree_->Branch("sigTrkMtdMva", &sigTrkMtdMva_);
  tree_->Branch("Track_pt", &Track_pt_);
  tree_->Branch("Track_eta", &Track_eta_);
  tree_->Branch("Track_phi", &Track_phi_);
  tree_->Branch("Track_Time", &Track_Time_);
  tree_->Branch("Track_TimeErr", &Track_TimeErr_);
  tree_->Branch("Track_MtdMva", &Track_MtdMva_);
  tree_->Branch("Track_dz", &Track_dz_);
  tree_->Branch("Track_dxy", &Track_dxy_);
  tree_->Branch("Trackster_dr", &Trackster_dr_);
  tree_->Branch("Trackster_isSeed", &Trackster_isSeed_);
  tree_->Branch("Trackster_pt", &Trackster_pt_);
  tree_->Branch("Trackster_eta", &Trackster_eta_);
  tree_->Branch("Trackster_phi", &Trackster_phi_);
  tree_->Branch("Trackster_Time", &Trackster_Time_);
  tree_->Branch("Trackster_TimeErr", &Trackster_TimeErr_);
  tree_->Branch("Trackster_TrackIdx", &Trackster_TrackIdx_);
  tree_->Branch("Trackster_Trackpt",  &Trackster_Trackpt_);
  tree_->Branch("Trackster_Tracketa", &Trackster_Tracketa_);
  tree_->Branch("Trackster_Trackphi", &Trackster_Trackphi_);

  tree_->Branch("HGCRecHit_dr",  &RecHit_dr_);
  tree_->Branch("HGCRecHit_eta", &RecHit_eta_);
  tree_->Branch("HGCRecHit_phi", &RecHit_phi_);
  tree_->Branch("HGCRecHit_layer", &RecHit_layer_);
  tree_->Branch("HGCRecHit_energy", &RecHit_energy_);
  tree_->Branch("HGCRecHit_Time",   &RecHit_time_);
  tree_->Branch("HGCRecHit_TimeError", &RecHit_time_error_);
  tree_->Branch("HGCRecHit_isSeed", &RecHit_isSeeded_);
}

void TrackerNtuplizer_Photon::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  nEvent_ = iEvent.id().event();
  nRun_   = iEvent.id().run();
  nLumi_  = iEvent.id().luminosityBlock();
  auto const& caloGeometry = iSetup.getData(caloGeometryToken_);
  recHitTools_.setGeometry(caloGeometry);

  auto photonHandle = iEvent.getHandle(photonProducer_);
  auto vertexHandle   = iEvent.getHandle(vertexProducer_);
  const edm::Handle <reco::TrackCollection> &trackHandle    = iEvent.getHandle(trackProducer_);
  auto beamPoint_      = iEvent.get(beamspotProducer_).position();
  auto genParticles    = iEvent.getHandle(genParticleProducer_);


  const edm::ValueMap<float> mtdt0_ = (const edm::ValueMap<float> &) (iEvent.get(mtdt0_H));
  const edm::ValueMap<float> mtdSigmat0_ = (const edm::ValueMap<float> &) (iEvent.get(mtdSigmat0_H));
  const edm::ValueMap<float> mtdTrkQualMVA_ = (const edm::ValueMap<float> &) (iEvent.get(mtdTrkQualMVA_H));

  const reco::TrackCollection *trackCollection_ = (trackHandle.product());
  const edm::Handle <reco::TrackCollection> trackCollectionH_ = trackHandle;

  const edm::Handle<HGCRecHitCollection> hitsEE_Handle_ = iEvent.getHandle(hitsEE_);
  const edm::Handle<HGCRecHitCollection> hitsFH_Handle_ = iEvent.getHandle(hitsFH_);
  const edm::Handle<HGCRecHitCollection> hitsBH_Handle_ = iEvent.getHandle(hitsBH_);


  std::vector <reco::Vertex> vertices = *vertexHandle;
  reco::Vertex pmVtx;
  int nGoodVertex = 0;
  for(int iVtx = 0; iVtx < (int) vertices.size(); iVtx++)
  {
     const reco::Vertex &vertex = vertices.at(iVtx);
     bool isGoodVertex = (
                !vertex.isFake() &&
                vertex.ndof() >= 4 );
     nGoodVertex += (int) isGoodVertex;
      if(nGoodVertex)
      {
         pmVtx = vertex;
         break;
      }
  }

  std::vector <reco::Photon> photons = *photonHandle;
  std::cout << "nPhoton: " << photons.size() << std::endl;
  for (unsigned int i = 0; i < photons.size(); ++i){

    // Get Gsf Electron track
    const reco::Photon* photon = &(photons.at(i));

   
    
    const reco::CaloClusterPtr& photon_seed = photon->superCluster()->seed();
    const std::vector<std::pair<DetId, float>>& seedHitsAndFractions = photon_seed->hitsAndFractions();

    std::cout << "------------------" << std::endl;


    Pho_pt_ = photon->pt();
    Pho_eta_ = photon->eta();
    Pho_phi_ = photon->phi();

    // GenMatch
    if (isMC_){
      matchedToGenPho_ = matchToTruth(*photon, *genParticles, 0.3);
    }

    int sigTrkIdx_ = -1;
    // Get signal track time information

    reco::TrackRef sigTrkRef;
    double sigTrkTime = -1;
    double sigTrkTimeErr = -1;
    double sigTrkMtdMva  = -1;

    if(sigTrkIdx_ > -1){
      sigTrkRef = reco::TrackRef(trackCollectionH_, sigTrkIdx_);
      sigTrkTime = mtdt0_[sigTrkRef];
      sigTrkTimeErr = mtdSigmat0_[sigTrkRef];
      sigTrkMtdMva  = mtdTrkQualMVA_[sigTrkRef];
    }
   
    sigTrkTime_     = sigTrkTime;
    sigTrkTimeErr_  = sigTrkTimeErr;
    sigTrkMtdMva_   = sigTrkMtdMva;



    PVTime_    = pmVtx.tError();
    PVTimeErr_ = pmVtx.t();

    // Loop tracks
    Track_pt_.clear();
    Track_eta_.clear();
    Track_phi_.clear();
    Track_Time_.clear();
    Track_TimeErr_.clear();
    Track_MtdMva_.clear();
    Track_dz_.clear();
    Track_dxy_.clear();

    int TrkIdx = -1;
    for(reco::TrackCollection::const_iterator itrTr = (*trackCollection_).begin(); itrTr != (*trackCollection_).end();++itrTr){
      TrkIdx ++;
      const reco::TrackRef trkRef(trackCollectionH_, TrkIdx);
      double this_pt = (*itrTr).pt();
      if(this_pt < ptMin_) continue; //Minimum pT cut

      // Isolation cone cut
      double dr = ROOT::Math::VectorUtil::DeltaR(itrTr->momentum(), photon->momentum());
      double deta = (*itrTr).eta() - Pho_eta_;
      bool isBarrel = std::abs(Pho_eta_) < 1.479;
      double intRadius = isBarrel ? intRadiusBarrel_ : intRadiusEndcap_;
      double strip = isBarrel ? stripBarrel_ : stripEndcap_;
      if (!(dr < extRadius_ && dr >= intRadius && std::abs(deta) >= strip)) continue; 

      double dz = 0.;
      switch(dzOption_){
        case egammaisolation::EgammaTrackSelector::dz:
          dz = fabs((*itrTr).dz() - photon->vertex().z());
          break;
        case egammaisolation::EgammaTrackSelector::vz:
          dz = fabs((*itrTr).vz() - photon->vertex().z());
          break;
        case egammaisolation::EgammaTrackSelector::bs:
          dz = fabs((*itrTr).dz(beamPoint_) - photon->vertex().z());
          break;
        case egammaisolation::EgammaTrackSelector::vtx:
          dz = fabs((*itrTr).dz(photon->vertex()));
          break;
        default:
          dz = fabs((*itrTr).vz() - photon->vertex().z());
          break;
      }

      double dxy = fabs((*itrTr).dxy(beamPoint_));

      double trkTime = mtdt0_[trkRef];
      double trkTimeErr = mtdSigmat0_[trkRef];
      double trkMtdMva = mtdTrkQualMVA_[trkRef];
      double eta     = (*itrTr).eta();
      double phi     = (*itrTr).phi();
      Track_pt_.push_back(this_pt);
      Track_eta_.push_back(eta);
      Track_phi_.push_back(phi);
      Track_Time_.push_back(trkTime);
      Track_TimeErr_.push_back(trkTimeErr);
      Track_MtdMva_.push_back(trkMtdMva);
      Track_dz_.push_back(dz);
      Track_dxy_.push_back(dxy);
    }

    iEvent.getByToken(tracksterToken_, trackstersH_);
    const auto& tracksters = *trackstersH_;
   
    iEvent.getByToken(layer_cluster_Token_, layer_cluster_Handle_);
    const auto& layerClusters = *layer_cluster_Handle_.product();

    // Find signal tracksters

    int sigTracksterIdx = -1;
    const ticl::Trackster* sigTracksterPtr = 0;
    Trackster_pt_.clear();
    Trackster_eta_.clear();
    Trackster_phi_.clear();
    Trackster_dr_.clear();
    Trackster_isSeed_.clear();
    Trackster_Time_.clear();
    Trackster_TimeErr_.clear();
    Trackster_TrackIdx_.clear();
    Trackster_Trackpt_.clear();
    Trackster_Tracketa_.clear();
    Trackster_Trackphi_.clear();

    int trackster_Index = -1;
    for (const auto& tst: tracksters){
      trackster_Index ++;
      if (tst.vertices().empty()) continue;
      auto momentum = tst.eigenvectors(0);
      double dr = ROOT::Math::VectorUtil::DeltaR(momentum, photon->momentum());
      bool isBarrel = std::abs(Pho_eta_) < 1.479;
      double intRadius = isBarrel ? intRadiusBarrel_ : intRadiusEndcap_;
      if (!(dr < extRadius_)) continue;
      Trackster_pt_.push_back(tst.raw_pt());
      Trackster_eta_.push_back(momentum.eta());
      Trackster_phi_.push_back(momentum.phi());
      Trackster_dr_.push_back(dr);
      Trackster_Time_.push_back(tst.time());
      Trackster_TimeErr_.push_back(tst.timeError());
      Trackster_TrackIdx_.push_back(tst.trackIdx());
      if(tst.trackIdx() != -1){
        const reco::TrackRef trkRef(trackCollectionH_, tst.trackIdx());
        const reco::Track& track = *trkRef;
        Trackster_Trackpt_.push_back(track.pt());
        Trackster_Tracketa_.push_back(track.eta());
        Trackster_Trackphi_.push_back(track.phi());
      }
      else{
        Trackster_Trackpt_.push_back(-1);
        Trackster_Tracketa_.push_back(-1);
        Trackster_Trackphi_.push_back(-1);
      }

      // Check if trackster is seeded or not.
      bool isSeedTrackster = false;
      for (const auto lcId : tst.vertices()){
        if (isSeedTrackster) break;
        const std::vector<std::pair<DetId, float>>& hit_and_fractions = layerClusters[lcId].hitsAndFractions();

        for(unsigned int seedhitId = 0; seedhitId < seedHitsAndFractions.size(); seedhitId++){
          if (isSeedTrackster) break;
          const auto seed_detId = seedHitsAndFractions[seedhitId].first;
          const auto seed_fraction = seedHitsAndFractions[seedhitId].second;
          for(unsigned int hitId = 0; hitId < hit_and_fractions.size(); hitId ++){
            if (isSeedTrackster) break;
            const auto hit_detId = hit_and_fractions[hitId].first;
            const auto hit_fraction = hit_and_fractions[hitId].second;
            if(hit_detId == seed_detId) isSeedTrackster = true;
          }
        }
      }
      Trackster_isSeed_.push_back(isSeedTrackster);
      if(isSeedTrackster) std::cout << "Find a seed trackster " << trackster_Index << std::endl;

    }

    auto checkAndFill = [this, &seedHitsAndFractions](const HGCRecHit& hit){

      const GlobalPoint position = recHitTools_.getPosition(hit.id());
      float eta = recHitTools_.getEta(position, 0);
      float phi = recHitTools_.getPhi(position);
      float deltar = reco::deltaR(Pho_eta_, Pho_phi_, eta, phi);
      bool isBarrel = std::abs(Pho_eta_) < 1.479;
      double intRadius = isBarrel ? intRadiusBarrel_ : intRadiusEndcap_;
      if (deltar < extRadius_){

        bool isSeeded = false;
        for(unsigned int seedhitId = 0; seedhitId < seedHitsAndFractions.size(); seedhitId++){
           if (isSeeded) break;
           if(seedHitsAndFractions[seedhitId].first == hit.id()) isSeeded = true;
        }
        int layer = recHitTools_.getLayerWithOffset(hit.id());
      


        RecHit_dr_.push_back(deltar);
        RecHit_eta_.push_back(eta);
        RecHit_phi_.push_back(phi);
        RecHit_layer_.push_back(layer);
        RecHit_energy_.push_back(hit.energy());
        RecHit_time_.push_back(hit.time());
        RecHit_time_error_.push_back(hit.timeError());
        RecHit_isSeeded_.push_back(isSeeded);
      }
    };


    RecHit_dr_.clear();
    RecHit_eta_.clear();
    RecHit_phi_.clear();
    RecHit_layer_.clear();
    RecHit_energy_.clear();
    RecHit_time_.clear();
    RecHit_time_error_.clear();
    RecHit_isSeeded_.clear();
    for (const auto& hit : *hitsEE_Handle_){
      checkAndFill(hit);
    }
    for (const auto& hit : *hitsFH_Handle_){
      checkAndFill(hit);
    }
    for (const auto& hit : *hitsBH_Handle_){
      checkAndFill(hit);
    }

    tree_->Fill();
  }
}

void TrackerNtuplizer_Photon::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("trackProducer", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("mtdt0", edm::InputTag("tofPID:t0"));
  desc.add<edm::InputTag>("mtdSigmat0", edm::InputTag("tofPID:sigmat0"));
  desc.add<edm::InputTag>("mtdTrkQualMVA", edm::InputTag("mtdTrackQualityMVA:mtdQualMVA"));
  desc.add<double>("gen_deltaR", 0.1);
}

int TrackerNtuplizer_Photon::matchToTruth(reco::Photon const &ph, edm::View<reco::GenParticle> const& genParticles, double deltaR) const{

  double dR = 999;
  reco::GenParticle const* closestPhoton = &genParticles[0];
  for (auto const& particle : genParticles){
    if (std::abs(particle.pdgId()) != 22 || particle.status() != 1) continue;
    double dRtmp = ROOT::Math::VectorUtil::DeltaR(ph.p4(), particle.p4());
    if (dRtmp < dR){
      dR = dRtmp;
      closestPhoton = &particle;
    }
  }

  if (dR < deltaR){
    if (closestPhoton->isPromptFinalState())
      return TRUE_PROMPT_PHOTON;
    else
      return TRUE_NON_PROMPT_PHOTON;
  }
  return FAKE_PHOTON;

}

DEFINE_FWK_MODULE(TrackerNtuplizer_Photon);

