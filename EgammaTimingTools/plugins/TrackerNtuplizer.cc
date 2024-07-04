#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaTrackSelector.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/GeometryVector/interface/PV3DBase.h"

#include <Math/VectorUtil.h>
#include <TTree.h>
#include <TFile.h>

enum ElectronMatchType{
  UNMATCHED,
  TRUE_PROMPT_ELECTRON,
  TRUE_ELECTRON_FROM_TAU,
  TRUE_NON_PROMPT_ELECTRON,
};


class TrackerNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit TrackerNtuplizer(const edm::ParameterSet&);
  int matchToTruth(reco::GsfElectron const& electron, edm::View<reco::GenParticle> const& genParticles) const;
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
  double sigTrkTime_;
  double sigTrkTimeErr_;
  double sigTrkMtdMva_;

  double Ele_pt_;
  double Ele_eta_;
  double Ele_phi_;

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

  // Tokens
  const edm::EDGetTokenT<reco::GsfElectronCollection> electronProducer_;
  const edm::EDGetTokenT<reco::TrackCollection> trackProducer_;
  const edm::EDGetTokenT<reco::BeamSpot> beamspotProducer_;
  const edm::EDGetTokenT<std::vector<reco::Vertex>> vertexProducer_;
  const edm::EDGetTokenT<edm::View<reco::GenParticle>> genParticleProducer_;

  edm::EDGetTokenT<std::vector<ticl::Trackster>> tracksterToken_;
  edm::Handle<std::vector<ticl::Trackster>>  trackstersH_;

  const edm::EDGetTokenT<edm::ValueMap<float>> mtdt0_H;
  const edm::EDGetTokenT<edm::ValueMap<float>> mtdSigmat0_H;
  const edm::EDGetTokenT<edm::ValueMap<float>> mtdTrkQualMVA_H;


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
  double matchedToGenEle_;
  TTree* tree_;
};

TrackerNtuplizer::TrackerNtuplizer(const edm::ParameterSet& config)
  : electronProducer_{consumes(config.getParameter<edm::InputTag>("electronProducer"))},
    trackProducer_{consumes(config.getParameter<edm::InputTag>("trackProducer"))},
    beamspotProducer_{consumes(config.getParameter<edm::InputTag>("BeamspotProducer"))},
    vertexProducer_{consumes(config.getParameter<edm::InputTag>("vertexProducer"))},
    genParticleProducer_{consumes<edm::View<reco::GenParticle>>(config.getParameter<edm::InputTag>("genParticles"))},
    tracksterToken_{consumes(config.getParameter<edm::InputTag>("tracksterSrc"))},
    mtdt0_H{consumes(config.getParameter<edm::InputTag>("mtdt0"))},
    mtdSigmat0_H{consumes(config.getParameter<edm::InputTag>("mtdSigmat0"))},
    mtdTrkQualMVA_H{consumes(config.getParameter<edm::InputTag>("mtdTrkQualMVA"))},    
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
  tree_->Branch("Ele_pt", &Ele_pt_);
  tree_->Branch("Ele_eta", &Ele_eta_);
  tree_->Branch("Ele_phi", &Ele_phi_);
  tree_->Branch("matchedToGenEle", &matchedToGenEle_);
  tree_->Branch("PV_Time", &PVTime_);
  tree_->Branch("PV_TimeErr", &PVTimeErr_);
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
  tree_->Branch("Trackster_pt", &Trackster_pt_);
  tree_->Branch("Trackster_eta", &Trackster_eta_);
  tree_->Branch("Trackster_phi", &Trackster_phi_);
  tree_->Branch("Trackster_Time", &Trackster_Time_);
  tree_->Branch("Trackster_TimeErr", &Trackster_TimeErr_);
  tree_->Branch("Trackster_TrackIdx", &Trackster_TrackIdx_);
}

void TrackerNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  nEvent_ = iEvent.id().event();
  nRun_   = iEvent.id().run();
  nLumi_  = iEvent.id().luminosityBlock();

  auto electronHandle = iEvent.getHandle(electronProducer_);
  auto vertexHandle   = iEvent.getHandle(vertexProducer_);
  const edm::Handle <reco::TrackCollection> &trackHandle    = iEvent.getHandle(trackProducer_);
  auto beamPoint_      = iEvent.get(beamspotProducer_).position();
  auto genParticles    = iEvent.getHandle(genParticleProducer_);


  const edm::ValueMap<float> mtdt0_ = (const edm::ValueMap<float> &) (iEvent.get(mtdt0_H));
  const edm::ValueMap<float> mtdSigmat0_ = (const edm::ValueMap<float> &) (iEvent.get(mtdSigmat0_H));
  const edm::ValueMap<float> mtdTrkQualMVA_ = (const edm::ValueMap<float> &) (iEvent.get(mtdTrkQualMVA_H));

  const reco::TrackCollection *trackCollection_ = (trackHandle.product());
  const edm::Handle <reco::TrackCollection> trackCollectionH_ = trackHandle;


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



  for (unsigned int i = 0; i < electronHandle->size(); ++i){

    // Get Gsf Electron track
    const reco::GsfElectron* electron = &(electronHandle->at(i));
    const reco::Track* tmpTrack = &(*(electron->gsfTrack()));
    math::XYZVector tmpElectronMomentumAtVtx = (*tmpTrack).momentum();
    double tmpElectronEtaAtVertex = (*tmpTrack).eta();

    int sigTrkIdx = -1;
    const reco::Track* sigTrkPtr = 0;
    
    Ele_pt_ = electron->pt();
    Ele_eta_ = electron->eta();
    Ele_phi_ = electron->phi();

    // GenMatch
    if (isMC_){
      matchedToGenEle_ = matchToTruth(*electron, *genParticles);
    }

    // Find signal Track
    int trkIdx = -1;
    for(reco::TrackCollection::const_iterator itrTr = (*trackCollection_).begin(); itrTr != (*trackCollection_).end(); ++itrTr) {
        trkIdx += 1;
        double dr = ROOT::Math::VectorUtil::DeltaR(itrTr->momentum(), tmpElectronMomentumAtVtx);
        bool isBarrel = std::abs(tmpElectronEtaAtVertex) < 1.479;
        double intRadius = isBarrel ? intRadiusBarrel_ : intRadiusEndcap_;
        if (dr < intRadius){
          if (!sigTrkPtr || (*itrTr).pt() > (*sigTrkPtr).pt()) {
            sigTrkPtr = &(*itrTr);
            sigTrkIdx = trkIdx;
          }
        }
    }

    // Get signal track time information

    reco::TrackRef sigTrkRef;
    double sigTrkTime = -1;
    double sigTrkTimeErr = -1;
    double sigTrkMtdMva  = -1;

    if(sigTrkIdx > -1){
      sigTrkRef = reco::TrackRef(trackCollectionH_, sigTrkIdx);
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
      double dr = ROOT::Math::VectorUtil::DeltaR(itrTr->momentum(), tmpElectronMomentumAtVtx);
      double deta = (*itrTr).eta() - tmpElectronEtaAtVertex;
      bool isBarrel = std::abs(tmpElectronEtaAtVertex) < 1.479;
      double intRadius = isBarrel ? intRadiusBarrel_ : intRadiusEndcap_;
      double strip = isBarrel ? stripBarrel_ : stripEndcap_;
      if (!(dr < extRadius_ && dr >= intRadius && std::abs(deta) >= strip)) continue; 

      double dz = 0.;
      switch(dzOption_){
        case egammaisolation::EgammaTrackSelector::dz:
          dz = fabs((*itrTr).dz() - (*tmpTrack).dz());
          break;
        case egammaisolation::EgammaTrackSelector::vz:
          dz = fabs((*itrTr).vz() - (*tmpTrack).vz());
          break;
        case egammaisolation::EgammaTrackSelector::bs:
          dz = fabs((*itrTr).dz(beamPoint_) - (*tmpTrack).dz(beamPoint_));
          break;
        case egammaisolation::EgammaTrackSelector::vtx:
          dz = fabs((*itrTr).dz(tmpTrack->vertex()));
          break;
        default:
          dz = fabs((*itrTr).vz() - (*tmpTrack).vz());
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

    // Find signal tracksters

    int sigTracksterIdx = -1;
    const ticl::Trackster* sigTracksterPtr = 0;
    Trackster_pt_.clear();
    Trackster_eta_.clear();
    Trackster_phi_.clear();
    Trackster_dr_.clear();
    Trackster_Time_.clear();
    Trackster_TimeErr_.clear();
    Trackster_TrackIdx_.clear();

    for (const auto& tst: tracksters){
      if (tst.vertices().empty()) continue;
      auto momentum = tst.eigenvectors(0);
      double dr = ROOT::Math::VectorUtil::DeltaR(momentum, tmpElectronMomentumAtVtx);
      bool isBarrel = std::abs(tmpElectronEtaAtVertex) < 1.479;
      double intRadius = isBarrel ? intRadiusBarrel_ : intRadiusEndcap_;
      if (!(dr < extRadius_)) continue;
      Trackster_pt_.push_back(tst.raw_pt());
      Trackster_eta_.push_back(momentum.eta());
      Trackster_phi_.push_back(momentum.phi());
      Trackster_dr_.push_back(dr);
      Trackster_Time_.push_back(tst.time());
      Trackster_TimeErr_.push_back(tst.timeError());
      Trackster_TrackIdx_.push_back(tst.trackIdx());
    }

    tree_->Fill();
  }
}

void TrackerNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("trackProducer", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("mtdt0", edm::InputTag("tofPID:t0"));
  desc.add<edm::InputTag>("mtdSigmat0", edm::InputTag("tofPID:sigmat0"));
  desc.add<edm::InputTag>("mtdTrkQualMVA", edm::InputTag("mtdTrackQualityMVA:mtdQualMVA"));
  desc.add<double>("gen_deltaR", 0.1);
}

int TrackerNtuplizer::matchToTruth(reco::GsfElectron const &electron, edm::View<reco::GenParticle> const& genParticles) const{

  double dR = 999;
  reco::GenParticle const* closestElectron = nullptr;
  for (auto const& particle : genParticles){
    if (std::abs(particle.pdgId()) != 11 || particle.status() != 1) continue;
    double dRtmp = ROOT::Math::VectorUtil::DeltaR(electron.p4(), particle.p4());
    if (dRtmp < dR){
      dR = dRtmp;
      closestElectron = &particle;
    }
  }

  if (closestElectron == nullptr || dR >= gen_deltaR_)
    return UNMATCHED;

  if (closestElectron->fromHardProcessFinalState())
    return TRUE_PROMPT_ELECTRON;

  if (closestElectron->isDirectHardProcessTauDecayProductFinalState())
    return TRUE_ELECTRON_FROM_TAU;

  return TRUE_NON_PROMPT_ELECTRON;

}
DEFINE_FWK_MODULE(TrackerNtuplizer);
