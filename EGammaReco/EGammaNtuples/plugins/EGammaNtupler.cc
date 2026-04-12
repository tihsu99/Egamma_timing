#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/isFinite.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HLTReco/interface/EgammaObject.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateIsolation.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"
#include "DataFormats/Common/interface/getRef.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"


#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"

#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "TTree.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"

#include "RecoEgamma/EgammaTools/interface/EgammaHGCALIDParamDefaults.h"
#include "RecoEgamma/EgammaTools/interface/HGCalShowerShapeHelper.h"
#include <Math/VectorUtil.h>
#include <unordered_map>


enum ElectronMatchType{
  UNMATCHED,
  TRUE_PROMPT_ELECTRON,
  TRUE_ELECTRON_FROM_TAU,
  TRUE_NON_PROMPT_ELECTRON,
};


using namespace reco;
namespace {

  typedef std::vector<HGCRecHit> HGCRecHitCollectionOld;

  //bool is if a valid dr was found, float is the dr
  std::pair<bool, float> getMaxDRNonSeedCluster(const reco::SuperCluster& sc) {
    float maxDR2 = 0.;
    const edm::Ptr<reco::CaloCluster>& seedClus = sc.seed();

    for (const auto& clus : sc.clusters()) {
      if (clus == seedClus) {
        continue;
      }

      // find cluster with max dR
      const double dr2 = reco::deltaR2(*clus, *seedClus);
      if (dr2 > maxDR2) {
        maxDR2 = dr2;
      }
    }
    return {sc.clustersSize() != 1, sc.clustersSize() != 1 ? std::sqrt(maxDR2) : 999.};
  }
  template <typename T>
  int countRecHits(const T& recHitHandle, float threshold) {
    int count = 0;
    if (recHitHandle.isValid()) {
      for (const auto& recHit : *recHitHandle) {
        if (recHit.energy() > threshold) {
          count++;
        }
      }
    }
    return count;
  }
}  // namespace


class EGammaNtuples : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit EGammaNtuples(const edm::ParameterSet&);
  ~EGammaNtuples() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
//  bool isEE(const reco::SuperCluster&);
//  bool isEB(const reco::SuperCluster&);
//  float cal_r9(const reco::SuperCluster&, const EcalRecHitCollection&, const CaloTopology&, const bool);
  void fillTree(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  int matchToTruth(const reco::RecoEcalCandidate &egamma,  const std::vector<reco::GenParticle>& genParticles) const;
  double getLayerClusterTime( const reco::CaloCluster &lc, const std::unordered_map<DetId, HGCRecHit> &recHitMap);

  // Utility: check overlap between two hit collections
  inline bool haveCommonHit(
    const std::vector<std::pair<DetId, float>>& hits1,
    const std::vector<std::pair<DetId, float>>& hits2)
  {
    for (const auto& h1 : hits1) {
        for (const auto& h2 : hits2) {
            if (h1.first == h2.first) return true;
        }
    }
    return false;
  }
//  std::vector<reco::GenParticle> getGenParticles(const std::vector<reco::GenParticle>& genParticles, std::string objType_);
//  reco::GenParticle getLastCopyPreFSR(reco::GenParticle part);
//  const reco::GenParticle* get_best_dr_match(const reco::SuperCluster& obj_to_match, const std::vector<reco::GenParticle>& genparts, double max_dr);

  // ----- member data -----
  edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticleToken_;
  edm::EDGetTokenT<std::vector<reco::SuperCluster>> scBarrelL1SeededToken_;
  edm::EDGetTokenT<std::vector<reco::SuperCluster>> scHGCalL1SeededToken_;
  edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>> ebRecHitsToken_;
  edm::EDGetTokenT<HGCRecHitCollection> HGCeeRecHitsToken_;
  edm::EDGetTokenT<HGCRecHitCollection> HGChebRecHitsToken_;
  edm::EDGetTokenT<HGCRecHitCollection> HGChefRecHitsToken_;
  edm::EDGetTokenT<reco::RecoEcalCandidateIsolationMap> sigmaIEtaIEtaToken_;
  edm::EDGetTokenT<reco::RecoEcalCandidateIsolationMap> sigmaIPhiIPhiToken_;
  edm::EDGetTokenT<reco::RecoEcalCandidateIsolationMap> sigmaIEtaIEtaNoiseCleanedToken_;
  edm::EDGetTokenT<reco::RecoEcalCandidateIsolationMap> sigmaIPhiIPhiNoiseCleanedToken_;
  edm::EDGetTokenT<int> nrHitsEB1GeVToken_;
  edm::EDGetTokenT<int> nrHitsEE1GeVToken_;
  edm::EDGetTokenT<std::vector<reco::RecoEcalCandidate>> eGammaCandidatesToken_;
  edm::EDGetTokenT<reco::VertexCollection> verticesToken_;
  edm::EDGetTokenT<reco::PFRecHitCollection> pfRecHitsHGCALToken_;
  edm::EDGetTokenT<std::vector<ticl::Trackster>> tracksterToken_;
  edm::EDGetTokenT<reco::CaloClusterCollection> layerClusterToken_;
  edm::EDGetTokenT<edm::ValueMap<std::pair<float, float>>> timelayerClusterToken_;

  std::string pType_;
  double ptThreshold_;
  double gen_deltaR_;
  double max_deltaR_;

  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;
  edm::ESGetToken<CaloTopology, CaloTopologyRecord> caloTopoToken_;

  const CaloTopology* caloTopology_;
  const CaloGeometry* caloGeometry_;
  std::map<DetId, const HGCRecHit*> hitMap;


  HGCalShowerShapeHelper hgcalShowerShapes_;
  float hgcalCylinderR_ = EgammaHGCALIDParamDefaults::kRCylinder;
  hgcal::RecHitTools recHitTools_;

  edm::Handle<reco::VertexCollection> vertices_;



  TFile *newfile = new TFile("output.root", "RECREATE", "", 207);
  TTree *tree_barrel   = new TTree("eg_barrel", "eg_barrel");
  TTree *tree_endcap   = new TTree("eg_endcap", "eg_endcap");



  float eg_energy;
  float eg_et;
  float eg_eta;
  float eg_phi;
  int   eg_status;

  float eg_sigmaIEtaIEta, eg_sigmaIPhiIPhi, eg_sigmaIEtaIEtaNoiseCleaned, eg_sigmaIPhiIPhiNoiseCleaned;
  float eg_sigma2uu, eg_sigma2vv, eg_sigma2ww;

  //HGCRecHit
  std::vector<double> HGCRecHit_energy_;
  std::vector<double> HGCRecHit_eta_;
  std::vector<double> HGCRecHit_phi_;
  std::vector<double> HGCRecHit_dr_;
  std::vector<double> HGCRecHit_time_;
  std::vector<double> HGCRecHit_timeErr_;
  std::vector<int>    HGCRecHit_layer_;
  std::vector<bool>   HGCRecHit_isSeed_;
  std::vector<bool>   HGCRecHit_inCluster_;


  //Trackster 
  std::vector<double> Trackster_pt_;
  std::vector<double> Trackster_eta_;
  std::vector<double> Trackster_phi_;
  std::vector<double> Trackster_dr_;
  std::vector<double> Trackster_time_;
  std::vector<double> Trackster_timeErr_;
  std::vector<bool>   Trackster_isSeed_;
  std::vector<bool>   Trackster_inCluster_;

  // Layer Cluster
  std::vector<double> LayerCluster_energy_;
  std::vector<double> LayerCluster_eta_;
  std::vector<double> LayerCluster_phi_;
  std::vector<double> LayerCluster_dr_;
  std::vector<double> LayerCluster_time_;
  std::vector<bool>   LayerCluster_isSeed_;
  std::vector<bool>   LayerCluster_inCluster_;

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif
};

EGammaNtuples::EGammaNtuples(const edm::ParameterSet& iConfig)
  : genParticleToken_(consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticles"))),
  scBarrelL1SeededToken_(consumes<std::vector<reco::SuperCluster>>(iConfig.getParameter<edm::InputTag>("scBarrelL1Seeded"))),
  scHGCalL1SeededToken_(consumes<std::vector<reco::SuperCluster>>(iConfig.getParameter<edm::InputTag>("scHGCalL1Seeded"))),
  ebRecHitsToken_(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(iConfig.getParameter<edm::InputTag>("ebRecHits"))),
  HGCeeRecHitsToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCeeRecHits"))),
  HGChebRecHitsToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGChebRecHits"))),
  HGChefRecHitsToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGChefRecHits"))),
  sigmaIEtaIEtaToken_(consumes<reco::RecoEcalCandidateIsolationMap> (iConfig.getParameter<edm::InputTag>("sigmaIEtaIEta"))),
  sigmaIPhiIPhiToken_(consumes<reco::RecoEcalCandidateIsolationMap> (iConfig.getParameter<edm::InputTag>("sigmaIPhiIPhi"))),
  sigmaIEtaIEtaNoiseCleanedToken_(consumes<reco::RecoEcalCandidateIsolationMap> (iConfig.getParameter<edm::InputTag>("sigmaIEtaIEtaNoiseCleaned"))),
  sigmaIPhiIPhiNoiseCleanedToken_(consumes<reco::RecoEcalCandidateIsolationMap> (iConfig.getParameter<edm::InputTag>("sigmaIPhiIPhiNoiseCleaned"))),
  nrHitsEB1GeVToken_(consumes<int>(iConfig.getParameter<edm::InputTag>("nrHitsEB1GeV"))),
  nrHitsEE1GeVToken_(consumes<int>(iConfig.getParameter<edm::InputTag>("nrHitsEE1GeV"))),
  eGammaCandidatesToken_(consumes<std::vector<reco::RecoEcalCandidate>>(iConfig.getParameter<edm::InputTag>("eGammaCandidates"))),
  verticesToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  pfRecHitsHGCALToken_(consumes<reco::PFRecHitCollection>(iConfig.getParameter<edm::InputTag>("pfHGCALRecHits"))),
  tracksterToken_(consumes<std::vector<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("trackster"))),
  layerClusterToken_(consumes<reco::CaloClusterCollection>(iConfig.getParameter<edm::InputTag>("layerClusters"))),
  timelayerClusterToken_(consumes<edm::ValueMap<std::pair<float, float>>>(iConfig.getParameter<edm::InputTag>("timelayerClusters"))),
  pType_(iConfig.getParameter<std::string>("pType")),
  ptThreshold_(iConfig.getParameter<double>("ptThreshold")),
  gen_deltaR_(iConfig.getParameter<double>("genMatchDR")),
  max_deltaR_(iConfig.getParameter<double>("maxDR")),
  caloGeomToken_(esConsumes<CaloGeometry, CaloGeometryRecord>()),
  caloTopoToken_(esConsumes<CaloTopology, CaloTopologyRecord>()){
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
    hgcalShowerShapes_.setTokens<edm::Transition::Event>(consumesCollector());

  }
  

EGammaNtuples::~EGammaNtuples() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}


void EGammaNtuples::fillTree(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace edm;

  edm::Handle<EcalRecHitCollection> ebRecHitsHandle;  
  iEvent.getByToken(ebRecHitsToken_, ebRecHitsHandle);

  edm::Handle<HGCRecHitCollection> HGCeeRecHitsHandle;
  edm::Handle<HGCRecHitCollection> HGChefRecHitsHandle;
  edm::Handle<HGCRecHitCollection> HGChebRecHitsHandle;
  edm::Handle<HGCRecHitCollection> HFNoseRecHitsHandle;

  iEvent.getByToken(HGCeeRecHitsToken_, HGCeeRecHitsHandle);
  iEvent.getByToken(HGChefRecHitsToken_, HGChefRecHitsHandle);
  iEvent.getByToken(HGChebRecHitsToken_, HGChebRecHitsHandle);

  std::vector<HGCRecHit> HGCeeRecHitsVec(HGCeeRecHitsHandle->begin(), HGCeeRecHitsHandle->end());
  std::vector<HGCRecHit> HGChefRecHitsVec(HGChefRecHitsHandle->begin(), HGChefRecHitsHandle->end());
  std::vector<HGCRecHit> HGChebRecHitsVec(HGChebRecHitsHandle->begin(), HGChebRecHitsHandle->end());

  std::vector<HGCRecHit> AllHGCRecHitsVec;
  AllHGCRecHitsVec.reserve(
    HGCeeRecHitsVec.size() + 
    HGChefRecHitsVec.size() + 
    HGChebRecHitsVec.size()  
  );

  AllHGCRecHitsVec.insert(AllHGCRecHitsVec.end(), HGCeeRecHitsVec.begin(), HGCeeRecHitsVec.end());
  AllHGCRecHitsVec.insert(AllHGCRecHitsVec.end(), HGChefRecHitsVec.begin(), HGChefRecHitsVec.end());
  AllHGCRecHitsVec.insert(AllHGCRecHitsVec.end(), HGChebRecHitsVec.begin(), HGChebRecHitsVec.end());


  edm::Handle<reco::RecoEcalCandidateIsolationMap> sigmaIEtaIEtaHandle;
  iEvent.getByToken(sigmaIEtaIEtaToken_,sigmaIEtaIEtaHandle);

  edm::Handle<reco::RecoEcalCandidateIsolationMap> sigmaIPhiIPhiHandle;
  iEvent.getByToken(sigmaIPhiIPhiToken_,sigmaIPhiIPhiHandle);

  edm::Handle<reco::RecoEcalCandidateIsolationMap> sigmaIEtaIEtaNoiseCleanedHandle;
  iEvent.getByToken(sigmaIEtaIEtaNoiseCleanedToken_,sigmaIEtaIEtaNoiseCleanedHandle);

  edm::Handle<reco::RecoEcalCandidateIsolationMap> sigmaIPhiIPhiNoiseCleanedHandle;
  iEvent.getByToken(sigmaIPhiIPhiNoiseCleanedToken_,sigmaIPhiIPhiNoiseCleanedHandle);

  edm::Handle<std::vector<reco::RecoEcalCandidate>> eGammaCandidatesHandle;
  iEvent.getByToken(eGammaCandidatesToken_,eGammaCandidatesHandle);

  edm::Handle<std::vector<ticl::Trackster>> tracksterHandle;
  iEvent.getByToken(tracksterToken_, tracksterHandle);
  const auto& tracksters = *tracksterHandle;

  edm::Handle<reco::CaloClusterCollection> clusterHandle;
  iEvent.getByToken(layerClusterToken_, clusterHandle);
  const std::vector<reco::CaloCluster> layerClusters = *(clusterHandle.product());


  edm::Handle<edm::ValueMap<std::pair<float, float>>> lcTimeHandle;
  iEvent.getByToken(timelayerClusterToken_, lcTimeHandle);



  // Objects
  const auto& gPs       = iEvent.get(genParticleToken_);


//  hgcalShowerShapes_.initPerEvent(iSetup, eeRecHitsVec);
  for (int i=0; i<static_cast<int>((*eGammaCandidatesHandle).size()); i++){
    const auto &cand = eGammaCandidatesHandle->at(i);
    reco::RecoEcalCandidateRef candidateRef = getRef(eGammaCandidatesHandle, i);
    reco::SuperClusterRef scRef = cand.superCluster();

    if (cand.pt() < ptThreshold_) continue;

    // EGamma Candidate kinematics
    eg_et = cand.et();
    eg_eta = cand.eta();
    eg_energy = cand.energy();
    eg_phi = cand.phi();
    eg_status =  matchToTruth(cand, gPs);


    if (std::abs(cand.eta()) < 1.479){ 
      eg_sigmaIEtaIEta = (*sigmaIEtaIEtaHandle)[candidateRef];
      eg_sigmaIPhiIPhi = (*sigmaIPhiIPhiHandle)[candidateRef];
      eg_sigmaIEtaIEtaNoiseCleaned = (*sigmaIEtaIEtaNoiseCleanedHandle)[candidateRef];
      eg_sigmaIPhiIPhiNoiseCleaned = (*sigmaIPhiIPhiNoiseCleanedHandle)[candidateRef];
      tree_barrel->Fill();
    }

    else{


      const reco::SuperCluster& hgcal_sc = *scRef;
      const std::vector<std::pair<DetId, float>>& seedHitsAndFractions = hgcal_sc.seed()->hitsAndFractions();
      const reco::CaloClusterPtrVector& eg_clusters = hgcal_sc.clusters();

      if (!(hgcal_sc.seed()->hitsAndFractions()[0].first.subdetId() == 0)) continue;

//      auto ssCalc = hgcalShowerShapes_.createCalc(hgcal_sc);
//      auto pcaWidths = ssCalc.getPCAWidths(hgcalCylinderR_);
//      auto energyHighestHits = ssCalc.getEnergyHighestHits(2);

 //     eg_sigma2uu = std::sqrt(pcaWidths.sigma2uu);
 //     eg_sigma2vv = std::sqrt(pcaWidths.sigma2vv);
 //     eg_sigma2ww = std::sqrt(pcaWidths.sigma2ww);


      HGCRecHit_dr_.clear();
      HGCRecHit_energy_.clear();
      HGCRecHit_eta_.clear();
      HGCRecHit_phi_.clear();
      HGCRecHit_time_.clear();
      HGCRecHit_timeErr_.clear();
      HGCRecHit_layer_.clear();
      HGCRecHit_isSeed_.clear();
      HGCRecHit_inCluster_.clear();


      std::unordered_map<DetId, HGCRecHit> recHitMap;
      for (auto const& rechit : AllHGCRecHitsVec) {
         recHitMap[rechit.detid()] = rechit;
         const GlobalPoint position = recHitTools_.getPosition(rechit.id());
         float eta = recHitTools_.getEta(position, 0);
         float phi = recHitTools_.getPhi(position);
         float deltar = reco::deltaR(eg_eta, eg_phi, eta, phi);
         if (deltar > max_deltaR_) continue;
         bool isSeeded = false;
         bool inCluster = false;
         int  layer = recHitTools_.getLayerWithOffset(rechit.id());
         for (const auto& h2 : seedHitsAndFractions) {
             if (rechit.id() == h2.first){
               isSeeded=true;
               break;
             }
         }
         for (const auto& bc: eg_clusters){
           for(const auto& h3 : bc->hitsAndFractions()){
             if (rechit.id() == h3.first){
               inCluster = true;
               break;
             }
           }
           if (inCluster) break;
         }

         HGCRecHit_dr_.push_back(deltar);
         HGCRecHit_energy_.push_back(rechit.energy());
         HGCRecHit_eta_.push_back(eta);
         HGCRecHit_phi_.push_back(phi);
         HGCRecHit_time_.push_back(rechit.time());
         HGCRecHit_timeErr_.push_back(rechit.timeError());
         HGCRecHit_layer_.push_back(layer);
         HGCRecHit_isSeed_.push_back(isSeeded);
         HGCRecHit_inCluster_.push_back(inCluster);
      }







      LayerCluster_dr_.clear();
      LayerCluster_energy_.clear();
      LayerCluster_eta_.clear();
      LayerCluster_phi_.clear();
      LayerCluster_time_.clear();
      LayerCluster_isSeed_.clear();
      LayerCluster_inCluster_.clear();

      std::vector<bool> LayerCluster_isSeed;
      std::vector<bool> LayerCluster_inCluster;
 
      for (size_t iclus = 0; iclus < clusterHandle->size(); ++iclus) {
        const auto& clus = (*clusterHandle)[iclus];
        edm::Ref<reco::CaloClusterCollection> clusRef(clusterHandle, iclus);

        bool isSeed = haveCommonHit(clus.hitsAndFractions(), seedHitsAndFractions);
        bool inCluster = false;
        for (const auto& bc : eg_clusters){
          if (haveCommonHit(clus.hitsAndFractions(), bc->hitsAndFractions())){
            inCluster = true;
            break;
          }
        }

        LayerCluster_isSeed.push_back(isSeed);
        LayerCluster_inCluster.push_back(inCluster);
        float dR = reco::deltaR(eg_eta, eg_phi, clus.eta(), clus.phi());
        if (dR > max_deltaR_) continue;
        LayerCluster_dr_.push_back(dR);
        LayerCluster_energy_.push_back(clus.energy());
        LayerCluster_eta_.push_back(clus.eta());
        LayerCluster_phi_.push_back(clus.phi());
//        LayerCluster_time_.push_back(getLayerClusterTime(clus, recHitMap));
        auto timePair = (*lcTimeHandle)[clusRef];  // pair<float, float>: {time, timeError}
        LayerCluster_time_.push_back(timePair.first);

        LayerCluster_isSeed_.push_back(isSeed);
        LayerCluster_inCluster_.push_back(inCluster);
      }

      Trackster_dr_.clear();
      Trackster_pt_.clear();
      Trackster_eta_.clear();
      Trackster_phi_.clear();
      Trackster_time_.clear();
      Trackster_timeErr_.clear();
      Trackster_isSeed_.clear();
      Trackster_inCluster_.clear();

      for (auto& tst: tracksters){
        if (tst.vertices().empty()) continue;
        auto momentum = tst.eigenvectors(0);
        float tst_dr = reco::deltaR(eg_eta, eg_phi, momentum.eta(), momentum.phi());
        if (tst_dr > max_deltaR_) continue;
        Trackster_dr_.push_back(tst_dr);
        Trackster_pt_.push_back(tst.raw_pt());
        Trackster_eta_.push_back(momentum.eta());
        Trackster_phi_.push_back(momentum.phi());
        Trackster_time_.push_back(tst.time());
        Trackster_timeErr_.push_back(tst.timeError());
       

        bool isSeedTrackster = false;
        bool TracksterInCluster = false;
        for (const auto lcId : tst.vertices()){
          if (isSeedTrackster && TracksterInCluster) break;
          isSeedTrackster = (isSeedTrackster | LayerCluster_isSeed[lcId]);
          TracksterInCluster = (TracksterInCluster | LayerCluster_inCluster[lcId]);
        }
        Trackster_isSeed_.push_back(isSeedTrackster);
        Trackster_inCluster_.push_back(TracksterInCluster);


      }




      tree_endcap->Fill();
    }


  }


}



void EGammaNtuples::beginJob() {


  tree_barrel->Branch("eg_et", &eg_et);
  tree_barrel->Branch("eg_energy", &eg_energy);
  tree_barrel->Branch("eg_eta", &eg_eta);
  tree_barrel->Branch("eg_phi", &eg_phi);
  tree_barrel->Branch("eg_status", &eg_status);
  tree_barrel->Branch("eg_sigmaIEtaIEta", &eg_sigmaIEtaIEta);
  tree_barrel->Branch("eg_sigmaIPhiIPhi", &eg_sigmaIPhiIPhi);
  tree_barrel->Branch("eg_sigmaIEtaIEtaNoiseCleaned", &eg_sigmaIEtaIEtaNoiseCleaned);
  tree_barrel->Branch("eg_sigmaIPhiIPhiNoiseCleaned", &eg_sigmaIPhiIPhiNoiseCleaned);


  tree_endcap->Branch("eg_et", &eg_et);
  tree_endcap->Branch("eg_energy", &eg_energy);
  tree_endcap->Branch("eg_eta", &eg_eta);
  tree_endcap->Branch("eg_phi", &eg_phi);
  tree_endcap->Branch("eg_status", &eg_status);
  tree_endcap->Branch("eg_sigma2uu", &eg_sigma2uu);
  tree_endcap->Branch("eg_sigma2ww", &eg_sigma2ww);
  tree_endcap->Branch("eg_sigma2vv", &eg_sigma2vv);


  tree_endcap->Branch("HGCRecHit_dr", &HGCRecHit_dr_);
  tree_endcap->Branch("HGCRecHit_energy", &HGCRecHit_energy_);
  tree_endcap->Branch("HGCRecHit_eta", &HGCRecHit_eta_);
  tree_endcap->Branch("HGCRecHit_phi", &HGCRecHit_phi_);
  tree_endcap->Branch("HGCRecHit_time", &HGCRecHit_time_);
  tree_endcap->Branch("HGCRecHit_timeErr", &HGCRecHit_timeErr_);
  tree_endcap->Branch("HGCRecHit_layer", &HGCRecHit_layer_);
  tree_endcap->Branch("HGCRecHit_isSeed", &HGCRecHit_isSeed_);
  tree_endcap->Branch("HGCRecHit_inCluster", &HGCRecHit_inCluster_);


  tree_endcap->Branch("layerCluster_dr", &LayerCluster_dr_);
  tree_endcap->Branch("layerCluster_energy", &LayerCluster_energy_);
  tree_endcap->Branch("layerCluster_eta", &LayerCluster_eta_);
  tree_endcap->Branch("layerCluster_phi", &LayerCluster_phi_);
  tree_endcap->Branch("layerCluster_time", &LayerCluster_time_);
  tree_endcap->Branch("layerCluster_isSeed", &LayerCluster_isSeed_);
  tree_endcap->Branch("layerCluster_inCluster", &LayerCluster_inCluster_);

  tree_endcap->Branch("trackster_dr", &Trackster_dr_);
  tree_endcap->Branch("trackster_pt", &Trackster_pt_);
  tree_endcap->Branch("trackster_eta", &Trackster_eta_);
  tree_endcap->Branch("trackster_phi", &Trackster_phi_);
  tree_endcap->Branch("trackster_time", &Trackster_time_);
  tree_endcap->Branch("trackster_timeErr", &Trackster_timeErr_);
  tree_endcap->Branch("trackster_isSeed", &Trackster_isSeed_);
  tree_endcap->Branch("trackster_inCluster", &Trackster_inCluster_);



}


void EGammaNtuples::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace edm;
  caloGeometry_ = &iSetup.getData(caloGeomToken_);
  caloTopology_ = &iSetup.getData(caloTopoToken_);

  recHitTools_.setGeometry(*caloGeometry_);

  fillTree(iEvent, iSetup);
  #ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
  #endif
}

// ------------ method called once each job just after ending the event loop  ------------
void EGammaNtuples::endJob() {
  // please remove this method if not needed
  newfile->cd();
  tree_barrel->Write();
  tree_endcap->Write();
  newfile->Close();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void EGammaNtuples::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("genParticles",edm::InputTag("genParticles"));
  desc.add<edm::InputTag>("scBarrelL1Seeded", edm::InputTag("hltParticleFlowSuperClusterECALUnseeded","particleFlowSuperClusterECALBarrel"));
  desc.add<edm::InputTag>("scHGCalL1Seeded", edm::InputTag("hltParticleFlowSuperClusterHGCalFromTICLUnseeded",""));
  desc.add<edm::InputTag>("ebRecHits", edm::InputTag("hltEgammaHLTExtra","EcalRecHitsEB"));
  desc.add<edm::InputTag>("HGCeeRecHits", edm::InputTag("HGCalRecHit","HGCEERecHits"));
  desc.add<edm::InputTag>("HGChebRecHits", edm::InputTag("HGCalRecHit","HGCEERecHits"));
  desc.add<edm::InputTag>("HGChefRecHits", edm::InputTag("HGCalRecHit","HGCEERecHits"));
  desc.add<edm::InputTag>("sigmaIEtaIEta", edm::InputTag("hltEgammaClusterShapeUnseeded","sigmaIEtaIEta5x5"));
  desc.add<edm::InputTag>("sigmaIPhiIPhi", edm::InputTag("hltEgammaClusterShapeUnseeded","sigmaIPhiIPhi5x5"));
  desc.add<edm::InputTag>("sigmaIEtaIEtaNoiseCleaned", edm::InputTag("hltEgammaClusterShapeUnseeded","sigmaIEtaIEta5x5NoiseCleaned"));
  desc.add<edm::InputTag>("sigmaIPhiIPhiNoiseCleaned", edm::InputTag("hltEgammaClusterShapeUnseeded","sigmaIPhiIPhi5x5NoiseCleaned"));
  desc.add<edm::InputTag>("nrHitsEB1GeV", edm::InputTag("hltEgammaHLTExtra","countEcalRecHitsEcalRecHitsEBThres1GeV"));
  desc.add<edm::InputTag>("nrHitsEE1GeV", edm::InputTag("hltEgammaHLTExtra","countEcalRecHitsHGCalRecHitsThres1GeV"));
  desc.add<edm::InputTag>("eGammaCandidates", edm::InputTag("hltEgammaCandidatesUnseeded"));
  desc.add<edm::InputTag>("vertices", edm::InputTag("offlinePrimaryVertices"));
  desc.add<edm::InputTag>("pfHGCALRecHits", edm::InputTag("particleFlowRecHitHGC"));
  desc.add<edm::InputTag>("trackster", edm::InputTag("hltTiclCandidate"));
  desc.add<edm::InputTag>("layerClusters", edm::InputTag("hltParticleFlowClusterECAL"));
  desc.add<edm::InputTag>("timelayerClusters", edm::InputTag("hltMergeLayerClusters", "timeLayerCluster"));
  desc.add<std::string>("pType","ele");
  desc.add<double>("ptThreshold", 10);
  desc.add<double>("genMatchDR", 0.1);
  desc.add<double>("maxDR", 0.4);
  descriptions.add("EGammaNtuples",desc);


  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}


int EGammaNtuples::matchToTruth(const reco::RecoEcalCandidate &egamma, const std::vector<reco::GenParticle>& genParticles) const{

  double dR = 999;
  reco::GenParticle const* closestElectron = nullptr;
  int pdgID = -1;
  if (pType_ == "ele") pdgID = 11;
  else if (pType_ == "pho") pdgID = 22;

  for (auto const& particle : genParticles){
    if (std::abs(particle.pdgId()) != pdgID || particle.status() != 1) continue;
    double dRtmp = ROOT::Math::VectorUtil::DeltaR(egamma.p4(), particle.p4());
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

double EGammaNtuples::getLayerClusterTime( const reco::CaloCluster &lc, const std::unordered_map<DetId, HGCRecHit> &recHitMap) {
    double timeSum = 0.0;
    double energySum = 0.0;

    for (auto const& hitFraction : lc.hitsAndFractions()) {
        const DetId& detid = hitFraction.first;
        double fraction = hitFraction.second;

        auto it = recHitMap.find(detid);
        if (it == recHitMap.end()) continue; // skip if hit not found

        const auto& recHit = it->second;

        double energy = recHit.energy() * fraction;
        double time   = recHit.time();

        if (time < -90) continue; // skip invalid times (CMS convention: -99 is invalid)


        timeSum += energy * time;
        energySum += energy;

    }

    return (energySum > 0) ? timeSum / energySum : -99;
}
//define this as a plug-in
DEFINE_FWK_MODULE(EGammaNtuples);

