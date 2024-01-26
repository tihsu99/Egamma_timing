//*****************************************************************************
// File:      EgammaElectronTkIsolationProducerNew.cc
// ----------------------------------------------------------------------------
// OrigAuth:  Matthias Mozer
// Institute: IIHE-VUB
//=============================================================================
//*****************************************************************************

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "ElectronTkIsolation.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

class EgammaElectronTkIsolationProducerNew : public edm::global::EDProducer<> {
public:
  explicit EgammaElectronTkIsolationProducerNew(const edm::ParameterSet&);

  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

private:
  const edm::EDGetTokenT<reco::GsfElectronCollection> electronProducer_;
  const edm::EDGetTokenT<reco::TrackCollection> trackProducer_;
  const edm::EDGetTokenT<reco::BeamSpot> beamspotProducer_;
  const edm::EDGetTokenT<std::vector<reco::Vertex>> vertexProducer_;
  
  const edm::EDGetTokenT<edm::ValueMap<float>> mtdt0_;
  const edm::EDGetTokenT<edm::ValueMap<float>> mtdSigmat0_;
  const edm::EDGetTokenT<edm::ValueMap<float>> mtdTrkQualMVA_;
  
  const double ptMin_;
  const double intRadiusBarrel_;
  const double intRadiusEndcap_;
  const double stripBarrel_;
  const double stripEndcap_;
  const double extRadius_;
  const double maxVtxDist_;
  const double drb_;
  
  const int dtRef_;
  const int dtType_;
  const double dtMax_;
  const double trkMtdMvaMin_;
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(EgammaElectronTkIsolationProducerNew);

EgammaElectronTkIsolationProducerNew::EgammaElectronTkIsolationProducerNew(const edm::ParameterSet& config)
    : electronProducer_{consumes(config.getParameter<edm::InputTag>("electronProducer"))},
      trackProducer_{consumes(config.getParameter<edm::InputTag>("trackProducer"))},
      beamspotProducer_{consumes(config.getParameter<edm::InputTag>("BeamspotProducer"))},
      vertexProducer_{consumes(config.getParameter<edm::InputTag>("vertexProducer"))},
      
      mtdt0_{consumes(config.getParameter<edm::InputTag>("mtdt0"))},
      mtdSigmat0_{consumes(config.getParameter<edm::InputTag>("mtdSigmat0"))},
      mtdTrkQualMVA_{consumes(config.getParameter<edm::InputTag>("mtdTrkQualMVA"))},
      
      ptMin_{config.getParameter<double>("ptMin")},
      intRadiusBarrel_{config.getParameter<double>("intRadiusBarrel")},
      intRadiusEndcap_{config.getParameter<double>("intRadiusEndcap")},
      stripBarrel_{config.getParameter<double>("stripBarrel")},
      stripEndcap_{config.getParameter<double>("stripEndcap")},
      extRadius_{config.getParameter<double>("extRadius")},
      maxVtxDist_{config.getParameter<double>("maxVtxDist")},
      drb_{config.getParameter<double>("maxVtxDistXY")},
      
      dtRef_{config.getParameter<int>("dtRef")},
      dtType_{config.getParameter<int>("dtType")},
      dtMax_{config.getParameter<double>("dtMax")},
      trkMtdMvaMin_{config.getParameter<double>("trkMtdMvaMin")}

{
  produces<edm::ValueMap<float>>();
}

void EgammaElectronTkIsolationProducerNew::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup&) const {
  // Get the  filtered objects
  auto electronHandle = iEvent.getHandle(electronProducer_);

  auto vertexHandle = iEvent.getHandle(vertexProducer_);
  std::vector <reco::Vertex> vertices = *vertexHandle;
  
  reco::Vertex pmVtx;
  int nGoodVertex = 0;
  
  for(int iVtx = 0; iVtx < (int) vertices.size(); iVtx++)
  {
      const reco::Vertex &vertex = vertices.at(iVtx);
      
      bool isGoodVertex = (
          !vertex.isFake() &&
          vertex.ndof() >= 4 //&&
          //fabs(vertices.z()) <= 24.0 &&
          //fabs(vertices.position().rho()) <= 2.0
      );
      
      nGoodVertex += (int) isGoodVertex;
      
      if(nGoodVertex)
      {
        pmVtx = vertex;
        break;
      }
  }

  //prepare product
  auto isoMap = std::make_unique<edm::ValueMap<float>>();
  edm::ValueMap<float>::Filler filler(*isoMap);
  std::vector<float> retV(electronHandle->size(), 0);

  //ElectronTkIsolation myTkIsolation(extRadius_,
  //                                  intRadiusBarrel_,
  //                                  intRadiusEndcap_,
  //                                  stripBarrel_,
  //                                  stripEndcap_,
  //                                  ptMin_,
  //                                  maxVtxDist_,
  //                                  drb_,
  //                                  &iEvent.get(trackProducer_),
  //                                  iEvent.get(beamspotProducer_).position()
  //                                  );
  //
  
  ElectronTkIsolation myTkIsolation(extRadius_,
                                    intRadiusBarrel_,
                                    intRadiusEndcap_,
                                    stripBarrel_,
                                    stripEndcap_,
                                    ptMin_,
                                    maxVtxDist_,
                                    drb_,
                                    dtRef_,
                                    dtType_,
                                    dtMax_,
                                    trkMtdMvaMin_,
                                    //&iEvent.get(trackProducer_),
                                    iEvent.getHandle(trackProducer_),
                                    iEvent.get(mtdt0_),
                                    iEvent.get(mtdSigmat0_),
                                    iEvent.get(mtdTrkQualMVA_),
                                    iEvent.get(beamspotProducer_).position(),
                                    pmVtx
                                    );

  for (unsigned int i = 0; i < electronHandle->size(); ++i) {
    double isoValue = myTkIsolation.getPtTracks(&(electronHandle->at(i)));
    retV[i] = isoValue;
  }

  //fill and insert valuemap
  filler.insert(electronHandle, retV.begin(), retV.end());
  filler.fill();
  iEvent.put(std::move(isoMap));
}
