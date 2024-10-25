#ifndef PhotonTkIsolation_h
#define PhotonTkIsolation_h

//*****************************************************************************
// File:      PhotonTkIsolation.h
// ----------------------------------------------------------------------------
// OrigAuth:  Matthias Mozer
// Institute: IIHE-VUB
//=============================================================================
//*****************************************************************************

//C++ includes
#include <string>

//CMSSW includes
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaTrackSelector.h"

class PhotonTkIsolation {
public:


  enum TRK_DT_TYPE {
      ABSOLUTE = 0,
      SIGNIFICANCE
  };
              
  enum TRK_DT_REF {
      PV = 0,
      SIG
  };

  //constructors
  PhotonTkIsolation(float extRadius,
                    float intRadius,
                    float etLow,
                    float lip,
                    float drb,
                    const reco::TrackCollection* trackCollection,
                    reco::TrackBase::Point beamPoint)
      : extRadius2_(extRadius * extRadius),
        intRadiusBarrel2_(intRadius * intRadius),
        intRadiusEndcap2_(intRadius * intRadius),
        stripBarrel_(0.0),
        stripEndcap_(0.0),
        etLow_(etLow),
        lip_(lip),
        drb_(drb),
        trackCollection_(trackCollection),
        beamPoint_(beamPoint) {
    setDzOption("vz");
  }

  PhotonTkIsolation(float extRadius,
                    float intRadius,
                    float strip,
                    float etLow,
                    float lip,
                    float drb,
                    const reco::TrackCollection* trackCollection,
                    reco::TrackBase::Point beamPoint)
      : extRadius2_(extRadius * extRadius),
        intRadiusBarrel2_(intRadius * intRadius),
        intRadiusEndcap2_(intRadius * intRadius),
        stripBarrel_(strip),
        stripEndcap_(strip),
        etLow_(etLow),
        lip_(lip),
        drb_(drb),
        trackCollection_(trackCollection),
        beamPoint_(beamPoint) {
    setDzOption("vz");
  }

  PhotonTkIsolation(float extRadius,
                    float intRadiusBarrel,
                    float intRadiusEndcap,
                    float stripBarrel,
                    float stripEndcap,
                    float etLow,
                    float lip,
                    float drb,
                    const reco::TrackCollection* trackCollection,
                    reco::TrackBase::Point beamPoint)
      : extRadius2_(extRadius * extRadius),
        intRadiusBarrel2_(intRadiusBarrel * intRadiusBarrel),
        intRadiusEndcap2_(intRadiusEndcap * intRadiusEndcap),
        stripBarrel_(stripBarrel),
        stripEndcap_(stripEndcap),
        etLow_(etLow),
        lip_(lip),
        drb_(drb),
        trackCollection_(trackCollection),
        beamPoint_(beamPoint) {
    setDzOption("vz");
  }

  PhotonTkIsolation(float extRadius,
                    float intRadiusBarrel,
                    float intRadiusEndcap,
                    float stripBarrel,
                    float stripEndcap,
                    float etLow,
                    float lip,
                    float drb,
                    int   dtRef,
                    int   dtType,
                    double dtMax,
                    double trkMtdMvaMin_,
//                    const reco::TrackCollection* trackCollection,
                    const edm::Handle <reco::TrackCollection> &trackCollectionH,
                    const edm::ValueMap<float> &mtdt0,
                    const edm::ValueMap<float> &mtdSigmat0,
                    const edm::ValueMap<float> &mtdTrkQualMVA,
                    reco::TrackBase::Point beamPoint,
                    const reco::Vertex &vertex)

      : extRadius2_(extRadius * extRadius),
        intRadiusBarrel2_(intRadiusBarrel * intRadiusBarrel),
        intRadiusEndcap2_(intRadiusEndcap * intRadiusEndcap),
        stripBarrel_(stripBarrel),
        stripEndcap_(stripEndcap),
        etLow_(etLow),
        lip_(lip),
        drb_(drb),
        dtRef_(dtRef),
        dtType_(dtType),
        dtMax_(dtMax),
        trkMtdMvaMin_(trkMtdMvaMin_),
        trackCollection_(trackCollectionH.product()),
        trackCollectionH_(trackCollectionH),
        mtdt0_(mtdt0),
        mtdSigmat0_(mtdSigmat0),
        mtdTrkQualMVA_(mtdTrkQualMVA),
        beamPoint_(beamPoint),
        vertex_(vertex){
    setDzOption("vz");
  }



  PhotonTkIsolation(float extRadius,
                    float intRadiusBarrel,
                    float intRadiusEndcap,
                    float stripBarrel,
                    float stripEndcap,
                    float etLow,
                    float lip,
                    float drb,
                    const reco::TrackCollection*,
                    reco::TrackBase::Point beamPoint,
                    const std::string&);

  //destructor
  ~PhotonTkIsolation();
  //methods

  std::pair<int, float> getIso(const reco::Candidate*) const;

  void setDzOption(const std::string& s);

private:
  float extRadius2_;
  float intRadiusBarrel2_;
  float intRadiusEndcap2_;
  float stripBarrel_;
  float stripEndcap_;
  float etLow_;
  float lip_;
  float drb_;

  int dtRef_;
  int dtType_;
  double dtMax_;
  double trkMtdMvaMin_;
  const reco::TrackCollection* trackCollection_;
  const edm::Handle <reco::TrackCollection> trackCollectionH_;
  const edm::ValueMap<float> mtdt0_;
  const edm::ValueMap<float> mtdSigmat0_;
  const edm::ValueMap<float> mtdTrkQualMVA_;
  reco::TrackBase::Point beamPoint_;
  reco::Vertex vertex_;
  int dzOption_;
};

#endif
