#ifndef EgammaIsolationAlgos_ElectronTkIsolation_h
#define EgammaIsolationAlgos_ElectronTkIsolation_h
//*****************************************************************************
// File:      ElectronTkIsolation.h
// ----------------------------------------------------------------------------
// OrigAuth:  Matthias Mozer
// Institute: IIHE-VUB
//=============================================================================
//*****************************************************************************

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaTrackSelector.h"

#include <string>
#include <vector>

class ElectronTkIsolation {
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
  ElectronTkIsolation(double extRadius,
                      double intRadius,
                      double ptLow,
                      double lip,
                      double drb,
                      const reco::TrackCollection* trackCollection,
                      reco::TrackBase::Point beamPoint)
      : extRadius_(extRadius),
        intRadiusBarrel_(intRadius),
        intRadiusEndcap_(intRadius),
        stripBarrel_(0.0),
        stripEndcap_(0.0),
        ptLow_(ptLow),
        lip_(lip),
        drb_(drb),
        trackCollection_(trackCollection),
        beamPoint_(beamPoint) {
    setAlgosToReject();
    setDzOption("vz");
  }

  ElectronTkIsolation(double extRadius,
                      double intRadiusBarrel,
                      double intRadiusEndcap,
                      double stripBarrel,
                      double stripEndcap,
                      double ptLow,
                      double lip,
                      double drb,
                      const reco::TrackCollection* trackCollection,
                      reco::TrackBase::Point beamPoint)
      : extRadius_(extRadius),
        intRadiusBarrel_(intRadiusBarrel),
        intRadiusEndcap_(intRadiusEndcap),
        stripBarrel_(stripBarrel),
        stripEndcap_(stripEndcap),
        ptLow_(ptLow),
        lip_(lip),
        drb_(drb),
        trackCollection_(trackCollection),
        beamPoint_(beamPoint) {
    setAlgosToReject();
    setDzOption("vz");
  }
  
  ElectronTkIsolation(double extRadius,
                      double intRadiusBarrel,
                      double intRadiusEndcap,
                      double stripBarrel,
                      double stripEndcap,
                      double ptLow,
                      double lip,
                      double drb,
                      int dtRef,
                      int dtType,
                      double dtMax,
                      double trkMtdMvaMin_,
                      //const reco::TrackCollection* trackCollection,
                      const edm::Handle <reco::TrackCollection> &trackCollectionH,
                      const edm::ValueMap<float> &mtdt0,
                      const edm::ValueMap<float> &mtdSigmat0,
                      const edm::ValueMap<float> &mtdTrkQualMVA,
                      reco::TrackBase::Point beamPoint,
                      const reco::Vertex &vertex
                      )
      : extRadius_(extRadius),
        intRadiusBarrel_(intRadiusBarrel),
        intRadiusEndcap_(intRadiusEndcap),
        stripBarrel_(stripBarrel),
        stripEndcap_(stripEndcap),
        ptLow_(ptLow),
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
        vertex_(vertex) {
    setAlgosToReject();
    setDzOption("vz");
  }

  ElectronTkIsolation(double extRadius,
                      double intRadiusBarrel,
                      double intRadiusEndcap,
                      double stripBarrel,
                      double stripEndcap,
                      double ptLow,
                      double lip,
                      double drb,
                      const reco::TrackCollection*,
                      reco::TrackBase::Point beamPoint,
                      const std::string&);

  //destructor
  ~ElectronTkIsolation();

  //methods

  void setDzOption(const std::string& s) {
    if (!s.compare("dz"))
      dzOption_ = egammaisolation::EgammaTrackSelector::dz;
    else if (!s.compare("vz"))
      dzOption_ = egammaisolation::EgammaTrackSelector::vz;
    else if (!s.compare("bs"))
      dzOption_ = egammaisolation::EgammaTrackSelector::bs;
    else if (!s.compare("vtx"))
      dzOption_ = egammaisolation::EgammaTrackSelector::vtx;
    else
      dzOption_ = egammaisolation::EgammaTrackSelector::dz;
  }

  int getNumberTracks(const reco::GsfElectron*) const;
  double getPtTracks(const reco::GsfElectron*) const;
  std::pair<int, double> getIso(const reco::GsfElectron*) const;
  std::pair<int, double> getIso(const reco::Track*) const;

private:
  bool passAlgo(const reco::TrackBase& trk) const;
  void setAlgosToReject();
  double extRadius_;
  double intRadiusBarrel_;
  double intRadiusEndcap_;
  double stripBarrel_;
  double stripEndcap_;
  double ptLow_;
  double lip_;
  double drb_;
  int dtRef_;
  int dtType_;
  double dtMax_;
  double trkMtdMvaMin_;
  std::vector<int> algosToReject_;  //vector is sorted
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
