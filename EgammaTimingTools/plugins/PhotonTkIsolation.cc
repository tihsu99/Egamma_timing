//*****************************************************************************
// File:      PhotonTkIsolation.cc
// ----------------------------------------------------------------------------
// OrigAuth:  Matthias Mozer
// Institute: IIHE-VUB
//=============================================================================
//*****************************************************************************

#include "PhotonTkIsolation.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaTrackSelector.h"

#include <Math/VectorUtil.h>

PhotonTkIsolation::PhotonTkIsolation(float extRadius,
                                     float intRadiusBarrel,
                                     float intRadiusEndcap,
                                     float stripBarrel,
                                     float stripEndcap,
                                     float etLow,
                                     float lip,
                                     float drb,
                                     const reco::TrackCollection* trackCollection,
                                     reco::TrackBase::Point beamPoint,
                                     const std::string& dzOptionString)
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
  setDzOption(dzOptionString);
}

void PhotonTkIsolation::setDzOption(const std::string& s) {
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

PhotonTkIsolation::~PhotonTkIsolation() {}

// unified acces to isolations
std::pair<int, float> PhotonTkIsolation::getIso(const reco::Candidate* photon) const {
  int counter = 0;
  float ptSum = 0.;

  //Take the photon position
  float photonEta = photon->eta();

  // Find the signal track
  int sigTrkIdx = -1;
  const reco::Track* sigTrkPtr = 0;
  if (dtRef_ == TRK_DT_REF::SIG) {
    int trkIdx = -1;
    for (reco::TrackCollection::const_iterator itrTr = (*trackCollection_).begin(); itrTr != (*trackCollection_).end(); ++itrTr){
      trkIdx++;
      double dr = ROOT::Math::VectorUtil::DeltaR(itrTr->momentum(), photon->momentum());
      bool   isBarrel = std::abs(photonEta) < 1.479;
      double intRadius = isBarrel ? intRadiusBarrel2_ : intRadiusEndcap2_;
      if (dr*dr < intRadius){
        if (!sigTrkPtr || (*itrTr).pt() > (*sigTrkPtr).pt()){
          sigTrkPtr = &(*itrTr);
          sigTrkIdx = trkIdx;
        }
      }
    }
  }



  // Collect signal track timing information
  reco::TrackRef sigTrkRef;
  double sigTrkTime = -1;
  double sigTrkTimeErr = -1;
  double sigTrkMtdMva = -1;
  if (sigTrkIdx >= 0){
    sigTrkRef = reco::TrackRef(trackCollectionH_, sigTrkIdx);
    sigTrkTime = mtdt0_[sigTrkRef];
    sigTrkMtdMva = mtdTrkQualMVA_[sigTrkRef];
    sigTrkTimeErr = (sigTrkMtdMva > trkMtdMvaMin_)? sigTrkTimeErr: -1;
  }

  int trkIdx = -1;
  //loop over tracks
  for (reco::TrackCollection::const_iterator trItr = trackCollection_->begin(); trItr != trackCollection_->end();
       ++trItr) {

    trkIdx ++; 
    const  reco::TrackRef trkRef(trackCollectionH_, trkIdx);

    //check z-distance of vertex
    float dzCut = 0;
    switch (dzOption_) {
      case egammaisolation::EgammaTrackSelector::dz:
        dzCut = fabs((*trItr).dz() - photon->vertex().z());
        break;
      case egammaisolation::EgammaTrackSelector::vz:
        dzCut = fabs((*trItr).vz() - photon->vertex().z());
        break;
      case egammaisolation::EgammaTrackSelector::bs:
        dzCut = fabs((*trItr).dz(beamPoint_) - photon->vertex().z());
        break;
      case egammaisolation::EgammaTrackSelector::vtx:
        dzCut = fabs((*trItr).dz(photon->vertex()));
        break;
      default:
        dzCut = fabs((*trItr).vz() - photon->vertex().z());
        break;
    }
    if (dzCut > lip_)
      continue;

    float this_pt = (*trItr).pt();
    if (this_pt < etLow_)
      continue;
    if (fabs((*trItr).dxy(beamPoint_)) > drb_)
      continue;  // only consider tracks from the main vertex

    // Timing variable cut
    double dt = 0;
    double trkTime = mtdt0_[trkRef];
    double trkTimeErr = mtdSigmat0_[trkRef];
    double trkMtdMva  = mtdTrkQualMVA_[trkRef];

    trkTimeErr = (trkMtdMva > trkMtdMvaMin_)? trkTimeErr: -1;
    if (dtRef_ == TRK_DT_REF::PV){
      double vtxTimeErr = vertex_.tError();
      if(dtMax_>0 && trkTimeErr > 0 && vtxTimeErr > 0)
      {
        dt = std::abs(trkTime - vertex_.t());
        if (dtType_ == TRK_DT_TYPE::SIGNIFICANCE){
          dt /= std::sqrt(trkTimeErr * trkTimeErr + vtxTimeErr * vtxTimeErr);
        }
        if (dt > dtMax_){
          continue;
        };
      }
    }

    else if (dtRef_ == TRK_DT_REF::SIG && sigTrkPtr){
      if (dtMax_>0 && trkTimeErr > 0 && sigTrkTimeErr > 0){
        dt = std::fabs(trkTime - sigTrkTime);
        if(dtType_ == TRK_DT_TYPE::SIGNIFICANCE){
          dt /= std::sqrt(trkTimeErr * trkTimeErr + sigTrkTimeErr * sigTrkTimeErr);
        }
        if (dt > dtMax_){
          continue;
        }
      }
    }



    float dr2 = reco::deltaR2(*trItr, *photon);
    float deta = (*trItr).eta() - photonEta;
    if (fabs(photonEta) < 1.479) {
      if (dr2 < extRadius2_ && dr2 >= intRadiusBarrel2_ && fabs(deta) >= stripBarrel_) {
        ++counter;
        ptSum += this_pt;
      }
    } else {
      if (dr2 < extRadius2_ && dr2 >= intRadiusEndcap2_ && fabs(deta) >= stripEndcap_) {
        ++counter;
        ptSum += this_pt;
      }
    }

  }  //end loop over tracks
  std::pair<int, float> retval;
  retval.first = counter;
  retval.second = ptSum;
  return retval;
}
