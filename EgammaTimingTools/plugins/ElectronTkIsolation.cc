//*****************************************************************************
// File:      ElectronTkIsolation.cc
// ----------------------------------------------------------------------------
// OrigAuth:  Matthias Mozer
// Institute: IIHE-VUB
//=============================================================================
//*****************************************************************************

#include "ElectronTkIsolation.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"

#include <Math/VectorUtil.h>

ElectronTkIsolation::ElectronTkIsolation(double extRadius,
                                         double intRadiusBarrel,
                                         double intRadiusEndcap,
                                         double stripBarrel,
                                         double stripEndcap,
                                         double ptLow,
                                         double lip,
                                         double drb,
                                         const reco::TrackCollection* trackCollection,
                                         reco::TrackBase::Point beamPoint,
                                         const std::string& dzOptionString)
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
  setDzOption(dzOptionString);
}

ElectronTkIsolation::~ElectronTkIsolation() {}

std::pair<int, double> ElectronTkIsolation::getIso(const reco::GsfElectron* electron) const {
  return getIso(&(*(electron->gsfTrack())));
}

// unified acces to isolations
std::pair<int, double> ElectronTkIsolation::getIso(const reco::Track* tmpTrack) const {
  int counter = 0;
  double ptSum = 0.;
  //Take the electron track
  math::XYZVector tmpElectronMomentumAtVtx = (*tmpTrack).momentum();
  double tmpElectronEtaAtVertex = (*tmpTrack).eta();

  int sigTrkIdx = -1;
  const reco::Track* sigTrkPtr = 0;
  
  // Find the signal track
  if (dtRef_ == TRK_DT_REF::SIG) {
    
    int trkIdx = -1;
    
    for (reco::TrackCollection::const_iterator itrTr = (*trackCollection_).begin(); itrTr != (*trackCollection_).end();
         ++itrTr) {
      trkIdx++;
      double dr = ROOT::Math::VectorUtil::DeltaR(itrTr->momentum(), tmpElectronMomentumAtVtx);
      //double deta = (*itrTr).eta() - tmpElectronEtaAtVertex;
      bool isBarrel = std::abs(tmpElectronEtaAtVertex) < 1.479;
      double intRadius = isBarrel ? intRadiusBarrel_ : intRadiusEndcap_;
      //double strip = isBarrel ? stripBarrel_ : stripEndcap_;
      if (dr < intRadius) {
        //printf("Sig trk. dr %f (%f), pt %f (tmpTrack pt %f), t0Error %f \n", dr, intRadius, (*itrTr).pt(), (*tmpTrack).pt(), (*itrTr).t0Error());
        if (!sigTrkPtr || (*itrTr).pt() > (*sigTrkPtr).pt()) {
          sigTrkPtr = &(*itrTr);
          sigTrkIdx = trkIdx;
        }
      }
    }
  }
  
  reco::TrackRef sigTrkRef;
  double sigTrkTime = -1;
  double sigTrkTimeErr = -1;
  double sigTrkMtdMva = -1;
  
  if (sigTrkIdx >= 0) {
    sigTrkRef = reco::TrackRef(trackCollectionH_, sigTrkIdx);
   // printf("Sig: trk %d, mtdt0 %f, mtdSigmat0 %f, mtdTrkQualMVA %f \n", sigTrkIdx, mtdt0_[sigTrkRef], mtdSigmat0_[sigTrkRef], mtdTrkQualMVA_[sigTrkRef]);
    
    sigTrkTime = mtdt0_[sigTrkRef];
    sigTrkTimeErr = mtdSigmat0_[sigTrkRef];
    sigTrkMtdMva = mtdTrkQualMVA_[sigTrkRef];
    
    sigTrkTimeErr = (sigTrkMtdMva > trkMtdMvaMin_)? sigTrkTimeErr: -1;
  }

  int trkIdx = -1;
  
  //printf("trackCollection_.size() %d, mtdt0_.size() %d, mtdSigmat0_.size() %d, mtdTrkQualMVA_.size() %d \n", (int) trackCollection_->size(), (int) mtdt0_.size(), (int) mtdSigmat0_.size(), (int) mtdTrkQualMVA_.size());
  
  for (reco::TrackCollection::const_iterator itrTr = (*trackCollection_).begin(); itrTr != (*trackCollection_).end();
       ++itrTr) {
    
    trkIdx++;
    const reco::TrackRef trkRef(trackCollectionH_, trkIdx);
    
    double this_pt = (*itrTr).pt();
    if (this_pt < ptLow_)
      continue;

    double dzCut = 0;
    switch (dzOption_) {
      case egammaisolation::EgammaTrackSelector::dz:
        dzCut = fabs((*itrTr).dz() - (*tmpTrack).dz());
        break;
      case egammaisolation::EgammaTrackSelector::vz:
        dzCut = fabs((*itrTr).vz() - (*tmpTrack).vz());
        break;
      case egammaisolation::EgammaTrackSelector::bs:
        dzCut = fabs((*itrTr).dz(beamPoint_) - (*tmpTrack).dz(beamPoint_));
        break;
      case egammaisolation::EgammaTrackSelector::vtx:
        dzCut = fabs((*itrTr).dz(tmpTrack->vertex()));
        break;
      default:
        dzCut = fabs((*itrTr).vz() - (*tmpTrack).vz());
        break;
    }
    if (dzCut > lip_)
      continue;
    if (fabs((*itrTr).dxy(beamPoint_)) > drb_)
      continue;
    
    double dt = 0;
    double trkTime = mtdt0_[trkRef];
    double trkTimeErr = mtdSigmat0_[trkRef];
    double trkMtdMva = mtdTrkQualMVA_[trkRef];
    //double trkTimeErr = (*itrTr).t0Error();
    
    trkTimeErr = (trkMtdMva > trkMtdMvaMin_)? trkTimeErr: -1;
    
    if (dtRef_ == TRK_DT_REF::PV)
    {
        double vtxTimeErr = vertex_.tError();
        
        if(dtMax_> 0 && trkTimeErr > 0 && vtxTimeErr > 0)
        {
            dt = std::fabs(trkTime - vertex_.t());
            
            if (dtType_ == TRK_DT_TYPE::SIGNIFICANCE)
            {
                dt /= std::sqrt(trkTimeErr*trkTimeErr + vtxTimeErr*vtxTimeErr);
            }
            
            //printf("dt (dtMax) %f (%f) \n", dt, dtMax_);
            
            if (dt > dtMax_) {
                //printf("Skipping track. \n");
                continue;
            };
        }
    }
    
    else if (dtRef_ == TRK_DT_REF::SIG && sigTrkPtr)
    {
        //printf("trkTimeErr %f, sigTrkTimeErr %f \n", trkTimeErr, sigTrkTimeErr);
        
        if(dtMax_> 0 && trkTimeErr > 0 && sigTrkTimeErr > 0)
        {
            dt = std::fabs(trkTime - sigTrkTime);
            
            if (dtType_ == TRK_DT_TYPE::SIGNIFICANCE)
            {
                dt /= std::sqrt(trkTimeErr*trkTimeErr + sigTrkTimeErr*sigTrkTimeErr);
            }
            
            //printf("dt (dtMax) %f (%f) \n", dt, dtMax_);
            
            if (dt > dtMax_) {
                //printf("Skipping track. \n");
                continue;
            };
        }
    }
    
    double dr = ROOT::Math::VectorUtil::DeltaR(itrTr->momentum(), tmpElectronMomentumAtVtx);
    double deta = (*itrTr).eta() - tmpElectronEtaAtVertex;
    bool isBarrel = std::abs(tmpElectronEtaAtVertex) < 1.479;
    double intRadius = isBarrel ? intRadiusBarrel_ : intRadiusEndcap_;
    double strip = isBarrel ? stripBarrel_ : stripEndcap_;
    if (dr < extRadius_ && dr >= intRadius && std::abs(deta) >= strip && passAlgo(*itrTr)) {
      ++counter;
      ptSum += this_pt;
    }

  }  //end loop over tracks

  //printf("Isolation %f", ptSum);
  
  std::pair<int, double> retval;
  retval.first = counter;
  retval.second = ptSum;

  return retval;
}

int ElectronTkIsolation::getNumberTracks(const reco::GsfElectron* electron) const {
  //counter for the tracks in the isolation cone
  return getIso(electron).first;
}

double ElectronTkIsolation::getPtTracks(const reco::GsfElectron* electron) const { return getIso(electron).second; }

bool ElectronTkIsolation::passAlgo(const reco::TrackBase& trk) const {
  int algo = trk.algo();
  bool rejAlgo = std::binary_search(algosToReject_.begin(), algosToReject_.end(), algo);
  return rejAlgo == false;
}

void ElectronTkIsolation::setAlgosToReject() {
  algosToReject_ = {reco::TrackBase::jetCoreRegionalStep};
  std::sort(algosToReject_.begin(), algosToReject_.end());
}
