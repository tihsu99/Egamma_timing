#ifndef EGAMMATIMINGTOOLS_COMPARISONNTUPLIZERUTILS_H
#define EGAMMATIMINGTOOLS_COMPARISONNTUPLIZERUTILS_H

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HLTReco/interface/EgammaObject.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include <cmath>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace egammaTiming {

struct GenMatchResult {
  int index = -1;
  double dr = 999.;
  int truthClass = 0;
  float pt = -1.f;
  float eta = -999.f;
  float phi = -999.f;
};

inline reco::Vertex getPrimaryVertex(const std::vector<reco::Vertex>& vertices, int& nGoodVertex) {
  reco::Vertex primaryVertex;
  nGoodVertex = 0;

  for (const auto& vertex : vertices) {
    const bool isGoodVertex = !vertex.isFake() && vertex.ndof() >= 4;
    nGoodVertex += static_cast<int>(isGoodVertex);
    if (isGoodVertex) {
      primaryVertex = vertex;
      break;
    }
  }

  return primaryVertex;
}

inline float getEgammaVar(const trigger::EgammaObject& object, const std::string& key, float fallback = -999.f) {
  int missing = 0;
  const float value = object.var(key, missing);
  return missing ? fallback : value;
}

template <typename HitCollectionA, typename HitCollectionB>
inline bool haveCommonHit(const HitCollectionA& hitsA, const HitCollectionB& hitsB) {
  for (const auto& hitA : hitsA) {
    for (const auto& hitB : hitsB) {
      if (hitA.first == hitB.first) {
        return true;
      }
    }
  }
  return false;
}

inline int classifyElectronTruth(const reco::GenParticle& particle) {
  if (particle.fromHardProcessFinalState()) {
    return 1;
  }
  if (particle.isDirectHardProcessTauDecayProductFinalState()) {
    return 2;
  }
  return 3;
}

inline int classifyPhotonTruth(const reco::GenParticle& particle) {
  return particle.isPromptFinalState() ? 1 : 2;
}

template <typename CandidateCollection>
inline std::pair<int, double> findBestRecoMatch(float eta,
                                                float phi,
                                                const CandidateCollection& collection,
                                                double maxDeltaR) {
  int bestIndex = -1;
  double bestDr2 = maxDeltaR * maxDeltaR;

  for (unsigned int index = 0; index < collection.size(); ++index) {
    const double dr2 = reco::deltaR2(eta, phi, collection[index].eta(), collection[index].phi());
    if (dr2 < bestDr2) {
      bestDr2 = dr2;
      bestIndex = static_cast<int>(index);
    }
  }

  return {bestIndex, bestIndex >= 0 ? std::sqrt(bestDr2) : 999.};
}

inline GenMatchResult findBestGenMatch(float eta,
                                       float phi,
                                       edm::View<reco::GenParticle> const& genParticles,
                                       int absPdgId,
                                       double maxDeltaR,
                                       int (*classify)(const reco::GenParticle&)) {
  GenMatchResult result;
  double bestDr2 = maxDeltaR * maxDeltaR;

  for (unsigned int index = 0; index < genParticles.size(); ++index) {
    const auto& particle = genParticles[index];
    if (std::abs(particle.pdgId()) != absPdgId || particle.status() != 1) {
      continue;
    }

    const double dr2 = reco::deltaR2(eta, phi, particle.eta(), particle.phi());
    if (dr2 < bestDr2) {
      bestDr2 = dr2;
      result.index = static_cast<int>(index);
      result.dr = std::sqrt(dr2);
      result.truthClass = classify(particle);
      result.pt = particle.pt();
      result.eta = particle.eta();
      result.phi = particle.phi();
    }
  }

  return result;
}

}  // namespace egammaTiming

#endif
