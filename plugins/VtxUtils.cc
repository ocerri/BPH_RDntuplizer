#ifndef __VTXUTILS__
#define __VTXUTILS__

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"

// Needed for Transient Tracks
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

// Needed for the kinematic fit
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"

#include <iostream>

using namespace std;

#define _KMass_ 0.493677
#define _KMassErr_ 0.000013
#define _PiMass_ 0.13957018
#define _PiMassErr_ 0.00000035

namespace vtxu {

  bool FitD0(const edm::EventSetup& iSetup, pat::PackedCandidate pi, pat::PackedCandidate K, int verbose = 0) {
    // Get transient track builder
    edm::ESHandle<TransientTrackBuilder> TTBuilder;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",TTBuilder);

    if (pi.charge()*K.charge() > 0) return false;

    reco::TransientTrack K_tk = TTBuilder->build(K.bestTrack());
    reco::TransientTrack pi_tk = TTBuilder->build(pi.bestTrack());

    KinematicParticleFactoryFromTransientTrack pFactory;
    KinematicParticleVertexFitter VtxFitter;

    std::vector<RefCountedKinematicParticle> parts;
    double chi = 0, ndf = 0;
    float mK = _KMass_, dmK = _KMass_;
    float mPi = _PiMass_, dmPi = _PiMass_;
    parts.push_back(pFactory.particle(K_tk, mK, chi, ndf, dmK));
    parts.push_back(pFactory.particle(pi_tk, mPi, chi, ndf, dmPi));
    RefCountedKinematicTree D0KinTree = VtxFitter.fit(parts);

    if (verbose) {
      cout << "D0KinTree->isValid(): " << D0KinTree->isValid() << endl;
    }

    D0KinTree->movePointerToTheTop();
    auto D0vtx = D0KinTree->currentDecayVertex();

    if (verbose) {
      cout << "D0vtx->chiSquared(): " << D0vtx->chiSquared() << endl;
    }

    return true;
  }
}

#endif
