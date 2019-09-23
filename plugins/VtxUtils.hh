#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicVertex.h"
#include <DataFormats/TrackReco/interface/Track.h>
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

// Pure ROOT import
#include <TLorentzVector.h>
#include <utility>

namespace vtxu {
  RefCountedKinematicTree FitD0(const edm::EventSetup&, pat::PackedCandidate, pat::PackedCandidate, bool, int);
  RefCountedKinematicTree FitKst_piK(const edm::EventSetup&, pat::PackedCandidate, pat::PackedCandidate, bool);
  RefCountedKinematicTree FitDst_fitD0wMassConstraint(const edm::EventSetup&, pat::PackedCandidate, pat::PackedCandidate, pat::PackedCandidate, bool, int);
  RefCountedKinematicTree FitJpsi_mumu(const edm::EventSetup&, pat::Muon, pat::Muon, bool);
  RefCountedKinematicTree FitDst(const edm::EventSetup&, pat::PackedCandidate, const RefCountedKinematicParticle, bool, int);
  RefCountedKinematicTree FitVtxMuDst(const edm::EventSetup&, const RefCountedKinematicParticle, pat::PackedCandidate, int);
  RefCountedKinematicTree FitVtxJpsiKst(const edm::EventSetup&, const RefCountedKinematicParticle, const RefCountedKinematicParticle, bool);
  RefCountedKinematicTree FitVtxDstK(const edm::EventSetup&, const RefCountedKinematicParticle, pat::PackedCandidate, int);
  RefCountedKinematicTree FitVtxDstPi(const edm::EventSetup&, const RefCountedKinematicParticle, pat::PackedCandidate, int);
  RefCountedKinematicTree FitVtxMuDstPi(const edm::EventSetup&, const RefCountedKinematicParticle, pat::PackedCandidate, pat::PackedCandidate, int);
  std::pair<double,double> computeDCA(const edm::EventSetup&, pat::PackedCandidate, GlobalPoint);
  std::pair<double,double> computeDCA(reco::TransientTrack, GlobalPoint);
  std::pair<double,double> vtxsDistance(reco::VertexRef, RefCountedKinematicVertex);
  std::pair<double,double> vtxsDistance(reco::Vertex, RefCountedKinematicVertex);
  TLorentzVector getTLVfromKinPart(const RefCountedKinematicParticle);
  TLorentzVector getTLVfromTrack(reco::Track, double);
  TLorentzVector getTLVfromCand(pat::PackedCandidate, double);
  double dPhi(double, double);
  double dR(double, double, double, double);
  double computeIP(reco::Candidate::Point, reco::Candidate::Point, reco::Candidate::Vector, bool);
}
