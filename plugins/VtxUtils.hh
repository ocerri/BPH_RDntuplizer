#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicVertex.h"
#include <DataFormats/TrackReco/interface/Track.h>
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

// Pure ROOT import
#include <TLorentzVector.h>
#include <utility>

namespace vtxu {
  struct kinFitResuts{
    bool isValid = false;
    double chi2 = -1;
    double dof = -1;
    double pval = -1;
    bool isGood = false;
  };

  /// error on dxy with respect to a user-given reference point + uncertainty (i.e. reco::Vertex position)
  double dxyError(const reco::TrackBase &tk, const reco::TrackBase::Point &vtx, math::Error<3>::type &vertexCov);

  // error on dxy with respect to a given beamspot
  double dxyError(const reco::TrackBase &tk, const reco::BeamSpot &theBeamSpot);

  RefCountedKinematicTree FitD0(const edm::EventSetup&, pat::PackedCandidate, pat::PackedCandidate, bool);
  RefCountedKinematicTree FitKst_piK(const edm::EventSetup&, pat::PackedCandidate, pat::PackedCandidate, bool);
  RefCountedKinematicTree FitPhi_KK(const edm::EventSetup&, pat::PackedCandidate, pat::PackedCandidate, bool);
  RefCountedKinematicTree FitDst_fitD0wMassConstraint(const edm::EventSetup&, pat::PackedCandidate, pat::PackedCandidate, pat::PackedCandidate, bool, int);
  RefCountedKinematicTree FitJpsi_mumu(const edm::EventSetup&, pat::Muon, pat::Muon, bool);
  RefCountedKinematicTree FitB_mumuK(const edm::EventSetup&, pat::Muon, pat::Muon, pat::PackedCandidate);
  RefCountedKinematicTree FitB_mumupiK(const edm::EventSetup&, pat::Muon, pat::Muon, pat::PackedCandidate, pat::PackedCandidate, bool, bool, bool, reco::Vertex* = nullptr);
  RefCountedKinematicTree FitB_mumupiK(const RefCountedKinematicParticle, const RefCountedKinematicParticle, const RefCountedKinematicParticle, const RefCountedKinematicParticle, bool, bool, bool, reco::Vertex* = nullptr);
  RefCountedKinematicTree FitB_mumupiK_tk(const edm::EventSetup&, pat::Muon, pat::Muon, pat::PackedCandidate, pat::PackedCandidate, pat::PackedCandidate, bool, bool, bool, reco::Vertex* = nullptr);
  RefCountedKinematicTree FitDst(const edm::EventSetup&, pat::PackedCandidate, const RefCountedKinematicParticle, bool);
  RefCountedKinematicTree FitD0pipi(const edm::EventSetup&, const RefCountedKinematicParticle, pat::PackedCandidate, pat::PackedCandidate);
  RefCountedKinematicTree FitD0pipipis(const edm::EventSetup&, const RefCountedKinematicParticle, pat::PackedCandidate, pat::PackedCandidate, pat::PackedCandidate);
  RefCountedKinematicTree FitD0_pions(const edm::EventSetup&, const RefCountedKinematicParticle, std::vector<pat::PackedCandidate>);
  RefCountedKinematicTree FitB_D0pismu(const edm::EventSetup&, const RefCountedKinematicParticle, pat::PackedCandidate, pat::Muon);
  RefCountedKinematicTree Fit_D0pismupi(const edm::EventSetup&, const RefCountedKinematicParticle, pat::PackedCandidate, pat::Muon, pat::PackedCandidate);
  RefCountedKinematicTree Fit_D0pihpis(const edm::EventSetup&, const RefCountedKinematicParticle, pat::PackedCandidate, pat::PackedCandidate);
  RefCountedKinematicTree FitVtxMuDst(const edm::EventSetup&, const RefCountedKinematicParticle, pat::Muon);
  RefCountedKinematicTree FitVtxJpsiKst(const edm::EventSetup&, const RefCountedKinematicParticle, const RefCountedKinematicParticle, bool);
  RefCountedKinematicTree FitVtxDstK(const edm::EventSetup&, const RefCountedKinematicParticle, pat::PackedCandidate, int);
  RefCountedKinematicTree FitVtxDstPi(const edm::EventSetup&, const RefCountedKinematicParticle, pat::PackedCandidate, int);
  RefCountedKinematicTree FitVtxMuDstPi(const edm::EventSetup&, const RefCountedKinematicParticle, pat::PackedCandidate, pat::PackedCandidate, int);
  std::pair<double,double> computeIP3D(reco::TransientTrack, GlobalVector, reco::Vertex);
  std::pair<double,double> computeIP3D(const edm::EventSetup&, pat::Muon, TLorentzVector, RefCountedKinematicVertex);
  std::pair<double,double> computeIP3D(const edm::EventSetup&, pat::PackedCandidate, TLorentzVector, RefCountedKinematicVertex);
  std::pair<double,double> computeDCA(reco::TransientTrack tk, GlobalPoint p, int kind=0);
  std::pair<double,double> computeDCA(const edm::EventSetup& iSetup, pat::PackedCandidate cand, GlobalPoint p, int kind=0);
  std::pair<double,double> computeDCA(const edm::EventSetup& iSetup, pat::Muon mu, GlobalPoint p, int kind=0);
  std::pair<double,double> vtxsDistance(reco::VertexRef, RefCountedKinematicVertex);
  std::pair<double,double> vtxsDistance(reco::Vertex, RefCountedKinematicVertex);
  std::pair<double,double> vtxsTransverseDistance(reco::Vertex, RefCountedKinematicVertex);
  kinFitResuts fitQuality(RefCountedKinematicTree, double = -1);
  double computePointingCos(reco::Vertex, const RefCountedKinematicVertex, const RefCountedKinematicParticle);
  double computePointingCosTransverse(reco::Vertex, const RefCountedKinematicVertex, const RefCountedKinematicParticle);
  TLorentzVector getTLVfromKinPart(const RefCountedKinematicParticle);
  TLorentzVector getTLVfromTrack(reco::Track, double);
  TLorentzVector getTLVfromCand(pat::PackedCandidate, double);
  TLorentzVector getTLVfromMuon(pat::Muon, double);
  double dPhi(double, double);
  double dR(double, double, double, double);
  double computeIP(reco::Candidate::Point, reco::Candidate::Point, reco::Candidate::Vector, bool);
  double computeDCA_linApprox(reco::Candidate::Point, reco::Candidate::Point, reco::Candidate::Vector, bool);
  double computeCTau(reco::GenParticle);
};
