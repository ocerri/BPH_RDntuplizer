#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"

// #include "RecoVertex/VertexTools/interface/VertexDistance3D.h"

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
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/MultiTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/MultiTrackPointingKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/CombinedKinematicConstraint.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"

// Needed for IP3D
#include "RecoTauTag/ImpactParameter/interface/ImpactParameterAlgorithm.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include <RecoBTag/BTagTools/interface/SignedImpactParameter3D.h>

#include <iostream>

#include "VtxUtils.hh"
#include "TMatrixDSym.h"
#include "TVectorD.h"

using namespace std;
using namespace vtxu;

#define _MuMass_ 0.1056583745
#define _MuMassErr_ 0.0000000024
#define _PiMass_ 0.13957018
#define _PiMassErr_ 0.00000035
#define _KMass_ 0.493677
#define _KMassErr_ 0.000013
#define _KstMass_ 0.89166
#define _KstMassErr_ 0.00011
#define _PhiMass_ 1.020
#define _D0Mass_ 1.86483
#define _D0MassErr_ 0.00017
#define _DstMass_ 2.01026
#define _DstMassErr_ 0.00017
#define _JpsiMass_ 3.096916
#define _JpsiMassErr_ 0.00001
#define _B0Mass_ 5.27963
#define _B0MassErr_ 0.00026

static int isMC = 0;

/* Returns a new vertex fit using the adaptive vertex fitter with tracks after
 * we fix the covariance matrix.
 *
 * See https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideOfflinePrimaryVertexProduction#Refitting_vertices_with_selected. */
reco::Vertex vtxu::refit_vertex(edm::Event& iEvent, const edm::EventSetup& iSetup, size_t ipv, bool beamSpotContrant, const std::vector<pat::PackedCandidate> &pfCandHandle)
{
    unsigned int i;
    std::vector<reco::TransientTrack> mytracks;

    edm::ESHandle<TransientTrackBuilder> TTBuilder;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",TTBuilder);

    edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
    iEvent.getByLabel("offlineBeamSpot", recoBeamSpotHandle);
    reco::BeamSpot vertexBeamSpot = *recoBeamSpotHandle;

    // edm::Handle<vector<reco::Vertex>> vtxHandle;
    // iEvent.getByLabel("offlineSlimmedPrimaryVertices", vtxHandle);

    // cout << "========> Vtx: " << ipv << endl;
    for (i = 0; i < pfCandHandle.size(); i++) {
        const pat::PackedCandidate &ptk = pfCandHandle[i];

        if (!ptk.hasTrackDetails()) continue;

        if (ptk.fromPV(ipv) < 3) continue;

        // Equivalent to the below
        // if (ptk.pvAssociationQuality() != pat::PackedCandidate::PVAssociationQuality::UsedInFitTight) continue;
        // auto vtxTK = ptk.vertexRef();
        // auto vtxMINI = (*vtxHandle)[ipv];
        // if (vtxTK->x() != vtxMINI.x() || vtxTK->y() != vtxMINI.y() || vtxTK->z() != vtxMINI.z()) continue;

        auto tk = ptk.bestTrack();
        reco::TransientTrack transientTrack = TTBuilder->build(fix_track(tk));
        transientTrack.setBeamSpot(vertexBeamSpot);
        mytracks.push_back(transientTrack);
    }

    if (mytracks.size() < 2) {
        // cout << "[WARNING]: Less than 2 tracks for vertex "<<ipv<<" fit" << endl;
        return reco::Vertex(); // Invalid vertex
    }

    AdaptiveVertexFitter theFitter;
    if (beamSpotContrant) {
      TransientVertex tmp = theFitter.vertex(mytracks, vertexBeamSpot);
      return tmp;
    }
    else {
      TransientVertex tmp = theFitter.vertex(mytracks);
      return tmp;
    }

}

void vtxu::set_isMC(int _isMC)
{
    isMC = _isMC;
}

/* Returns the impact parameter uncertainty with respect to a given vertex. */
double vtxu::dxyError(const reco::TrackBase &tk, const reco::Vertex &vtx) {
    return dxyError(tk, vtx.position(), vtx.covariance());
}

/* Returns the impact parameter uncertainty with respect to a given point and
 * covariance matrix. Taken from a newer version of cmssw.
 *
 * See https://github.com/cms-sw/cmssw/blob/29f5fc15b34591745c5cd3c2c6eb9793aa6f371b/DataFormats/TrackReco/src/TrackBase.cc#L148. */
double vtxu::dxyError(const reco::TrackBase &tk, const reco::TrackBase::Point &vtx, const math::Error<3>::type &vertexCov) {
  // Gradient of TrackBase::dxy(const Point &myBeamSpot) with respect to track parameters. Using unrolled expressions to avoid calling for higher dimension matrices
  // ( 0, 0, x_vert * cos(phi) + y_vert * sin(phi), 1, 0 )
  // Gradient with respect to point parameters
  // ( sin(phi), -cos(phi))
  // Propagate covariance assuming cross-terms of the covariance between track and vertex parameters are 0
  return std::sqrt((vtx.x() * tk.px() + vtx.y() * tk.py()) * (vtx.x() * tk.px() + vtx.y() * tk.py()) / (tk.pt() * tk.pt()) *
                    tk.covariance()(reco::TrackBase::i_phi, reco::TrackBase::i_phi) +
                2 * (vtx.x() * tk.px() + vtx.y() * tk.py()) / tk.pt() * tk.covariance()(reco::TrackBase::i_phi, reco::TrackBase::i_dxy) + tk.covariance()(reco::TrackBase::i_dxy, reco::TrackBase::i_dxy) +
                tk.py() * tk.py() / (tk.pt() * tk.pt()) * vertexCov(0, 0) - 2 * tk.py() * tk.px() / (tk.pt() * tk.pt()) * vertexCov(0, 1) +
                tk.px() * tk.px() / (tk.pt() * tk.pt()) * vertexCov(1, 1));
}

/* Returns the impact parameter uncertainty with respect to a given beamspot.
 * Also taken from a newer version of cmssw. */
double vtxu::dxyError(const reco::TrackBase &tk, const reco::BeamSpot &theBeamSpot) {
    return dxyError(tk, theBeamSpot.position(tk.vz()), theBeamSpot.rotatedCovariance3D());
}

reco::Track vtxu::fix_track(const reco::TrackRef& tk)
{
    reco::Track t = reco::Track(*tk);
    return fix_track(&t);
}

/* Check for a not positive definite covariance matrix. If the covariance
 * matrix is not positive definite, we force it to be positive definite by
 * adding the minimum eigenvalue to the diagonal of the covariance matrix plus
 * `delta`.
 *
 * See https://nhigham.com/2020/12/22/what-is-a-modified-cholesky-factorization/.
 *
 * Note: There may be better ways of doing this, but since most of these tracks
 * don't end up in the final sample, this is probably good enough. */
reco::Track vtxu::fix_track(const reco::Track *tk, double delta)
{
    unsigned int i, j;
    double min_eig = 1;

    /* Get the original covariance matrix. */
    reco::TrackBase::CovarianceMatrix cov = tk->covariance();

    /* Convert it from an SMatrix to a TMatrixD so we can get the eigenvalues. */
    TMatrixDSym new_cov(cov.kRows);
    for (i = 0; i < cov.kRows; i++) {
        for (j = 0; j < cov.kRows; j++) {
            /* Need to check for nan or inf, because for some reason these
             * cause a segfault when calling Eigenvectors().
             *
             * No idea what to do here or why this happens. */
            if (std::isnan(cov(i,j)) || std::isinf(cov(i,j)))
                cov(i,j) = 1e-6;
            new_cov(i,j) = cov(i,j);
        }
    }

    /* Get the eigenvalues. */
    TVectorD eig(cov.kRows);
    new_cov.EigenVectors(eig);
    for (i = 0; i < cov.kRows; i++)
        if (eig(i) < min_eig)
            min_eig = eig(i);

    /* If the minimum eigenvalue is less than zero, then subtract it from the
     * diagonal and add `delta`. */
    if (min_eig < 0) {
        for (i = 0; i < cov.kRows; i++)
            cov(i,i) -= min_eig - delta;
    }

    /* Correct for the fact that the impact parameter error is significantly
     * different in data and MC. Technically we should be *decreasing* the data
     * uncertainty since according to Marco Musich:
     *
     * > we suspect that the track uncertainty in data is over-estimated, while
     * > it should be well calibrated in MC.  To corroborate this, here is a
     * > plot from the recent alignment paper which shows the pull of the pixel
     * > residuals for various reconstruction passes and MC:
     * > for an ideally calibrated detector you expect a peak around 1, but for
     * > the "alignment during data-taking" which is what was used for the
     * > BParking Prompt reconstruction you can see it's pretty shifted to
     * > lower values, meaning the error is overestimated.
     * > So I would recommend to scale down the data uncertainties.
     *
     * However, we can't do this because this would cause the impact parameter
     * significance to increase which means that some events which didn't
     * trigger should have triggered. Therefore, we *increase* the error on the
     * MC.
     *
     * This number comes from looking at the impact parameter error for the
     * trigger muon for the B0 -> J/psi K* events. I tested both a scaling and
     * an offset and the offset seemed to match really well, whereas the
     * scaling didn't work. */
    if (isMC)
        cov(reco::TrackBase::i_dxy, reco::TrackBase::i_dxy) = pow(sqrt(cov(reco::TrackBase::i_dxy, reco::TrackBase::i_dxy)) + 0.0003,2);

    return reco::Track(tk->chi2(), tk->ndof(), tk->referencePoint(), tk->momentum(), tk->charge(), cov, tk->algo(), (reco::TrackBase::TrackQuality) tk->qualityMask());
}

RefCountedKinematicTree vtxu::FitD0(const edm::EventSetup& iSetup, pat::PackedCandidate pi, pat::PackedCandidate K, bool mass_constraint) {
  // Get transient track builder
  edm::ESHandle<TransientTrackBuilder> TTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",TTBuilder);

  reco::TransientTrack pi_tk = TTBuilder->build(fix_track(pi.bestTrack()));
  reco::TransientTrack K_tk = TTBuilder->build(fix_track(K.bestTrack()));

  KinematicParticleFactoryFromTransientTrack pFactory;

  std::vector<RefCountedKinematicParticle> parts;
  double chi = 0, ndf = 0;
  float mK = _KMass_, dmK = _KMassErr_;
  float mPi = _PiMass_, dmPi = _PiMassErr_;
  parts.push_back(pFactory.particle(K_tk, mK, chi, ndf, dmK));
  parts.push_back(pFactory.particle(pi_tk, mPi, chi, ndf, dmPi));

  if (!mass_constraint) {
    KinematicParticleVertexFitter VtxFitter;
    RefCountedKinematicTree D0KinTree = VtxFitter.fit(parts);
    return D0KinTree;
  }
  else {
    ParticleMass D0mass = _D0Mass_;
    MultiTrackKinematicConstraint * D0mass_c = new TwoTrackMassKinematicConstraint(D0mass);
    KinematicConstrainedVertexFitter kcVtxFitter;
    RefCountedKinematicTree D0KinTree = kcVtxFitter.fit(parts, D0mass_c);
    return D0KinTree;
  }
}


RefCountedKinematicTree vtxu::FitKst_piK(const edm::EventSetup& iSetup, pat::PackedCandidate pi, pat::PackedCandidate K, bool mass_constraint) {
  // Get transient track builder
  edm::ESHandle<TransientTrackBuilder> TTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",TTBuilder);

  reco::TransientTrack pi_tk = TTBuilder->build(fix_track(pi.bestTrack()));
  reco::TransientTrack K_tk = TTBuilder->build(fix_track(K.bestTrack()));

  KinematicParticleFactoryFromTransientTrack pFactory;

  std::vector<RefCountedKinematicParticle> parts;
  double chi = 0, ndf = 0;
  float mK = _KMass_, dmK = _KMassErr_;
  float mPi = _PiMass_, dmPi = _PiMassErr_;
  parts.push_back(pFactory.particle(pi_tk, mPi, chi, ndf, dmPi));
  parts.push_back(pFactory.particle(K_tk, mK, chi, ndf, dmK));

  if (!mass_constraint) {
    KinematicParticleVertexFitter VtxFitter;
    RefCountedKinematicTree kinTree = VtxFitter.fit(parts);
    return kinTree;
  }
  else {
    ParticleMass mass = _KstMass_;
    MultiTrackKinematicConstraint * mass_c = new TwoTrackMassKinematicConstraint(mass);
    KinematicConstrainedVertexFitter kcVtxFitter;
    RefCountedKinematicTree kinTree = kcVtxFitter.fit(parts, mass_c);
    return kinTree;
  }
}

RefCountedKinematicTree vtxu::FitPhi_KK(const edm::EventSetup& iSetup, pat::PackedCandidate K1, pat::PackedCandidate K2, bool mass_constraint) {
  // Get transient track builder
  edm::ESHandle<TransientTrackBuilder> TTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",TTBuilder);

  reco::TransientTrack K1_tk = TTBuilder->build(fix_track(K1.bestTrack()));
  reco::TransientTrack K2_tk = TTBuilder->build(fix_track(K2.bestTrack()));

  KinematicParticleFactoryFromTransientTrack pFactory;

  std::vector<RefCountedKinematicParticle> parts;
  double chi = 0, ndf = 0;
  float mK = _KMass_, dmK = _KMassErr_;
  parts.push_back(pFactory.particle(K1_tk, mK, chi, ndf, dmK));
  parts.push_back(pFactory.particle(K2_tk, mK, chi, ndf, dmK));

  if (!mass_constraint) {
    KinematicParticleVertexFitter VtxFitter;
    RefCountedKinematicTree kinTree = VtxFitter.fit(parts);
    return kinTree;
  }
  else {
    ParticleMass mass = _PhiMass_;
    MultiTrackKinematicConstraint * mass_c = new TwoTrackMassKinematicConstraint(mass);
    KinematicConstrainedVertexFitter kcVtxFitter;
    RefCountedKinematicTree kinTree = kcVtxFitter.fit(parts, mass_c);
    return kinTree;
  }
}


RefCountedKinematicTree vtxu::FitDst_fitD0wMassConstraint(const edm::EventSetup& iSetup, pat::PackedCandidate pisoft, pat::PackedCandidate pi, pat::PackedCandidate K, bool mass_constraint, int verbose = 0) {
  // Get the mass constrained D0
  auto D0KinTree = vtxu::FitD0(iSetup, pi, K, true);
  if(!D0KinTree->isValid()) return D0KinTree;
  D0KinTree->movePointerToTheTop();
  auto D0_reco = D0KinTree->currentParticle();

  // Get transient track builder
  edm::ESHandle<TransientTrackBuilder> TTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",TTBuilder);

  reco::TransientTrack pisoft_tk = TTBuilder->build(fix_track(pisoft.bestTrack()));

  KinematicParticleFactoryFromTransientTrack pFactory;

  std::vector<RefCountedKinematicParticle> parts;
  double chi = 0, ndf = 0;
  float mPi = _PiMass_, dmPi = _PiMassErr_;
  parts.push_back(D0_reco);
  parts.push_back(pFactory.particle(pisoft_tk, mPi, chi, ndf, dmPi));

  if (!mass_constraint) {
    KinematicParticleVertexFitter VtxFitter;
    RefCountedKinematicTree DstKinTree = VtxFitter.fit(parts);
    return DstKinTree;
  }
  else {
    ParticleMass Dstmass = _DstMass_;
    MultiTrackKinematicConstraint * Dstmass_c = new TwoTrackMassKinematicConstraint(Dstmass);
    KinematicConstrainedVertexFitter kcVtxFitter;
    RefCountedKinematicTree DstKinTree = kcVtxFitter.fit(parts, Dstmass_c);
    return DstKinTree;
  }
}

RefCountedKinematicTree vtxu::FitJpsi_mumu(const edm::EventSetup& iSetup, pat::Muon m1, pat::Muon m2, bool mass_constraint) {
  // Get transient track builder
  edm::ESHandle<TransientTrackBuilder> TTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",TTBuilder);

  reco::TransientTrack m1_tk = TTBuilder->build(fix_track(m1.muonBestTrack()));
  reco::TransientTrack m2_tk = TTBuilder->build(fix_track(m2.muonBestTrack()));

  std::vector<RefCountedKinematicParticle> parts;
  KinematicParticleFactoryFromTransientTrack pFactory;
  double chi = 0, ndf = 0;
  float mMu = _MuMass_, dmMu = _MuMassErr_;
  parts.push_back(pFactory.particle(m1_tk, mMu, chi, ndf, dmMu));
  parts.push_back(pFactory.particle(m2_tk, mMu, chi, ndf, dmMu));

  if (!mass_constraint) {
    KinematicParticleVertexFitter VtxFitter;
    RefCountedKinematicTree KinTree = VtxFitter.fit(parts);
    return KinTree;
  }
  else {
    ParticleMass mass = _JpsiMass_;
    MultiTrackKinematicConstraint * mass_c = new TwoTrackMassKinematicConstraint(mass);
    KinematicConstrainedVertexFitter kcVtxFitter;
    RefCountedKinematicTree KinTree = kcVtxFitter.fit(parts, mass_c);
    return KinTree;
  }
}

RefCountedKinematicTree vtxu::FitB_mumuK(const edm::EventSetup& iSetup, pat::Muon m1, pat::Muon m2, pat::PackedCandidate K) {
  // Get transient track builder
  edm::ESHandle<TransientTrackBuilder> TTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",TTBuilder);

  reco::TransientTrack m1_tk = TTBuilder->build(fix_track(m1.muonBestTrack()));
  reco::TransientTrack m2_tk = TTBuilder->build(fix_track(m2.muonBestTrack()));
  reco::TransientTrack K_tk = TTBuilder->build(fix_track(K.bestTrack()));

  std::vector<RefCountedKinematicParticle> parts;
  KinematicParticleFactoryFromTransientTrack pFactory;
  double chi = 0, ndf = 0;
  float mMu = _MuMass_, dmMu = _MuMassErr_;
  float mK = _KMass_, dmK = _KMassErr_;
  parts.push_back(pFactory.particle(m1_tk, mMu, chi, ndf, dmMu));
  parts.push_back(pFactory.particle(m2_tk, mMu, chi, ndf, dmMu));
  parts.push_back(pFactory.particle(K_tk, mK, chi, ndf, dmK));

  KinematicParticleVertexFitter VtxFitter;
  RefCountedKinematicTree KinTree = VtxFitter.fit(parts);
  return KinTree;
}

RefCountedKinematicTree vtxu::FitB_mumupiK(const edm::EventSetup& iSetup, pat::Muon m1, pat::Muon m2, pat::PackedCandidate pi, pat::PackedCandidate K, bool JpsiKstmass_constraint, bool Bmass_constraint, bool pointing_constraint, reco::Vertex* ptrPV) {
  // Get transient track builder
  edm::ESHandle<TransientTrackBuilder> TTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",TTBuilder);

  reco::TransientTrack m1_tk = TTBuilder->build(fix_track(m1.muonBestTrack()));
  reco::TransientTrack m2_tk = TTBuilder->build(fix_track(m2.muonBestTrack()));
  reco::TransientTrack pi_tk = TTBuilder->build(fix_track(pi.bestTrack()));
  reco::TransientTrack K_tk = TTBuilder->build(fix_track(K.bestTrack()));

  std::vector<RefCountedKinematicParticle> parts;
  KinematicParticleFactoryFromTransientTrack pFactory;
  double chi = 0, ndf = 0;
  float mMu = _MuMass_, dmMu = _MuMassErr_;
  float mK = _KMass_, dmK = _KMassErr_;
  float mPi = _PiMass_, dmPi = _PiMassErr_;
  parts.push_back(pFactory.particle(m1_tk, mMu, chi, ndf, dmMu));
  parts.push_back(pFactory.particle(m2_tk, mMu, chi, ndf, dmMu));
  parts.push_back(pFactory.particle(pi_tk, mPi, chi, ndf, dmPi));
  parts.push_back(pFactory.particle(K_tk, mK, chi, ndf, dmK));

  if (!Bmass_constraint && !JpsiKstmass_constraint && !pointing_constraint) {
    KinematicParticleVertexFitter VtxFitter;
    RefCountedKinematicTree KinTree = VtxFitter.fit(parts);
    return KinTree;
  }
  else {
    vector<MultiTrackKinematicConstraint* > kinConstraints;

    if(JpsiKstmass_constraint) {
      ParticleMass mass_Jpsi = _JpsiMass_;
      MultiTrackKinematicConstraint * c_mass_Jpsi = new MultiTrackMassKinematicConstraint(mass_Jpsi, 2);
      kinConstraints.push_back(c_mass_Jpsi);

      ParticleMass mass_Kst = _KstMass_;
      MultiTrackKinematicConstraint * c_mass_Kst = new MultiTrackMassKinematicConstraint(mass_Kst, 2);
      kinConstraints.push_back(c_mass_Kst);
    }

    if(Bmass_constraint){
      ParticleMass mass_B = _B0Mass_;
      MultiTrackKinematicConstraint * c_mass_B = new MultiTrackMassKinematicConstraint(mass_B, 4);
      kinConstraints.push_back(c_mass_B);
    }

    if(pointing_constraint) {
      if (ptrPV == nullptr) {
        cerr << "[ERROR]: Pointing constraint add but no vertex provided" << endl;
      }
      else {
        GlobalPoint vtxPosition(ptrPV->x(), ptrPV->y(), ptrPV->z());
        MultiTrackKinematicConstraint * c_pointing_PV = new MultiTrackPointingKinematicConstraint (vtxPosition);
        kinConstraints.push_back(c_pointing_PV);
      }
    }

    MultiTrackKinematicConstraint * combConstraints = new CombinedKinematicConstraint(kinConstraints);

    KinematicConstrainedVertexFitter kcVtxFitter;
    RefCountedKinematicTree KinTree = kcVtxFitter.fit(parts, combConstraints);
    return KinTree;
  }
}

RefCountedKinematicTree vtxu::FitB_mumupiK(const RefCountedKinematicParticle m1, const RefCountedKinematicParticle m2, const RefCountedKinematicParticle pi, const RefCountedKinematicParticle K, bool JpsiKstmass_constraint, bool Bmass_constraint, bool pointing_constraint, reco::Vertex* ptrPV) {

  reco::TransientTrack m1_tk = m1->refittedTransientTrack();
  reco::TransientTrack m2_tk = m2->refittedTransientTrack();
  reco::TransientTrack pi_tk = pi->refittedTransientTrack();
  reco::TransientTrack K_tk = K->refittedTransientTrack();

  std::vector<RefCountedKinematicParticle> parts;
  KinematicParticleFactoryFromTransientTrack pFactory;
  double chi = 0, ndf = 0;
  float mMu = _MuMass_, dmMu = _MuMassErr_;
  float mK = _KMass_, dmK = _KMassErr_;
  float mPi = _PiMass_, dmPi = _PiMassErr_;
  parts.push_back(pFactory.particle(m1_tk, mMu, chi, ndf, dmMu));
  parts.push_back(pFactory.particle(m2_tk, mMu, chi, ndf, dmMu));
  parts.push_back(pFactory.particle(K_tk, mK, chi, ndf, dmK));
  parts.push_back(pFactory.particle(pi_tk, mPi, chi, ndf, dmPi));

  if (!Bmass_constraint && !JpsiKstmass_constraint && !pointing_constraint) {
    KinematicParticleVertexFitter VtxFitter;
    RefCountedKinematicTree KinTree = VtxFitter.fit(parts);
    return KinTree;
  }
  else {
    vector<MultiTrackKinematicConstraint* > kinConstraints;

    if(JpsiKstmass_constraint) {
      ParticleMass mass_Jpsi = _JpsiMass_;
      MultiTrackKinematicConstraint * c_mass_Jpsi = new MultiTrackMassKinematicConstraint(mass_Jpsi, 2);
      kinConstraints.push_back(c_mass_Jpsi);

      ParticleMass mass_Kst = _KstMass_;
      MultiTrackKinematicConstraint * c_mass_Kst = new MultiTrackMassKinematicConstraint(mass_Kst, 2);
      kinConstraints.push_back(c_mass_Kst);
    }

    if(Bmass_constraint){
      ParticleMass mass_B = _B0Mass_;
      MultiTrackKinematicConstraint * c_mass_B = new MultiTrackMassKinematicConstraint(mass_B, 4);
      kinConstraints.push_back(c_mass_B);
    }

    if(pointing_constraint) {
      if (ptrPV == nullptr) {
        cerr << "[ERROR]: Pointing constraint add but no vertex provided" << endl;
      }
      else {
        GlobalPoint vtxPosition(ptrPV->x(), ptrPV->y(), ptrPV->z());
        MultiTrackKinematicConstraint * c_pointing_PV = new MultiTrackPointingKinematicConstraint (vtxPosition);
        kinConstraints.push_back(c_pointing_PV);
      }
    }

    MultiTrackKinematicConstraint * combConstraints = new CombinedKinematicConstraint(kinConstraints);

    KinematicConstrainedVertexFitter kcVtxFitter;
    RefCountedKinematicTree KinTree = kcVtxFitter.fit(parts, combConstraints);
    return KinTree;
  }
}

RefCountedKinematicTree vtxu::FitB_mumupiK_tk(const edm::EventSetup& iSetup, pat::Muon m1, pat::Muon m2, pat::PackedCandidate pi, pat::PackedCandidate K, pat::PackedCandidate addTk, bool JpsiKstmass_constraint, bool Bmass_constraint, bool pointing_constraint, reco::Vertex* ptrPV) {
  // Get transient track builder
  edm::ESHandle<TransientTrackBuilder> TTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",TTBuilder);

  reco::TransientTrack m1_tk = TTBuilder->build(fix_track(m1.muonBestTrack()));
  reco::TransientTrack m2_tk = TTBuilder->build(fix_track(m2.muonBestTrack()));
  reco::TransientTrack pi_tk = TTBuilder->build(fix_track(pi.bestTrack()));
  reco::TransientTrack K_tk = TTBuilder->build(fix_track(K.bestTrack()));
  reco::TransientTrack tk_tk = TTBuilder->build(fix_track(addTk.bestTrack()));

  std::vector<RefCountedKinematicParticle> parts;
  KinematicParticleFactoryFromTransientTrack pFactory;
  double chi = 0, ndf = 0;
  float mMu = _MuMass_, dmMu = _MuMassErr_;
  float mK = _KMass_, dmK = _KMassErr_;
  float mPi = _PiMass_, dmPi = _PiMassErr_;
  parts.push_back(pFactory.particle(m1_tk, mMu, chi, ndf, dmMu));
  parts.push_back(pFactory.particle(m2_tk, mMu, chi, ndf, dmMu));
  parts.push_back(pFactory.particle(pi_tk, mPi, chi, ndf, dmPi));
  parts.push_back(pFactory.particle(K_tk, mK, chi, ndf, dmK));
  parts.push_back(pFactory.particle(tk_tk, mPi, chi, ndf, dmPi));

  if (!Bmass_constraint && !JpsiKstmass_constraint && !pointing_constraint) {
    KinematicParticleVertexFitter VtxFitter;
    RefCountedKinematicTree KinTree = VtxFitter.fit(parts);
    return KinTree;
  }
  else {
    vector<MultiTrackKinematicConstraint* > kinConstraints;

    if(JpsiKstmass_constraint) {
      ParticleMass mass_Jpsi = _JpsiMass_;
      MultiTrackKinematicConstraint * c_mass_Jpsi = new MultiTrackMassKinematicConstraint(mass_Jpsi, 2);
      kinConstraints.push_back(c_mass_Jpsi);

      ParticleMass mass_Kst = _KstMass_;
      MultiTrackKinematicConstraint * c_mass_Kst = new MultiTrackMassKinematicConstraint(mass_Kst, 2);
      kinConstraints.push_back(c_mass_Kst);
    }

    if(Bmass_constraint){
      ParticleMass mass_B = _B0Mass_;
      MultiTrackKinematicConstraint * c_mass_B = new MultiTrackMassKinematicConstraint(mass_B, 4);
      kinConstraints.push_back(c_mass_B);
    }

    if(pointing_constraint) {
      if (ptrPV == nullptr) {
        cerr << "[ERROR]: Pointing constraint add but no vertex provided" << endl;
      }
      else {
        GlobalPoint vtxPosition(ptrPV->x(), ptrPV->y(), ptrPV->z());
        MultiTrackKinematicConstraint * c_pointing_PV = new MultiTrackPointingKinematicConstraint (vtxPosition);
        kinConstraints.push_back(c_pointing_PV);
      }
    }

    MultiTrackKinematicConstraint * combConstraints = new CombinedKinematicConstraint(kinConstraints);

    KinematicConstrainedVertexFitter kcVtxFitter;
    RefCountedKinematicTree KinTree = kcVtxFitter.fit(parts, combConstraints);
    return KinTree;
  }
}

RefCountedKinematicTree vtxu::FitDst(const edm::EventSetup& iSetup, pat::PackedCandidate pisoft, const RefCountedKinematicParticle D0, bool mass_constraint) {
// Get transient track builder
  edm::ESHandle<TransientTrackBuilder> TTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",TTBuilder);

  reco::TransientTrack pisoft_tk = TTBuilder->build(fix_track(pisoft.bestTrack()));
  reco::TransientTrack D0_tk = D0->refittedTransientTrack();

  KinematicParticleFactoryFromTransientTrack pFactory;

  std::vector<RefCountedKinematicParticle> parts;
  double chi = 0, ndf = 0;
  float mPi = _PiMass_, dmPi = _PiMassErr_;
  parts.push_back(pFactory.particle(pisoft_tk, mPi, chi, ndf, dmPi));
  float mD0 = D0->currentState().mass();
  float dmD0 = sqrt(D0->currentState().kinematicParametersError().matrix()(6,6));
  parts.push_back(pFactory.particle(D0_tk, mD0, chi, ndf, dmD0));
  if (!mass_constraint) {
    KinematicParticleVertexFitter VtxFitter;
    RefCountedKinematicTree DstKinTree = VtxFitter.fit(parts);
    return DstKinTree;
  }
  else {
    ParticleMass Dstmass = _DstMass_;
    MultiTrackKinematicConstraint * Dstmass_c = new TwoTrackMassKinematicConstraint(Dstmass);
    KinematicConstrainedVertexFitter kcVtxFitter;
    RefCountedKinematicTree DstKinTree = kcVtxFitter.fit(parts, Dstmass_c);
    return DstKinTree;
  }
}

RefCountedKinematicTree vtxu::FitD0pipi(const edm::EventSetup& iSetup, const RefCountedKinematicParticle D0, pat::PackedCandidate pi1, pat::PackedCandidate pi2) {
// Get transient track builder
  edm::ESHandle<TransientTrackBuilder> TTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",TTBuilder);

  reco::TransientTrack D0_tk = D0->refittedTransientTrack();
  reco::TransientTrack pi1_tk = TTBuilder->build(fix_track(pi1.bestTrack()));
  reco::TransientTrack pi2_tk = TTBuilder->build(fix_track(pi2.bestTrack()));

  KinematicParticleFactoryFromTransientTrack pFactory;

  std::vector<RefCountedKinematicParticle> parts;
  double chi = 0, ndf = 0;
  float mD0 = D0->currentState().mass();
  float dmD0 = sqrt(D0->currentState().kinematicParametersError().matrix()(6,6));
  parts.push_back(pFactory.particle(D0_tk, mD0, chi, ndf, dmD0));
  float mPi = _PiMass_, dmPi = _PiMassErr_;
  parts.push_back(pFactory.particle(pi1_tk, mPi, chi, ndf, dmPi));
  parts.push_back(pFactory.particle(pi2_tk, mPi, chi, ndf, dmPi));

  KinematicParticleVertexFitter VtxFitter;
  RefCountedKinematicTree kinTree = VtxFitter.fit(parts);
  return kinTree;
}

RefCountedKinematicTree vtxu::FitD0pipipis(const edm::EventSetup& iSetup, const RefCountedKinematicParticle D0, pat::PackedCandidate pi1, pat::PackedCandidate pi2, pat::PackedCandidate pis) {
// Get transient track builder
  edm::ESHandle<TransientTrackBuilder> TTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",TTBuilder);

  reco::TransientTrack D0_tk = D0->refittedTransientTrack();
  reco::TransientTrack pi1_tk = TTBuilder->build(fix_track(pi1.bestTrack()));
  reco::TransientTrack pi2_tk = TTBuilder->build(fix_track(pi2.bestTrack()));
  reco::TransientTrack pis_tk = TTBuilder->build(fix_track(pis.bestTrack()));

  KinematicParticleFactoryFromTransientTrack pFactory;

  std::vector<RefCountedKinematicParticle> parts;
  double chi = 0, ndf = 0;
  float mD0 = D0->currentState().mass();
  float dmD0 = sqrt(D0->currentState().kinematicParametersError().matrix()(6,6));
  parts.push_back(pFactory.particle(D0_tk, mD0, chi, ndf, dmD0));
  float mPi = _PiMass_, dmPi = _PiMassErr_;
  parts.push_back(pFactory.particle(pi1_tk, mPi, chi, ndf, dmPi));
  parts.push_back(pFactory.particle(pi2_tk, mPi, chi, ndf, dmPi));
  parts.push_back(pFactory.particle(pis_tk, mPi, chi, ndf, dmPi));

  KinematicParticleVertexFitter VtxFitter;
  RefCountedKinematicTree kinTree = VtxFitter.fit(parts);
  return kinTree;
}

RefCountedKinematicTree vtxu::FitD0_pions(const edm::EventSetup& iSetup, const RefCountedKinematicParticle D0, vector<pat::PackedCandidate> pions) {
// Get transient track builder
  edm::ESHandle<TransientTrackBuilder> TTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",TTBuilder);

  reco::TransientTrack D0_tk = D0->refittedTransientTrack();

  KinematicParticleFactoryFromTransientTrack pFactory;
  std::vector<RefCountedKinematicParticle> parts;
  double chi = 0, ndf = 0;
  float mD0 = D0->currentState().mass();
  float dmD0 = sqrt(D0->currentState().kinematicParametersError().matrix()(6,6));
  parts.push_back(pFactory.particle(D0_tk, mD0, chi, ndf, dmD0));
  float mPi = _PiMass_, dmPi = _PiMassErr_;

  for(auto pi : pions) {
    reco::TransientTrack pi_tk = TTBuilder->build(fix_track(pi.bestTrack()));
    parts.push_back(pFactory.particle(pi_tk, mPi, chi, ndf, dmPi));
  }

  KinematicParticleVertexFitter VtxFitter;
  RefCountedKinematicTree kinTree = VtxFitter.fit(parts);
  return kinTree;
}

RefCountedKinematicTree vtxu::FitB_D0pismu(const edm::EventSetup& iSetup, const RefCountedKinematicParticle D0, pat::PackedCandidate pis, pat::Muon mu){
  edm::ESHandle<TransientTrackBuilder> TTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",TTBuilder);

  reco::TransientTrack D0_tk = D0->refittedTransientTrack();
  reco::TransientTrack pis_tk = TTBuilder->build(fix_track(pis.bestTrack()));
  reco::TransientTrack mu_tk = TTBuilder->build(fix_track(mu.muonBestTrack()));

  KinematicParticleFactoryFromTransientTrack pFactory;
  std::vector<RefCountedKinematicParticle> parts;
  double chi = 0, ndf = 0;

  float mD0 = D0->currentState().mass();
  float dmD0 = sqrt(D0->currentState().kinematicParametersError().matrix()(6,6));
  parts.push_back(pFactory.particle(D0_tk, mD0, chi, ndf, dmD0));
  float mPi = _PiMass_, dmPi = _PiMassErr_;
  parts.push_back(pFactory.particle(pis_tk, mPi, chi, ndf, dmPi));
  float mMu = _MuMass_, dmMu = _MuMassErr_;
  parts.push_back(pFactory.particle(mu_tk, mMu, chi, ndf, dmMu));

  KinematicParticleVertexFitter VtxFitter;
  RefCountedKinematicTree BKinTree = VtxFitter.fit(parts);
  return BKinTree;
}

RefCountedKinematicTree vtxu::Fit_D0pismupi(const edm::EventSetup& iSetup, const RefCountedKinematicParticle D0, pat::PackedCandidate pis, pat::Muon mu, pat::PackedCandidate pi){
  edm::ESHandle<TransientTrackBuilder> TTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",TTBuilder);

  reco::TransientTrack D0_tk = D0->refittedTransientTrack();
  reco::TransientTrack pis_tk = TTBuilder->build(fix_track(pis.bestTrack()));
  reco::TransientTrack mu_tk = TTBuilder->build(fix_track(mu.muonBestTrack()));
  reco::TransientTrack pi_tk = TTBuilder->build(fix_track(pi.bestTrack()));

  KinematicParticleFactoryFromTransientTrack pFactory;
  std::vector<RefCountedKinematicParticle> parts;
  double chi = 0, ndf = 0;

  float mD0 = D0->currentState().mass();
  float dmD0 = sqrt(D0->currentState().kinematicParametersError().matrix()(6,6));
  parts.push_back(pFactory.particle(D0_tk, mD0, chi, ndf, dmD0));
  float mPi = _PiMass_, dmPi = _PiMassErr_;
  parts.push_back(pFactory.particle(pis_tk, mPi, chi, ndf, dmPi));
  float mMu = _MuMass_, dmMu = _MuMassErr_;
  parts.push_back(pFactory.particle(mu_tk, mMu, chi, ndf, dmMu));
  parts.push_back(pFactory.particle(pi_tk, mPi, chi, ndf, dmPi));

  KinematicParticleVertexFitter VtxFitter;
  RefCountedKinematicTree BKinTree = VtxFitter.fit(parts);
  return BKinTree;
}

RefCountedKinematicTree vtxu::Fit_D0pihpis(const edm::EventSetup& iSetup, const RefCountedKinematicParticle D0, pat::PackedCandidate pih, pat::PackedCandidate pis){
  edm::ESHandle<TransientTrackBuilder> TTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",TTBuilder);

  reco::TransientTrack D0_tk = D0->refittedTransientTrack();
  reco::TransientTrack pih_tk = TTBuilder->build(fix_track(pih.bestTrack()));
  reco::TransientTrack pis_tk = TTBuilder->build(fix_track(pis.bestTrack()));

  KinematicParticleFactoryFromTransientTrack pFactory;
  std::vector<RefCountedKinematicParticle> parts;
  double chi = 0, ndf = 0;

  float mD0 = D0->currentState().mass();
  float dmD0 = sqrt(D0->currentState().kinematicParametersError().matrix()(6,6));
  parts.push_back(pFactory.particle(D0_tk, mD0, chi, ndf, dmD0));
  float mPi = _PiMass_, dmPi = _PiMassErr_;
  parts.push_back(pFactory.particle(pih_tk, mPi, chi, ndf, dmPi));
  parts.push_back(pFactory.particle(pis_tk, mPi, chi, ndf, dmPi));

  KinematicParticleVertexFitter VtxFitter;
  RefCountedKinematicTree BKinTree = VtxFitter.fit(parts);
  return BKinTree;
}


RefCountedKinematicTree vtxu::FitVtxMuDst(const edm::EventSetup& iSetup, const RefCountedKinematicParticle Dst, pat::Muon mu) {
  // Get transient track builder
  edm::ESHandle<TransientTrackBuilder> TTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",TTBuilder);

  reco::TransientTrack Mu_tk = TTBuilder->build(fix_track(mu.muonBestTrack()));
  reco::TransientTrack Dst_tk = Dst->refittedTransientTrack();

  KinematicParticleFactoryFromTransientTrack pFactory;
  std::vector<RefCountedKinematicParticle> parts;
  double chi = 0, ndf = 0;
  float mMu = _MuMass_, dmMu = _MuMassErr_;
  parts.push_back(pFactory.particle(Mu_tk, mMu, chi, ndf, dmMu));
  float mDst = Dst->currentState().mass();
  float dmDst = sqrt(Dst->currentState().kinematicParametersError().matrix()(6,6));
  parts.push_back(pFactory.particle(Dst_tk, mDst, chi, ndf, dmDst));

  KinematicParticleVertexFitter VtxFitter;
  RefCountedKinematicTree BKinTree = VtxFitter.fit(parts);
  return BKinTree;
}

RefCountedKinematicTree vtxu::FitVtxJpsiKst(const edm::EventSetup& iSetup, const RefCountedKinematicParticle Jpsi, const RefCountedKinematicParticle Kst, bool mass_constraint) {
  reco::TransientTrack Jpsi_tk = Jpsi->refittedTransientTrack();
  reco::TransientTrack Kst_tk = Kst->refittedTransientTrack();

  std::vector<RefCountedKinematicParticle> parts;// = {Jpsi, Kst};
  KinematicParticleFactoryFromTransientTrack pFactory;
  double chi = 0, ndf = 0;
  float mJpsi = Jpsi->currentState().mass();
  float dmJpsi = sqrt(Jpsi->currentState().kinematicParametersError().matrix()(6,6));
  parts.push_back(pFactory.particle(Jpsi_tk, mJpsi, chi, ndf, dmJpsi));
  float mKst = Kst->currentState().mass();
  float dmKst = sqrt(Kst->currentState().kinematicParametersError().matrix()(6,6));
  parts.push_back(pFactory.particle(Kst_tk, mKst, chi, ndf, dmKst));

  if (!mass_constraint) {
    KinematicParticleVertexFitter VtxFitter;
    RefCountedKinematicTree KinTree = VtxFitter.fit(parts);
    return KinTree;
  }
  else {
    ParticleMass mass = _B0Mass_;
    MultiTrackKinematicConstraint * mass_c = new TwoTrackMassKinematicConstraint(mass);
    KinematicConstrainedVertexFitter kcVtxFitter;
    RefCountedKinematicTree KinTree = kcVtxFitter.fit(parts, mass_c);
    return KinTree;
  }
}

RefCountedKinematicTree vtxu::FitVtxDstK(const edm::EventSetup& iSetup, const RefCountedKinematicParticle Dst, pat::PackedCandidate K, int verbose = 0) {
  // Get transient track builder
  edm::ESHandle<TransientTrackBuilder> TTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",TTBuilder);

  reco::TransientTrack K_tk = TTBuilder->build(fix_track(K.bestTrack()));
  reco::TransientTrack Dst_tk = Dst->refittedTransientTrack();

  KinematicParticleFactoryFromTransientTrack pFactory;
  std::vector<RefCountedKinematicParticle> parts;
  double chi = 0, ndf = 0;
  float mK = _KMass_, dmK = _KMassErr_;
  parts.push_back(pFactory.particle(K_tk, mK, chi, ndf, dmK));
  float mDst = Dst->currentState().mass();
  float dmDst = sqrt(Dst->currentState().kinematicParametersError().matrix()(6,6));
  parts.push_back(pFactory.particle(Dst_tk, mDst, chi, ndf, dmDst));

  KinematicParticleVertexFitter VtxFitter;
  RefCountedKinematicTree BKinTree = VtxFitter.fit(parts);
  return BKinTree;
}

RefCountedKinematicTree vtxu::FitVtxDstPi(const edm::EventSetup& iSetup, const RefCountedKinematicParticle Dst, pat::PackedCandidate pi, int verbose = 0) {
  // Get transient track builder
  edm::ESHandle<TransientTrackBuilder> TTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",TTBuilder);

  reco::TransientTrack Pi_tk = TTBuilder->build(fix_track(pi.bestTrack()));
  reco::TransientTrack Dst_tk = Dst->refittedTransientTrack();

  KinematicParticleFactoryFromTransientTrack pFactory;
  std::vector<RefCountedKinematicParticle> parts;
  double chi = 0, ndf = 0;
  float mPi = _PiMass_, dmPi = _PiMassErr_;
  parts.push_back(pFactory.particle(Pi_tk, mPi, chi, ndf, dmPi));
  float mDst = Dst->currentState().mass();
  float dmDst = sqrt(Dst->currentState().kinematicParametersError().matrix()(6,6));
  parts.push_back(pFactory.particle(Dst_tk, mDst, chi, ndf, dmDst));

  KinematicParticleVertexFitter VtxFitter;
  RefCountedKinematicTree KinTree = VtxFitter.fit(parts);
  return KinTree;
}

RefCountedKinematicTree vtxu::FitVtxMuDstPi(const edm::EventSetup& iSetup, const RefCountedKinematicParticle Dst, pat::PackedCandidate mu, pat::PackedCandidate pi, int verbose = 0) {
  // Get transient track builder
  edm::ESHandle<TransientTrackBuilder> TTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",TTBuilder);

  reco::TransientTrack Dst_tk = Dst->refittedTransientTrack();
  reco::TransientTrack Mu_tk = TTBuilder->build(fix_track(mu.bestTrack()));
  reco::TransientTrack Pi_tk = TTBuilder->build(fix_track(pi.bestTrack()));

  KinematicParticleFactoryFromTransientTrack pFactory;
  std::vector<RefCountedKinematicParticle> parts;
  double chi = 0, ndf = 0;
  float mDst = Dst->currentState().mass();
  float dmDst = sqrt(Dst->currentState().kinematicParametersError().matrix()(6,6));
  parts.push_back(pFactory.particle(Dst_tk, mDst, chi, ndf, dmDst));
  float mMu = _MuMass_, dmMu = _MuMassErr_;
  parts.push_back(pFactory.particle(Mu_tk, mMu, chi, ndf, dmMu));
  float mPi = _PiMass_, dmPi = _PiMassErr_;
  parts.push_back(pFactory.particle(Pi_tk, mPi, chi, ndf, dmPi));

  KinematicParticleVertexFitter VtxFitter;
  RefCountedKinematicTree BKinTree = VtxFitter.fit(parts);
  return BKinTree;
}

vtxu::kinFitResuts vtxu::fitQuality(RefCountedKinematicTree t, double pval_thr){
  vtxu::kinFitResuts out;
  if(t->isValid()) {
    out.isValid = true;
    t->movePointerToTheTop();
    out.chi2 = t->currentDecayVertex()->chiSquared();
    out.dof = t->currentDecayVertex()->degreesOfFreedom();
    if (out.chi2 < 0) {
      fprintf(stdout, "[WARNING] Vertex fit chi2 is less than zero!\n");
      out.isGood = false;
      out.pval = -1;
    }
    else {
      out.pval = ChiSquaredProbability(out.chi2, out.dof);
      if (pval_thr > 0) {
        out.isGood = out.pval > pval_thr;
      }
    }
  }
  return out;
}

std::pair<double,double> vtxu::computeIP3D(const reco::TransientTrack tk, const GlobalVector direction, const reco::Vertex vtx) {
  SignedImpactParameter3D signed_ip3D;
  Measurement1D ip3D = signed_ip3D.apply(tk, direction, vtx).second;
  return make_pair(ip3D.value(), ip3D.error());
}

std::pair<double,double> vtxu::computeIP3D(const edm::EventSetup& iSetup, pat::Muon mu, TLorentzVector p4_direction , RefCountedKinematicVertex vtx) {
  edm::ESHandle<TransientTrackBuilder> TTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",TTBuilder);
  reco::TransientTrack tk = TTBuilder->build(fix_track(mu.muonBestTrack()));

  GlobalVector direction(p4_direction.X(), p4_direction.Y(), p4_direction.Z());

  reco::Vertex::Point vtxPosition(vtx->position().x(),vtx->position().y(),vtx->position().z());
  reco::Vertex::Error vtxError;
  for (uint i=0; i<3; i++) {
    for (uint j=0; j<3; j++) {vtxError(i,j) = vtx->error().matrix()(i,j);}
  }
  float chi2 = vtx->chiSquared();
  float dof = vtx->degreesOfFreedom();
  const reco::Vertex recoVtx(vtxPosition, vtxError, chi2, dof, 2);

  return vtxu::computeIP3D(tk, direction, recoVtx);
}

std::pair<double,double> vtxu::computeIP3D(const edm::EventSetup& iSetup, pat::PackedCandidate cand, TLorentzVector p4_direction , RefCountedKinematicVertex vtx) {
  edm::ESHandle<TransientTrackBuilder> TTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",TTBuilder);
  reco::TransientTrack tk = TTBuilder->build(fix_track(cand.bestTrack()));

  GlobalVector direction(p4_direction.X(), p4_direction.Y(), p4_direction.Z());

  reco::Vertex::Point vtxPosition(vtx->position().x(),vtx->position().y(),vtx->position().z());
  reco::Vertex::Error vtxError;
  for (uint i=0; i<3; i++) {
    for (uint j=0; j<3; j++) {vtxError(i,j) = vtx->error().matrix()(i,j);}
  }
  float chi2 = vtx->chiSquared();
  float dof = vtx->degreesOfFreedom();
  const reco::Vertex recoVtx(vtxPosition, vtxError, chi2, dof, 2);

  return vtxu::computeIP3D(tk, direction, recoVtx);
}

pair<double,double> vtxu::computeDCA(reco::TransientTrack tk, GlobalPoint p, int kind) {
  TrajectoryStateClosestToPoint stateCA = tk.trajectoryStateClosestToPoint(p);
  double dT = stateCA.perigeeParameters().transverseImpactParameter();
  double dL = stateCA.perigeeParameters().longitudinalImpactParameter();
  double dCA = hypot(dL, dT);

  double EdT = stateCA.perigeeError().transverseImpactParameterError();
  double EdL = stateCA.perigeeError().longitudinalImpactParameterError();
  double EdCA = hypot(dT*EdT, dL*EdL)/dCA;
  if (kind == 0) return make_pair(dCA,EdCA);
  else if (kind == 1) return make_pair(dT,EdT);
  else if (kind == 2) return make_pair(dL,EdL);
  else assert(false);
}

pair<double,double> vtxu::computeDCA(const edm::EventSetup& iSetup, pat::PackedCandidate cand, GlobalPoint p, int kind) {
  // Get transient track builder
  edm::ESHandle<TransientTrackBuilder> TTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",TTBuilder);

  reco::TransientTrack tk = TTBuilder->build(fix_track(cand.bestTrack()));
  return vtxu::computeDCA(tk, p, kind);
}

pair<double,double> vtxu::computeDCA(const edm::EventSetup& iSetup, pat::Muon mu, GlobalPoint p, int kind) {
  // Get transient track builder
  edm::ESHandle<TransientTrackBuilder> TTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",TTBuilder);

  reco::TransientTrack tk = TTBuilder->build(fix_track(mu.muonBestTrack()));
  return vtxu::computeDCA(tk, p, kind);
}

TLorentzVector vtxu::getTLVfromKinPart(const RefCountedKinematicParticle p) {
  auto pvec = p->currentState().globalMomentum();
  auto mass = p->currentState().mass();
  TLorentzVector out;
  out.SetXYZM(pvec.x(), pvec.y(), pvec.z(), mass);
  return out;
}

TLorentzVector vtxu::getTLVfromTrack(reco::Track t, double mass) {
  TLorentzVector out;
  auto p3 = t.momentum();
  auto pt = sqrt(p3.perp2());
  out.SetPtEtaPhiM(pt, p3.eta(), p3.phi(), mass);
  return out;
}

TLorentzVector vtxu::getTLVfromCand(pat::PackedCandidate p, double mass) {
  TLorentzVector out;
  out.SetPtEtaPhiM(p.pt(), p.eta(), p.phi(), mass);
  return out;
}

TLorentzVector vtxu::getTLVfromMuon(pat::Muon p, double mass) {
  TLorentzVector out;
  out.SetPtEtaPhiM(p.pt(), p.eta(), p.phi(), mass);
  return out;
}

double vtxu::dPhi(double p1, double p2) {
  double dPhi = p1 - p2;
  double pi = 3.14159265358979323846;
  while (fabs(dPhi) > pi) {
    int sgn = dPhi > 0? 1 : -1;
    dPhi -= sgn*2*pi;
  }
  return dPhi;
}

double vtxu::dR(double p1, double p2, double e1, double e2) {
  return hypot(dPhi(p1, p2), e1 - e2);
}

std::pair<double,double> vtxu::vtxsDistance(reco::VertexRef v1, RefCountedKinematicVertex v2){
  double dx = v2->position().x() - v1->x();
  double dy = v2->position().y() - v1->y();
  double dz = v2->position().z() - v1->z();

  double d = sqrt(dx*dx + dy*dy + dz*dz);

  double dd[3] = {dx/d, dy/d, dz/d};
  auto e1 = v1->covariance();
  auto e2 = v2->error().matrix();

  double Ed2 = 0;
  for(uint i=0; i<3; ++i){
    for(uint j=0; j<3; ++j){
      Ed2 += dd[i]* ( e1.At(i,j)+e2.At(i,j) ) *dd[j];
    }
  }

  return make_pair(d,sqrt(Ed2));
}

std::pair<double,double> vtxu::vtxsDistance(reco::Vertex v1, RefCountedKinematicVertex v2){
  double dx = v2->position().x() - v1.x();
  double dy = v2->position().y() - v1.y();
  double dz = v2->position().z() - v1.z();

  double d = sqrt(dx*dx + dy*dy + dz*dz);

  double dd[3] = {dx/d, dy/d, dz/d};
  auto e1 = v1.covariance();
  auto e2 = v2->error().matrix();

  double Ed2 = 0;
  for(uint i=0; i<3; ++i){
    for(uint j=0; j<3; ++j){
      Ed2 += dd[i]* ( e1.At(i,j)+e2.At(i,j) ) *dd[j];
    }
  }

  return make_pair(d,sqrt(Ed2));
}

// std::pair<double,double> vtxu::vtxsDistance3D(reco::Vertex v1, RefCountedKinematicVertex v2){
//   GlobalPoint v1_p(v1.x(), v1.y(), v1.z());
//   GlobalPoint v2_p(v2->position().x(), v2->position().y(), v2->position().z());
//
//   auto m1 = v1.covariance();
//   GlobalError v1_e(m1.At(1,0), m1.At(1,0), m1.At(1,1), m1.At(2,0), m1.At(2,1), m1.At(2,2));
//
//   auto m2 = v2->error().matrix();
//   GlobalError v2_e(m2.At(1,0), m2.At(1,0), m2.At(1,1), m2.At(2,0), m2.At(2,1), m2.At(2,2));
//
//   VertexDistance3D vtxTool;
//   auto val = vtxTool.distance(v1_p, v1_e, v2_p, v2_e).value();
//   auto err = vtxTool.distance(v1_p, v1_e, v2_p, v2_e).error();
//
//   return make_pair(val, err);
// }

std::pair<double,double> vtxu::vtxsTransverseDistanceFromBeamSpot(const reco::BeamSpot &theBeamSpot, RefCountedKinematicVertex v2)
{
  double dx = v2->position().x() - theBeamSpot.x0();
  double dy = v2->position().y() - theBeamSpot.y0();

  double d = sqrt(dx*dx + dy*dy);

  double dd[2] = {dx/d, dy/d};
  auto e1 = theBeamSpot.covariance();
  auto e2 = v2->error().matrix();

  double Ed2 = 0;
  for(uint i=0; i<2; ++i){
    for(uint j=0; j<2; ++j){
      Ed2 += dd[i]* ( e1.At(i,j)+e2.At(i,j) ) *dd[j];
    }
  }

  return make_pair(d,sqrt(Ed2));
}

std::pair<double,double> vtxu::vtxsTransverseDistance(reco::Vertex v1, RefCountedKinematicVertex v2){
  double dx = v2->position().x() - v1.x();
  double dy = v2->position().y() - v1.y();

  double d = sqrt(dx*dx + dy*dy);

  double dd[2] = {dx/d, dy/d};
  auto e1 = v1.covariance();
  auto e2 = v2->error().matrix();

  double Ed2 = 0;
  for(uint i=0; i<2; ++i){
    for(uint j=0; j<2; ++j){
      Ed2 += dd[i]* ( e1.At(i,j)+e2.At(i,j) ) *dd[j];
    }
  }

  return make_pair(d,sqrt(Ed2));
}


double vtxu::computePointingCos(reco::Vertex vtxP, const RefCountedKinematicVertex vtxKinPartDecay, const RefCountedKinematicParticle p) {
  TVector3 dvtx(vtxKinPartDecay->position().x() - vtxP.position().x(),
                vtxKinPartDecay->position().y() - vtxP.position().y(),
                vtxKinPartDecay->position().z() - vtxP.position().z()
               );

  TVector3 p3(p->currentState().globalMomentum().x(),
              p->currentState().globalMomentum().y(),
              p->currentState().globalMomentum().z()
            );

  double dalpha = dvtx.Angle(p3);
  return cos(dalpha);
}

double vtxu::computePointingCosTransverse(reco::Vertex vtxP, const RefCountedKinematicVertex vtxKinPartDecay, const RefCountedKinematicParticle p) {
  TVector3 dvtx(vtxKinPartDecay->position().x() - vtxP.position().x(),
                vtxKinPartDecay->position().y() - vtxP.position().y(),
                0
               );

  TVector3 p3(p->currentState().globalMomentum().x(),
              p->currentState().globalMomentum().y(),
              0
            );

  double dalpha = dvtx.Angle(p3);
  return cos(dalpha);
}

double vtxu::computeIP(reco::Candidate::Point axis, reco::Candidate::Point trajPoint, reco::Candidate::Vector trajVector, bool linearApprox){
  if(linearApprox) {
    double nx = trajVector.x() / hypot(trajVector.x(), trajVector.y());
    double ny = trajVector.y() / hypot(trajVector.x(), trajVector.y());

    double amp_x = axis.x() - trajPoint.x();
    double amp_y = axis.y() - trajPoint.y();
    double dx = amp_x - nx * (amp_x * nx + amp_y * ny);
    double dy = amp_y - ny * (amp_x * nx + amp_y * ny);

    return hypot(dx, dy);
  }
  else {
    cout << "Not yet implemented" << endl;
    return -1;
  }
}

double vtxu::computeDCA_linApprox(reco::Candidate::Point axis, reco::Candidate::Point trajPoint, reco::Candidate::Vector trajVector, bool transverseOnly){
  double tV_z = transverseOnly ? 0 : trajVector.z();

  double norm = sqrt(trajVector.x()*trajVector.x() + trajVector.y()*trajVector.y() + tV_z*tV_z);
  double nx = trajVector.x() / norm;
  double ny = trajVector.y() / norm;
  double nz = tV_z / norm;

  double amp_x = axis.x() - trajPoint.x();
  double amp_y = axis.y() - trajPoint.y();
  double amp_z = transverseOnly ? 0 : axis.z() - trajPoint.z();

  double amp_dot_n = amp_x * nx + amp_y * ny + amp_z * nz;
  double dx = amp_x - nx * amp_dot_n;
  double dy = amp_y - ny * amp_dot_n;
  double dz = amp_z - nz * amp_dot_n;

  return sqrt(dx*dx + dy*dy + dz*dz);
}

double vtxu::computeCTau(reco::GenParticle p){
  // cout << endl << "Computing lifetime for " << p.pdgId() << endl;
  if(p.numberOfDaughters() == 0) {
    // cout << "Stable particle" << endl;
    return 1e6-1;
  }

  // auto strVtx = [](reco::Candidate::Point v) {
  //   return Form("[%.3f, %.3f, %.3f]", v.x(), v.y(), v.z());
  // };
  auto distanceVtx = [](reco::Candidate::Point v1, reco::Candidate::Point v2) {
    double dd = (v1.x() - v2.x())*(v1.x() - v2.x());
    dd += (v1.y() - v2.y())*(v1.y() - v2.y());
    dd += (v1.z() - v2.z())*(v1.z() - v2.z());
    return sqrt(dd);
  };

  auto decayVtx = p.daughterRefVector()[0]->vertex();
  // cout << "Decay vtx " << strVtx(decayVtx) << endl;
  auto productionVtx = p.vertex();
  // cout << Form("Particle: P=%.3f, ", p.p()) << strVtx(productionVtx) << Form(", Delta=%.4f", distanceVtx(decayVtx, productionVtx)) << endl;

  auto mother = p.mother();
  while(mother->numberOfDaughters() == 1) {
    // cout << "Backward step " << mother->pdgId() << " with only 1 daughter" << endl;
    // cout << Form("Particle: P=%.3f, ", mother->p()) << strVtx(mother->vertex()) << Form(", Delta=%.4f", distanceVtx(mother->vertex(), productionVtx)) << endl;
    productionVtx = mother->vertex();
    mother = mother->mother();
  }
  // cout << "Mother: " << mother->pdgId() << endl;

  auto dl = distanceVtx(decayVtx, productionVtx);
  auto beta = p.p()/p.energy();
  auto gamma = 1/ sqrt(1 - beta*beta);
  auto ctau = dl / (beta*gamma);
  // cout << Form("ctau: %.3f mm", ctau*10) << endl;
  // cout << endl;
  return ctau;
}
