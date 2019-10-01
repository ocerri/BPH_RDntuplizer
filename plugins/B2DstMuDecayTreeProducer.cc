#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include <iostream>
#include <string>
#include <map>
#include <TMath.h>

#include "VtxUtils.hh"

#define __pThad_min__ 0.5 // loose cut
#define __dzMax__ 1.0
#define __dRMax__ 2.0
#define __sigIPpfCand_min__ 2. // loose cut
#define __PvalChi2Vtx_max__ 0.99 // Very loose cut
#define __dmD0_max__ 0.065 // 5*0.013 GeV(i.e. 3 sigma) form the D0 mass
#define __cos_D0_PV_min__ 0.97 // Optimized in D0 fitting
#define __sigd_vtxD0_PV_min__ 2.0 // Optimized in D0 fitting
#define __dmDst_max__ 0.0040 // 5*0.8 MeV (i.e 3 inflated sigma) from Dst mass
#define __mass_MuDst_max__ 7.0 // Some reasonable cut on the mass
#define __cos_MuDst_PV_min__ 0.8 // Some loose cut tuned on MC
#define __PvalChi2FakeVtx_min__ 0.90 // Very loose cut


using namespace std;

class B2DstMuDecayTreeProducer : public edm::EDProducer {

public:

    explicit B2DstMuDecayTreeProducer(const edm::ParameterSet &iConfig);
    // void DeclareTLVToOut(string, map<string, vector<float>>*);
    void AddTLVToOut(TLorentzVector, string, map<string, vector<float>>*);

    ~B2DstMuDecayTreeProducer() override {};

private:
    virtual void produce(edm::Event&, const edm::EventSetup&);

    // ----------member data ---------------------------
    edm::EDGetTokenT<vector<pat::PackedCandidate>> PFCandSrc_;
    edm::EDGetTokenT<vector<reco::Vertex>> vtxSrc_;
    edm::EDGetTokenT<vector<pat::Muon>> TrgMuonSrc_;

    double mass_Mu  = 0.10565;
    double mass_Pi  = 0.13957;
    double mass_K   = 0.49367;
    double mass_D0  = 1.86483;
    double mass_Dst = 2.01026;
    double mass_B0  = 5.27963;

    int verbose = 0;
};



B2DstMuDecayTreeProducer::B2DstMuDecayTreeProducer(const edm::ParameterSet &iConfig)
{
    verbose = iConfig.getParameter<int>( "verbose" );
    PFCandSrc_ = consumes<vector<pat::PackedCandidate>>(edm::InputTag("packedPFCandidates"));
    vtxSrc_ = consumes<vector<reco::Vertex>>(edm::InputTag("offlineSlimmedPrimaryVertices"));
    TrgMuonSrc_ = consumes<vector<pat::Muon>> ( iConfig.getParameter<edm::InputTag>( "trgMuons" ) );

    produces<map<string, float>>("outputNtuplizer");
    produces<map<string, vector<float>>>("outputVecNtuplizer");
}


void B2DstMuDecayTreeProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    edm::Handle<vector<pat::PackedCandidate>> pfCandHandle;
    iEvent.getByToken(PFCandSrc_, pfCandHandle);
    unsigned int N_pfCand = pfCandHandle->size();

    edm::Handle<vector<reco::Vertex>> vtxHandle;
    iEvent.getByToken(vtxSrc_, vtxHandle);
    auto primaryVtx = (*vtxHandle)[0];

    edm::Handle<vector<pat::Muon>> trgMuonsHandle;
    iEvent.getByToken(TrgMuonSrc_, trgMuonsHandle);
    auto trgMu = (*trgMuonsHandle)[0];

    // Output collection
    unique_ptr<map<string, float>> outputNtuplizer(new map<string, float>);
    unique_ptr<map<string, vector<float>>> outputVecNtuplizer(new map<string, vector<float>>);

    // auto mu = (*trgMuonsHandle)[0];
    // pat::PackedCandidate trgMu;
    // uint i_trgMu = 0;
    // for (i_trgMu = 0; i_trgMu < N_pfCand; ++i_trgMu) {
    //   const pat::PackedCandidate & p = (*pfCandHandle)[i_trgMu];
    //   if (fabs(mu.pt() - p.pt())/mu.pt() > 5e-3) continue;
    //   if (fabs(mu.eta() - p.eta()) > 5e-3) continue;
    //   if (fabs(vtxu::dPhi(mu.phi(), p.phi())) > 5e-3) continue;
    //   trgMu = p;
    //   break;
    // }
    //
    // auto vtxMu = trgMu.vertexRef();
    // GlobalPoint p_vtxMu(vtxMu->x(), vtxMu->y(), vtxMu->z());


    int n_K = 0, n_pi = 0, n_D0 = 0, n_pis = 0, n_Dst = 0, n_B = 0;

    // (*outputVecNtuplizer)["chi2_Kpi"] = {};
    // (*outputVecNtuplizer)["mass_Kpi"] = {};
    // (*outputVecNtuplizer)["dca_Kpi_vtxMu"] = {};
    // (*outputVecNtuplizer)["sigdca_Kpi_vtxMu"] = {};
    // (*outputVecNtuplizer)["cos_Kpi_vtxMu"] = {};
    // (*outputVecNtuplizer)["d_vtxKpi_vtxMu"] = {};
    // (*outputVecNtuplizer)["sigd_vtxKpi_vtxMu"] = {};
    //
    // (*outputVecNtuplizer)["chi2_D0pis"] = {};
    // (*outputVecNtuplizer)["mass_D0pis"] = {};
    // (*outputVecNtuplizer)["dca_D0pis_vtxMu"] = {};
    // (*outputVecNtuplizer)["sigdca_D0pis_vtxMu"] = {};
    // (*outputVecNtuplizer)["cos_D0pis_vtxMu"] = {};
    // (*outputVecNtuplizer)["d_vtxD0pis_vtxMu"] = {};
    // (*outputVecNtuplizer)["sigd_vtxD0pis_vtxMu"] = {};
    //
    // (*outputVecNtuplizer)["chi2_MuDst"] = {};
    // (*outputVecNtuplizer)["mass_MuDst"] = {};
    // (*outputVecNtuplizer)["cos_MuDst_vtxBest"] = {};
    //
    // DeclareTLVToOut("K", &(*outputVecNtuplizer));
    // DeclareTLVToOut("pi", &(*outputVecNtuplizer));
    // DeclareTLVToOut("pis", &(*outputVecNtuplizer));
    // DeclareTLVToOut("D0", &(*outputVecNtuplizer));
    // DeclareTLVToOut("Dst", &(*outputVecNtuplizer));
    // DeclareTLVToOut("B", &(*outputVecNtuplizer));
    //
    // (*outputVecNtuplizer)["M2_miss"] = {};
    // (*outputVecNtuplizer)["q2"] = {};
    // (*outputVecNtuplizer)["Est_mu"] = {};

    if (verbose) {cout <<"-------------------- Evt -----------------------\n";}
    /*
    ############################################################################
                              Look for the K+ (321)
    ############################################################################
    */
    for(uint i_K = 0; i_K < N_pfCand; ++i_K) {
      // if(i_K==i_trgMu) continue;

      const pat::PackedCandidate & K = (*pfCandHandle)[i_K];
      if (!K.hasTrackDetails()) continue;
      //Require a positive charged hadron
      if (K.pdgId() != 211 ) continue;
      if (K.pt() < __pThad_min__) continue;
      // Require to be close to the trigger muon;
      auto K_tk = K.bestTrack();
      if (fabs(K_tk->dz(primaryVtx.position()) - trgMu.muonBestTrack()->dz(primaryVtx.position())) > __dzMax__) continue;
      if (vtxu::dR(K.phi(), trgMu.phi(), K.eta(), trgMu.eta()) > __dRMax__) continue;
      // Require significant impact parameter
      auto dxy = K_tk->dxy(primaryVtx.position());
      auto sigdxy_K_PV = fabs(dxy)/K_tk->dxyError();
      auto K_norm_chi2 = K_tk->normalizedChi2();
      auto K_N_valid_hits = K_tk->numberOfValidHits();
      if (sigdxy_K_PV < __sigIPpfCand_min__) continue;

      n_K++;

      /*
      ############################################################################
                               Look for the pi- (-211)
      ############################################################################
      */
      for(uint i_pi = 0; i_pi < N_pfCand; ++i_pi) {
        if(i_pi==i_K) continue;

        const pat::PackedCandidate & pi = (*pfCandHandle)[i_pi];
        if (!pi.hasTrackDetails()) continue;
        //Require a negative charged hadron
        if (pi.pdgId() != -211 ) continue;
        //Require a minimum pt
        if(pi.pt() < __pThad_min__) continue;
        // Require to be close to the trigger muon;
        auto pi_tk = pi.bestTrack();
        if (fabs(pi_tk->dz(primaryVtx.position()) - trgMu.muonBestTrack()->dz(primaryVtx.position())) > __dzMax__) continue;
        if (vtxu::dR(pi.phi(), trgMu.phi(), pi.eta(), trgMu.eta()) > __dRMax__) continue;
        // Require significant impact parameter
        auto dxy = pi_tk->dxy(primaryVtx.position());
        auto sigdxy_pi_PV = fabs(dxy)/pi.dxyError();
        auto pi_norm_chi2 = pi_tk->normalizedChi2();
        auto pi_N_valid_hits = pi_tk->numberOfValidHits();
        if (sigdxy_pi_PV < __sigIPpfCand_min__) continue;

        n_pi++;

        //Fit the vertex
        auto D0KinTree = vtxu::FitD0(iSetup, pi, K, false);
        bool accept_D0 = false;
        double chi2_Kpi;
        if(D0KinTree->isValid()) {
          D0KinTree->movePointerToTheTop();
          chi2_Kpi = D0KinTree->currentDecayVertex()->chiSquared();
          accept_D0 = chi2_Kpi > 0 && chi2_Kpi < TMath::ChisquareQuantile(__PvalChi2Vtx_max__, D0KinTree->currentDecayVertex()->degreesOfFreedom());
        }
        if(!accept_D0) continue;

        auto mass_Kpi = D0KinTree->currentParticle()->currentState().mass();
        if (fabs(mass_Kpi - mass_D0) > __dmD0_max__) continue;

        // Refit with mass constraint
        D0KinTree = vtxu::FitD0(iSetup, pi, K, true);
        if(!D0KinTree->isValid()) continue;
        D0KinTree->movePointerToTheTop();
        auto D0 = D0KinTree->currentParticle();
        auto vtxD0 = D0KinTree->currentDecayVertex();

        auto cos_D0_PV = vtxu::computePointingCos(primaryVtx, vtxD0, D0);
        auto d_vtxD0_PV = vtxu::vtxsDistance(primaryVtx, vtxD0);
        auto sigd_vtxD0_PV = d_vtxD0_PV.first/d_vtxD0_PV.second;

        // Make selection cut to get good D0
        accept_D0 = cos_D0_PV > __cos_D0_PV_min__;
        accept_D0 &= sigd_vtxD0_PV > __sigd_vtxD0_PV_min__;
        if(!accept_D0) continue;
        n_D0++;

        /*
        ############################################################################
                         Look for the soft pion to make a Dst-
        ############################################################################
        */
        for(uint i_pis = 0; i_pis < N_pfCand; ++i_pis) {
          if(i_pis==i_K || i_pis==i_pi) continue;

          const pat::PackedCandidate & pis = (*pfCandHandle)[i_pis];
          if (!pis.hasTrackDetails()) continue;
          //Require a negative charged hadron
          if (pis.pdgId() != -211 ) continue;

          // Require to be close to the trigger muon;
          auto pis_tk = pis.bestTrack();
          if (fabs(pis_tk->dz(primaryVtx.position()) - trgMu.muonBestTrack()->dz(primaryVtx.position())) > __dzMax__) continue;
          if (vtxu::dR(pis.phi(), trgMu.phi(), pis.eta(), trgMu.eta()) > __dRMax__) continue;
          // Require significant impact parameter
          auto dxy = pis_tk->dxy(primaryVtx.position());
          auto sigdxy_pis_PV = fabs(dxy)/pis.dxyError();
          auto pis_norm_chi2 = pis_tk->normalizedChi2();
          auto pis_N_valid_hits = pis_tk->numberOfValidHits();
          if (sigdxy_pis_PV < __sigIPpfCand_min__) continue;

          n_pis++;

          // Fit the Dst vertex
          auto DstKinTree = vtxu::FitDst(iSetup, pis, D0, false);
          bool accept_Dst = false;
          double chi2_D0pis;
          if(DstKinTree->isValid()) {
            DstKinTree->movePointerToTheTop();
            chi2_D0pis = DstKinTree->currentDecayVertex()->chiSquared();
            accept_Dst = chi2_D0pis > 0 && chi2_D0pis < TMath::ChisquareQuantile(__PvalChi2Vtx_max__, DstKinTree->currentDecayVertex()->degreesOfFreedom());
          }
          if(!accept_Dst) continue;

          auto mass_D0pis = DstKinTree->currentParticle()->currentState().mass();
          if (fabs(mass_D0pis - mass_Dst) > __dmDst_max__) continue;

          // Refit with mass constraint
          DstKinTree = vtxu::FitDst(iSetup, pis, D0, true);
          if(!DstKinTree->isValid()) continue;
          DstKinTree->movePointerToTheTop();
          auto Dst = DstKinTree->currentParticle();
          auto vtxDst = DstKinTree->currentDecayVertex();

          auto cos_Dst_PV = vtxu::computePointingCos(primaryVtx, vtxDst, Dst);
          auto d_vtxDst_PV = vtxu::vtxsDistance(primaryVtx, vtxDst);
          auto sigd_vtxDst_PV = d_vtxDst_PV.first/d_vtxDst_PV.second;

          n_Dst++;
          if (verbose) {cout << "D* found\n";}

          //  Look for the Dst- and trigger muon to make a vertex
          auto BKinTree = vtxu::FitVtxMuDst(iSetup, Dst, trgMu);
          bool accept_MuDst = false;
          double chi2_MuDst;
          if(BKinTree->isValid()) {
            BKinTree->movePointerToTheTop();
            chi2_MuDst = BKinTree->currentDecayVertex()->chiSquared();
            auto max_chi2 = TMath::ChisquareQuantile(__PvalChi2Vtx_max__, BKinTree->currentDecayVertex()->degreesOfFreedom());
            accept_MuDst = chi2_MuDst > 0 && chi2_MuDst < max_chi2;
          }
          if(!accept_MuDst) continue;

          auto MuDst = BKinTree->currentParticle();
          auto vtxB = BKinTree->currentDecayVertex();

          BKinTree->movePointerToTheFirstChild();
          auto refit_Mu = BKinTree->currentParticle();
          BKinTree->movePointerToTheNextChild();
          auto refit_Dst = BKinTree->currentParticle();

          auto mass_MuDst = MuDst->currentState().mass();
          auto cos_MuDst_PV = vtxu::computePointingCos(primaryVtx, vtxB, MuDst);
          auto d_vtxB_PV = vtxu::vtxsDistance(primaryVtx, vtxB);
          auto sigd_vtxB_PV = d_vtxB_PV.first/d_vtxB_PV.second;

          bool accept_B = mass_MuDst < __mass_MuDst_max__;
          accept_B &= cos_MuDst_PV > __cos_MuDst_PV_min__;
          if(!accept_B) continue;
          if (verbose) {cout << "B->D* mu found\n";}

          /*
          ############################################################################
                  Veto the presence of additional tracks in the D* mu vertex
          ############################################################################
          */
          // int N_compatible_tk = 0;
          // for(uint i_tk = 0; i_tk < N_pfCand; ++i_tk) {
          //   // Pion different from the pion from D0
          //   if( i_tk==i_trgMu || i_tk==i_K || i_tk==i_pi || i_tk==i_pis) continue;
          //
          //   const pat::PackedCandidate & ptk = (*pfCandHandle)[i_tk];
          //   //Require a charged candidate
          //   if ( ptk.charge() == 0 ) continue;
          //   if (!ptk.hasTrackDetails()) continue;
          //
          //   // Require to be close to the trigger muon;
          //   if (fabs(ptk.dz() - trgMu.dz()) > __dzMax__) continue;
          //   if (vtxu::dR(ptk.phi(), trgMu.phi(), ptk.eta(), trgMu.eta()) > __dRMax__) continue;
          //
          //   //  Look for the Dst- and trigger muon to make a vertex
          //   auto VtxKinTree = vtxu::FitVtxDstPi(iSetup, Dst, ptk, 0);
          //   // auto VtxKinTree = vtxu::FitVtxMuDstPi(iSetup, Dst, trgMu, ptk, 0);
          //   bool accept_Vtx = false;
          //   double chi2_vtx;
          //   if(VtxKinTree->isValid()) {
          //     VtxKinTree->movePointerToTheTop();
          //     auto vtx = VtxKinTree->currentDecayVertex();
          //     chi2_vtx = vtx->chiSquared();
          //     auto max_chi2 = TMath::ChisquareQuantile(__PvalChi2FakeVtx_min__, vtx->degreesOfFreedom());
          //     if (chi2_vtx > 0 && chi2_vtx < max_chi2) accept_Vtx = true;
          //   }
          //   if(!accept_Vtx) continue;
          //
          //   auto particle = VtxKinTree->currentParticle();
          //   auto mass_particle = particle->currentState().mass();
          //   if(verbose) {cout << Form("Tk pt=%.1f eta=%.1f phi=%.1f - Vtx chi2=%.1f Mass=%.3f ", ptk.pt(), ptk.eta(), ptk.phi(), chi2_vtx, mass_particle) << endl;}
          //
          //   // accept_Vtx &=  mass_particle >
          //   N_compatible_tk++;
          // }
          // if(verbose) {cout << "Compatible tracks: " << N_compatible_tk << endl;}


          n_B++;

          /*
          ############################################################################
                                Compute analysis variables
          ############################################################################
          */
          auto p4_vis = vtxu::getTLVfromKinPart(MuDst);
          auto p4_Dst = vtxu::getTLVfromKinPart(refit_Dst);
          auto p4_mu = vtxu::getTLVfromKinPart(refit_Mu);

          // // ------------- Transverse Approx ------------------
          double pt_B_reco = p4_vis.Pt() * mass_B0/ p4_vis.M();
          TVector3 flightB(vtxB->position().x() - primaryVtx.position().x(),
                           vtxB->position().y() - primaryVtx.position().y(),
                           vtxB->position().z() - primaryVtx.position().z()
                          );
          auto B_vect = flightB * ( pt_B_reco / flightB.Perp() );
          TLorentzVector p4_B;
          p4_B.SetVectM(B_vect, mass_B0);

          auto M2_miss = (p4_B - p4_vis).M2();
          auto q2 = (p4_B - p4_Dst).M2();

          TLorentzVector p4st_mu(p4_mu);
          p4st_mu.Boost(-1*p4_B.BoostVector());
          auto Est_mu = p4st_mu.E();

          (*outputVecNtuplizer)["sigdxy_K_PV"].push_back(sigdxy_K_PV);
          (*outputVecNtuplizer)["K_norm_chi2"].push_back(K_norm_chi2);
          (*outputVecNtuplizer)["K_N_valid_hits"].push_back(K_N_valid_hits);
          (*outputVecNtuplizer)["sigdxy_pi_PV"].push_back(sigdxy_pi_PV);
          (*outputVecNtuplizer)["pi_norm_chi2"].push_back(pi_norm_chi2);
          (*outputVecNtuplizer)["pi_N_valid_hits"].push_back(pi_N_valid_hits);

          (*outputVecNtuplizer)["chi2_Kpi"].push_back(chi2_Kpi);
          (*outputVecNtuplizer)["mass_Kpi"].push_back(mass_Kpi);
          D0KinTree->movePointerToTheFirstChild();
          auto refit_K = D0KinTree->currentParticle();
          AddTLVToOut(vtxu::getTLVfromKinPart(refit_K), string("K"), &(*outputVecNtuplizer));
          D0KinTree->movePointerToTheNextChild();
          auto refit_pi = D0KinTree->currentParticle();
          AddTLVToOut(vtxu::getTLVfromKinPart(refit_pi), string("pi"), &(*outputVecNtuplizer));
          (*outputVecNtuplizer)["cos_D0_PV"].push_back(cos_D0_PV);
          (*outputVecNtuplizer)["d_vtxD0_PV"].push_back(d_vtxD0_PV.first);
          (*outputVecNtuplizer)["sigd_vtxD0_PV"].push_back(sigd_vtxD0_PV);

          (*outputVecNtuplizer)["sigdxy_pis_PV"].push_back(sigdxy_pis_PV);
          (*outputVecNtuplizer)["pis_norm_chi2"].push_back(pis_norm_chi2);
          (*outputVecNtuplizer)["pis_N_valid_hits"].push_back(pis_N_valid_hits);

          (*outputVecNtuplizer)["chi2_D0pis"].push_back(chi2_D0pis);
          (*outputVecNtuplizer)["mass_D0pis"].push_back(mass_D0pis);
          DstKinTree->movePointerToTheFirstChild();
          auto refit_pis = DstKinTree->currentParticle();
          AddTLVToOut(vtxu::getTLVfromKinPart(refit_pis), string("pis"), &(*outputVecNtuplizer));
          DstKinTree->movePointerToTheNextChild();
          auto refit_D0 = DstKinTree->currentParticle();
          AddTLVToOut(vtxu::getTLVfromKinPart(refit_D0), string("D0"), &(*outputVecNtuplizer));
          (*outputVecNtuplizer)["cos_Dst_PV"].push_back(cos_Dst_PV);
          (*outputVecNtuplizer)["d_vtxDst_PV"].push_back(d_vtxDst_PV.first);
          (*outputVecNtuplizer)["sigd_vtxDst_PV"].push_back(sigd_vtxDst_PV);

          (*outputVecNtuplizer)["chi2_MuDst"].push_back(chi2_MuDst);
          (*outputVecNtuplizer)["mass_MuDst"].push_back(mass_MuDst);
          (*outputVecNtuplizer)["cos_MuDst_PV"].push_back(cos_MuDst_PV);
          (*outputVecNtuplizer)["d_vtxB_PV"].push_back(d_vtxB_PV.first);
          (*outputVecNtuplizer)["sigd_vtxB_PV"].push_back(sigd_vtxB_PV);

          AddTLVToOut(p4_Dst, string("Dst"), &(*outputVecNtuplizer));
          AddTLVToOut(p4_mu, string("mu"), &(*outputVecNtuplizer));
          AddTLVToOut(p4_B, string("B"), &(*outputVecNtuplizer));

          (*outputVecNtuplizer)["M2_miss"].push_back(M2_miss);
          (*outputVecNtuplizer)["q2"].push_back(q2);
          (*outputVecNtuplizer)["Est_mu"].push_back(Est_mu);
        }


      }

    }

    (*outputNtuplizer)["Run"] = iEvent.run();
    (*outputNtuplizer)["LumiBlock"] = iEvent.luminosityBlock();
    (*outputNtuplizer)["eventNumber"] = iEvent.id().event();

    (*outputNtuplizer)["n_K"] = n_K;
    (*outputNtuplizer)["n_pi"] = n_pi;
    (*outputNtuplizer)["n_D0"] = n_D0;
    (*outputNtuplizer)["n_pis"] = n_pis;
    (*outputNtuplizer)["n_Dst"] = n_Dst;
    (*outputNtuplizer)["n_B"] = n_B;


    iEvent.put(move(outputNtuplizer), "outputNtuplizer");
    iEvent.put(move(outputVecNtuplizer), "outputVecNtuplizer");
    return;
}


// void B2DstMuDecayTreeProducer::DeclareTLVToOut(string n, map<string, vector<float>>* outv) {
//   (*outv)[n+"_pt"] = {};
//   (*outv)[n+"_pz"] = {};
//   (*outv)[n+"_eta"] = {};
//   (*outv)[n+"_phi"] = {};
//   (*outv)[n+"_P"] = {};
//   (*outv)[n+"_E"] = {};
//   return;
// }


void B2DstMuDecayTreeProducer::AddTLVToOut(TLorentzVector v, string n, map<string, vector<float>>* outv) {
  (*outv)[n+"_pt"].push_back(v.Pt());
  (*outv)[n+"_pz"].push_back(v.Pz());
  (*outv)[n+"_eta"].push_back(v.Eta());
  (*outv)[n+"_phi"].push_back(v.Phi());
  (*outv)[n+"_P"].push_back(v.P());
  (*outv)[n+"_E"].push_back(v.E());
  return;
}

DEFINE_FWK_MODULE(B2DstMuDecayTreeProducer);
