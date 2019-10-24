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
#define __PvalChi2Vtx_min__ 0.05 // Very loose cut
#define __dmD0_max__ 0.1 // loose cut
#define __sigdxy_vtx_PV_min__ 2.0 // Optimized in D0 fitting
#define __dmDst_max__ 0.050 // Loose cut
#define __mass_D0pismu_max__ 7.0 // Some reasonable cut on the mass
// #define __cos_MuDst_PV_min__ 0.8 // Some loose cut tuned on MC
// #define __PvalChi2FakeVtx_min__ 0.90 // Very loose cut


using namespace std;

class B2DstMuDecayTreeProducer : public edm::EDProducer {

public:

    explicit B2DstMuDecayTreeProducer(const edm::ParameterSet &iConfig);
    void AddTLVToOut(TLorentzVector, string, map<string, vector<float>>*);
    bool qualityMuonID(pat::Muon, reco::Vertex);

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

    int n_K = 0, n_pi = 0, n_D0 = 0, n_pis = 0, n_Dst = 0, n_B = 0;

    if (verbose) {cout <<"-------------------- Evt -----------------------\n";}
    //######## Require muon quality ##################
    if(!qualityMuonID(trgMu, primaryVtx)) {
      (*outputNtuplizer)["n_B"] = 0;
      iEvent.put(move(outputNtuplizer), "outputNtuplizer");
      iEvent.put(move(outputVecNtuplizer), "outputVecNtuplizer");
      return;
    }

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
        auto res_piK = vtxu::fitQuality(D0KinTree, __PvalChi2Vtx_min__);
        if(!res_piK.isGood) continue;

        D0KinTree->movePointerToTheTop();
        auto mass_piK = D0KinTree->currentParticle()->currentState().mass();
        if (fabs(mass_piK - mass_D0) > __dmD0_max__) continue;

        auto D0 = D0KinTree->currentParticle();
        auto vtxD0 = D0KinTree->currentDecayVertex();

        auto cos_D0_PV = vtxu::computePointingCos(primaryVtx, vtxD0, D0);
        auto cosT_D0_PV = vtxu::computePointingCosTransverse(primaryVtx, vtxD0, D0);
        auto d_vtxD0_PV = vtxu::vtxsDistance(primaryVtx, vtxD0);
        auto sigd_vtxD0_PV = d_vtxD0_PV.first/d_vtxD0_PV.second;
        auto dxy_vtxD0_PV = vtxu::vtxsTransverseDistance(primaryVtx, vtxD0);
        auto sigdxy_vtxD0_PV = dxy_vtxD0_PV.first/dxy_vtxD0_PV.second;

        if(sigdxy_vtxD0_PV < __sigdxy_vtx_PV_min__) continue;
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
          auto res_D0pis = vtxu::fitQuality(DstKinTree, __PvalChi2Vtx_min__);
          if(!res_D0pis.isGood) continue;

          DstKinTree->movePointerToTheTop();
          auto mass_D0pis = DstKinTree->currentParticle()->currentState().mass();
          if (fabs(mass_D0pis - mass_Dst) > __dmDst_max__) continue;

          auto Dst = DstKinTree->currentParticle();
          auto vtxDst = DstKinTree->currentDecayVertex();

          auto cos_Dst_PV = vtxu::computePointingCos(primaryVtx, vtxDst, Dst);
          auto cosT_Dst_PV = vtxu::computePointingCosTransverse(primaryVtx, vtxDst, Dst);
          auto d_vtxDst_PV = vtxu::vtxsDistance(primaryVtx, vtxDst);
          auto sigd_vtxDst_PV = d_vtxDst_PV.first/d_vtxDst_PV.second;
          auto dxy_vtxDst_PV = vtxu::vtxsTransverseDistance(primaryVtx, vtxDst);
          auto sigdxy_vtxDst_PV = dxy_vtxDst_PV.first/dxy_vtxDst_PV.second;

          n_Dst++;
          if (verbose) {cout << "D* found\n";}

          // Vertex fitting from B and Dst products
          auto BKinTree = vtxu::FitB_D0pismu(iSetup, D0, pis, trgMu);
          auto res = vtxu::fitQuality(BKinTree, __PvalChi2Vtx_min__);
          if(!res.isGood) continue;
          BKinTree->movePointerToTheTop();
          auto mass_D0pismu = BKinTree->currentParticle()->currentState().mass();
          if (mass_D0pismu > __mass_D0pismu_max__) continue; // Last cut! from now on always save
          if (verbose) {cout << "B->D* mu found\n";}

          (*outputVecNtuplizer)["chi2_D0pismu"].push_back(res.chi2);
          (*outputVecNtuplizer)["dof_D0pismu"].push_back(res.dof);
          (*outputVecNtuplizer)["pval_D0pismu"].push_back(res.pval);
          (*outputVecNtuplizer)["mass_D0pismu"].push_back(mass_D0pismu);

          auto D0pismu = BKinTree->currentParticle();
          auto vtxB = BKinTree->currentDecayVertex();

          BKinTree->movePointerToTheFirstChild();
          auto refit_D0 = BKinTree->currentParticle();
          AddTLVToOut(vtxu::getTLVfromKinPart(refit_D0), string("D0_refitD0pismu"), &(*outputVecNtuplizer));
          BKinTree->movePointerToTheNextChild();
          auto refit_pis = BKinTree->currentParticle();
          AddTLVToOut(vtxu::getTLVfromKinPart(refit_pis), string("pis_refitD0pismu"), &(*outputVecNtuplizer));
          BKinTree->movePointerToTheNextChild();
          auto refit_mu = BKinTree->currentParticle();
          AddTLVToOut(vtxu::getTLVfromKinPart(refit_mu), string("mu_refitD0pismu"), &(*outputVecNtuplizer));

          auto p4_Dst_refitD0pismu = vtxu::getTLVfromKinPart(refit_pis) + vtxu::getTLVfromKinPart(refit_D0);
          (*outputVecNtuplizer)["mass_D0pis_refitD0pismu"].push_back(p4_Dst_refitD0pismu.M());
          AddTLVToOut(p4_Dst_refitD0pismu, string("Dst_refitD0pismu"), &(*outputVecNtuplizer));

          auto cos_D0pismu_PV = vtxu::computePointingCos(primaryVtx, vtxB, D0pismu);
          auto cosT_D0pismu_PV = vtxu::computePointingCosTransverse(primaryVtx, vtxB, D0pismu);
          auto d_vtxD0pismu_PV = vtxu::vtxsDistance(primaryVtx, vtxB);
          auto sigd_vtxD0pismu_PV = d_vtxD0pismu_PV.first/d_vtxD0pismu_PV.second;
          auto dxy_vtxD0pismu_PV = vtxu::vtxsDistance(primaryVtx, vtxB);
          auto sigdxy_vtxD0pismu_PV = dxy_vtxD0pismu_PV.first/dxy_vtxD0pismu_PV.second;
          (*outputVecNtuplizer)["cos_D0pismu_PV"].push_back(cos_D0pismu_PV);
          (*outputVecNtuplizer)["cosT_D0pismu_PV"].push_back(cosT_D0pismu_PV);
          (*outputVecNtuplizer)["d_vtxD0pismu_PV"].push_back(d_vtxD0pismu_PV.first);
          (*outputVecNtuplizer)["sigd_vtxD0pismu_PV"].push_back(sigd_vtxD0pismu_PV);
          (*outputVecNtuplizer)["dxy_vtxD0pismu_PV"].push_back(dxy_vtxD0pismu_PV.first);
          (*outputVecNtuplizer)["sigdxy_vtxD0pismu_PV"].push_back(sigdxy_vtxD0pismu_PV);

          /*
          ############################################################################
                                Compute analysis variables
          ############################################################################
          */
          auto p4_vis = vtxu::getTLVfromKinPart(D0pismu);
          TLorentzVector p4_Dst = p4_Dst_refitD0pismu;
          auto p4_mu = vtxu::getTLVfromKinPart(refit_mu);

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

          AddTLVToOut(p4_B, string("B_D0pismu"), &(*outputVecNtuplizer));
          (*outputVecNtuplizer)["M2_miss_D0pismu"].push_back(M2_miss);
          (*outputVecNtuplizer)["q2_D0pismu"].push_back(q2);
          (*outputVecNtuplizer)["Est_mu_D0pismu"].push_back(Est_mu);
          /*
          ############################################################################
          */


          //  Look for the Dst- and trigger muon to make a vertex
          BKinTree = vtxu::FitVtxMuDst(iSetup, Dst, trgMu);
          res = vtxu::fitQuality(BKinTree, __PvalChi2Vtx_min__);
          if(res.isValid) {
            (*outputVecNtuplizer)["isGood_Dstmu"].push_back(res.isGood);
            (*outputVecNtuplizer)["chi2_Dstmu"].push_back(res.chi2);
            (*outputVecNtuplizer)["dof_Dstmu"].push_back(res.dof);
            (*outputVecNtuplizer)["pval_Dstmu"].push_back(res.pval);

            BKinTree->movePointerToTheTop();
            auto Dstmu = BKinTree->currentParticle();
            (*outputVecNtuplizer)["mass_Dstmu"].push_back(Dstmu->currentState().mass());

            auto vtxB = BKinTree->currentDecayVertex();

            auto cos_Dstmu_PV = vtxu::computePointingCos(primaryVtx, vtxB, Dstmu);
            auto cosT_Dstmu_PV = vtxu::computePointingCosTransverse(primaryVtx, vtxB, Dstmu);
            (*outputVecNtuplizer)["cos_Dstmu_PV"].push_back(cos_Dstmu_PV);
            (*outputVecNtuplizer)["cosT_Dstmu_PV"].push_back(cosT_Dstmu_PV);

            BKinTree->movePointerToTheFirstChild();
            auto refit_mu = BKinTree->currentParticle();
            BKinTree->movePointerToTheNextChild();
            auto refit_Dst = BKinTree->currentParticle();

            auto p4_vis = vtxu::getTLVfromKinPart(Dstmu);
            auto p4_Dst = vtxu::getTLVfromKinPart(refit_Dst);
            auto p4_mu = vtxu::getTLVfromKinPart(refit_mu);

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

            AddTLVToOut(p4_B, string("B_Dstmu"), &(*outputVecNtuplizer));
            (*outputVecNtuplizer)["M2_miss_Dstmu"].push_back(M2_miss);
            (*outputVecNtuplizer)["q2_Dstmu"].push_back(q2);
            (*outputVecNtuplizer)["Est_mu_Dstmu"].push_back(Est_mu);
          }
          else {
            (*outputVecNtuplizer)["isGood_Dstmu"].push_back(-1);
            (*outputVecNtuplizer)["chi2_Dstmu"].push_back(-1);
            (*outputVecNtuplizer)["dof_Dstmu"].push_back(-1);
            (*outputVecNtuplizer)["pval_Dstmu"].push_back(-1);
            (*outputVecNtuplizer)["mass_Dstmu"].push_back(-1);
            (*outputVecNtuplizer)["cos_Dstmu_PV"].push_back(-1);
            (*outputVecNtuplizer)["cosT_Dstmu_PV"].push_back(-1);
            TLorentzVector void_p4;
            AddTLVToOut(void_p4, string("B_Dstmu"), &(*outputVecNtuplizer));
            (*outputVecNtuplizer)["M2_miss_Dstmu"].push_back(-1);
            (*outputVecNtuplizer)["q2_Dstmu"].push_back(-1);
            (*outputVecNtuplizer)["Est_mu_Dstmu"].push_back(-1);
          }

          // bool accept_B = mass_MuDst < __mass_MuDst_max__;
          // accept_B &= cos_MuDst_PV > __cos_MuDst_PV_min__;
          // if(!accept_B) continue;

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

          (*outputVecNtuplizer)["sigdxy_K_PV"].push_back(sigdxy_K_PV);
          (*outputVecNtuplizer)["K_norm_chi2"].push_back(K_norm_chi2);
          (*outputVecNtuplizer)["K_N_valid_hits"].push_back(K_N_valid_hits);
          (*outputVecNtuplizer)["sigdxy_pi_PV"].push_back(sigdxy_pi_PV);
          (*outputVecNtuplizer)["pi_norm_chi2"].push_back(pi_norm_chi2);
          (*outputVecNtuplizer)["pi_N_valid_hits"].push_back(pi_N_valid_hits);

          (*outputVecNtuplizer)["chi2_piK"].push_back(res_piK.chi2);
          (*outputVecNtuplizer)["dof_piK"].push_back(res_piK.dof);
          (*outputVecNtuplizer)["pval_piK"].push_back(res_piK.pval);
          (*outputVecNtuplizer)["mass_piK"].push_back(mass_piK);
          D0KinTree->movePointerToTheFirstChild();
          auto refit_K = D0KinTree->currentParticle();
          AddTLVToOut(vtxu::getTLVfromKinPart(refit_K), string("K"), &(*outputVecNtuplizer));
          D0KinTree->movePointerToTheNextChild();
          auto refit_pi = D0KinTree->currentParticle();
          AddTLVToOut(vtxu::getTLVfromKinPart(refit_pi), string("pi"), &(*outputVecNtuplizer));
          (*outputVecNtuplizer)["cos_D0_PV"].push_back(cos_D0_PV);
          (*outputVecNtuplizer)["cosT_D0_PV"].push_back(cosT_D0_PV);
          (*outputVecNtuplizer)["d_vtxD0_PV"].push_back(d_vtxD0_PV.first);
          (*outputVecNtuplizer)["sigd_vtxD0_PV"].push_back(sigd_vtxD0_PV);
          (*outputVecNtuplizer)["dxy_vtxD0_PV"].push_back(dxy_vtxD0_PV.first);
          (*outputVecNtuplizer)["sigdxy_vtxD0_PV"].push_back(sigdxy_vtxD0_PV);

          (*outputVecNtuplizer)["sigdxy_pis_PV"].push_back(sigdxy_pis_PV);
          (*outputVecNtuplizer)["pis_norm_chi2"].push_back(pis_norm_chi2);
          (*outputVecNtuplizer)["pis_N_valid_hits"].push_back(pis_N_valid_hits);

          (*outputVecNtuplizer)["chi2_D0pis"].push_back(res_D0pis.chi2);
          (*outputVecNtuplizer)["dof_D0pis"].push_back(res_D0pis.dof);
          (*outputVecNtuplizer)["pval_D0pis"].push_back(res_D0pis.pval);
          (*outputVecNtuplizer)["mass_D0pis"].push_back(mass_D0pis);
          (*outputVecNtuplizer)["cos_Dst_PV"].push_back(cos_Dst_PV);
          (*outputVecNtuplizer)["cosT_Dst_PV"].push_back(cosT_Dst_PV);
          (*outputVecNtuplizer)["d_vtxDst_PV"].push_back(d_vtxDst_PV.first);
          (*outputVecNtuplizer)["sigd_vtxDst_PV"].push_back(sigd_vtxDst_PV);
          (*outputVecNtuplizer)["dxy_vtxDst_PV"].push_back(dxy_vtxDst_PV.first);
          (*outputVecNtuplizer)["sigdxy_vtxDst_PV"].push_back(sigdxy_vtxDst_PV);

          DstKinTree->movePointerToTheFirstChild();
          auto refitD0pis_pis = DstKinTree->currentParticle();
          AddTLVToOut(vtxu::getTLVfromKinPart(refitD0pis_pis), string("pis_refitD0pis"), &(*outputVecNtuplizer));
          DstKinTree->movePointerToTheNextChild();
          auto refitD0pis_D0 = DstKinTree->currentParticle();
          AddTLVToOut(vtxu::getTLVfromKinPart(refitD0pis_D0), string("D0_refitD0pis"), &(*outputVecNtuplizer));

          // DstKinTree->movePointerToTheFirstChild();
          // auto refit_pis = DstKinTree->currentParticle();
          // AddTLVToOut(vtxu::getTLVfromKinPart(refit_pis), string("pis"), &(*outputVecNtuplizer));
          // DstKinTree->movePointerToTheNextChild();
          // auto refit_D0 = DstKinTree->currentParticle();
          // AddTLVToOut(vtxu::getTLVfromKinPart(refit_D0), string("D0"), &(*outputVecNtuplizer));
          //
          // (*outputVecNtuplizer)["chi2_MuDst"].push_back(chi2_MuDst);
          // (*outputVecNtuplizer)["mass_MuDst"].push_back(mass_MuDst);
          // (*outputVecNtuplizer)["cos_MuDst_PV"].push_back(cos_MuDst_PV);
          // (*outputVecNtuplizer)["d_vtxB_PV"].push_back(d_vtxB_PV.first);
          // (*outputVecNtuplizer)["sigd_vtxB_PV"].push_back(sigd_vtxB_PV);
          //
          // AddTLVToOut(p4_Dst, string("Dst"), &(*outputVecNtuplizer));
          // AddTLVToOut(p4_mu, string("mu"), &(*outputVecNtuplizer));
          // AddTLVToOut(p4_B, string("B"), &(*outputVecNtuplizer));
          //
          // (*outputVecNtuplizer)["M2_miss"].push_back(M2_miss);
          // (*outputVecNtuplizer)["q2"].push_back(q2);
          // (*outputVecNtuplizer)["Est_mu"].push_back(Est_mu);
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


void B2DstMuDecayTreeProducer::AddTLVToOut(TLorentzVector v, string n, map<string, vector<float>>* outv) {
  (*outv)[n+"_pt"].push_back(v.Pt());
  (*outv)[n+"_pz"].push_back(v.Pz());
  (*outv)[n+"_eta"].push_back(v.Eta());
  (*outv)[n+"_phi"].push_back(v.Phi());
  return;
}

bool B2DstMuDecayTreeProducer::qualityMuonID(pat::Muon m, reco::Vertex pVtx) {
  if(m.innerTrack().isNull()) return false;

  if(m.innerTrack()->hitPattern().pixelLayersWithMeasurement() < 2) return false;
  if(!m.innerTrack()->quality(reco::TrackBase::highPurity)) return false;
  if(!m.isGood("TMOneStationTight")) return false;
  if(m.innerTrack()->normalizedChi2() > 1.8) return false;

  double dxy = m.innerTrack()->dxy(pVtx.position());
  float sigdxy = fabs(dxy)/m.innerTrack()->dxyError();
  if (sigdxy < 2) return false;

  return true;
}

DEFINE_FWK_MODULE(B2DstMuDecayTreeProducer);
