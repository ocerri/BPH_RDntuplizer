#include <iostream>
#include <string>
#include <map>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <TMath.h>
#include <TLorentzVector.h>

#include "VtxUtils.hh"

#define __pThad_min__ 0.5 // loose cut
#define __dzMax__ 1.0
#define __dRMax__ 2.0
#define __sigIPpfCand_min__ 2. // loose cut
#define __PvalChi2Vtx_min__ 0.05 // loose cut
#define __dmD0_max__ 0.05 // loose cut
#define __sigdxy_vtx_PV_min__ 5.0 // loose cut
#define __dmDst_max__ 0.050 // loose cut
#define __mass_D0pismu_max__ 10.0 // Some reasonable cut on the mass
#define __pTaddTracks_min__ 0.3 // loose cut
#define __mass_D0pismupi_max__ 10.0 // Some reasonable cut on the mass


using namespace std;

class TagAndProbeBp2DmstPipPipProducer : public edm::EDProducer {

public:

    explicit TagAndProbeBp2DmstPipPipProducer(const edm::ParameterSet &iConfig);
    void AddTLVToOut(TLorentzVector, string, map<string, vector<float>>*);
    bool qualityMuonID(pat::Muon, reco::Vertex);
    void updateCounter(int, vector<bool>&);

    ~TagAndProbeBp2DmstPipPipProducer() {
      cout << Form("TagAndProbeBp2DmstPipPipProducer counters:%d\n", (int)counters.size());
      for(auto v : counters) {
        cout << v << endl;
      }
      cout << "\n\n";

      edm::Service<TFileService> fs;
      TH1I* hCounters;
      int nC = counters.size();
      hCounters = fs->make<TH1I>("hB2DstPiDecayCounters", "hB2DstPiDecayCounters", nC, 0, nC);
      for(int i=0; i < nC; i++) hCounters->SetBinContent(i+1, counters[i]);
    };

private:
    virtual void produce(edm::Event&, const edm::EventSetup&);

    // ----------member data ---------------------------
    edm::EDGetTokenT<vector<pat::PackedCandidate>> PFCandSrc_;
    edm::EDGetTokenT<vector<reco::Vertex>> vtxSrc_;
    edm::EDGetTokenT<vector<pat::Muon>> TrgMuonSrc_;

    int charge_K = +1;
    int charge_pi = -1;
    int charge_pi1 = +1;
    int charge_pi2 = +1;
    int charge_pis = -1;

    double mass_Mu  = 0.10565;
    double mass_pi  = 0.13957;
    double mass_K   = 0.49367;
    double mass_D0  = 1.86483;
    double mass_Dst = 2.01026;
    double mass_B0  = 5.27963;

    int verbose = 0;

    vector<uint> counters;
};



TagAndProbeBp2DmstPipPipProducer::TagAndProbeBp2DmstPipPipProducer(const edm::ParameterSet &iConfig)
{
    PFCandSrc_ = consumes<vector<pat::PackedCandidate>>(edm::InputTag("packedPFCandidates"));
    vtxSrc_ = consumes<vector<reco::Vertex>>(edm::InputTag("offlineSlimmedPrimaryVertices"));
    TrgMuonSrc_ = consumes<vector<pat::Muon>> ( iConfig.getParameter<edm::InputTag>( "trgMuons" ) );

    charge_K = iConfig.getParameter<int>( "charge_K" );
    charge_pi = iConfig.getParameter<int>( "charge_pi" );
    charge_pi1 = iConfig.getParameter<int>( "charge_pi1" );
    charge_pi2 = iConfig.getParameter<int>( "charge_pi2" );
    charge_pis = iConfig.getParameter<int>( "charge_pis" );

    verbose = iConfig.getParameter<int>( "verbose" );

    produces<map<string, float>>("outputNtuplizer");
    produces<map<string, vector<float>>>("outputVecNtuplizer");
    for(uint i=0; i<10; i++) counters.push_back(0);
}


void TagAndProbeBp2DmstPipPipProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
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

    int n_K = 0, n_pi = 0, n_D0 = 0, n_pi1 = 0, n_pi2 = 0, n_B = 0;
    int n_pis = 0, n_Dst = 0, n_fullB = 0;

    if (verbose) {cout <<"-------------------- Evt -----------------------\n";}
    vector<bool> countersFlag(counters.size(), false);
    counters[0]++;

    vector<reco::Vertex> possibleVtxs = {};

    for(uint i_vtx = 0; i_vtx<vtxHandle->size(); i_vtx++) {
      auto vtx = (*vtxHandle)[i_vtx];
      if(vtx.normalizedChi2() > 1.5) continue;
      bool isSoft = trgMu.isSoftMuon(vtx);
      if(isSoft){
        possibleVtxs.push_back(vtx);
        // if (i_vtx==0) softMuForPV = true;
      }
    }
    if (possibleVtxs.size() == 0) {
      (*outputNtuplizer)["n_B"] = 0;
      iEvent.put(move(outputNtuplizer), "outputNtuplizer");
      iEvent.put(move(outputVecNtuplizer), "outputVecNtuplizer");
      return;
    }
    updateCounter(1, countersFlag);


    /*
    ############################################################################
                              Look for the K+ (321)
    ############################################################################
    */
    for(uint i_K = 0; i_K < N_pfCand; ++i_K) {
      const pat::PackedCandidate & K = (*pfCandHandle)[i_K];
      if (!K.hasTrackDetails()) continue;
      //Require a positive charged hadron
      if (K.pdgId() != charge_K*211 ) continue;
      if (K.pt() < __pThad_min__) continue;
      auto K_tk = K.bestTrack();
      // Require to be close to the triggering muon
      if (fabs(K_tk->dz(primaryVtx.position()) - trgMu.muonBestTrack()->dz(primaryVtx.position())) > __dzMax__) continue;
      // Require significant impact parameter
      auto dxy = K_tk->dxy(primaryVtx.position());
      auto sigdxy_K_PV = fabs(dxy)/K_tk->dxyError();
      auto K_norm_chi2 = K_tk->normalizedChi2();
      auto K_N_valid_hits = K_tk->numberOfValidHits();
      if (sigdxy_K_PV < __sigIPpfCand_min__) continue;
      updateCounter(2, countersFlag);

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
        if (pi.pdgId() != charge_pi*211 ) continue;
        //Require a minimum pt
        if(pi.pt() < __pThad_min__) continue;
        // Require to be close to the trigger muon;
        auto pi_tk = pi.bestTrack();
        if (fabs(pi_tk->dz(primaryVtx.position()) - K_tk->dz(primaryVtx.position())) > __dzMax__) continue;
        if (vtxu::dR(pi.phi(), K.phi(), pi.eta(), K.eta()) > __dRMax__) continue;
        // Require significant impact parameter
        auto dxy = pi_tk->dxy(primaryVtx.position());
        auto sigdxy_pi_PV = fabs(dxy)/pi.dxyError();
        auto pi_norm_chi2 = pi_tk->normalizedChi2();
        auto pi_N_valid_hits = pi_tk->numberOfValidHits();
        if (sigdxy_pi_PV < __sigIPpfCand_min__) continue;
        updateCounter(3, countersFlag);

        n_pi++;

        //Fit the vertex
        auto D0KinTree = vtxu::FitD0(iSetup, pi, K, false);
        auto res_piK = vtxu::fitQuality(D0KinTree, __PvalChi2Vtx_min__);
        if(!res_piK.isGood) continue;
        updateCounter(4, countersFlag);

        D0KinTree->movePointerToTheTop();
        auto mass_piK = D0KinTree->currentParticle()->currentState().mass();
        if (fabs(mass_piK - mass_D0) > __dmD0_max__) continue;
        updateCounter(5, countersFlag);

        auto D0 = D0KinTree->currentParticle();
        auto vtxD0 = D0KinTree->currentDecayVertex();

        auto cos_D0_PV = vtxu::computePointingCos(primaryVtx, vtxD0, D0);
        auto cosT_D0_PV = vtxu::computePointingCosTransverse(primaryVtx, vtxD0, D0);
        auto d_vtxD0_PV = vtxu::vtxsDistance(primaryVtx, vtxD0);
        auto sigd_vtxD0_PV = d_vtxD0_PV.first/d_vtxD0_PV.second;
        auto dxy_vtxD0_PV = vtxu::vtxsTransverseDistance(primaryVtx, vtxD0);
        auto sigdxy_vtxD0_PV = dxy_vtxD0_PV.first/dxy_vtxD0_PV.second;

        if(cosT_D0_PV < 0.95) continue;
        if(sigdxy_vtxD0_PV < __sigdxy_vtx_PV_min__) continue;
        updateCounter(6, countersFlag);

        n_D0++;

        /*
        ############################################################################
                         Look for the first pi+ to make a B+
        ############################################################################
        */
        for(uint i_pi1 = 0; i_pi1 < N_pfCand; ++i_pi1) {
          if(i_pi1==i_K || i_pi1==i_pi) continue;

          const pat::PackedCandidate & pi1 = (*pfCandHandle)[i_pi1];
          if (!pi1.hasTrackDetails()) continue;
          //Require a positive charged hadron
          if (pi1.pdgId() != charge_pi1*211 ) continue;
          //Require a minimum pt
          if(pi1.pt() < 0.8) continue;

          auto pi1_tk = pi1.bestTrack();
          if (fabs(pi1_tk->dz(primaryVtx.position()) - K_tk->dz(primaryVtx.position())) > __dzMax__) continue;
          if (vtxu::dR(pi1.phi(), K.phi(), pi1.eta(), K.eta()) > __dRMax__) continue;
          // Require significant impact parameter
          auto dxy = pi1_tk->dxy(primaryVtx.position());
          auto sigdxy_pi1_PV = fabs(dxy)/pi.dxyError();
          auto pi1_norm_chi2 = pi_tk->normalizedChi2();
          auto pi1_N_valid_hits = pi_tk->numberOfValidHits();
          if (sigdxy_pi1_PV < __sigIPpfCand_min__) continue;
          updateCounter(7, countersFlag);

          n_pi1++;

          // Fit the B+ vertex
          auto D0pi1KinTree = vtxu::FitDst(iSetup, pi1, D0, false);
          auto res_D0pi1 = vtxu::fitQuality(D0pi1KinTree, __PvalChi2Vtx_min__);
          if(!res_D0pi1.isGood) continue;
          D0pi1KinTree->movePointerToTheTop();
          auto mass_D0pi1 = D0pi1KinTree->currentParticle()->currentState().mass();
          if (mass_D0pi1 > 7.) continue;

          /*
          ############################################################################
                           Look for the second pi+ to make a B+
          ############################################################################
          */
          for(uint i_pi2 = i_pi1+1; i_pi2 < N_pfCand; ++i_pi2) {
            if(i_pi2==i_K || i_pi2==i_pi) continue;

            const pat::PackedCandidate & pi2 = (*pfCandHandle)[i_pi2];
            if (!pi2.hasTrackDetails()) continue;
            //Require a positive charged hadron
            if (pi2.pdgId() != charge_pi2*211 ) continue;
            //Require a minimum pt
            if(pi2.pt() < 0.8) continue;

            auto pi2_tk = pi2.bestTrack();
            if (fabs(pi2_tk->dz(primaryVtx.position()) - K_tk->dz(primaryVtx.position())) > __dzMax__) continue;
            if (vtxu::dR(pi2.phi(), K.phi(), pi2.eta(), K.eta()) > __dRMax__) continue;

            // Require significant impact parameter
            auto dxy = pi2_tk->dxy(primaryVtx.position());
            auto sigdxy_pi2_PV = fabs(dxy)/pi.dxyError();
            auto pi2_norm_chi2 = pi_tk->normalizedChi2();
            auto pi2_N_valid_hits = pi_tk->numberOfValidHits();
            if (sigdxy_pi2_PV < __sigIPpfCand_min__) continue;

            n_pi2++;

            // Fit the B+ vertex
            auto D0pipiKinTree = vtxu::FitD0pipi(iSetup, D0, pi1, pi2);
            auto res_D0pipi = vtxu::fitQuality(D0pipiKinTree, __PvalChi2Vtx_min__);
            if(!res_D0pipi.isGood) continue;

            D0pipiKinTree->movePointerToTheTop();
            auto D0pipi = D0pipiKinTree->currentParticle();
            auto mass_D0pipi = D0pipi->currentState().mass();
            if (mass_D0pipi < 4.) continue;
            if (mass_D0pipi > 7.) continue;

            auto vtxD0pipi = D0pipiKinTree->currentDecayVertex();

            D0pipiKinTree->movePointerToTheFirstChild();

            auto refit_D0 = D0pipiKinTree->currentParticle();
            auto p4_refit_D0 = vtxu::getTLVfromKinPart(refit_D0);
            D0pipiKinTree->movePointerToTheNextChild();

            auto refit_pi1 = D0pipiKinTree->currentParticle();
            auto p4_refit_pi1 = vtxu::getTLVfromKinPart(refit_pi1);
            D0pipiKinTree->movePointerToTheNextChild();

            auto refit_pi2 = D0pipiKinTree->currentParticle();
            auto p4_refit_pi2 = vtxu::getTLVfromKinPart(refit_pi2);

            double exp_pt = (mass_Dst/mass_D0)*p4_refit_D0.Pt();
            double exp_M = mass_piK + (mass_Dst - mass_D0);
            TLorentzVector p4_exp_Dst;
            p4_exp_Dst.SetPtEtaPhiM(exp_pt, p4_refit_D0.Eta(), p4_refit_D0.Phi(), exp_M);

            auto p4_D0pipi_scaledDst = p4_exp_Dst + p4_refit_pi1 + p4_refit_pi2;
            auto mass_D0pipi_scaledDst = p4_D0pipi_scaledDst.M();

            // Looking for the best vertex to associate it with
            uint i_best = 0;
            auto maxCos = vtxu::computePointingCos(primaryVtx, vtxD0pipi, D0pipi);
            for(uint i_vtx = 1; i_vtx < possibleVtxs.size(); i_vtx++) {
              auto auxCos = vtxu::computePointingCos(possibleVtxs[i_vtx], vtxD0pipi, D0pipi);
              if(auxCos > maxCos) {maxCos = auxCos; i_best = i_vtx;}
            }
            auto bestVtx = possibleVtxs[i_best];

            auto cos_D0pipi_PV = vtxu::computePointingCos(bestVtx, vtxD0pipi, D0pipi);
            auto cosT_D0pipi_PV = vtxu::computePointingCosTransverse(bestVtx, vtxD0pipi, D0pipi);
            auto d_vtxD0pipi_PV = vtxu::vtxsDistance(bestVtx, vtxD0pipi);
            auto sigd_vtxD0pipi_PV = d_vtxD0pipi_PV.first/d_vtxD0pipi_PV.second;
            auto dxy_vtxD0pipi_PV = vtxu::vtxsDistance(bestVtx, vtxD0pipi);
            auto sigdxy_vtxD0pipi_PV = dxy_vtxD0pipi_PV.first/dxy_vtxD0pipi_PV.second;

            if(cos_D0pipi_PV < 0.95) continue;
            n_B++;

            if (verbose) {cout << "B->D*- pi+ pi+ candidate found\n";}

            (*outputVecNtuplizer)["sigdxy_K_PV"].push_back(sigdxy_K_PV);
            (*outputVecNtuplizer)["K_norm_chi2"].push_back(K_norm_chi2);
            (*outputVecNtuplizer)["K_N_valid_hits"].push_back(K_N_valid_hits);
            AddTLVToOut(vtxu::getTLVfromCand(K, mass_K), string("K"), &(*outputVecNtuplizer));
            (*outputVecNtuplizer)["sigdxy_pi_PV"].push_back(sigdxy_pi_PV);
            (*outputVecNtuplizer)["pi_norm_chi2"].push_back(pi_norm_chi2);
            (*outputVecNtuplizer)["pi_N_valid_hits"].push_back(pi_N_valid_hits);
            AddTLVToOut(vtxu::getTLVfromCand(pi, mass_pi), string("pi"), &(*outputVecNtuplizer));
            double m = (vtxu::getTLVfromCand(pi, mass_pi) + vtxu::getTLVfromCand(K, mass_K)).M();
            (*outputVecNtuplizer)["massb_piK"].push_back(m);
            m = (vtxu::getTLVfromCand(pi, mass_K) + vtxu::getTLVfromCand(K, mass_K)).M();
            (*outputVecNtuplizer)["massb_KK"].push_back(m);
            m = (vtxu::getTLVfromCand(pi, mass_K) + vtxu::getTLVfromCand(K, mass_pi)).M();
            (*outputVecNtuplizer)["massb_Kpi_conj"].push_back(m);

            (*outputVecNtuplizer)["chi2_piK"].push_back(res_piK.chi2);
            (*outputVecNtuplizer)["dof_piK"].push_back(res_piK.dof);
            (*outputVecNtuplizer)["pval_piK"].push_back(res_piK.pval);
            (*outputVecNtuplizer)["mass_piK"].push_back(mass_piK);
            AddTLVToOut(vtxu::getTLVfromKinPart(D0), string("D0"), &(*outputVecNtuplizer));
            D0KinTree->movePointerToTheFirstChild();
            auto refit_K = D0KinTree->currentParticle();
            AddTLVToOut(vtxu::getTLVfromKinPart(refit_K), string("K_refitpiK"), &(*outputVecNtuplizer));
            D0KinTree->movePointerToTheNextChild();
            auto refit_pi = D0KinTree->currentParticle();
            AddTLVToOut(vtxu::getTLVfromKinPart(refit_pi), string("pi_refitpiK"), &(*outputVecNtuplizer));
            (*outputVecNtuplizer)["cos_D0_PV"].push_back(cos_D0_PV);
            (*outputVecNtuplizer)["cosT_D0_PV"].push_back(cosT_D0_PV);
            (*outputVecNtuplizer)["d_vtxD0_PV"].push_back(d_vtxD0_PV.first);
            (*outputVecNtuplizer)["sigd_vtxD0_PV"].push_back(sigd_vtxD0_PV);
            (*outputVecNtuplizer)["dxy_vtxD0_PV"].push_back(dxy_vtxD0_PV.first);
            (*outputVecNtuplizer)["sigdxy_vtxD0_PV"].push_back(sigdxy_vtxD0_PV);

            (*outputVecNtuplizer)["sigdxy_pi1_PV"].push_back(sigdxy_pi1_PV);
            (*outputVecNtuplizer)["pi1_norm_chi2"].push_back(pi1_norm_chi2);
            (*outputVecNtuplizer)["pi1_N_valid_hits"].push_back(pi1_N_valid_hits);
            AddTLVToOut(vtxu::getTLVfromCand(pi1, mass_pi), string("pi1"), &(*outputVecNtuplizer));

            (*outputVecNtuplizer)["sigdxy_pi2_PV"].push_back(sigdxy_pi2_PV);
            (*outputVecNtuplizer)["pi2_norm_chi2"].push_back(pi2_norm_chi2);
            (*outputVecNtuplizer)["pi2_N_valid_hits"].push_back(pi2_N_valid_hits);
            AddTLVToOut(vtxu::getTLVfromCand(pi2, mass_pi), string("pi2"), &(*outputVecNtuplizer));

            (*outputVecNtuplizer)["chi2_D0pipi"].push_back(res_D0pipi.chi2);
            (*outputVecNtuplizer)["dof_D0pipi"].push_back(res_D0pipi.dof);
            (*outputVecNtuplizer)["pval_D0pipi"].push_back(res_D0pipi.pval);
            (*outputVecNtuplizer)["mass_D0pipi"].push_back(mass_D0pipi);
            (*outputVecNtuplizer)["mass_D0pipi_scaledDst"].push_back(mass_D0pipi_scaledDst);
            (*outputVecNtuplizer)["cos_D0pipi_PV"].push_back(cos_D0pipi_PV);
            (*outputVecNtuplizer)["cosT_D0pipi_PV"].push_back(cosT_D0pipi_PV);
            (*outputVecNtuplizer)["d_vtxD0pipi_PV"].push_back(d_vtxD0pipi_PV.first);
            (*outputVecNtuplizer)["sigd_vtxD0pipi_PV"].push_back(sigd_vtxD0pipi_PV);
            (*outputVecNtuplizer)["dxy_vtxD0pipi_PV"].push_back(dxy_vtxD0pipi_PV.first);
            (*outputVecNtuplizer)["sigdxy_vtxD0pipi_PV"].push_back(sigdxy_vtxD0pipi_PV);
            AddTLVToOut(p4_D0pipi_scaledDst, string("D0pipi_sDst"), &(*outputVecNtuplizer));
            AddTLVToOut(vtxu::getTLVfromKinPart(D0pipi), string("D0pipi"), &(*outputVecNtuplizer));
            (*outputVecNtuplizer)["dphi_trgMu_D0pipi_sDst"].push_back(vtxu::dPhi(trgMu.phi(), p4_D0pipi_scaledDst.Phi()));

            /*
            ############################################################################
                             Look for the soft pion to make a Dst-
            ############################################################################
            */
            (*outputVecNtuplizer)["sigdxy_pis_PV"] = {};
            (*outputVecNtuplizer)["pis_norm_chi2"] = {};
            (*outputVecNtuplizer)["pis_N_valid_hits"] = {};
            (*outputVecNtuplizer)["pis_pt"] = {};
            (*outputVecNtuplizer)["pis_eta"] = {};
            (*outputVecNtuplizer)["pis_phi"] = {};
            (*outputVecNtuplizer)["chi2_D0pis"] = {};
            (*outputVecNtuplizer)["dof_D0pis"] = {};
            (*outputVecNtuplizer)["pval_D0pis"] = {};
            (*outputVecNtuplizer)["mass_D0pis"] = {};
            (*outputVecNtuplizer)["cos_Dst_PV"] = {};
            (*outputVecNtuplizer)["pis_refitD0pis_pt"] = {};
            (*outputVecNtuplizer)["pis_refitD0pis_eta"] = {};
            (*outputVecNtuplizer)["pis_refitD0pis_phi"] = {};
            (*outputVecNtuplizer)["chi2_D0pihpis"] = {};
            (*outputVecNtuplizer)["dof_D0pihpis"] = {};
            (*outputVecNtuplizer)["pval_D0pihpis"] = {};
            (*outputVecNtuplizer)["mass_D0pihpis"] = {};
            for(uint i_pis = 0; i_pis < N_pfCand; ++i_pis) {
              if(i_pis==i_K || i_pis==i_pi || i_pis==i_pi1 || i_pis==i_pi2) continue;

              const pat::PackedCandidate & pis = (*pfCandHandle)[i_pis];
              if (!pis.hasTrackDetails()) continue;
              if (pis.pdgId() != charge_pis*211 ) continue;

              // Require to be close to the trigger muon;
              auto pis_tk = pis.bestTrack();
              if (fabs(pis_tk->dz(primaryVtx.position()) - K_tk->dz(primaryVtx.position())) > __dzMax__) continue;
              if (vtxu::dR(pis.phi(), K.phi(), pis.eta(), K.eta()) > __dRMax__) continue;
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
              n_Dst++;
              if (verbose) {cout << "D* found\n";}

              auto Dst = DstKinTree->currentParticle();
              auto vtxDst = DstKinTree->currentDecayVertex();

              auto cos_Dst_PV = vtxu::computePointingCos(primaryVtx, vtxDst, Dst);

              (*outputVecNtuplizer)["sigdxy_pis_PV"].push_back(sigdxy_pis_PV);
              (*outputVecNtuplizer)["pis_norm_chi2"].push_back(pis_norm_chi2);
              (*outputVecNtuplizer)["pis_N_valid_hits"].push_back(pis_N_valid_hits);
              AddTLVToOut(vtxu::getTLVfromCand(pis, mass_pi), string("pis"), &(*outputVecNtuplizer));

              (*outputVecNtuplizer)["chi2_D0pis"].push_back(res_D0pis.chi2);
              (*outputVecNtuplizer)["dof_D0pis"].push_back(res_D0pis.dof);
              (*outputVecNtuplizer)["pval_D0pis"].push_back(res_D0pis.pval);
              (*outputVecNtuplizer)["mass_D0pis"].push_back(mass_D0pis);
              (*outputVecNtuplizer)["cos_Dst_PV"].push_back(cos_Dst_PV);

              DstKinTree->movePointerToTheFirstChild();
              auto refitD0pis_pis = DstKinTree->currentParticle();
              AddTLVToOut(vtxu::getTLVfromKinPart(refitD0pis_pis), string("pis_refitD0pis"), &(*outputVecNtuplizer));


              auto fullBKinTree = vtxu::FitD0pipipis(iSetup, D0, pi1, pi2, pis);
              auto res_D0pipipis = vtxu::fitQuality(fullBKinTree, __PvalChi2Vtx_min__);
              if(res_D0pipipis.isGood){
                n_fullB++;
                (*outputVecNtuplizer)["chi2_D0pihpis"].push_back(res_D0pipipis.chi2);
                (*outputVecNtuplizer)["dof_D0pihpis"].push_back(res_D0pipipis.dof);
                (*outputVecNtuplizer)["pval_D0pihpis"].push_back(res_D0pipipis.pval);

                fullBKinTree->movePointerToTheTop();
                auto mass_D0pipipis = fullBKinTree->currentParticle()->currentState().mass();
                (*outputVecNtuplizer)["mass_D0pipipis"].push_back(mass_D0pipipis);
              }
              else{
                (*outputVecNtuplizer)["chi2_D0pipipis"].push_back(-1);
                (*outputVecNtuplizer)["dof_D0pipipis"].push_back(-1);
                (*outputVecNtuplizer)["pval_D0pipipis"].push_back(-1);
                (*outputVecNtuplizer)["mass_D0pipipis"].push_back(-1);
              }
            }

            if(n_B > 30) break;
          }
          if(n_B > 30) break;
        }
        if(n_B > 30) break;
      }
      if(n_B > 30) break;
    }

    (*outputNtuplizer)["Run"] = iEvent.run();
    (*outputNtuplizer)["LumiBlock"] = iEvent.luminosityBlock();
    (*outputNtuplizer)["eventNumber"] = iEvent.id().event();

    (*outputNtuplizer)["n_K"] = n_K;
    (*outputNtuplizer)["n_pi"] = n_pi;
    (*outputNtuplizer)["n_D0"] = n_D0;
    (*outputNtuplizer)["n_pi1"] = n_pi1;
    (*outputNtuplizer)["n_pi2"] = n_pi2;
    (*outputNtuplizer)["n_B"] = n_B;
    (*outputNtuplizer)["n_pis"] = n_pis;
    (*outputNtuplizer)["n_Dst"] = n_Dst;
    (*outputNtuplizer)["n_fullB"] = n_fullB;

    iEvent.put(move(outputNtuplizer), "outputNtuplizer");
    iEvent.put(move(outputVecNtuplizer), "outputVecNtuplizer");
    return;
}


void TagAndProbeBp2DmstPipPipProducer::AddTLVToOut(TLorentzVector v, string n, map<string, vector<float>>* outv) {
  (*outv)[n+"_pt"].push_back(v.Pt());
  (*outv)[n+"_eta"].push_back(v.Eta());
  (*outv)[n+"_phi"].push_back(v.Phi());
  return;
}

bool TagAndProbeBp2DmstPipPipProducer::qualityMuonID(pat::Muon m, reco::Vertex pVtx) {
  if(m.innerTrack().isNull()) return false;

  if (!m.isSoftMuon(pVtx)) return false;
  // if(m.innerTrack()->hitPattern().pixelLayersWithMeasurement() < 2) return false;
  // if(!m.innerTrack()->quality(reco::TrackBase::highPurity)) return false;
  // if(!m.isGood("TMOneStationTight")) return false;
  if(m.innerTrack()->normalizedChi2() > 1.8) return false;

  double dxy = m.innerTrack()->dxy(pVtx.position());
  float sigdxy = fabs(dxy)/m.innerTrack()->dxyError();
  if (sigdxy < 2) return false;

  return true;
}

void TagAndProbeBp2DmstPipPipProducer::updateCounter(int idx, vector<bool> &cF) {
  if(!cF[idx]) {
    cF[idx] = true;
    counters[idx]++;
  }
  return;
}

DEFINE_FWK_MODULE(TagAndProbeBp2DmstPipPipProducer);
