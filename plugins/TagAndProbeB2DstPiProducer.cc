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
#define __sigdxy_vtx_PV_min__ 2.0 // loose cut
#define __dmDst_max__ 0.050 // loose cut
#define __mass_D0pismu_max__ 10.0 // Some reasonable cut on the mass
#define __pTaddTracks_min__ 0.3 // loose cut
#define __mass_D0pismupi_max__ 10.0 // Some reasonable cut on the mass


using namespace std;

class TagAndProbeB2DstPiProducer : public edm::EDProducer {

public:

    explicit TagAndProbeB2DstPiProducer(const edm::ParameterSet &iConfig);
    void AddTLVToOut(TLorentzVector, string, map<string, vector<float>>*);
    bool qualityMuonID(pat::Muon, reco::Vertex);
    void updateCounter(int, vector<bool>&);

    ~TagAndProbeB2DstPiProducer() {
      cout << Form("TagAndProbeB2DstPiProducer counters:%d\n", (int)counters.size());
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
    int charge_pih = +1;
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



TagAndProbeB2DstPiProducer::TagAndProbeB2DstPiProducer(const edm::ParameterSet &iConfig)
{
    PFCandSrc_ = consumes<vector<pat::PackedCandidate>>(edm::InputTag("packedPFCandidates"));
    vtxSrc_ = consumes<vector<reco::Vertex>>(edm::InputTag("offlineSlimmedPrimaryVertices"));
    TrgMuonSrc_ = consumes<vector<pat::Muon>> ( iConfig.getParameter<edm::InputTag>( "trgMuons" ) );

    charge_K = iConfig.getParameter<int>( "charge_K" );
    charge_pi = iConfig.getParameter<int>( "charge_pi" );
    charge_pih = iConfig.getParameter<int>( "charge_pih" );
    charge_pis = iConfig.getParameter<int>( "charge_pis" );

    verbose = iConfig.getParameter<int>( "verbose" );

    produces<map<string, float>>("outputNtuplizer");
    produces<map<string, vector<float>>>("outputVecNtuplizer");
    for(uint i=0; i<10; i++) counters.push_back(0);
}


void TagAndProbeB2DstPiProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
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

    int n_K = 0, n_pi = 0, n_D0 = 0, n_pih = 0, n_B = 0;
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

        if(cosT_D0_PV < 0.9) continue;
        if(sigdxy_vtxD0_PV < __sigdxy_vtx_PV_min__) continue;
        updateCounter(6, countersFlag);

        n_D0++;

        /*
        ############################################################################
                         Look for the hard pi+ to make a B0
        ############################################################################
        */
        for(uint i_pih = 0; i_pih < N_pfCand; ++i_pih) {
          if(i_pih==i_K || i_pih==i_pi) continue;

          const pat::PackedCandidate & pih = (*pfCandHandle)[i_pih];
          if (!pih.hasTrackDetails()) continue;
          //Require a positive charged hadron
          if (pih.pdgId() != charge_pih*211 ) continue;
          //Require a minimum pt
          if(pih.pt() < 6.5) continue;

          auto pih_tk = pih.bestTrack();
          if (fabs(pih_tk->dz(primaryVtx.position()) - K_tk->dz(primaryVtx.position())) > __dzMax__) continue;
          if (vtxu::dR(pih.phi(), K.phi(), pih.eta(), K.eta()) > __dRMax__) continue;
          // Require significant impact parameter
          auto dxy = pih_tk->dxy(primaryVtx.position());
          auto sigdxy_pih_PV = fabs(dxy)/pi.dxyError();
          auto pih_norm_chi2 = pi_tk->normalizedChi2();
          auto pih_N_valid_hits = pi_tk->numberOfValidHits();
          if (sigdxy_pih_PV < __sigIPpfCand_min__) continue;
          updateCounter(7, countersFlag);

          n_pih++;

          // Fit the B0 vertex
          auto D0pihKinTree = vtxu::FitDst(iSetup, pih, D0, false);
          auto res_D0pih = vtxu::fitQuality(D0pihKinTree, __PvalChi2Vtx_min__);
          if(!res_D0pih.isGood) continue;
          updateCounter(8, countersFlag);

          D0pihKinTree->movePointerToTheFirstChild();
          auto refit_pih = D0pihKinTree->currentParticle();
          auto p4_refit_pih = vtxu::getTLVfromKinPart(refit_pih);
          D0pihKinTree->movePointerToTheNextChild();
          auto refit_D0 = D0pihKinTree->currentParticle();
          auto p4_refit_D0 = vtxu::getTLVfromKinPart(refit_D0);

          double exp_pt = (mass_Dst/mass_D0)*p4_refit_D0.Pt();
          double exp_M = mass_piK + (mass_Dst - mass_D0);
          TLorentzVector p4_exp_Dst;
          p4_exp_Dst.SetPtEtaPhiM(exp_pt, p4_refit_D0.Eta(), p4_refit_D0.Phi(), exp_M);

          auto p4_D0pih_scaledDst = p4_exp_Dst + p4_refit_pih;
          auto mass_D0pih_scaledDst = p4_D0pih_scaledDst.M();
          if (mass_D0pih_scaledDst < 4.) continue;
          if (mass_D0pih_scaledDst > 7.) continue;
          updateCounter(9, countersFlag);

          D0pihKinTree->movePointerToTheTop();
          auto D0pih = D0pihKinTree->currentParticle();
          auto mass_D0pih = D0pih->currentState().mass();
          auto vtxD0pih = D0pihKinTree->currentDecayVertex();

          // Looking for the best vertex to associate it with
          uint i_best = 0;
          auto maxCos = vtxu::computePointingCos(primaryVtx, vtxD0pih, D0pih);
          for(uint i_vtx = 1; i_vtx < possibleVtxs.size(); i_vtx++) {
            auto auxCos = vtxu::computePointingCos(possibleVtxs[i_vtx], vtxD0pih, D0pih);
            if(auxCos > maxCos) {maxCos = auxCos; i_best = i_vtx;}
          }
          auto bestVtx = possibleVtxs[i_best];

          auto cos_D0pih_PV = vtxu::computePointingCos(bestVtx, vtxD0pih, D0pih);
          auto cosT_D0pih_PV = vtxu::computePointingCosTransverse(bestVtx, vtxD0pih, D0pih);
          auto d_vtxD0pih_PV = vtxu::vtxsDistance(bestVtx, vtxD0pih);
          auto sigd_vtxD0pih_PV = d_vtxD0pih_PV.first/d_vtxD0pih_PV.second;
          auto dxy_vtxD0pih_PV = vtxu::vtxsDistance(bestVtx, vtxD0pih);
          auto sigdxy_vtxD0pih_PV = dxy_vtxD0pih_PV.first/dxy_vtxD0pih_PV.second;

          if(cos_D0pih_PV < 0.9) continue;
          n_B++;
          if (verbose) {cout << "B->D* pi candidate found\n";}
          updateCounter(10, countersFlag);

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

          (*outputVecNtuplizer)["sigdxy_pih_PV"].push_back(sigdxy_pih_PV);
          (*outputVecNtuplizer)["pih_norm_chi2"].push_back(pih_norm_chi2);
          (*outputVecNtuplizer)["pih_N_valid_hits"].push_back(pih_N_valid_hits);
          AddTLVToOut(vtxu::getTLVfromCand(pih, mass_pi), string("pih"), &(*outputVecNtuplizer));

          (*outputVecNtuplizer)["chi2_D0pih"].push_back(res_D0pih.chi2);
          (*outputVecNtuplizer)["dof_D0pih"].push_back(res_D0pih.dof);
          (*outputVecNtuplizer)["pval_D0pih"].push_back(res_D0pih.pval);
          (*outputVecNtuplizer)["mass_D0pih_scaledDst"].push_back(mass_D0pih_scaledDst);
          (*outputVecNtuplizer)["mass_D0pih"].push_back(mass_D0pih);
          (*outputVecNtuplizer)["cos_D0pih_PV"].push_back(cos_D0pih_PV);
          (*outputVecNtuplizer)["cosT_D0pih_PV"].push_back(cosT_D0pih_PV);
          (*outputVecNtuplizer)["d_vtxD0pih_PV"].push_back(d_vtxD0pih_PV.first);
          (*outputVecNtuplizer)["sigd_vtxD0pih_PV"].push_back(sigd_vtxD0pih_PV);
          (*outputVecNtuplizer)["dxy_vtxD0pih_PV"].push_back(dxy_vtxD0pih_PV.first);
          (*outputVecNtuplizer)["sigdxy_vtxD0pih_PV"].push_back(sigdxy_vtxD0pih_PV);
          AddTLVToOut(p4_D0pih_scaledDst, string("D0pih_sDst"), &(*outputVecNtuplizer));
          AddTLVToOut(vtxu::getTLVfromKinPart(D0pih), string("D0pih"), &(*outputVecNtuplizer));
          (*outputVecNtuplizer)["dphi_trgMu_D0pih_sDst"].push_back(vtxu::dPhi(trgMu.phi(), p4_D0pih_scaledDst.Phi()));

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
            if(i_pis==i_K || i_pis==i_pi || i_pis==i_pih) continue;

            const pat::PackedCandidate & pis = (*pfCandHandle)[i_pis];
            if (!pis.hasTrackDetails()) continue;
            //Require a negative charged hadron
            if (pis.pdgId() != -211 ) continue;

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


            auto fullBKinTree = vtxu::Fit_D0pihpis(iSetup, D0, pih, pis);
            auto res_D0pihpis = vtxu::fitQuality(fullBKinTree, __PvalChi2Vtx_min__);
            if(res_D0pihpis.isGood){
              n_fullB++;
              (*outputVecNtuplizer)["chi2_D0pihpis"].push_back(res_D0pihpis.chi2);
              (*outputVecNtuplizer)["dof_D0pihpis"].push_back(res_D0pihpis.dof);
              (*outputVecNtuplizer)["pval_D0pihpis"].push_back(res_D0pihpis.pval);

              fullBKinTree->movePointerToTheTop();
              auto mass_D0pihpis = fullBKinTree->currentParticle()->currentState().mass();
              (*outputVecNtuplizer)["mass_D0pihpis"].push_back(mass_D0pihpis);
            }
            else{
              (*outputVecNtuplizer)["chi2_D0pihpis"].push_back(-1);
              (*outputVecNtuplizer)["dof_D0pihpis"].push_back(-1);
              (*outputVecNtuplizer)["pval_D0pihpis"].push_back(-1);
              (*outputVecNtuplizer)["mass_D0pihpis"].push_back(-1);
            }
          }

          if(n_B > 50) break;
        }
        if(n_B > 50) break;
      }
      if(n_B > 50) break;
    }

    (*outputNtuplizer)["Run"] = iEvent.run();
    (*outputNtuplizer)["LumiBlock"] = iEvent.luminosityBlock();
    (*outputNtuplizer)["eventNumber"] = iEvent.id().event();

    (*outputNtuplizer)["n_K"] = n_K;
    (*outputNtuplizer)["n_pi"] = n_pi;
    (*outputNtuplizer)["n_D0"] = n_D0;
    (*outputNtuplizer)["n_pih"] = n_pih;
    (*outputNtuplizer)["n_B"] = n_B;
    (*outputNtuplizer)["n_pis"] = n_pis;
    (*outputNtuplizer)["n_Dst"] = n_Dst;
    (*outputNtuplizer)["n_fullB"] = n_fullB;

    iEvent.put(move(outputNtuplizer), "outputNtuplizer");
    iEvent.put(move(outputVecNtuplizer), "outputVecNtuplizer");
    return;
}


void TagAndProbeB2DstPiProducer::AddTLVToOut(TLorentzVector v, string n, map<string, vector<float>>* outv) {
  (*outv)[n+"_pt"].push_back(v.Pt());
  (*outv)[n+"_eta"].push_back(v.Eta());
  (*outv)[n+"_phi"].push_back(v.Phi());
  return;
}

bool TagAndProbeB2DstPiProducer::qualityMuonID(pat::Muon m, reco::Vertex pVtx) {
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

void TagAndProbeB2DstPiProducer::updateCounter(int idx, vector<bool> &cF) {
  if(!cF[idx]) {
    cF[idx] = true;
    counters[idx]++;
  }
  return;
}

DEFINE_FWK_MODULE(TagAndProbeB2DstPiProducer);
