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

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <iostream>
#include <string>
#include <map>
#include <TMath.h>

#include "VtxUtils.hh"

#define __pThad_min__ 0.5 // loose cut
#define __dzMax__ 1.0
#define __dRMax__ 2.0
#define __sigIPpfCand_min__ 2. // loose cut
#define __PvalChi2Vtx_min__ 0.05 // loose cut
#define __dmD0_max__ 0.15 // loose cut
#define __sigdxy_vtx_PV_min__ 4.0 // loose cut
#define __dmDst_max__ 0.15 // loose cut
#define __mass_D0pipmu_max__ 10. // Some reasonable cut on the mass
#define __pTaddTracks_min__ 0.3 // loose cut
// #define __mass_D0pismupi_max__ 10. // Some reasonable cut on the mass


using namespace std;

class TagAndProbeBp2DststMuProducer : public edm::EDProducer {

public:

    explicit TagAndProbeBp2DststMuProducer(const edm::ParameterSet &iConfig);
    void AddTLVToOut(TLorentzVector, string, map<string, vector<float>>*);
    bool qualityMuonID(pat::Muon, reco::Vertex);
    void updateCounter(int, vector<bool>&);
    double Mass_varM(TLorentzVector, double, TLorentzVector, double);

    ~TagAndProbeBp2DststMuProducer() {
      cout << Form("Total number of fit crashed %u", fitCrash) << endl;
      cout << Form("TagAndProbeBp2DststMuProducer counters:%d\n", (int)counters.size());
      for(auto v : counters) {
        cout << v << endl;
      }
      cout << "\n\n";

      edm::Service<TFileService> fs;
      TH1I* hCounters;
      int nC = counters.size();
      hCounters = fs->make<TH1I>("hB2DstMuDecayCounters", "hB2DstMuDecayCounters", nC, 0, nC);
      for(int i=0; i < nC; i++) hCounters->SetBinContent(i+1, counters[i]);
    };

private:
    virtual void produce(edm::Event&, const edm::EventSetup&);

    // ----------member data ---------------------------
    edm::EDGetTokenT<vector<pat::PackedCandidate>> PFCandSrc_;
    edm::EDGetTokenT<vector<reco::Vertex>> vtxSrc_;
    edm::EDGetTokenT<vector<pat::Muon>> TrgMuonSrc_;

    int charge_muon = +1;
    int charge_K = +1;
    int charge_pi = -1;
    int charge_pip = +1;

    double mass_Mu  = 0.105658;
    double mass_pi  = 0.139570;
    double mass_K   = 0.493677;
    double mass_D0  = 1.86483;
    double mass_Dst = 2.01026;
    double mass_B0  = 5.27963;

    int verbose = 0;

    vector<uint> counters;
    uint fitCrash = 0;
};



TagAndProbeBp2DststMuProducer::TagAndProbeBp2DststMuProducer(const edm::ParameterSet &iConfig)
{
    PFCandSrc_ = consumes<vector<pat::PackedCandidate>>(edm::InputTag("packedPFCandidates"));
    vtxSrc_ = consumes<vector<reco::Vertex>>(edm::InputTag("offlineSlimmedPrimaryVertices"));
    TrgMuonSrc_ = consumes<vector<pat::Muon>> ( iConfig.getParameter<edm::InputTag>( "trgMuons" ) );

    charge_muon = iConfig.getParameter<int>( "charge_muon" );
    charge_K = iConfig.getParameter<int>( "charge_K" );
    charge_pi = iConfig.getParameter<int>( "charge_pi" );
    charge_pip = iConfig.getParameter<int>( "charge_pip" );

    verbose = iConfig.getParameter<int>( "verbose" );

    produces<map<string, float>>("outputNtuplizer");
    produces<map<string, vector<float>>>("outputVecNtuplizer");
    for(uint i=0; i<13; i++) counters.push_back(0);
}


void TagAndProbeBp2DststMuProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    edm::Handle<vector<pat::PackedCandidate>> pfCandHandle;
    iEvent.getByToken(PFCandSrc_, pfCandHandle);
    unsigned int N_pfCand = pfCandHandle->size();

    edm::Handle<vector<reco::Vertex>> vtxHandle;
    iEvent.getByToken(vtxSrc_, vtxHandle);
    auto primaryVtx = (*vtxHandle)[0];

    edm::Handle<vector<pat::Muon>> trgMuonsHandle;
    iEvent.getByToken(TrgMuonSrc_, trgMuonsHandle);

    // Output collection
    unique_ptr<map<string, float>> outputNtuplizer(new map<string, float>);
    unique_ptr<map<string, vector<float>>> outputVecNtuplizer(new map<string, vector<float>>);

    int n_mu = 0, n_K = 0, n_pi = 0, n_D0 = 0, n_pip = 0, n_Dst = 0, n_B = 0;

    if (verbose) {cout <<"-------------------- Evt -----------------------\n";}
    vector<bool> countersFlag(counters.size(), false);
    counters[0]++;
    for(uint i_trgMu = 0; i_trgMu < trgMuonsHandle->size(); i_trgMu++) {
      //######## Require muon quality ##################
      auto trgMu = (*trgMuonsHandle)[i_trgMu];
      if (trgMu.charge() != charge_muon) continue;
      if (trgMu.innerTrack().isNull()) continue;
      if (fabs(trgMu.eta()) > 1.5) continue;
      updateCounter(1, countersFlag);

      vector<reco::Vertex> possibleVtxs = {};
      // bool softMuForPV = false;
      for(uint i_vtx = 0; i_vtx<vtxHandle->size(); i_vtx++) {
        auto vtx = (*vtxHandle)[i_vtx];
        if(vtx.normalizedChi2() > 1.5) continue;
        bool isSoft = trgMu.isSoftMuon(vtx);
        if(isSoft){
          possibleVtxs.push_back(vtx);
          // if (i_vtx==0) softMuForPV = true;
        }
      }
      if (possibleVtxs.size() == 0) continue;
      updateCounter(2, countersFlag);

      // if (!trgMu.isSoftMuon(primaryVtx)) continue;
      n_mu++;
      /*
      ############################################################################
                                Look for the K+ (321)
      ############################################################################
      */
      for(uint i_K = 0; i_K < N_pfCand; ++i_K) {
        const pat::PackedCandidate & K = (*pfCandHandle)[i_K];
        if (!K.hasTrackDetails()) continue;
        if (K.pdgId() != charge_K*211 ) continue;
        if (K.pt() < __pThad_min__) continue;
        if (fabs(K.eta()) > 2.4) continue;
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
        updateCounter(3, countersFlag);
        if (verbose) {cout << Form("K%s found at %u\n", charge_K>0?"+":"-", i_K);}
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
          if (pi.pdgId() != charge_pi*211 ) continue;
          if (fabs(pi.eta()) > 2.4) continue;
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
          updateCounter(4, countersFlag);
          if (verbose) {cout << Form("pi%s found at %u\n", charge_pi>0?"+":"-", i_pi);}
          n_pi++;

          //Fit the vertex
          auto D0KinTree = vtxu::FitD0(iSetup, pi, K, false);
          auto res_piK = vtxu::fitQuality(D0KinTree, __PvalChi2Vtx_min__);
          if(!res_piK.isGood) continue;
          if (verbose) {cout << "pi-K vertex fit good\n";}
          updateCounter(5, countersFlag);

          D0KinTree->movePointerToTheTop();
          auto mass_piK = D0KinTree->currentParticle()->currentState().mass();
          if (fabs(mass_piK - mass_D0) > __dmD0_max__) continue;
          if (verbose) {cout << "pi-K mass cut passed\n";}
          updateCounter(6, countersFlag);

          auto D0 = D0KinTree->currentParticle();
          auto vtxD0 = D0KinTree->currentDecayVertex();

          auto cos_D0_PV = vtxu::computePointingCos(primaryVtx, vtxD0, D0);
          auto cosT_D0_PV = vtxu::computePointingCosTransverse(primaryVtx, vtxD0, D0);
          auto d_vtxD0_PV = vtxu::vtxsDistance(primaryVtx, vtxD0);
          auto sigd_vtxD0_PV = d_vtxD0_PV.first/d_vtxD0_PV.second;
          auto dxy_vtxD0_PV = vtxu::vtxsTransverseDistance(primaryVtx, vtxD0);
          auto sigdxy_vtxD0_PV = dxy_vtxD0_PV.first/dxy_vtxD0_PV.second;

          if(cosT_D0_PV < 0.5) continue;
          if(sigdxy_vtxD0_PV < __sigdxy_vtx_PV_min__) continue;
          if (verbose) {cout << "pi-K vertex displacement passed\n";}
          updateCounter(7, countersFlag);
          if (verbose) {cout << "D0 candidate\n";}
          n_D0++;

          /*
          ############################################################################
                           Look for the soft pion to make a D**0
          ############################################################################
          */
          for(uint i_pip = 0; i_pip < N_pfCand; ++i_pip) {
            if(i_pip==i_K || i_pip==i_pi) continue;

            const pat::PackedCandidate & pip = (*pfCandHandle)[i_pip];
            if (!pip.hasTrackDetails()) continue;
            if (pip.pdgId() != charge_pip*211 ) continue;
            if (fabs(pip.eta()) > 2.4) continue;
            // Require to be close to the trigger muon;
            auto pip_tk = pip.bestTrack();
            if (fabs(pip_tk->dz(primaryVtx.position()) - trgMu.muonBestTrack()->dz(primaryVtx.position())) > __dzMax__) continue;
            if (vtxu::dR(pip.phi(), trgMu.phi(), pip.eta(), trgMu.eta()) > __dRMax__) continue;
            // Require significant impact parameter
            auto dxy = pip_tk->dxy(primaryVtx.position());
            auto sigdxy_pip_PV = fabs(dxy)/pip.dxyError();
            auto pip_norm_chi2 = pip_tk->normalizedChi2();
            auto pip_N_valid_hits = pip_tk->numberOfValidHits();
            if (sigdxy_pip_PV < __sigIPpfCand_min__) continue;
            updateCounter(8, countersFlag);
            if (verbose) {cout << Form("pip%s found at %u\n", charge_pip>0?"+":"-", i_pip);}
            n_pip++;

            // Vertex fitting from B and Dst products
            auto BKinTree = vtxu::FitB_D0pismu(iSetup, D0, pip, trgMu);
            auto res = vtxu::fitQuality(BKinTree, __PvalChi2Vtx_min__);
            if(!res.isGood) continue;
            updateCounter(11, countersFlag);

            BKinTree->movePointerToTheTop();
            auto mass_D0pipmu = BKinTree->currentParticle()->currentState().mass();
            if (mass_D0pipmu > __mass_D0pipmu_max__) continue;
            updateCounter(12, countersFlag);
            if (verbose) {cout << "Candidate B->D** mu found\n";}

            auto D0pipmu = BKinTree->currentParticle();
            auto vtxB = BKinTree->currentDecayVertex();

            BKinTree->movePointerToTheFirstChild();
            auto refit_D0 = BKinTree->currentParticle();
            auto p4_refit_D0 = vtxu::getTLVfromKinPart(refit_D0);
            AddTLVToOut(p4_refit_D0, string("D0_refitD0pipmu"), &(*outputVecNtuplizer));
            BKinTree->movePointerToTheNextChild();
            auto refit_pip = BKinTree->currentParticle();
            auto p4_refit_pip = vtxu::getTLVfromKinPart(refit_pip);
            AddTLVToOut(p4_refit_pip, string("pip_refitD0pipmu"), &(*outputVecNtuplizer));
            BKinTree->movePointerToTheNextChild();
            auto refit_mu = BKinTree->currentParticle();
            AddTLVToOut(vtxu::getTLVfromKinPart(refit_mu), string("mu_refitD0pipmu"), &(*outputVecNtuplizer));

            double exp_pt = p4_refit_D0.Pt()*(mass_Dst/mass_D0);
            double exp_M = mass_piK + (mass_Dst - mass_D0);
            TLorentzVector p4_exp_Dst;
            p4_exp_Dst.SetPtEtaPhiM(exp_pt, p4_refit_D0.Eta(), p4_refit_D0.Phi(), exp_M);
            auto p4_exp_pis = p4_exp_Dst - p4_refit_D0;
            auto p4_exp_Dstst = p4_exp_Dst + p4_refit_pip;

            auto mass_D0pip = (p4_refit_D0 + p4_refit_pip).M();
            auto mass_expDstpip = (p4_exp_Dst + p4_refit_pip).M();
            if (mass_expDstpip > 3.) continue;

            TLorentzVector p4st_pip(p4_refit_pip);
            p4st_pip.Boost(-1*p4_exp_Dstst.BoostVector());
            (*outputVecNtuplizer)["Est_pip"].push_back(p4st_pip.E());
            auto CosThetaSt_pip = cos(p4st_pip.Angle(p4_exp_Dstst.BoostVector()));
            (*outputVecNtuplizer)["CosThetaSt_pip"].push_back(CosThetaSt_pip);

            (*outputVecNtuplizer)["mass_D0pip"].push_back(mass_D0pip);
            (*outputVecNtuplizer)["mass_expDstpip"].push_back(mass_expDstpip);
            (*outputVecNtuplizer)["dm_expDstpip_pik"].push_back(mass_expDstpip-mass_piK);

            (*outputVecNtuplizer)["chi2_D0pipmu"].push_back(res.chi2);
            (*outputVecNtuplizer)["dof_D0pipmu"].push_back(res.dof);
            (*outputVecNtuplizer)["pval_D0pipmu"].push_back(res.pval);
            (*outputVecNtuplizer)["mass_D0pipmu"].push_back(mass_D0pipmu);

            AddTLVToOut(p4_exp_Dst+p4_refit_pip, string("Dstst_expD0pip"), &(*outputVecNtuplizer));
            AddTLVToOut(p4_exp_pis, string("pis_expD0pipmu"), &(*outputVecNtuplizer));


            // Looking for the best vertex to associate it with
            uint i_best = 0;
            auto maxCos = vtxu::computePointingCos(possibleVtxs[0], vtxB, D0pipmu);
            for(uint i_vtx = 1; i_vtx < possibleVtxs.size(); i_vtx++) {
              auto auxCos = vtxu::computePointingCos(possibleVtxs[i_vtx], vtxB, D0pipmu);
              if(auxCos > maxCos) {
                maxCos = auxCos;
                i_best = i_vtx;
              }
            }
            auto bestVtx = possibleVtxs[i_best];

            auto cos_D0pipmu_PV = vtxu::computePointingCos(bestVtx, vtxB, D0pipmu);
            auto cosT_D0pipmu_PV = vtxu::computePointingCosTransverse(bestVtx, vtxB, D0pipmu);
            auto d_vtxD0pipmu_PV = vtxu::vtxsDistance(bestVtx, vtxB);
            auto sigd_vtxD0pipmu_PV = d_vtxD0pipmu_PV.first/d_vtxD0pipmu_PV.second;
            auto dxy_vtxD0pipmu_PV = vtxu::vtxsDistance(bestVtx, vtxB);
            auto sigdxy_vtxD0pipmu_PV = dxy_vtxD0pipmu_PV.first/dxy_vtxD0pipmu_PV.second;
            (*outputVecNtuplizer)["cos_D0pipmu_PV"].push_back(cos_D0pipmu_PV);
            (*outputVecNtuplizer)["cosT_D0pipmu_PV"].push_back(cosT_D0pipmu_PV);
            (*outputVecNtuplizer)["d_vtxD0pipmu_PV"].push_back(d_vtxD0pipmu_PV.first);
            (*outputVecNtuplizer)["sigd_vtxD0pipmu_PV"].push_back(sigd_vtxD0pipmu_PV);
            (*outputVecNtuplizer)["dxy_vtxD0pipmu_PV"].push_back(dxy_vtxD0pipmu_PV.first);
            (*outputVecNtuplizer)["sigdxy_vtxD0pipmu_PV"].push_back(sigdxy_vtxD0pipmu_PV);

            float localVertexDensity = 0;
            for(auto vtx : (*vtxHandle)) {
              float dz = bestVtx.position().z() - vtx.position().z();
              if(fabs(dz) < 1.5*__dzMax__) localVertexDensity++;
            }
            (*outputVecNtuplizer)["localVertexDensity"].push_back(localVertexDensity/(2*1.5*__dzMax__));

            n_B++;

            (*outputVecNtuplizer)["mu_trgMu_idx"].push_back(i_trgMu);
            AddTLVToOut(vtxu::getTLVfromMuon(trgMu, mass_Mu), string("mu"), &(*outputVecNtuplizer));

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
            auto refit_K = vtxu::getTLVfromKinPart(D0KinTree->currentParticle());
            AddTLVToOut(refit_K, string("K_refitpiK"), &(*outputVecNtuplizer));
            D0KinTree->movePointerToTheNextChild();
            auto refit_pi = vtxu::getTLVfromKinPart(D0KinTree->currentParticle());
            AddTLVToOut(refit_pi, string("pi_refitpiK"), &(*outputVecNtuplizer));
            (*outputVecNtuplizer)["mass_piK_hKK"].push_back(Mass_varM(refit_pi, mass_K, refit_K, mass_K));
            (*outputVecNtuplizer)["mass_piK_hpipi"].push_back(Mass_varM(refit_pi, mass_pi, refit_K, mass_pi));
            (*outputVecNtuplizer)["mass_piK_hKpi"].push_back(Mass_varM(refit_pi, mass_K, refit_K, mass_pi));

            (*outputVecNtuplizer)["cos_D0_PV"].push_back(cos_D0_PV);
            (*outputVecNtuplizer)["cosT_D0_PV"].push_back(cosT_D0_PV);
            (*outputVecNtuplizer)["d_vtxD0_PV"].push_back(d_vtxD0_PV.first);
            (*outputVecNtuplizer)["sigd_vtxD0_PV"].push_back(sigd_vtxD0_PV);
            (*outputVecNtuplizer)["dxy_vtxD0_PV"].push_back(dxy_vtxD0_PV.first);
            (*outputVecNtuplizer)["sigdxy_vtxD0_PV"].push_back(sigdxy_vtxD0_PV);

            (*outputVecNtuplizer)["sigdxy_pip_PV"].push_back(sigdxy_pip_PV);
            (*outputVecNtuplizer)["pip_norm_chi2"].push_back(pip_norm_chi2);
            (*outputVecNtuplizer)["pip_N_valid_hits"].push_back(pip_N_valid_hits);
            AddTLVToOut(vtxu::getTLVfromCand(pip, mass_pi), string("pip"), &(*outputVecNtuplizer));

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
            (*outputVecNtuplizer)["chi2_D0pismu"] = {};
            (*outputVecNtuplizer)["dof_D0pis"] = {};
            (*outputVecNtuplizer)["dof_D0pismu"] = {};
            (*outputVecNtuplizer)["pval_D0pis"] = {};
            (*outputVecNtuplizer)["pval_D0pismu"] = {};
            (*outputVecNtuplizer)["mass_D0pis"] = {};
            (*outputVecNtuplizer)["mass_D0pismu"] = {};
            (*outputVecNtuplizer)["dm_D0pis_piK"] = {};
            (*outputVecNtuplizer)["cos_Dst_PV"] = {};
            (*outputVecNtuplizer)["pis_refitD0pismu_pt"] = {};
            (*outputVecNtuplizer)["pis_refitD0pismu_eta"] = {};
            (*outputVecNtuplizer)["pis_refitD0pismu_phi"] = {};
            for(uint i_pis = 0; i_pis < N_pfCand; ++i_pis) {
              if(i_pis==i_pi || i_pis==i_K || i_pis==i_pip) continue;

              const pat::PackedCandidate & pis = (*pfCandHandle)[i_pis];
              if (!pis.hasTrackDetails()) continue;
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

              // Fit the Dst vertex
              auto DstKinTree = vtxu::FitDst(iSetup, pis, D0, false);
              auto res_D0pis = vtxu::fitQuality(DstKinTree, __PvalChi2Vtx_min__);
              if(!res_D0pis.isGood) continue;

              DstKinTree->movePointerToTheTop();
              auto mass_D0pis = DstKinTree->currentParticle()->currentState().mass();
              if (fabs(mass_D0pis - mass_Dst) > __dmDst_max__) continue;

              if (verbose) {cout << "D* found\n";}

              auto BKinTree = vtxu::FitB_D0pismu(iSetup, D0, pis, trgMu);
              auto res = vtxu::fitQuality(BKinTree, __PvalChi2Vtx_min__);
              if(!res.isGood) continue;
              updateCounter(11, countersFlag);

              BKinTree->movePointerToTheTop();
              auto mass_D0pismu = BKinTree->currentParticle()->currentState().mass();
              if (mass_D0pismu > 7.) continue; // Last cut! from now on always save
              updateCounter(12, countersFlag);
              if (verbose) {cout << "B->D* mu found\n";}

              (*outputVecNtuplizer)["chi2_D0pismu"].push_back(res.chi2);
              (*outputVecNtuplizer)["dof_D0pismu"].push_back(res.dof);
              (*outputVecNtuplizer)["pval_D0pismu"].push_back(res.pval);
              (*outputVecNtuplizer)["mass_D0pismu"].push_back(mass_D0pismu);

              auto D0pismu = BKinTree->currentParticle();
              auto vtxB = BKinTree->currentDecayVertex();

              BKinTree->movePointerToTheFirstChild();
              auto refit_D0 = BKinTree->currentParticle();
              auto p4_refiD0pismu_D0 = vtxu::getTLVfromKinPart(refit_D0);
              BKinTree->movePointerToTheNextChild();
              auto refit_pis = BKinTree->currentParticle();
              auto p4_refitD0pismu_pis = vtxu::getTLVfromKinPart(refit_pis);
              AddTLVToOut(p4_refitD0pismu_pis, string("pis_refitD0pismu"), &(*outputVecNtuplizer));

              auto dm = (p4_refiD0pismu_D0 + p4_refitD0pismu_pis).M() - mass_piK;

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
              (*outputVecNtuplizer)["dm_D0pis_piK"].push_back(dm);
              (*outputVecNtuplizer)["cos_Dst_PV"].push_back(cos_Dst_PV);
            }

            if(n_B > 100) break;
          }
          if(n_B > 100) break;
        }
        if(n_B > 100) break;
      }
      if(n_B > 100) break;
    }

    (*outputNtuplizer)["Run"] = iEvent.run();
    (*outputNtuplizer)["LumiBlock"] = iEvent.luminosityBlock();
    (*outputNtuplizer)["eventNumber"] = iEvent.id().event();

    (*outputNtuplizer)["n_mu"] = n_mu;
    (*outputNtuplizer)["n_K"] = n_K;
    (*outputNtuplizer)["n_pi"] = n_pi;
    (*outputNtuplizer)["n_D0"] = n_D0;
    (*outputNtuplizer)["n_pis"] = n_pip;
    (*outputNtuplizer)["n_Dst"] = n_Dst;
    (*outputNtuplizer)["n_B"] = n_B;


    iEvent.put(move(outputNtuplizer), "outputNtuplizer");
    iEvent.put(move(outputVecNtuplizer), "outputVecNtuplizer");
    return;
}


void TagAndProbeBp2DststMuProducer::AddTLVToOut(TLorentzVector v, string n, map<string, vector<float>>* outv) {
  (*outv)[n+"_pt"].push_back(v.Pt());
  (*outv)[n+"_pz"].push_back(v.Pz());
  (*outv)[n+"_eta"].push_back(v.Eta());
  (*outv)[n+"_phi"].push_back(v.Phi());
  return;
}

bool TagAndProbeBp2DststMuProducer::qualityMuonID(pat::Muon m, reco::Vertex pVtx) {
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

void TagAndProbeBp2DststMuProducer::updateCounter(int idx, vector<bool> &cF) {
  if(!cF[idx]) {
    cF[idx] = true;
    counters[idx]++;
  }
  return;
}

double TagAndProbeBp2DststMuProducer::Mass_varM(TLorentzVector p1, double m1, TLorentzVector p2, double m2){
    double E1 = hypot(m1, p1.P());
    double E2 = hypot(m2, p2.P());
    double p1p2 = p1.Pt()*p2.Pt();
    p1p2 *= cos(p1.Phi() - p2.Phi()) + sinh(p1.Eta())*sinh(p2.Eta());

    double M = m1*m1 + m2*m2 + 2*(E1*E2 - p1p2);
    return sqrt(M);
}

DEFINE_FWK_MODULE(TagAndProbeBp2DststMuProducer);
