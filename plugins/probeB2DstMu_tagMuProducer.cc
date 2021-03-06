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

#include "VtxUtils.hh"

#define __pThad_min__ 0.5 // loose cut
#define __dzMax__ 1.0
#define __dRMax__ 2.0
#define __sigIPpfCand_min__ 2. // loose cut
#define __PvalChi2Vtx_min__ 0.05 // loose cut
#define __dmD0_max__ 0.1 // loose cut
#define __sigdxy_vtx_PV_min__ 2.0 // loose cut
#define __dmDst_max__ 0.050 // loose cut
#define __mass_D0pismu_max__ 10.0 // Some reasonable cut on the mass
#define __pTaddTracks_min__ 0.3 // loose cut
#define __mass_D0pismupi_max__ 10.0 // Some reasonable cut on the mass


using namespace std;

class probeB2DstMu_tagMuProducer : public edm::EDProducer {

public:

    explicit probeB2DstMu_tagMuProducer(const edm::ParameterSet &iConfig);
    void AddTLVToOut(TLorentzVector, string, map<string, vector<float>>*);
    bool qualityMuonID(pat::Muon, reco::Vertex);
    void updateCounter(int, vector<bool>&);

    ~probeB2DstMu_tagMuProducer() {
      cout << Form("probeB2DstMu_tagMuProducer counters:%d\n", (int)counters.size());
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
    edm::EDGetTokenT<vector<pat::Muon>> muonSrc_;
    edm::EDGetTokenT<vector<pat::Muon>> TrgMuonSrc_;

    double mass_Mu  = 0.10565;
    double mass_pi  = 0.13957;
    double mass_K   = 0.49367;
    double mass_D0  = 1.86483;
    double mass_Dst = 2.01026;
    double mass_B0  = 5.27963;

    int verbose = 0;

    vector<uint> counters;
};



probeB2DstMu_tagMuProducer::probeB2DstMu_tagMuProducer(const edm::ParameterSet &iConfig)
{
    verbose = iConfig.getParameter<int>( "verbose" );
    PFCandSrc_ = consumes<vector<pat::PackedCandidate>>(edm::InputTag("packedPFCandidates"));
    vtxSrc_ = consumes<vector<reco::Vertex>>(edm::InputTag("offlineSlimmedPrimaryVertices"));
    muonSrc_ = consumes<vector<pat::Muon>> ( edm::InputTag("slimmedMuons") );
    TrgMuonSrc_ = consumes<vector<pat::Muon>> ( iConfig.getParameter<edm::InputTag>( "trgMuons" ) );

    produces<map<string, float>>("outputNtuplizer");
    produces<map<string, vector<float>>>("outputVecNtuplizer");
    for(uint i=0; i<13; i++) counters.push_back(0);
}


void probeB2DstMu_tagMuProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    edm::Handle<vector<pat::PackedCandidate>> pfCandHandle;
    iEvent.getByToken(PFCandSrc_, pfCandHandle);
    unsigned int N_pfCand = pfCandHandle->size();

    edm::Handle<vector<reco::Vertex>> vtxHandle;
    iEvent.getByToken(vtxSrc_, vtxHandle);
    auto primaryVtx = (*vtxHandle)[0];

    edm::Handle<vector<pat::Muon>> muonsHandle;
    iEvent.getByToken(muonSrc_, muonsHandle);

    edm::Handle<vector<pat::Muon>> trgMuonsHandle;
    iEvent.getByToken(TrgMuonSrc_, trgMuonsHandle);
    auto nTrgMu = trgMuonsHandle->size();

    // Output collection
    unique_ptr<map<string, float>> outputNtuplizer(new map<string, float>);
    unique_ptr<map<string, vector<float>>> outputVecNtuplizer(new map<string, vector<float>>);

    int n_mu = 0, n_K = 0, n_pi = 0, n_D0 = 0, n_pis = 0, n_Dst = 0, n_B = 0;

    if (verbose) {cout <<"-------------------- Evt -----------------------\n";}
    vector<bool> countersFlag(counters.size(), false);
    counters[0]++;
    for(uint i_mu = 0; i_mu < muonsHandle->size(); i_mu++) {
      auto mu = (*muonsHandle)[i_mu];
      //######## Require probe side ###################
      if (mu.triggered("HLT_Mu*_IP*_part*_v*") && nTrgMu ==  1) continue;

      //######## Require muon quality ##################
      if (mu.charge() != 1) continue;
      if (mu.innerTrack().isNull()) continue;
      if (fabs(mu.eta()) > 2.2) continue;
      if (mu.pt() < 3.25) continue;
      updateCounter(1, countersFlag);

      vector<reco::Vertex> possibleVtxs = {};
      for(uint i_vtx = 0; i_vtx<vtxHandle->size(); i_vtx++) {
        auto vtx = (*vtxHandle)[i_vtx];
        if(vtx.normalizedChi2() > 1.5) continue;
        bool isSoft = mu.isSoftMuon(vtx);
        if(isSoft) possibleVtxs.push_back(vtx);
      }
      if (possibleVtxs.size() == 0) continue;
      updateCounter(2, countersFlag);

      n_mu++;
      /*
      ############################################################################
                                Look for the K+ (321)
      ############################################################################
      */
      for(uint i_K = 0; i_K < N_pfCand; ++i_K) {
        const pat::PackedCandidate & K = (*pfCandHandle)[i_K];
        if (!K.hasTrackDetails()) continue;
        //Require a positive charged hadron
        if (K.pdgId() != 211 ) continue;
        if (K.pt() < __pThad_min__) continue;
        // Require to be close to the trigger muon;
        auto K_tk = K.bestTrack();
        if (fabs(K_tk->dz(primaryVtx.position()) - mu.muonBestTrack()->dz(primaryVtx.position())) > __dzMax__) continue;
        if (vtxu::dR(K.phi(), mu.phi(), K.eta(), mu.eta()) > __dRMax__) continue;
        // Require significant impact parameter
        auto dxy = K_tk->dxy(primaryVtx.position());
        auto sigdxy_K_PV = fabs(dxy)/K_tk->dxyError();
        auto K_norm_chi2 = K_tk->normalizedChi2();
        auto K_N_valid_hits = K_tk->numberOfValidHits();
        if (sigdxy_K_PV < __sigIPpfCand_min__) continue;
        updateCounter(3, countersFlag);

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
          if (fabs(pi_tk->dz(primaryVtx.position()) - mu.muonBestTrack()->dz(primaryVtx.position())) > __dzMax__) continue;
          if (vtxu::dR(pi.phi(), mu.phi(), pi.eta(), mu.eta()) > __dRMax__) continue;
          // Require significant impact parameter
          auto dxy = pi_tk->dxy(primaryVtx.position());
          auto sigdxy_pi_PV = fabs(dxy)/pi.dxyError();
          auto pi_norm_chi2 = pi_tk->normalizedChi2();
          auto pi_N_valid_hits = pi_tk->numberOfValidHits();
          if (sigdxy_pi_PV < __sigIPpfCand_min__) continue;
          updateCounter(4, countersFlag);

          n_pi++;

          //Fit the vertex
          auto D0KinTree = vtxu::FitD0(iSetup, pi, K, false);
          auto res_piK = vtxu::fitQuality(D0KinTree, __PvalChi2Vtx_min__);
          if(!res_piK.isGood) continue;
          updateCounter(5, countersFlag);

          D0KinTree->movePointerToTheTop();
          auto mass_piK = D0KinTree->currentParticle()->currentState().mass();
          if (fabs(mass_piK - mass_D0) > __dmD0_max__) continue;
          updateCounter(6, countersFlag);

          auto D0 = D0KinTree->currentParticle();
          auto vtxD0 = D0KinTree->currentDecayVertex();

          auto cos_D0_PV = vtxu::computePointingCos(primaryVtx, vtxD0, D0);
          auto cosT_D0_PV = vtxu::computePointingCosTransverse(primaryVtx, vtxD0, D0);
          auto d_vtxD0_PV = vtxu::vtxsDistance(primaryVtx, vtxD0);
          auto sigd_vtxD0_PV = d_vtxD0_PV.first/d_vtxD0_PV.second;
          auto dxy_vtxD0_PV = vtxu::vtxsTransverseDistance(primaryVtx, vtxD0);
          auto sigdxy_vtxD0_PV = dxy_vtxD0_PV.first/dxy_vtxD0_PV.second;

          if(sigdxy_vtxD0_PV < __sigdxy_vtx_PV_min__) continue;
          updateCounter(7, countersFlag);

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
            if (fabs(pis_tk->dz(primaryVtx.position()) - mu.muonBestTrack()->dz(primaryVtx.position())) > __dzMax__) continue;
            if (vtxu::dR(pis.phi(), mu.phi(), pis.eta(), mu.eta()) > __dRMax__) continue;
            // Require significant impact parameter
            auto dxy = pis_tk->dxy(primaryVtx.position());
            auto sigdxy_pis_PV = fabs(dxy)/pis.dxyError();
            auto pis_norm_chi2 = pis_tk->normalizedChi2();
            auto pis_N_valid_hits = pis_tk->numberOfValidHits();
            if (sigdxy_pis_PV < __sigIPpfCand_min__) continue;
            updateCounter(8, countersFlag);

            n_pis++;

            // Fit the Dst vertex
            auto DstKinTree = vtxu::FitDst(iSetup, pis, D0, false);
            auto res_D0pis = vtxu::fitQuality(DstKinTree, __PvalChi2Vtx_min__);
            if(!res_D0pis.isGood) continue;
            updateCounter(9, countersFlag);

            DstKinTree->movePointerToTheTop();
            auto mass_D0pis = DstKinTree->currentParticle()->currentState().mass();
            if (fabs(mass_D0pis - mass_Dst) > __dmDst_max__) continue;
            updateCounter(10, countersFlag);

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
            auto BKinTree = vtxu::FitB_D0pismu(iSetup, D0, pis, mu);
            auto res = vtxu::fitQuality(BKinTree, __PvalChi2Vtx_min__);
            if(!res.isGood) continue;
            updateCounter(11, countersFlag);

            BKinTree->movePointerToTheTop();
            auto mass_D0pismu = BKinTree->currentParticle()->currentState().mass();
            if (mass_D0pismu > __mass_D0pismu_max__) continue; // Last cut! from now on always save
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

            // Looking for the best vertex to associate it with
            uint i_best = 0;
            auto maxCos = vtxu::computePointingCos(possibleVtxs[0], vtxB, D0pismu);
            for(uint i_vtx = 1; i_vtx < possibleVtxs.size(); i_vtx++) {
              auto auxCos = vtxu::computePointingCos(possibleVtxs[i_vtx], vtxB, D0pismu);
              if(auxCos > maxCos) {
                maxCos = auxCos;
                i_best = i_vtx;
              }
            }
            auto bestVtx = possibleVtxs[i_best];

            auto cos_D0pismu_PV = vtxu::computePointingCos(bestVtx, vtxB, D0pismu);
            auto cosT_D0pismu_PV = vtxu::computePointingCosTransverse(bestVtx, vtxB, D0pismu);
            auto d_vtxD0pismu_PV = vtxu::vtxsDistance(bestVtx, vtxB);
            auto sigd_vtxD0pismu_PV = d_vtxD0pismu_PV.first/d_vtxD0pismu_PV.second;
            auto dxy_vtxD0pismu_PV = vtxu::vtxsDistance(bestVtx, vtxB);
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

            TVector3 flightB(vtxB->position().x() - bestVtx.position().x(),
                             vtxB->position().y() - bestVtx.position().y(),
                             vtxB->position().z() - bestVtx.position().z()
                            );

            // // ------------- Transverse Approx ------------------
            double pt_B_reco = p4_vis.Pt() * mass_B0/ p4_vis.M();

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
                    Veto the presence of additional tracks in the D* mu vertex
            ############################################################################
            */
            uint N_compatible_tk = 0;
            vector<double> tksAdd_massHad = {};
            vector<double> tksAdd_massVis = {};
            vector<double> tksAdd_pval = {};
            vector<double> tksAdd_pt = {};
            vector<double> tksAdd_sigdca_vtxB = {};
            vector<double> tksAdd_cos_PV = {};
            if(verbose) {cout << "Looking for additional tracks" << endl;}
            for(uint i_tk = 0; i_tk < N_pfCand; ++i_tk) {
              // Pion different from the pion from D0
              if( i_tk==i_K || i_tk==i_pi || i_tk==i_pis) continue;

              const pat::PackedCandidate & ptk = (*pfCandHandle)[i_tk];
              if (!ptk.hasTrackDetails()) continue;
              //Require a positive charged hadron
              if (abs(ptk.pdgId()) != 211 ) continue;
              if (ptk.pt() < __pTaddTracks_min__) continue;
              // Require to be close to the trigger muon;
              auto tk = ptk.bestTrack();
              if (fabs(tk->dz(primaryVtx.position()) - mu.muonBestTrack()->dz(primaryVtx.position())) > __dzMax__) continue;
              if (vtxu::dR(ptk.phi(), mu.phi(), ptk.eta(), mu.eta()) > __dRMax__) continue;

              GlobalPoint auxp(vtxB->position().x(), vtxB->position().y(), vtxB->position().z());
              auto dca = vtxu::computeDCA(iSetup, ptk, auxp);

              // Check if it makes a consistent vertex
              auto kinTree = vtxu::Fit_D0pismupi(iSetup, D0, pis, mu, ptk);
              auto res = vtxu::fitQuality(kinTree, __PvalChi2Vtx_min__);
              if(!res.isGood) continue;
              if (verbose) {cout << Form("Trk pars: pt=%.2f eta=%.2f phi=%.2f\n", ptk.pt(), ptk.eta(), ptk.phi());}

              kinTree->movePointerToTheTop();
              auto p_vis = kinTree->currentParticle();
              auto m_vis = p_vis->currentState().mass();
              if (m_vis > __mass_D0pismupi_max__) continue;

              kinTree->movePointerToTheFirstChild();
              auto refit_D0 = kinTree->currentParticle();
              kinTree->movePointerToTheNextChild();
              auto refit_pis = kinTree->currentParticle();
              kinTree->movePointerToTheNextChild();
              auto refit_mu = kinTree->currentParticle();
              kinTree->movePointerToTheNextChild();
              auto refit_pi = kinTree->currentParticle();

              auto p4_D0pispi = vtxu::getTLVfromKinPart(refit_pis) + vtxu::getTLVfromKinPart(refit_D0) + vtxu::getTLVfromKinPart(refit_pi);
              auto m_D0pispi = p4_D0pispi.M();



              tksAdd_massVis.push_back(m_vis);
              tksAdd_massHad.push_back(m_D0pispi);
              tksAdd_pval.push_back(res.pval);
              tksAdd_pt.push_back(ptk.pt());
              tksAdd_sigdca_vtxB.push_back(fabs(dca.first)/dca.second);
              tksAdd_cos_PV.push_back(vtxu::computePointingCos(bestVtx, vtxB, refit_pi));
              N_compatible_tk++;
            }
            if ( (*outputVecNtuplizer)["nTksAdd"].size() == 0 ) {
              (*outputVecNtuplizer)["tksAdd_massVis"] = {};
              (*outputVecNtuplizer)["tksAdd_massHad"] = {};
              (*outputVecNtuplizer)["tksAdd_pval"] = {};
              (*outputVecNtuplizer)["tksAdd_pt"] = {};
              (*outputVecNtuplizer)["tksAdd_sigdca_vtxB"] = {};
              (*outputVecNtuplizer)["tksAdd_cos_PV"] = {};
            }
            (*outputVecNtuplizer)["nTksAdd"].push_back(N_compatible_tk);
            if(verbose) {cout << "Compatible tracks: " << N_compatible_tk << endl;}
            for(uint i = 0; i < N_compatible_tk; i++){
              (*outputVecNtuplizer)["tksAdd_massVis"].push_back(tksAdd_massVis[i]);
              (*outputVecNtuplizer)["tksAdd_massHad"].push_back(tksAdd_massHad[i]);
              (*outputVecNtuplizer)["tksAdd_pval"].push_back(tksAdd_pval[i]);
              (*outputVecNtuplizer)["tksAdd_pt"].push_back(tksAdd_pt[i]);
              (*outputVecNtuplizer)["tksAdd_sigdca_vtxB"].push_back(tksAdd_sigdca_vtxB[i]);
              (*outputVecNtuplizer)["tksAdd_cos_PV"].push_back(tksAdd_cos_PV[i]);
            }

            int auxC = 0;
            for(auto auxN : (*outputVecNtuplizer)["nTksAdd"]) auxC += auxN;
            if( auxC != int ((*outputVecNtuplizer)["tksAdd_massVis"].size()) ){
              cout << "Number of tracks and tracks details lenght not matching" << endl;
              assert(false);
            }


            n_B++;


            AddTLVToOut(vtxu::getTLVfromMuon(mu, mass_Mu), string("mu"), &(*outputVecNtuplizer));
            GlobalPoint auxp(vtxDst->position().x(), vtxDst->position().y(), vtxDst->position().z());
            auto dca = vtxu::computeDCA(iSetup, mu, auxp);
            (*outputVecNtuplizer)["mu_dca_vtxDst"].push_back(dca.first);
            auto tkMu = mu.innerTrack();
            auto dxyMu = tkMu->dxy(bestVtx.position());
            (*outputVecNtuplizer)["mu_sigdxy_PV"].push_back(fabs(dxyMu)/tkMu->dxyError());

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

            (*outputVecNtuplizer)["sigdxy_pis_PV"].push_back(sigdxy_pis_PV);
            (*outputVecNtuplizer)["pis_norm_chi2"].push_back(pis_norm_chi2);
            (*outputVecNtuplizer)["pis_N_valid_hits"].push_back(pis_N_valid_hits);
            AddTLVToOut(vtxu::getTLVfromCand(pis, mass_pi), string("pis"), &(*outputVecNtuplizer));

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
    (*outputNtuplizer)["n_pis"] = n_pis;
    (*outputNtuplizer)["n_Dst"] = n_Dst;
    (*outputNtuplizer)["n_B"] = n_B;


    iEvent.put(move(outputNtuplizer), "outputNtuplizer");
    iEvent.put(move(outputVecNtuplizer), "outputVecNtuplizer");
    return;
}


void probeB2DstMu_tagMuProducer::AddTLVToOut(TLorentzVector v, string n, map<string, vector<float>>* outv) {
  (*outv)[n+"_pt"].push_back(v.Pt());
  (*outv)[n+"_pz"].push_back(v.Pz());
  (*outv)[n+"_eta"].push_back(v.Eta());
  (*outv)[n+"_phi"].push_back(v.Phi());
  return;
}

bool probeB2DstMu_tagMuProducer::qualityMuonID(pat::Muon m, reco::Vertex pVtx) {
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

void probeB2DstMu_tagMuProducer::updateCounter(int idx, vector<bool> &cF) {
  if(!cF[idx]) {
    cF[idx] = true;
    counters[idx]++;
  }
  return;
}

DEFINE_FWK_MODULE(probeB2DstMu_tagMuProducer);
