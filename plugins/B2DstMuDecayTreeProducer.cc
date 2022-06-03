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

#define __pThad_min__ 0.6 // loose cut
#define __dzMax__ 1.0
#define __dRMax__ 2.0
#define __sigIPpfCand_min__ 2. // loose cut
#define __PvalChi2Vtx_min__ 0.05 // loose cut
#define __dmD0_max__ 0.05 // loose cut
#define __sigdxy_vtx_PV_min__ 2.0 // loose cut
#define __dmDst_max__ 0.15 // loose cut
#define __dm_DstMiunsD0_max__ 0.003
#define __mass_D0pismu_max__ 8. // Some reasonable cut on the mass
#define __pTaddTracks_min__ 0.3 // loose cut
#define __mass_D0pismupi_max__ 10. // Some reasonable cut on the mass


using namespace std;

class B2DstMuDecayTreeProducer : public edm::EDProducer {

public:

    explicit B2DstMuDecayTreeProducer(const edm::ParameterSet &iConfig);
    void AddTLVToOut(TLorentzVector, string, map<string, vector<float>>*);
    void updateCounter(int, vector<bool>&);
    double Mass_varM(TLorentzVector, double, TLorentzVector, double);

    ~B2DstMuDecayTreeProducer() {
      cout << Form("Total number of fit crashed %u", fitCrash) << endl;
      cout << Form("B2DstMuDecayTreeProducer counters:%d\n", (int)counters.size());
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
    edm::EDGetTokenT<reco::BeamSpot> beamSpotSrc_;

    int charge_muon = +1;
    int charge_K = +1;
    int charge_pi = -1;
    int charge_pis = -1;

    double mass_Mu  = 0.105658;
    double mass_pi  = 0.139570;
    double mass_K   = 0.493677;
    double mass_D0  = 1.86483;
    double mass_Dst = 2.01026;
    double mass_B0  = 5.27963;

    int verbose = 0;
    int isMC = 0;

    vector<uint> counters;
    uint fitCrash = 0;
};



B2DstMuDecayTreeProducer::B2DstMuDecayTreeProducer(const edm::ParameterSet &iConfig) :
    beamSpotSrc_( consumes<reco::BeamSpot> ( edm::InputTag("offlineBeamSpot") ) )
{
    PFCandSrc_ = consumes<vector<pat::PackedCandidate>>(edm::InputTag("packedPFCandidates"));
    vtxSrc_ = consumes<vector<reco::Vertex>>(edm::InputTag("offlineSlimmedPrimaryVertices"));
    TrgMuonSrc_ = consumes<vector<pat::Muon>> ( iConfig.getParameter<edm::InputTag>( "trgMuons" ) );

    // If charge_muon == 0 accept both cahrges
    charge_muon = iConfig.getParameter<int>( "charge_muon" );

    // Charges relative to muon charge
    charge_K = iConfig.getParameter<int>( "charge_K" );
    charge_pi = iConfig.getParameter<int>( "charge_pi" );
    charge_pis = iConfig.getParameter<int>( "charge_pis" );

    verbose = iConfig.getParameter<int>( "verbose" );

    isMC = iConfig.getParameter<int>("isMC");
    vtxu::set_isMC(isMC);

    produces<map<string, float>>("outputNtuplizer");
    produces<map<string, vector<float>>>("outputVecNtuplizer");
    for(uint i=0; i<13; i++) counters.push_back(0);
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

    edm::Handle<reco::BeamSpot> beamSpotHandle;
    iEvent.getByToken(beamSpotSrc_, beamSpotHandle);

    // Output collection
    unique_ptr<map<string, float>> outputNtuplizer(new map<string, float>);
    unique_ptr<map<string, vector<float>>> outputVecNtuplizer(new map<string, vector<float>>);

    int n_mu = 0, n_K = 0, n_pi = 0, n_D0 = 0, n_pis = 0, n_Dst = 0, n_B = 0;

    vector<reco::Vertex> possibleVtxs = {};
    for(uint i_vtx = 0; i_vtx<vtxHandle->size(); i_vtx++) {
      auto vtx = (*vtxHandle)[i_vtx];
      /* Cut on the number of degrees of freedom. The number of degrees of
       * freedom is calculated as:
       *
       *     ndof = -3 + 2*(w_1 + w_2 + ...)
       *
       *  where w_1, w_2, etc. are weights for each track associated with the
       *  vertex. Quoting from https://arxiv.org/pdf/1405.6569.pdf:
       *
       *  > In the adaptive vertex fit, each track in the vertex is assigned
       *  > a weight between 0 and 1, which reflects the likelihood that it
       *  > genuinely > belongs to the vertex.
       *
       *  Therefore, the number of degrees of freedom is strongly correlated
       *  with the number of tracks. We cut here requiring more than 4 which
       *  translates roughly to at least 4 "good" tracks. This is because for
       *  ndof <= 4 the data and MC do not seem to agree, and there is an
       *  excess of vertices with low ndof in data. */
      if (vtx.ndof() <= 4) continue;
      reco::Vertex tmp = vtxu::refit_vertex(iEvent, iSetup, i_vtx, 1, *pfCandHandle);
      if (tmp.isValid())
        possibleVtxs.push_back(tmp);
      else {
        cout << "[ERROR] Invalid vertex refit for " << i_vtx << endl;
        // possibleVtxs.push_back(vtx);
      }
    }

    if (verbose) {cout <<"-------------------- Evt -----------------------\n";}
    vector<bool> countersFlag(counters.size(), false);
    counters[0]++;
    for(uint i_trgMu = 0; i_trgMu < trgMuonsHandle->size(); i_trgMu++) {
      //######## Require muon quality ##################
      if (possibleVtxs.size() == 0) continue;
      auto trgMu = (*trgMuonsHandle)[i_trgMu];
      if (charge_muon != 0 && trgMu.charge() != charge_muon) continue;
      if (trgMu.innerTrack().isNull()) continue;
      if (fabs(trgMu.eta()) > 1.5) continue;
      updateCounter(1, countersFlag);

      if (!trgMu.isMediumMuon()) continue;
      updateCounter(2, countersFlag);

      n_mu++;

      int *valid = (int *) malloc(sizeof(int)*N_pfCand);

      /* Create an array which tells us whether each track is potentially
       * valid. Here we check dz and dR. This way we don't have to recompute
       * this quantity for each for loop. */
      for (uint i = 0; i < N_pfCand; i++) {
        /* Set it to not valid by default. */
        valid[i] = 0;
        const pat::PackedCandidate &ptk = (*pfCandHandle)[i];
        if (!ptk.hasTrackDetails()) continue;
        if (ptk.pt() < __pThad_min__) continue;
        // Require to be close to the trigger muon;
        auto tk = ptk.bestTrack();
        if (fabs(tk->dz(primaryVtx.position()) - trgMu.muonBestTrack()->dz(primaryVtx.position())) > __dzMax__) continue;
        if (vtxu::dR(ptk.phi(), trgMu.phi(), ptk.eta(), trgMu.eta()) > __dRMax__) continue;
        valid[i] = 1;
      }
      /*
      ############################################################################
                                Look for the K+ (321)
      ############################################################################
      */
      for(uint i_K = 0; i_K < N_pfCand; ++i_K) {
        if (!valid[i_K]) continue;
        const pat::PackedCandidate & K = (*pfCandHandle)[i_K];
        if (!K.hasTrackDetails()) continue;
        if (K.pdgId() != charge_K*trgMu.charge()*211 ) continue;
        if (K.pt() < __pThad_min__) continue;
        // Require to be close to the trigger muon;
        auto K_tk = vtxu::fix_track(K.bestTrack());
        if (fabs(K_tk.dz(primaryVtx.position()) - trgMu.muonBestTrack()->dz(primaryVtx.position())) > __dzMax__) continue;
        if (vtxu::dR(K.phi(), trgMu.phi(), K.eta(), trgMu.eta()) > __dRMax__) continue;
        // Require significant impact parameter
        auto dxy = K_tk.dxy(*beamSpotHandle);
        double dxyErr = vtxu::dxyError(K_tk,*beamSpotHandle);
        auto sigdxy_K_PV = fabs(dxy)/dxyErr;
        auto K_norm_chi2 = K_tk.normalizedChi2();
        auto K_N_valid_hits = K_tk.numberOfValidHits();
        if (sigdxy_K_PV < __sigIPpfCand_min__) continue;
        updateCounter(3, countersFlag);
        if (verbose) {cout << Form("K%s found at %u\n", K.pdgId() > 0?"+":"-", i_K);}
        n_K++;

        /*
        ############################################################################
                                 Look for the pi- (-211)
        ############################################################################
        */
        for(uint i_pi = 0; i_pi < N_pfCand; ++i_pi) {
          if (!valid[i_pi]) continue;
          if(i_pi==i_K) continue;

          const pat::PackedCandidate & pi = (*pfCandHandle)[i_pi];
          if (!pi.hasTrackDetails()) continue;
          if (pi.pdgId() != charge_pi*trgMu.charge()*211 ) continue;
          //Require a minimum pt
          if(pi.pt() < __pThad_min__) continue;
          // Require to be close to the trigger muon;
          auto pi_tk = vtxu::fix_track(pi.bestTrack());
          if (fabs(pi_tk.dz(primaryVtx.position()) - trgMu.muonBestTrack()->dz(primaryVtx.position())) > __dzMax__) continue;
          if (vtxu::dR(pi.phi(), trgMu.phi(), pi.eta(), trgMu.eta()) > __dRMax__) continue;
          // Require significant impact parameter
          auto dxy = pi_tk.dxy(*beamSpotHandle);
          double dxyErr = vtxu::dxyError(pi_tk,*beamSpotHandle);
          auto sigdxy_pi_PV = fabs(dxy)/dxyErr;
          auto pi_norm_chi2 = pi_tk.normalizedChi2();
          auto pi_N_valid_hits = pi_tk.numberOfValidHits();
          if (sigdxy_pi_PV < __sigIPpfCand_min__) continue;
          updateCounter(4, countersFlag);
          if (verbose) {cout << Form("pi%s found at %u\n", pi.pdgId()>0?"+":"-", i_pi);}
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
          auto dxy_vtxD0_PV = vtxu::vtxsTransverseDistanceFromBeamSpot(*beamSpotHandle, vtxD0);
          auto sigdxy_vtxD0_PV = dxy_vtxD0_PV.first/dxy_vtxD0_PV.second;

          if(sigdxy_vtxD0_PV < __sigdxy_vtx_PV_min__) continue;
          if (verbose) {cout << "pi-K vertex displacement passed\n";}
          updateCounter(7, countersFlag);
          if (verbose) {cout << "D0 candidate\n";}
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
            if (pis.pdgId() != charge_pis*trgMu.charge()*211 ) continue;
            // Require to be close to the trigger muon;
            auto pis_tk = vtxu::fix_track(pis.bestTrack());
            if (fabs(pis_tk.dz(primaryVtx.position()) - trgMu.muonBestTrack()->dz(primaryVtx.position())) > __dzMax__) continue;
            if (vtxu::dR(pis.phi(), trgMu.phi(), pis.eta(), trgMu.eta()) > __dRMax__) continue;
            // Require significant impact parameter
            auto dxy = pis_tk.dxy(*beamSpotHandle);
            double dxyErr = vtxu::dxyError(pis_tk,*beamSpotHandle);
            auto sigdxy_pis_PV = fabs(dxy)/dxyErr;
            auto pis_norm_chi2 = pis_tk.normalizedChi2();
            auto pis_N_valid_hits = pis_tk.numberOfValidHits();
            if (sigdxy_pis_PV < __sigIPpfCand_min__) continue;
            updateCounter(8, countersFlag);
            if (verbose) {cout << Form("pis%s found at %u\n", pis.pdgId()>0?"+":"-", i_pis);}
            n_pis++;

            // Fit the Dst vertex
            RefCountedKinematicTree DstKinTree;
            // try {
            DstKinTree = vtxu::FitDst(iSetup, pis, D0, false);
            // }
            // catch(...) {
            //   if (verbose) {cout << "Fit crushed...\n";}
            //   fitCrash++;
            //   continue;
            // }
            auto res_D0pis = vtxu::fitQuality(DstKinTree, __PvalChi2Vtx_min__);
            if(!res_D0pis.isGood) continue;
            if (verbose) {cout << "D0-pis vertex fit good\n";}
            updateCounter(9, countersFlag);

            DstKinTree->movePointerToTheTop();
            auto mass_D0pis = DstKinTree->currentParticle()->currentState().mass();
            if (fabs(mass_D0pis - mass_Dst) > __dmDst_max__) continue;
            if (fabs(mass_D0pis - mass_piK - (mass_Dst - mass_D0) ) > __dm_DstMiunsD0_max__) continue;
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
            auto BKinTree = vtxu::FitB_D0pismu(iSetup, D0, pis, trgMu);
            auto res = vtxu::fitQuality(BKinTree, __PvalChi2Vtx_min__);
            // cout << "[DEBUG]: pval(D0pimu) = " << res.pval << endl;
            if(!res.isValid || res.chi2 > 50 || res.pval < 0.001) continue;
            updateCounter(11, countersFlag);

            BKinTree->movePointerToTheTop();
            auto mass_D0pismu = BKinTree->currentParticle()->currentState().mass();
            // cout << "[DEBUG]: m(D0pimu) = " << mass_D0pismu << endl;
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

            (*outputVecNtuplizer)["N_goodVtx"].push_back(possibleVtxs.size());

            (*outputVecNtuplizer)["PV_chi2"].push_back(bestVtx.chi2());
            (*outputVecNtuplizer)["PV_ndof"].push_back(bestVtx.ndof());
            (*outputVecNtuplizer)["PV_nTracks"].push_back(bestVtx.nTracks());

            (*outputVecNtuplizer)["vtx_PV_x"].push_back(bestVtx.position().x());
            (*outputVecNtuplizer)["vtx_PV_y"].push_back(bestVtx.position().y());
            (*outputVecNtuplizer)["vtx_PV_z"].push_back(bestVtx.position().z());

            (*outputVecNtuplizer)["vtx_B_decay_x"].push_back(vtxB->position().x());
            (*outputVecNtuplizer)["vtx_B_decay_y"].push_back(vtxB->position().y());
            (*outputVecNtuplizer)["vtx_B_decay_z"].push_back(vtxB->position().z());

            auto cos_D0pismu_PV = vtxu::computePointingCos(bestVtx, vtxB, D0pismu);
            auto cosT_D0pismu_PV = vtxu::computePointingCosTransverse(bestVtx, vtxB, D0pismu);
            auto d_vtxD0pismu_PV = vtxu::vtxsDistance(bestVtx, vtxB);
            auto sigd_vtxD0pismu_PV = d_vtxD0pismu_PV.first/d_vtxD0pismu_PV.second;
            auto dxy_vtxD0pismu_PV = vtxu::vtxsTransverseDistance(bestVtx, vtxB);
            auto sigdxy_vtxD0pismu_PV = dxy_vtxD0pismu_PV.first/dxy_vtxD0pismu_PV.second;
            (*outputVecNtuplizer)["cos_D0pismu_PV"].push_back(cos_D0pismu_PV);
            (*outputVecNtuplizer)["cosT_D0pismu_PV"].push_back(cosT_D0pismu_PV);
            (*outputVecNtuplizer)["d_vtxD0pismu_PV"].push_back(d_vtxD0pismu_PV.first);
            (*outputVecNtuplizer)["sigd_vtxD0pismu_PV"].push_back(sigd_vtxD0pismu_PV);
            (*outputVecNtuplizer)["dxy_vtxD0pismu_PV"].push_back(dxy_vtxD0pismu_PV.first);
            (*outputVecNtuplizer)["sigdxy_vtxD0pismu_PV"].push_back(sigdxy_vtxD0pismu_PV);

            float localVertexDensity = 0;
            float localVertexDensity_30mm = 0;
            float localVertexDensity_10mm = 0;
            float localVertexDensity_5mm = 0;
            float localVertexDensity_2mm = 0;
            float localVertexDensity_1mm = 0;
            float localVertexDensity_cos800 = 0;
            float localVertexDensity_cos950 = 0;
            float localVertexDensity_cos990 = 0;
            float localVertexDensity_cos995 = 0;
            float localVertexDensity_cos999 = 0;
            for (uint i_vtx = 0; i_vtx < possibleVtxs.size(); i_vtx++) {
              auto vtx = possibleVtxs[i_vtx];
              float dz = bestVtx.position().z() - vtx.position().z();
              if(fabs(dz) < 1.5*__dzMax__) localVertexDensity++;
              if(fabs(dz) < 1.5) localVertexDensity_30mm++;
              if(fabs(dz) < 0.5) localVertexDensity_10mm++;
              if(fabs(dz) < 0.25) localVertexDensity_5mm++;
              if(fabs(dz) < 0.1) localVertexDensity_2mm++;
              if(fabs(dz) < 0.05) localVertexDensity_1mm++;
              auto cos = vtxu::computePointingCos(vtx, vtxB, D0pismu);
              if (cos > 0.800) localVertexDensity_cos800++;
              if (cos > 0.950) localVertexDensity_cos950++;
              if (cos > 0.990) localVertexDensity_cos990++;
              if (cos > 0.995) localVertexDensity_cos995++;
              if (cos > 0.999) localVertexDensity_cos999++;
            }
            (*outputVecNtuplizer)["localVertexDensity"].push_back(localVertexDensity/(2*1.5*__dzMax__));
            (*outputVecNtuplizer)["localVertexDensity_30mm"].push_back(localVertexDensity_30mm);
            (*outputVecNtuplizer)["localVertexDensity_10mm"].push_back(localVertexDensity_10mm);
            (*outputVecNtuplizer)["localVertexDensity_5mm"].push_back(localVertexDensity_5mm);
            (*outputVecNtuplizer)["localVertexDensity_2mm"].push_back(localVertexDensity_2mm);
            (*outputVecNtuplizer)["localVertexDensity_1mm"].push_back(localVertexDensity_1mm);
            (*outputVecNtuplizer)["localVertexDensity_cos800"].push_back(localVertexDensity_cos800);
            (*outputVecNtuplizer)["localVertexDensity_cos950"].push_back(localVertexDensity_cos950);
            (*outputVecNtuplizer)["localVertexDensity_cos990"].push_back(localVertexDensity_cos990);
            (*outputVecNtuplizer)["localVertexDensity_cos995"].push_back(localVertexDensity_cos995);
            (*outputVecNtuplizer)["localVertexDensity_cos999"].push_back(localVertexDensity_cos999);

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
            vector<double> tksAdd_massMuTk = {};
            vector<double> tksAdd_pval = {};
            vector<double> tksAdd_charge = {};
            vector<double> tksAdd_pdgId = {};
            vector<double> tksAdd_pt = {};
            vector<double> tksAdd_ptError = {};
            vector<double> tksAdd_eta = {};
            vector<double> tksAdd_etaError = {};
            vector<double> tksAdd_phi = {};
            vector<double> tksAdd_phiError = {};
            vector<double> tksAdd_dz = {};
            vector<double> tksAdd_lostInnerHits = {};
            vector<double> tksAdd_sigdca_vtxB = {};
            vector<double> tksAdd_cos_PV = {};
            if(verbose) {cout << "Looking for additional tracks compatible with D0pismu vertex" << endl;}

            double dR_max = vtxu::dR(K.phi(), trgMu.phi(), K.eta(), trgMu.eta());
            dR_max = max(dR_max, vtxu::dR(pi.phi(), trgMu.phi(), pi.eta(), trgMu.eta()));
            dR_max = max(dR_max, vtxu::dR(pis.phi(), trgMu.phi(), pis.eta(), trgMu.eta()));
            dR_max = min(1., dR_max*1.5);
            if(verbose) {cout << Form("Looking for additional neutrals within dR=%.3f from muon", dR_max) << endl;}
            uint N_compatible_neu = 0;
            vector<double> neu_pdgId = {};
            vector<double> neu_pt = {};
            vector<double> neu_eta = {};
            vector<double> neu_phi = {};
            vector<double> neu_energy = {};
            vector<double> neu_et2 = {};
            vector<double> neu_dR_fromMu = {};

            int i_muon_PFcand = -1;

            for(uint i_tk = 0; i_tk < N_pfCand; ++i_tk) {
              // PF candidate different from K, pi and pis
              if( i_tk==i_K || i_tk==i_pi || i_tk==i_pis) continue;

              const pat::PackedCandidate & ptk = (*pfCandHandle)[i_tk];

              double dR_fromMu = vtxu::dR(ptk.phi(), trgMu.phi(), ptk.eta(), trgMu.eta());
              // PF candidate not matching with the muon
              if ( fabs(ptk.pt() - trgMu.pt())/trgMu.pt() < 0.01  && dR_fromMu < 0.01 ) {
                if (ptk.pdgId() == trgMu.pdgId()) {
                    i_muon_PFcand = i_tk;
                    continue;
                }
              }

              if (ptk.charge() != 0 ) {
                if (!ptk.hasTrackDetails()) continue;
                if (ptk.pt() < __pTaddTracks_min__) continue;
                // Require to be close to the trigger muon;
                auto tk = vtxu::fix_track(ptk.bestTrack());
                if (fabs(tk.dz(primaryVtx.position()) - trgMu.muonBestTrack()->dz(primaryVtx.position())) > __dzMax__) continue;
                if (dR_fromMu > __dRMax__) continue;

                // GlobalPoint auxp(vtxB->position().x(), vtxB->position().y(), vtxB->position().z());
                // auto dca = vtxu::computeDCA(iSetup, ptk, auxp);
                auto dca = vtxu::computeIP3D(iSetup, ptk, p4_vis , vtxB);

                // Check if it makes a consistent vertex
                auto kinTree = vtxu::Fit_D0pismupi(iSetup, D0, pis, trgMu, ptk);
                auto res = vtxu::fitQuality(kinTree, __PvalChi2Vtx_min__);
                if(!res.isGood) continue;
                if (verbose) {cout << Form("Trk pars: pt=%.2f eta=%.2f phi=%.2f\n", ptk.pt(), ptk.eta(), ptk.phi());}

                kinTree->movePointerToTheTop();
                auto p_vis = kinTree->currentParticle();
                auto m_vis = p_vis->currentState().mass();
                if (m_vis > __mass_D0pismupi_max__) continue;

                kinTree->movePointerToTheFirstChild();
                auto refit_D0 = vtxu::getTLVfromKinPart(kinTree->currentParticle());
                kinTree->movePointerToTheNextChild();
                auto refit_pis = vtxu::getTLVfromKinPart(kinTree->currentParticle());
                kinTree->movePointerToTheNextChild();
                auto refit_mu = vtxu::getTLVfromKinPart(kinTree->currentParticle());
                kinTree->movePointerToTheNextChild();
                auto refit_pi = kinTree->currentParticle();
                auto refit_pi_p4 = vtxu::getTLVfromKinPart(refit_pi);

                auto m_D0pispi = (refit_pis + refit_D0 + refit_pi_p4).M();
                auto m_MuTk = (refit_pi_p4 + refit_mu).M();

                tksAdd_massVis.push_back(m_vis);
                tksAdd_massHad.push_back(m_D0pispi);
                tksAdd_massMuTk.push_back(m_MuTk);
                tksAdd_pval.push_back(res.pval);
                tksAdd_charge.push_back(ptk.charge());
                tksAdd_pdgId.push_back(ptk.pdgId());
                tksAdd_pt.push_back(refit_pi_p4.Pt());
                tksAdd_ptError.push_back(tk.ptError());
                tksAdd_eta.push_back(refit_pi_p4.Eta());
                tksAdd_etaError.push_back(tk.etaError());
                tksAdd_phi.push_back(refit_pi_p4.Phi());
                tksAdd_phiError.push_back(tk.phiError());
                tksAdd_dz.push_back(tk.dz(bestVtx.position()));
                tksAdd_lostInnerHits.push_back(ptk.lostInnerHits());
                tksAdd_sigdca_vtxB.push_back(dca.first/dca.second);
                tksAdd_cos_PV.push_back(vtxu::computePointingCos(bestVtx, vtxB, refit_pi));
                N_compatible_tk++;
              }
              else {
                if (dR_fromMu > dR_max) continue;
                if (ptk.energy() < 0.5) continue;

                neu_pdgId.push_back(ptk.pdgId());
                neu_pt.push_back(ptk.pt());
                neu_eta.push_back(ptk.eta());
                neu_phi.push_back(ptk.phi());
                neu_energy.push_back(ptk.energy());
                neu_et2.push_back(ptk.et2());
                neu_dR_fromMu.push_back(dR_fromMu);

                N_compatible_neu++;
              }
            }

            if ( (*outputVecNtuplizer)["nTksAdd"].size() == 0 ) {
              (*outputVecNtuplizer)["tksAdd_massVis"] = {};
              (*outputVecNtuplizer)["tksAdd_massHad"] = {};
              (*outputVecNtuplizer)["tksAdd_massMuTk"] = {};
              (*outputVecNtuplizer)["tksAdd_pval"] = {};
              (*outputVecNtuplizer)["tksAdd_charge"] = {};
              (*outputVecNtuplizer)["tksAdd_pdgId"] = {};
              (*outputVecNtuplizer)["tksAdd_pt"] = {};
              (*outputVecNtuplizer)["tksAdd_ptError"] = {};
              (*outputVecNtuplizer)["tksAdd_eta"] = {};
              (*outputVecNtuplizer)["tksAdd_etaError"] = {};
              (*outputVecNtuplizer)["tksAdd_phi"] = {};
              (*outputVecNtuplizer)["tksAdd_phiError"] = {};
              (*outputVecNtuplizer)["tksAdd_dz"] = {};
              (*outputVecNtuplizer)["tksAdd_lostInnerHits"] = {};
              (*outputVecNtuplizer)["tksAdd_sigIP3D_vtxB"] = {};
              (*outputVecNtuplizer)["tksAdd_cos_PV"] = {};
            }
            (*outputVecNtuplizer)["nTksAdd"].push_back(N_compatible_tk);
            if(verbose) {cout << "Compatible tracks: " << N_compatible_tk << endl;}
            for(uint i = 0; i < N_compatible_tk; i++){
              (*outputVecNtuplizer)["tksAdd_massVis"].push_back(tksAdd_massVis[i]);
              (*outputVecNtuplizer)["tksAdd_massHad"].push_back(tksAdd_massHad[i]);
              (*outputVecNtuplizer)["tksAdd_massMuTk"].push_back(tksAdd_massMuTk[i]);
              (*outputVecNtuplizer)["tksAdd_pval"].push_back(tksAdd_pval[i]);
              (*outputVecNtuplizer)["tksAdd_charge"].push_back(tksAdd_charge[i]);
              (*outputVecNtuplizer)["tksAdd_pdgId"].push_back(tksAdd_pdgId[i]);
              (*outputVecNtuplizer)["tksAdd_pt"].push_back(tksAdd_pt[i]);
              (*outputVecNtuplizer)["tksAdd_ptError"].push_back(tksAdd_ptError[i]);
              (*outputVecNtuplizer)["tksAdd_eta"].push_back(tksAdd_eta[i]);
              (*outputVecNtuplizer)["tksAdd_etaError"].push_back(tksAdd_etaError[i]);
              (*outputVecNtuplizer)["tksAdd_phi"].push_back(tksAdd_phi[i]);
              (*outputVecNtuplizer)["tksAdd_phiError"].push_back(tksAdd_phiError[i]);
              (*outputVecNtuplizer)["tksAdd_dz"].push_back(tksAdd_dz[i]);
              (*outputVecNtuplizer)["tksAdd_lostInnerHits"].push_back(tksAdd_lostInnerHits[i]);
              (*outputVecNtuplizer)["tksAdd_sigIP3D_vtxB"].push_back(tksAdd_sigdca_vtxB[i]);
              (*outputVecNtuplizer)["tksAdd_cos_PV"].push_back(tksAdd_cos_PV[i]);
            }

            int auxC = 0;
            for(auto auxN : (*outputVecNtuplizer)["nTksAdd"]) auxC += auxN;
            if( auxC != int ((*outputVecNtuplizer)["tksAdd_massVis"].size()) ){
              cout << "Number of tracks and tracks details lenght not matching" << endl;
              assert(false);
            }


            if ( (*outputVecNtuplizer)["nNeuAdd"].size() == 0 ) {
              (*outputVecNtuplizer)["neuAdd_pdgId"] = {};
              (*outputVecNtuplizer)["neuAdd_pt"] = {};
              (*outputVecNtuplizer)["neuAdd_eta"] = {};
              (*outputVecNtuplizer)["neuAdd_phi"] = {};
              (*outputVecNtuplizer)["neuAdd_energy"] = {};
              (*outputVecNtuplizer)["neuAdd_et2"] = {};
              (*outputVecNtuplizer)["neuAdd_dR_fromMu"] = {};
            }
            (*outputVecNtuplizer)["nNeuAdd"].push_back(N_compatible_neu);
            if(verbose) {cout << "Compatible neutrals: " << N_compatible_neu << endl;}
            for(uint i = 0; i < N_compatible_neu; i++){
              (*outputVecNtuplizer)["neuAdd_pdgId"].push_back(neu_pdgId[i]);
              (*outputVecNtuplizer)["neuAdd_pt"].push_back(neu_pt[i]);
              (*outputVecNtuplizer)["neuAdd_eta"].push_back(neu_eta[i]);
              (*outputVecNtuplizer)["neuAdd_phi"].push_back(neu_phi[i]);
              (*outputVecNtuplizer)["neuAdd_energy"].push_back(neu_energy[i]);
              (*outputVecNtuplizer)["neuAdd_et2"].push_back(neu_et2[i]);
              (*outputVecNtuplizer)["neuAdd_dR_fromMu"].push_back(neu_dR_fromMu[i]);
            }

            auxC = 0;
            for(auto auxN : (*outputVecNtuplizer)["nNeuAdd"]) auxC += auxN;
            if( auxC != int ((*outputVecNtuplizer)["neuAdd_pdgId"].size()) ){
              cout << "Number of neutrals and neutrals details lenght not matching" << endl;
              assert(false);
            }

            n_B++;

            (*outputVecNtuplizer)["mu_trgMu_idx"].push_back(i_trgMu);
            (*outputVecNtuplizer)["mu_charge"].push_back(trgMu.charge());
            AddTLVToOut(vtxu::getTLVfromMuon(trgMu, mass_Mu), string("mu"), &(*outputVecNtuplizer));
            GlobalPoint auxp_PV(bestVtx.position().x(), bestVtx.position().y(), bestVtx.position().z());
            auto dcaPV = vtxu::computeDCA(iSetup, trgMu, auxp_PV);
            (*outputVecNtuplizer)["mu_dca_PV"].push_back(dcaPV.first);
            (*outputVecNtuplizer)["mu_sigdca_PV"].push_back(dcaPV.first/dcaPV.second);


            GlobalPoint auxp(vtxDst->position().x(), vtxDst->position().y(), vtxDst->position().z());
            auto dca = vtxu::computeDCA(iSetup, trgMu, auxp);
            (*outputVecNtuplizer)["mu_dca_vtxDst"].push_back(dca.first);
            (*outputVecNtuplizer)["mu_sigdca_vtxDst"].push_back(dca.first/dca.second);

            auto ip3D = vtxu::computeIP3D(iSetup, trgMu, vtxu::getTLVfromKinPart(Dst) , vtxDst);
            (*outputVecNtuplizer)["mu_IP3D_vtxDst"].push_back(ip3D.first);
            (*outputVecNtuplizer)["mu_sigIP3D_vtxDst"].push_back(ip3D.first/ip3D.second);

            auto dcaT = vtxu::computeDCA(iSetup, trgMu, auxp, 1);
            (*outputVecNtuplizer)["mu_dcaT_vtxDst"].push_back(dcaT.first);
            (*outputVecNtuplizer)["mu_sigdcaT_vtxDst"].push_back(fabs(dcaT.first)/dcaT.second);
            auto dcaL = vtxu::computeDCA(iSetup, trgMu, auxp, 2);
            (*outputVecNtuplizer)["mu_dcaL_vtxDst"].push_back(dcaL.first);
            (*outputVecNtuplizer)["mu_sigdcaL_vtxDst"].push_back(fabs(dcaL.first)/dcaL.second);

            (*outputVecNtuplizer)["mu_isGlobalMuon"].push_back(trgMu.isGlobalMuon());
            (*outputVecNtuplizer)["mu_trackerStandalonePosLocalChi2"].push_back(trgMu.combinedQuality().chi2LocalPosition);
            if ( !trgMu.globalTrack().isNull() )
            {
              (*outputVecNtuplizer)["mu_globalTkNormChi2"].push_back(trgMu.globalTrack()->normalizedChi2());
            }
            else {
              (*outputVecNtuplizer)["mu_globalTkNormChi2"].push_back( 99999. );
            }
            (*outputVecNtuplizer)["mu_kickFinder"].push_back(trgMu.combinedQuality().trkKink);
            (*outputVecNtuplizer)["mu_segmentCompatibility"].push_back(	trgMu.segmentCompatibility() );
            if (i_muon_PFcand == -1) (*outputVecNtuplizer)["mu_lostInnerHits"].push_back( -10 );
            else (*outputVecNtuplizer)["mu_lostInnerHits"].push_back( (*pfCandHandle)[i_muon_PFcand].lostInnerHits() );
            (*outputVecNtuplizer)["mu_looseId"].push_back(trgMu.isLooseMuon());
            (*outputVecNtuplizer)["mu_mediumId"].push_back(trgMu.isMediumMuon());
            (*outputVecNtuplizer)["mu_tightId"].push_back(trgMu.isTightMuon(bestVtx));
            (*outputVecNtuplizer)["mu_softIdBestVtx"].push_back(trgMu.isSoftMuon(bestVtx));

            auto pfIsoR04 = trgMu.pfIsolationR04();
            auto iso04 = pfIsoR04.sumChargedHadronPt + max(pfIsoR04.sumNeutralHadronEt + pfIsoR04.sumPhotonEt - pfIsoR04.sumPUPt/2.,0.0);
            (*outputVecNtuplizer)["mu_db_iso04"].push_back(iso04);
            auto iso04_corrdDst = iso04 - K.pt() - pi.pt() - pis.pt();
            (*outputVecNtuplizer)["mu_db_corr_iso04"].push_back(iso04_corrdDst);

            auto pfIsoR03 = trgMu.pfIsolationR03();
            auto iso03 = pfIsoR03.sumChargedHadronPt + max(pfIsoR03.sumNeutralHadronEt + pfIsoR03.sumPhotonEt - pfIsoR03.sumPUPt/2.,0.0);
            (*outputVecNtuplizer)["mu_db_iso03"].push_back(iso03);

            GlobalPoint auxp2(vtxB->position().x(), vtxB->position().y(), vtxB->position().z());
            dca = vtxu::computeDCA(iSetup, trgMu, auxp2);
            (*outputVecNtuplizer)["mu_dca_vtxDstMu"].push_back(dca.first);
            (*outputVecNtuplizer)["mu_sigdca_vtxDstMu"].push_back(dca.first/dca.second);
            dcaT = vtxu::computeDCA(iSetup, trgMu, auxp2, 1);
            (*outputVecNtuplizer)["mu_dcaT_vtxDstMu"].push_back(dcaT.first);
            (*outputVecNtuplizer)["mu_sigdcaT_vtxDstMu"].push_back(fabs(dcaT.first)/dcaT.second);
            dcaL = vtxu::computeDCA(iSetup, trgMu, auxp2, 2);
            (*outputVecNtuplizer)["mu_dcaL_vtxDstMu"].push_back(dcaL.first);
            (*outputVecNtuplizer)["mu_sigdcaL_vtxDstMu"].push_back(fabs(dcaL.first)/dcaL.second);


            (*outputVecNtuplizer)["sigdxy_K_PV"].push_back(sigdxy_K_PV);
            (*outputVecNtuplizer)["K_norm_chi2"].push_back(K_norm_chi2);
            (*outputVecNtuplizer)["K_N_valid_hits"].push_back(K_N_valid_hits);
            (*outputVecNtuplizer)["K_lostInnerHits"].push_back(K.lostInnerHits());
            AddTLVToOut(vtxu::getTLVfromCand(K, mass_K), string("K"), &(*outputVecNtuplizer));
            (*outputVecNtuplizer)["sigdxy_pi_PV"].push_back(sigdxy_pi_PV);
            (*outputVecNtuplizer)["pi_norm_chi2"].push_back(pi_norm_chi2);
            (*outputVecNtuplizer)["pi_N_valid_hits"].push_back(pi_N_valid_hits);
            (*outputVecNtuplizer)["pi_lostInnerHits"].push_back(pi.lostInnerHits());
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

            (*outputVecNtuplizer)["sigdxy_pis_PV"].push_back(sigdxy_pis_PV);
            (*outputVecNtuplizer)["pis_norm_chi2"].push_back(pis_norm_chi2);
            (*outputVecNtuplizer)["pis_N_valid_hits"].push_back(pis_N_valid_hits);
            (*outputVecNtuplizer)["pis_lostInnerHits"].push_back(pis.lostInnerHits());
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
      free(valid);
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


void B2DstMuDecayTreeProducer::AddTLVToOut(TLorentzVector v, string n, map<string, vector<float>>* outv) {
  (*outv)[n+"_pt"].push_back(v.Pt());
  (*outv)[n+"_pz"].push_back(v.Pz());
  (*outv)[n+"_eta"].push_back(v.Eta());
  (*outv)[n+"_phi"].push_back(v.Phi());
  return;
}

void B2DstMuDecayTreeProducer::updateCounter(int idx, vector<bool> &cF) {
  if(!cF[idx]) {
    cF[idx] = true;
    counters[idx]++;
  }
  return;
}

double B2DstMuDecayTreeProducer::Mass_varM(TLorentzVector p1, double m1, TLorentzVector p2, double m2){
    double E1 = hypot(m1, p1.P());
    double E2 = hypot(m2, p2.P());
    double p1p2 = p1.Pt()*p2.Pt();
    p1p2 *= cos(p1.Phi() - p2.Phi()) + sinh(p1.Eta())*sinh(p2.Eta());

    double M = m1*m1 + m2*m2 + 2*(E1*E2 - p1p2);
    return sqrt(M);
}

DEFINE_FWK_MODULE(B2DstMuDecayTreeProducer);
