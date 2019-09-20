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

#define __PvalChi2Vtx_max__ 0.99 // Very loose cut
#define __dzMax__ 1.0
#define __dRMax__ 2.0
#define __cos_Kpi_vtxMu_min__ 0.97 // Optimized in D0 fitting
#define __d_vtxKpi_vtxMu_min__ 0.02 // Optimized in D0 fitting
#define __dmD0_max__ 0.039 // 3*0.013 GeV(i.e. 3 sigma) form the D0 mass
#define __dmD0pis_max__ 0.0024 // 3*0.8 MeV (i.e 3 inflated sigma) from Dst mass
#define __pT_Ks_min__ 0.3 // Studied in nb to reduce fake
#define __mDstK_max__ 6.8 // Some reasonable cut on the mass
#define __mDstK_min__ 3.8 // Some reasonable cut on the mass
#define __cos_DstK_vtxMu_min__ 0.8 // Some loose cut tuned on MC


using namespace std;

class B2DstKDecayTreeProducer : public edm::EDProducer {

public:

    explicit B2DstKDecayTreeProducer(const edm::ParameterSet &iConfig);
    void DeclareTLVToOut(string, map<string, vector<float>>*);
    void AddTLVToOut(TLorentzVector, string, map<string, vector<float>>*);

    ~B2DstKDecayTreeProducer() override {};

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



B2DstKDecayTreeProducer::B2DstKDecayTreeProducer(const edm::ParameterSet &iConfig)
{
    verbose = iConfig.getParameter<int>( "verbose" );
    PFCandSrc_ = consumes<vector<pat::PackedCandidate>>(edm::InputTag("packedPFCandidates"));
    vtxSrc_ = consumes<vector<reco::Vertex>>(edm::InputTag("offlineSlimmedPrimaryVertices"));
    TrgMuonSrc_ = consumes<vector<pat::Muon>> ( iConfig.getParameter<edm::InputTag>( "trgMuons" ) );

    produces<map<string, float>>("outputNtuplizer");
    produces<map<string, vector<float>>>("outputVecNtuplizer");
}


void B2DstKDecayTreeProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    edm::Handle<vector<pat::PackedCandidate>> pfCandHandle;
    iEvent.getByToken(PFCandSrc_, pfCandHandle);
    unsigned int N_pfCand = pfCandHandle->size();

    edm::Handle<vector<reco::Vertex>> vtxHandle;
    iEvent.getByToken(vtxSrc_, vtxHandle);


    edm::Handle<vector<pat::Muon>> trgMuonsHandle;
    iEvent.getByToken(TrgMuonSrc_, trgMuonsHandle);

    // Output collection
    unique_ptr<map<string, float>> outputNtuplizer(new map<string, float>);
    unique_ptr<map<string, vector<float>>> outputVecNtuplizer(new map<string, vector<float>>);

    auto mu = (*trgMuonsHandle)[0];
    pat::PackedCandidate trgMu;
    uint i_trgMu = 0;
    for (i_trgMu = 0; i_trgMu < N_pfCand; ++i_trgMu) {
      const pat::PackedCandidate & p = (*pfCandHandle)[i_trgMu];
      if (fabs(mu.pt() - p.pt())/mu.pt() > 5e-3) continue;
      if (fabs(mu.eta() - p.eta()) > 5e-3) continue;
      if (fabs(vtxu::dPhi(mu.phi(), p.phi())) > 5e-3) continue;
      trgMu = p;
      break;
    }

    auto vtxMu = trgMu.vertexRef();
    GlobalPoint p_vtxMu(vtxMu->x(), vtxMu->y(), vtxMu->z());


    int n_K = 0, n_pi = 0, n_D0 = 0, n_pis = 0, n_Dst = 0, n_Ks = 0, n_B = 0;

    (*outputVecNtuplizer)["chi2_Kpi"] = {};
    (*outputVecNtuplizer)["mass_Kpi"] = {};
    (*outputVecNtuplizer)["dca_Kpi_vtxMu"] = {};
    (*outputVecNtuplizer)["sigdca_Kpi_vtxMu"] = {};
    (*outputVecNtuplizer)["cos_Kpi_vtxMu"] = {};
    (*outputVecNtuplizer)["d_vtxKpi_vtxMu"] = {};
    (*outputVecNtuplizer)["sigd_vtxKpi_vtxMu"] = {};

    (*outputVecNtuplizer)["chi2_D0pis"] = {};
    (*outputVecNtuplizer)["mass_D0pis"] = {};
    (*outputVecNtuplizer)["dca_D0pis_vtxMu"] = {};
    (*outputVecNtuplizer)["sigdca_D0pis_vtxMu"] = {};
    (*outputVecNtuplizer)["cos_D0pis_vtxMu"] = {};
    (*outputVecNtuplizer)["d_vtxD0pis_vtxMu"] = {};
    (*outputVecNtuplizer)["sigd_vtxD0pis_vtxMu"] = {};

    (*outputVecNtuplizer)["d_Ks_PV"] = {};
    (*outputVecNtuplizer)["sigd_Ks_PV"] = {};
    (*outputVecNtuplizer)["chi2_DstK"] = {};
    (*outputVecNtuplizer)["mass_DstKs"] = {};
    (*outputVecNtuplizer)["cos_DstK_vtxBest"] = {};
    (*outputVecNtuplizer)["d_vtxDstK_vtxBest"] = {};
    (*outputVecNtuplizer)["sigd_vtxDstK_vtxBest"] = {};

    DeclareTLVToOut("B", &(*outputVecNtuplizer));
    DeclareTLVToOut("Ks", &(*outputVecNtuplizer));
    DeclareTLVToOut("Dst", &(*outputVecNtuplizer));
    DeclareTLVToOut("D0", &(*outputVecNtuplizer));
    DeclareTLVToOut("pis", &(*outputVecNtuplizer));
    DeclareTLVToOut("K", &(*outputVecNtuplizer));
    DeclareTLVToOut("pi", &(*outputVecNtuplizer));

    if (verbose) {cout <<"-------------------- Evt -----------------------\n";}
    /*
    ############################################################################
                              Look for the K+ (321)
    ############################################################################
    */
    for(uint i_k = 0; i_k < N_pfCand; ++i_k) {
      if(i_k==i_trgMu) continue;

      const pat::PackedCandidate & k = (*pfCandHandle)[i_k];
      if (!k.hasTrackDetails()) continue;
      //Require a positive charged hadron
      if (k.pdgId() != 211 ) continue;
      // Require to be close to the trigger muon;
      if (fabs(k.dz() - trgMu.dz()) > __dzMax__) continue;
      // Require significant displacement from PV?
      // if (verbose) {cout << "###-### K number " << n_K << endl;}

      n_K++;

      /*
      ############################################################################
                               Look for the pi- (-211)
      ############################################################################
      */
      for(uint i_pi = 0; i_pi < N_pfCand; ++i_pi) {
        if(i_pi==i_trgMu || i_pi==i_k) continue;

        const pat::PackedCandidate & pi = (*pfCandHandle)[i_pi];
        if (!pi.hasTrackDetails()) continue;
        //Require a negative charged hadron
        if (pi.pdgId() != -211 ) continue;
        // Require to be close to the trigger muon;
        if (fabs(pi.dz() - trgMu.dz()) > __dzMax__) continue;

        // Require to be close to the K;
        if (fabs(pi.dz() - k.dz()) > __dzMax__) continue;
        if (vtxu::dR(pi.phi(), k.phi(), pi.eta(), k.eta()) > __dRMax__) continue;

        //Fit the vertex
        auto D0KinTree = vtxu::FitD0(iSetup, pi, k, false, 0);
        bool accept_pi = false;
        double chi2_Kpi;
        if(D0KinTree->isValid()) {
          D0KinTree->movePointerToTheTop();
          auto D0vtx = D0KinTree->currentDecayVertex();
          chi2_Kpi = D0vtx->chiSquared();
          auto max_chi2 = TMath::ChisquareQuantile(__PvalChi2Vtx_max__, D0vtx->degreesOfFreedom());
          if (chi2_Kpi > 0 && chi2_Kpi < max_chi2) accept_pi = true;
        }
        if(!accept_pi) continue;
        n_pi++;

        auto D0 = D0KinTree->currentParticle();
        auto D0vtx = D0KinTree->currentDecayVertex();

        auto mass_Kpi = D0->currentState().mass();

        auto dca = vtxu::computeDCA(D0->refittedTransientTrack(), p_vtxMu);
        auto dca_Kpi_vtxMu = dca.first;
        auto sigdca_Kpi_vtxMu = fabs(dca.first)/dca.second;

        TVector3 dvtx(D0vtx->position().x() - vtxMu->x(),
                      D0vtx->position().y() - vtxMu->y(),
                      D0vtx->position().z() - vtxMu->z()
                     );
        TVector3 pD0(D0->currentState().globalMomentum().x(),
                     D0->currentState().globalMomentum().y(),
                     D0->currentState().globalMomentum().z()
                     );
        double dalpha = dvtx.Angle(pD0);
        auto cos_Kpi_vtxMu = cos(dalpha);

        auto dvtx_KpiPV = vtxu::vtxsDistance(vtxMu, D0vtx);
        auto d_vtxKpi_vtxMu = dvtx_KpiPV.first;
        auto sigd_vtxKpi_vtxMu = fabs(d_vtxKpi_vtxMu/dvtx_KpiPV.second);

        // Make selection cut to get good D0
        bool accept_D0 = cos_Kpi_vtxMu > __cos_Kpi_vtxMu_min__;
        accept_D0 &= d_vtxKpi_vtxMu > __d_vtxKpi_vtxMu_min__;
        accept_D0 &= fabs(mass_Kpi - mass_D0) < __dmD0_max__;
        if(!accept_D0) continue;
        n_D0++;

        /*
        ############################################################################
                         Look for the soft pion to make a Dst-
        ############################################################################
        */
        for(uint i_pis = 0; i_pis < N_pfCand; ++i_pis) {
          // Pion different from the pion from D0
          if(i_pis==i_trgMu || i_pis==i_k || i_pis==i_pi) continue;

          const pat::PackedCandidate & pis = (*pfCandHandle)[i_pis];
          if (!pis.hasTrackDetails()) continue;
          //Require a negative charged hadron
          if (pis.pdgId() != -211 ) continue;

          // Require to be close to the trigger muon;
          if (fabs(pis.dz() - trgMu.dz()) > __dzMax__) continue;

          auto d_pis_PV = pis.dxy();
          auto sigd_pis_PV = fabs(d_pis_PV/pis.dxyError());
          auto pis_tk = pis.bestTrack();
          auto pis_norm_chi2 = pis_tk->normalizedChi2();
          auto pis_N_valid_hits = pis_tk->numberOfValidHits();

          // Fit the Dst vertex
          //  Refitting the D0 with mass constrint
          auto DstKinTree = vtxu::FitDst_fitD0wMassConstraint(iSetup, pis, pi, k, false, 0);
          //  Without refitting the D0 with mass constrint
          // auto DstKinTree = vtxu::FitDst(iSetup, pis, D0, false, 0);
          bool accept_pis = false;
          double chi2_D0pis;
          if(DstKinTree->isValid()) {
            DstKinTree->movePointerToTheTop();
            auto Dstvtx = DstKinTree->currentDecayVertex();
            chi2_D0pis = Dstvtx->chiSquared();
            auto max_chi2 = TMath::ChisquareQuantile(__PvalChi2Vtx_max__, Dstvtx->degreesOfFreedom());
            if (chi2_D0pis > 0 && chi2_D0pis < max_chi2) accept_pis = true;
          }
          if(!accept_pis) continue;
          n_pis++;

          auto Dst = DstKinTree->currentParticle();
          auto Dstvtx = DstKinTree->currentDecayVertex();

          auto mass_D0pis = Dst->currentState().mass();

          auto dca = vtxu::computeDCA(Dst->refittedTransientTrack(), p_vtxMu);
          auto dca_D0pis_vtxMu = dca.first;
          auto sigdca_D0pis_vtxMu = fabs(dca.first)/dca.second;

          TVector3 dvtx(Dstvtx->position().x() - vtxMu->x(),
                        Dstvtx->position().y() - vtxMu->y(),
                        Dstvtx->position().z() - vtxMu->z()
                       );
          TVector3 pDst(Dst->currentState().globalMomentum().x(),
                       Dst->currentState().globalMomentum().y(),
                       Dst->currentState().globalMomentum().z()
                       );
          double dalpha = dvtx.Angle(pDst);
          auto cos_D0pis_vtxMu = cos(dalpha);

          auto dvtx_D0pisPV = vtxu::vtxsDistance(vtxMu, Dstvtx);
          auto d_vtxD0pis_vtxMu = dvtx_D0pisPV.first;
          auto sigd_vtxD0pis_vtxMu = fabs(d_vtxD0pis_vtxMu/dvtx_D0pisPV.second);

          // Make selection cut to get good Dst
          bool accept_Dst = fabs(mass_D0pis - mass_Dst) < __dmD0pis_max__;
          if (!accept_Dst) continue;
          n_Dst++;
          if (verbose) {cout << "D* found\n";}
          // Refi the Dst with its mass constraint
          DstKinTree = vtxu::FitDst_fitD0wMassConstraint(iSetup, pis, pi, k, true, 0);
          if(DstKinTree->isValid()) Dst = DstKinTree->currentParticle();
          else continue;

          /*
          ############################################################################
                           Look for the soft K to make a B0
          ############################################################################
          */
          for(uint i_Ks = 0; i_Ks < N_pfCand; ++i_Ks) {
            // Kaon different from the previous condidates
            if(i_Ks==i_trgMu || i_Ks==i_k || i_Ks==i_pi || i_Ks==i_pis) continue;
            const pat::PackedCandidate & Ks = (*pfCandHandle)[i_Ks];
            if (!Ks.hasTrackDetails()) continue;
            //Require a positive charged hadron
            if (Ks.pdgId() != 211 ) continue;

            // Require to be close to the trigger muon;
            if (fabs(Ks.dz() - trgMu.dz()) > __dzMax__) continue;

            if (Ks.pt() < __pT_Ks_min__ ) continue;

            auto d_Ks_PV = Ks.dxy();
            auto sigd_Ks_PV = fabs(d_Ks_PV/Ks.dxyError());
            auto Ks_tk = Ks.bestTrack();
            auto Ks_norm_chi2 = Ks_tk->normalizedChi2();
            auto Ks_N_valid_hits = Ks_tk->numberOfValidHits();

            // Fit the B vertex
            auto BKinTree = vtxu::FitVtxDstK(iSetup, Dst, Ks, 0);
            bool accept_DstK = false;
            double chi2_DstK;
            if(BKinTree->isValid()) {
              BKinTree->movePointerToTheTop();
              auto Bvtx = BKinTree->currentDecayVertex();
              chi2_DstK = Bvtx->chiSquared();
              auto max_chi2 = TMath::ChisquareQuantile(__PvalChi2Vtx_max__, Bvtx->degreesOfFreedom());
              if (chi2_DstK > 0 && chi2_DstK < max_chi2) accept_DstK = true;
            }
            if(!accept_DstK) continue;
            n_Ks++;

            auto DstK = BKinTree->currentParticle();
            auto Bvtx = BKinTree->currentDecayVertex();

            BKinTree->movePointerToTheFirstChild();
            auto refit_K = BKinTree->currentParticle();
            BKinTree->movePointerToTheNextChild();
            auto refit_Dst = BKinTree->currentParticle();

            auto mass_DstK = DstK->currentState().mass();

            TVector3 dvtxB(Bvtx->position().x() - vtxMu->x(),
                               Bvtx->position().y() - vtxMu->y(),
                               Bvtx->position().z() - vtxMu->z()
                              );
            TVector3 pDstK(DstK->currentState().globalMomentum().x(),
                           DstK->currentState().globalMomentum().y(),
                           DstK->currentState().globalMomentum().z()
                          );
            double dalphaB = dvtxB.Angle(pDstK);
            auto cos_DstK_vtxBest = cos(dalphaB);
            auto dvtx_DstK_vtxBest = vtxu::vtxsDistance(vtxMu, Bvtx);

            // Look for the best vertex
            if (verbose) {cout << Form("Muon vtx: d = %.2f  cos = %.3e", dvtxB.Mag(), 1-cos_DstK_vtxBest) << endl;}
            for(auto v : *vtxHandle) {
              if ( fabs(vtxMu->z()-v.z()) < 2 ) {
                TVector3 df(Bvtx->position().x() - v.x(),
                            Bvtx->position().y() - v.y(),
                            Bvtx->position().z() - v.z()
                            );

                double da = df.Angle(pDstK);
                auto cosda = cos(da);
                if (cosda > cos_DstK_vtxBest) {
                  dvtxB = df;
                  cos_DstK_vtxBest = cosda;
                  dvtx_DstK_vtxBest = vtxu::vtxsDistance(v, Bvtx);
                }
              }
            }
            if (verbose) {cout << Form("Best vtx: d = %.2f  cos = %.3e", dvtxB.Mag(), 1-cos_DstK_vtxBest) << endl;}
            auto d_vtxDstK_vtxBest = dvtx_DstK_vtxBest.first;
            auto sigd_vtxDstK_vtxBest = fabs(d_vtxDstK_vtxBest/dvtx_DstK_vtxBest.second);


            bool accept_B = mass_DstK < __mDstK_max__ && mass_DstK > __mDstK_min__;
            accept_B &= cos_DstK_vtxBest > __cos_DstK_vtxMu_min__;
            if(!accept_B) continue;
            if (verbose) {cout << "B -> D*- K+ found\n";}
            n_B++;

            /*
            ############################################################################
                                  Compute analysis variables
            ############################################################################
            */
            (*outputVecNtuplizer)["chi2_Kpi"].push_back(chi2_Kpi);
            (*outputVecNtuplizer)["mass_Kpi"].push_back(mass_Kpi);
            (*outputVecNtuplizer)["dca_Kpi_vtxMu"].push_back(dca_Kpi_vtxMu);
            (*outputVecNtuplizer)["sigdca_Kpi_vtxMu"].push_back(sigdca_Kpi_vtxMu);
            (*outputVecNtuplizer)["cos_Kpi_vtxMu"].push_back(cos_Kpi_vtxMu);
            (*outputVecNtuplizer)["d_vtxKpi_vtxMu"].push_back(d_vtxKpi_vtxMu);
            (*outputVecNtuplizer)["sigd_vtxKpi_vtxMu"].push_back(sigd_vtxKpi_vtxMu);

            (*outputVecNtuplizer)["d_pis_PV"].push_back(d_pis_PV);
            (*outputVecNtuplizer)["sigd_pis_PV"].push_back(sigd_pis_PV);
            (*outputVecNtuplizer)["pis_norm_chi2"].push_back(pis_norm_chi2);
            (*outputVecNtuplizer)["pis_N_valid_hits"].push_back(pis_N_valid_hits);
            (*outputVecNtuplizer)["chi2_D0pis"].push_back(chi2_D0pis);
            (*outputVecNtuplizer)["mass_D0pis"].push_back(mass_D0pis);
            (*outputVecNtuplizer)["dca_D0pis_vtxMu"].push_back(dca_D0pis_vtxMu);
            (*outputVecNtuplizer)["sigdca_D0pis_vtxMu"].push_back(sigdca_D0pis_vtxMu);
            (*outputVecNtuplizer)["cos_D0pis_vtxMu"].push_back(cos_D0pis_vtxMu);
            (*outputVecNtuplizer)["d_vtxD0pis_vtxMu"].push_back(d_vtxD0pis_vtxMu);
            (*outputVecNtuplizer)["sigd_vtxD0pis_vtxMu"].push_back(sigd_vtxD0pis_vtxMu);

            (*outputVecNtuplizer)["d_Ks_PV"].push_back(d_Ks_PV);
            (*outputVecNtuplizer)["sigd_Ks_PV"].push_back(sigd_Ks_PV);
            (*outputVecNtuplizer)["Ks_norm_chi2"].push_back(Ks_norm_chi2);
            (*outputVecNtuplizer)["Ks_N_valid_hits"].push_back(Ks_N_valid_hits);
            (*outputVecNtuplizer)["chi2_DstKs"].push_back(chi2_DstK);
            (*outputVecNtuplizer)["mass_DstKs"].push_back(mass_DstK);
            (*outputVecNtuplizer)["cos_DstKs_vtxBest"].push_back(cos_DstK_vtxBest);
            (*outputVecNtuplizer)["d_vtxDstK_vtxBest"].push_back(d_vtxDstK_vtxBest);
            (*outputVecNtuplizer)["sigd_vtxDstK_vtxBest"].push_back(sigd_vtxDstK_vtxBest);

            AddTLVToOut(vtxu::getTLVfromCand(k, mass_K), string("K"), &(*outputVecNtuplizer));
            AddTLVToOut(vtxu::getTLVfromCand(pi, mass_Pi), string("pi"), &(*outputVecNtuplizer));
            AddTLVToOut(vtxu::getTLVfromCand(pis, mass_Pi), string("pis"), &(*outputVecNtuplizer));
            AddTLVToOut(vtxu::getTLVfromKinPart(D0), string("D0"), &(*outputVecNtuplizer));
            AddTLVToOut(vtxu::getTLVfromKinPart(DstK), string("B"), &(*outputVecNtuplizer));
            AddTLVToOut(vtxu::getTLVfromKinPart(refit_Dst), string("Dst"), &(*outputVecNtuplizer));
            AddTLVToOut(vtxu::getTLVfromKinPart(refit_K), string("Ks"), &(*outputVecNtuplizer));
          }
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
    (*outputNtuplizer)["n_Ks"] = n_Ks;
    (*outputNtuplizer)["n_B"] = n_B;


    iEvent.put(move(outputNtuplizer), "outputNtuplizer");
    iEvent.put(move(outputVecNtuplizer), "outputVecNtuplizer");
    return;
}


void B2DstKDecayTreeProducer::DeclareTLVToOut(string n, map<string, vector<float>>* outv) {
  (*outv)[n+"_pt"] = {};
  (*outv)[n+"_pz"] = {};
  (*outv)[n+"_eta"] = {};
  (*outv)[n+"_phi"] = {};
  (*outv)[n+"_P"] = {};
  (*outv)[n+"_E"] = {};
  return;
}


void B2DstKDecayTreeProducer::AddTLVToOut(TLorentzVector v, string n, map<string, vector<float>>* outv) {
  (*outv)[n+"_pt"].push_back(v.Pt());
  (*outv)[n+"_pz"].push_back(v.Pz());
  (*outv)[n+"_eta"].push_back(v.Eta());
  (*outv)[n+"_phi"].push_back(v.Phi());
  (*outv)[n+"_P"].push_back(v.P());
  (*outv)[n+"_E"].push_back(v.E());
  return;
}

DEFINE_FWK_MODULE(B2DstKDecayTreeProducer);
