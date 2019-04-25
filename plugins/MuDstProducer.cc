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

#define __ProbChi2VtxCut__ 0.99
#define __dzMax__ 1.0
#define __dRMax__ 1.0
#define __cos_kpi_vtxMu_min__ 0.97
#define __d_vtxkpi_vtxMu_min__ 0.02
#define __dmD0_max__ 0.039 // 3*0.013 (i.e. 3 sigma)
#define __cos_D0pis_vtxMu_min__ 0.95
#define __dmKpi_D0pis_max__ 0.010 //Just reasonable for the moment


using namespace std;

class MuDstProducer : public edm::EDProducer {

public:

    explicit MuDstProducer(const edm::ParameterSet &iConfig);

    ~MuDstProducer() override {};

private:

    virtual void produce(edm::Event&, const edm::EventSetup&);

    // ----------member data ---------------------------

    edm::EDGetTokenT<vector<pat::PackedCandidate>> PFCandSrc_;
    edm::EDGetTokenT<vector<reco::Vertex>> vtxSrc_;
    edm::EDGetTokenT<vector<pat::Muon>> TrgMuonSrc_;

    double mass_D0 = 1.86483;

    int verbose = 0;
};



MuDstProducer::MuDstProducer(const edm::ParameterSet &iConfig)
{
    verbose = iConfig.getParameter<int>( "verbose" );
    PFCandSrc_ = consumes<vector<pat::PackedCandidate>>(edm::InputTag("packedPFCandidates"));
    vtxSrc_ = consumes<vector<reco::Vertex>>(edm::InputTag("offlineSlimmedPrimaryVertices"));
    // TrgMuonSrc_ = consumes<vector<pat::Muon>>(edm::InputTag("trgMuonsMatched", "BPHTriggerPathProducer"));
    TrgMuonSrc_ = consumes<vector<pat::Muon>> ( iConfig.getParameter<edm::InputTag>( "trgMuons" ) );

    // produces<vector<pat::PackedCandidate>>("recoCandMatched");
    // produces<vector<string>>("recoCandMatchedNames");
    produces<map<string, float>>("outputNtuplizer");
    produces<map<string, vector<float>>>("outputVecNtuplizer");
}


void MuDstProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
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
    // cout << Form("TrgMu vtx: %.2f %.2f %.2f\n", mu.vx(), mu.vy(), mu.vz());
    pat::PackedCandidate trgMu;
    for (auto p : *pfCandHandle) {
      if (fabs(mu.pt() - p.pt())/mu.pt() > 5e-3) continue;
      if (fabs(mu.eta() - p.eta()) > 5e-3) continue;
      if (fabs(vtxu::dPhi(mu.phi(), p.phi())) > 5e-3) continue;
      trgMu = p;
      break;
    }

    auto vtxMu = trgMu.vertexRef();
    GlobalPoint p_vtxMu(vtxMu->x(), vtxMu->y(), vtxMu->z());


    int n_K = 0, n_pi = 0, n_D0 = 0, n_pis = 0, n_Dst = 0;
    (*outputVecNtuplizer)["chi2_kpi"] = {};
    (*outputVecNtuplizer)["mass_kpi"] = {};
    (*outputVecNtuplizer)["dca_kpi_vtxMu"] = {};
    (*outputVecNtuplizer)["sigdca_kpi_vtxMu"] = {};
    (*outputVecNtuplizer)["cos_kpi_vtxMu"] = {};
    (*outputVecNtuplizer)["d_vtxkpi_vtxMu"] = {};
    (*outputVecNtuplizer)["sigd_vtxkpi_vtxMu"] = {};
    (*outputVecNtuplizer)["chi2_D0pis"] = {};
    (*outputVecNtuplizer)["mass_D0pis"] = {};
    (*outputVecNtuplizer)["dca_D0pis_vtxMu"] = {};
    (*outputVecNtuplizer)["sigdca_D0pis_vtxMu"] = {};
    (*outputVecNtuplizer)["cos_D0pis_vtxMu"] = {};
    (*outputVecNtuplizer)["d_vtxD0pis_vtxMu"] = {};
    (*outputVecNtuplizer)["sigd_vtxD0pis_vtxMu"] = {};

    if (verbose) {cout <<"-------------------- Evt -----------------------\n";}
    // Look for the K+ (321)
    for(auto k : *pfCandHandle) {
      if (!k.hasTrackDetails()) continue;
      //Require a positive charged hadron
      if (k.pdgId() != 211 ) continue;
      // Require to be close to the trigger muon;
      if (fabs(k.dz() - trgMu.dz()) > __dzMax__) continue;
      if (vtxu::dR(k.phi(), trgMu.phi(), k.eta(), trgMu.eta()) > __dRMax__) continue;
      // Require significant displacement from PV
      if (verbose) {cout << "###-### K number " << n_K << endl;}

      n_K++;
      for(uint i_pi = 0; i_pi < N_pfCand; ++i_pi) {
        const pat::PackedCandidate & pi = (*pfCandHandle)[i_pi];
        if (!pi.hasTrackDetails()) continue;
        //Require a negative charged hadron
        if (pi.pdgId() != -211 ) continue;
        // Require to be close to the trigger muon;
        if (fabs(pi.dz() - trgMu.dz()) > __dzMax__) continue;
        if (vtxu::dR(pi.phi(), trgMu.phi(), pi.eta(), trgMu.eta()) > __dRMax__) continue;

        // Require to be close to the K;
        if (fabs(pi.dz() - k.dz()) > __dzMax__) continue;
        if (vtxu::dR(pi.phi(), k.phi(), pi.eta(), k.eta()) > __dRMax__) continue;

        //Fit the vertex
        auto D0KinTree = vtxu::FitD0(iSetup, pi, k, false, 0);
        bool accept_pi = false;
        double chi2_kpi;
        if(D0KinTree->isValid()) {
          D0KinTree->movePointerToTheTop();
          auto D0vtx = D0KinTree->currentDecayVertex();
          chi2_kpi = D0vtx->chiSquared();
          auto max_chi2 = TMath::ChisquareQuantile(__ProbChi2VtxCut__, D0vtx->degreesOfFreedom());
          if (chi2_kpi < max_chi2) accept_pi = true;
        }
        if(!accept_pi) continue;
        n_pi++;

        auto D0_reco = D0KinTree->currentParticle();
        auto D0vtx = D0KinTree->currentDecayVertex();

        auto mass_kpi = D0_reco->currentState().mass();

        auto dca = vtxu::computeDCA(D0_reco->refittedTransientTrack(), p_vtxMu);
        auto dca_kpi_vtxMu = dca.first;
        auto sigdca_kpi_vtxMu = fabs(dca.first)/dca.second;

        TVector3 dvtx(D0vtx->position().x() - vtxMu->x(),
                      D0vtx->position().y() - vtxMu->y(),
                      D0vtx->position().z() - vtxMu->z()
                     );
        TVector3 pD0(D0_reco->currentState().globalMomentum().x(),
                     D0_reco->currentState().globalMomentum().y(),
                     D0_reco->currentState().globalMomentum().z()
                     );
        double dalpha = dvtx.Angle(pD0);
        auto cos_kpi_vtxMu = cos(dalpha);

        auto dvtx_kpiPV = vtxu::vtxsDistance(vtxMu, D0vtx);
        auto d_vtxkpi_vtxMu = dvtx_kpiPV.first;
        auto sigd_vtxkpi_vtxMu = d_vtxkpi_vtxMu/dvtx_kpiPV.second;

        // Make selection cut to get good D0
        bool accept_D0 = cos_kpi_vtxMu > __cos_kpi_vtxMu_min__;
        accept_D0 &= d_vtxkpi_vtxMu > __d_vtxkpi_vtxMu_min__;
        accept_D0 &= fabs(mass_kpi - mass_D0) < __dmD0_max__;
        if(!accept_D0) continue;
        n_D0++;

        // Look for the soft pion to make a Dst-
        for(uint i_pis = 0; i_pis < N_pfCand; ++i_pis) {
          // Pion different from the pion from D0
          if(i_pis==i_pi) continue;

          const pat::PackedCandidate & pis = (*pfCandHandle)[i_pis];
          if (!pis.hasTrackDetails()) continue;
          //Require a negative charged hadron
          if (pis.pdgId() != -211 ) continue;

          // Require to be close to the trigger muon;
          if (fabs(pis.dz() - trgMu.dz()) > __dzMax__) continue;
          if (vtxu::dR(pis.phi(), trgMu.phi(), pis.eta(), trgMu.eta()) > __dRMax__) continue;

          // Fit the Dst vertex
          //  Refitting the D0 with mass constrint
          auto DstKinTree = vtxu::FitDst_fitD0wMassConstraint(iSetup, pis, pi, k, false, 0);

          //  Without refitting the D0 with mass constrint
          // Does not work for unknown reasons
          // auto DstKinTree = vtxu::FitDst(iSetup, pis, D0KinTree->currentParticle(), false, 0);

          bool accept_pis = false;
          double chi2_D0pis;
          if(DstKinTree->isValid()) {
            DstKinTree->movePointerToTheTop();
            auto Dstvtx = DstKinTree->currentDecayVertex();
            chi2_D0pis = Dstvtx->chiSquared();
            auto max_chi2 = TMath::ChisquareQuantile(__ProbChi2VtxCut__, Dstvtx->degreesOfFreedom());
            if (chi2_D0pis < max_chi2) accept_pis = true;
          }
          if(!accept_pis) continue;
          n_pis++;

          auto Dst_reco = DstKinTree->currentParticle();
          auto Dstvtx = DstKinTree->currentDecayVertex();

          auto mass_D0pis = Dst_reco->currentState().mass();

          auto dca = vtxu::computeDCA(Dst_reco->refittedTransientTrack(), p_vtxMu);
          auto dca_D0pis_vtxMu = dca.first;
          auto sigdca_D0pis_vtxMu = fabs(dca.first)/dca.second;

          TVector3 dvtx(Dstvtx->position().x() - vtxMu->x(),
                        Dstvtx->position().y() - vtxMu->y(),
                        Dstvtx->position().z() - vtxMu->z()
                       );
          TVector3 pDst(Dst_reco->currentState().globalMomentum().x(),
                       Dst_reco->currentState().globalMomentum().y(),
                       Dst_reco->currentState().globalMomentum().z()
                       );
          double dalpha = dvtx.Angle(pDst);
          auto cos_D0pis_vtxMu = cos(dalpha);

          auto dvtx_D0pisPV = vtxu::vtxsDistance(vtxMu, Dstvtx);
          auto d_vtxD0pis_vtxMu = dvtx_D0pisPV.first;
          auto sigd_vtxD0pis_vtxMu = d_vtxD0pis_vtxMu/dvtx_D0pisPV.second;

          // Make selection cut to get good Dst
          bool accept_Dst = cos_D0pis_vtxMu > __cos_D0pis_vtxMu_min__;
          accept_Dst &= fabs(mass_D0pis - mass_kpi) < __dmKpi_D0pis_max__;

          if (accept_Dst) {
            n_Dst++;
          }

          (*outputVecNtuplizer)["chi2_kpi"].push_back(chi2_kpi);
          (*outputVecNtuplizer)["mass_kpi"].push_back(mass_kpi);
          (*outputVecNtuplizer)["dca_kpi_vtxMu"].push_back(dca_kpi_vtxMu);
          (*outputVecNtuplizer)["sigdca_kpi_vtxMu"].push_back(sigdca_kpi_vtxMu);
          (*outputVecNtuplizer)["cos_kpi_vtxMu"].push_back(cos_kpi_vtxMu);
          (*outputVecNtuplizer)["d_vtxkpi_vtxMu"].push_back(d_vtxkpi_vtxMu);
          (*outputVecNtuplizer)["sigd_vtxkpi_vtxMu"].push_back(sigd_vtxkpi_vtxMu);
          (*outputVecNtuplizer)["chi2_D0pis"].push_back(chi2_D0pis);
          (*outputVecNtuplizer)["mass_D0pis"].push_back(mass_D0pis);
          (*outputVecNtuplizer)["dca_D0pis_vtxMu"].push_back(dca_D0pis_vtxMu);
          (*outputVecNtuplizer)["sigdca_D0pis_vtxMu"].push_back(sigdca_D0pis_vtxMu);
          (*outputVecNtuplizer)["cos_D0pis_vtxMu"].push_back(cos_D0pis_vtxMu);
          (*outputVecNtuplizer)["d_vtxD0pis_vtxMu"].push_back(d_vtxD0pis_vtxMu);
          (*outputVecNtuplizer)["sigd_vtxD0pis_vtxMu"].push_back(sigd_vtxD0pis_vtxMu);
        }


      }

    }

    (*outputNtuplizer)["n_K"] = n_K;
    (*outputNtuplizer)["n_pi"] = n_pi;
    (*outputNtuplizer)["n_D0"] = n_D0;
    (*outputNtuplizer)["n_pis"] = n_pis;
    (*outputNtuplizer)["n_Dst"] = n_Dst;


    iEvent.put(move(outputNtuplizer), "outputNtuplizer");
    iEvent.put(move(outputVecNtuplizer), "outputVecNtuplizer");
    return;
}


DEFINE_FWK_MODULE(MuDstProducer);
