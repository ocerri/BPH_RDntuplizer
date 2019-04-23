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

using namespace std;

class MuD0Producer : public edm::EDProducer {

public:

    explicit MuD0Producer(const edm::ParameterSet &iConfig);

    ~MuD0Producer() override {};

private:

    virtual void produce(edm::Event&, const edm::EventSetup&);

    // ----------member data ---------------------------

    edm::EDGetTokenT<vector<pat::PackedCandidate>> PFCandSrc_;
    edm::EDGetTokenT<vector<reco::Vertex>> vtxSrc_;
    edm::EDGetTokenT<vector<pat::Muon>> TrgMuonSrc_;


    // double mass_B0 = 5.27961;
    // double mass_mu = 0.1056583745;
    // double mass_K = 0.493677;
    // double mass_pi = 0.13957018;

    int verbose = 0;
};



MuD0Producer::MuD0Producer(const edm::ParameterSet &iConfig)
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


void MuD0Producer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
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
    // cout << Form("TrgMuPF vtx: %.2f %.2f %.2f\n", trgMu.vx(), trgMu.vy(), trgMu.vz());
    auto vtxMu = trgMu.vertexRef();
    GlobalPoint p_vtxMu(vtxMu->x(), vtxMu->y(), vtxMu->z());
    // auto dca = vtxu::computeDCA(iSetup, trgMu, p_vtxMu);
    // cout << Form("trgMu dca %.3f +/- %.3f\n", dca.first, dca.second);
    // cout << Form("TrgMuPF vtxRef: %.2f %.2f %.2f\n", vtxtrgmu->x(), vtxtrgmu->y(), vtxtrgmu->z());

    // reco::Vertex PV;
    // double DCA_trgMuPV = -1;
    // int N_PV = 0;
    // for(auto v : *vtxHandle){
    //   if(mu.isTightMuon(v)) {
    //     PV = v;
    //     N_PV++;
    //     cout << Form("vtx: %.2f %.2f %.2f", v.x(), v.y(), v.z());
    //     // cout << " PV" << flush;
    //     GlobalPoint vp(v.x(), v.y(), v.z());
    //     auto dca = vtxu::computeDCA(iSetup, trgMu, vp);
    //     cout << Form("  dca %.3f +/- %.3f\n", dca.first, dca.second);
    //   }
    //   // cout << endl;
    // }
    // if(N_PV>1 && verbose) {cout << "[Warning]: " << N_PV << " PV reco\n";}


    int n_K = 0;
    vector<float> n_pi = {};
    vector<float> mass_kpi = {};
    vector<float> dca_kpi_vtxMu = {};
    vector<float> sigdca_kpi_vtxMu = {};
    vector<float> cos_kpi_vtxMu = {};
    vector<float> d_vtxkpi_vtxMu = {};
    vector<float> sigd_vtxkpi_vtxMu = {};

    if (verbose) {cout <<"-------------------- Evt -----------------------\n";}
    // Look for the K+ (321)
    for(auto k : *pfCandHandle) {
      if (!k.hasTrackDetails()) continue;
      //Require a positive charged hadron
      if (k.pdgId() != 211 ) continue;
      // Require to be close to the trigger muon;
      if (fabs(k.dz() - trgMu.dz()) > 1.) continue;
      if (vtxu::dR(k.phi(), trgMu.phi(), k.eta(), trgMu.eta()) > 1.) continue;
      // Require significant displacement from PV
      n_K++;
      if (verbose) {cout << "###-### K number " << n_K << endl;}

      int n_pi_aux = 0;
      for(uint i_pi = 0; i_pi < N_pfCand; ++i_pi) {
        const pat::PackedCandidate & pi = (*pfCandHandle)[i_pi];
        if (!pi.hasTrackDetails()) continue;
        //Require a negative charged hadron
        if (pi.pdgId() != -211 ) continue;
        // Require to be close to the trigger muon;
        if (fabs(pi.dz() - trgMu.dz()) > 1.) continue;
        if (vtxu::dR(pi.phi(), trgMu.phi(), pi.eta(), trgMu.eta()) > 1.) continue;

        // Require to be close to the K;
        if (fabs(pi.dz() - k.dz()) > 1.) continue;
        if (vtxu::dR(pi.phi(), k.phi(), pi.eta(), k.eta()) > 1.) continue;

        //Fit the vertex
        auto D0KinTree = vtxu::FitD0(iSetup, pi, k, false, 0);
        bool accept_pi = false;
        if(D0KinTree->isValid()) {
          D0KinTree->movePointerToTheTop();
          auto D0vtx = D0KinTree->currentDecayVertex();
          auto chi2 = D0vtx->chiSquared();
          auto max_chi2 = TMath::ChisquareQuantile(0.95, D0vtx->degreesOfFreedom());
          if (chi2 < max_chi2) accept_pi = true;
        }
        if(!accept_pi) continue;

        n_pi_aux++;
        if (verbose) {cout << "---> pi number " << n_pi_aux << endl;}
        auto D0_reco = D0KinTree->currentParticle();
        auto D0vtx = D0KinTree->currentDecayVertex();

        auto mass = D0_reco->currentState().mass();
        // auto mass_err = sqrt(D0_reco.kinematicParametersError().matrix()(6,6));
        mass_kpi.push_back(mass);
        auto dca = vtxu::computeDCA(D0_reco->refittedTransientTrack(), p_vtxMu);
        dca_kpi_vtxMu.push_back(dca.first);
        sigdca_kpi_vtxMu.push_back(fabs(dca.first)/dca.second);

        TVector3 dvtx(D0vtx->position().x() - vtxMu->x(),
                      D0vtx->position().y() - vtxMu->y(),
                      D0vtx->position().z() - vtxMu->z()
                     );
        TVector3 pD0(D0_reco->currentState().globalMomentum().x(),
                     D0_reco->currentState().globalMomentum().y(),
                     D0_reco->currentState().globalMomentum().z()
                     );
        double dalpha = dvtx.Angle(pD0);
        cos_kpi_vtxMu.push_back(cos(dalpha));

        auto dvtx_kpiPV = vtxu::vtxsDistance(vtxMu, D0vtx);
        d_vtxkpi_vtxMu.push_back(dvtx_kpiPV.first);
        sigd_vtxkpi_vtxMu.push_back(dvtx_kpiPV.second);

        if (verbose) {
          cout << "Mass: " << mass << endl;
          cout << Form("kpiPV dca: %.3f +/- %.3f\n", dca.first, dca.second);
          cout << Form("dalpha: %.2f (%.2f)\n", dalpha, cos(dalpha));
          cout << Form("dvtx kpi PV: %.3f +/- %.3f\n", dvtx_kpiPV.first, dvtx_kpiPV.second);
        }
      }

      n_pi.push_back(n_pi_aux);
    }

    (*outputNtuplizer)["n_K"] = n_K;

    (*outputVecNtuplizer)["n_pi"] = n_pi;
    (*outputVecNtuplizer)["mass_kpi"] = mass_kpi;
    (*outputVecNtuplizer)["dca_kpi_vtxMu"] = dca_kpi_vtxMu;
    (*outputVecNtuplizer)["sigdca_kpi_vtxMu"] = sigdca_kpi_vtxMu;
    (*outputVecNtuplizer)["d_vtxkpi_vtxMu"] = d_vtxkpi_vtxMu;
    (*outputVecNtuplizer)["sigd_vtxkpi_vtxMu"] = sigd_vtxkpi_vtxMu;
    (*outputVecNtuplizer)["cos_kpi_vtxMu"] = cos_kpi_vtxMu;


    iEvent.put(move(outputNtuplizer), "outputNtuplizer");
    iEvent.put(move(outputVecNtuplizer), "outputVecNtuplizer");
    return;
}


DEFINE_FWK_MODULE(MuD0Producer);
