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
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include <iostream>
#include <string>
#include <map>

#include "VtxUtils.hh"

using namespace std;

class MuonImpactParameterStudyProducer : public edm::EDProducer {

public:

    explicit MuonImpactParameterStudyProducer(const edm::ParameterSet &iConfig);
    void AddTLVToOut(TLorentzVector, string, map<string, float>*);

    ~MuonImpactParameterStudyProducer() override {};

private:
    virtual void produce(edm::Event&, const edm::EventSetup&);

    // ----------member data ---------------------------
    edm::EDGetTokenT<vector<reco::GenParticle>> PrunedParticlesSrc_;
    edm::EDGetTokenT<vector<pat::Muon>> muonSrc_;
    edm::EDGetTokenT<vector<reco::Vertex>> vtxSrc_;

    int verbose = 0;
};



MuonImpactParameterStudyProducer::MuonImpactParameterStudyProducer(const edm::ParameterSet &iConfig)
{
    verbose = iConfig.getParameter<int>( "verbose" );
    PrunedParticlesSrc_ = consumes<vector<reco::GenParticle>>(edm::InputTag("prunedGenParticles"));
    muonSrc_ = consumes<vector<pat::Muon>>(edm::InputTag("slimmedMuons"));
    vtxSrc_ = consumes<vector<reco::Vertex>> (edm::InputTag("offlineSlimmedPrimaryVertices"));

    produces<map<string, float>>("outputNtuplizer");
    produces<map<string, vector<float>>>("outputVecNtuplizer");
}


void MuonImpactParameterStudyProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    if (verbose) {cout << "-----------  IP Study ----------\n";}

    // Get prunedGenParticles
    edm::Handle<std::vector<reco::GenParticle>> PrunedGenParticlesHandle;
    iEvent.getByToken(PrunedParticlesSrc_, PrunedGenParticlesHandle);
    unsigned int N_PrunedGenParticles = PrunedGenParticlesHandle->size();

    edm::Handle<vector<pat::Muon>> muonsHandle;
    iEvent.getByToken(muonSrc_, muonsHandle);
    unsigned int N_recoMuons = muonsHandle->size();

    edm::Handle<vector<reco::Vertex>> vtxHandle;
    iEvent.getByToken(vtxSrc_, vtxHandle);
    auto primaryVtx = (*vtxHandle)[0];


    // Output collection
    unique_ptr<map<string, float>> outputNtuplizer(new map<string, float>);
    unique_ptr<map<string, vector<float>>> outputVecNtuplizer(new map<string, vector<float>>);
    (*outputNtuplizer)["n_mu"] = 0;

    if(N_recoMuons < 2) {
      iEvent.put(move(outputNtuplizer), "outputNtuplizer");
      iEvent.put(move(outputVecNtuplizer), "outputVecNtuplizer");
      return;
    }

    // Looking for the MC muons
    reco::Candidate::Point interactionPoint(-9999999999, -999, -999);
    vector<reco::GenParticle> muonsMC;
    for(unsigned int i = 0; i < N_PrunedGenParticles; i++) {
      auto p = (*PrunedGenParticlesHandle)[i];
      if(interactionPoint.x() == -9999999999 && p.isHardProcess()) {
        interactionPoint = p.vertex();
        if(verbose) {cout << "[Hard process] " << p.pdgId() << ": " << p.vx() << ", " << p.vy() << endl;}
      }

      if (abs(p.pdgId()) != 13 || p.numberOfDaughters()!=0) continue;
      muonsMC.push_back(p);
    }
    if(verbose) {cout << "MC muons found: " << muonsMC.size() << endl;}

    vector<pat::Muon> good_muonsRECO;
    for(auto muon : (*muonsHandle)) {
      if(muon.innerTrack().isNull()) continue;
      if (!muon.isSoftMuon(primaryVtx)) continue;
      good_muonsRECO.push_back(muon);
    }
    if(verbose) {cout << "Good muons reco found: " << good_muonsRECO.size() << endl;}

    unsigned int n_mu = 0;
    for(auto mu_mc : muonsMC) {
      int i_recoMatch = -1;
      double best_dR = 1e9;
      for(uint i=0; i<good_muonsRECO.size(); i++) {
        auto mu_reco = good_muonsRECO[i];
        double dPhi = vtxu::dPhi(mu_mc.phi(), mu_reco.phi());
        double dEta = mu_mc.eta() - mu_reco.eta();
        double dPt_rel = (mu_mc.pt() - mu_reco.pt())/mu_mc.pt();

        if (fabs(dPhi) < 0.01 && fabs(dEta) < 0.01 && fabs(dPt_rel) < 0.1 ) {
          auto dR = hypot(dPhi, dEta);
          if(dR < best_dR) {
            best_dR = dR;
            i_recoMatch = i;
          }
        }
      }

      if(i_recoMatch != -1) {
        n_mu++;
        auto mu_reco = good_muonsRECO[i_recoMatch];
        (*outputVecNtuplizer)["recoMu_pt"].push_back(mu_reco.pt());
        (*outputVecNtuplizer)["recoMu_eta"].push_back(mu_reco.eta());
        (*outputVecNtuplizer)["recoMu_phi"].push_back(mu_reco.phi());
        (*outputVecNtuplizer)["recoMu_charge"].push_back(mu_reco.charge());
        auto tk = mu_reco.innerTrack();
        auto dxy = tk->dxy(primaryVtx.position());
        (*outputVecNtuplizer)["recoMu_dz"].push_back(tk->dz(primaryVtx.position()));
        (*outputVecNtuplizer)["recoMu_dxy"].push_back(dxy);
        (*outputVecNtuplizer)["recoMu_sigdxy"].push_back(fabs(dxy)/tk->dxyError());


        (*outputVecNtuplizer)["mcMu_pt"].push_back(mu_mc.pt());
        (*outputVecNtuplizer)["mcMu_eta"].push_back(mu_mc.eta());
        (*outputVecNtuplizer)["mcMu_phi"].push_back(mu_mc.phi());
        (*outputVecNtuplizer)["mcMu_charge"].push_back(mu_mc.charge());
        auto mu_mc_IP = vtxu::computeIP(interactionPoint, mu_mc.vertex(), mu_mc.momentum(), true);
        (*outputVecNtuplizer)["mcMu_IP"].push_back(mu_mc_IP);

        double dPhi = vtxu::dPhi(mu_mc.phi(), mu_reco.phi());
        double dEta = mu_mc.eta() - mu_reco.eta();
        double dPt_rel = (mu_mc.pt() - mu_reco.pt())/mu_mc.pt();
        (*outputVecNtuplizer)["delta_pt_rel"].push_back(dPt_rel);
        (*outputVecNtuplizer)["delta_eta"].push_back(dEta);
        (*outputVecNtuplizer)["delta_phi"].push_back(dPhi);
      }
    }
    if(verbose) {cout << "Matches found: " << n_mu << endl;}

    (*outputNtuplizer)["n_mu"] = n_mu;
    iEvent.put(move(outputNtuplizer), "outputNtuplizer");
    iEvent.put(move(outputVecNtuplizer), "outputVecNtuplizer");
    return;
}


DEFINE_FWK_MODULE(MuonImpactParameterStudyProducer);
