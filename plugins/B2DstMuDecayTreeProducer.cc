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

#include "VtxUtils.hh"

using namespace std;

class B2DstMuDecayTreeProducer : public edm::EDProducer {

public:

    explicit B2DstMuDecayTreeProducer(const edm::ParameterSet &iConfig);

    ~B2DstMuDecayTreeProducer() override {};

private:

    virtual void produce(edm::Event&, const edm::EventSetup&);

    // ----------member data ---------------------------

    edm::EDGetTokenT<vector<pat::PackedCandidate>> PFCandSrc_;
    edm::EDGetTokenT<vector<pat::Muon>> TrgMuonSrc_;


    double mass_B0 = 5.27961;
    double mass_mu = 0.1056583745;
    double mass_K = 0.493677;
    double mass_pi = 0.13957018;

    int verbose = 0;
};



B2DstMuDecayTreeProducer::B2DstMuDecayTreeProducer(const edm::ParameterSet &iConfig)
{
    verbose = iConfig.getParameter<int>( "verbose" );
    PFCandSrc_ = consumes<vector<pat::PackedCandidate>>(edm::InputTag("packedPFCandidates"));
    // TrgMuonSrc_ = consumes<vector<pat::Muon>>(edm::InputTag("trgMuonsMatched", "BPHTriggerPathProducer"));
    TrgMuonSrc_ = consumes<vector<pat::Muon>> ( iConfig.getParameter<edm::InputTag>( "trgMuons" ) );

    // produces<vector<pat::PackedCandidate>>("recoCandMatched");
    // produces<vector<string>>("recoCandMatchedNames");
    produces<map<string, float>>("outputNtuplizer");
    produces<map<string, vector<float>>>("outputVecNtuplizer");
}


void B2DstMuDecayTreeProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    edm::Handle<vector<pat::PackedCandidate>> pfCandHandle;
    iEvent.getByToken(PFCandSrc_, pfCandHandle);

    edm::Handle<vector<pat::Muon>> trgMuonsHandle;
    iEvent.getByToken(TrgMuonSrc_, trgMuonsHandle);

    // Output collection
    unique_ptr<map<string, float>> outputNtuplizer(new map<string, float>);
    unique_ptr<map<string, vector<float>>> outputVecNtuplizer(new map<string, vector<float>>);
    // unique_ptr<vector<pat::PackedCandidate>> RECO_MCmatch( new vector<pat::PackedCandidate> );
    // unique_ptr<vector<string>> RECO_MCmatchNames( new vector<string> );

    auto mu = (*trgMuonsHandle)[0];
    pat::PackedCandidate trgMu;
    for (auto p : *pfCandHandle) {
      if (fabs(mu.pt() - p.pt())/mu.pt() > 5e-3) continue;
      if (fabs(mu.eta() - p.eta()) > 5e-3) continue;
      if (fabs(vtxu::dPhi(mu.phi(), p.phi())) > 5e-3) continue;
      trgMu = p;
      break;
    }

    int n_K = 0;
    vector<float> n_pi;

    // Look for the K+ (321)
    for(auto k : *pfCandHandle) {
      //Require a positive charged hadron
      if (k.pdgId() != 211 ) continue;
      // Require to be close to the trigger muon;
      if (fabs(k.dz() - trgMu.dz()) > 1.) continue;
      if (vtxu::dR(k.phi(), trgMu.phi(), k.eta(), trgMu.eta()) > 1.) continue;
      n_K++;

      int n_pi_aux = 0;
      for(auto pi : *pfCandHandle) {
        //Require a negative charged hadron
        if (pi.pdgId() != -211 ) continue;
        // Require to be close to the trigger muon;
        if (fabs(pi.dz() - trgMu.dz()) > 1.) continue;
        if (vtxu::dR(pi.phi(), trgMu.phi(), pi.eta(), trgMu.eta()) > 1.) continue;
        n_pi_aux++;
      }
      n_pi.push_back(n_pi_aux);
    }

    (*outputNtuplizer)["n_K"] = n_K;
    cout << n_K << endl;
    (*outputVecNtuplizer)["n_pi"] = n_pi;
    cout << n_pi.size() << endl;
    iEvent.put(move(outputNtuplizer), "outputNtuplizer");
    iEvent.put(move(outputVecNtuplizer), "outputVecNtuplizer");
    return;
}


DEFINE_FWK_MODULE(B2DstMuDecayTreeProducer);
