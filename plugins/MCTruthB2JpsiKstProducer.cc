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
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include <iostream>
#include <string>
#include <map>

#include "VtxUtils.hh"

using namespace std;

class MCTruthB2JpsiKstProducer : public edm::EDProducer {

public:

    explicit MCTruthB2JpsiKstProducer(const edm::ParameterSet &iConfig);
    void AddTLVToOut(TLorentzVector, string, map<string, float>*);

    ~MCTruthB2JpsiKstProducer() override {};

private:
    virtual void produce(edm::Event&, const edm::EventSetup&);

    // ----------member data ---------------------------
    edm::EDGetTokenT<vector<reco::GenParticle>> PrunedParticlesSrc_;
    edm::EDGetTokenT<vector<pat::PackedGenParticle>> PackedParticlesSrc_;
    edm::EDGetTokenT<vector<pat::Muon>> TrgMuonSrc_;

    edm::EDGetTokenT<vector<pat::PackedCandidate>> PFCandSrc_;

    int verbose = 0;
};



MCTruthB2JpsiKstProducer::MCTruthB2JpsiKstProducer(const edm::ParameterSet &iConfig)
{
    verbose = iConfig.getParameter<int>( "verbose" );
    PrunedParticlesSrc_ = consumes<vector<reco::GenParticle>>(edm::InputTag("prunedGenParticles"));
    PackedParticlesSrc_ = consumes<vector<pat::PackedGenParticle>>(edm::InputTag("packedGenParticles"));
    TrgMuonSrc_ = consumes<vector<pat::Muon>> ( iConfig.getParameter<edm::InputTag>( "trgMuons" ) );

    PFCandSrc_ = consumes<vector<pat::PackedCandidate>>(edm::InputTag("packedPFCandidates"));

    produces<map<string, float>>("outputNtuplizer");
    produces<int>("indexBmc");
}


void MCTruthB2JpsiKstProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    if (verbose) {cout << "-----------  MC Truth ----------\n";}

    // Get prunedGenParticles
    edm::Handle<std::vector<reco::GenParticle>> PrunedGenParticlesHandle;
    iEvent.getByToken(PrunedParticlesSrc_, PrunedGenParticlesHandle);
    unsigned int N_PrunedGenParticles = PrunedGenParticlesHandle->size();

    // Get packedGenParticles
    edm::Handle<std::vector<pat::PackedGenParticle>> PackedGenParticlesHandle;
    iEvent.getByToken(PackedParticlesSrc_, PackedGenParticlesHandle);
    // unsigned int N_packedGenParticles = PackedGenParticlesHandle->size();

    // Get trigger muon
    edm::Handle<vector<pat::Muon>> trgMuonsHandle;
    iEvent.getByToken(TrgMuonSrc_, trgMuonsHandle);
    auto trgMu = (*trgMuonsHandle)[0];

    // Get PF cand
    edm::Handle<vector<pat::PackedCandidate>> pfCandHandle;
    iEvent.getByToken(PFCandSrc_, pfCandHandle);

    // Output collection
    unique_ptr<map<string, float>> outputNtuplizer(new map<string, float>);
    unique_ptr<int> indexBmc(new int);


    if(verbose) {
      uint nb = 0;
      for(auto p : *PrunedGenParticlesHandle) {
        if (abs(p.pdgId()) == 10423) {
          cout << "Found: " << p.pdgId() << " (-->[" << p.numberOfDaughters() << "]" << flush;
          for(auto d : p.daughterRefVector()) {
            cout << " " << d->pdgId() << flush;
          }
          cout << ')' << endl;
        }
        if ((abs(p.pdgId()) == 511 || abs(p.pdgId()) == 521  || abs(p.pdgId()) == 531  || abs(p.pdgId()) == 541) && p.numberOfDaughters()>1) {
          nb++;
          cout << "Found: " << p.pdgId() << " (-->" << flush;
          for(auto d : p.daughterRefVector()) {
            cout << " " << d->pdgId() << flush;
          }
          cout << ')' << endl;
        }
      }
      cout << "Number of B found: " << nb << endl;
    }

    // Looking for the B0 -> J/psi K*0
    map<string, TLorentzVector> p4;
    p4["mu"] = TLorentzVector();
    double mu_impactParam = -1;
    int i_B = -1;
    reco::Candidate::Point interactionPoint(-9999999999, -999, -999);
    for(unsigned int i = 0; i < N_PrunedGenParticles; i++) {
      auto p = (*PrunedGenParticlesHandle)[i];
      if(p.isHardProcess() && interactionPoint.x() == -9999999999) {
        interactionPoint = p.vertex();
        if(verbose) {cout << "[Hard process] " << p.pdgId() << ": " << p.vx() << ", " << p.vy() << endl;}
      }

      if (p.pdgId() == 511 && p.numberOfDaughters()>1) {
        bool Jpsi_found = false;
        bool Kst_found = false;
        for(auto d : p.daughterRefVector()) {
          if(d->pdgId() == 443) Jpsi_found = true;
          else if(d->pdgId() == 313) Kst_found = true;
        }
        if(Jpsi_found && Kst_found) i_B = i;
      }

      bool semilepDecay = p.pdgId() == -511 || abs(p.pdgId()) == 521 || abs(p.pdgId()) == 531;
      if (semilepDecay && p.numberOfDaughters()>1) {
        for(auto d : p.daughterRefVector()) {
          if(abs(d->pdgId()) == 13) {
            p4["mu"].SetPtEtaPhiM(d->pt(), d->eta(), d->phi(), d->mass());
            mu_impactParam = vtxu::computeIP(interactionPoint, d->vertex(), d->momentum(), true);
            if(verbose) {
              cout << "Muon from " << p.pdgId() << Form(" found (pt: %.1f, eta: %.1f, phi: %.1f, ip: %f)", d->pt(), d->eta(), d->phi(), mu_impactParam) << endl;
            }
          }
        }
      }
    }
    (*indexBmc) = i_B;

    p4["B"] = TLorentzVector();
    p4["Jpsi"] = TLorentzVector();
    p4["mup"] = TLorentzVector();
    p4["mum"] = TLorentzVector();
    p4["Kst"] = TLorentzVector();
    p4["pi"] = TLorentzVector();
    p4["K"] = TLorentzVector();

    if(i_B >= 0){
      auto p = (*PrunedGenParticlesHandle)[i_B];
      p4["B"].SetPtEtaPhiM(p.pt(), p.eta(), p.phi(), p.mass());
      for(auto d : p.daughterRefVector()) {
        if(d->pdgId() == 313) {
          p4["Kst"].SetPtEtaPhiM(d->pt(), d->eta(), d->phi(), d->mass());
          for (auto dd : (*PackedGenParticlesHandle)) {
            if (dd.numberOfMothers() == 0) continue;
            auto m = dd.mother(0);
            if(m->pdgId() == 313) {
              bool match = d->pt()==m->pt() && d->eta()==m->eta() && d->phi()==m->phi();
              if (!match) continue;
              if(dd.pdgId() == -211) {
                p4["pi"].SetPtEtaPhiM(dd.pt(), dd.eta(), dd.phi(), dd.mass());
              }
              else if(dd.pdgId() == 321) {
                p4["K"].SetPtEtaPhiM(dd.pt(), dd.eta(), dd.phi(), dd.mass());
              }
            }
          }
        }
        else if(d->pdgId() == 443) {
          p4["Jpsi"].SetPtEtaPhiM(d->pt(), d->eta(), d->phi(), d->mass());
          for(auto dd : d->daughterRefVector()) {
            if(dd->pdgId() == 13) {
              p4["mum"].SetPtEtaPhiM(dd->pt(), dd->eta(), dd->phi(), dd->mass());
            }
            else if(dd->pdgId() == -13) {
              p4["mup"].SetPtEtaPhiM(dd->pt(), dd->eta(), dd->phi(), dd->mass());
            }
          }
        }
      }
    }

    for(auto kv : p4) {
      AddTLVToOut(kv.second, "MC_"+kv.first, &(*outputNtuplizer));
      if(verbose) {
        bool matched = false;
        size_t i_c = -1;
        if (kv.first=="K" || kv.first=="pi") {
          for (size_t ip=0; ip < pfCandHandle->size(); ip++) {
            auto p = (*pfCandHandle)[ip];
            matched = fabs(p.pt() - kv.second.Pt())/p.pt() < 0.2;
            matched &= vtxu::dR(p.phi(), kv.second.Phi(), p.eta(), kv.second.Eta()) < 0.2;
            if(matched) {
              i_c = ip;
              break;
            }
          }
        }
        cout << kv.first << Form(": %.2f %.2f %.2f", kv.second.Pt(), kv.second.Eta(), kv.second.Phi());
        if (kv.first=="K" || kv.first=="pi") {
          cout << " (" << matched << ")";
          if (matched) {
            auto p = (*pfCandHandle)[i_c];
            cout << endl << "\t" << i_c << Form(": %.2f %.2f %.2f", p.pt(), p.eta(), p.phi()) << endl;
            bool antiMuID = !p.isTrackerMuon() && !p.isStandAloneMuon();
            cout << "\t   hasTrack: " << p.hasTrackDetails() << " antiMuID:" << antiMuID << " charge:" << p.charge();
          }
        }
        cout << endl;
      }
    }
    (*outputNtuplizer)["MC_mu_IP"] = mu_impactParam;



    iEvent.put(move(outputNtuplizer), "outputNtuplizer");
    iEvent.put(move(indexBmc), "indexBmc");
    return;
}


void MCTruthB2JpsiKstProducer::AddTLVToOut(TLorentzVector v, string n, map<string, float>* outv) {
  (*outv)[n+"_pt"] = v.Pt();
  (*outv)[n+"_eta"] = v.Eta();
  (*outv)[n+"_phi"] = v.Phi();
  (*outv)[n+"_P"] = v.P();
  return;
}


DEFINE_FWK_MODULE(MCTruthB2JpsiKstProducer);
