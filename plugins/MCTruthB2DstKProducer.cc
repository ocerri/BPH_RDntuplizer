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

class MCTruthB2DstKProducer : public edm::EDProducer {

public:

    explicit MCTruthB2DstKProducer(const edm::ParameterSet &iConfig);
    void AddTLVToOut(TLorentzVector, string, map<string, float>*);

    ~MCTruthB2DstKProducer() override {};

private:
    virtual void produce(edm::Event&, const edm::EventSetup&);

    // ----------member data ---------------------------
    edm::EDGetTokenT<vector<reco::GenParticle>> PrunedParticlesSrc_;
    edm::EDGetTokenT<vector<pat::Muon>> TrgMuonSrc_;

    double mass_B0 = 5.27961;
    double mass_mu = 0.1056583745;
    double mass_K = 0.493677;
    double mass_pi = 0.13957018;

    int verbose = 0;
};



MCTruthB2DstKProducer::MCTruthB2DstKProducer(const edm::ParameterSet &iConfig)
{
    verbose = iConfig.getParameter<int>( "verbose" );
    PrunedParticlesSrc_ = consumes<vector<reco::GenParticle>>(edm::InputTag("prunedGenParticles"));
    TrgMuonSrc_ = consumes<vector<pat::Muon>> ( iConfig.getParameter<edm::InputTag>( "trgMuons" ) );

    produces<map<string, float>>("outputNtuplizer");
    produces<int>("indexBmc");
}


void MCTruthB2DstKProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    if (verbose) {cout << "-----------  MC Truth ----------\n";}

    // Get prunedGenParticles
    edm::Handle<std::vector<reco::GenParticle>> PrunedGenParticlesHandle;
    iEvent.getByToken(PrunedParticlesSrc_, PrunedGenParticlesHandle);
    unsigned int N_PrunedGenParticles = PrunedGenParticlesHandle->size();

    // Get trigger muon
    edm::Handle<vector<pat::Muon>> trgMuonsHandle;
    iEvent.getByToken(TrgMuonSrc_, trgMuonsHandle);
    auto trgMu = (*trgMuonsHandle)[0];

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

    // Looking for the B0 -> D*-K+
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
        bool Dst_found = false;
        bool K_found = false;
        for(auto d : p.daughterRefVector()) {
          if(d->pdgId() == -413) Dst_found = true;
          else if(d->pdgId() == 321) K_found = true;
        }
        if(Dst_found && K_found) i_B = i;
      }

      bool semilepDecay = p.pdgId() == -511 || abs(p.pdgId()) == 521 || abs(p.pdgId()) == 531;
      if (semilepDecay && p.numberOfDaughters()>1) {
        for(auto d : p.daughterRefVector()) {
          if(abs(d->pdgId()) == 13) {
            double dR = vtxu::dR(trgMu.phi(), d->phi(), trgMu.eta(), d->eta());
            double dPt_rel_mu = fabs(trgMu.pt() - d->pt())/trgMu.pt();
            if (fabs(dR) < 0.05 && fabs(dPt_rel_mu) < 0.1 ) {
              p4["mu"].SetPtEtaPhiM(d->pt(), d->eta(), d->phi(), d->mass());
              mu_impactParam = vtxu::computeIP(interactionPoint, d->vertex(), d->momentum(), true);
              if(verbose) {
                cout << Form("Muon found in B decay: %.2f %.2f %.2f", d->pt(), d->eta(), d->phi()) << endl;
                auto vtx = p.vertex();
                cout << Form("B vertex: %.2f %.2f %.2f", vtx.x(), vtx.y(), vtx.z()) << endl;
                auto muVtx = d->vertex();
                cout << Form("Muon vertex: %.2f %.2f %.2f", muVtx.x(), muVtx.y(), muVtx.z()) << endl;
                cout << "IP: " << mu_impactParam << endl;
              }
            }
          }
        }
      }
    }
    (*indexBmc) = i_B;

    p4["B"] = TLorentzVector();
    p4["Ks"] = TLorentzVector();
    p4["Dst"] = TLorentzVector();
    p4["pis"] = TLorentzVector();
    p4["D0"] = TLorentzVector();

    if(i_B >= 0){
      auto p = (*PrunedGenParticlesHandle)[i_B];
      p4["B"].SetPtEtaPhiM(p.pt(), p.eta(), p.phi(), p.mass());

      for(auto d : p.daughterRefVector()) {
        if(d->pdgId() == 321) {
          p4["Ks"].SetPtEtaPhiM(d->pt(), d->eta(), d->phi(), d->mass());
        }
        else if(d->pdgId() == -413) {
          p4["Dst"].SetPtEtaPhiM(d->pt(), d->eta(), d->phi(), d->mass());

          for(auto dd : d->daughterRefVector()) {
            if(dd->pdgId() == -211) {
              p4["pis"].SetPtEtaPhiM(dd->pt(), dd->eta(), dd->phi(), dd->mass());
            }
            else if(abs(dd->pdgId()) == 421) {
              p4["D0"].SetPtEtaPhiM(dd->pt(), dd->eta(), dd->phi(), dd->mass());
              // for(auto ddd : dd->daughterRefVector()) {
              //   if (ddd->pdgId() == -211 || ddd->pdgId() == 321)
              //   {
              //     string name = ddd->pdgId() == -211 ? "pi" : "K";
              //   }
              // }
            }
          }
        }
      }
    }

    for(auto kv : p4) {
      AddTLVToOut(kv.second, "MC_"+kv.first, &(*outputNtuplizer));
      if(verbose) {
        cout << kv.first << Form(": %.2f %.2f %.2f", kv.second.Pt(), kv.second.Eta(), kv.second.Phi()) << endl;
      }
    }
    (*outputNtuplizer)["MC_mu_IP"] = mu_impactParam;



    iEvent.put(move(outputNtuplizer), "outputNtuplizer");
    iEvent.put(move(indexBmc), "indexBmc");
    return;
}


void MCTruthB2DstKProducer::AddTLVToOut(TLorentzVector v, string n, map<string, float>* outv) {
  (*outv)[n+"_pt"] = v.Pt();
  (*outv)[n+"_eta"] = v.Eta();
  (*outv)[n+"_phi"] = v.Phi();
  (*outv)[n+"_P"] = v.P();
  return;
}


DEFINE_FWK_MODULE(MCTruthB2DstKProducer);
