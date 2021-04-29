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
    edm::EDGetTokenT<vector<pat::PackedCandidate>> PFCandSrc_;

    edm::EDGetTokenT<map<string, vector<float>>> decayTreeOutSrc_;
    edm::EDGetTokenT<vector<pat::Muon>> trgMuonSrc_;
    int verbose = 0;
};



MCTruthB2JpsiKstProducer::MCTruthB2JpsiKstProducer(const edm::ParameterSet &iConfig)
{
    verbose = iConfig.getParameter<int>( "verbose" );
    PrunedParticlesSrc_ = consumes<vector<reco::GenParticle>>(edm::InputTag("prunedGenParticles"));
    PackedParticlesSrc_ = consumes<vector<pat::PackedGenParticle>>(edm::InputTag("packedGenParticles"));
    PFCandSrc_ = consumes<vector<pat::PackedCandidate>>(edm::InputTag("packedPFCandidates"));
    decayTreeOutSrc_ = consumes<map<string, vector<float>>>(iConfig.getParameter<edm::InputTag>( "decayTreeVecOut" ));
    trgMuonSrc_ = consumes<vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>( "triggerMuons" ));

    produces<int>("indexBmc");
    produces<map<string, float>>("outputNtuplizer");
    produces<map<string, vector<float>>>("outputVecNtuplizer");
}


void MCTruthB2JpsiKstProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    if (verbose) {cout << "-----------  MC Truth ----------\n";}

    // Get prunedGenParticles
    edm::Handle<vector<reco::GenParticle>> PrunedGenParticlesHandle;
    iEvent.getByToken(PrunedParticlesSrc_, PrunedGenParticlesHandle);
    unsigned int N_PrunedGenParticles = PrunedGenParticlesHandle->size();

    // Get packedGenParticles
    edm::Handle<vector<pat::PackedGenParticle>> PackedGenParticlesHandle;
    iEvent.getByToken(PackedParticlesSrc_, PackedGenParticlesHandle);

    // Get PF cand
    edm::Handle<vector<pat::PackedCandidate>> pfCandHandle;
    iEvent.getByToken(PFCandSrc_, pfCandHandle);

    // Get output from decay tree producer
    edm::Handle<map<string, vector<float>>> outMapHandle;
    iEvent.getByToken(decayTreeOutSrc_, outMapHandle);
    vector<float> reco_B_pt, reco_B_eta, reco_B_phi;
    for( auto const& kv : (*outMapHandle) ) {
      if (kv.first == "B_mumupiK_pt") reco_B_pt = kv.second;
      if (kv.first == "B_mumupiK_eta") reco_B_eta = kv.second;
      if (kv.first == "B_mumupiK_phi") reco_B_phi = kv.second;
    }
    if(verbose) {
      cout << "Reco B candidates" << endl;
      for(uint j = 0; j < reco_B_pt.size(); j++) {
        cout << Form("%d: pt:%.2f eta:%.2f phi:%.2f", j, reco_B_pt[j], reco_B_eta[j], reco_B_phi[j]) << endl;
      }
      cout << endl;
    }

    // Get trigger muons
    edm::Handle<vector<pat::Muon>> trgMuHandle;
    iEvent.getByToken(trgMuonSrc_, trgMuHandle);

    // Output collection
    unique_ptr<int> indexBmc(new int);
    unique_ptr<map<string, float>> outputNtuplizer(new map<string, float>);
    unique_ptr<map<string, vector<float>>> outputVecNtuplizer(new map<string, vector<float>>);

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
    vector<reco::GenParticle> genMuons;  //Used later to match trigger muons
    map<string, TLorentzVector> p4;
    int i_B = -1;
    int i_cand = -1;
    int nB02JpsiKst = 0;
    float best_dR = 1e9, best_dPt = 1e9;
    reco::Candidate::Point interactionPoint(-9999999999, -999, -999);
    for(unsigned int i = 0; i < N_PrunedGenParticles; i++) {
      auto p = (*PrunedGenParticlesHandle)[i];
      if(p.isHardProcess() && interactionPoint.x() == -9999999999) {
        interactionPoint = p.vertex();
        if(verbose) {cout << "[Hard process] " << p.pdgId() << ": " << p.vx() << ", " << p.vy() << endl;}
      }

      if (abs(p.pdgId()) == 511 && p.numberOfDaughters()>1) {
        bool Jpsi_found = false;
        bool Kst_found = false;
        for(auto d : p.daughterRefVector()) {
          if(abs(d->pdgId()) == 443) Jpsi_found = true;
          else if(abs(d->pdgId()) == 313) Kst_found = true;
        }
        if(Jpsi_found && Kst_found) nB02JpsiKst++;
        else continue;

        if(verbose) {
          cout << Form("%d (%sB0->JpsiKst): pt:%.2f eta:%.2f phi:%.2f", i, p.pdgId()>0? "" : "anti-", p.pt(), p.eta(), p.phi()) << endl;
        }

        for(uint j = 0; j < reco_B_pt.size(); j++) {
          float dR = vtxu::dR(p.phi(), reco_B_phi[j], p.eta(), reco_B_eta[j]);
          float dPt = fabs(reco_B_pt[j] - p.pt())/p.pt();
          if(verbose) {cout << Form("%.2e %.2e", dR, dPt) << endl;}
          if(dR  < best_dR && dPt < best_dPt) {
            best_dR = dR;
            best_dPt = dPt;
            i_B = i;
            i_cand = j;
            if(verbose) {cout << "selected" << endl;}
          }
        }
      }
      else if (abs(p.pdgId()) ==  13) genMuons.push_back(p);
    }

    (*indexBmc) = i_B;
    (*outputNtuplizer)["MC_nBd_to_JpsiKst"] = nB02JpsiKst;
    (*outputNtuplizer)["MC_idxCand"] = i_cand;

    (*outputNtuplizer)["MC_d_vtxB"] = -1;
    (*outputNtuplizer)["MC_dxy_vtxB"] = -1;

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

      // Charge conjugation sign
      int ccSign = p.pdgId() > 0? 1: -1;

      bool displSet = false;
      for(auto d : p.daughterRefVector()) {
        if (!displSet) {
          auto decayVtx = d->vertex();

          auto dxy = hypot(decayVtx.x() - interactionPoint.x(), decayVtx.y() - interactionPoint.y());
          auto d = hypot(decayVtx.z() - interactionPoint.z(), dxy);
          (*outputNtuplizer)["MC_d_vtxB"] = d;
          (*outputNtuplizer)["MC_dxy_vtxB"] = dxy;
          displSet = true;
        }

        if(d->pdgId()*ccSign == 313) {
          p4["Kst"].SetPtEtaPhiM(d->pt(), d->eta(), d->phi(), d->mass());
          for (auto dd : (*PackedGenParticlesHandle)) {
            if (dd.numberOfMothers() == 0) continue;
            auto m = dd.mother(0);
            if(m->pdgId()*ccSign == 313) {
              bool match = d->pt()==m->pt() && d->eta()==m->eta() && d->phi()==m->phi();
              if (!match) continue;
              if(dd.pdgId()*ccSign == -211) {
                p4["pi"].SetPtEtaPhiM(dd.pt(), dd.eta(), dd.phi(), dd.mass());
              }
              else if(dd.pdgId()*ccSign == 321) {
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

    // Retrieve the cloaset MC mu to the trigger muons
    for(pat::Muon recoMu : (*trgMuHandle)) {
      float best_dR = 1e9, best_dPt = 1e9;
      bool assigned = false;
      reco::GenParticle bestMu;
      for(reco::GenParticle mcMu : genMuons) {
        float dR = vtxu::dR(recoMu.phi(), mcMu.phi(), recoMu.eta(), mcMu.eta());
        float dPt = fabs(recoMu.pt() - mcMu.pt())/mcMu.pt();
        if(dR < best_dR && dPt < best_dPt) {
          best_dR = dR;
          best_dPt = dPt;
          bestMu = mcMu;
          assigned = true;
        }
      }
      if(assigned) {
        (*outputVecNtuplizer)["MC_trgMu_pt"].push_back(bestMu.pt());
        (*outputVecNtuplizer)["MC_trgMu_eta"].push_back(bestMu.eta());
        (*outputVecNtuplizer)["MC_trgMu_phi"].push_back(bestMu.phi());
        float impactParam = vtxu::computeIP(interactionPoint, bestMu.vertex(), bestMu.momentum(), true);
        (*outputVecNtuplizer)["MC_trgMu_IP"].push_back(impactParam);
      }
      else {
        (*outputVecNtuplizer)["MC_trgMu_pt"].push_back(-1);
        (*outputVecNtuplizer)["MC_trgMu_eta"].push_back(-999);
        (*outputVecNtuplizer)["MC_trgMu_phi"].push_back(-999);
        (*outputVecNtuplizer)["MC_trgMu_IP"].push_back(-999);
      }
    }

    iEvent.put(move(indexBmc), "indexBmc");
    iEvent.put(move(outputNtuplizer), "outputNtuplizer");
    iEvent.put(move(outputVecNtuplizer), "outputVecNtuplizer");
    return;
}


void MCTruthB2JpsiKstProducer::AddTLVToOut(TLorentzVector v, string n, map<string, float>* outv) {
  (*outv)[n+"_pt"] = v.Pt();
  (*outv)[n+"_eta"] = v.Eta();
  (*outv)[n+"_phi"] = v.Phi();
  return;
}


DEFINE_FWK_MODULE(MCTruthB2JpsiKstProducer);
