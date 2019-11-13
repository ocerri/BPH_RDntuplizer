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

bool auxIsAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle){
  if(ancestor->pt() == particle->pt() && ancestor->eta() == particle->eta() && ancestor->phi() == particle->phi() && ancestor->pdgId() == particle->pdgId()){
    return true;
  }
  for(size_t i=0;i< particle->numberOfMothers();i++)
  {
    if(auxIsAncestor(ancestor,particle->mother(i))) return true;
  }
  return false;
}

void auxPrintDau(const reco::Candidate* p, int kk=0) {
  auto Nd = p->numberOfDaughters();
  if(kk != 0 && Nd == 0) return;
  for(auto i = 0; i < kk; i++) {cout << "\t";}
  if(kk == 0) {cout << "\033[1;36;m" << p->pdgId() << "\033[0m -->";}
  else {cout << p->pdgId() << " -->";}
  for(uint i_d = 0; i_d < Nd; i_d++) {
    const auto d = p->daughter(i_d);
    if((d->pdgId() == -13 || d->pdgId() == -413) && kk >= 0) {
      cout << " \033[1;32m" << d->pdgId() << "\033[0m";
    }
    else if(abs(d->pdgId()) > 10 && abs(d->pdgId()) < 17) {
      cout << " \033[1;35m" << d->pdgId() << "\033[0m";
    }
    else {cout << " " << d->pdgId();}
  }
  cout << endl;
  if(kk < 0) return;
  for(uint i_d = 0; i_d < Nd; i_d++) auxPrintDau(p->daughter(i_d), kk+1);
}

class MCTruthB2DstMuProducer : public edm::EDProducer {

public:

    explicit MCTruthB2DstMuProducer(const edm::ParameterSet &iConfig);
    void AddTLVToOut(TLorentzVector, string, map<string, float>*);

    ~MCTruthB2DstMuProducer() override {};

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



MCTruthB2DstMuProducer::MCTruthB2DstMuProducer(const edm::ParameterSet &iConfig)
{
    verbose = iConfig.getParameter<int>( "verbose" );
    PrunedParticlesSrc_ = consumes<vector<reco::GenParticle>>(edm::InputTag("prunedGenParticles"));
    TrgMuonSrc_ = consumes<vector<pat::Muon>> ( iConfig.getParameter<edm::InputTag>( "trgMuons" ) );

    produces<map<string, float>>("outputNtuplizer");
    produces<map<string, vector<float>>>("outputVecNtuplizer");
    produces<int>("indexBmc");
}


void MCTruthB2DstMuProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    if (verbose) {cout << "-----------  MC Truth ----------\n";}

    // Get prunedGenParticles
    edm::Handle<std::vector<reco::GenParticle>> PrunedGenParticlesHandle;
    iEvent.getByToken(PrunedParticlesSrc_, PrunedGenParticlesHandle);
    unsigned int N_PrunedGenParticles = PrunedGenParticlesHandle->size();

    edm::Handle<vector<pat::Muon>> trgMuonsHandle;
    iEvent.getByToken(TrgMuonSrc_, trgMuonsHandle);
    auto trgMu = (*trgMuonsHandle)[0];
    TLorentzVector trgMu_TLV;
    trgMu_TLV.SetPtEtaPhiM(trgMu.pt(), trgMu.eta(), trgMu.phi(), mass_mu);


    // Output collection
    unique_ptr<map<string, float>> outputNtuplizer(new map<string, float>);
    unique_ptr<map<string, vector<float>>> outputVecNtuplizer(new map<string, vector<float>>);
    unique_ptr<int> indexBmc(new int);


    if(verbose) {
      uint nb = 0;
      for(unsigned int i = 0; i < N_PrunedGenParticles; i++) {
        auto p = (*PrunedGenParticlesHandle)[i];
        int pId = abs(p.pdgId());
        bool sel = pId == 511 || pId == 521;
        sel = sel || pId == 531 || pId == 541;
        sel = sel || pId == 10423;

        if(sel) {
          int kk = 0;
          if(p.numberOfDaughters() < 2) kk = -1;
          cout << "idx: " << i << endl;
          auxPrintDau(&p, kk);
          cout << endl;
          nb ++;
        }
      }
      cout << "Number of B found: " << nb << endl;
    }

    // Looking for the B -> Mu/Tau + D* + X with the cloasest muon to the trigger muon
    int i_B = -1;
    float mInv_min = 9999999999999999;
    for(unsigned int i = 0; i < N_PrunedGenParticles; i++) {
      auto p = (*PrunedGenParticlesHandle)[i];
      if (p.pdgId() != 511 || p.numberOfDaughters()<=1) continue;

      TLorentzVector p_Mu;
      bool mu_found = false;
      bool Dst_found = false;
      for(auto d : p.daughterRefVector()) {
        if(d->pdgId() == -13) {
          p_Mu.SetPtEtaPhiM(d->pt(), d->eta(), d->phi(), d->mass());
          mu_found = true;
        }
        else if(d->pdgId() == -15) {
          for (auto dd : d->daughterRefVector()) {
            if (dd->pdgId() == -13) {
              p_Mu.SetPtEtaPhiM(dd->pt(), dd->eta(), dd->phi(), dd->mass());
              mu_found = true;
            }
          }
        }
        else if(d->pdgId() == -413) Dst_found = true;
      }

      if(!mu_found || !Dst_found) continue;
      auto mInv = (p_Mu + trgMu_TLV).M();
      if(mInv < mInv_min) { mInv_min = mInv; i_B = i;}
    }

    // Look for the bkg process
    if(i_B == -1) {
      for(unsigned int i = 0; i < N_PrunedGenParticles; i++) {
        auto p = (*PrunedGenParticlesHandle)[i];
        if(p.numberOfDaughters()<=1) continue; // Needs to decay and not oscillate into a different one
        int pId = abs(p.pdgId());
        if(pId != 511 && pId != 521 && pId != 531 && pId != 541) continue; // Require a B or a B_s meson

        TLorentzVector p_Mu;
        bool mu_found = false;
        bool Dst_found = false;

        for(auto d : *PrunedGenParticlesHandle) {
          if(d.pdgId() == -13 && auxIsAncestor(&p, &d)) {
            p_Mu.SetPtEtaPhiM(d.pt(), d.eta(), d.phi(), d.mass());
            mu_found = true;
          }
          else if(d.pdgId() == -413 && auxIsAncestor(&p, &d)) Dst_found = true;
        }
        if(!mu_found || !Dst_found) continue;
        auto mInv = (p_Mu + trgMu_TLV).M();
        if(mInv < mInv_min) { mInv_min = mInv; i_B = i;}

        if(verbose) {
          cout << "Recognized B --> D*Mu +(X): " << endl;
          auxPrintDau(&p, -1);
          cout << endl << endl;
        }

      }
    }

    (*indexBmc) = i_B;

    float trgMuon_match_BMuon = 0;
    map<string, TLorentzVector> p4;
    p4["B"] = TLorentzVector();
    p4["mu"] = TLorentzVector();
    p4["Dst"] = TLorentzVector();
    p4["pis"] = TLorentzVector();
    p4["D0"] = TLorentzVector();
    p4["pi"] = TLorentzVector();
    p4["K"] = TLorentzVector();

    map<string, reco::Candidate::Point> vtx;
    vtx["B"] = reco::Candidate::Point();
    vtx["Dst"] = reco::Candidate::Point();
    vtx["D0"] = reco::Candidate::Point();
    vtx["K"] = reco::Candidate::Point();



    if(i_B >= 0){
      auto p = (*PrunedGenParticlesHandle)[i_B];
      vtx["B"] = p.vertex();
      p4["B"].SetPtEtaPhiM(p.pt(), p.eta(), p.phi(), p.mass());
      (*outputVecNtuplizer)["MC_decay"].push_back(p.pdgId());

      for(auto d : p.daughterRefVector()) {
        (*outputVecNtuplizer)["MC_decay"].push_back(d->pdgId());
      }

      for(auto d : *PrunedGenParticlesHandle) {
        if(d.pdgId() == -13 && auxIsAncestor(&p, &d)) {
          p4["mu"].SetPtEtaPhiM(d.pt(), d.eta(), d.phi(), d.mass());
        }
        else if(d.pdgId() == -15 && auxIsAncestor(&p, &d)) {
          for (auto dd : d.daughterRefVector()) {
            if (dd->pdgId() == -13) {
              p4["mu"].SetPtEtaPhiM(dd->pt(), dd->eta(), dd->phi(), dd->mass());
            }
          }
        }
        else if(d.pdgId() == -413 && auxIsAncestor(&p, &d)) {
          p4["Dst"].SetPtEtaPhiM(d.pt(), d.eta(), d.phi(), d.mass());
          vtx["Dst"] = d.vertex();
          for(auto dd : d.daughterRefVector()) {
            if(dd->pdgId() == -211) {
              p4["pis"].SetPtEtaPhiM(dd->pt(), dd->eta(), dd->phi(), dd->mass());
            }
            else if(abs(dd->pdgId()) == 421) {
              p4["D0"].SetPtEtaPhiM(dd->pt(), dd->eta(), dd->phi(), dd->mass());
              vtx["D0"] = dd->vertex();
              for(auto ddd : dd->daughterRefVector()) {
                if (ddd->pdgId() == -211 || ddd->pdgId() == 321)
                {
                  string name = ddd->pdgId() == -211 ? "pi" : "K";
                  p4[name].SetPtEtaPhiM(ddd->pt(), ddd->eta(), ddd->phi(), ddd->mass());
                  if (ddd->pdgId() == 321) vtx["K"] = ddd->vertex();
                }
              }
            }
          }
        }
      }

      double dPhi_mu = vtxu::dPhi(trgMu.phi(), p4["mu"].Phi());
      double dEta_mu = trgMu.eta() - p4["mu"].Eta();
      double dPt_rel_mu = (trgMu.pt() - p4["mu"].Pt())/trgMu.pt();
      if (p4["mu"].Pt() > 0 && fabs(dPhi_mu) < 0.01 && fabs(dEta_mu) < 0.01 && fabs(dPt_rel_mu) < 0.1 ) {
        trgMuon_match_BMuon = 1.;
      }
    }
    else {
      (*outputVecNtuplizer)["MC_decay"].push_back(-1);
    }

    (*outputNtuplizer)["trgMuon_match_BMuon"] = trgMuon_match_BMuon;

    for(auto kv : p4) {
      AddTLVToOut(kv.second, "MC_"+kv.first, &(*outputNtuplizer));
    }

    for(auto kv : vtx) {
      auto n = kv.first;
      auto v = kv.second;
      (*outputNtuplizer)["MC_prodVtx_"+n+"_dxy"] = hypot(v.x(), v.y());
      (*outputNtuplizer)["MC_prodVtx_"+n+"_phi"] = atan2(v.y(), v.x());
      (*outputNtuplizer)["MC_prodVtx_"+n+"_z"] = v.z();
    }

    if(verbose) {
      cout << i_B << ": ";
      for(auto v : (*outputVecNtuplizer)["MC_decay"]) {cout << v << " ";} cout << endl;
      cout << "Muon matched: " << trgMuon_match_BMuon << endl;
      cout << Form("Muon pt:%.1f eta:%.1f", p4["mu"].Pt(), p4["mu"].Eta()) << endl;
      if(p4["mu"].Pt() > 0) {
        cout << "Distance between MC muon and trigger muon" << endl;
        cout << "dPhi: " << vtxu::dPhi(trgMu.phi(), p4["mu"].Phi()) << endl;
        cout << "dEta: " << trgMu.eta() - p4["mu"].Eta() << endl;
        cout << "dPt: " << (trgMu.pt() - p4["mu"].Pt())/trgMu.pt() << endl;
      }
    }

    (*outputNtuplizer)["MC_M2_miss"] = (p4["B"] - p4["Dst"] - p4["mu"]).M2();
    (*outputNtuplizer)["MC_q2"] = (p4["B"] - p4["Dst"]).M2();

    TLorentzVector p4st_mu(p4["mu"]);
    p4st_mu.Boost(-1*p4["B"].BoostVector());
    (*outputNtuplizer)["MC_Est_mu"] = p4st_mu.E();


    iEvent.put(move(outputNtuplizer), "outputNtuplizer");
    iEvent.put(move(outputVecNtuplizer), "outputVecNtuplizer");
    iEvent.put(move(indexBmc), "indexBmc");
    return;
}


void MCTruthB2DstMuProducer::AddTLVToOut(TLorentzVector v, string n, map<string, float>* outv) {
  (*outv)[n+"_pt"] = v.Pt();
  (*outv)[n+"_eta"] = v.Eta();
  (*outv)[n+"_phi"] = v.Phi();
  (*outv)[n+"_P"] = v.P();
  return;
}


DEFINE_FWK_MODULE(MCTruthB2DstMuProducer);
