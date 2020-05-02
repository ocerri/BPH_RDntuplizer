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
    edm::EDGetTokenT<map<string, vector<float>>> decayTreeOutSrc_;

    double mass_mu = 0.1056583745;
    double mass_pi = 0.13957018;
    double mass_K = 0.493677;
    double mass_Dst = 2.01026;
    double mass_B0 = 5.27961;

    int verbose = 0;
};



MCTruthB2DstMuProducer::MCTruthB2DstMuProducer(const edm::ParameterSet &iConfig)
{
    verbose = iConfig.getParameter<int>( "verbose" );
    PrunedParticlesSrc_ = consumes<vector<reco::GenParticle>>(edm::InputTag("prunedGenParticles"));
    decayTreeOutSrc_ = consumes<map<string, vector<float>>>(iConfig.getParameter<edm::InputTag>( "decayTreeVecOut" ));

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

    edm::Handle<map<string, vector<float>>> outMapHandle;
    iEvent.getByToken(decayTreeOutSrc_, outMapHandle);
    vector<float> muPt, muEta, muPhi;
    vector<float> DstPt, DstEta, DstPhi;
    for( auto const& kv : (*outMapHandle) ) {
      if (kv.first == "Dst_refitD0pismu_pt") DstPt = kv.second;
      if (kv.first == "Dst_refitD0pismu_eta") DstEta = kv.second;
      if (kv.first == "Dst_refitD0pismu_phi") DstPhi = kv.second;
      if (kv.first == "mu_refitD0pismu_pt") muPt = kv.second;
      if (kv.first == "mu_refitD0pismu_eta") muEta = kv.second;
      if (kv.first == "mu_refitD0pismu_phi") muPhi = kv.second;
    }
    vector<TLorentzVector> reco_muTLVs, reco_DstTLVs;
    for(uint i = 0; i < muPt.size(); i ++) {
      TLorentzVector m, dst;
      dst.SetPtEtaPhiM(DstPt[i], DstEta[i], DstPhi[i], mass_Dst);
      reco_DstTLVs.push_back(dst);
      m.SetPtEtaPhiM(muPt[i], muEta[i], muPhi[i], mass_mu);
      reco_muTLVs.push_back(m);
    }


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
    int i_cand = -1;
    float best_mInvMu = 1e9, best_mInvDst = 1e9;
    TLorentzVector bestMu, bestDst;
    int nB02DstMuX = 0;

    reco::Candidate::Point interactionPoint(-9999999999, -999, -999);
    for(uint i = 0; i < N_PrunedGenParticles; i++) {
      auto p = (*PrunedGenParticlesHandle)[i];
      if(p.isHardProcess() && interactionPoint.x() == -9999999999) {
        interactionPoint = p.vertex();
        if(verbose) {cout << "[Hard process] " << p.pdgId() << ": " << p.vx() << ", " << p.vy() << endl;}
      }

      if (p.pdgId() != 511 || p.numberOfDaughters()<=1) continue;

      TLorentzVector pMu, pDst;
      bool mu_found = false;
      bool Dst_found = false;
      for(auto d : p.daughterRefVector()) {
        if(d->pdgId() == -13) {
          pMu.SetPtEtaPhiM(d->pt(), d->eta(), d->phi(), d->mass());
          mu_found = true;
        }
        else if(d->pdgId() == -15) {
          for (auto dd : d->daughterRefVector()) {
            if (dd->pdgId() == -13) {
              pMu.SetPtEtaPhiM(dd->pt(), dd->eta(), dd->phi(), dd->mass());
              mu_found = true;
            }
          }
        }
        else if(d->pdgId() == -413) {
          pDst.SetPtEtaPhiM(d->pt(), d->eta(), d->phi(), d->mass());
          Dst_found = true;
        }
      }

      if(!mu_found || !Dst_found) continue;
      nB02DstMuX++;
      if(verbose) {
        cout << "Recognized B0 --> D*MuNu: " << endl;
        auxPrintDau(&p, -1);
        cout << endl << endl;
      }

      for(uint j = 0; j < reco_muTLVs.size(); j++) {
        double mMu = (pMu + reco_muTLVs[j]).M();
        double mDst = (pDst + reco_DstTLVs[j]).M();
        if(mMu < best_mInvMu && mDst < best_mInvDst) {
          best_mInvMu = mMu;
          best_mInvDst = mDst;
          bestMu = reco_muTLVs[j];
          bestDst = reco_DstTLVs[j];
          i_B = i;
          i_cand = j;
        }
      }
    }

    // Look for the bkg process
    if(i_B == -1) {
      for(unsigned int i = 0; i < N_PrunedGenParticles; i++) {
        auto p = (*PrunedGenParticlesHandle)[i];
        if(p.numberOfDaughters()<=1) continue; // Needs to decay and not oscillate into a different one
        int pId = abs(p.pdgId());
        if(pId != 511 && pId != 521 && pId != 531 && pId != 541) continue; // Require a B or a B_s meson

        TLorentzVector pMu, pDst;
        bool mu_found = false;
        bool Dst_found = false;

        for(auto d : *PrunedGenParticlesHandle) {
          if(d.pdgId() == -13 && auxIsAncestor(&p, &d)) {
            pMu.SetPtEtaPhiM(d.pt(), d.eta(), d.phi(), d.mass());
            mu_found = true;
          }
          else if(d.pdgId() == -413 && auxIsAncestor(&p, &d)) {
            pDst.SetPtEtaPhiM(d.pt(), d.eta(), d.phi(), d.mass());
            Dst_found = true;
          }
        }

        if(!mu_found || !Dst_found) continue;
        nB02DstMuX++;
        if(verbose) {
          cout << "Recognized B --> D*Mu +(X): " << endl;
          auxPrintDau(&p, -1);
          cout << endl << endl;
        }

        for(uint j = 0; j < reco_muTLVs.size(); j++) {
          double mMu = (pMu + reco_muTLVs[j]).M();
          double mDst = (pDst + reco_DstTLVs[j]).M();
          if(mMu < best_mInvMu && mDst < best_mInvDst) {
            best_mInvMu = mMu;
            best_mInvDst = mDst;
            bestMu = reco_muTLVs[j];
            bestDst = reco_DstTLVs[j];
            i_B = i;
            i_cand = j;
          }
        }
      }
    }

    (*outputNtuplizer)["MC_nB02DstMuX"] = nB02DstMuX;
    (*outputNtuplizer)["MC_idxCand"] = i_cand;
    (*indexBmc) = i_B;

    // assert(i_B >= 0);

    float recoMuon_match_BMuon = 0;
    map<string, TLorentzVector> p4;
    p4["B"] = TLorentzVector();
    p4["mu"] = TLorentzVector();
    p4["Dst"] = TLorentzVector();
    p4["pis"] = TLorentzVector();
    p4["D0"] = TLorentzVector();
    p4["pi"] = TLorentzVector();
    p4["K"] = TLorentzVector();
    double mu_impactParam = -1;

    map<string, reco::Candidate::Point> vtx;
    vtx["B"] = reco::Candidate::Point();
    vtx["Dst"] = reco::Candidate::Point();
    vtx["D0"] = reco::Candidate::Point();
    vtx["K"] = reco::Candidate::Point();

    // (*outputNtuplizer)["MC_muMotherPdgId"] = 0;
    // (*outputNtuplizer)["MC_muMotherNdaughters"] = 0;
    (*outputNtuplizer)["MC_munuSisterPdgId"] = 0;
    (*outputNtuplizer)["MC_DstSisPdgId_light"] = 0;
    (*outputNtuplizer)["MC_DstSisPdgId_heavy"] = 0;


    if(i_B >= 0){
      auto p = (*PrunedGenParticlesHandle)[i_B];
      vtx["B"] = p.vertex();
      p4["B"].SetPtEtaPhiM(p.pt(), p.eta(), p.phi(), p.mass());
      (*outputVecNtuplizer)["MC_decay"].push_back(p.mother()->pdgId());
      (*outputVecNtuplizer)["MC_decay"].push_back(p.pdgId());

      for(auto d : p.daughterRefVector()) {
        (*outputVecNtuplizer)["MC_decay"].push_back(d->pdgId());
      }

      for(auto d : *PrunedGenParticlesHandle) {
        if(d.pdgId() == -13 && auxIsAncestor(&p, &d)) {
          p4["mu"].SetPtEtaPhiM(d.pt(), d.eta(), d.phi(), d.mass());
          mu_impactParam = vtxu::computeIP(interactionPoint, d.vertex(), d.momentum(), true);
          auto muMother = d.mother();
          (*outputNtuplizer)["MC_muMotherPdgId"] = muMother->pdgId();
          (*outputNtuplizer)["MC_muMotherNdaughters"] = muMother->numberOfDaughters();
          for(uint iss = 0; iss <  muMother->numberOfDaughters(); iss++) {
            auto muSis = muMother->daughter(iss);
            if (muSis->pdgId() != 22 && muSis->pdgId() != 14 && muSis->pdgId() != -13) (*outputNtuplizer)["MC_munuSisterPdgId"] = muSis->pdgId();
          }
        }
        else if(d.pdgId() == -413 && auxIsAncestor(&p, &d)) {
          p4["Dst"].SetPtEtaPhiM(d.pt(), d.eta(), d.phi(), d.mass());
          vtx["Dst"] = d.vertex();
          auto DstMother = d.mother();
          (*outputNtuplizer)["MC_DstMotherPdgId"] = DstMother->pdgId();
          (*outputNtuplizer)["MC_DstMotherNdaughters"] = DstMother->numberOfDaughters();
          for(uint iss = 0; iss <  DstMother->numberOfDaughters(); iss++) {
            auto DstSis = DstMother->daughter(iss);
            auto pdgId = DstSis->pdgId();
            if (pdgId == 22 || pdgId == 14 || pdgId == -413) continue;
            if (abs(pdgId) > 400) (*outputNtuplizer)["MC_DstSisPdgId_heavy"] = pdgId;
            else (*outputNtuplizer)["MC_DstSisPdgId_light"] = pdgId;
          }
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

      double dPhi_mu = vtxu::dPhi(bestMu.Phi(), p4["mu"].Phi());
      double dEta_mu = bestMu.Eta() - p4["mu"].Eta();
      double dPt_rel_mu = (bestMu.Pt() - p4["mu"].Pt())/bestMu.Pt();
      if (p4["mu"].Pt() > 0 && fabs(dPhi_mu) < 0.01 && fabs(dEta_mu) < 0.01 && fabs(dPt_rel_mu) < 0.1 ) {
        recoMuon_match_BMuon = 1.;
      }
    }
    else {
      (*outputVecNtuplizer)["MC_decay"].push_back(-1);
    }

    (*outputNtuplizer)["recoMuon_match_BMuon"] = recoMuon_match_BMuon;

    for(auto kv : p4) {
      AddTLVToOut(kv.second, "MC_"+kv.first, &(*outputNtuplizer));
    }
    (*outputNtuplizer)["MC_mu_IP"] = mu_impactParam;

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
      cout << "Muon matched: " << recoMuon_match_BMuon << endl;
      cout << Form("Muon pt:%.1f eta:%.1f IP:%1.2e", p4["mu"].Pt(), p4["mu"].Eta(), mu_impactParam) << endl;
      if(p4["mu"].Pt() > 0) {
        cout << "Distance between MC muon and trigger muon" << endl;
        cout << "dPhi: " << vtxu::dPhi(bestMu.Phi(), p4["mu"].Phi()) << endl;
        cout << "dEta: " << bestMu.Eta() - p4["mu"].Eta() << endl;
        cout << "dPt: " << (bestMu.Pt() - p4["mu"].Pt())/bestMu.Pt() << endl;
      }
    }

    (*outputNtuplizer)["MC_M_vis"] = (p4["Dst"] + p4["mu"]).M();
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
