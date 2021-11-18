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

class MCTruthB2DstMuProducer_standalone : public edm::EDProducer {

public:

    explicit MCTruthB2DstMuProducer_standalone(const edm::ParameterSet &iConfig);
    void AddTLVToOut(TLorentzVector, string, map<string, float>*);

    ~MCTruthB2DstMuProducer_standalone() override {};

private:
    virtual void produce(edm::Event&, const edm::EventSetup&);

    bool auxIsAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle){
      if(ancestor->pt() == particle->pt() && ancestor->eta() == particle->eta() && ancestor->phi() == particle->phi() && ancestor->pdgId() == particle->pdgId()){
        return true;
      }
      for(size_t i=0;i< particle->numberOfMothers();i++)
      {
        if(auxIsAncestor(ancestor,particle->mother(i))) return true;
      }
      return false;
    };

    void auxPrintDau(const reco::Candidate* p, int kk=0) {
      auto Nd = p->numberOfDaughters();
      if(kk != 0 && Nd == 0) return;
      for(auto i = 0; i < kk; i++) {cout << "\t";}
      if(kk == 0) {cout << "\033[1;36;m" << p->pdgId() << "\033[0m -->";}
      else {cout << p->pdgId() << " -->";}
      for(uint i_d = 0; i_d < Nd; i_d++) {
        const auto d = p->daughter(i_d);
        if((abs(d->pdgId()) == 13 || abs(d->pdgId()) == 413)) {
          cout << " \033[1;32m" << d->pdgId() << "\033[0m";
        }
        else if(abs(d->pdgId()) > 10 && abs(d->pdgId()) < 17) {
          cout << " \033[" << abs(d->pdgId())%2 << ";35m" << d->pdgId() << "\033[0m";
        }
        else {cout << " " << d->pdgId();}
      }
      cout << endl;
      if(kk < 0) return;
      for(uint i_d = 0; i_d < Nd; i_d++) auxPrintDau(p->daughter(i_d), kk+1);
    };

    // ----------member data ---------------------------
    edm::EDGetTokenT<vector<reco::GenParticle>> GenParticlesSrc_;
    edm::EDGetTokenT<vector<pat::PackedGenParticle>> PackedParticlesSrc_;


    double mass_mu = 0.1056583745;
    double mass_pi = 0.13957018;
    double mass_K = 0.493677;
    double mass_Dst = 2.01026;
    double mass_B0 = 5.27961;

    int verbose = 0;
};



MCTruthB2DstMuProducer_standalone::MCTruthB2DstMuProducer_standalone(const edm::ParameterSet &iConfig)
{
    verbose = iConfig.getParameter<int>( "verbose" );
    GenParticlesSrc_ = consumes<vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>( "genParticlesCollection" ));
    PackedParticlesSrc_ = consumes<vector<pat::PackedGenParticle>>(edm::InputTag("packedGenParticles"));

    produces<map<string, float>>("outputNtuplizer");
    produces<map<string, vector<float>>>("outputVecNtuplizer");
}


void MCTruthB2DstMuProducer_standalone::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    if (verbose) {cout << "------------------  MC Truth -----------------\n";}

    // Get prunedGenParticles
    edm::Handle<std::vector<reco::GenParticle>> GenParticlesHandle;
    iEvent.getByToken(GenParticlesSrc_, GenParticlesHandle);
    unsigned int N_PrunedGenParticles = GenParticlesHandle->size();

    // Get packedGenParticles
    edm::Handle<std::vector<pat::PackedGenParticle>> PackedGenParticlesHandle;
    iEvent.getByToken(PackedParticlesSrc_, PackedGenParticlesHandle);
    unsigned int N_PackedGenParticles = PackedGenParticlesHandle->size();


    // Output collection
    unique_ptr<map<string, float>> outputNtuplizer(new map<string, float>);
    unique_ptr<map<string, vector<float>>> outputVecNtuplizer(new map<string, vector<float>>);


    if(verbose) {
      cout << "Bottom hadrons MC history" << endl;
      for(unsigned int i = 0; i < N_PrunedGenParticles; i++) {
        auto p = (*GenParticlesHandle)[i];
        if (p.numberOfDaughters() <= 1) continue;
        int pId = abs(p.pdgId());

        if ( (pId/100)%10 != 5 && (pId/1000)%10 != 5 ) continue;
        cout << "idx: " << i << ", mother: " << p.mother()->pdgId() <<endl;
        auxPrintDau(&p, -1);
        cout << endl;
      }
      cout << endl << endl;
    }

    // Looking for the B(s) -> (mu D*)^0 + X cloasest to one of the candidates
    reco::Candidate::Point interactionPoint(-9999999999, -999, -999);
    int i_B = -1;
    int j_mu = -1, j_Dst = -1;
    int nB2DstMuX = 0;
    for(uint i = 0; i < N_PrunedGenParticles; i++) {
      auto p = (*GenParticlesHandle)[i];
      if(p.isHardProcess() && interactionPoint.x() == -9999999999) {
        interactionPoint = p.vertex();
        if(verbose) {cout << "[Hard process] " << p.pdgId() << ": " << p.vx() << ", " << p.vy() << endl;}
      }


      if(p.numberOfDaughters()<=1) continue; // Needs to decay and not oscillate into a different one
      int pId = abs(p.pdgId());
      if(pId != 511 && pId != 521 && pId != 531) continue; // Require B or B_s meson

      vector<int> idx_mu;
      vector<int> idx_Dst;

      for(uint j = 0; j < N_PrunedGenParticles; j++) {
        auto d = (*GenParticlesHandle)[j];
        int dID = abs(d.pdgId());
        if(dID != 13 && dID != 413) continue; // If not mu or D*, skip
        if (!auxIsAncestor(&p, &d)) continue;

        if (dID == 13) idx_mu.push_back(j);
        else if (dID == 413) {
          bool pis_found = false;
          bool DuToKpi_found = false;
          for(auto dd : d.daughterRefVector()) {
            if(abs(dd->pdgId()) == 211) pis_found = true;
            else if(abs(dd->pdgId()) == 421) {

              bool pi_found = false;
              bool K_found = false;
              bool other_nonGamma = false;
              for(auto ddd : dd->daughterRefVector()) {
                if (abs(ddd->pdgId()) == 211) pi_found = true;
                else if ( abs(ddd->pdgId()) == 321 ) K_found = true;
                else if ( ddd->pdgId() != 22 ) other_nonGamma = true; // Allow only PHOTOS decays
              }

              if (pi_found && K_found && !other_nonGamma) DuToKpi_found = true;
            }
          }
          if (pis_found && DuToKpi_found) idx_Dst.push_back(j);
        }
      }

      vector<pair<int, int>> idxMuDstPairs;
      for(auto jMu : idx_mu) {
        auto mu = (*GenParticlesHandle)[jMu];
        for (auto jDst : idx_Dst) {
          auto Dst = (*GenParticlesHandle)[jDst];
          // mu+(-13)Dst-(-413) or mu-(13)Dst+(413)
          if (mu.pdgId()*Dst.pdgId() > 0)  idxMuDstPairs.push_back(make_pair(jMu, jDst));
        }
      }

      if (idxMuDstPairs.size() == 0) continue;
      nB2DstMuX++;
      if(verbose) {
        cout << "Recognized B(s) -> D* mu X in MC particle " << i << ": " << endl;
        auxPrintDau(&p, 0);
        cout << endl << endl;
      }
      i_B = i;
      j_mu = idxMuDstPairs[0].first;
      j_Dst = idxMuDstPairs[0].second;

    }


    vector<uint> bOgIdx;
    int bSelAncestorIdx = -1;
    for(unsigned int i = 0; i < N_PrunedGenParticles && i_B >=0; i++) {
      auto p = (*GenParticlesHandle)[i];
      int pId = abs(p.pdgId());

      if ( (pId/100)%10 != 5 && (pId/1000)%10 != 5 ) continue;

      bool bMother = false;
      for (uint j=0; j < p.numberOfMothers(); j++) {
        int mID = abs(p.mother(j)->pdgId());
        if ( (mID/100)%10 == 5 || (mID/1000)%10 == 5 ) bMother = true;
      }
      if (bMother) continue;
      if (verbose) {
        cout << "idx: " << i << ", mother: " << p.mother()->pdgId() <<endl;
        auxPrintDau(&p, -1);
      }
      if (auxIsAncestor(&p, &((*GenParticlesHandle)[i_B])   )) {
        bSelAncestorIdx = i;
        if (verbose) {cout << "Anchestor of selected B" << endl;}
      }
      else bOgIdx.push_back(i);
      if (verbose) {cout << endl;}
    }

    float BB_dR = 1e9;
    float BB_dphi = 1e9;
    float BB_mass = -1;
    if (bSelAncestorIdx != -1 && bOgIdx.size()>0 ) {
      auto bAncestor = (*GenParticlesHandle)[bSelAncestorIdx];
      TLorentzVector lAnc;
      lAnc.SetPtEtaPhiM(bAncestor.pt(), bAncestor.eta(), bAncestor.phi(), bAncestor.mass());
      for (auto idx : bOgIdx) {
        auto p = (*GenParticlesHandle)[idx];
        auto dR = vtxu::dR(p.phi(), bAncestor.phi(), p.eta(), bAncestor.eta());
        if (dR < BB_dR) {
          if (verbose) {cout << "Chosen " << idx << endl;}
          BB_dR = dR;
          BB_dphi = vtxu::dPhi(p.phi(), bAncestor.phi());
          TLorentzVector lp;
          lp.SetPtEtaPhiM(p.pt(), p.eta(), p.phi(), p.mass());
          BB_mass = (lp + lAnc).M();
        }
      }
    }
    // assert(bOgIdx.size() == 2);

    (*outputNtuplizer)["MC_nB2DstMuX"] = nB2DstMuX;

    (*outputNtuplizer)["MC_nAddOgB"] = bOgIdx.size();
    (*outputNtuplizer)["MC_bestBB_dR"] = BB_dR;
    (*outputNtuplizer)["MC_bestBB_dphi"] = BB_dphi;
    (*outputNtuplizer)["MC_bestBB_mass"] = BB_mass;

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
    vtx["mu"] = reco::Candidate::Point();
    vtx["Dst"] = reco::Candidate::Point();
    vtx["D0"] = reco::Candidate::Point();
    vtx["pi"] = reco::Candidate::Point();

    (*outputVecNtuplizer)["MC_decay"] = {};
    (*outputNtuplizer)["MC_muMotherPdgId"] = 0;
    (*outputVecNtuplizer)["MC_muSistersPdgId"] = {};
    (*outputNtuplizer)["MC_DstMotherPdgId"] = 0;
    (*outputVecNtuplizer)["MC_DstSistersPdgId"] = {};
    (*outputNtuplizer)["MC_CharmedDstSisPdgId"] = 0;
    (*outputNtuplizer)["MC_StrangeDstSisPdgId"] = 0;
    (*outputNtuplizer)["MC_B_ctau"] = -1;
    (*outputNtuplizer)["MC_MassCharmedBDaugther"] = -1;

    (*outputNtuplizer)["MC_mu_TransvIP_PV"] = 0;
    (*outputNtuplizer)["MC_mu_TransvIP_vtxDst"] = 0;
    (*outputNtuplizer)["MC_mu_IP_vtxDst"] = 0;

    if(i_B >= 0){
      auto p = (*GenParticlesHandle)[i_B];

      // auto bPart = p.mother();
      // cout << p.pdgId() << " <-- " << bPart->pdgId() << flush;
      // while (bPart->numberOfMothers() > 0) {
      //   int mId = abs(bPart->mother()->pdgId());
      //   if (mId == 5 || (mId/ 100)%10 == 5 || (mId/1000)%10 == 5) {
      //     bPart = bPart->mother();
      //     cout <<  " <-- " << bPart->pdgId() << flush;
      //   }
      //   else break;
      // }
      // cout << endl;
      //
      // auto oPart = bPart->numberOfMothers()>0 ? bPart->mother() : bPart;
      // cout << oPart->pdgId() << " -->" << flush;
      // for(uint nDau=0; nDau < oPart->numberOfDaughters(); nDau++) {
      //   cout << " " << oPart->daughter(nDau)->pdgId() << flush;
      // }
      // cout << endl;


      vtx["B"] = p.vertex();
      p4["B"].SetPtEtaPhiM(p.pt(), p.eta(), p.phi(), p.mass());
      (*outputVecNtuplizer)["MC_decay"].push_back(p.mother()->pdgId());
      (*outputVecNtuplizer)["MC_decay"].push_back(0);
      (*outputVecNtuplizer)["MC_decay"].push_back(p.pdgId());
      (*outputVecNtuplizer)["MC_decay"].push_back(0);
      (*outputNtuplizer)["MC_B_ctau"] = vtxu::computeCTau(p);

      for(auto d : p.daughterRefVector()) {
        (*outputVecNtuplizer)["MC_decay"].push_back(d->pdgId());
        int pdgId = (int) abs(d->pdgId());
        if ((pdgId/100)%10 == 4) {
          (*outputNtuplizer)["MC_MassCharmedBDaugther"] = d->mass();
          if (verbose) {
            cout << d->pdgId() << Form(" mass: %.4f GeV", d->mass()) << endl;
          }
        }
      }

      auto mu = (*GenParticlesHandle)[j_mu];
      p4["mu"].SetPtEtaPhiM(mu.pt(), mu.eta(), mu.phi(), mu.mass());
      vtx["mu"] = mu.vertex();
      (*outputNtuplizer)["MC_mu_TransvIP_PV"] = vtxu::computeIP(interactionPoint, mu.vertex(), mu.momentum(), true);
      auto muMother = mu.mother();
      (*outputNtuplizer)["MC_muMotherPdgId"] = muMother->pdgId();
      for(uint iss = 0; iss <  muMother->numberOfDaughters(); iss++) {
        auto muSis = muMother->daughter(iss);
        if ( muSis->pdgId() != mu.pdgId() ) (*outputVecNtuplizer)["MC_muSistersPdgId"].push_back(muSis->pdgId());
      }

      auto Dst = (*GenParticlesHandle)[j_Dst];
      p4["Dst"].SetPtEtaPhiM(Dst.pt(), Dst.eta(), Dst.phi(), Dst.mass());
      vtx["Dst"] = Dst.vertex();
      auto DstMother = Dst.mother();
      (*outputNtuplizer)["MC_DstMotherPdgId"] = DstMother->pdgId();
      for(uint iss = 0; iss <  DstMother->numberOfDaughters(); iss++) {
        auto DstSis = DstMother->daughter(iss);
        int pdgId = DstSis->pdgId();
        if ( pdgId == Dst.pdgId() ) continue;

        (*outputVecNtuplizer)["MC_DstSistersPdgId"].push_back(DstSis->pdgId());
        int heavyQuark = (abs(pdgId)/100) %10;
        if ( heavyQuark == 4 ) (*outputNtuplizer)["MC_CharmedDstSisPdgId"] = pdgId;
        else if ( heavyQuark == 3 ) (*outputNtuplizer)["MC_StrangeDstSisPdgId"] = pdgId;
      }

      for(auto dd : Dst.daughterRefVector()) {
        if(abs(dd->pdgId()) == 211) {
          p4["pis"].SetPtEtaPhiM(dd->pt(), dd->eta(), dd->phi(), dd->mass());
        }
        else if(abs(dd->pdgId()) == 421) {
          p4["D0"].SetPtEtaPhiM(dd->pt(), dd->eta(), dd->phi(), dd->mass());
          vtx["D0"] = dd->vertex();
          for(auto ddd : dd->daughterRefVector()) {
            if (abs(ddd->pdgId()) == 211) {
              p4["pi"].SetPtEtaPhiM(ddd->pt(), ddd->eta(), ddd->phi(), ddd->mass());
              vtx["pi"] = ddd->vertex();
            }
            else if (abs(ddd->pdgId()) == 321) p4["K"].SetPtEtaPhiM(ddd->pt(), ddd->eta(), ddd->phi(), ddd->mass());
          }
        }
      }

      (*outputNtuplizer)["MC_mu_TransvIP_vtxDst"] = vtxu::computeDCA_linApprox(vtx["Dst"], vtx["mu"], mu.momentum(), true);
      (*outputNtuplizer)["MC_mu_IP_vtxDst"] = vtxu::computeDCA_linApprox(vtx["Dst"], vtx["mu"], mu.momentum(), false);


      if (verbose) {cout << "Looking for additional particle from the selected B -> D*Mu+X decay"<< endl;}
      for(int j = 0; ((uint)j) < N_PackedGenParticles; j++) {
        // if (j == j_Dst || j == j_mu || j == i_B) continue;
        auto pp = (*PackedGenParticlesHandle)[j];
        if (!auxIsAncestor(&p, &pp)) continue;
        if (auxIsAncestor(&Dst, &pp)) continue;
        if (fabs(mu.pt() - pp.pt())/mu.pt() < 0.001 && fabs(mu.eta() - pp.eta()) < 0.01) continue;
        if (abs(pp.pdgId()) == 12 || abs(pp.pdgId()) == 14 || abs(pp.pdgId()) == 14) continue;

        if (verbose) {cout << j << ": " << pp.pdgId() << endl;}
      }
      if (verbose) {cout << "--------------\n\n" << endl;}
    }


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
      cout << Form("MC muon pt=%.1f eta=%.2f IP=%1.2e", p4["mu"].Pt(), p4["mu"].Eta(), (*outputNtuplizer)["MC_mu_TransvIP_PV"]) << endl;
    }

    (*outputNtuplizer)["MC_M_vis"] = (p4["Dst"] + p4["mu"]).M();
    (*outputNtuplizer)["MC_M2_miss"] = (p4["B"] - p4["Dst"] - p4["mu"]).M2();
    (*outputNtuplizer)["MC_q2"] = (p4["B"] - p4["Dst"]).M2();

    TLorentzVector p4st_mu(p4["mu"]);
    p4st_mu.Boost(-1*p4["B"].BoostVector());
    (*outputNtuplizer)["MC_Est_mu"] = p4st_mu.E();


    iEvent.put(move(outputNtuplizer), "outputNtuplizer");
    iEvent.put(move(outputVecNtuplizer), "outputVecNtuplizer");
    return;
}


void MCTruthB2DstMuProducer_standalone::AddTLVToOut(TLorentzVector v, string n, map<string, float>* outv) {
  (*outv)[n+"_pt"] = v.Pt();
  (*outv)[n+"_eta"] = v.Eta();
  (*outv)[n+"_phi"] = v.Phi();
  (*outv)[n+"_P"] = v.P();
  return;
}


DEFINE_FWK_MODULE(MCTruthB2DstMuProducer_standalone);
