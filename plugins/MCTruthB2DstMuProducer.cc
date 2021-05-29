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
    edm::EDGetTokenT<vector<pat::PackedGenParticle>> PackedParticlesSrc_;
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
    PackedParticlesSrc_ = consumes<vector<pat::PackedGenParticle>>(edm::InputTag("packedGenParticles"));
    decayTreeOutSrc_ = consumes<map<string, vector<float>>>(iConfig.getParameter<edm::InputTag>( "decayTreeVecOut" ));

    produces<map<string, float>>("outputNtuplizer");
    produces<map<string, vector<float>>>("outputVecNtuplizer");
    produces<int>("indexBmc");
}


void MCTruthB2DstMuProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    if (verbose) {cout << "------------------  MC Truth -----------------\n";}

    // Get prunedGenParticles
    edm::Handle<std::vector<reco::GenParticle>> PrunedGenParticlesHandle;
    iEvent.getByToken(PrunedParticlesSrc_, PrunedGenParticlesHandle);
    unsigned int N_PrunedGenParticles = PrunedGenParticlesHandle->size();

    edm::Handle<map<string, vector<float>>> outMapHandle;
    iEvent.getByToken(decayTreeOutSrc_, outMapHandle);
    vector<float> muCharge;
    vector<float> muPt, muEta, muPhi;
    vector<float> DstPt, DstEta, DstPhi;
    vector<float> AddTkCharge, AddTkPt, AddTkEta, AddTkPhi;
    for( auto const& kv : (*outMapHandle) ) {
      if (kv.first == "mu_charge") muCharge = kv.second;
      else if (kv.first == "Dst_refitD0pismu_pt") DstPt = kv.second;
      else if (kv.first == "Dst_refitD0pismu_eta") DstEta = kv.second;
      else if (kv.first == "Dst_refitD0pismu_phi") DstPhi = kv.second;
      else if (kv.first == "mu_refitD0pismu_pt") muPt = kv.second;
      else if (kv.first == "mu_refitD0pismu_eta") muEta = kv.second;
      else if (kv.first == "mu_refitD0pismu_phi") muPhi = kv.second;
      else if (kv.first == "tksAdd_charge") AddTkCharge = kv.second;
      else if (kv.first == "tksAdd_pt") AddTkPt = kv.second;
      else if (kv.first == "tksAdd_eta") AddTkEta = kv.second;
      else if (kv.first == "tksAdd_phi") AddTkPhi = kv.second;
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

    // Get the pdgId of the cloasest MC particle to each reconstructed muon
    (*outputVecNtuplizer)["MC_mu_cloasestMCpart_pdgId"] = {};
    (*outputVecNtuplizer)["MC_mu_cloasestMCpart_dR"] = {};
    (*outputVecNtuplizer)["MC_mu_cloasestMCpart_dpt"] = {};

    for (auto muTLV : reco_muTLVs) {
      float best_dR = 1e9;
      float best_dpt = 1e9;
      int best_pdgId = 0;
      for(auto p : (*PrunedGenParticlesHandle)) {
        if (p.charge() == 0 || p.pt() < 5) continue;
        auto dR = vtxu::dR(muTLV.Phi(), p.phi(), muTLV.Eta(), p.eta());
        auto dpt = fabs( p.pt() - muTLV.Pt() ) / p.pt();

        if ( dR/0.01 + dpt/0.05 < best_dR/0.01 + best_dpt/0.05) {
          best_dR = dR;
          best_dpt = dpt;
          best_pdgId = p.pdgId();
        }
      }
      (*outputVecNtuplizer)["MC_mu_cloasestMCpart_pdgId"].push_back(best_pdgId);
      (*outputVecNtuplizer)["MC_mu_cloasestMCpart_dR"].push_back(best_dR);
      (*outputVecNtuplizer)["MC_mu_cloasestMCpart_dpt"].push_back(best_dpt);

      if (verbose) {
        cout << Form("Cand %.2f -> pdgId match %d: dR=%.4f dpt=%.4f", muTLV.Pt(), best_pdgId, best_dR, best_dpt) << endl;
      }
    }


    if(verbose) {
      cout << "Bottom hadrons MC history" << endl;
      for(unsigned int i = 0; i < N_PrunedGenParticles; i++) {
        auto p = (*PrunedGenParticlesHandle)[i];
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
    int i_cand = -1;
    float best_distance = 1e9;
    TLorentzVector bestMu, bestDst;
    int nB2DstMuX = 0;
    for(uint i = 0; i < N_PrunedGenParticles; i++) {
      auto p = (*PrunedGenParticlesHandle)[i];
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
        auto d = (*PrunedGenParticlesHandle)[j];
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
        auto mu = (*PrunedGenParticlesHandle)[jMu];
        for (auto jDst : idx_Dst) {
          auto Dst = (*PrunedGenParticlesHandle)[jDst];
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

      for (auto idxMuDst : idxMuDstPairs) {
        auto mu = (*PrunedGenParticlesHandle)[idxMuDst.first];
        auto Dst = (*PrunedGenParticlesHandle)[idxMuDst.second];

        for(uint j = 0; j < reco_muTLVs.size(); j++) {
          if ( mu.charge() != muCharge[j] ) continue; // If they do not have the right sign

          // delta R resolution of 0.01 for muon and double for D*, similar for pt
          float distance = vtxu::dR(reco_muTLVs[j].Phi(), mu.phi(), reco_muTLVs[j].Eta(), mu.eta()) / 0.01;
          distance += vtxu::dR(reco_DstTLVs[j].Phi(), Dst.phi(), reco_DstTLVs[j].Eta(), Dst.eta()) / 0.03;
          distance += ( fabs(mu.pt() - reco_muTLVs[j].Pt()) / mu.pt() ) / 0.03;
          distance += ( fabs(Dst.pt() - reco_DstTLVs[j].Pt()) / Dst.pt() ) / 0.10;

          if(distance < best_distance) {
            best_distance = distance;
            i_B = i;
            j_mu = idxMuDst.first;
            j_Dst = idxMuDst.second;
            i_cand = j;
            bestDst = reco_DstTLVs[j];
            bestMu = reco_muTLVs[j];

            if (verbose) {
              cout << Form("Best pair updated to mu:%d Dst:%d with candidate %d", j_mu, j_Dst, i_cand) << endl;
              cout << Form ( "Mu dR=%.4f", vtxu::dR(reco_muTLVs[j].Phi(), mu.phi(), reco_muTLVs[j].Eta(), mu.eta()) ) << endl;
              cout << Form ( "Mu dpt/pt=%.4f", ( fabs(mu.pt() - reco_muTLVs[j].Pt()) / mu.pt() ) ) << endl;
              cout << Form ( "Dst dR=%.4f", vtxu::dR(reco_DstTLVs[j].Phi(), Dst.phi(), reco_DstTLVs[j].Eta(), Dst.eta()) ) << endl;
              cout << Form ( "Dst dpt/pt=%.4f", ( fabs(Dst.pt() - reco_DstTLVs[j].Pt()) / Dst.pt() ) ) << endl;
            }
          }
        }
      }

    }

    (*outputNtuplizer)["MC_nB2DstMuX"] = nB2DstMuX;
    (*outputNtuplizer)["MC_idxCand"] = i_cand;
    (*indexBmc) = i_B;


    float recoMuon_match_BMuon = 0;
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
      auto p = (*PrunedGenParticlesHandle)[i_B];
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

      auto mu = (*PrunedGenParticlesHandle)[j_mu];
      p4["mu"].SetPtEtaPhiM(mu.pt(), mu.eta(), mu.phi(), mu.mass());
      vtx["mu"] = mu.vertex();
      (*outputNtuplizer)["MC_mu_TransvIP_PV"] = vtxu::computeIP(interactionPoint, mu.vertex(), mu.momentum(), true);
      auto muMother = mu.mother();
      (*outputNtuplizer)["MC_muMotherPdgId"] = muMother->pdgId();
      for(uint iss = 0; iss <  muMother->numberOfDaughters(); iss++) {
        auto muSis = muMother->daughter(iss);
        if ( muSis->pdgId() != mu.pdgId() ) (*outputVecNtuplizer)["MC_muSistersPdgId"].push_back(muSis->pdgId());
      }

      auto Dst = (*PrunedGenParticlesHandle)[j_Dst];
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

      double dR_mu = vtxu::dR(bestMu.Phi(), mu.phi(), bestMu.Eta(), mu.eta());
      double dPt_rel_mu = (bestMu.Pt() - mu.pt())/mu.pt();
      if (dR_mu < 0.02 && fabs(dPt_rel_mu) < 0.1 ) recoMuon_match_BMuon = 1.;

      (*outputNtuplizer)["MC_mu_TransvIP_vtxDst"] = vtxu::computeDCA_linApprox(vtx["Dst"], vtx["mu"], mu.momentum(), true);
      (*outputNtuplizer)["MC_mu_IP_vtxDst"] = vtxu::computeDCA_linApprox(vtx["Dst"], vtx["mu"], mu.momentum(), false);
    }

    (*outputNtuplizer)["recoMuon_match_BMuon"] = recoMuon_match_BMuon;

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
      cout << "Muon matched: " << recoMuon_match_BMuon << endl;
      cout << Form("MC muon pt=%.1f eta=%.2f IP=%1.2e", p4["mu"].Pt(), p4["mu"].Eta(), (*outputNtuplizer)["MC_mu_TransvIP_PV"]) << endl;
    }

    (*outputNtuplizer)["MC_M_vis"] = (p4["Dst"] + p4["mu"]).M();
    (*outputNtuplizer)["MC_M2_miss"] = (p4["B"] - p4["Dst"] - p4["mu"]).M2();
    (*outputNtuplizer)["MC_q2"] = (p4["B"] - p4["Dst"]).M2();

    TLorentzVector p4st_mu(p4["mu"]);
    p4st_mu.Boost(-1*p4["B"].BoostVector());
    (*outputNtuplizer)["MC_Est_mu"] = p4st_mu.E();


    (*outputVecNtuplizer)["MC_addTkFlag"] = {};
    (*outputVecNtuplizer)["MC_addTk_dEta"] = {};
    (*outputVecNtuplizer)["MC_addTk_dPhi"] = {};
    (*outputVecNtuplizer)["MC_addTk_dPt"] = {};
    (*outputVecNtuplizer)["MC_addTk_dz"] = {};
    (*outputVecNtuplizer)["MC_addTk_dxy"] = {};
    (*outputVecNtuplizer)["MC_addTk_pdgId"] = {};
    (*outputVecNtuplizer)["MC_addTk_pdgIdMother"] = {};
    (*outputVecNtuplizer)["MC_addTk_pdgIdMotherMother"] = {};
    for (uint i = 0; i < AddTkCharge.size(); i++) {
        (*outputVecNtuplizer)["MC_addTkFlag"].push_back(-1);
        (*outputVecNtuplizer)["MC_addTk_dEta"].push_back(0);
        (*outputVecNtuplizer)["MC_addTk_dPhi"].push_back(0);
        (*outputVecNtuplizer)["MC_addTk_dPt"].push_back(-1);
        (*outputVecNtuplizer)["MC_addTk_dz"].push_back(0);
        (*outputVecNtuplizer)["MC_addTk_dxy"].push_back(0);
        (*outputVecNtuplizer)["MC_addTk_pdgId"].push_back(0);
        (*outputVecNtuplizer)["MC_addTk_pdgIdMother"].push_back(0);
        (*outputVecNtuplizer)["MC_addTk_pdgIdMotherMother"].push_back(0);
    }
    if (AddTkCharge.size() > 0) {
      // Get packedGenParticles
      edm::Handle<std::vector<pat::PackedGenParticle>> PackedGenParticlesHandle;
      iEvent.getByToken(PackedParticlesSrc_, PackedGenParticlesHandle);
      // unsigned int N_PackedGenParticles = PackedGenParticlesHandle->size();
      for(auto packGenP : (*PackedGenParticlesHandle)) {
        if (packGenP.charge() == 0) continue;
        if (packGenP.pt() < 0.2) continue;
        for (uint i = 0; i < AddTkCharge.size(); i++) {
          if ((*outputVecNtuplizer)["MC_addTkFlag"][i] != -1) continue;
          if (packGenP.charge() != AddTkCharge[i]) continue;

          float dEta = fabs(packGenP.eta() - AddTkEta[i]);
          float dPhi = fabs(vtxu::dPhi(packGenP.phi(), AddTkPhi[i]));
          float dPt = fabs(packGenP.pt() - AddTkPt[i])/packGenP.pt();
          if (hypot(dEta, dPhi) < 0.002 && dPt < 0.03) {
            (*outputVecNtuplizer)["MC_addTkFlag"][i] = 1;
            (*outputVecNtuplizer)["MC_addTk_dEta"][i] = dEta;
            (*outputVecNtuplizer)["MC_addTk_dPhi"][i] = dPhi;
            (*outputVecNtuplizer)["MC_addTk_dPt"][i] = dPt;
            (*outputVecNtuplizer)["MC_addTk_dz"][i] = packGenP.dz();
            (*outputVecNtuplizer)["MC_addTk_dxy"][i] = packGenP.dxy();
            (*outputVecNtuplizer)["MC_addTk_pdgId"][i] = packGenP.pdgId();
            (*outputVecNtuplizer)["MC_addTk_pdgIdMother"][i] = packGenP.mother(0)->pdgId();
            if (packGenP.mother(0)->mother(0)) {
              (*outputVecNtuplizer)["MC_addTk_pdgIdMotherMother"][i] = packGenP.mother(0)->mother(0)->pdgId();
            }
          }
        }
      }
    }


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
