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
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include <iostream>
#include <string>
#include <map>

#include "VtxUtils.hh"

using namespace std;

class RECOMCmatchDecayRecoProducer : public edm::EDProducer {

public:

    explicit RECOMCmatchDecayRecoProducer(const edm::ParameterSet &iConfig);

    ~RECOMCmatchDecayRecoProducer() override {};

    tuple<pat::PackedCandidate, double, double> RecoMC_matching(reco::GenParticle, vector<pat::PackedCandidate>);

private:

    virtual void produce(edm::Event&, const edm::EventSetup&);

    // ----------member data ---------------------------

    edm::EDGetTokenT<vector<reco::GenParticle>> PrunedParticlesSrc_;
    edm::EDGetTokenT<vector<pat::PackedCandidate>> PFCandSrc_;

    int verbose = 0;
};



RECOMCmatchDecayRecoProducer::RECOMCmatchDecayRecoProducer(const edm::ParameterSet &iConfig)
{
    verbose = iConfig.getParameter<int>( "verbose" );
    PFCandSrc_ = consumes<vector<pat::PackedCandidate>>(edm::InputTag("packedPFCandidates"));
    PrunedParticlesSrc_ = consumes<vector<reco::GenParticle>>(edm::InputTag("prunedGenParticles"));
    // produces<vector<pat::PackedCandidate>>("recoCandMatched");
    // produces<vector<string>>("recoCandMatchedNames");
    produces<map<string, float>>("outputNtuplizer");
}


void RECOMCmatchDecayRecoProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    if (verbose) {cout << "-----------  MC matching ----------\n";}

    // Get prunedGenParticles
    edm::Handle<std::vector<reco::GenParticle>> PrunedGenParticlesHandle;
    iEvent.getByToken(PrunedParticlesSrc_, PrunedGenParticlesHandle);

    edm::Handle<vector<pat::PackedCandidate>> pfCandHandle;
    iEvent.getByToken(PFCandSrc_, pfCandHandle);

    // Output collection
    unique_ptr<map<string, float>> outputNtuplizer(new map<string, float>);
    // unique_ptr<vector<pat::PackedCandidate>> RECO_MCmatch( new vector<pat::PackedCandidate> );
    // unique_ptr<vector<string>> RECO_MCmatchNames( new vector<string> );

    map<string, pat::PackedCandidate> matchedPart;
    map<string, pair<double, double>> matchedPartDist;

    for(auto p : *PrunedGenParticlesHandle) {
      if (p.pdgId() == 511 && p.numberOfDaughters()>1) {
        if(verbose) {
          cout << "Found: " << p.pdgId() << " (-->" << flush;
          for(auto d : p.daughterRefVector()) {
            cout << " " << d->pdgId() << flush;
          }
          cout << ')' << endl;
          cout << Form("Generated at: {%.4f, %.4f, %.4f}\n", p.vx(), p.vy(), p.vz());
        }
        for(auto d : p.daughterRefVector()) {
          if(d->pdgId() == -13) {
            auto match_res = RecoMC_matching(*d, *pfCandHandle);
            auto dR = get<1>(match_res);
            auto dpt = get<2>(match_res);
            if (verbose) {
              cout << Form("Mu generated at: {%.4f, %.4f, %.4f}\n", d->vx(), d->vy(), d->vz());
              cout << "Mu dR: " << dR << " dpt: " << dpt << endl;
            }
            if(dR > 0) {
              auto p_reco = get<0>(match_res);
              if(p_reco.isStandAloneMuon() && dR < 0.015 && fabs(dpt) < 0.05) {
                matchedPart["mu"] = p_reco;
              }
            }
          }
          else if(d->pdgId() == -413) {
            if (verbose) {
              cout << Form("D*- generated at: {%.4f, %.4f, %.4f}\n", d->vx(), d->vy(), d->vz());
            }
            for(auto dd : d->daughterRefVector()) {
              if(dd->pdgId() == -211) {
                auto match_res = RecoMC_matching(*dd, *pfCandHandle);
                auto dR = get<1>(match_res);
                auto dpt = get<2>(match_res);
                if (verbose) {
                  cout << Form("Softpi generated at: {%.4f, %.4f, %.4f}\n", d->vx(), d->vy(), d->vz());
                  cout << "Softpi dR: " << dR << " dpt: " << dpt << endl;

                }
                if(dR > 0) {
                  auto p_reco = get<0>(match_res);
                  if(dR < 0.02 && fabs(dpt) < 0.05) {
                    matchedPart["pisoft"] = p_reco;
                  }
                }
              }
              else if(abs(dd->pdgId()) == 421) {
                for(auto ddd : dd->daughterRefVector()) {
                  if (ddd->pdgId() == -211 || ddd->pdgId() == 321)
                  {
                    string name = ddd->pdgId() == -211 ? "pi" : "K";
                    auto match_res = RecoMC_matching(*ddd, *pfCandHandle);
                    auto dR = get<1>(match_res);
                    auto dpt = get<2>(match_res);
                    if (verbose) {
                      cout << name.c_str() << " dR: " << dR << " dpt: " << dpt << endl;
                    }
                    if(dR > 0) {
                      auto p_reco = get<0>(match_res);
                      if(dR < 0.02 && fabs(dpt) < 0.05) {
                        matchedPart[name] = p_reco;
                      }
                    }
                  }
                }
              }
            }
          }
        }

        break;
      }
    }

    bool matched_K = false;
    bool matched_pi = false;
    bool matched_pisoft = false;
    bool matched_mu = false;
    // cout << "Evt" << endl;
    // cout << matchedPart.size() << ":";
    for(auto kv : matchedPart) {
      // cout << " " << kv.first;
      if(kv.first == "K") {
        if(kv.second.bestTrack() != 0) {
          matched_K =  true;
        }
      }
      else if(kv.first == "pi") {
        if(kv.second.bestTrack() != 0) {
          matched_pi =  true;
        }
      }
      else if(kv.first == "pisoft") {
        if(kv.second.bestTrack() != 0) {
          matched_pisoft =  true;
        }
      }
      else if(kv.first == "mu") {
        if(kv.second.bestTrack() != 0) {
          matched_mu =  true;
        }
      }

      if (verbose) {
        cout << kv.first << endl;
        cout << "PV assoc: " << kv.second.fromPV() << endl;
        cout << Form("part.v at: {%.4f, %.4f, %.4f}\n", kv.second.vx(), kv.second.vy(), kv.second.vz());
        auto vtx = kv.second.vertex();
        cout << Form("part.v at: {%.4f, %.4f, %.4f}\n", vtx.x(), vtx.y(), vtx.z());
        auto vtx_ptr = kv.second.vertexRef();
        cout << Form("vtx_ptr at: {%.4f, %.4f, %.4f}\n", vtx_ptr->x(), vtx_ptr->y(), vtx_ptr->z());
      }
    }
    // cout << endl;
    // if(matched_pi && !(matchedPart["pi"].bestTrack())) matched_pi = false;
    // cout << matchedPart["pi"].bestTrack() << endl;
    // if(matched_K && !(matchedPart["K"].bestTrack())) matched_K = false;
    // cout << matchedPart["K"].bestTrack() << endl;


    if (matched_pi && matched_K) {
      // cout << "Matched pi and K\n";
      (*outputNtuplizer)["D0prefit_mass"] = (matchedPart["pi"].p4() + matchedPart["K"].p4()).M();

      auto D0KinTree = vtxu::FitD0(iSetup, matchedPart["pi"], matchedPart["K"], false, 0);
      if(!D0KinTree->isValid()) return;

      D0KinTree->movePointerToTheTop();
      auto D0vtx = D0KinTree->currentDecayVertex();

      (*outputNtuplizer)["D0vtx_chi2"] = D0vtx->chiSquared();
      if (verbose) {
        cout << "D0vtx->chiSquared(): " << D0vtx->chiSquared() << endl;
      }

      auto D0_reco = D0KinTree->currentParticle()->currentState();
      (*outputNtuplizer)["D0reco_mass"] = D0_reco.mass();
      (*outputNtuplizer)["D0reco_massErr"] = sqrt(D0_reco.kinematicParametersError().matrix()(6,6));

      D0KinTree->movePointerToTheTop();
      D0KinTree->movePointerToTheFirstChild();
      auto K_refitD0 = D0KinTree->currentParticle()->refittedTransientTrack().track();
      D0KinTree->movePointerToTheNextChild();
      auto pi_refitD0 = D0KinTree->currentParticle()->refittedTransientTrack().track();

      (*outputNtuplizer)["PiKrefitD0_mass"] = (vtxu::getTLVfromTrack(pi_refitD0, 0.13957018) + vtxu::getTLVfromTrack(K_refitD0, 0.493677)).M();

      cout << "---DEBUG---\n";
      auto v1 = vtxu::getTLVfromTrack(pi_refitD0, 0.13957018);
      auto v2 = vtxu::getTLVfromKinPart(D0KinTree->currentParticle());
      v1.Print();
      v2.Print();

      if(matched_pisoft) {
        (*outputNtuplizer)["Dstprefit_mass"] = (matchedPart["pisoft"].p4() + matchedPart["pi"].p4() + matchedPart["K"].p4()).M();

        auto DstKinTree = vtxu::FitDst(iSetup, matchedPart["pisoft"], matchedPart["pi"], matchedPart["K"], false, 0);
        if(!DstKinTree->isValid()) return;
        DstKinTree->movePointerToTheTop();
        auto Dst_reco = D0KinTree->currentParticle()->currentState();
        (*outputNtuplizer)["Dstreco_mass"] = Dst_reco.mass();
        (*outputNtuplizer)["Dstreco_massErr"] = sqrt(Dst_reco.kinematicParametersError().matrix()(6,6));

        if (matched_mu) {
          // Refit with D* mass constraint
          auto DstKinTree = vtxu::FitDst(iSetup, matchedPart["pisoft"], matchedPart["pi"], matchedPart["K"], true, 0);
          if(!DstKinTree->isValid()) return;
          DstKinTree->movePointerToTheTop();

          // Get B flight direction
          auto Dstvtx = DstKinTree->currentDecayVertex()->position();
          auto pVtx = matchedPart["mu"].vertexRef();
          TVector3 d_flightB(Dstvtx.x()-pVtx->x(), Dstvtx.y()-pVtx->y(), Dstvtx.z()-pVtx->z());

          // auto Dst_reco = D0KinTree->currentParticle()->currentState();
        }

      }
    }

    // This was a debug
    // (*outputNtuplizer)["matchMu"] = matchedPart.count("mu");
    // (*outputNtuplizer)["matchPisoft"] = matchedPart.count("pisoft");
    // (*outputNtuplizer)["matchPi"] = matchedPart.count("pi");
    // (*outputNtuplizer)["matchK"] = matchedPart.count("K");

    // Needed to filter out the events with no matched particles
    (*outputNtuplizer)["NMatchedPart"] = matchedPart.size() * (matched_pi && matched_K && matched_pisoft);
    // cout << (*outputNtuplizer)["NMatchedPart"] << " " << matchedPart.size() << endl;
    // cout << " ------------- \n";

    iEvent.put(move(outputNtuplizer), "outputNtuplizer");
    return;
}

tuple<pat::PackedCandidate, double, double> RECOMCmatchDecayRecoProducer::RecoMC_matching(reco::GenParticle p_MC, vector<pat::PackedCandidate> reco_coll) {
  pat::PackedCandidate p_best;
  double dR_best = -1, dpt_best = 0;

  for(auto p : reco_coll) {
    if(p.charge() == p_MC.charge()) {
      double dR = hypot(p_MC.phi()-p.phi(), p_MC.eta()-p.eta());
      double dpt = (p.pt() - p_MC.pt()) / p_MC.pt();
      if(dR_best == -1 || dR < dR_best || (fabs(dR-dR_best) < 1e-3 && dpt < dpt_best)) {
        p_best = p;
        dR_best = dR;
        dpt_best = dpt;
      }
    }
  }

  tuple<pat::PackedCandidate, double, double> out = {p_best, dR_best, dpt_best};
  return out;
}


DEFINE_FWK_MODULE(RECOMCmatchDecayRecoProducer);
