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

#include "VtxUtils.cc"

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
        }
        for(auto d : p.daughterRefVector()) {
          if(d->pdgId() == -13) {
            auto match_res = RecoMC_matching(*d, *pfCandHandle);
            auto dR = get<1>(match_res);
            auto dpt = get<2>(match_res);
            if (verbose) {
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
            for(auto dd : d->daughterRefVector()) {
              if(dd->pdgId() == -211) {
                auto match_res = RecoMC_matching(*dd, *pfCandHandle);
                auto dR = get<1>(match_res);
                auto dpt = get<2>(match_res);
                if (verbose) {
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

    if (matchedPart.size() == 4) {
      vtxu::FitD0(iSetup, matchedPart["pi"], matchedPart["K"], 1);
    }

    unique_ptr<map<string, float>> outputNtuplizer(new map<string, float>);

    // This was a debug
    // (*outputNtuplizer)["matchMu"] = matchedPart.count("mu");
    // (*outputNtuplizer)["matchPisoft"] = matchedPart.count("pisoft");
    // (*outputNtuplizer)["matchPi"] = matchedPart.count("pi");
    // (*outputNtuplizer)["matchK"] = matchedPart.count("K");

    // Needed to filter out the events with no matched particles
    (*outputNtuplizer)["NMatchedPart"] = matchedPart.size();

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
