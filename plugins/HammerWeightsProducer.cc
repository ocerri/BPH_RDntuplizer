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

#include <iostream>
#include <string>
#include <map>

#include "Hammer/Hammer.hh"
#include "Hammer/Process.hh"
#include "Hammer/Particle.hh"

using namespace std;

class HammerWeightsProducer : public edm::EDProducer {

public:

    explicit HammerWeightsProducer(const edm::ParameterSet &iConfig);

    ~HammerWeightsProducer() override {};

private:

    virtual void produce(edm::Event&, const edm::EventSetup&);

    // ----------member data ---------------------------

    edm::EDGetTokenT<vector<reco::GenParticle>> PrunedParticlesSrc_;
    edm::EDGetTokenT<int> indexBmcSrc_;

    Hammer::Hammer hammer;

    double mass_B0 = 5.27961;
    double mass_mu = 0.1056583745;
    double mass_K = 0.493677;
    double mass_pi = 0.13957018;

    int verbose = 0;
};



HammerWeightsProducer::HammerWeightsProducer(const edm::ParameterSet &iConfig)
{
    PrunedParticlesSrc_ = consumes<vector<reco::GenParticle>>(edm::InputTag("prunedGenParticles"));
    indexBmcSrc_ = consumes<int> (edm::InputTag("MCpart", "indexBmc"));

    verbose = iConfig.getParameter<int>( "verbose" );

    auto decayOfInterest = iConfig.getParameter<vector<string>>( "decayOfInterest" );
    for(auto s : decayOfInterest) {
      if(verbose) {cout << "[Hammer]: Including decay " << s << endl;}
      hammer.includeDecay(s);
    }

    auto inputFFScheme_ = iConfig.getParameter<vector<string>>("inputFFScheme");
    if(verbose) {cout << "[Hammer]: Input scheme" << endl;}
    map<string, string> inputFFScheme;
    for(uint i = 0; i < inputFFScheme_.size(); i++) {
      if(i%2 == 1) continue;
      inputFFScheme[inputFFScheme_[i]] = inputFFScheme_[i+1];
      if(verbose){cout << "\t" << inputFFScheme_[i] << ": " << inputFFScheme_[i+1] << endl;}
    }
    hammer.setFFInputScheme(inputFFScheme);

    hammer.addFFScheme("Central", {{"BD", "BGLVar"}, {"BD*", "BGLVar"}});

    hammer.initRun();

    produces<map<string, float>>("outputNtuplizer");
}


void HammerWeightsProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    if (verbose) {cout << "-----------  Hammer weights ----------\n";}

    // Get prunedGenParticles
    edm::Handle<std::vector<reco::GenParticle>> PrunedGenParticlesHandle;
    iEvent.getByToken(PrunedParticlesSrc_, PrunedGenParticlesHandle);

    // Get the B index
    edm::Handle<int> indexBmcHandle;
    iEvent.getByToken(indexBmcSrc_, indexBmcHandle);
    int i_B = (*indexBmcHandle);
    if(i_B <0){
      cout << "Invalid B idx (i.e. no B MC set)" << endl;
      assert(false);
    }
    // cout << "i_B retieved: " << i_B << endl;

    // Output collection
    unique_ptr<map<string, float>> outputNtuplizer(new map<string, float>);

    // Initialize the Hammer event
    hammer.initEvent();
    Hammer::Process B2DstLNu_Dst2DPi;
    vector<size_t> Bvtx_idxs;
    int idxTau = -1;
    vector<size_t> Tauvtx_idxs;
    int idxDst = -1;
    vector<size_t> Dstvtx_idxs;

    auto p = (*PrunedGenParticlesHandle)[i_B];
    Hammer::Particle pB({p.energy(), p.px(), p.py(), p.pz()}, p.pdgId());
    auto idxB = B2DstLNu_Dst2DPi.addParticle(pB);
    for(auto d : p.daughterRefVector()) {
      Hammer::Particle B_dau({d->energy(), d->px(), d->py(), d->pz()}, d->pdgId());
      auto idx_d = B2DstLNu_Dst2DPi.addParticle(B_dau);
      Bvtx_idxs.push_back(idx_d);

      if(d->pdgId() == -15) {
        idxTau = idx_d;
        for (auto dd : d->daughterRefVector()) {
          Hammer::Particle Tau_dau({dd->energy(), dd->px(), dd->py(), dd->pz()}, dd->pdgId());
          auto idx_dd = B2DstLNu_Dst2DPi.addParticle(Tau_dau);
          Tauvtx_idxs.push_back(idx_dd);
        }
      }
      else if(d->pdgId() == -413) {
        idxDst = idx_d;
        for (auto dd : d->daughterRefVector()) {
          Hammer::Particle Dst_dau({dd->energy(), dd->px(), dd->py(), dd->pz()}, dd->pdgId());
          auto idx_dd = B2DstLNu_Dst2DPi.addParticle(Dst_dau);
          Dstvtx_idxs.push_back(idx_dd);
        }
      }
    }

    B2DstLNu_Dst2DPi.addVertex(idxB, Bvtx_idxs);
    if(idxTau != -1) {
      B2DstLNu_Dst2DPi.addVertex(idxTau, Tauvtx_idxs);
    }
    if(idxDst != -1) {
      B2DstLNu_Dst2DPi.addVertex(idxDst, Dstvtx_idxs);
    }

    hammer.addProcess(B2DstLNu_Dst2DPi);
    hammer.processEvent();

    auto weights = hammer.getWeights("Central");
    if(verbose) {
      cout << "Central:" << endl;
      for(auto elem: weights) {
        cout << "\t" << elem.first << " " << elem.second << endl;
      }
    }

    (*outputNtuplizer)["test"] = 0.;

    iEvent.put(move(outputNtuplizer), "outputNtuplizer");
    return;
}


DEFINE_FWK_MODULE(HammerWeightsProducer);
