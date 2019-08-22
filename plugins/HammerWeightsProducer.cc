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

    double mass_B0 = 5.27961;
    double mass_mu = 0.1056583745;
    double mass_K = 0.493677;
    double mass_pi = 0.13957018;

    int verbose = 0;
};



HammerWeightsProducer::HammerWeightsProducer(const edm::ParameterSet &iConfig)
{
    verbose = iConfig.getParameter<int>( "verbose" );
    PrunedParticlesSrc_ = consumes<vector<reco::GenParticle>>(edm::InputTag("prunedGenParticles"));
    indexBmcSrc_ = consumes<int> (edm::InputTag("MCpart", "indexBmc"));

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

    // Output collection
    unique_ptr<map<string, float>> outputNtuplizer(new map<string, float>);

    cout << "i_B retieved: " << i_B << endl;
    // if(i_B >= 0){
    //   auto p = (*PrunedGenParticlesHandle)[i_B];
    //
    //   for(auto d : p.daughterRefVector()) {
    //     if(d->pdgId() == -13) {
    //       p4["mu"].SetPtEtaPhiM(d->pt(), d->eta(), d->phi(), d->mass());
    //     }
    //     if(d->pdgId() == -15) {
    //       for (auto dd : d->daughterRefVector()) {
    //         if (dd->pdgId() == -13) {
    //           p4["mu"].SetPtEtaPhiM(dd->pt(), dd->eta(), dd->phi(), dd->mass());
    //         }
    //       }
    //     }
    //     else if(d->pdgId() == -413) {
    //       p4["Dst"].SetPtEtaPhiM(d->pt(), d->eta(), d->phi(), d->mass());
    //
    //       for(auto dd : d->daughterRefVector()) {
    //         if(dd->pdgId() == -211) {
    //           p4["pis"].SetPtEtaPhiM(dd->pt(), dd->eta(), dd->phi(), dd->mass());
    //         }
    //         else if(abs(dd->pdgId()) == 421) {
    //           p4["D0"].SetPtEtaPhiM(dd->pt(), dd->eta(), dd->phi(), dd->mass());
    //         }
    //       }
    //     }
    //   }
    // }

    (*outputNtuplizer)["test"] = 0.;

    iEvent.put(move(outputNtuplizer), "outputNtuplizer");
    return;
}


DEFINE_FWK_MODULE(HammerWeightsProducer);
