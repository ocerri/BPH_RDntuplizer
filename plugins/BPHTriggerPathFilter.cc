// system include files
#include <memory>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/PatCandidates/interface/Muon.h"

using namespace std;

class BPHTriggerPathFilter : public edm::stream::EDFilter<> {
   public:
      explicit BPHTriggerPathFilter(const edm::ParameterSet&);
      ~BPHTriggerPathFilter() {
        cout << Form("BPHTrigger filter efficiency: %d/%d = %1.2e", N_passed_events, N_analyzed_events, (double)N_passed_events/N_analyzed_events) << endl;
      };

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;
      // ----------member data ---------------------------

      edm::EDGetTokenT<vector<pat::Muon>> muonSrc_;
      int N_analyzed_events = 0;
      int N_passed_events = 0;
      int verbose = 0;
};

BPHTriggerPathFilter::BPHTriggerPathFilter(const edm::ParameterSet& iConfig):
  muonSrc_( consumes<vector<pat::Muon>> ( iConfig.getParameter<edm::InputTag>( "trgMuons" ) ) )
{
   //now do what ever initialization is needed
}

// ------------ method called on each new Event  ------------
bool
BPHTriggerPathFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  N_analyzed_events++;

  edm::Handle<std::vector<pat::Muon>> muonHandle;
  iEvent.getByToken(muonSrc_, muonHandle);
  unsigned int muonNumber = muonHandle->size();

  if(muonNumber > 0) {
    N_passed_events++;
    return true;
  }
  else return false;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
BPHTriggerPathFilter::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
BPHTriggerPathFilter::endStream() {
}

DEFINE_FWK_MODULE(BPHTriggerPathFilter);
