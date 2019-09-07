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

using namespace std;

class B2DstKDecayTreeFilter : public edm::stream::EDFilter<> {
   public:
      explicit B2DstKDecayTreeFilter(const edm::ParameterSet&);
      ~B2DstKDecayTreeFilter() {
        cout << Form("B2DstKDecay filter efficiency: %d/%d = %1.2e", N_passed_events, N_analyzed_events, (double)N_passed_events/N_analyzed_events) << endl;
      };

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<map<string, float>> B2DstKDecayTreeOutSrc_;
      int N_analyzed_events = 0;
      int N_passed_events = 0;
      int verbose = 0;
};

B2DstKDecayTreeFilter::B2DstKDecayTreeFilter(const edm::ParameterSet& iConfig)
{
  verbose = iConfig.getParameter<int>( "verbose" );
  B2DstKDecayTreeOutSrc_ = consumes<map<string, float>>(edm::InputTag("B2DstKDT", "outputNtuplizer"));
}

// ------------ method called on each new Event  ------------
bool B2DstKDecayTreeFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  N_analyzed_events++;

  edm::Handle<map<string, float>> outMapHandle;
  iEvent.getByToken(B2DstKDecayTreeOutSrc_, outMapHandle);

  for( auto const& kv : (*outMapHandle) ) {
    if (kv.first == "n_B") {
        if(verbose) {cout << Form("Filter found %.0f B candidates\n", kv.second);}
        if (kv.second > 0) {
          N_passed_events++;
          return true;
        }
        else return false;
    }
  }

  return false;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
B2DstKDecayTreeFilter::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
B2DstKDecayTreeFilter::endStream() {
}

DEFINE_FWK_MODULE(B2DstKDecayTreeFilter);
