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

class B2DstMuDecayTreeFilter : public edm::stream::EDFilter<> {
   public:
      explicit B2DstMuDecayTreeFilter(const edm::ParameterSet&);
      ~B2DstMuDecayTreeFilter() {};

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;
      // ----------member data ---------------------------

      edm::EDGetTokenT<map<string,float>> muonSrc_;
      int verbose = 0;
};

B2DstMuDecayTreeFilter::B2DstMuDecayTreeFilter(const edm::ParameterSet& iConfig)
{
  verbose = iConfig.getParameter<int>( "verbose" );
  consumesMany<map<string, float>>();
}

// ------------ method called on each new Event  ------------
bool B2DstMuDecayTreeFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  vector< edm::Handle<map<string, float>> > outMapHandle;
  iEvent.getManyByType(outMapHandle);
  for( auto h : outMapHandle ) {
    for( auto& kv : *(h.product()) ) {
      if (kv.first == "n_B") {
          if(verbose) {cout << Form("Filter found %.0f B candidates\n", kv.second);}
          if (kv.second > 0) return true;
          else return false;
      }
    }
  }

  return false;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
B2DstMuDecayTreeFilter::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
B2DstMuDecayTreeFilter::endStream() {
}

DEFINE_FWK_MODULE(B2DstMuDecayTreeFilter);
