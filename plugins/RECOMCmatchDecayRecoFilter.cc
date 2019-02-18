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

class RECOMCmatchDecayRecoFilter : public edm::stream::EDFilter<> {
   public:
      explicit RECOMCmatchDecayRecoFilter(const edm::ParameterSet&);
      ~RECOMCmatchDecayRecoFilter() {};

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;
      // ----------member data ---------------------------

      int verbose = 0;
};

RECOMCmatchDecayRecoFilter::RECOMCmatchDecayRecoFilter(const edm::ParameterSet& iConfig)
{
  verbose = iConfig.getParameter<int>( "verbose" );
  consumesMany<map<string, float>>();
}

// ------------ method called on each new Event  ------------
bool
RECOMCmatchDecayRecoFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  vector< edm::Handle<map<string, float>> > outMapHandle;
  iEvent.getManyByType(outMapHandle);
  for( auto h : outMapHandle ) {
    for( auto& kv : *(h.product()) ) {
        if (kv.first == "NMatchedPart") {
           if (kv.second == 4) return true;
           else return false;
        }
    }
  }

  return false;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
RECOMCmatchDecayRecoFilter::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
RECOMCmatchDecayRecoFilter::endStream() {
}

DEFINE_FWK_MODULE(RECOMCmatchDecayRecoFilter);
