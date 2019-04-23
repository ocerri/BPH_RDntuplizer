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

class MuDstDecayTreeFilter : public edm::stream::EDFilter<> {
   public:
      explicit MuDstDecayTreeFilter(const edm::ParameterSet&);
      ~MuDstDecayTreeFilter() {};

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;
      // ----------member data ---------------------------

      edm::EDGetTokenT<map<string,float>> muonSrc_;
      int verbose = 0;
};

MuDstDecayTreeFilter::MuDstDecayTreeFilter(const edm::ParameterSet& iConfig)
{
  consumesMany<map<string, float>>();
}

// ------------ method called on each new Event  ------------
bool MuDstDecayTreeFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  vector< edm::Handle<map<string, float>> > outMapHandle;
  iEvent.getManyByType(outMapHandle);
  for( auto h : outMapHandle ) {
    for( auto& kv : *(h.product()) ) {
      if (kv.first == "n_pis") {
          if(verbose) {cout << Form("Filter found %.0f soft pions vertexing with a D0\n", kv.second);}
          if (kv.second > 0) return true;
          else return false;
      }
    }
  }

  return false;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
MuDstDecayTreeFilter::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
MuDstDecayTreeFilter::endStream() {
}

DEFINE_FWK_MODULE(MuDstDecayTreeFilter);
