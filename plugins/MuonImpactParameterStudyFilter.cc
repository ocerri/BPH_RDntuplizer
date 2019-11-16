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

class MuonImpactParameterStudyFilter : public edm::stream::EDFilter<> {
   public:
      explicit MuonImpactParameterStudyFilter(const edm::ParameterSet&);
      ~MuonImpactParameterStudyFilter() {
        cout << Form("MuonImpactParameterStudyFilter filter efficiency: %d/%d = %1.2e", N_passed_events, N_analyzed_events, (double)N_passed_events/N_analyzed_events) << endl;
      };

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<map<string, float>> B2DstMuDecayTreeOutSrc_;
      int N_analyzed_events = 0;
      int N_passed_events = 0;
      int verbose = 0;
};

MuonImpactParameterStudyFilter::MuonImpactParameterStudyFilter(const edm::ParameterSet& iConfig)
{
  verbose = iConfig.getParameter<int>( "verbose" );
  B2DstMuDecayTreeOutSrc_ = consumes<map<string, float>>(edm::InputTag("MuProd", "outputNtuplizer"));
}

// ------------ method called on each new Event  ------------
bool MuonImpactParameterStudyFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  N_analyzed_events++;

  edm::Handle<map<string, float>> outMapHandle;
  iEvent.getByToken(B2DstMuDecayTreeOutSrc_, outMapHandle);

  for( auto const& kv : (*outMapHandle) ) {
    if (kv.first == "n_mu") {
        if(verbose) {cout << Form("Filter found %.0f muons\n", kv.second);}
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
MuonImpactParameterStudyFilter::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
MuonImpactParameterStudyFilter::endStream() {
}

DEFINE_FWK_MODULE(MuonImpactParameterStudyFilter);
