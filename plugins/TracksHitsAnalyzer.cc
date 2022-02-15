#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

using namespace std;

class TracksHitsAnalyzer : public edm::EDAnalyzer {

public:

    explicit TracksHitsAnalyzer(const edm::ParameterSet &);
    ~TracksHitsAnalyzer() {};

private:
    virtual void beginJob(const edm::EventSetup&) ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() ;

    // ----------member data ---------------------------
    edm::EDGetTokenT<vector<reco::Track>> trackSrc_;
    int verbose = 0;
};



TracksHitsAnalyzer::TracksHitsAnalyzer(const edm::ParameterSet &iConfig)
{
    trackSrc_ = consumes<vector<reco::Track>>(edm::InputTag("generalTracks"));
    verbose = iConfig.getParameter<int>( "verbose" );
}


void TracksHitsAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    edm::Handle<vector<reco::Track>> trackHandle;
    iEvent.getByToken(trackSrc_, trackHandle);

    edm::ESHandle<TransientTrackingRecHitBuilder> theTrackerRecHitBuilder;
    // iSetup.get<TransientRecHitRecord>().get("WithTrackAngle",theTrackerRecHitBuilder);
    iSetup.get<TransientRecHitRecord>().get("WithoutRefit",theTrackerRecHitBuilder);

    for(uint i_tk = 0; i_tk<trackHandle->size(); i_tk++) {
      auto tk = (*trackHandle)[i_tk];
      cout << "Track " << i_tk << ", nHits: " << tk.recHitsSize() << endl;
      for(uint i_hit = 0; i_hit < tk.recHitsSize(); i_hit++) {
        auto hit = tk.recHit(i_hit);
        if (hit->isValid()) {
          cout << "Building transient hit " << i_hit << "..." << flush;
          TransientTrackingRecHit::RecHitPointer tthit = theTrackerRecHitBuilder->build(&*hit);
          GlobalPoint gPosition =  tthit->globalPosition();
          cout << "valid hit found with global position = "<< gPosition << endl;
        }
      }
    }

}

void TracksHitsAnalyzer::beginJob(const edm::EventSetup&) { }
void TracksHitsAnalyzer::endJob() { }


DEFINE_FWK_MODULE(TracksHitsAnalyzer);
