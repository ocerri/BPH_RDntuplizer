#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
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
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
// #include "Geometry/CommonDetUnit/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetType.h"



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
    // edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> trackGeomSrc_;
    int verbose = 0;
};



TracksHitsAnalyzer::TracksHitsAnalyzer(const edm::ParameterSet &iConfig)
// : trackGeomSrc_(esConsumes())
{
    trackSrc_ = consumes<vector<reco::Track>>(edm::InputTag("standAloneMuons"));
    verbose = iConfig.getParameter<int>( "verbose" );
}


void TracksHitsAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    edm::Handle<vector<reco::Track>> trackHandle;
    iEvent.getByToken(trackSrc_, trackHandle);

    edm::ESHandle<TransientTrackingRecHitBuilder> theTrackerRecHitBuilder;
    // iSetup.get<TransientRecHitRecord>().get("WithTrackAngle",theTrackerRecHitBuilder);
    iSetup.get<TransientRecHitRecord>().get("WithoutRefit",theTrackerRecHitBuilder);

    // edm::ESHandle<TrackerGeometry> theGeometry;
    // iSetup.get<TrackerGeometryRecord>(theGeometry);
    // auto const& pDD = iSetup.getData(trackGeomSrc_);

    for(uint i_tk = 0; i_tk<trackHandle->size(); i_tk++) {
      auto tk = (*trackHandle)[i_tk];
      cout << "Track " << i_tk << ", nHits: " << tk.recHitsSize() << endl;
      for(uint i_hit = 0; i_hit < tk.recHitsSize(); i_hit++) {
        auto hit = tk.recHit(i_hit);
        if (hit->isValid()) {
          if (hit->getType() != 0) continue;
          cout << "Type: " << hit->getType() << ", localPosition " << hit->localPosition() << endl;
          cout << hit->geographicalId().rawId() << endl;
          cout << hit->geographicalId().subdetId() << endl;

          // auto module = hit->det();
          // auto gPos = module->toGlobal(hit->localPosition());
          // auto r = hypot(gPos.x(), gPos.y())
          // cout << "valid hit found with global position = "<< gPos << ", r=" << r << endl;

          cout << "Building transient hit " << i_hit << "..." << flush;
          TransientTrackingRecHit::RecHitPointer tthit = theTrackerRecHitBuilder->build(&*hit);
          cout << "built..." << flush;
          if (tthit->isValid()) {
            cout << tthit->det() << endl;
            cout << "local Position" << flush;
            cout << tthit->localPosition() << endl;
            try {
              GlobalPoint gPosition =  tthit->globalPosition();
              cout << "valid hit found with global position = "<< gPosition << endl;
            }
            catch (...) {cout << "no globalPosition" << endl;}
          }
          else { cout << "not valid." << endl;}
        }
      }
    }

}

void TracksHitsAnalyzer::beginJob(const edm::EventSetup&) { }
void TracksHitsAnalyzer::endJob() { }


DEFINE_FWK_MODULE(TracksHitsAnalyzer);
