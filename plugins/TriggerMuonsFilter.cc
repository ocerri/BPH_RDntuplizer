// system include files
#include <memory>
#include <iostream>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

using namespace std;

class TriggerMuonsFilter : public edm::stream::EDFilter<> {
   public:
      explicit TriggerMuonsFilter(const edm::ParameterSet&);
      ~TriggerMuonsFilter() {
        cout << Form("Muons trigger filter efficiency: %d/%d = %1.2e", N_passed_events, N_analyzed_events, (double)N_passed_events/N_analyzed_events) << endl;

        fs->make<TNamed>("TriggerMuonsFilterEfficiency", Form("%d/%d", N_passed_events, N_analyzed_events));
      };

   private:
      // void beginJob(const edm::EventSetup&) {};
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      // void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<vector<pat::Muon>> muonSrc_;
      edm::EDGetTokenT<vector<reco::Vertex>> vtxSrc_;

      edm::Service<TFileService> fs;
      TH1I* hAllNvts;
      TH1I* hAllVtxZ;

      int N_analyzed_events = 0;
      int N_passed_events = 0;
      int muonCharge_ = 0;
      int verbose = 0;
};


TriggerMuonsFilter::TriggerMuonsFilter(const edm::ParameterSet& iConfig):
  muonSrc_( consumes<vector<pat::Muon>> ( edm::InputTag("slimmedMuons") ) ),
  vtxSrc_( consumes<vector<reco::Vertex>> ( edm::InputTag("offlineSlimmedPrimaryVertices") ) ),
  muonCharge_( iConfig.getParameter<int>( "muon_charge" ) ),
  verbose( iConfig.getParameter<int>( "verbose" ) )
{
  produces<vector<pat::Muon>>("trgMuonsMatched");
  produces<map<string, float>>("outputNtuplizer");
  produces<map<string, vector<float>>>("outputVecNtuplizer");
  hAllNvts = fs->make<TH1I>("hAllNvts", "Number of vertexes from all the MINIAOD events", 101, -0.5, 100.5);
  hAllVtxZ = fs->make<TH1I>("hAllVtxZ", "Z coordinate of vertexes from all the MINIAOD events", 100, -25, 25);
}

bool TriggerMuonsFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  N_analyzed_events++;
  if (verbose) {cout << "Event " << N_analyzed_events << endl;}

  edm::Handle<vector<pat::Muon>> muonHandle;
  iEvent.getByToken(muonSrc_, muonHandle);

  edm::Handle<vector<reco::Vertex>> vtxHandle;
  iEvent.getByToken(vtxSrc_, vtxHandle);
  auto primaryVtx = (*vtxHandle)[0];
  hAllNvts->Fill((int)vtxHandle->size());
  for(auto vtx : (*vtxHandle)) hAllVtxZ->Fill(vtx.position().z());

  // Output collection
  unique_ptr<vector<pat::Muon>> trgMuonsMatched( new vector<pat::Muon> );
  unique_ptr<map<string, float>> outputNtuplizer(new map<string, float>);
  unique_ptr<map<string, vector<float>>> outputVecNtuplizer(new map<string, vector<float>>);

  vector<string> triggerTag = {"Mu12_IP6", "Mu9_IP5", "Mu7_IP4", "Mu9_IP4", "Mu8_IP5", "Mu8_IP6", "Mu9_IP6", "Mu8_IP3"};

  if (verbose) {cout << "\n MUONS LIST" << endl;}
  for (auto muon : (*muonHandle)) {
    if (verbose) {
      cout << "\t" << Form("id:%d  pt:%.1f  eta:%.1f  phi:%.1f softID:%d", muon.pdgId(), muon.pt(), muon.eta(), muon.phi(), muon.isSoftMuon(primaryVtx)) << endl;
    }
    if(muonCharge_ != 0 && muon.charge() != muonCharge_) continue;
    if(!muon.triggered("HLT_Mu*_IP*_part*_v*")) continue;
    trgMuonsMatched->push_back(muon);

    for(auto tag : triggerTag) {
      string trgPath = "HLT_" + tag + "_part*_v*";
      (*outputVecNtuplizer)["trgMu_HLT_" + tag].push_back(muon.triggered(trgPath.c_str()));
      if (verbose && muon.triggered(trgPath.c_str())) {cout << "\t\t" << "HLT_" + tag << endl;}
    }
    (*outputVecNtuplizer)["trgMu_pt"].push_back(muon.pt());
    (*outputVecNtuplizer)["trgMu_eta"].push_back(muon.eta());
    (*outputVecNtuplizer)["trgMu_phi"].push_back(muon.phi());
    (*outputVecNtuplizer)["trgMu_charge"].push_back(muon.charge());
    if(!muon.innerTrack().isNull()) {
      auto tk = muon.innerTrack();
      auto dxy = tk->dxy(primaryVtx.position());
      (*outputVecNtuplizer)["trgMu_dz"].push_back(tk->dz(primaryVtx.position()));
      (*outputVecNtuplizer)["trgMu_dxy"].push_back(dxy);
      (*outputVecNtuplizer)["trgMu_sigdxy"].push_back(fabs(dxy)/tk->dxyError());
    }
    else {
      (*outputVecNtuplizer)["trgMu_dz"].push_back(-999);
      (*outputVecNtuplizer)["trgMu_dxy"].push_back(-999);
      (*outputVecNtuplizer)["trgMu_sigdxy"].push_back(-999);
    }
  }
  bool acceptEvent = trgMuonsMatched->size() > 0;
  (*outputNtuplizer)["N_trgMu"] = trgMuonsMatched->size();
  if(verbose) {cout << "\nTriggered muons: " << trgMuonsMatched->size() << endl;}

  (*outputNtuplizer)["N_vertexes"] = vtxHandle->size();
  if(verbose) {cout << "\nNumber of vertices: " << vtxHandle->size() << endl;}
  (*outputNtuplizer)["primaryVtx_x"] = primaryVtx.x();
  (*outputNtuplizer)["primaryVtx_y"] = primaryVtx.y();
  (*outputNtuplizer)["primaryVtx_z"] = primaryVtx.z();
  (*outputNtuplizer)["primaryVtx_sig_xx"] = primaryVtx.covariance(0, 0);
  (*outputNtuplizer)["primaryVtx_sig_xy"] = primaryVtx.covariance(0, 1);
  (*outputNtuplizer)["primaryVtx_sig_xz"] = primaryVtx.covariance(0, 2);
  (*outputNtuplizer)["primaryVtx_sig_yy"] = primaryVtx.covariance(1, 1);
  (*outputNtuplizer)["primaryVtx_sig_yz"] = primaryVtx.covariance(1, 2);
  (*outputNtuplizer)["primaryVtx_sig_zz"] = primaryVtx.covariance(2, 2);

  iEvent.put(move(trgMuonsMatched), "trgMuonsMatched");
  iEvent.put(move(outputNtuplizer), "outputNtuplizer");
  iEvent.put(move(outputVecNtuplizer), "outputVecNtuplizer");
  if (verbose) {cout << "======================== " << endl;}

  if (acceptEvent) {
    N_passed_events++;
    return true;
  }
  else return false;
}

DEFINE_FWK_MODULE(TriggerMuonsFilter);
