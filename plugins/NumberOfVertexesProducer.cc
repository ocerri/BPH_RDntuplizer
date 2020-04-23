#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include <iostream>
#include <string>
#include <regex>

using namespace std;

class NumberOfVertexesProducer : public edm::EDProducer {
   public:
      explicit NumberOfVertexesProducer(const edm::ParameterSet& iConfig);
      ~NumberOfVertexesProducer() {};

   private:
      void beginJob(const edm::EventSetup&) {};
      void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<edm::TriggerResults> triggerBitsSrc_;
      edm::EDGetTokenT <pat::PackedTriggerPrescales> triggerPrescalesSrc_;
      edm::EDGetTokenT<vector<reco::Vertex>> vtxSrc_;

      edm::Service<TFileService> fs;
      map<string, TH1D*> hNvtx;
      map<string, TH1D*> hNvtxPassed;

      vector<string> triggerTags = {"Mu12_IP6", "Mu9_IP6", "Mu7_IP4"};
      int verbose = 0;
};


NumberOfVertexesProducer::NumberOfVertexesProducer(const edm::ParameterSet& iConfig):
  triggerBitsSrc_( consumes<edm::TriggerResults> ( iConfig.getParameter<edm::InputTag>("triggerBits") ) ),
  triggerPrescalesSrc_(consumes<pat::PackedTriggerPrescales>(edm::InputTag("patTrigger"))),

  vtxSrc_( consumes<vector<reco::Vertex>> ( edm::InputTag("offlineSlimmedPrimaryVertices") ) ),
  verbose( iConfig.getParameter<int>( "verbose" ) )
{
  produces<map<string, float>>("outputNtuplizer");
  for (auto trgTag: triggerTags) {
    hNvtx[trgTag] = fs->make<TH1D>(("hNvtx"+trgTag).c_str(), ("Number of vertexes from events with active "+trgTag).c_str(), 81, -0.5, 80.5);
    hNvtxPassed[trgTag] = fs->make<TH1D>(("hNvtxPassed"+trgTag).c_str(), ("Number of vertexes from events with passed "+trgTag).c_str(), 81, -0.5, 80.5);
  }
}

void NumberOfVertexesProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBitsSrc_, triggerBits);

  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByToken(triggerPrescalesSrc_, triggerPrescales);

  edm::Handle<vector<reco::Vertex>> vtxHandle;
  iEvent.getByToken(vtxSrc_, vtxHandle);
  auto primaryVtx = (*vtxHandle)[0];

  // Output collection
  unique_ptr<map<string, float>> outputNtuplizer(new map<string, float>);

  //BPH trigger footprint
  regex txt_regex_path("HLT_Mu[0-9]+_IP[0-9]_part[0-9].*");
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  if (verbose) {cout << "\n ==== TRIGGER PATHS ==== " << endl;}

  for (auto trgTag : triggerTags){(*outputNtuplizer)["prescale_" + trgTag] = 0;}

  map<string, bool> trgActive;
  map<string, bool> trgPassed;
  for (auto trgTag : triggerTags) {
    trgActive[trgTag] = false;
    trgPassed[trgTag] = false;
  }

  for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
    auto trgName = names.triggerName(i);
    if (!regex_match(trgName, txt_regex_path)) continue;
    if (verbose) {
      cout << "Trigger " << trgName << ", prescale " << triggerPrescales->getPrescaleForIndex(i) << endl;
    }

    for (auto trgTag : triggerTags){
      bool match = trgName.substr(4, trgTag.size()) == trgTag.c_str();
      if (match && triggerPrescales->getPrescaleForIndex(i) > 0) trgActive[trgTag] = true;
      if (match && triggerBits->accept(i)) trgPassed[trgTag] = true;

      // for (int part = 0; part <= 5; part++) {
      //   regex rule(Form("HLT_%s_part%d.*", trgTag.c_str(), part));
      //   if (regex_match(trgName, rule)) {
      //     // (*outputNtuplizer)[Form("prescale_%s_part%d", trgTag.c_str(), part)] = triggerPrescales->getPrescaleForIndex(i);
      //     if (triggerPrescales->getPrescaleForIndex(i) > 0) (*outputNtuplizer)["prescale_" + trgTag]++;
      //   }
      // }

    }

  }

  auto Nvtx = vtxHandle->size();
  for (auto kv : hNvtx) {
    if(trgActive[kv.first]) {
      if (verbose) {cout << "Filling active " << kv.first << endl;}
      kv.second->Fill(Nvtx);
    }

    if(trgPassed[kv.first]) {
      if (verbose) {cout << "Filling passed " << kv.first << endl;}
      hNvtxPassed[kv.first]->Fill(Nvtx);
    }
  }

  (*outputNtuplizer)["N_vertexes"] = vtxHandle->size();

  iEvent.put(move(outputNtuplizer), "outputNtuplizer");

  if (verbose) {cout << "======================== " << endl;}
  return;
}

DEFINE_FWK_MODULE(NumberOfVertexesProducer);
