#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"

#include <iostream>
#include <string>


using namespace std;

class FlatTreeWriter : public edm::EDAnalyzer {
public:
  explicit FlatTreeWriter( const edm::ParameterSet & );
  ~FlatTreeWriter() {};

private:
  virtual void analyze( const edm::Event& , const edm::EventSetup& ) override;

  // edm::InputTag src_;
  string cmssw;
  // vector<string> output_soruce_modules;

  TTree* tree;
  map<string, float> out_map;
  map<string, vector<float>> outv_map;

  bool first_call = true;

  int verbose = 0;
};

FlatTreeWriter::FlatTreeWriter( const edm::ParameterSet & cfg ) :
  cmssw( cfg.getParameter<string>("cmssw") )
  // output_soruce_modules( cfg.getParameter<vector<string>>("output_soruce_modules") )
{
  verbose = cfg.getParameter<int>( "verbose" );

  edm::Service<TFileService> fs;
  tree = fs->make<TTree>( "Tevts", "Events Tree from BPHRDntuplizer");

  fs->make<TNamed>("cmssw",cmssw.c_str() );

  consumesMany<map<string, float>>();
  consumesMany<map<string, vector<float>>>();
}

void FlatTreeWriter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  vector< edm::Handle<map<string, float>> > outMapHandle;
  iEvent.getManyByType(outMapHandle);
  for( auto h : outMapHandle ) {
    for( auto& kv : *(h.product()) ) {
      out_map[kv.first] = kv.second;
      if (verbose) {cout << Form("%s: %.2f", kv.first.c_str(), kv.second) << endl;}
    }
  }

  vector< edm::Handle<map<string, vector<float>>> > outVMapHandle;
  iEvent.getManyByType(outVMapHandle);
  for( auto h : outVMapHandle ) {
    for( auto& kv : *(h.product()) ) {
      outv_map[kv.first] = kv.second;
      if (verbose) {cout << kv.first.c_str() << ": "<< kv.second.size() << endl;}
    }
  }


  if (first_call) {
    // Creating the branches in the output tree
    for(auto& kv : out_map) {
      auto k = kv.first;
      tree->Branch(k.c_str(), &(out_map[k]));
    }
    first_call = false;
  }

  tree->Fill();
}

DEFINE_FWK_MODULE(FlatTreeWriter);
