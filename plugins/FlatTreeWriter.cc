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
  string cfg_name = "";
  string commit_hash = "";
  // vector<string> output_soruce_modules;

  TTree* tree;
  map<string, float> out_map;
  map<string, vector<float>> outv_map;

  bool first_call = true;

  int verbose = 0;
};

FlatTreeWriter::FlatTreeWriter( const edm::ParameterSet & cfg ) :
  cmssw( cfg.getParameter<string>("cmssw") ),
  cfg_name( cfg.getParameter<string>("cfg_name") ),
  commit_hash( cfg.getParameter<string>("commit_hash") )
  // output_soruce_modules( cfg.getParameter<vector<string>>("output_soruce_modules") )
{
  verbose = cfg.getParameter<int>( "verbose" );

  edm::Service<TFileService> fs;
  tree = fs->make<TTree>( "Tevts", "Events Tree from BPHRDntuplizer");

  fs->make<TNamed>("cmssw", cmssw.c_str() );
  fs->make<TNamed>("cfgName", cfg_name.c_str() );
  fs->make<TNamed>("commitHash", commit_hash.c_str() );

  consumesMany<map<string, float>>();
  consumesMany<map<string, vector<float>>>();
}

void FlatTreeWriter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  cout << "--------------------- Tree Wirter -----------------------" << endl;

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
      if (verbose) {
        cout << kv.first.c_str() << " (" << kv.second.size() << "):";
        for(float val : kv.second) {cout << Form(" %.2f", val);}
        cout << endl;
      }
    }
  }


  if (first_call) {
    if(verbose) {cout << "\nCreating the branches in the output tree:\n";}
    for(auto& kv : out_map) {
      auto k = kv.first;
      if(verbose) {cout << k << endl;}
      tree->Branch(k.c_str(), &(out_map[k]));
    }

    for(auto& kv : outv_map) {
      auto k = kv.first;
      if(verbose) {cout << k << endl;}
      tree->Branch(k.c_str(), &(outv_map[k]));
    }
    first_call = false;
    if(verbose) {cout << '\n';}
  }

  tree->Fill();

  if(verbose) {
    cout << "--------------------------------------------------------" << endl;
    cout << "--------------------- Event over -----------------------" << endl;
    cout << "--------------------------------------------------------" << endl;
    cout << endl << endl << endl;
  }
}

DEFINE_FWK_MODULE(FlatTreeWriter);
