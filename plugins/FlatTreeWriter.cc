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
  string cmssw_;
  TTree* tree;
  map<string, float> out_map;
};

FlatTreeWriter::FlatTreeWriter( const edm::ParameterSet & cfg ) :
  // src_( cfg.getUntrackedParameter<edm::InputTag>( "src" ) ),
  cmssw_( cfg.getParameter<string>("cmssw") )
{
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>( "Tevts", "Events Tree from BPHRDntuplizer");

  fs->make<TNamed>("cmssw",cmssw_.c_str() );
}

void FlatTreeWriter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  out_map["test"] = 7;

  if (tree->GetNbranches() == 0) {
    // Creating the branches in the output tree
    for(auto& kv : out_map) {
      auto k = kv.first;
      tree->Branch(k.c_str(), &(out_map[k]));
    }
  }

  tree->Fill();
}

DEFINE_FWK_MODULE(FlatTreeWriter);
