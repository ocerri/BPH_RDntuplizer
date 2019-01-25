#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include <iostream>
#include <string>
#include <regex>

using namespace std;

class BPHTriggerPathFilter : public edm::EDFilter {
   public:
      explicit BPHTriggerPathFilter(const edm::ParameterSet& iConfig);
      ~BPHTriggerPathFilter();

   private:
      bool filter(edm::Event& iEvent, const edm::EventSetup& iSetup) override;

      vector<pat::Muon> TriggerObj_matching(edm::Handle<vector<pat::Muon>>, pat::TriggerObjectStandAlone);


      // ----------member data ---------------------------
      edm::EDGetTokenT<edm::TriggerResults> triggerBitsSrc_;
      int muon_charge_ = 0;
      edm::EDGetTokenT<vector<pat::TriggerObjectStandAlone>> triggerObjectsSrc_;
      edm::EDGetTokenT<std::vector<pat::Muon>> muonSrc_;

      int verbose = 1;
};


BPHTriggerPathFilter::BPHTriggerPathFilter(const edm::ParameterSet& iConfig):
  triggerBitsSrc_( consumes<edm::TriggerResults> ( iConfig.getParameter<edm::InputTag>("triggerBits") ) ),
  muon_charge_( iConfig.getParameter<int>( "muon_charge" ) ),
  triggerObjectsSrc_(consumes<vector<pat::TriggerObjectStandAlone>> ( iConfig.getParameter<edm::InputTag>("triggerObjects") ) ),
  muonSrc_( consumes<std::vector<pat::Muon>> ( iConfig.getParameter<edm::InputTag>( "muonCollection" ) ) )
{
  produces<std::vector<pat::Muon>>("TrgMu").setBranchAlias( "TrgMu");
}

BPHTriggerPathFilter::~BPHTriggerPathFilter(){};

bool BPHTriggerPathFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBitsSrc_, triggerBits);

  edm::Handle<vector<pat::TriggerObjectStandAlone>> triggerObjects;
  iEvent.getByToken(triggerObjectsSrc_, triggerObjects);

  edm::Handle<std::vector<pat::Muon>> muonHandle;
  iEvent.getByToken(muonSrc_, muonHandle);
  unsigned int muonNumber = muonHandle->size();

  // Output collection
  std::unique_ptr<std::vector<pat::Muon>> result( new std::vector<pat::Muon> );

  //BPH trigger footprint
  std::regex txt_regex_path("HLT_Mu[0-9]+_IP[0-9]+.*");
  if (verbose) {std::cout << "\n == BPH TRIGGER OBJ == " << std::endl;}

  for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
    if(muon_charge_ && obj.charge() != muon_charge_) continue;

    obj.unpackNamesAndLabels(iEvent, *triggerBits);
    vector pathNamesAll = obj.pathNames(false);
    vector pathNamesLast = obj.pathNames(true);

    unsigned int obj_BPH_path = 0;
    for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
      if (regex_match(pathNamesAll[h], txt_regex_path)) obj_BPH_path++;
    }

    std::regex txt_regex_coll("hlt.*MuonCandidates::HLT");
    bool HLT_muon = regex_match(obj.collection(), txt_regex_coll);

    if (obj_BPH_path>0 && HLT_muon){
      if (verbose) {
        cout << "\t\t" << obj.charge() << endl;
        cout << "\tTriggered mu" << obj.charge() << ":  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << endl;
        // Print trigger object collection and type
        cout << "\tCollection: " << obj.collection() << endl;

        // Print associated trigger paths
        cout << "\tPaths (" << pathNamesAll.size()<<"/"<<pathNamesLast.size()<<"):"<< endl;
        for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
          cout << "\t\t" << pathNamesAll[h] << endl;
        }
      }

      auto matching_muons = TriggerObj_matching(muonHandle, obj);
      if (matching_muons.size()>0) result->push_back(matching_muons[0]);

    }
  }

  if (verbose) {
    cout << "\n MUONS LIST" << endl;
    for (unsigned int k = 0; k < muonNumber; ++k) {
      const pat::Muon & muon = (*muonHandle)[k];
      cout << "\t" << " " << muon.pdgId() << "  " << muon.pt() << "  " << muon.eta() << "  " << muon.phi() << endl;
    }
  }

  iEvent.put(std::move(result));

  if (verbose) {std::cout << "======================== " << std::endl;}
  return result->size()>0;
}



vector<pat::Muon> BPHTriggerPathFilter::TriggerObj_matching(edm::Handle<vector<pat::Muon>> muon_list, pat::TriggerObjectStandAlone obj) {
  double max_DeltaR = 0.005;
  double max_Delta_pt_rel = 0.01;

  double bestM_DeltaR = max_DeltaR;

  vector<pat::Muon> out;

  for( auto muon : *muon_list) {
    double dEta = muon.eta() - obj.eta();
    double dPhi = muon.phi() - obj.phi();
    double deltaR = sqrt(dEta*dEta + dPhi*dPhi);

    double dpt_rel = abs(muon.pt() - obj.pt())/obj.pt();

    if (dpt_rel < max_Delta_pt_rel && deltaR < max_DeltaR) {
      if(verbose) {cout << "\t\tMuon matched with deltaR=" << deltaR << " and dpt_rel=" << dpt_rel << endl;}
      if (deltaR <= bestM_DeltaR) {
        bestM_DeltaR = deltaR;
        out.insert(out.begin(), muon);
      }
      else out.push_back(muon);
    }
  }

  return out;
}

DEFINE_FWK_MODULE(BPHTriggerPathFilter);
