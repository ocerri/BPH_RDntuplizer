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
#include "DataFormats/VertexReco/interface/Vertex.h"

#include <iostream>
#include <string>
#include <regex>

using namespace std;

class BPHTriggerPathProducer : public edm::EDProducer {
   public:
      explicit BPHTriggerPathProducer(const edm::ParameterSet& iConfig);
      ~BPHTriggerPathProducer() {};

   private:
      void beginJob(const edm::EventSetup&) {};
      void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;

      vector<pat::Muon> TriggerObj_matching(edm::Handle<vector<pat::Muon>>, pat::TriggerObjectStandAlone, edm::Handle<std::vector<reco::Vertex>>, int);


      // ----------member data ---------------------------
      edm::EDGetTokenT<edm::TriggerResults> triggerBitsSrc_;
      int muonCharge_ = 0;
      edm::EDGetTokenT<vector<pat::TriggerObjectStandAlone>> triggerObjectsSrc_;
      edm::EDGetTokenT<vector<pat::Muon>> muonSrc_;
      edm::EDGetTokenT<vector<reco::Vertex>> vtxSrc_;
      int verbose = 0;
};


BPHTriggerPathProducer::BPHTriggerPathProducer(const edm::ParameterSet& iConfig):
  triggerBitsSrc_( consumes<edm::TriggerResults> ( iConfig.getParameter<edm::InputTag>("triggerBits") ) ),
  muonCharge_( iConfig.getParameter<int>( "muon_charge" ) ),
  triggerObjectsSrc_(consumes<vector<pat::TriggerObjectStandAlone>> ( iConfig.getParameter<edm::InputTag>("triggerObjects") ) ),
  muonSrc_( consumes<vector<pat::Muon>> ( iConfig.getParameter<edm::InputTag>( "muonCollection" ) ) ),
  vtxSrc_( consumes<vector<reco::Vertex>> ( iConfig.getParameter<edm::InputTag>( "vertexCollection" ) ) ),
  verbose( iConfig.getParameter<int>( "verbose" ) )
{
  produces<vector<pat::Muon>>("trgMuonsMatched");
  produces<map<string, float>>("outputNtuplizer");
}

void BPHTriggerPathProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBitsSrc_, triggerBits);

  edm::Handle<vector<pat::TriggerObjectStandAlone>> triggerObjects;
  iEvent.getByToken(triggerObjectsSrc_, triggerObjects);

  edm::Handle<std::vector<pat::Muon>> muonHandle;
  iEvent.getByToken(muonSrc_, muonHandle);
  unsigned int muonNumber = muonHandle->size();

  edm::Handle<std::vector<reco::Vertex>> vtxHandle;
  iEvent.getByToken(vtxSrc_, vtxHandle);

  // Output collection
  unique_ptr<vector<pat::Muon>> trgMuonsMatched( new vector<pat::Muon> );

  //BPH trigger footprint
  std::regex txt_regex_path("HLT_Mu[0-9]+_IP[0-9]+.*");
  if (verbose) {std::cout << "\n == BPH TRIGGER OBJ == " << std::endl;}

  for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
    obj.unpackNamesAndLabels(iEvent, *triggerBits);
    vector pathNamesLast = obj.pathNames(true);
    unsigned int obj_BPH_path = 0;
    for (unsigned h = 0, n = pathNamesLast.size(); h < n; ++h) {
      if (regex_match(pathNamesLast[h], txt_regex_path)) obj_BPH_path++;
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
        cout << "\tBPH Paths:"<< endl;
        for (unsigned h = 0, n = pathNamesLast.size(); h < n; ++h) {
          if (regex_match(pathNamesLast[h], txt_regex_path)) {
            cout << "\t\t" << pathNamesLast[h] << endl;
          }
        }
      }

      auto matching_muons = TriggerObj_matching(muonHandle, obj, vtxHandle, muonCharge_);
      if (matching_muons.size()>0) trgMuonsMatched->push_back(matching_muons[0]);

    }
  }

  if (verbose) {
    cout << "\n MUONS LIST" << endl;
    for (unsigned int k = 0; k < muonNumber; ++k) {
      const pat::Muon & muon = (*muonHandle)[k];
      bool isTightMuon = false;
      for (auto vtx : *vtxHandle) {
        if(muon.isTightMuon(vtx)) isTightMuon = true;
      }
      if (isTightMuon) {
        cout << "\t" << " " << muon.pdgId() << "  " << muon.pt() << "  " << muon.eta() << "  " << muon.phi() << endl;
      }
    }
  }

  unique_ptr<map<string, float>> outputNtuplizer(new map<string, float>);
  if (trgMuonsMatched->size()) {
    auto m = (*trgMuonsMatched)[0];
    (*outputNtuplizer)["trgMu_pt"] = m.pt();
    (*outputNtuplizer)["trgMu_eta"] = m.eta();
    (*outputNtuplizer)["trgMu_phi"] = m.phi();
    (*outputNtuplizer)["trgMu_charge"] = m.charge();
  }
  else {
    (*outputNtuplizer)["trgMu_pt"] = -1;
    (*outputNtuplizer)["trgMu_eta"] = 0;
    (*outputNtuplizer)["trgMu_phi"] = 0;
    (*outputNtuplizer)["trgMu_charge"] = 0;
  }

  iEvent.put(move(trgMuonsMatched), "trgMuonsMatched");
  iEvent.put(move(outputNtuplizer), "outputNtuplizer");

  if (verbose) {std::cout << "======================== " << std::endl;}
  return;
}



vector<pat::Muon> BPHTriggerPathProducer::TriggerObj_matching(edm::Handle<vector<pat::Muon>> muon_list,
                                                              pat::TriggerObjectStandAlone obj,
                                                              edm::Handle<std::vector<reco::Vertex>> vtx_list,
                                                              int mu_charge=0)
{
  double max_DeltaR = 0.005;
  double max_Delta_pt_rel = 0.01;

  double bestM_DeltaR = max_DeltaR;

  vector<pat::Muon> out;

  for( auto muon : *muon_list) {
    if(mu_charge && mu_charge!=muon.charge()) continue;

    double dEta = muon.eta() - obj.eta();
    double dPhi = muon.phi() - obj.phi();
    double pi = 3.14159265358979323846;
    while (fabs(dPhi) > pi) {
      int sgn = dPhi > 0? 1 : -1;
      dPhi -= sgn*2*pi;
    }

    double deltaR = sqrt(dEta*dEta + dPhi*dPhi);

    double dpt_rel = abs(muon.pt() - obj.pt())/obj.pt();

    if (dpt_rel < max_Delta_pt_rel && deltaR < max_DeltaR) {
      bool isTightMuon = false;
      for (auto vtx : *vtx_list) {
        if(muon.isTightMuon(vtx)) isTightMuon = true;
      }
      if (isTightMuon) {
        if(verbose) {cout << "\t\tMuon matched with deltaR=" << deltaR << " and dpt_rel=" << dpt_rel << endl;}
        if (deltaR <= bestM_DeltaR) {
          bestM_DeltaR = deltaR;
          out.insert(out.begin(), muon);
        }
        else out.push_back(muon);
      }
    }
  }

  return out;
}

DEFINE_FWK_MODULE(BPHTriggerPathProducer);
