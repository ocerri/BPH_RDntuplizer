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
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

// Needed for Transient Tracks
// #include "TrackingTools/TransientTrack/interface/TransientTrack.h"
// #include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
// #include "TrackingTools/Records/interface/TransientTrackRecord.h"
// #include "MagneticField/Engine/interface/MagneticField.h"

#include "VtxUtils.hh"

#include <iostream>
#include <string>
#include <regex>

using namespace std;

class TagAndProbeProducer : public edm::EDProducer {
   public:
      explicit TagAndProbeProducer(const edm::ParameterSet& iConfig);
      ~TagAndProbeProducer() {};

   private:
      void beginJob(const edm::EventSetup&) {};
      void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;

      vector<pat::Muon> TriggerObj_matching(edm::Handle<vector<pat::Muon>>, pat::TriggerObjectStandAlone, reco::Vertex, int);


      // ----------member data ---------------------------
      edm::EDGetTokenT<edm::TriggerResults> triggerBitsSrc_;
      edm::EDGetTokenT <pat::PackedTriggerPrescales> triggerPrescalesSrc_;
      edm::EDGetTokenT<vector<pat::TriggerObjectStandAlone>> triggerObjectsSrc_;

      edm::EDGetTokenT<vector<pat::Muon>> muonSrc_;
      edm::EDGetTokenT<vector<reco::Vertex>> vtxSrc_;
      int muonCharge_ = 0;
      int verbose = 0;
};


TagAndProbeProducer::TagAndProbeProducer(const edm::ParameterSet& iConfig):
  triggerBitsSrc_( consumes<edm::TriggerResults> ( edm::InputTag("TriggerResults","","HLT") ) ),
  triggerPrescalesSrc_(consumes<pat::PackedTriggerPrescales>(edm::InputTag("patTrigger"))),
  triggerObjectsSrc_(consumes<vector<pat::TriggerObjectStandAlone>> ( edm::InputTag("slimmedPatTrigger") ) ),
  muonSrc_( consumes<vector<pat::Muon>> ( edm::InputTag("slimmedMuons") ) ),
  vtxSrc_( consumes<vector<reco::Vertex>> ( edm::InputTag("offlineSlimmedPrimaryVertices") ) ),
  muonCharge_( iConfig.getParameter<int>( "muon_charge" ) ),
  verbose( iConfig.getParameter<int>( "verbose" ) )
{
  produces<vector<pat::Muon>>("trgMuonsMatched");
  produces<map<string, float>>("outputNtuplizer");
}

void TagAndProbeProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBitsSrc_, triggerBits);

  edm::Handle<vector<pat::TriggerObjectStandAlone>> triggerObjects;
  iEvent.getByToken(triggerObjectsSrc_, triggerObjects);

  edm::Handle<vector<pat::Muon>> muonHandle;
  iEvent.getByToken(muonSrc_, muonHandle);

  edm::Handle<vector<reco::Vertex>> vtxHandle;
  iEvent.getByToken(vtxSrc_, vtxHandle);
  auto primaryVtx = (*vtxHandle)[0];

  // Output collection
  unique_ptr<vector<pat::Muon>> trgMuonsMatched( new vector<pat::Muon> );

  //BPH trigger footprint
  regex txt_regex_path("HLT_Mu[0-9]+_IP[0-9]_part[0-9].*");

  if (verbose) {
    edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
    iEvent.getByToken(triggerPrescalesSrc_, triggerPrescales);
    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
    cout << "\n == TRIGGER PATHS= " << endl;
    for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
      auto trgName = names.triggerName(i);
      if (!regex_match(trgName, txt_regex_path)) continue;
      cout << "Trigger " << trgName << ", prescale " << triggerPrescales->getPrescaleForIndex(i);
      cout << ": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)") << endl;
    }
  }

  vector<string> trigger_paths = {"Mu12_IP6", "Mu9_IP5", "Mu7_IP4", "Mu9_IP4", "Mu8_IP5", "Mu8_IP6", "Mu9_IP6", "Mu8_IP3"};
  map<string, bool> trigger_bits;
  for(auto n : trigger_paths) trigger_bits[n] = false;

  if (verbose) {cout << "\n == BPH TRIGGER OBJ == " << endl;}
  for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
    regex txt_regex_coll("hlt.*MuonCandidates::HLT");
    if (!regex_match(obj.collection(), txt_regex_coll)) continue;

    obj.unpackNamesAndLabels(iEvent, *triggerBits);
    vector pathNamesLast = obj.pathNames(true);
    unsigned int obj_BPH_path = 0;
    for (unsigned h = 0, n = pathNamesLast.size(); h < n; ++h) {
      if (regex_match(pathNamesLast[h], txt_regex_path)){
        obj_BPH_path++;
        if (trgMuonsMatched->size() == 0) {
          for(auto n : trigger_paths) {
            TString s_aux = pathNamesLast[h];
            if ( s_aux.BeginsWith("HLT_" + n) ) trigger_bits[n] = true;
          }
        }
      }
    }

    if (obj_BPH_path>0){
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

      auto matching_muons = TriggerObj_matching(muonHandle, obj, primaryVtx, muonCharge_);
      if (matching_muons.size()>0) trgMuonsMatched->push_back(matching_muons[0]);

    }
  }

  if (verbose) {
    cout << "\n MUONS LIST" << endl;
    for (auto muon : (*muonHandle)) {
      cout << "\t" << Form("id:%d  pt:%.1f  eta:%.1f  phi:%.1f softID:%d trgBPH:%d", muon.pdgId(), muon.pt(), muon.eta(), muon.phi(), muon.isSoftMuon(primaryVtx), muon.triggered("HLT_Mu*_IP*_part*_v*")) << endl;
    }
  }

  unique_ptr<map<string, float>> outputNtuplizer(new map<string, float>);
  (*outputNtuplizer)["N_trgMu"] = trgMuonsMatched->size();
  if (trgMuonsMatched->size()) {
    for(auto n : trigger_paths) {
      (*outputNtuplizer)["trgMu_HLT_" + n] = trigger_bits[n];
    }
    auto m = (*trgMuonsMatched)[0];
    (*outputNtuplizer)["trgMu_pt"] = m.pt();
    (*outputNtuplizer)["trgMu_eta"] = m.eta();
    (*outputNtuplizer)["trgMu_phi"] = m.phi();
    (*outputNtuplizer)["trgMu_charge"] = m.charge();
    auto tk = m.innerTrack();
    auto dxy = tk->dxy(primaryVtx.position());
    (*outputNtuplizer)["trgMu_dz"] = tk->dz(primaryVtx.position());
    (*outputNtuplizer)["trgMu_dxy"] = dxy;
    (*outputNtuplizer)["trgMu_sigdxy"] = fabs(dxy)/tk->dxyError();

    if(verbose) {
      cout << "\nTriggered muons: " << trgMuonsMatched->size() << endl;
      cout << Form("Muon dxy: %.4f +/- %.4f (sig = %.1f)", dxy, tk->dxyError(), fabs(dxy)/tk->dxyError());
      cout << Form(", dz: %.3f", tk->dz(primaryVtx.position())) << endl;
    }
  }


  (*outputNtuplizer)["N_vertexes"] = vtxHandle->size();
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

  if (verbose) {cout << "======================== " << endl;}
  return;
}



vector<pat::Muon> TagAndProbeProducer::TriggerObj_matching(edm::Handle<vector<pat::Muon>> muon_list,
                                                              pat::TriggerObjectStandAlone obj,
                                                              reco::Vertex vtx,
                                                              int mu_charge=0)
{
  double max_DeltaR = 0.01;
  double max_Delta_pt_rel = 0.05;

  double bestM_DeltaR = max_DeltaR;

  vector<pat::Muon> out;

  for( auto muon : *muon_list) {
    if(mu_charge && mu_charge!=muon.charge()) continue;
    if(muon.innerTrack().isNull()) continue;
    if (!muon.isSoftMuon(vtx)) continue;

    double deltaR = vtxu::dR(muon.phi(), obj.phi(), muon.eta(), obj.eta());
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

DEFINE_FWK_MODULE(TagAndProbeProducer);
