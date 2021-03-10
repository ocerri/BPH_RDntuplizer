// system include files
#include <memory>
#include <iostream>
#include <string>
#include <regex>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "VtxUtils.hh"

using namespace std;

class TriggerObjectsFilter : public edm::stream::EDFilter<> {
   public:
      explicit TriggerObjectsFilter(const edm::ParameterSet&);
      ~TriggerObjectsFilter() {
        cout << Form("Trigger filter efficiency: %d/%d = %1.2e", N_passed_events, N_analyzed_events, (double)N_passed_events/N_analyzed_events) << endl;

        fs->make<TNamed>("TriggerObjectsFilterEfficiency", Form("%d/%d", N_passed_events, N_analyzed_events));
      };

   private:
      // void beginJob(const edm::EventSetup&) {};
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      // void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;

      vector<pat::PackedCandidate> TriggerObj_matching(edm::Handle<vector<pat::PackedCandidate>>, pat::TriggerObjectStandAlone, reco::Vertex, int);

      // ----------member data ---------------------------
      edm::EDGetTokenT<edm::TriggerResults> triggerBitsSrc_;
      edm::EDGetTokenT <pat::PackedTriggerPrescales> triggerPrescalesSrc_;
      edm::EDGetTokenT<vector<pat::TriggerObjectStandAlone>> triggerObjectsSrc_;
      edm::EDGetTokenT<vector<pat::PackedCandidate>> PFCandSrc_;
      edm::EDGetTokenT<vector<reco::Vertex>> vtxSrc_;
      edm::EDGetTokenT<reco::BeamSpot> beamSpotSrc_;


      edm::Service<TFileService> fs;
      TH1I* hAllNvts;
      TH1I* hAllVtxZ;

      int N_analyzed_events = 0;
      int N_passed_events = 0;
      int trgObjectCharge_ = 0;
      int verbose = 0;
};


TriggerObjectsFilter::TriggerObjectsFilter(const edm::ParameterSet& iConfig)
{
  triggerBitsSrc_ = consumes<edm::TriggerResults> ( edm::InputTag("TriggerResults","","HLT") ) ;
  triggerPrescalesSrc_ =consumes<pat::PackedTriggerPrescales>(edm::InputTag("patTrigger"));
  triggerObjectsSrc_ =consumes<vector<pat::TriggerObjectStandAlone>> ( edm::InputTag("slimmedPatTrigger") ) ;

  PFCandSrc_ = consumes<vector<pat::PackedCandidate>>(edm::InputTag("packedPFCandidates"));
  vtxSrc_ = consumes<vector<reco::Vertex>> ( edm::InputTag("offlineSlimmedPrimaryVertices") );
  beamSpotSrc_ = consumes<reco::BeamSpot> ( edm::InputTag("offlineBeamSpot") );

  trgObjectCharge_ = iConfig.getParameter<int>( "object_charge" );
  verbose = iConfig.getParameter<int>( "verbose" );

  produces<vector<pat::PackedCandidate>>("trgCandidatesMatched");
  produces<map<string, float>>("outputNtuplizer");
  produces<map<string, vector<float>>>("outputVecNtuplizer");
  hAllNvts = fs->make<TH1I>("hAllNvts", "Number of vertexes from all the MINIAOD events", 101, -0.5, 100.5);
  hAllVtxZ = fs->make<TH1I>("hAllVtxZ", "Z coordinate of vertexes from all the MINIAOD events", 100, -25, 25);
}

bool TriggerObjectsFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  N_analyzed_events++;
  if (verbose) {cout << "Event " << N_analyzed_events << endl;}

  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBitsSrc_, triggerBits);

  edm::Handle<vector<pat::TriggerObjectStandAlone>> triggerObjects;
  iEvent.getByToken(triggerObjectsSrc_, triggerObjects);

  edm::Handle<vector<pat::PackedCandidate>> pfCandHandle;
  iEvent.getByToken(PFCandSrc_, pfCandHandle);

  edm::Handle<vector<reco::Vertex>> vtxHandle;
  iEvent.getByToken(vtxSrc_, vtxHandle);
  auto primaryVtx = (*vtxHandle)[0];
  hAllNvts->Fill((int)vtxHandle->size());
  for(auto vtx : (*vtxHandle)) hAllVtxZ->Fill(vtx.position().z());

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamSpotSrc_, beamSpotHandle);

  // Output collection
  unique_ptr<vector<pat::PackedCandidate>> trgCandidatesMatched( new vector<pat::PackedCandidate> );
  unique_ptr<map<string, float>> outputNtuplizer(new map<string, float>);
  unique_ptr<map<string, vector<float>>> outputVecNtuplizer(new map<string, vector<float>>);

  //BPH trigger footprint
  regex txt_regex_path("HLT_Mu[1279]+_IP[46]_part[0-9].*");
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

  vector<string> triggerTag = {"Mu12_IP6", "Mu7_IP4", "Mu9_IP6"};
  map<string, bool> trigger_bits;
  for(auto n : triggerTag) trigger_bits[n] = false;

  bool acceptEvent = false;

  if (verbose) {cout << "\n == BPH TRIGGER OBJ == " << endl;}
  for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
    regex txt_regex_coll("hlt.*MuonCandidates::HLT");
    if (!regex_match(obj.collection(), txt_regex_coll)) continue;
    if(trgObjectCharge_ != 0 && obj.charge() != trgObjectCharge_) continue;

    obj.unpackNamesAndLabels(iEvent, *triggerBits);
    vector pathNamesLast = obj.pathNames(true);

    map<string, bool> triggersAux;
    for(auto n : triggerTag) triggersAux[n] = false;
    bool acceptObject = false;
    for (unsigned h = 0, n = pathNamesLast.size(); h < n; ++h) {
      if (!regex_match(pathNamesLast[h], txt_regex_path)) continue;
      for(auto n : triggerTag) {
        TString s_aux = pathNamesLast[h];
        if ( s_aux.BeginsWith("HLT_" + n) ) {
          acceptEvent = true;
          trigger_bits[n] = true;
          acceptObject = true;
          triggersAux[n] = true;
        }
      }
    }

    if (acceptObject){
      if (verbose) {
        cout << "\t\t" << obj.charge() << endl;
        cout << "\tTrigger object: charge" << obj.charge() << ",  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << endl;
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


      for(auto tag : triggerTag) {
        (*outputVecNtuplizer)["trgObj_HLT_" + tag].push_back(triggersAux[tag]);
        if (verbose && triggersAux[tag]) {cout << "\t\t\t" << "HLT_" + tag << endl;}
      }

      auto matching_cand = TriggerObj_matching(pfCandHandle, obj, primaryVtx, trgObjectCharge_);
      if (matching_cand.size()>0) {
        if (verbose) {
          cout << "Trigger object matched with " << matching_cand.size() << " PF candidates" << endl;
          cout << "pdgID:";
          for (auto c : matching_cand) {
            cout << " " << c.pdgId();
          }
          cout << endl;
        }
        auto cand = matching_cand[0];
        trgCandidatesMatched->push_back(cand);

        (*outputVecNtuplizer)["trgCand_pt"].push_back(cand.pt());
        (*outputVecNtuplizer)["trgCand_eta"].push_back(cand.eta());
        (*outputVecNtuplizer)["trgCand_phi"].push_back(cand.phi());
        (*outputVecNtuplizer)["trgCand_charge"].push_back(cand.charge());
        (*outputVecNtuplizer)["trgCand_pdgId"].push_back(cand.pdgId());
        (*outputVecNtuplizer)["trgCand_isCaloMuon"].push_back(cand.isCaloMuon());
        (*outputVecNtuplizer)["trgCand_isConvertedPhoton"].push_back(cand.isConvertedPhoton());
        (*outputVecNtuplizer)["trgCand_isElectron"].push_back(cand.isElectron());
        (*outputVecNtuplizer)["trgCand_isGlobalMuon"].push_back(cand.isGlobalMuon());
        (*outputVecNtuplizer)["trgCand_isGoodEgamma"].push_back(cand.isGoodEgamma());
        (*outputVecNtuplizer)["trgCand_isIsolatedChargedHadron"].push_back(cand.isIsolatedChargedHadron());
        (*outputVecNtuplizer)["trgCand_isMuon"].push_back(cand.isMuon());
        (*outputVecNtuplizer)["trgCand_isPhoton"].push_back(cand.isPhoton());
        (*outputVecNtuplizer)["trgCand_isStandAloneMuon"].push_back(cand.isStandAloneMuon());
        (*outputVecNtuplizer)["trgCand_isTrackerMuon"].push_back(cand.isTrackerMuon());

        if(cand.hasTrackDetails()) {
          auto tk = cand.bestTrack();
          (*outputVecNtuplizer)["trgMu_dz"].push_back(tk->dz(primaryVtx.position()));
          auto dxyUnc = tk->dxyError();

          auto dxy_PV = fabs(tk->dxy(primaryVtx.position()));
          (*outputVecNtuplizer)["trgMu_dxy_PV"].push_back(dxy_PV);
          (*outputVecNtuplizer)["trgMu_sigdxy_PV"].push_back( dxy_PV/dxyUnc );

          auto dxy_BS = fabs(tk->dxy((*beamSpotHandle)));
          (*outputVecNtuplizer)["trgMu_dxy_BS"].push_back(dxy_BS);
          (*outputVecNtuplizer)["trgMu_sigdxy_BS"].push_back( dxy_BS/dxyUnc );
        }
        else {
        (*outputVecNtuplizer)["trgMu_dz"].push_back(-9999);
        (*outputVecNtuplizer)["trgMu_dxy_PV"].push_back(-9999);
        (*outputVecNtuplizer)["trgMu_sigdxy_PV"].push_back(-9999);
        (*outputVecNtuplizer)["trgMu_dxy_BS"].push_back(-9999);
        (*outputVecNtuplizer)["trgMu_sigdxy_BS"].push_back(-9999);
      }
      }
      else {
        (*outputVecNtuplizer)["trgCand_pt"].push_back(-1);
        (*outputVecNtuplizer)["trgCand_eta"].push_back(-9999);
        (*outputVecNtuplizer)["trgCand_phi"].push_back(-9999);
        (*outputVecNtuplizer)["trgCand_charge"].push_back(0);
        (*outputVecNtuplizer)["trgCand_pdgId"].push_back(0);
        (*outputVecNtuplizer)["trgCand_isCaloMuon"].push_back(0);
        (*outputVecNtuplizer)["trgCand_isConvertedPhoton"].push_back(0);
        (*outputVecNtuplizer)["trgCand_isElectron"].push_back(0);
        (*outputVecNtuplizer)["trgCand_isGlobalMuon"].push_back(0);
        (*outputVecNtuplizer)["trgCand_isGoodEgamma"].push_back(0);
        (*outputVecNtuplizer)["trgCand_isIsolatedChargedHadron"].push_back(0);
        (*outputVecNtuplizer)["trgCand_isMuon"].push_back(0);
        (*outputVecNtuplizer)["trgCand_isPhoton"].push_back(0);
        (*outputVecNtuplizer)["trgCand_isStandAloneMuon"].push_back(0);
        (*outputVecNtuplizer)["trgCand_isTrackerMuon"].push_back(0);
        (*outputVecNtuplizer)["trgMu_dz"].push_back(-9999);
        (*outputVecNtuplizer)["trgMu_dxy_PV"].push_back(-9999);
        (*outputVecNtuplizer)["trgMu_sigdxy_PV"].push_back(-9999);
        (*outputVecNtuplizer)["trgMu_dxy_BS"].push_back(-9999);
        (*outputVecNtuplizer)["trgMu_sigdxy_BS"].push_back(-9999);
      }
    }
  }

  for(auto tag : triggerTag) {
    (*outputNtuplizer)["trgPassed_HLT_" + tag] = trigger_bits[tag];
    if (verbose && trigger_bits[tag]) {cout << "Event passed HLT_" + tag << endl;}
  }

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

  iEvent.put(move(trgCandidatesMatched), "trgCandidatesMatched");
  iEvent.put(move(outputNtuplizer), "outputNtuplizer");
  iEvent.put(move(outputVecNtuplizer), "outputVecNtuplizer");
  if (verbose) {cout << "======================== " << endl;}

  if (acceptEvent) {
    N_passed_events++;
    return true;
  }
  else return false;
}

vector<pat::PackedCandidate> TriggerObjectsFilter::TriggerObj_matching(edm::Handle<vector<pat::PackedCandidate>> candidates_list,
                                                              pat::TriggerObjectStandAlone obj,
                                                              reco::Vertex vtx,
                                                              int charge=0)
{
  double max_DeltaR = 0.01;
  double max_Delta_pt_rel = 0.05;

  double bestM_DeltaR = max_DeltaR;

  vector<pat::PackedCandidate> out;

  for( auto cand : *candidates_list) {
    if(charge && charge!=cand.charge()) continue;

    double deltaR = vtxu::dR(cand.phi(), obj.phi(), cand.eta(), obj.eta());
    double dpt_rel = abs(cand.pt() - obj.pt())/obj.pt();
    if (dpt_rel < max_Delta_pt_rel && deltaR < max_DeltaR) {
      if(verbose) {cout << "\t\tMatched with deltaR=" << deltaR << " and dpt_rel=" << dpt_rel << endl;}
      if (deltaR <= bestM_DeltaR) {
        bestM_DeltaR = deltaR;
        out.insert(out.begin(), cand);
      }
      else out.push_back(cand);
    }
  }
  return out;
}

DEFINE_FWK_MODULE(TriggerObjectsFilter);
