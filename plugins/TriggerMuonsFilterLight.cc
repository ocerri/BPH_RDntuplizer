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
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

// L1 trigger
#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "VtxUtils.hh"

// Pileup info
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

using namespace std;

class TriggerMuonsFilterLight : public edm::stream::EDFilter<> {
   public:
      explicit TriggerMuonsFilterLight(const edm::ParameterSet&);
      ~TriggerMuonsFilterLight() {
        cout << Form("Muons trigger filter efficiency: %d/%d = %1.2e", N_passed_events, N_analyzed_events, (double)N_passed_events/N_analyzed_events) << endl;

        fs->make<TNamed>("TriggerMuonsFilterLightEfficiency", Form("%d/%d", N_passed_events, N_analyzed_events));
      };

   private:
      // void beginJob(const edm::EventSetup&) {};
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      tuple<uint, float, float, float> matchL1Muon(pat::Muon muReco, BXVector<l1t::Muon> muonsL1, uint skipIdx=9999);

      // ----------member data ---------------------------
      edm::EDGetTokenT<BXVector<l1t::Muon>> l1MuonSrc_;
      edm::EDGetTokenT<vector<pat::Muon>> muonSrc_;
      edm::EDGetTokenT<vector<reco::Vertex>> vtxSrc_;
      edm::EDGetTokenT<reco::BeamSpot> beamSpotSrc_;

      edm::EDGetTokenT<vector<PileupSummaryInfo>> pileupMCSrc_;

      edm::Service<TFileService> fs;
      TH1I* hAllNvts;
      TH1I* hAllNTrueIntMC;
      TH1I* hAllVtxZ;

      int N_analyzed_events = 0;
      int N_passed_events = 0;
      int muonCharge_ = 0;
      bool isMC_ = 0;
      int verbose = 0;
};


TriggerMuonsFilterLight::TriggerMuonsFilterLight(const edm::ParameterSet& iConfig):
  l1MuonSrc_( consumes<BXVector<l1t::Muon>> ( edm::InputTag("gmtStage2Digis", "Muon", "RECO") ) ),
  muonSrc_( consumes<vector<pat::Muon>> ( edm::InputTag("slimmedMuons") ) ),
  vtxSrc_( consumes<vector<reco::Vertex>> ( edm::InputTag("offlineSlimmedPrimaryVertices") ) ),
  beamSpotSrc_( consumes<reco::BeamSpot> ( edm::InputTag("offlineBeamSpot") ) ),
  muonCharge_( iConfig.getParameter<int>( "muon_charge" ) ),
  isMC_( iConfig.getParameter<int>( "isMC" ) ),
  verbose( iConfig.getParameter<int>( "verbose" ) )
{
  produces<vector<pat::Muon>>("trgMuonsMatched");
  produces<map<string, float>>("outputNtuplizer");
  produces<map<string, vector<float>>>("outputVecNtuplizer");
  hAllNvts = fs->make<TH1I>("hAllNvts", "Number of vertexes from all the MINIAOD events", 101, -0.5, 100.5);
  hAllVtxZ = fs->make<TH1I>("hAllVtxZ", "Z coordinate of vertexes from all the MINIAOD events", 100, -25, 25);
  if(isMC_) {
    pileupMCSrc_ = consumes<vector<PileupSummaryInfo>> ( edm::InputTag("slimmedAddPileupInfo") );
    hAllNTrueIntMC = fs->make<TH1I>("hAllNTrueIntMC", "Number of true interactions generated in MC", 101, -0.5, 100.5);
  }

  std::srand(std::time(nullptr));
}

bool TriggerMuonsFilterLight::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Random dropout
  if (std::rand()%10 != 0) return false;

  N_analyzed_events++;
  if (verbose) {cout << "Event " << N_analyzed_events << endl;}

  edm::Handle<BXVector<l1t::Muon>> l1MuonHandle;
  iEvent.getByToken(l1MuonSrc_, l1MuonHandle);

  edm::Handle<vector<pat::Muon>> muonHandle;
  iEvent.getByToken(muonSrc_, muonHandle);

  edm::Handle<vector<reco::Vertex>> vtxHandle;
  iEvent.getByToken(vtxSrc_, vtxHandle);
  auto primaryVtx = (*vtxHandle)[0];
  hAllNvts->Fill((int)vtxHandle->size());
  for(auto vtx : (*vtxHandle)) hAllVtxZ->Fill(vtx.position().z());

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamSpotSrc_, beamSpotHandle);

  // Output collection
  unique_ptr<vector<pat::Muon>> trgMuonsMatched( new vector<pat::Muon> );
  unique_ptr<map<string, float>> outputNtuplizer(new map<string, float>);
  unique_ptr<map<string, vector<float>>> outputVecNtuplizer(new map<string, vector<float>>);

  // Save MC pileup information
  if (isMC_) {
    edm::Handle<vector<PileupSummaryInfo>> pileupMCHandle;
    iEvent.getByToken(pileupMCSrc_, pileupMCHandle);
    uint idxBX0 = -1;
    for (uint i=0; i < pileupMCHandle->size(); i++) {
      if (pileupMCHandle->at(i).getBunchCrossing() == 0) {
        idxBX0 = i;
        break;
      }
    }
    float puMC = pileupMCHandle->at(idxBX0).getTrueNumInteractions();
    (*outputNtuplizer)["nTrueIntMC"] = puMC;
    hAllNTrueIntMC->Fill(puMC);
  }

  vector<string> triggerTag = {"Mu12_IP6", "Mu9_IP6", "Mu7_IP4"};

  if (verbose) {cout << "\n MUONS LIST" << endl;}
  for (auto muon : (*muonHandle)) {
    if (verbose) {
      cout << "\t" << Form("id:%d  pt:%.1f  eta:%.1f  phi:%.1f softID:%d", muon.pdgId(), muon.pt(), muon.eta(), muon.phi(), muon.isSoftMuon(primaryVtx)) << endl;
    }
    if(muonCharge_ != 0 && muon.charge() != muonCharge_) continue;
    if(!muon.triggered("HLT_Mu*_IP*_part*_v*")) continue;
    trgMuonsMatched->push_back(muon);

    auto out = matchL1Muon(muon, *l1MuonHandle);
    (*outputVecNtuplizer)["trgMu_L1_pt"].push_back(get<3>(out));
    if (get<0>(out) == 9999) (*outputVecNtuplizer)["trgMu_L1_eta"].push_back(-999);
    else (*outputVecNtuplizer)["trgMu_L1_eta"].push_back( l1MuonHandle->at(0,get<0>(out)).eta() );
    (*outputVecNtuplizer)["trgMu_L1_dR"].push_back( get<1>(out) );

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
      (*outputVecNtuplizer)["trgMu_dz"].push_back(tk->dz(primaryVtx.position()));
      auto dxyUnc = tk->dxyError();
      auto dxy_BS = fabs(tk->dxy((*beamSpotHandle)));
      (*outputVecNtuplizer)["trgMu_dxy_BS"].push_back(dxy_BS);
      (*outputVecNtuplizer)["trgMu_sigdxy_BS"].push_back( dxy_BS/dxyUnc );
    }
    else {
      (*outputVecNtuplizer)["trgMu_dz"].push_back(-999);
      (*outputVecNtuplizer)["trgMu_dxy_BS"].push_back(-999);
      (*outputVecNtuplizer)["trgMu_sigdxy_BS"].push_back(-999);
    }
  }
  bool acceptEvent = trgMuonsMatched->size() > 0;
  (*outputNtuplizer)["N_trgMu"] = trgMuonsMatched->size();
  if(verbose) {cout << "\nTriggered muons: " << trgMuonsMatched->size() << endl;}

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

tuple<uint, float, float, float> TriggerMuonsFilterLight::matchL1Muon(pat::Muon muReco, BXVector<l1t::Muon> muonsL1, uint skipIdx) {
  uint idxMatch = 9999;
  float best_dR = 1e6;
  float best_dpt = 1e6;
  float best_pt = -1;

  // Approximately radius 2.5 m in a 3.8T magnetic field
  float dphi = 2 * TMath::ASin(2.5 * 0.3 * 3.8 / (2 * muReco.pt()) );
  float phiProp = muReco.phi() + TMath::Sign(1, muReco.pdgId()) * dphi;

  for (uint i=0; i < muonsL1.size(0); i++) {
    if (i == skipIdx) continue;
    auto m = muonsL1.at(0,i);
    if (m.hwQual() < 12) continue;
    float dR = vtxu::dR(m.phi(), phiProp, m.eta(), muReco.eta());
    float dpt = fabs(muReco.pt() - m.pt())/muReco.pt();
    if ((dR < best_dR && dpt < best_dpt) || (dpt + dR < best_dpt + best_dR)) {
      best_dR = dR;
      best_dpt = dpt;
      idxMatch = i;
      best_pt = m.pt();
    }
  }

  tuple<uint, float, float, float> out(idxMatch, best_dR, best_dpt, best_pt);
  return out;
}

DEFINE_FWK_MODULE(TriggerMuonsFilterLight);
