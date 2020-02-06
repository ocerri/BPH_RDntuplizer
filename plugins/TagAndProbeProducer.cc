// system include files
#include <memory>
#include <iostream>
#include <string>
#include <regex>

#include "TLorentzVector.h"
#include "TTree.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "VtxUtils.hh"

using namespace std;

class TagAndProbeProducer : public edm::stream::EDFilter<> {
   public:
      explicit TagAndProbeProducer(const edm::ParameterSet&);
      ~TagAndProbeProducer() {};

   private:
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      void addToTree();

      // ----------member data ---------------------------
      edm::EDGetTokenT<edm::TriggerResults> triggerBitsSrc_;
      edm::EDGetTokenT <pat::PackedTriggerPrescales> triggerPrescales_;
      edm::EDGetTokenT<vector<pat::Muon>> muonSrc_;
      edm::EDGetTokenT<vector<reco::Vertex>> vtxSrc_;

      TTree* tree;
      map<string, float> outMap;
      bool treeDeclared = false;

      int verbose = 0;

      double massMu  = 0.10565;
      double massJpsi  = 3.09691;
};


TagAndProbeProducer::TagAndProbeProducer(const edm::ParameterSet& iConfig):
  triggerBitsSrc_(consumes<edm::TriggerResults>( edm::InputTag("TriggerResults","","HLT") )),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>( edm::InputTag("patTrigger") )),
  muonSrc_( consumes<vector<pat::Muon>> ( edm::InputTag("slimmedMuons") ) ),
  vtxSrc_( consumes<vector<reco::Vertex>> ( edm::InputTag("offlineSlimmedPrimaryVertices") ) ),
  verbose( iConfig.getParameter<int>( "verbose" ) )
{
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>( "T", "Events Tree from TAG AND PROBE");
}

bool TagAndProbeProducer::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBitsSrc_, triggerBits);

  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByToken(triggerPrescales_, triggerPrescales);

  edm::Handle<vector<pat::Muon>> muonHandle;
  iEvent.getByToken(muonSrc_, muonHandle);
  uint nMuons = muonHandle->size();
  if(nMuons < 2) return false;

  edm::Handle<vector<reco::Vertex>> vtxHandle;
  iEvent.getByToken(vtxSrc_, vtxHandle);
  auto primaryVtx = (*vtxHandle)[0];

  if (verbose) {
    cout << "\n MUONS LIST" << endl;
    for (auto muon : (*muonHandle)) {
        cout << Form("id:%d  pt:%.1f  eta:%.1f  phi:%.1f", muon.pdgId(), muon.pt(), muon.eta(), muon.phi());
        cout << Form(" tightID:%d softID:%d trgBPH:%d", muon.isTightMuon(primaryVtx), muon.isSoftMuon(primaryVtx), muon.triggered("HLT_Mu*_IP*_part*_v*"));
        cout << endl;
      }
  }

  vector<string> triggerTag = {"Mu12_IP6", "Mu9_IP5", "Mu7_IP4", "Mu9_IP4", "Mu8_IP5", "Mu8_IP6", "Mu9_IP6", "Mu8_IP3"};
  for(auto tag : triggerTag) outMap["prescale" + tag] = -1;

  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  regex txt_regex_path("HLT_Mu[0-9]+_IP[0-9]_part[0-9]_v[0-9]");
  if (verbose) { cout << "\n == TRIGGER PATHS == \n";}
  for (unsigned int i = 0; i < triggerBits->size(); ++i) {
      auto n =  names.triggerName(i);
      if(!regex_match(n, txt_regex_path)) continue;
      for(auto tag : triggerTag) {
        if(n.substr(4, tag.size()) == tag) {
          outMap["prescale" + tag] = triggerPrescales->getPrescaleForIndex(i);
          if(verbose) {cout << tag << "\t" << n << "\t" << triggerBits->wasrun(i) << "\t" << triggerPrescales->getPrescaleForIndex(i) << endl;}
        }
      }
  }

  for(uint i=0; i < nMuons-1; i++) {
    auto mTag = (*muonHandle)[i];
    if(!mTag.isTightMuon(primaryVtx)) continue;
    for(uint j=i+1; j < nMuons; j++) {
        auto mProbe = (*muonHandle)[j];
        if(mTag.charge()*mProbe.charge() != -1) continue;

        TLorentzVector pTag, pProbe;
        pTag.SetPtEtaPhiM(mTag.pt(), mTag.eta(), mTag.phi(), massMu);
        pProbe.SetPtEtaPhiM(mProbe.pt(), mProbe.eta(), mProbe.phi(), massMu);

        if( fabs((pTag + pProbe).M() - massJpsi) > 0.5 ) continue;

        outMap["mTag_pt"] = mTag.pt();
        outMap["mTag_eta"] = mTag.eta();
        outMap["mTag_phi"] = mTag.phi();
        outMap["mTag_hasInnerTrk"] = !mProbe.innerTrack().isNull();
        outMap["mProbe_pt"] = mProbe.pt();
        outMap["mProbe_eta"] = mProbe.eta();
        outMap["mProbe_phi"] = mProbe.phi();
        outMap["massInv"] = (pTag + pProbe).M();
        for(auto tag : triggerTag) {
          string trgPath = "HLT_" + tag + "_part*_v*";
          outMap["mProbe_HLT_" + tag] = mProbe.triggered(trgPath.c_str());
        }
        outMap["mProbe_hasInnerTrk"] = !mProbe.innerTrack().isNull();
        outMap["mProbe_tightID"] = mProbe.isTightMuon(primaryVtx);
        outMap["mProbe_softID"] = mProbe.isSoftMuon(primaryVtx);
        if(!mProbe.innerTrack().isNull()) {
          auto tk = mProbe.innerTrack();
          auto dxy = tk->dxy(primaryVtx.position());
          outMap["mProbe_dz"] = tk->dz(primaryVtx.position());
          outMap["mProbe_dxy"] = dxy;
          outMap["mProbe_sigdxy"] = fabs(dxy)/tk->dxyError();

          if(!mProbe.innerTrack().isNull()) {
            auto kinTree = vtxu::FitJpsi_mumu(iSetup, mTag, mProbe, false);
            auto res = vtxu::fitQuality(kinTree, 0.1);
            outMap["vtx_isValid"] = res.isValid;
            outMap["vtx_chi2"] = res.chi2;
            outMap["vtx_dof"] = res.dof;
            outMap["vtx_pval"] = res.pval;
            outMap["vtx_isGood"] = res.isGood;
          }
          else {
            outMap["vtx_isValid"] = 0;
            outMap["vtx_chi2"] = -1;
            outMap["vtx_dof"] = -1;
            outMap["vtx_pval"] = -1;
            outMap["vtx_isGood"] = 0;
          }
        }
        else {
          outMap["mProbe_dz"] = 0;
          outMap["mProbe_dxy"] = 0;
          outMap["mProbe_sigdxy"] = 0;
          outMap["vtx_isValid"] = 0;
          outMap["vtx_chi2"] = -1;
          outMap["vtx_dof"] = -1;
          outMap["vtx_pval"] = -1;
          outMap["vtx_isGood"] = 0;
        }
        addToTree();
    }
  }
  if (verbose) {cout << "======================== " << endl;}
  return true;
}

void TagAndProbeProducer::addToTree() {
  if (!treeDeclared) {
    if(verbose) {cout << "\nCreating the branches in the output tree:\n";}
    for(auto& kv : outMap) {
      auto k = kv.first;
      if(verbose) {cout << "\t" << k;}
      tree->Branch(k.c_str(), &(outMap[k]));
    }
    treeDeclared = true;
    if(verbose) {cout << "\n\n";}
  }

  tree->Fill();
}

DEFINE_FWK_MODULE(TagAndProbeProducer);
