// system include files
#include <memory>
#include <iostream>
#include <string>
#include <regex>

#include "TLorentzVector.h"
#include "TTree.h"
#include "TH2D.h"

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

      edm::Service<TFileService> fs;

      TH1I* hAllNvts;
      TH1I* hAllVtxZ;
      TTree* tree;
      map<string, float> outMap;
      bool treeDeclared = false;

      int muonIDScaleFactors = false;
      int verbose = 0;

      double massMu  = 0.10565;
      double massJpsi  = 3.09691;

      TH2D* hMuonIdSF;
};


TagAndProbeProducer::TagAndProbeProducer(const edm::ParameterSet& iConfig):
  triggerBitsSrc_(consumes<edm::TriggerResults>( edm::InputTag("TriggerResults","","HLT") )),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>( edm::InputTag("patTrigger") )),
  muonSrc_( consumes<vector<pat::Muon>> ( edm::InputTag("slimmedMuons") ) ),
  vtxSrc_( consumes<vector<reco::Vertex>> ( edm::InputTag("offlineSlimmedPrimaryVertices") ) ),
  muonIDScaleFactors( iConfig.getParameter<int>( "muonIDScaleFactors" ) ),
  verbose( iConfig.getParameter<int>( "verbose" ) )
{
  hAllNvts = fs->make<TH1I>("hAllNvts", "Number of vertexes from all the MINIAOD events", 101, -0.5, 100.5);
  hAllVtxZ = fs->make<TH1I>("hAllVtxZ", "Z coordinate of vertexes from all the MINIAOD events", 100, -25, 25);
  tree = fs->make<TTree>( "T", "Events Tree from TAG AND PROBE");

  if (muonIDScaleFactors) {
    TFile fAux = TFile("/storage/user/ocerri/BPhysics/data/calibration/muonIDscaleFactors/Run2018ABCD_SF_MuonID_Jpsi.root", "READ");
    hMuonIdSF = (TH2D*) fAux.Get("NUM_SoftID_DEN_genTracks_pt_abseta");
  }

}

bool TagAndProbeProducer::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBitsSrc_, triggerBits);

  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByToken(triggerPrescales_, triggerPrescales);


  edm::Handle<vector<reco::Vertex>> vtxHandle;
  iEvent.getByToken(vtxSrc_, vtxHandle);
  auto primaryVtx = (*vtxHandle)[0];
  auto nRecoVtx = vtxHandle->size();
  hAllNvts->Fill((int)nRecoVtx);
  for(auto vtx : (*vtxHandle)) hAllVtxZ->Fill(vtx.position().z());

  edm::Handle<vector<pat::Muon>> muonHandle;
  iEvent.getByToken(muonSrc_, muonHandle);
  uint nMuons = muonHandle->size();
  if(nMuons < 2) return false;

  if (verbose) {
    cout << "\n MUONS LIST" << endl;
    for (auto muon : (*muonHandle)) {
        cout << Form("id:%d  pt:%.1f  eta:%.1f  phi:%.1f", muon.pdgId(), muon.pt(), muon.eta(), muon.phi());
        cout << Form(" tightID:%d softID:%d trgBPH:%d", muon.isTightMuon(primaryVtx), muon.isSoftMuon(primaryVtx), muon.triggered("HLT_Mu*_IP*_part*_v*"));
        cout << endl;
      }
  }

  vector<string> triggerTag = {"Mu12_IP6", "Mu9_IP5", "Mu7_IP4", "Mu9_IP4", "Mu9_IP6"};
  for(auto tag : triggerTag) outMap["prescale_" + tag] = -1;

  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  regex txt_regex_path("HLT_Mu[0-9]+_IP[0-9]_part[0-9]_v[0-9]");
  if (verbose) { cout << "\n == TRIGGER PATHS == \n";}
  for (unsigned int i = 0; i < triggerBits->size(); ++i) {
      auto n =  names.triggerName(i);
      if (verbose && triggerBits->accept(i)) {cout << n << ": PASS" << endl;}
      if(!regex_match(n, txt_regex_path)) continue;
      for(auto tag : triggerTag) {
        if(n.substr(4, tag.size()) == tag) {
          outMap["prescale" + tag] = triggerPrescales->getPrescaleForIndex(i);
          if(verbose) {cout << tag << "\t" << n << "\t" << triggerBits->wasrun(i) << "\t" << triggerPrescales->getPrescaleForIndex(i) << endl;}
        }
      }
  }

  vector<uint> idxTriggeringMuons;
  for(uint i=0; i < nMuons; i++) {
    auto m = (*muonHandle)[i];
    if(m.triggered("HLT_Mu*_IP*")) idxTriggeringMuons.push_back(i);
  }
  if(idxTriggeringMuons.size() == 0) return false;


  for(uint j=0; j < nMuons; j++) {
    if(idxTriggeringMuons.size() == 1 && idxTriggeringMuons[0] == j) continue;
    auto mProbe = (*muonHandle)[j];
    if ( fabs(mProbe.eta()) > 1.6 ) continue;
    if ( mProbe.pt() < 5. ) continue;
    TLorentzVector pTag, pProbe;
    pProbe.SetPtEtaPhiM(mProbe.pt(), mProbe.eta(), mProbe.phi(), massMu);

    pat::Muon mTag;
    double ptTagMu = -1;
    for(uint i=0; i < nMuons; i++) {
      if(i==j) continue;
      auto m = (*muonHandle)[i];
      if(!m.isSoftMuon(primaryVtx)) continue;
      if(m.charge()*mProbe.charge() != -1) continue;
      pTag.SetPtEtaPhiM(m.pt(), m.eta(), m.phi(), massMu);
      if( fabs((pTag + pProbe).M() - massJpsi) > 0.5 ) continue;
      if(m.pt() > ptTagMu){
        ptTagMu = m.pt();
        mTag = m;
      }
    }
    if(ptTagMu == -1) pTag.SetPtEtaPhiM(0, -99999999, -9999999999, massMu);
    else pTag.SetPtEtaPhiM(mTag.pt(), mTag.eta(), mTag.phi(), massMu);

    outMap["nVtx"] = nRecoVtx;
    outMap["mTag_pt"] = mTag.pt();
    outMap["mTag_eta"] = mTag.eta();
    outMap["mTag_phi"] = mTag.phi();
    if(ptTagMu == -1) outMap["mTag_hasInnerTrk"] = 0;
    else outMap["mTag_hasInnerTrk"] = !mTag.innerTrack().isNull();
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
    if (muonIDScaleFactors) {
      if (outMap["mProbe_softID"]) {
        auto ix = hMuonIdSF->GetXaxis()->FindBin(min(39.9,mProbe.pt()));
        auto iy = hMuonIdSF->GetYaxis()->FindBin(min(2.4, fabs(mProbe.eta())));
        outMap["sfMuonID"] = hMuonIdSF->GetBinContent(ix, iy);
      }
      else outMap["sfMuonID"] = 1.;
    }
    if (!outMap["mProbe_softID"]) continue;
    if(!mProbe.innerTrack().isNull()) {
      auto tk = mProbe.innerTrack();
      auto dxy = tk->dxy(primaryVtx.position());
      outMap["mProbe_dz"] = tk->dz(primaryVtx.position());
      outMap["mProbe_dxy"] = dxy;
      outMap["mProbe_sigdxy"] = fabs(dxy)/tk->dxyError();

      if(outMap["mProbe_hasInnerTrk"] && outMap["mTag_hasInnerTrk"]) {
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
