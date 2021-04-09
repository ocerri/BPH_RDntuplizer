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
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

// L1 trigger
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1Trigger/interface/Muon.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "VtxUtils.hh"

using namespace std;

class TagAndProbeProducer_MC : public edm::stream::EDFilter<> {
   public:
      explicit TagAndProbeProducer_MC(const edm::ParameterSet&);
      ~TagAndProbeProducer_MC() {
        cout << "Total events in output tree: " << tree->GetEntries() << endl;
      };

   private:
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      tuple<uint, float, float, float> matchL1Muon(pat::Muon, BXVector<l1t::Muon>);
      void addToTree();

      // ----------member data ---------------------------
      edm::EDGetTokenT<BXVector<GlobalAlgBlk>> L1triggerBitsSrc_;
      edm::EDGetTokenT<edm::TriggerResults> L1triggerResultsSrc_;
      edm::EDGetTokenT<edm::TriggerResults> triggerBitsSrc_;
      edm::EDGetTokenT <pat::PackedTriggerPrescales> triggerPrescales_;
      edm::EDGetTokenT<BXVector<l1t::Muon>> l1MuonSrc_;
      edm::EDGetTokenT<vector<pat::Muon>> muonSrc_;
      edm::EDGetTokenT<vector<reco::Vertex>> vtxSrc_;
      edm::EDGetTokenT<reco::BeamSpot> beamSpotSrc_;

      edm::EDGetTokenT<vector<PileupSummaryInfo>> pileupMCSrc_;

      edm::Service<TFileService> fs;

      TH1I* hAllNvts;
      TH1I* hAllNTrueIntMC;
      TH1I* hAllVtxZ;
      TTree* tree;
      map<string, float> outMap;
      bool treeDeclared = false;

      bool isRealData;
      unsigned int runNum;
      unsigned int lumiNum;
      unsigned long long eventNum;

      int muonIDScaleFactors = false;
      int requireTag = 1;
      int verbose = 0;

      double massMu  = 0.10565;
      double massJpsi  = 3.09691;

      TH2D* hMuonIdSF;
};


TagAndProbeProducer_MC::TagAndProbeProducer_MC(const edm::ParameterSet& iConfig):
  L1triggerBitsSrc_(consumes<BXVector<GlobalAlgBlk>>( edm::InputTag("gtStage2Digis") )),
  L1triggerResultsSrc_(consumes<edm::TriggerResults>( edm::InputTag("l1bits") )),
  triggerBitsSrc_(consumes<edm::TriggerResults>( edm::InputTag("TriggerResults","","HLT") )),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>( edm::InputTag("patTrigger") )),
  l1MuonSrc_( consumes<BXVector<l1t::Muon>> ( edm::InputTag("gmtStage2Digis", "Muon", "RECO") ) ),
  muonSrc_( consumes<vector<pat::Muon>> ( edm::InputTag("slimmedMuons") ) ),
  vtxSrc_( consumes<vector<reco::Vertex>> ( edm::InputTag("offlineSlimmedPrimaryVertices") ) ),
  beamSpotSrc_( consumes<reco::BeamSpot> ( edm::InputTag("offlineBeamSpot") ) ),
  pileupMCSrc_( consumes<vector<PileupSummaryInfo>> ( edm::InputTag("slimmedAddPileupInfo") ) ),
  muonIDScaleFactors( iConfig.getParameter<int>( "muonIDScaleFactors" ) ),
  requireTag( iConfig.getParameter<int>( "requireTag" ) ),
  verbose( iConfig.getParameter<int>( "verbose" ) )
{
  hAllNvts = fs->make<TH1I>("hAllNvts", "Number of vertexes from all the MINIAOD events", 101, -0.5, 100.5);
  hAllNTrueIntMC = fs->make<TH1I>("hAllNTrueIntMC", "Number of true interactions generated in MC", 101, -0.5, 100.5);
  hAllVtxZ = fs->make<TH1I>("hAllVtxZ", "Z coordinate of vertexes from all the MINIAOD events", 100, -25, 25);
  tree = fs->make<TTree>( "T", "Events Tree from TAG AND PROBE");

  if (muonIDScaleFactors) {
    TFile fAux = TFile("/storage/user/ocerri/BPhysics/data/calibration/muonIDscaleFactors/Run2018ABCD_SF_MuonID_Jpsi.root", "READ");
    hMuonIdSF = (TH2D*) fAux.Get("NUM_SoftID_DEN_genTracks_pt_abseta");
  }

}

bool TagAndProbeProducer_MC::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  isRealData = iEvent.isRealData() ? 1 : 0 ;
  runNum     = iEvent.id().run();
  lumiNum    = iEvent.luminosityBlock();
  eventNum   = iEvent.id().event();

  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBitsSrc_, triggerBits);

  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByToken(triggerPrescales_, triggerPrescales);


  edm::Handle<vector<reco::Vertex>> vtxHandle;
  iEvent.getByToken(vtxSrc_, vtxHandle);
  auto primaryVtx = (*vtxHandle)[0];
  auto nRecoVtx = vtxHandle->size();
  outMap["nVtx"] = nRecoVtx;
  hAllNvts->Fill((int)nRecoVtx);
  for(auto vtx : (*vtxHandle)) hAllVtxZ->Fill(vtx.position().z());

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
  outMap["nTrueIntMC"] = puMC;
  hAllNTrueIntMC->Fill((int)puMC);

  edm::Handle<vector<pat::Muon>> muonHandle;
  iEvent.getByToken(muonSrc_, muonHandle);
  uint nMuons = muonHandle->size();
  if(nMuons < 2) return false;
  if (verbose) {cout << "\n\n =================  Event " << eventNum << " =================== " << endl;}


  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamSpotSrc_, beamSpotHandle);

  // L1 triggers NanoAOD
  // edm::Handle<edm::TriggerResults> L1triggerResults;
  // iEvent.getByToken(L1triggerResultsSrc_, L1triggerResults);

  // if (verbose) {cout << "======== L1 trigger NanoAOD ======== : " << L1triggerResults->size() << endl;}
  //
  // // regex regex_MuNano(".*SingleMu[0-9]+.*");
  //
  // for (unsigned int i = 0; i < L1triggerResults->size(); ++i) {
  //     auto n =  L1triggerResults->name(i);
  //     auto result = L1triggerResults->accept(i);
  //     if (verbose) {
  //       cout << n << ": " << result << endl;
  //     }
  // }
  //
  // if (verbose) {cout << "====================================" << endl;}
  //
  // // L1 trigger
  //
  // edm::ESHandle<L1GtTriggerMenu> L1TriggerMenuHandle;
  // iSetup.get<L1GtTriggerMenuRcd>().get(L1TriggerMenuHandle);
  //
  // edm::Handle<BXVector<GlobalAlgBlk>> L1TriggerBitsHandle;
  // iEvent.getByToken(L1triggerBitsSrc_, L1TriggerBitsHandle);
  // const std::vector<bool>* l1FinalDecision = &L1TriggerBitsHandle->at(0, 0).getAlgoDecisionFinal();
  //
  // auto prescaleColumn = L1TriggerBitsHandle->at(0, 0).getPreScColumn();
  // outMap["l1PrescaleColumn"] = prescaleColumn;
  // if (verbose) {cout << "--> L1 prescale column: " << prescaleColumn << endl;}
  //
  // vector<string> consideredL1Seeds = {"22", "25", "18er1p5", "14er1p5", "12er1p5", "10er1p5", "9er1p5", "8er1p5", "7er1p5", "6er1p5"};
  // for (auto ptThr : consideredL1Seeds) {
  //   string name = "L1_SingleMu"+ptThr;
  //   outMap[name] = -1;
  // }
  //
  // regex regex_Mu(".*Mu.*");
  // if (verbose) {cout << "======== L1 trigger results ========" << endl;}
  //
  // for (auto const& keyval : L1TriggerMenuHandle->gtAlgorithmAliasMap()) {
  //   auto name = keyval.first;
  //   auto idx  = keyval.second.algoBitNumber();
  //   bool result = (*l1FinalDecision)[idx];
  //   if (verbose && regex_match(name, regex_Mu)) {
  //     cout << name << ": " << (int)result << endl;
  //   }
  //   for (auto ptThr : consideredL1Seeds) {
  //     regex auxRegex("L1_SingleMu"+ptThr);
  //     if (regex_match(name, auxRegex)) {
  //       outMap[name] = result;
  //       break;
  //     }
  //   }
  //   // if (regex_match(name, regex_Mu)) outMap[name] = result;
  // }
  // if (verbose) {cout << "====================================" << endl;}

  edm::Handle<BXVector<l1t::Muon>> l1MuonHandle;
  iEvent.getByToken(l1MuonSrc_, l1MuonHandle);
  uint nL1Muons = (uint)l1MuonHandle->size(0);

  if (verbose) {
    cout << "============= L1 muons =============" << endl;
    for (uint i=0; i < nL1Muons; i++) {
      auto m = l1MuonHandle->at(0,i);
      cout << i << ": " << Form("pt=%.2f, eta=%.2f, phi=%.2f, quality=%d", m.pt(), m.eta(), m.phi(), m.hwQual()) << endl;
    }
    cout << "====================================" << endl;
  }

  // HLT trigger

  vector<string> triggerTag = {"Mu12_IP6", "Mu9_IP5", "Mu7_IP4", "Mu9_IP4", "Mu9_IP6"};
  for(auto tag : triggerTag) outMap["prescale" + tag] = -1;

  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  regex txt_regex_path("HLT_Mu[0-9]+_IP[0-9]_part[0-9]_v[0-9]");
  if (verbose) { cout << "\n == TRIGGER PATHS == \n";}
  for (unsigned int i = 0; i < triggerBits->size(); ++i) {
      auto n =  names.triggerName(i);
      if (verbose && triggerBits->accept(i) && false) {
        cout << n << ": PASS" << endl;
      }
      if(!regex_match(n, txt_regex_path)) continue;
      for(auto tag : triggerTag) {
        if(n.substr(4, tag.size()) == tag) {
          outMap["prescale" + tag] = triggerPrescales->getPrescaleForIndex(i);
          // if(verbose) {cout << tag << "\t" << n << "\t" << triggerBits->wasrun(i) << "\t" << triggerPrescales->getPrescaleForIndex(i) << endl;}
        }
      }
  }

  if (verbose) {
    cout << "\n MUONS LIST" << endl;
    for (auto muon : (*muonHandle)) {
        cout << Form("id:%d  pt:%.1f  eta:%.1f  phi:%.1f", muon.pdgId(), muon.pt(), muon.eta(), muon.phi());
        cout << Form(" tightID:%d softID:%d trgBPH:%d", muon.isTightMuon(primaryVtx), muon.isSoftMuon(primaryVtx), muon.triggered("HLT_Mu*_IP*_part*_v*"));
        if(muon.pt() > 5.5) {
          auto out = matchL1Muon(muon, *l1MuonHandle);
          cout << ", L1 match " << get<0>(out) << Form(": dR=%.4f, dpt=%.3f", get<1>(out), get<2>(out));
        }
        cout << endl;
      }
  }

  vector<uint> idxTriggeringMuons;
  for(uint i=0; i < nMuons; i++) {
    auto m = (*muonHandle)[i];
    if(m.triggered("HLT_Mu*_IP*")) idxTriggeringMuons.push_back(i);
  }
  if(idxTriggeringMuons.size() == 0 && requireTag) return false;


  for(uint j=0; j < nMuons; j++) {
    if(idxTriggeringMuons.size() == 1 && idxTriggeringMuons[0] == j && requireTag) continue;

    auto mProbe = (*muonHandle)[j];
    if ( fabs(mProbe.eta()) > 1.6 ) continue;
    if ( mProbe.pt() < 5.5 ) continue;
    if (mProbe.innerTrack().isNull()) continue;
    auto tkProbe = mProbe.innerTrack();
    auto dxyProbe_PV = fabs(tkProbe->dxy(primaryVtx.position()));
    auto dxyProbe_BS = fabs(tkProbe->dxy((*beamSpotHandle)));

    auto dxyProbeUnc = tkProbe->dxyError();

    if (dxyProbe_BS/dxyProbeUnc < 1) continue;
    TLorentzVector pProbe;
    pProbe.SetPtEtaPhiM(mProbe.pt(), mProbe.eta(), mProbe.phi(), massMu);

    pat::Muon mTag;
    TLorentzVector pTag;
    double ptTagMu = -1;
    for(auto i : idxTriggeringMuons) {
      if(i==j) continue;
      auto m = (*muonHandle)[i];
      if(m.charge()*mProbe.charge() != -1) continue;
      if(m.innerTrack().isNull()) continue;
      TLorentzVector pAux;
      pAux.SetPtEtaPhiM(m.pt(), m.eta(), m.phi(), massMu);
      if( fabs((pAux + pProbe).M() - massJpsi) > 0.2 ) continue;
      if(m.pt() > ptTagMu){
        ptTagMu = m.pt();
        mTag = m;
        pTag = pAux;
      }
    }
    if (requireTag && ptTagMu == -1) continue;

    if(ptTagMu != -1) {
      outMap["mTag_pt"] = mTag.pt();
      outMap["mTag_eta"] = mTag.eta();
      outMap["mTag_phi"] = mTag.phi();

      auto tkTag = mTag.innerTrack();
      auto dxyTag_PV = fabs(tkTag->dxy(primaryVtx.position()));
      auto dxyTag_BS = fabs(tkTag->dxy((*beamSpotHandle)));
      outMap["mTag_dxy_PV"] = dxyTag_PV;
      outMap["mTag_sigdxy_PV"] = dxyTag_PV/tkTag->dxyError();
      outMap["mTag_dxy_BS"] = dxyTag_BS;
      outMap["mTag_sigdxy_BS"] = dxyTag_BS/tkTag->dxyError();

      for(auto tag : triggerTag) {
        string trgPath = "HLT_" + tag + "_part*_v*";
        outMap["mTag_HLT_" + tag] = mTag.triggered(trgPath.c_str());
      }
      outMap["mTag_tightID"] = mTag.isTightMuon(primaryVtx);
      outMap["mTag_softID"] = mTag.isSoftMuon(primaryVtx);

      auto out = matchL1Muon(mTag, *l1MuonHandle);
      outMap["mTag_L1_pt"] = get<3>(out);
      if (get<0>(out) == 999) outMap["mTag_L1_eta"] = -999;
      else outMap["mTag_L1_eta"] = l1MuonHandle->at(0,get<0>(out)).eta();
      outMap["mTag_L1_dR"] = get<1>(out);
    }
    else {
      outMap["mTag_pt"] = -1;
      outMap["mTag_eta"] = -9999;
      outMap["mTag_phi"] = -9999;
      outMap["mTag_dxy_PV"] = -1;
      outMap["mTag_sigdxy_PV"] = -1;
      outMap["mTag_dxy_BS"] = -1;
      outMap["mTag_sigdxy_BS"] = -1;
      for(auto tag : triggerTag) outMap["mTag_HLT_" + tag] = -1;
      outMap["mTag_tightID"] = -1;
      outMap["mTag_softID"] = -1;
      outMap["mTag_L1_pt"] = -1;
      outMap["mTag_L1_eta"] = -999;
      outMap["mTag_L1_dR"] = -1;
    }


    outMap["mProbe_pt"] = mProbe.pt();
    outMap["mProbe_eta"] = mProbe.eta();
    outMap["mProbe_phi"] = mProbe.phi();
    outMap["mProbe_dz"] = tkProbe->dz(primaryVtx.position());
    outMap["mProbe_dxy_PV"] = dxyProbe_PV;
    outMap["mProbe_sigdxy_PV"] = dxyProbe_PV/dxyProbeUnc;
    outMap["mProbe_dxy_BS"] = dxyProbe_BS;
    outMap["mProbe_sigdxy_BS"] = dxyProbe_BS/dxyProbeUnc;

    auto out = matchL1Muon(mProbe, *l1MuonHandle);
    outMap["mProbe_L1_pt"] = get<3>(out);
    if (get<0>(out) == 999) outMap["mProbe_L1_eta"] = -999;
    else outMap["mProbe_L1_eta"] = l1MuonHandle->at(0,get<0>(out)).eta();
    outMap["mProbe_L1_dR"] = get<1>(out);

    for(auto tag : triggerTag) {
      string trgPath = "HLT_" + tag + "_part*_v*";
      outMap["mProbe_HLT_" + tag] = mProbe.triggered(trgPath.c_str());
    }
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

    if(ptTagMu != -1) {
      outMap["massMuMu"] = (pTag + pProbe).M();
      outMap["deltaR_tagProbe"] = vtxu::dR(pTag.Phi(), mProbe.phi(), pTag.Eta(), mProbe.eta());
      auto kinTree = vtxu::FitJpsi_mumu(iSetup, mTag, mProbe, false);
      auto res = vtxu::fitQuality(kinTree, 0.05);
      outMap["vtx_isValid"] = res.isValid;
      outMap["vtx_chi2"] = res.chi2;
      outMap["vtx_dof"] = res.dof;
      outMap["vtx_pval"] = res.pval;
      outMap["vtx_isGood"] = res.isGood;
      if (res.isValid) {
        kinTree->movePointerToTheTop();
        auto mass = kinTree->currentParticle()->currentState().mass();
        outMap["massMuMu_refit"] = mass;
      }
      else outMap["massMuMu_refit"] = -1;
    }
    else {
      outMap["massMuMu"] = -1;
      outMap["deltaR_tagProbe"] = -1;
      outMap["vtx_isValid"] = 0;
      outMap["vtx_chi2"] = -1;
      outMap["vtx_dof"] = -1;
      outMap["vtx_pval"] = -1;
      outMap["vtx_isGood"] = 0;
      outMap["massMuMu_refit"] = -1;
    }

    addToTree();
  }
  if (verbose) {cout << "======================== " << endl;}
  return true;
}

tuple<uint, float, float, float> TagAndProbeProducer_MC::matchL1Muon(pat::Muon muReco, BXVector<l1t::Muon> muonsL1) {
  uint idxMatch = 999;
  float best_dR = 1e6;
  float best_dpt = 1e6;
  float best_pt = -1;

  // Approximately radius 2.5 m in a 3.8T magnetic field
  float dphi = 2 * TMath::ASin(2.5 * 0.3 * 3.8 / (2 * muReco.pt()) );
  float phiProp = muReco.phi() + TMath::Sign(1, muReco.pdgId()) * dphi;

  for (uint i=0; i < muonsL1.size(0); i++) {
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

void TagAndProbeProducer_MC::addToTree() {
  if (!treeDeclared) {
    if(verbose) {cout << "\nCreating the branches in the output tree:\n";}
    tree->Branch("isRealData", &isRealData);
    tree->Branch("runNum", &runNum);
    tree->Branch("lumiNum", &lumiNum);
    tree->Branch("eventNum", &eventNum);

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

DEFINE_FWK_MODULE(TagAndProbeProducer_MC);
