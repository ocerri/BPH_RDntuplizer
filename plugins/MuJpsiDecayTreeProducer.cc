#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include <iostream>
#include <string>
#include <map>
#include <TMath.h>

#include "VtxUtils.hh"

#define __PvalChi2Vtx_max__ 0.95 // Very loose cut
#define __dzMax__ 3.0
#define __dmJpsi_max__ 1.0 // loose cut
#define __sigIPpfCand_min__ 2. // loose cut
#define __pThad_min__ 0.5 // loose cut
#define __dmKst_max__ 0.5 // loose cut
#define __dmB0_max__ 0.5 // loose cut

using namespace std;

class MuJpsiDecayTreeProducer : public edm::EDProducer {

public:

    explicit MuJpsiDecayTreeProducer(const edm::ParameterSet &iConfig);
    void AddTLVToOut(TLorentzVector, string, map<string, vector<float>>*);
    bool isMuonFromJpsiID(pat::Muon, reco::Vertex, pat::Muon);

    ~MuJpsiDecayTreeProducer() override {};

private:
    virtual void produce(edm::Event&, const edm::EventSetup&);

    // ----------member data ---------------------------
    edm::EDGetTokenT<vector<pat::PackedCandidate>> PFCandSrc_;
    edm::EDGetTokenT<vector<reco::Vertex>> vtxSrc_;
    edm::EDGetTokenT<vector<pat::Muon>> TrgMuonSrc_;
    edm::EDGetTokenT<vector<pat::Muon>> muonSrc_;

    double mass_Mu   = 0.10565;
    double mass_Pi   = 0.13957;
    double mass_K    = 0.49367;
    double mass_Jpsi = 3.09691;
    double mass_Kst  = 0.89166;
    double mass_B0   = 5.27963;

    int verbose = 0;
};



MuJpsiDecayTreeProducer::MuJpsiDecayTreeProducer(const edm::ParameterSet &iConfig)
{
    verbose = iConfig.getParameter<int>( "verbose" );
    PFCandSrc_ = consumes<vector<pat::PackedCandidate>>(edm::InputTag("packedPFCandidates"));
    muonSrc_ = consumes<vector<pat::Muon>>(edm::InputTag("slimmedMuons"));
    vtxSrc_ = consumes<vector<reco::Vertex>>(edm::InputTag("offlineSlimmedPrimaryVertices"));
    TrgMuonSrc_ = consumes<vector<pat::Muon>> ( iConfig.getParameter<edm::InputTag>( "trgMuons" ) );

    produces<map<string, float>>("outputNtuplizer");
    produces<map<string, vector<float>>>("outputVecNtuplizer");
}


void MuJpsiDecayTreeProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    if (verbose) {cout <<"-------------------- Evt -----------------------\n";}

    // Output collection
    unique_ptr<map<string, float>> outputNtuplizer(new map<string, float>);
    (*outputNtuplizer)["n_B"] = 0;
    unique_ptr<map<string, vector<float>>> outputVecNtuplizer(new map<string, vector<float>>);

    edm::Handle<vector<pat::Muon>> muonHandle;
    iEvent.getByToken(muonSrc_, muonHandle);
    unsigned int nMu = muonHandle->size();
    if(nMu < 3) {
      iEvent.put(move(outputNtuplizer), "outputNtuplizer");
      iEvent.put(move(outputVecNtuplizer), "outputVecNtuplizer");
      return;
    }

    edm::Handle<vector<reco::Vertex>> vtxHandle;
    iEvent.getByToken(vtxSrc_, vtxHandle);
    auto primaryVtx = (*vtxHandle)[0];
    if (verbose) {
      cout << "Event with " << nMu << " muons" << endl;
      int idxMu = 0;
      for(auto m : *muonHandle ) {
        cout << Form("i: %d, pt: %.1f, eta: %.1f, phi: %.1f, SoftId: ", idxMu, m.pt(), m.eta(), m.phi()) << m.isSoftMuon(primaryVtx);
        cout << Form(", dz: %.2f, hasInnerTrack: ", m.muonBestTrack()->dz(primaryVtx.position())) << !m.innerTrack().isNull() << endl;
        idxMu++;
      }
    }

    edm::Handle<vector<pat::Muon>> trgMuonsHandle;
    iEvent.getByToken(TrgMuonSrc_, trgMuonsHandle);
    auto trgMu = (*trgMuonsHandle)[0];

    int n_Jpsi = 0;
    for(uint i1 = 0; i1 < nMu -1; i1++){
      auto m1 = (*muonHandle)[i1];
      if (!isMuonFromJpsiID(m1, primaryVtx, trgMu)) continue;

      for(uint i2 = i1+1; i2 < nMu; i2++){
        auto m2 = (*muonHandle)[i2];
        if (!isMuonFromJpsiID(m2, primaryVtx, trgMu)) continue;
        if(m1.charge() * m2.charge() != -1) continue;

        auto mup = m1;
        auto mum = m2;
        if(m1.charge() < 0) {
          mum = m1;
          mup = m2;
        }
        if(verbose){cout << "Fitting the J/Psi from: " << Form("%d %d", i1, i2) << endl;}
        auto kinTree = vtxu::FitJpsi_mumu(iSetup, mup, mum, false);
        bool accept = false;
        double chi2;
        if(kinTree->isValid()) {
          kinTree->movePointerToTheTop();
          auto vtx = kinTree->currentDecayVertex();
          chi2 = vtx->chiSquared();
          if(verbose){cout << Form("chi2 geom: %.2f (max %.2f)", chi2, TMath::ChisquareQuantile(__PvalChi2Vtx_max__, vtx->degreesOfFreedom())) << endl;}
          accept = chi2 > 0 && chi2 < TMath::ChisquareQuantile(__PvalChi2Vtx_max__, vtx->degreesOfFreedom());
        }
        if(!accept) continue;
        double massMuPair = kinTree->currentParticle()->currentState().mass();
        accept &= fabs(massMuPair - mass_Jpsi) < __dmJpsi_max__;
        if(!accept) continue;

        kinTree = vtxu::FitJpsi_mumu(iSetup, mup, mum, true);
        accept = false;
        double chi2_mass;
        if(kinTree->isValid()) {
          kinTree->movePointerToTheTop();
          auto vtx = kinTree->currentDecayVertex();
          chi2_mass = vtx->chiSquared();
          if(verbose){cout << Form("chi2 mass: %.2f (max %.2f)", chi2_mass, TMath::ChisquareQuantile(__PvalChi2Vtx_max__, vtx->degreesOfFreedom())) << endl;}
          // accept = chi2_mass > 0 && chi2_mass < TMath::ChisquareQuantile(__PvalChi2Vtx_max__, vtx->degreesOfFreedom());
          accept = true;
        }
        if(!accept) continue;

        kinTree->movePointerToTheTop();
        auto Jpsi = kinTree->currentParticle();
        auto vtxJpsi = kinTree->currentDecayVertex();

        auto cos_Jpsi_vtxMu = vtxu::computePointingCos(primaryVtx, vtxJpsi, Jpsi);
        auto d_vtxJpsi_PV = vtxu::vtxsDistance(primaryVtx, vtxJpsi);
        auto sigd_vtxJpsi_PV = fabs(d_vtxJpsi_PV.first/d_vtxJpsi_PV.second);

        // Save variables in the output
        n_Jpsi++;
        kinTree->movePointerToTheFirstChild();
        auto refit_MuPlus = kinTree->currentParticle();
        kinTree->movePointerToTheNextChild();
        auto refit_MuMinus = kinTree->currentParticle();

        AddTLVToOut(vtxu::getTLVfromKinPart(refit_MuPlus), string("mup"), &(*outputVecNtuplizer));
        auto dxy_mup = mup.innerTrack()->dxy(primaryVtx.position());
        (*outputVecNtuplizer)["mup_dxy"].push_back(dxy_mup);
        (*outputVecNtuplizer)["mup_sigdxy"].push_back(fabs(dxy_mup)/mup.innerTrack()->dxyError());
        float mup_isTrg = trgMu.pt()==mup.pt() && trgMu.phi()==mup.phi() && trgMu.eta()==mup.eta();
        (*outputVecNtuplizer)["mup_isTrg"].push_back(mup_isTrg);

        AddTLVToOut(vtxu::getTLVfromKinPart(refit_MuMinus), string("mum"), &(*outputVecNtuplizer));
        auto dxy_mum = mum.innerTrack()->dxy(primaryVtx.position());
        (*outputVecNtuplizer)["mum_dxy"].push_back(dxy_mum);
        (*outputVecNtuplizer)["mum_sigdxy"].push_back(fabs(dxy_mum)/mum.innerTrack()->dxyError());
        float mum_isTrg = trgMu.pt()==mum.pt() && trgMu.phi()==mum.phi() && trgMu.eta()==mum.eta();
        (*outputVecNtuplizer)["mum_isTrg"].push_back(mum_isTrg);

        (*outputVecNtuplizer)["chi2_mumu"].push_back(chi2);
        (*outputVecNtuplizer)["mass_mumu"].push_back(massMuPair);
        (*outputVecNtuplizer)["chi2_Jpsi"].push_back(chi2_mass);
        (*outputVecNtuplizer)["cos_Jpsi_vtxMu"].push_back(cos_Jpsi_vtxMu);
        (*outputVecNtuplizer)["d_vtxJpsi_PV"].push_back(d_vtxJpsi_PV.first);
        (*outputVecNtuplizer)["sigd_vtxJpsi_PV"].push_back(sigd_vtxJpsi_PV);
        AddTLVToOut(vtxu::getTLVfromKinPart(Jpsi), string("Jpsi"), &(*outputVecNtuplizer));

      }
    }

    if(verbose) {cout << n_Jpsi << " Jpsi candidate found" << endl;}

    (*outputNtuplizer)["Run"] = iEvent.run();
    (*outputNtuplizer)["LumiBlock"] = iEvent.luminosityBlock();
    (*outputNtuplizer)["eventNumber"] = iEvent.id().event();
    (*outputNtuplizer)["n_B"] = n_Jpsi;


    iEvent.put(move(outputNtuplizer), "outputNtuplizer");
    iEvent.put(move(outputVecNtuplizer), "outputVecNtuplizer");
    return;
}


void MuJpsiDecayTreeProducer::AddTLVToOut(TLorentzVector v, string n, map<string, vector<float>>* outv) {
  (*outv)[n+"_pt"].push_back(v.Pt());
  (*outv)[n+"_pz"].push_back(v.Pz());
  (*outv)[n+"_eta"].push_back(v.Eta());
  (*outv)[n+"_phi"].push_back(v.Phi());
  return;
}


bool MuJpsiDecayTreeProducer::isMuonFromJpsiID(pat::Muon m, reco::Vertex pVtx, pat::Muon trgMu) {
  if(m.innerTrack().isNull()) return false;
  if (fabs(m.innerTrack()->dz(pVtx.position()) - trgMu.innerTrack()->dz(pVtx.position())) > __dzMax__) return false;
  // if(trgMu.pt()==m.pt() && trgMu.phi()==m.phi() && trgMu.eta()==m.eta()) return false;
  // if(!m.isSoftMuon(pVtx)) return false;
  if(m.innerTrack()->hitPattern().pixelLayersWithMeasurement() < 2) return false;
  if(!m.innerTrack()->quality(reco::TrackBase::highPurity)) return false;
  if(!m.isGood("TMOneStationTight")) return false;
  if(m.innerTrack()->normalizedChi2() > 1.8) return false;
  return true;
}

DEFINE_FWK_MODULE(MuJpsiDecayTreeProducer);
