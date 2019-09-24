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
#define __dzMax__ 2.0
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

    vector<RefCountedKinematicTree> vecJpsiKinTree;
    vector<double> vecMassMuMu;
    vector<double> vecChi2MuMu;
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

        if(accept) {
          vecJpsiKinTree.push_back(kinTree);
          vecMassMuMu.push_back(massMuPair);
          vecChi2MuMu.push_back(chi2);
        }
      }
    }

    if(vecJpsiKinTree.size() < 1) {
      iEvent.put(move(outputNtuplizer), "outputNtuplizer");
      iEvent.put(move(outputVecNtuplizer), "outputVecNtuplizer");
      return;
    }
    else if(verbose) {cout << vecJpsiKinTree.size() << " Jpsi candidate found" << endl;}

    edm::Handle<vector<pat::PackedCandidate>> pfCandHandle;
    iEvent.getByToken(PFCandSrc_, pfCandHandle);
    unsigned int N_pfCand = pfCandHandle->size();

    int n_K = 0, n_pi = 0, n_Kst = 0, n_B = 0;

    /*
    ############################################################################
                              Look for the K+ (321)
    ############################################################################
    */
    for(uint i_k = 0; i_k < N_pfCand; ++i_k) {
      const pat::PackedCandidate & K = (*pfCandHandle)[i_k];
      if (!K.hasTrackDetails()) continue;
      //Require a positive charged hadron
      if (K.isTrackerMuon() || K.isStandAloneMuon()) continue;
      if (K.charge() < 0) continue;
      //Require a minimum pt
      if(K.pt() < __pThad_min__) continue;
      // Require to be close to the trigger muon;
      auto K_tk = K.bestTrack();
      if (fabs(K_tk->dz(primaryVtx.position()) - trgMu.muonBestTrack()->dz(primaryVtx.position())) > __dzMax__) continue;
      // Require significant impact parameter
      auto dxy = K_tk->dxy(primaryVtx.position());
      auto sigd_K_PV = fabs(dxy)/K_tk->dxyError();
      auto K_norm_chi2 = K_tk->normalizedChi2();
      auto K_N_valid_hits = K_tk->numberOfValidHits();
      if (sigd_K_PV < __sigIPpfCand_min__) continue;

      n_K++;
      if(verbose){cout << Form("K cand found i:%d, pt:%.2f, eta:%.2f, phi:%.2f ", i_k, K.pt(), K.eta(), K.phi()) << endl;}
      /*
      ############################################################################
                               Look for the pi- (-211)
      ############################################################################
      */
      for(uint i_pi = 0; i_pi < N_pfCand; ++i_pi) {
        const pat::PackedCandidate & pi = (*pfCandHandle)[i_pi];
        if (!pi.hasTrackDetails()) continue;
        //Require a negative charged hadron
        if (pi.isTrackerMuon() || pi.isStandAloneMuon()) continue;
        if (pi.charge() > 0) continue;
        //Require a minimum pt
        if(pi.pt() < __pThad_min__) continue;
        // Require to be close to the trigger muon;
        auto pi_tk = pi.bestTrack();
        if (fabs(pi_tk->dz(primaryVtx.position()) - trgMu.muonBestTrack()->dz(primaryVtx.position())) > __dzMax__) continue;
        // Require significant impact parameter
        auto dxy = pi_tk->dxy(primaryVtx.position());
        auto sigd_pi_PV = fabs(dxy)/pi.dxyError();
        auto pi_norm_chi2 = pi_tk->normalizedChi2();
        auto pi_N_valid_hits = pi_tk->numberOfValidHits();
        if (sigd_pi_PV < __sigIPpfCand_min__) continue;

        n_pi++;
        if(verbose){cout << Form("pi cand found i:%d, pt:%.2f, eta:%.2f, phi:%.2f ", i_pi, pi.pt(), pi.eta(), pi.phi()) << endl;}

        //Fit the vertex w/o mass constraint
        auto KstKinTree = vtxu::FitKst_piK(iSetup, pi, K, false);
        bool accept = false;
        double chi2_Kpi;
        if(KstKinTree->isValid()) {
          KstKinTree->movePointerToTheTop();
          auto vtxKst = KstKinTree->currentDecayVertex();
          chi2_Kpi = vtxKst->chiSquared();
          if(verbose){cout << "chi2 geom: " << chi2_Kpi << endl;}
          accept = chi2_Kpi > 0 && chi2_Kpi < TMath::ChisquareQuantile(__PvalChi2Vtx_max__, vtxKst->degreesOfFreedom());
        }
        if(!accept) continue;
        auto Kst = KstKinTree->currentParticle();
        auto mass_Kpi = Kst->currentState().mass();
        bool accept_Kst = fabs(mass_Kpi - mass_Kst) < __dmKst_max__;
        if(!accept_Kst) continue;

        //Fit the vertex WITH mass constraint
        KstKinTree = vtxu::FitKst_piK(iSetup, pi, K, true);
        accept = false;
        double chi2_KpiMass;
        if(KstKinTree->isValid()) {
          KstKinTree->movePointerToTheTop();
          auto vtxKst = KstKinTree->currentDecayVertex();
          chi2_KpiMass = vtxKst->chiSquared();
          if(verbose){cout << "chi2 mass: " << chi2_KpiMass << endl;}
          // accept = chi2_KpiMass > 0 && chi2_KpiMass < TMath::ChisquareQuantile(__PvalChi2Vtx_max__, vtxKst->degreesOfFreedom());
          accept = true;
        }
        if(!accept) continue;
        Kst = KstKinTree->currentParticle();
        auto vtxKst = KstKinTree->currentDecayVertex();

        auto cos_Kst_PV = vtxu::computePointingCos(primaryVtx, vtxKst, Kst);
        auto d_vtxKst_PV = vtxu::vtxsDistance(primaryVtx, vtxKst);
        auto sigd_vtxKst_PV = d_vtxKst_PV.first/d_vtxKst_PV.second;

        if (verbose) {cout << "K*_0 -> K+ pi- found\n";}
        n_Kst++;

        /*
        ############################################################################
                  Make a B0 -> J/psi K* using the candidates found above
        ############################################################################
        */
        for(uint i_J = 0; i_J < vecJpsiKinTree.size(); ++i_J) {
          auto JpsiKinTree = vecJpsiKinTree[i_J];
          JpsiKinTree->movePointerToTheTop();
          auto Jpsi = JpsiKinTree->currentParticle();
          auto vtxJpsi = JpsiKinTree->currentDecayVertex();

          auto cos_Jpsi_vtxMu = vtxu::computePointingCos(primaryVtx, vtxJpsi, Jpsi);
          auto d_vtxJpsi_PV = vtxu::vtxsDistance(primaryVtx, vtxJpsi);
          auto sigd_vtxJpsi_PV = fabs(d_vtxJpsi_PV.first/d_vtxJpsi_PV.second);

          // Fit the B vertex
          RefCountedKinematicTree BKinTree;
          try{BKinTree = vtxu::FitVtxJpsiKst(iSetup, Jpsi, Kst, false);}
          catch(...) {cout << "Fitting B -> J/psi K* failed" << endl; continue;}
          bool accept = false;
          double chi2_JpsiKst;
          if(BKinTree->isValid()) {
            BKinTree->movePointerToTheTop();
            auto vtxB = BKinTree->currentDecayVertex();
            chi2_JpsiKst = vtxB->chiSquared();
            accept = chi2_JpsiKst > 0 && chi2_JpsiKst < TMath::ChisquareQuantile(__PvalChi2Vtx_max__, vtxB->degreesOfFreedom());
          }
          if(!accept) continue;
          auto pairJpsiKst = BKinTree->currentParticle();
          auto mass_JpsiKst = pairJpsiKst->currentState().mass();
          accept = fabs(mass_JpsiKst - mass_B0) < __dmB0_max__;
          if(!accept) continue;
          else if(verbose) {
            cout << Form("Candidate B0 mass (%.2f GeV) check passed", mass_JpsiKst) << endl;
          }

          BKinTree = vtxu::FitVtxJpsiKst(iSetup, Jpsi, Kst, true);
          accept = false;
          double chi2_B;
          if(BKinTree->isValid()) {
            BKinTree->movePointerToTheTop();
            auto vtxB = BKinTree->currentDecayVertex();
            chi2_B = vtxB->chiSquared();
            // accept = chi2_B > 0 && chi2_B < TMath::ChisquareQuantile(__PvalChi2Vtx_max__, vtxB->degreesOfFreedom());
            accept = true;
          }
          if(!accept) continue;
          auto B = BKinTree->currentParticle();
          auto vtxB = BKinTree->currentDecayVertex();

          double cos_B_PV = vtxu::computePointingCos(primaryVtx, vtxB, B);
          auto d_vtxB_PV = vtxu::vtxsDistance(primaryVtx, vtxB);
          double sigd_vtxB_PV = d_vtxB_PV.first/d_vtxB_PV.second;
          if (verbose) {
            cout << Form("PV: sigd = %.2f  cos = %.3e", sigd_vtxB_PV, 1-cos_B_PV) << endl;
            cout << "B0 -> J/psi K* found\n";
          }
          n_B++;

          /*
          ############################################################################
                                Compute analysis variables
          ############################################################################
          */
          KstKinTree->movePointerToTheFirstChild();
          auto refit_pi = KstKinTree->currentParticle();
          AddTLVToOut(vtxu::getTLVfromKinPart(refit_pi), string("pi"), &(*outputVecNtuplizer));
          (*outputVecNtuplizer)["sigd_pi_PV"].push_back(sigd_pi_PV);
          (*outputVecNtuplizer)["pi_norm_chi2"].push_back(pi_norm_chi2);
          (*outputVecNtuplizer)["pi_N_valid_hits"].push_back(pi_N_valid_hits);
          KstKinTree->movePointerToTheNextChild();
          auto refit_K = KstKinTree->currentParticle();
          AddTLVToOut(vtxu::getTLVfromKinPart(refit_K), string("K"), &(*outputVecNtuplizer));
          (*outputVecNtuplizer)["sigd_K_PV"].push_back(sigd_K_PV);
          (*outputVecNtuplizer)["K_norm_chi2"].push_back(K_norm_chi2);
          (*outputVecNtuplizer)["K_N_valid_hits"].push_back(K_N_valid_hits);
          (*outputVecNtuplizer)["chi2_Kpi"].push_back(chi2_Kpi);
          (*outputVecNtuplizer)["mass_Kpi"].push_back(mass_Kpi);
          (*outputVecNtuplizer)["chi2_Kst"].push_back(chi2_KpiMass);
          (*outputVecNtuplizer)["cos_Kst_PV"].push_back(cos_Kst_PV);
          (*outputVecNtuplizer)["d_vtxKst_PV"].push_back(d_vtxKst_PV.first);
          (*outputVecNtuplizer)["sigd_vtxKst_PV"].push_back(sigd_vtxKst_PV);


          JpsiKinTree->movePointerToTheFirstChild();
          auto refit_Mu1 = JpsiKinTree->currentParticle();
          JpsiKinTree->movePointerToTheNextChild();
          auto refit_Mu2 = JpsiKinTree->currentParticle();
          AddTLVToOut(vtxu::getTLVfromKinPart(refit_Mu1), string("mup"), &(*outputVecNtuplizer));
          AddTLVToOut(vtxu::getTLVfromKinPart(refit_Mu2), string("mum"), &(*outputVecNtuplizer));
          (*outputVecNtuplizer)["chi2_mumu"].push_back(vecChi2MuMu[i_J]);
          (*outputVecNtuplizer)["mass_mumu"].push_back(vecMassMuMu[i_J]);
          (*outputVecNtuplizer)["chi2_Jpsi"].push_back(vtxJpsi->chiSquared());
          (*outputVecNtuplizer)["cos_Jpsi_vtxMu"].push_back(cos_Jpsi_vtxMu);
          (*outputVecNtuplizer)["d_vtxJpsi_PV"].push_back(d_vtxJpsi_PV.first);
          (*outputVecNtuplizer)["sigd_vtxJpsi_PV"].push_back(sigd_vtxJpsi_PV);

          (*outputVecNtuplizer)["chi2_JpsiKst"].push_back(chi2_JpsiKst);
          (*outputVecNtuplizer)["mass_JpsiKst"].push_back(mass_JpsiKst);
          (*outputVecNtuplizer)["chi2_B"].push_back(chi2_B);
          (*outputVecNtuplizer)["cos_B_PV"].push_back(cos_B_PV);
          (*outputVecNtuplizer)["d_vtxB_PV"].push_back(d_vtxB_PV.first);
          (*outputVecNtuplizer)["sigd_vtxB_PV"].push_back(sigd_vtxB_PV);

          AddTLVToOut(vtxu::getTLVfromKinPart(B), string("B"), &(*outputVecNtuplizer));
          BKinTree->movePointerToTheFirstChild();
          auto refit_Jpsi = BKinTree->currentParticle();
          AddTLVToOut(vtxu::getTLVfromKinPart(refit_Jpsi), string("Jpsi"), &(*outputVecNtuplizer));
          BKinTree->movePointerToTheNextChild();
          auto refit_Kst = BKinTree->currentParticle();
          AddTLVToOut(vtxu::getTLVfromKinPart(refit_Kst), string("Kst"), &(*outputVecNtuplizer));
        }
      }
    }

    (*outputNtuplizer)["Run"] = iEvent.run();
    (*outputNtuplizer)["LumiBlock"] = iEvent.luminosityBlock();
    (*outputNtuplizer)["eventNumber"] = iEvent.id().event();

    (*outputNtuplizer)["n_K"] = n_K;
    (*outputNtuplizer)["n_pi"] = n_pi;
    (*outputNtuplizer)["n_Kst"] = n_Kst;
    (*outputNtuplizer)["n_Jpsi"] = vecJpsiKinTree.size();
    (*outputNtuplizer)["n_B"] = n_B;


    iEvent.put(move(outputNtuplizer), "outputNtuplizer");
    iEvent.put(move(outputVecNtuplizer), "outputVecNtuplizer");
    return;
}


void MuJpsiDecayTreeProducer::AddTLVToOut(TLorentzVector v, string n, map<string, vector<float>>* outv) {
  (*outv)[n+"_pt"].push_back(v.Pt());
  (*outv)[n+"_pz"].push_back(v.Pz());
  (*outv)[n+"_eta"].push_back(v.Eta());
  (*outv)[n+"_phi"].push_back(v.Phi());
  (*outv)[n+"_P"].push_back(v.P());
  (*outv)[n+"_E"].push_back(v.E());
  return;
}


bool MuJpsiDecayTreeProducer::isMuonFromJpsiID(pat::Muon m, reco::Vertex pVtx, pat::Muon trgMu) {
  if(m.innerTrack().isNull()) return false;
  if (fabs(m.innerTrack()->dz(pVtx.position()) - trgMu.innerTrack()->dz(pVtx.position())) > __dzMax__) return false;
  if(trgMu.pt()==m.pt() && trgMu.phi()==m.phi() && trgMu.eta()==m.eta()) return false;
  if (!m.isSoftMuon(pVtx)) return false;
  if(m.innerTrack()->hitPattern().pixelLayersWithMeasurement() < 2) return false;
  if(m.innerTrack()->normalizedChi2() > 1.8) return false;
  return true;
}

DEFINE_FWK_MODULE(MuJpsiDecayTreeProducer);
