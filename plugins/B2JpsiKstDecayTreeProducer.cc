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

#define __PvalChi2Vtx_max__ 0.99 // Very loose cut
#define __dzMax__ 1.0
#define __dmJpsi_max__ 1.0 // loose cut
#define __sigIPpfCand_min__ 1.5 // loose cut
#define __dRMax__ 3.0
#define __cos_Kpi_vtxMu_min__ 0.97
#define __d_vtxKpi_vtxMu_min__ 0.02
#define __dmKst_max__ 0.5 // loose cut
#define __dmB0_max__ 0.5 // loose cut
#define __cos_JpsiKst_vtxMu_min__ 0.95 // loose cut

using namespace std;

class B2JpsiKstDecayTreeProducer : public edm::EDProducer {

public:

    explicit B2JpsiKstDecayTreeProducer(const edm::ParameterSet &iConfig);
    void DeclareTLVToOut(string, map<string, vector<float>>*);
    void AddTLVToOut(TLorentzVector, string, map<string, vector<float>>*);

    ~B2JpsiKstDecayTreeProducer() override {};

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



B2JpsiKstDecayTreeProducer::B2JpsiKstDecayTreeProducer(const edm::ParameterSet &iConfig)
{
    verbose = iConfig.getParameter<int>( "verbose" );
    PFCandSrc_ = consumes<vector<pat::PackedCandidate>>(edm::InputTag("packedPFCandidates"));
    muonSrc_ = consumes<vector<pat::Muon>>(edm::InputTag("slimmedMuons"));
    vtxSrc_ = consumes<vector<reco::Vertex>>(edm::InputTag("offlineSlimmedPrimaryVertices"));
    TrgMuonSrc_ = consumes<vector<pat::Muon>> ( iConfig.getParameter<edm::InputTag>( "trgMuons" ) );

    produces<map<string, float>>("outputNtuplizer");
    produces<map<string, vector<float>>>("outputVecNtuplizer");
}


void B2JpsiKstDecayTreeProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
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
    // else if (verbose) {
    //   cout << "Event with " << nMu << " muons" << endl;
    //   int idxMu = 0;
    //   for(auto m : *muonHandle ) {
    //     cout << Form("i: %d, pt: %.1f, eta: %.1f, phi: %.1f, LooseId: ", idxMu, m.pt(), m.eta(), m.phi()) << m.isLooseMuon();
    //     cout << Form(", dz: %.2f", m.muonBestTrack()->dz()) << endl;
    //     idxMu++;
    //   }
    // }

    edm::Handle<vector<pat::Muon>> trgMuonsHandle;
    iEvent.getByToken(TrgMuonSrc_, trgMuonsHandle);
    auto trgMu = (*trgMuonsHandle)[0];

    vector<RefCountedKinematicTree> vecJpsiKinTree;
    vector<double> vecMassMuMu;
    vector<double> vecChi2MuMu;
    for(uint i1 = 0; i1 < nMu -1; i1++){
      auto m1 = (*muonHandle)[i1];
      if (!m1.isLooseMuon()) continue;
      if (fabs(m1.muonBestTrack()->dz() - trgMu.muonBestTrack()->dz()) > __dzMax__) continue;
      if(trgMu.pt()==m1.pt() && trgMu.phi()==m1.phi() && trgMu.eta()==m1.eta()) continue;

      for(uint i2 = i1+1; i2 < nMu; i2++){
        auto m2 = (*muonHandle)[i2];
        if (!m2.isLooseMuon()) continue;
        if (fabs(m2.muonBestTrack()->dz() - trgMu.muonBestTrack()->dz()) > __dzMax__) continue;
        if(trgMu.pt()==m2.pt() && trgMu.phi()==m2.phi() && trgMu.eta()==m2.eta()) continue;
        if(m1.charge() * m2.charge() != -1) continue;

        auto mup = m1;
        auto mum = m2;
        if(m1.charge() < 0) {
          mum = m1;
          mup = m2;
        }
        auto kinTree = vtxu::FitJpsi_mumu(iSetup, mup, mum, false);
        bool accept = false;
        double chi2;
        if(kinTree->isValid()) {
          kinTree->movePointerToTheTop();
          auto vtx = kinTree->currentDecayVertex();
          chi2 = vtx->chiSquared();
          accept = chi2 > 0 && chi2 < TMath::ChisquareQuantile(__PvalChi2Vtx_max__, vtx->degreesOfFreedom());
        }
        if(!accept) continue;
        auto massMuPair = kinTree->currentParticle()->currentState().mass();
        accept &= fabs(massMuPair - mass_Jpsi) < __dmJpsi_max__;
        if(!accept) continue;

        kinTree = vtxu::FitJpsi_mumu(iSetup, mup, mum, true);
        accept = false;
        double chi2_mass;
        if(kinTree->isValid()) {
          kinTree->movePointerToTheTop();
          auto vtx = kinTree->currentDecayVertex();
          chi2_mass = vtx->chiSquared();
          accept = chi2_mass > 0 && chi2_mass < TMath::ChisquareQuantile(__PvalChi2Vtx_max__, vtx->degreesOfFreedom());
        }
        if(!accept) continue;
        auto muPair = kinTree->currentParticle();

        // Add J/psi vertex displacement
        // TODO: add here

        if(accept) {
          vecJpsiKinTree.push_back(kinTree);
          vecMassMuMu.push_back(massMuPair);
          vecChi2MuMu.push_back(chi2);
          if(verbose) {
            cout << "Mass Jpsi: " << massMuPair << endl;
          }
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

    edm::Handle<vector<reco::Vertex>> vtxHandle;
    iEvent.getByToken(vtxSrc_, vtxHandle);


    int n_K = 0, n_pi = 0, n_Kst = 0, n_B = 0;

    pat::PackedCandidate trgMu_pf;
    uint i_trgMu = 0;
    for (i_trgMu = 0; i_trgMu < N_pfCand; ++i_trgMu) {
      const pat::PackedCandidate & p = (*pfCandHandle)[i_trgMu];
      if (fabs(trgMu.pt() - p.pt())/trgMu.pt() > 5e-3) continue;
      if (fabs(trgMu.eta() - p.eta()) > 5e-3) continue;
      if (fabs(vtxu::dPhi(trgMu.phi(), p.phi())) > 5e-3) continue;
      trgMu_pf = p;
      break;
    }

    auto vtxMu = trgMu_pf.vertexRef();
    GlobalPoint p_vtxMu(vtxMu->x(), vtxMu->y(), vtxMu->z());

    // (*outputVecNtuplizer)["chi2_Kpi"] = {};
    // (*outputVecNtuplizer)["mass_Kpi"] = {};
    // (*outputVecNtuplizer)["cos_Kpi_vtxMu"] = {};
    // (*outputVecNtuplizer)["d_vtxKpi_vtxMu"] = {};
    // (*outputVecNtuplizer)["sigd_vtxKpi_vtxMu"] = {};
    //
    // (*outputVecNtuplizer)["chi2_D0pis"] = {};
    // (*outputVecNtuplizer)["mass_D0pis"] = {};
    // (*outputVecNtuplizer)["cos_D0pis_vtxMu"] = {};
    // (*outputVecNtuplizer)["d_vtxD0pis_vtxMu"] = {};
    // (*outputVecNtuplizer)["sigd_vtxD0pis_vtxMu"] = {};
    //
    // (*outputVecNtuplizer)["d_Ks_PV"] = {};
    // (*outputVecNtuplizer)["sigd_Ks_PV"] = {};
    // (*outputVecNtuplizer)["chi2_JpsiKst"] = {};
    // (*outputVecNtuplizer)["mass_JpsiKsts"] = {};
    // (*outputVecNtuplizer)["cos_JpsiKst_vtxBest"] = {};
    // (*outputVecNtuplizer)["d_vtxJpsiKst_vtxBest"] = {};
    // (*outputVecNtuplizer)["sigd_vtxJpsiKst_vtxBest"] = {};
    //
    // DeclareTLVToOut("B", &(*outputVecNtuplizer));
    // DeclareTLVToOut("Ks", &(*outputVecNtuplizer));
    // DeclareTLVToOut("Dst", &(*outputVecNtuplizer));
    // DeclareTLVToOut("D0", &(*outputVecNtuplizer));
    // DeclareTLVToOut("pis", &(*outputVecNtuplizer));
    // DeclareTLVToOut("K", &(*outputVecNtuplizer));
    // DeclareTLVToOut("pi", &(*outputVecNtuplizer));

    /*
    ############################################################################
                              Look for the K+ (321)
    ############################################################################
    */
    for(uint i_k = 0; i_k < N_pfCand; ++i_k) {
      const pat::PackedCandidate & K = (*pfCandHandle)[i_k];
      if (!K.hasTrackDetails()) continue;
      //Require a positive charged hadron
      if (K.pdgId() != 211 ) continue;
      // Require to be close to the trigger muon;
      if (fabs(K.bestTrack()->dz() - trgMu.muonBestTrack()->dz()) > __dzMax__) continue;
      // Require significant impact parameter
      auto sigd_K_PV = fabs(K.dxy()/K.dxyError());
      auto K_tk = K.bestTrack();
      auto K_norm_chi2 = K_tk->normalizedChi2();
      auto K_N_valid_hits = K_tk->numberOfValidHits();
      if (sigd_K_PV < __sigIPpfCand_min__) continue;

      n_K++;

      /*
      ############################################################################
                               Look for the pi- (-211)
      ############################################################################
      */
      for(uint i_pi = 0; i_pi < N_pfCand; ++i_pi) {
        const pat::PackedCandidate & pi = (*pfCandHandle)[i_pi];
        if (!pi.hasTrackDetails()) continue;
        //Require a negative charged hadron
        if (pi.pdgId() != -211 ) continue;
        // Require to be close to the trigger muon;
        if (fabs(pi.bestTrack()->dz() - trgMu.muonBestTrack()->dz()) > __dzMax__) continue;
        // Require to be close to the K;
        if (fabs(pi.bestTrack()->dz() - K.bestTrack()->dz()) > __dzMax__) continue;
        if (vtxu::dR(pi.phi(), K.phi(), pi.eta(), K.eta()) > __dRMax__) continue;
        // Require significant impact parameter
        auto sigd_pi_PV = fabs(pi.dxy()/pi.dxyError());
        auto pi_tk = pi.bestTrack();
        auto pi_norm_chi2 = pi_tk->normalizedChi2();
        auto pi_N_valid_hits = pi_tk->numberOfValidHits();
        if (sigd_pi_PV < __sigIPpfCand_min__) continue;

        n_pi++;

        //Fit the vertex w/o mass constraint
        auto KstKinTree = vtxu::FitKst_piK(iSetup, pi, K, false);
        bool accept = false;
        double chi2_Kpi;
        if(KstKinTree->isValid()) {
          KstKinTree->movePointerToTheTop();
          auto Kstvtx = KstKinTree->currentDecayVertex();
          chi2_Kpi = Kstvtx->chiSquared();
          accept = chi2_Kpi > 0 && chi2_Kpi < TMath::ChisquareQuantile(__PvalChi2Vtx_max__, Kstvtx->degreesOfFreedom());
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
          auto Kstvtx = KstKinTree->currentDecayVertex();
          chi2_KpiMass = Kstvtx->chiSquared();
          accept = chi2_KpiMass > 0 && chi2_KpiMass < TMath::ChisquareQuantile(__PvalChi2Vtx_max__, Kstvtx->degreesOfFreedom());
        }
        if(!accept) continue;
        Kst = KstKinTree->currentParticle();
        auto Kstvtx = KstKinTree->currentDecayVertex();

        TVector3 dvtx(Kstvtx->position().x() - vtxMu->x(),
                      Kstvtx->position().y() - vtxMu->y(),
                      Kstvtx->position().z() - vtxMu->z()
                     );
        TVector3 pKst(Kst->currentState().globalMomentum().x(),
                     Kst->currentState().globalMomentum().y(),
                     Kst->currentState().globalMomentum().z()
                     );
        double dalpha = dvtx.Angle(pKst);
        auto cos_Kpi_vtxMu = cos(dalpha);

        auto dvtx_KpiPV = vtxu::vtxsDistance(vtxMu, Kstvtx);
        auto d_vtxKpi_vtxMu = dvtx_KpiPV.first;
        auto sigd_vtxKpi_vtxMu = fabs(d_vtxKpi_vtxMu/dvtx_KpiPV.second);

        // Make selection cut to get good Kst
        accept_Kst &= cos_Kpi_vtxMu > __cos_Kpi_vtxMu_min__;
        accept_Kst &= d_vtxKpi_vtxMu > __d_vtxKpi_vtxMu_min__;

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
          auto Jpsivtx = JpsiKinTree->currentDecayVertex();
          TVector3 dvtx(Jpsivtx->position().x() - vtxMu->x(),
                        Jpsivtx->position().y() - vtxMu->y(),
                        Jpsivtx->position().z() - vtxMu->z()
                       );
          TVector3 pJpsi(Jpsi->currentState().globalMomentum().x(),
                       Jpsi->currentState().globalMomentum().y(),
                       Jpsi->currentState().globalMomentum().z()
                       );
          double dalpha = dvtx.Angle(pJpsi);
          auto cos_mumu_vtxMu = cos(dalpha);

          auto dvtx_mumuPV = vtxu::vtxsDistance(vtxMu, Jpsivtx);
          auto d_vtxmumu_vtxMu = dvtx_mumuPV.first;
          auto sigd_vtxmumu_vtxMu = fabs(d_vtxmumu_vtxMu/dvtx_mumuPV.second);

          // Fit the B vertex
          RefCountedKinematicTree BKinTree;
          try{BKinTree = vtxu::FitVtxJpsiKst(iSetup, Jpsi, Kst, false);}
          catch(...) {cout << "Fitting B -> J/psi K* failed" << endl; continue;}
          bool accept = false;
          double chi2_JpsiKst;
          if(BKinTree->isValid()) {
            BKinTree->movePointerToTheTop();
            auto Bvtx = BKinTree->currentDecayVertex();
            chi2_JpsiKst = Bvtx->chiSquared();
            accept = chi2_JpsiKst > 0 && chi2_JpsiKst < TMath::ChisquareQuantile(__PvalChi2Vtx_max__, Bvtx->degreesOfFreedom());
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
            auto Bvtx = BKinTree->currentDecayVertex();
            chi2_B = Bvtx->chiSquared();
            accept = chi2_B > 0 && chi2_B < TMath::ChisquareQuantile(__PvalChi2Vtx_max__, Bvtx->degreesOfFreedom());
          }
          if(!accept) continue;
          auto B = BKinTree->currentParticle();
          auto Bvtx = BKinTree->currentDecayVertex();

          TVector3 dvtxB(Bvtx->position().x() - vtxMu->x(),
                             Bvtx->position().y() - vtxMu->y(),
                             Bvtx->position().z() - vtxMu->z()
                            );
          TVector3 pB(B->currentState().globalMomentum().x(),
                         B->currentState().globalMomentum().y(),
                         B->currentState().globalMomentum().z()
                        );
          double dalphaB = dvtxB.Angle(pB);
          auto cos_JpsiKst_vtxBest = cos(dalphaB);
          auto dvtx_JpsiKst_vtxBest = vtxu::vtxsDistance(vtxMu, Bvtx);

          // Look for the best vertex
          if (verbose) {cout << Form("Muon vtx: d = %.2f  cos = %.3e", dvtxB.Mag(), 1-cos_JpsiKst_vtxBest) << endl;}
          for(auto v : *vtxHandle) {
            if ( fabs(vtxMu->z()-v.z()) < 2 ) {
              TVector3 df(Bvtx->position().x() - v.x(),
                          Bvtx->position().y() - v.y(),
                          Bvtx->position().z() - v.z()
                          );

              double da = df.Angle(pB);
              auto cosda = cos(da);
              if (cosda > cos_JpsiKst_vtxBest) {
                dvtxB = df;
                cos_JpsiKst_vtxBest = cosda;
                dvtx_JpsiKst_vtxBest = vtxu::vtxsDistance(v, Bvtx);
              }
            }
          }
          if (verbose) {cout << Form("Best vtx: d = %.2f  cos = %.3e", dvtxB.Mag(), 1-cos_JpsiKst_vtxBest) << endl;}
          auto d_vtxJpsiKst_vtxBest = dvtx_JpsiKst_vtxBest.first;
          auto sigd_vtxJpsiKst_vtxBest = fabs(d_vtxJpsiKst_vtxBest/dvtx_JpsiKst_vtxBest.second);

          accept &= cos_JpsiKst_vtxBest > __cos_JpsiKst_vtxMu_min__;
          if(!accept) continue;

          if (verbose) {cout << "B0 -> J/psi K* found\n";}
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
          (*outputVecNtuplizer)["cos_Kpi_vtxMu"].push_back(cos_Kpi_vtxMu);
          (*outputVecNtuplizer)["d_vtxKpi_vtxMu"].push_back(d_vtxKpi_vtxMu);
          (*outputVecNtuplizer)["sigd_vtxKpi_vtxMu"].push_back(sigd_vtxKpi_vtxMu);


          JpsiKinTree->movePointerToTheFirstChild();
          auto refit_Mu1 = JpsiKinTree->currentParticle();
          JpsiKinTree->movePointerToTheNextChild();
          auto refit_Mu2 = JpsiKinTree->currentParticle();
          AddTLVToOut(vtxu::getTLVfromKinPart(refit_Mu1), string("mup"), &(*outputVecNtuplizer));
          AddTLVToOut(vtxu::getTLVfromKinPart(refit_Mu2), string("mum"), &(*outputVecNtuplizer));
          (*outputVecNtuplizer)["chi2_mumu"].push_back(vecMassMuMu[i_J]);
          (*outputVecNtuplizer)["mass_mumu"].push_back(vecChi2MuMu[i_J]);
          auto vtxJpsi = JpsiKinTree->currentDecayVertex();
          auto chi2_JpsiMass = vtxJpsi->chiSquared();
          (*outputVecNtuplizer)["chi2_Jpsi"].push_back(chi2_JpsiMass);
          (*outputVecNtuplizer)["cos_mumu_vtxMu"].push_back(cos_mumu_vtxMu);
          (*outputVecNtuplizer)["d_vtxmumu_vtxMu"].push_back(d_vtxmumu_vtxMu);
          (*outputVecNtuplizer)["sigd_vtxmumu_vtxMu"].push_back(sigd_vtxmumu_vtxMu);

          (*outputVecNtuplizer)["chi2_JpsiKst"].push_back(chi2_JpsiKst);
          (*outputVecNtuplizer)["mass_JpsiKst"].push_back(mass_JpsiKst);
          (*outputVecNtuplizer)["chi2_B"].push_back(chi2_B);
          (*outputVecNtuplizer)["cos_JpsiKsts_vtxBest"].push_back(cos_JpsiKst_vtxBest);
          (*outputVecNtuplizer)["d_vtxJpsiKst_vtxBest"].push_back(d_vtxJpsiKst_vtxBest);
          (*outputVecNtuplizer)["sigd_vtxJpsiKst_vtxBest"].push_back(sigd_vtxJpsiKst_vtxBest);

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


void B2JpsiKstDecayTreeProducer::DeclareTLVToOut(string n, map<string, vector<float>>* outv) {
  (*outv)[n+"_pt"] = {};
  (*outv)[n+"_pz"] = {};
  (*outv)[n+"_eta"] = {};
  (*outv)[n+"_phi"] = {};
  (*outv)[n+"_P"] = {};
  (*outv)[n+"_E"] = {};
  return;
}


void B2JpsiKstDecayTreeProducer::AddTLVToOut(TLorentzVector v, string n, map<string, vector<float>>* outv) {
  (*outv)[n+"_pt"].push_back(v.Pt());
  (*outv)[n+"_pz"].push_back(v.Pz());
  (*outv)[n+"_eta"].push_back(v.Eta());
  (*outv)[n+"_phi"].push_back(v.Phi());
  (*outv)[n+"_P"].push_back(v.P());
  (*outv)[n+"_E"].push_back(v.E());
  return;
}

DEFINE_FWK_MODULE(B2JpsiKstDecayTreeProducer);
