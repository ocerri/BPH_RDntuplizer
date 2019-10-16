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

#define __PvalChi2Vtx_min__ 0.05 // Very loose cut
#define __dzMax__ 2.0
#define __dmJpsi_max__ 1.0 // loose cut
#define __sigIPpfCand_min__ 2. // loose cut
#define __pThad_min__ 0.5 // loose cut
#define __dmKst_max__ 0.5 // loose cut
#define __dmB0_max__ 1.0 // loose cut

using namespace std;

class B2JpsiKstDecayTreeProducer : public edm::EDProducer {

public:

    explicit B2JpsiKstDecayTreeProducer(const edm::ParameterSet &iConfig);
    void AddTLVToOut(TLorentzVector, string, map<string, vector<float>>*);
    bool isMuonFromJpsiID(pat::Muon, reco::Vertex, pat::Muon);

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
    if(nMu < 2) {
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
    vector<double> vecMass_mumu;
    vector<vtxu::kinFitResuts> vecRes_mumu;
    vector<pair<pat::Muon, pat::Muon>> vecJpsiMuons;
    for(uint i1 = 0; i1 < nMu -1; i1++){
      auto m1 = (*muonHandle)[i1];
      if (!isMuonFromJpsiID(m1, primaryVtx, trgMu)) continue;

      for(uint i2 = i1+1; i2 < nMu; i2++){
        auto m2 = (*muonHandle)[i2];
        if (!isMuonFromJpsiID(m2, primaryVtx, trgMu)) continue;
        if(m1.charge() * m2.charge() != -1) continue;

        auto mup = m1;
        auto mum = m2;
        if(m1.charge() < 0) {mum = m1; mup = m2;}

        if(verbose){cout << "Fitting the J/Psi from: " << Form("%d %d", i1, i2) << endl;}
        auto kinTree = vtxu::FitJpsi_mumu(iSetup, mup, mum, false);
        auto res = vtxu::fitQuality(kinTree, __PvalChi2Vtx_min__);
        if(!res.isGood) continue;

        kinTree->movePointerToTheTop();
        double massMuPair = kinTree->currentParticle()->currentState().mass();
        if (fabs(massMuPair - mass_Jpsi) > __dmJpsi_max__) continue;
        if(verbose) {cout << Form("chi2: %.2f (%.2f) - pval: %.2f - mass: %.2f", res.chi2, res.dof, res.pval, massMuPair) << endl;}

        // kinTree = vtxu::FitJpsi_mumu(iSetup, mup, mum, true);
        // auto resC = vtxu::fitQuality(kinTree, __PvalChi2Vtx_min__);
        // if(!resC.isValid) continue;
        // if(verbose) {cout << Form("chi2 mass: %.2f (%.2f) - pval: %.2f", resC.chi2, resC.dof, resC.pval) << endl;}

        vecJpsiKinTree.push_back(kinTree);
        vecMass_mumu.push_back(massMuPair);
        vecRes_mumu.push_back(res);
        vecJpsiMuons.push_back(make_pair(mup, mum));
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
      auto K_sigdxy_PV = fabs(dxy)/K_tk->dxyError();
      auto K_norm_chi2 = K_tk->normalizedChi2();
      auto K_N_valid_hits = K_tk->numberOfValidHits();
      if (K_sigdxy_PV < __sigIPpfCand_min__) continue;

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
        auto pi_sigdxy_PV = fabs(dxy)/pi.dxyError();
        auto pi_norm_chi2 = pi_tk->normalizedChi2();
        auto pi_N_valid_hits = pi_tk->numberOfValidHits();
        if (pi_sigdxy_PV < __sigIPpfCand_min__) continue;

        n_pi++;
        if(verbose){cout << Form("pi cand found i:%d, pt:%.2f, eta:%.2f, phi:%.2f ", i_pi, pi.pt(), pi.eta(), pi.phi()) << endl;}

        //Fit the vertex w/o mass constraint
        auto KstKinTree = vtxu::FitKst_piK(iSetup, pi, K, false);
        auto res_piK = vtxu::fitQuality(KstKinTree, __PvalChi2Vtx_min__);
        if(!res_piK.isGood) continue;

        KstKinTree->movePointerToTheTop();
        auto mass_piK = KstKinTree->currentParticle()->currentState().mass();
        if (fabs(mass_piK - mass_Kst) > __dmKst_max__) continue;

        //Fit the vertex WITH mass constraint
        // KstKinTree = vtxu::FitKst_piK(iSetup, pi, K, true);
        // auto res_piK_cKstMass = vtxu::fitQuality(KstKinTree, __PvalChi2Vtx_min__);
        // if(!res_piK_cKstMass.isValid) continue;

        KstKinTree->movePointerToTheTop();
        auto Kst = KstKinTree->currentParticle();
        auto vtxKst = KstKinTree->currentDecayVertex();

        auto cos_Kst_PV = vtxu::computePointingCos(primaryVtx, vtxKst, Kst);
        auto cosT_Kst_PV = vtxu::computePointingCosTransverse(primaryVtx, vtxKst, Kst);
        auto d_vtxKst_PV = vtxu::vtxsDistance(primaryVtx, vtxKst);
        auto sigd_vtxKst_PV = d_vtxKst_PV.first/d_vtxKst_PV.second;
        auto dxy_vtxKst_PV = vtxu::vtxsTransverseDistance(primaryVtx, vtxKst);

        auto PhiKinTree = vtxu::FitPhi_KK(iSetup, pi, K, false);
        float mass_KK = -1;
        if(PhiKinTree->isValid()) {
          PhiKinTree->movePointerToTheTop();
          mass_KK = PhiKinTree->currentParticle()->currentState().mass();
        }

        auto KstbarKinTree = vtxu::FitKst_piK(iSetup, K, pi, false);
        float mass_piK_CPconj = -1;
        if(KstbarKinTree->isValid()) {
          KstbarKinTree->movePointerToTheTop();
          mass_piK_CPconj = KstbarKinTree->currentParticle()->currentState().mass();
        }

        if (verbose) {cout << "K*_0 -> K+ pi- found\n";}
        n_Kst++;

        /*
        ############################################################################
                  Make a B0 -> J/psi K* using the candidates found above
        ############################################################################
        */
        for(uint i_J = 0; i_J < vecJpsiKinTree.size(); ++i_J) {
          auto JpsiKinTree = vecJpsiKinTree[i_J];
          // auto res_mumu_cJpsiMass = vtxu::fitQuality(JpsiKinTree, __PvalChi2Vtx_min__);
          JpsiKinTree->movePointerToTheTop();
          auto Jpsi = JpsiKinTree->currentParticle();
          auto vtxJpsi = JpsiKinTree->currentDecayVertex();

          auto cos_Jpsi_PV = vtxu::computePointingCos(primaryVtx, vtxJpsi, Jpsi);
          auto cosT_Jpsi_PV = vtxu::computePointingCosTransverse(primaryVtx, vtxJpsi, Jpsi);
          auto d_vtxJpsi_PV = vtxu::vtxsDistance(primaryVtx, vtxJpsi);
          auto sigd_vtxJpsi_PV = fabs(d_vtxJpsi_PV.first/d_vtxJpsi_PV.second);
          auto dxy_vtxJpsi_PV = vtxu::vtxsTransverseDistance(primaryVtx, vtxJpsi);

          // auto debugPrint = [](string tag, vtxu::kinFitResuts r, double m) {
          //   cout << tag << endl;
          //   cout << "Valid: " << r.isValid << ", Good: " << r.isGood << endl;
          //   cout << Form("chi2: %.2f (%.2f), pval: %.2f, mass: %.2f", r.chi2, r.dof, r.pval, m) << endl << endl;
          // };

          auto addToOutBFit = [](string tag, vtxu::kinFitResuts r, double m, RefCountedKinematicTree kinTree, reco::Vertex primaryVtx, map<string, vector<float>>* outv) {
            (*outv)["isValid_"+tag].push_back(r.isValid);
            (*outv)["isGood_"+tag].push_back(r.isGood);
            (*outv)["chi2_"+tag].push_back(r.chi2);
            (*outv)["dof_"+tag].push_back(r.dof);
            (*outv)["pval_"+tag].push_back(r.pval);
            (*outv)["mass_"+tag].push_back(m);

            if(r.isValid) {
              kinTree->movePointerToTheTop();
              auto B = kinTree->currentParticle();
              auto vtxB = kinTree->currentDecayVertex();

              double cos_B_PV = vtxu::computePointingCos(primaryVtx, vtxB, B);
              (*outv)["cos_B_PV_"+tag].push_back(cos_B_PV);
              auto d_vtxB_PV = vtxu::vtxsDistance(primaryVtx, vtxB);
              double sigd_vtxB_PV = d_vtxB_PV.first/d_vtxB_PV.second;
              (*outv)["sigd_vtxB_PV_"+tag].push_back(sigd_vtxB_PV);

              // B2JpsiKstDecayTreeProducer::AddTLVToOut(vtxu::getTLVfromKinPart(B), string("B_"+tag), outv);
              auto pvec = B->currentState().globalMomentum();
              TLorentzVector out;
              out.SetXYZM(pvec.x(), pvec.y(), pvec.z(), m);
              (*outv)["B_"+tag+"_pt"].push_back(out.Pt());
              (*outv)["B_"+tag+"_eta"].push_back(out.Eta());
              (*outv)["B_"+tag+"_phi"].push_back(out.Phi());
            }
            else{
              (*outv)["cos_B_PV_"+tag].push_back(0);
              (*outv)["sigd_vtxB_PV_"+tag].push_back(0);
              (*outv)["B_"+tag+"_pt"].push_back(0);
              (*outv)["B_"+tag+"_eta"].push_back(0);
              (*outv)["B_"+tag+"_phi"].push_back(0);
            }
          };


          RefCountedKinematicTree BKinTree;
          BKinTree = vtxu::FitB_mumupiK(iSetup, vecJpsiMuons[i_J].first, vecJpsiMuons[i_J].second, pi, K, false, false, false);
          auto res = vtxu::fitQuality(BKinTree, __PvalChi2Vtx_min__);
          if(!res.isValid) continue;
          auto mass = BKinTree->currentParticle()->currentState().mass();
          if(fabs(mass - mass_B0) > __dmB0_max__) continue;
          // debugPrint("mumupiK", res, mass);
          addToOutBFit("mumupiK", res, mass, BKinTree, primaryVtx, &(*outputVecNtuplizer));



          BKinTree = vtxu::FitVtxJpsiKst(iSetup, Jpsi, Kst, false);
          res = vtxu::fitQuality(BKinTree, __PvalChi2Vtx_min__);
          if(res.isValid) mass = BKinTree->currentParticle()->currentState().mass();
          else mass = 0;
          // debugPrint("JpsiKst", res, mass);
          addToOutBFit("JpsiKst", res, mass, BKinTree, primaryVtx, &(*outputVecNtuplizer));

          BKinTree = vtxu::FitVtxJpsiKst(iSetup, Jpsi, Kst, true);
          res = vtxu::fitQuality(BKinTree, __PvalChi2Vtx_min__);
          if(res.isValid) mass = BKinTree->currentParticle()->currentState().mass();
          else mass = 0;
          // debugPrint("JpsiKst mass constraint", res, mass);
          addToOutBFit("JpsiKst_cBMass", res, mass, BKinTree, primaryVtx, &(*outputVecNtuplizer));

          // BKinTree = vtxu::FitB_mumupiK(iSetup, vecJpsiMuons[i_J].first, vecJpsiMuons[i_J].second, pi, K, true, false, false);
          // res = vtxu::fitQuality(BKinTree, __PvalChi2Vtx_min__);
          // if(res.isValid) mass = BKinTree->currentParticle()->currentState().mass();
          // else mass = 0;
          // // debugPrint("mumupiK JpsiKst mass", res, mass);
          // addToOutBFit("mumupiK_cJpsiKstMass", res, mass, BKinTree, primaryVtx, &(*outputVecNtuplizer));

          BKinTree = vtxu::FitB_mumupiK(iSetup, vecJpsiMuons[i_J].first, vecJpsiMuons[i_J].second, pi, K, false, true, false);
          res = vtxu::fitQuality(BKinTree, __PvalChi2Vtx_min__);
          if(res.isValid) mass = BKinTree->currentParticle()->currentState().mass();
          else mass = 0;
          // debugPrint("mumupiK B mass", res, mass);
          addToOutBFit("mumupiK_cBMass", res, mass, BKinTree, primaryVtx, &(*outputVecNtuplizer));

          BKinTree = vtxu::FitB_mumupiK(iSetup, vecJpsiMuons[i_J].first, vecJpsiMuons[i_J].second, pi, K, false, false, true, &primaryVtx);
          res = vtxu::fitQuality(BKinTree, __PvalChi2Vtx_min__);
          if(res.isValid) mass = BKinTree->currentParticle()->currentState().mass();
          else mass = 0;
          // debugPrint("mumupiK pointing", res, mass);
          addToOutBFit("mumupiK_cPointing", res, mass, BKinTree, primaryVtx, &(*outputVecNtuplizer));


          // BKinTree = vtxu::FitB_mumupiK(iSetup, vecJpsiMuons[i_J].first, vecJpsiMuons[i_J].second, pi, K, true, false, true, &primaryVtx);
          // res = vtxu::fitQuality(BKinTree, __PvalChi2Vtx_min__);
          // if(res.isValid) mass = BKinTree->currentParticle()->currentState().mass();
          // else mass = 0;
          // // debugPrint("mumupiK JpsiKst mass + pointing", res, mass);
          // addToOutBFit("mumupiK_cJpsiKstMassPointing", res, mass, BKinTree, primaryVtx, &(*outputVecNtuplizer));


          JpsiKinTree->movePointerToTheFirstChild();
          auto refit_Mu1 = JpsiKinTree->currentParticle();
          JpsiKinTree->movePointerToTheNextChild();
          auto refit_Mu2 = JpsiKinTree->currentParticle();
          KstKinTree->movePointerToTheFirstChild();
          auto refit_pi = KstKinTree->currentParticle();
          KstKinTree->movePointerToTheNextChild();
          auto refit_K = KstKinTree->currentParticle();

          BKinTree = vtxu::FitB_mumupiK(refit_Mu1, refit_Mu2, refit_pi, refit_K, false, false, false);
          res = vtxu::fitQuality(BKinTree, __PvalChi2Vtx_min__);
          if(res.isValid) mass = BKinTree->currentParticle()->currentState().mass();
          else mass = 0;
          // debugPrint("mumupiK", res, mass);
          addToOutBFit("mumupiKrefit", res, mass, BKinTree, primaryVtx, &(*outputVecNtuplizer));

          BKinTree = vtxu::FitB_mumupiK(refit_Mu1, refit_Mu2, refit_pi, refit_K, false, true, false);
          res = vtxu::fitQuality(BKinTree, __PvalChi2Vtx_min__);
          if(res.isValid) mass = BKinTree->currentParticle()->currentState().mass();
          else mass = 0;
          // debugPrint("mumupiK B mass", res, mass);
          addToOutBFit("mumupiKrefit_cBmass", res, mass, BKinTree, primaryVtx, &(*outputVecNtuplizer));

          n_B++;

          /*
          ############################################################################
                                Compute analysis variables
          ############################################################################
          */
          AddTLVToOut(vtxu::getTLVfromCand(pi, mass_Pi), string("pi"), &(*outputVecNtuplizer));
          (*outputVecNtuplizer)["pi_sigdxy_PV"].push_back(pi_sigdxy_PV);
          (*outputVecNtuplizer)["pi_norm_chi2"].push_back(pi_norm_chi2);
          (*outputVecNtuplizer)["pi_N_valid_hits"].push_back(pi_N_valid_hits);
          AddTLVToOut(vtxu::getTLVfromCand(K, mass_K), string("K"), &(*outputVecNtuplizer));
          (*outputVecNtuplizer)["K_sigdxy_PV"].push_back(K_sigdxy_PV);
          (*outputVecNtuplizer)["K_norm_chi2"].push_back(K_norm_chi2);
          (*outputVecNtuplizer)["K_N_valid_hits"].push_back(K_N_valid_hits);

          (*outputVecNtuplizer)["chi2_piK"].push_back(res_piK.chi2);
          (*outputVecNtuplizer)["dof_piK"].push_back(res_piK.dof);
          (*outputVecNtuplizer)["pval_piK"].push_back(res_piK.pval);
          (*outputVecNtuplizer)["mass_piK"].push_back(mass_piK);
          // (*outputVecNtuplizer)["chi2_piK_cKstMass"].push_back(res_piK_cKstMass.chi2);
          // (*outputVecNtuplizer)["dof_piK_cKstMass"].push_back(res_piK_cKstMass.dof);
          // (*outputVecNtuplizer)["pval_piK_cKstMass"].push_back(res_piK_cKstMass.pval);
          AddTLVToOut(vtxu::getTLVfromKinPart(Kst), string("Kst"), &(*outputVecNtuplizer));
          // KstKinTree->movePointerToTheFirstChild();
          // auto refit_pi = KstKinTree->currentParticle();
          AddTLVToOut(vtxu::getTLVfromKinPart(refit_pi), string("piRefit"), &(*outputVecNtuplizer));
          // KstKinTree->movePointerToTheNextChild();
          // auto refit_K = KstKinTree->currentParticle();
          AddTLVToOut(vtxu::getTLVfromKinPart(refit_K), string("KRefit"), &(*outputVecNtuplizer));
          (*outputVecNtuplizer)["cos_Kst_PV"].push_back(cos_Kst_PV);
          (*outputVecNtuplizer)["cosT_Kst_PV"].push_back(cosT_Kst_PV);
          (*outputVecNtuplizer)["d_vtxKst_PV"].push_back(d_vtxKst_PV.first);
          (*outputVecNtuplizer)["sigd_vtxKst_PV"].push_back(sigd_vtxKst_PV);
          (*outputVecNtuplizer)["sigdxy_vtxKst_PV"].push_back(dxy_vtxKst_PV.first/dxy_vtxKst_PV.second);
          (*outputVecNtuplizer)["mass_KK"].push_back(mass_KK);
          (*outputVecNtuplizer)["mass_piK_CPconj"].push_back(mass_piK_CPconj);

          auto mup = vecJpsiMuons[i_J].first;
          AddTLVToOut(vtxu::getTLVfromTrack(*(mup.bestTrack()), mass_Mu), string("mup"), &(*outputVecNtuplizer));
          auto dxy_mup = mup.innerTrack()->dxy(primaryVtx.position());
          (*outputVecNtuplizer)["mup_dxy"].push_back(dxy_mup);
          (*outputVecNtuplizer)["mup_sigdxy_PV"].push_back(fabs(dxy_mup)/mup.innerTrack()->dxyError());
          float mup_isTrg = trgMu.pt()==mup.pt() && trgMu.phi()==mup.phi() && trgMu.eta()==mup.eta();
          (*outputVecNtuplizer)["mup_isTrg"].push_back(mup_isTrg);

          auto mum = vecJpsiMuons[i_J].second;
          AddTLVToOut(vtxu::getTLVfromTrack(*(mum.bestTrack()), mass_Mu), string("mum"), &(*outputVecNtuplizer));
          auto dxy_mum = mum.innerTrack()->dxy(primaryVtx.position());
          (*outputVecNtuplizer)["mum_dxy"].push_back(dxy_mum);
          (*outputVecNtuplizer)["mum_sigdxy_PV"].push_back(fabs(dxy_mum)/mum.innerTrack()->dxyError());
          float mum_isTrg = trgMu.pt()==mum.pt() && trgMu.phi()==mum.phi() && trgMu.eta()==mum.eta();
          (*outputVecNtuplizer)["mum_isTrg"].push_back(mum_isTrg);

          (*outputVecNtuplizer)["chi2_mumu"].push_back(vecRes_mumu[i_J].chi2);
          (*outputVecNtuplizer)["dof_mumu"].push_back(vecRes_mumu[i_J].dof);
          (*outputVecNtuplizer)["pval_mumu"].push_back(vecRes_mumu[i_J].pval);
          (*outputVecNtuplizer)["mass_mumu"].push_back(vecMass_mumu[i_J]);
          // (*outputVecNtuplizer)["chi2_mumu_cJpsiMass"].push_back(res_mumu_cJpsiMass.chi2);
          // (*outputVecNtuplizer)["dof_mumu_cJpsiMass"].push_back(res_mumu_cJpsiMass.dof);
          // (*outputVecNtuplizer)["pval_mumu_cJpsiMass"].push_back(res_mumu_cJpsiMass.pval);
          AddTLVToOut(vtxu::getTLVfromKinPart(Jpsi), string("Jpsi"), &(*outputVecNtuplizer));
          // JpsiKinTree->movePointerToTheFirstChild();
          // auto refit_Mu1 = JpsiKinTree->currentParticle();
          AddTLVToOut(vtxu::getTLVfromKinPart(refit_Mu1), string("mupRefit"), &(*outputVecNtuplizer));
          // JpsiKinTree->movePointerToTheNextChild();
          // auto refit_Mu2 = JpsiKinTree->currentParticle();
          AddTLVToOut(vtxu::getTLVfromKinPart(refit_Mu2), string("mumRefit"), &(*outputVecNtuplizer));
          (*outputVecNtuplizer)["cos_Jpsi_PV"].push_back(cos_Jpsi_PV);
          (*outputVecNtuplizer)["cosT_Jpsi_PV"].push_back(cosT_Jpsi_PV);
          (*outputVecNtuplizer)["d_vtxJpsi_PV"].push_back(d_vtxJpsi_PV.first);
          (*outputVecNtuplizer)["sigd_vtxJpsi_PV"].push_back(sigd_vtxJpsi_PV);
          (*outputVecNtuplizer)["sigdxy_vtxJpsi_PV"].push_back(dxy_vtxJpsi_PV.first/dxy_vtxJpsi_PV.second);

          // AddTLVToOut(vtxu::getTLVfromKinPart(B), string("B"), &(*outputVecNtuplizer));
          // BKinTree->movePointerToTheFirstChild();
          // auto refit_Jpsi = BKinTree->currentParticle();
          // AddTLVToOut(vtxu::getTLVfromKinPart(refit_Jpsi), string("Jpsi"), &(*outputVecNtuplizer));
          // BKinTree->movePointerToTheNextChild();
          // auto refit_Kst = BKinTree->currentParticle();
          // AddTLVToOut(vtxu::getTLVfromKinPart(refit_Kst), string("Kst"), &(*outputVecNtuplizer));
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


void B2JpsiKstDecayTreeProducer::AddTLVToOut(TLorentzVector v, string n, map<string, vector<float>>* outv) {
  (*outv)[n+"_pt"].push_back(v.Pt());
  // (*outv)[n+"_pz"].push_back(v.Pz());
  (*outv)[n+"_eta"].push_back(v.Eta());
  (*outv)[n+"_phi"].push_back(v.Phi());
  // (*outv)[n+"_P"].push_back(v.P());
  // (*outv)[n+"_E"].push_back(v.E());
  return;
}


bool B2JpsiKstDecayTreeProducer::isMuonFromJpsiID(pat::Muon m, reco::Vertex pVtx, pat::Muon trgMu) {
  if(m.innerTrack().isNull()) return false;
  if (fabs(m.innerTrack()->dz(pVtx.position()) - trgMu.innerTrack()->dz(pVtx.position())) > __dzMax__) return false;

  if(m.innerTrack()->hitPattern().pixelLayersWithMeasurement() < 2) return false;
  if(!m.innerTrack()->quality(reco::TrackBase::highPurity)) return false;
  if(!m.isGood("TMOneStationTight")) return false;
  if(m.innerTrack()->normalizedChi2() > 1.8) return false;

  double dxy = m.innerTrack()->dxy(pVtx.position());
  float sigdxy = fabs(dxy)/m.innerTrack()->dxyError();
  if (sigdxy < 2) return false;

  return true;
}

DEFINE_FWK_MODULE(B2JpsiKstDecayTreeProducer);
