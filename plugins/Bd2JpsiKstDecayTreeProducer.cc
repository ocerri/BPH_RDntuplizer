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
#define __dzMax__ 3.0
#define __dmJpsi_max__ 0.4 // loose cut
#define __sigIPpfCand_min__ 2. // loose cut
#define __pThad_min__ 0.5 // loose cut
#define __dmKst_max__ 0.3 // loose cut
#define __dmB0_max__ 0.5 // loose cut

using namespace std;

class Bd2JpsiKstDecayTreeProducer : public edm::EDProducer {

public:

    explicit Bd2JpsiKstDecayTreeProducer(const edm::ParameterSet &iConfig);
    void AddTLVToOut(TLorentzVector, string, map<string, vector<float>>*);
    bool isMuonFromJpsiID(pat::Muon, reco::Vertex, pat::Muon);
    int trgMuIdx(pat::Muon, vector<pat::Muon>);

    ~Bd2JpsiKstDecayTreeProducer() override {};

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



Bd2JpsiKstDecayTreeProducer::Bd2JpsiKstDecayTreeProducer(const edm::ParameterSet &iConfig)
{
    verbose = iConfig.getParameter<int>( "verbose" );
    PFCandSrc_ = consumes<vector<pat::PackedCandidate>>(edm::InputTag("packedPFCandidates"));
    muonSrc_ = consumes<vector<pat::Muon>>(edm::InputTag("slimmedMuons"));
    vtxSrc_ = consumes<vector<reco::Vertex>>(edm::InputTag("offlineSlimmedPrimaryVertices"));
    TrgMuonSrc_ = consumes<vector<pat::Muon>> ( iConfig.getParameter<edm::InputTag>( "trgMuons" ) );

    produces<map<string, float>>("outputNtuplizer");
    produces<map<string, vector<float>>>("outputVecNtuplizer");
}


void Bd2JpsiKstDecayTreeProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
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

    vector<RefCountedKinematicTree> vecJpsiKinTree;
    vector<double> vecMass_mumu;
    vector<vtxu::kinFitResuts> vecRes_mumu;
    vector<pair<pat::Muon, pat::Muon>> vecJpsiMuons;
    for(uint i1 = 0; i1 < nMu -1; i1++){
      auto m1 = (*muonHandle)[i1];
      if (m1.innerTrack().isNull()) continue;
      // if (!m1.isSoftMuon(primaryVtx)) continue;
      if (!m1.isMediumMuon()) continue;

      for(uint i2 = i1+1; i2 < nMu; i2++){
        auto m2 = (*muonHandle)[i2];
        if (m2.innerTrack().isNull()) continue;
        // if (!m2.isSoftMuon(primaryVtx)) continue;
        if (!m2.isMediumMuon()) continue;
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
      //Require a charged hadron
      if (K.charge() == 0) continue;
      if (!K.hasTrackDetails()) continue;
      if (K.isTrackerMuon() || K.isStandAloneMuon()) continue;
      //Require a minimum pt
      if(K.pt() < __pThad_min__) continue;
      // Require to be close to the trigger muon;
      auto K_tk = K.bestTrack();
      // if (fabs(K_tk->dz(primaryVtx.position()) - trgMu.muonBestTrack()->dz(primaryVtx.position())) > __dzMax__) continue;
      if (fabs( K_tk->dz(primaryVtx.position()) ) > __dzMax__) continue;
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
        if (pi.charge() + K.charge() != 0) continue;
        //Require a minimum pt
        if(pi.pt() < __pThad_min__) continue;
        // Require to be close to the trigger muon;
        auto pi_tk = pi.bestTrack();
        // if (fabs(pi_tk->dz(primaryVtx.position()) - trgMu.muonBestTrack()->dz(primaryVtx.position())) > __dzMax__) continue;
        if (fabs( pi_tk->dz(primaryVtx.position()) ) > __dzMax__) continue;
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

        KstKinTree->movePointerToTheTop();
        auto Kst = KstKinTree->currentParticle();
        auto vtxKst = KstKinTree->currentDecayVertex();

        auto cos_Kst_PV = vtxu::computePointingCos(primaryVtx, vtxKst, Kst);
        auto cosT_Kst_PV = vtxu::computePointingCosTransverse(primaryVtx, vtxKst, Kst);
        auto d_vtxKst_PV = vtxu::vtxsDistance(primaryVtx, vtxKst);
        auto sigd_vtxKst_PV = d_vtxKst_PV.first/d_vtxKst_PV.second;
        auto dxy_vtxKst_PV = vtxu::vtxsTransverseDistance(primaryVtx, vtxKst);
        auto sigdxy_vtxKst_PV = fabs(dxy_vtxKst_PV.first)/dxy_vtxKst_PV.second;
        if (sigdxy_vtxKst_PV < 3) continue;

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

              (*outv)["vtx_PV_x"].push_back(primaryVtx.position().x());
              (*outv)["vtx_PV_y"].push_back(primaryVtx.position().y());
              (*outv)["vtx_PV_z"].push_back(primaryVtx.position().z());

              (*outv)["vtx_B_decay_x"].push_back(vtxB->position().x());
              (*outv)["vtx_B_decay_y"].push_back(vtxB->position().y());
              (*outv)["vtx_B_decay_z"].push_back(vtxB->position().z());

              double cos_B_PV = vtxu::computePointingCos(primaryVtx, vtxB, B);
              (*outv)["cos_B_PV_"+tag].push_back(cos_B_PV);
              auto cosT_B_PV = vtxu::computePointingCosTransverse(primaryVtx, vtxB, B);
              (*outv)["cosT_B_PV_"+tag].push_back(cosT_B_PV);
              auto d_vtxB_PV = vtxu::vtxsDistance(primaryVtx, vtxB);
              double sigd_vtxB_PV = d_vtxB_PV.first/d_vtxB_PV.second;
              (*outv)["d_vtxB_PV_"+tag].push_back(d_vtxB_PV.first);
              (*outv)["sigd_vtxB_PV_"+tag].push_back(sigd_vtxB_PV);
              auto dxy_vtxB_PV = vtxu::vtxsTransverseDistance(primaryVtx, vtxB);
              (*outv)["dxy_vtxB_PV"].push_back(dxy_vtxB_PV.first);
              (*outv)["sigdxy_vtxB_PV"].push_back(dxy_vtxB_PV.first/dxy_vtxB_PV.second);

              // Bd2JpsiKstDecayTreeProducer::AddTLVToOut(vtxu::getTLVfromKinPart(B), string("B_"+tag), outv);
              auto pvec = B->currentState().globalMomentum();
              TLorentzVector out;
              out.SetXYZM(pvec.x(), pvec.y(), pvec.z(), m);
              (*outv)["B_"+tag+"_pt"].push_back(out.Pt());
              (*outv)["B_"+tag+"_eta"].push_back(out.Eta());
              (*outv)["B_"+tag+"_phi"].push_back(out.Phi());
            }
            else{
              (*outv)["vtx_PV_x"].push_back(-999);
              (*outv)["vtx_PV_y"].push_back(-999);
              (*outv)["vtx_PV_z"].push_back(-999);
              (*outv)["vtx_B_decay_x"].push_back(-999);
              (*outv)["vtx_B_decay_y"].push_back(-999);
              (*outv)["vtx_B_decay_z"].push_back(-999);
              (*outv)["cos_B_PV_"+tag].push_back(0);
              (*outv)["cosT_B_PV_"+tag].push_back(0);
              (*outv)["d_vtxB_PV_"+tag].push_back(0);
              (*outv)["sigd_vtxB_PV_"+tag].push_back(0);
              (*outv)["dxy_vtxB_PV_"+tag].push_back(0);
              (*outv)["sigdxy_vtxB_PV_"+tag].push_back(0);
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

          BKinTree->movePointerToTheTop();
          auto B = BKinTree->currentParticle();
          auto vtxB = BKinTree->currentDecayVertex();
          // Looking for the best vertex to associate it with
          uint i_best = 0;
          auto maxCos = vtxu::computePointingCos(primaryVtx, vtxB, B);
          for(uint i_vtx = 1; i_vtx < vtxHandle->size(); i_vtx++) {
            if ((*vtxHandle)[i_vtx].ndof() < 4) continue;
            auto auxCos = vtxu::computePointingCos((*vtxHandle)[i_vtx], vtxB, B);
            if(auxCos > maxCos) {
              maxCos = auxCos;
              i_best = i_vtx;
            }
          }
          auto bestVtx = (*vtxHandle)[i_best];

          addToOutBFit("mumupiK", res, mass, BKinTree, bestVtx, &(*outputVecNtuplizer));

          BKinTree->movePointerToTheFirstChild();
          auto refit_Mu1 = BKinTree->currentParticle();
          auto p4_refit_Mu1 = vtxu::getTLVfromKinPart(refit_Mu1);
          BKinTree->movePointerToTheNextChild();
          auto refit_Mu2 = BKinTree->currentParticle();
          auto p4_refit_Mu2 = vtxu::getTLVfromKinPart(refit_Mu2);
          BKinTree->movePointerToTheNextChild();
          auto refit_pi = BKinTree->currentParticle();
          auto p4_refit_pi = vtxu::getTLVfromKinPart(refit_pi);
          BKinTree->movePointerToTheNextChild();
          auto refit_K = BKinTree->currentParticle();
          auto p4_refit_K = vtxu::getTLVfromKinPart(refit_K);

          n_B++;

          /*
          ############################################################################
                                Compute analysis variables
          ############################################################################
          */
          AddTLVToOut(vtxu::getTLVfromCand(pi, mass_Pi), string("pi"), &(*outputVecNtuplizer));
          (*outputVecNtuplizer)["pi_charge"].push_back(pi.charge());
          (*outputVecNtuplizer)["pi_sigdxy_PV"].push_back(pi_sigdxy_PV);
          (*outputVecNtuplizer)["pi_norm_chi2"].push_back(pi_norm_chi2);
          (*outputVecNtuplizer)["pi_N_valid_hits"].push_back(pi_N_valid_hits);
          AddTLVToOut(vtxu::getTLVfromCand(K, mass_K), string("K"), &(*outputVecNtuplizer));
          (*outputVecNtuplizer)["K_charge"].push_back(K.charge());
          (*outputVecNtuplizer)["K_sigdxy_PV"].push_back(K_sigdxy_PV);
          (*outputVecNtuplizer)["K_norm_chi2"].push_back(K_norm_chi2);
          (*outputVecNtuplizer)["K_N_valid_hits"].push_back(K_N_valid_hits);

          (*outputVecNtuplizer)["chi2_piK"].push_back(res_piK.chi2);
          (*outputVecNtuplizer)["dof_piK"].push_back(res_piK.dof);
          (*outputVecNtuplizer)["pval_piK"].push_back(res_piK.pval);
          (*outputVecNtuplizer)["mass_piK"].push_back(mass_piK);
          (*outputVecNtuplizer)["cos_Kst_PV"].push_back(cos_Kst_PV);
          (*outputVecNtuplizer)["cosT_Kst_PV"].push_back(cosT_Kst_PV);
          (*outputVecNtuplizer)["d_vtxKst_PV"].push_back(d_vtxKst_PV.first);
          (*outputVecNtuplizer)["sigd_vtxKst_PV"].push_back(sigd_vtxKst_PV);
          (*outputVecNtuplizer)["sigdxy_vtxKst_PV"].push_back(sigdxy_vtxKst_PV);
          (*outputVecNtuplizer)["mass_KK"].push_back(mass_KK);
          (*outputVecNtuplizer)["mass_piK_CPconj"].push_back(mass_piK_CPconj);

          auto mup = vecJpsiMuons[i_J].first;
          AddTLVToOut(vtxu::getTLVfromTrack(*(mup.bestTrack()), mass_Mu), string("mup"), &(*outputVecNtuplizer));
          auto dxy_mup = mup.innerTrack()->dxy(primaryVtx.position());
          (*outputVecNtuplizer)["mup_dxy_PV"].push_back(dxy_mup);
          (*outputVecNtuplizer)["mup_sigdxy_PV"].push_back(fabs(dxy_mup)/mup.innerTrack()->dxyError());
          (*outputVecNtuplizer)["mup_isTrg"].push_back(trgMuIdx(mup, (*trgMuonsHandle) ));

          auto mum = vecJpsiMuons[i_J].second;
          AddTLVToOut(vtxu::getTLVfromTrack(*(mum.bestTrack()), mass_Mu), string("mum"), &(*outputVecNtuplizer));
          auto dxy_mum = mum.innerTrack()->dxy(primaryVtx.position());
          (*outputVecNtuplizer)["mum_dxy_PV"].push_back(dxy_mum);
          (*outputVecNtuplizer)["mum_sigdxy_PV"].push_back(fabs(dxy_mum)/mum.innerTrack()->dxyError());
          (*outputVecNtuplizer)["mum_isTrg"].push_back(trgMuIdx(mum, (*trgMuonsHandle) ));

          (*outputVecNtuplizer)["chi2_mumu"].push_back(vecRes_mumu[i_J].chi2);
          (*outputVecNtuplizer)["dof_mumu"].push_back(vecRes_mumu[i_J].dof);
          (*outputVecNtuplizer)["pval_mumu"].push_back(vecRes_mumu[i_J].pval);
          (*outputVecNtuplizer)["mass_mumu"].push_back(vecMass_mumu[i_J]);
          // (*outputVecNtuplizer)["chi2_mumu_cJpsiMass"].push_back(res_mumu_cJpsiMass.chi2);
          // (*outputVecNtuplizer)["dof_mumu_cJpsiMass"].push_back(res_mumu_cJpsiMass.dof);
          // (*outputVecNtuplizer)["pval_mumu_cJpsiMass"].push_back(res_mumu_cJpsiMass.pval);
          (*outputVecNtuplizer)["cos_Jpsi_PV"].push_back(cos_Jpsi_PV);
          (*outputVecNtuplizer)["cosT_Jpsi_PV"].push_back(cosT_Jpsi_PV);
          (*outputVecNtuplizer)["d_vtxJpsi_PV"].push_back(d_vtxJpsi_PV.first);
          (*outputVecNtuplizer)["sigd_vtxJpsi_PV"].push_back(sigd_vtxJpsi_PV);
          (*outputVecNtuplizer)["sigdxy_vtxJpsi_PV"].push_back(dxy_vtxJpsi_PV.first/dxy_vtxJpsi_PV.second);

          AddTLVToOut(p4_refit_Mu1, string("mupRefit"), &(*outputVecNtuplizer));
          AddTLVToOut(p4_refit_Mu2, string("mumRefit"), &(*outputVecNtuplizer));
          AddTLVToOut(p4_refit_Mu1+p4_refit_Mu1, string("JpsiRefit"), &(*outputVecNtuplizer));
          AddTLVToOut(p4_refit_pi, string("piRefit"), &(*outputVecNtuplizer));
          AddTLVToOut(p4_refit_K, string("KRefit"), &(*outputVecNtuplizer));
          AddTLVToOut(p4_refit_pi+p4_refit_K, string("KstRefit"), &(*outputVecNtuplizer));

          /*
          ############################################################################
                  Study the presence of additional tracks in the D* mu vertex
          ############################################################################
          */
          uint N_compatible_tk = 0;
          vector<double> tksAdd_massHad = {};
          vector<double> tksAdd_massVis = {};
          vector<double> tksAdd_pval = {};
          vector<double> tksAdd_charge = {};
          vector<double> tksAdd_pdgId = {};
          vector<double> tksAdd_pt = {};
          vector<double> tksAdd_eta = {};
          vector<double> tksAdd_phi = {};
          vector<double> tksAdd_dz = {};
          vector<double> tksAdd_lostInnerHits = {};
          vector<double> tksAdd_sigdca_vtxB = {};
          vector<double> tksAdd_cos_PV = {};
          if(verbose) {cout << "Looking for additional tracks compatible with D0pismu vertex" << endl;}

          double dR_max = 1.;
          if(verbose) {cout << Form("Looking for additional neutrals within dR=%.3f from muon", dR_max) << endl;}
          uint N_compatible_neu = 0;
          vector<double> neu_pdgId = {};
          vector<double> neu_pt = {};
          vector<double> neu_eta = {};
          vector<double> neu_phi = {};
          vector<double> neu_energy = {};
          vector<double> neu_et2 = {};
          vector<double> neu_dR_fromMu = {};

          int i_mup_PFcand = -1;
          int i_mum_PFcand = -1;

          for(uint i_tk = 0; i_tk < N_pfCand; ++i_tk) {
            // PF candidate different from K, pi and pis
            if( i_tk==i_k || i_tk==i_pi) continue;

            const pat::PackedCandidate & ptk = (*pfCandHandle)[i_tk];

            // PF candidate not matching with the muon
            double dR_fromMup = hypot(ptk.eta() - mup.eta(), ptk.phi() - mup.phi());
            if ( fabs(ptk.pt() - mup.pt())/mup.pt() < 0.01  && dR_fromMup < 0.01 ) {
              if (ptk.pdgId() == mup.pdgId()) i_mup_PFcand = i_tk;
              continue;
            }
            double dR_fromMum = hypot(ptk.eta() - mum.eta(), ptk.phi() - mum.phi());
            if ( fabs(ptk.pt() - mum.pt())/mum.pt() < 0.01  && dR_fromMum < 0.01 ) {
              if (ptk.pdgId() == mum.pdgId()) i_mum_PFcand = i_tk;
              continue;
            }

            if (ptk.charge() != 0 ) {
              if (!ptk.hasTrackDetails()) continue;
              if (ptk.pt() < 0.3) continue;
              // Require to be close to the trigger muon;
              auto tk = ptk.bestTrack();
              if (fabs(tk->dz(primaryVtx.position()) - mum.muonBestTrack()->dz(primaryVtx.position())) > __dzMax__) continue;
              if (fabs(tk->dz(primaryVtx.position()) - mup.muonBestTrack()->dz(primaryVtx.position())) > __dzMax__) continue;
              if (vtxu::dR(ptk.phi(), mum.phi(), ptk.eta(), mum.eta()) > 3) continue;
              if (vtxu::dR(ptk.phi(), mup.phi(), ptk.eta(), mup.eta()) > 3) continue;

              GlobalPoint auxp(vtxB->position().x(), vtxB->position().y(), vtxB->position().z());
              auto dca = vtxu::computeDCA(iSetup, ptk, auxp);

              // Check if it makes a consistent vertex
              auto kinTree = vtxu::FitB_mumupiK_tk(iSetup, mup, mum, pi, K, ptk, false, false, false);
              auto res = vtxu::fitQuality(kinTree, __PvalChi2Vtx_min__);
              if(!res.isGood) continue;
              if (verbose) {cout << Form("Trk pars: pt=%.2f eta=%.2f phi=%.2f\n", ptk.pt(), ptk.eta(), ptk.phi());}

              kinTree->movePointerToTheTop();
              auto p_vis = kinTree->currentParticle();
              auto m_vis = p_vis->currentState().mass();

              kinTree->movePointerToTheFirstChild(); // mu+
              kinTree->movePointerToTheNextChild();  // mu-
              kinTree->movePointerToTheNextChild();  // pi
              auto refit_pi = kinTree->currentParticle();
              auto refit_pi_p4 = vtxu::getTLVfromKinPart(refit_pi);
              kinTree->movePointerToTheNextChild();  // K
              auto refit_K = kinTree->currentParticle();
              auto refit_K_p4 = vtxu::getTLVfromKinPart(refit_K);
              kinTree->movePointerToTheNextChild();  // addTk
              auto refit_tk = kinTree->currentParticle();
              auto refit_tk_p4 = vtxu::getTLVfromKinPart(refit_tk);

              auto m_had = (refit_pi_p4 + refit_K_p4 + refit_tk_p4).M();

              tksAdd_massVis.push_back(m_vis);
              tksAdd_massHad.push_back(m_had);
              tksAdd_pval.push_back(res.pval);
              tksAdd_charge.push_back(ptk.charge());
              tksAdd_pdgId.push_back(ptk.pdgId());
              tksAdd_pt.push_back(refit_tk_p4.Pt());
              tksAdd_eta.push_back(refit_tk_p4.Eta());
              tksAdd_phi.push_back(refit_tk_p4.Phi());
              tksAdd_dz.push_back(tk->dz(bestVtx.position()));
              tksAdd_lostInnerHits.push_back(ptk.lostInnerHits());
              tksAdd_sigdca_vtxB.push_back(fabs(dca.first)/dca.second);
              tksAdd_cos_PV.push_back(vtxu::computePointingCos(bestVtx, vtxB, refit_tk));
              N_compatible_tk++;
            }
            else {
              if (min(dR_fromMum, dR_fromMup) > dR_max) continue;
              if (ptk.energy() < 0.5) continue;

              neu_pdgId.push_back(ptk.pdgId());
              neu_pt.push_back(ptk.pt());
              neu_eta.push_back(ptk.eta());
              neu_phi.push_back(ptk.phi());
              neu_energy.push_back(ptk.energy());
              neu_et2.push_back(ptk.et2());
              neu_dR_fromMu.push_back(min(dR_fromMup, dR_fromMum));

              N_compatible_neu++;
            }
          }

          if ( (*outputVecNtuplizer)["nTksAdd"].size() == 0 ) {
            (*outputVecNtuplizer)["tksAdd_massVis"] = {};
            (*outputVecNtuplizer)["tksAdd_massHad"] = {};
            (*outputVecNtuplizer)["tksAdd_pval"] = {};
            (*outputVecNtuplizer)["tksAdd_charge"] = {};
            (*outputVecNtuplizer)["tksAdd_pdgId"] = {};
            (*outputVecNtuplizer)["tksAdd_pt"] = {};
            (*outputVecNtuplizer)["tksAdd_eta"] = {};
            (*outputVecNtuplizer)["tksAdd_phi"] = {};
            (*outputVecNtuplizer)["tksAdd_dz"] = {};
            (*outputVecNtuplizer)["tksAdd_lostInnerHits"] = {};
            (*outputVecNtuplizer)["tksAdd_sigdca_vtxB"] = {};
            (*outputVecNtuplizer)["tksAdd_cos_PV"] = {};
          }
          (*outputVecNtuplizer)["nTksAdd"].push_back(N_compatible_tk);
          if(verbose) {cout << "Compatible tracks: " << N_compatible_tk << endl;}
          for(uint i = 0; i < N_compatible_tk; i++){
            (*outputVecNtuplizer)["tksAdd_massVis"].push_back(tksAdd_massVis[i]);
            (*outputVecNtuplizer)["tksAdd_massHad"].push_back(tksAdd_massHad[i]);
            (*outputVecNtuplizer)["tksAdd_pval"].push_back(tksAdd_pval[i]);
            (*outputVecNtuplizer)["tksAdd_charge"].push_back(tksAdd_charge[i]);
            (*outputVecNtuplizer)["tksAdd_pdgId"].push_back(tksAdd_pdgId[i]);
            (*outputVecNtuplizer)["tksAdd_pt"].push_back(tksAdd_pt[i]);
            (*outputVecNtuplizer)["tksAdd_eta"].push_back(tksAdd_eta[i]);
            (*outputVecNtuplizer)["tksAdd_phi"].push_back(tksAdd_phi[i]);
            (*outputVecNtuplizer)["tksAdd_dz"].push_back(tksAdd_dz[i]);
            (*outputVecNtuplizer)["tksAdd_lostInnerHits"].push_back(tksAdd_lostInnerHits[i]);
            (*outputVecNtuplizer)["tksAdd_sigdca_vtxB"].push_back(tksAdd_sigdca_vtxB[i]);
            (*outputVecNtuplizer)["tksAdd_cos_PV"].push_back(tksAdd_cos_PV[i]);
          }

          int auxC = 0;
          for(auto auxN : (*outputVecNtuplizer)["nTksAdd"]) auxC += auxN;
          if( auxC != int ((*outputVecNtuplizer)["tksAdd_massVis"].size()) ){
            cout << "Number of tracks and tracks details lenght not matching" << endl;
            assert(false);
          }


          if ( (*outputVecNtuplizer)["nNeuAdd"].size() == 0 ) {
            (*outputVecNtuplizer)["neuAdd_pdgId"] = {};
            (*outputVecNtuplizer)["neuAdd_pt"] = {};
            (*outputVecNtuplizer)["neuAdd_eta"] = {};
            (*outputVecNtuplizer)["neuAdd_phi"] = {};
            (*outputVecNtuplizer)["neuAdd_energy"] = {};
            (*outputVecNtuplizer)["neuAdd_et2"] = {};
            (*outputVecNtuplizer)["neuAdd_dR_fromMu"] = {};
          }
          (*outputVecNtuplizer)["nNeuAdd"].push_back(N_compatible_neu);
          if(verbose) {cout << "Compatible neutrals: " << N_compatible_neu << endl;}
          for(uint i = 0; i < N_compatible_neu; i++){
            (*outputVecNtuplizer)["neuAdd_pdgId"].push_back(neu_pdgId[i]);
            (*outputVecNtuplizer)["neuAdd_pt"].push_back(neu_pt[i]);
            (*outputVecNtuplizer)["neuAdd_eta"].push_back(neu_eta[i]);
            (*outputVecNtuplizer)["neuAdd_phi"].push_back(neu_phi[i]);
            (*outputVecNtuplizer)["neuAdd_energy"].push_back(neu_energy[i]);
            (*outputVecNtuplizer)["neuAdd_et2"].push_back(neu_et2[i]);
            (*outputVecNtuplizer)["neuAdd_dR_fromMu"].push_back(neu_dR_fromMu[i]);
          }

          auxC = 0;
          for(auto auxN : (*outputVecNtuplizer)["nNeuAdd"]) auxC += auxN;
          if( auxC != int ((*outputVecNtuplizer)["neuAdd_pdgId"].size()) ){
            cout << "Number of neutrals and neutrals details lenght not matching" << endl;
            assert(false);
          }

          if(n_B >= 100) break;
        }
        if(n_B >= 100) break;
      }
      if(n_B >= 100) break;
    }

    (*outputNtuplizer)["n_K"] = n_K;
    (*outputNtuplizer)["n_pi"] = n_pi;
    (*outputNtuplizer)["n_Kst"] = n_Kst;
    (*outputNtuplizer)["n_Jpsi"] = vecJpsiKinTree.size();
    (*outputNtuplizer)["n_B"] = n_B;


    iEvent.put(move(outputNtuplizer), "outputNtuplizer");
    iEvent.put(move(outputVecNtuplizer), "outputVecNtuplizer");
    return;
}


void Bd2JpsiKstDecayTreeProducer::AddTLVToOut(TLorentzVector v, string n, map<string, vector<float>>* outv) {
  (*outv)[n+"_pt"].push_back(v.Pt());
  // (*outv)[n+"_pz"].push_back(v.Pz());
  (*outv)[n+"_eta"].push_back(v.Eta());
  (*outv)[n+"_phi"].push_back(v.Phi());
  // (*outv)[n+"_P"].push_back(v.P());
  // (*outv)[n+"_E"].push_back(v.E());
  return;
}

int Bd2JpsiKstDecayTreeProducer::trgMuIdx(pat::Muon m, vector<pat::Muon> trgMuons) {
  int idx = -1;
  int Nmu = trgMuons.size();
  for(int i = 0; i < Nmu; i++) {
    bool match = trgMuons[i].pt() == m.pt();
    match &= trgMuons[i].eta() == m.eta();
    match &= trgMuons[i].phi() == m.phi();
    match &= trgMuons[i].charge() == m.charge();
    if(match) {
      idx = i;
      break;
    }
  }
  return idx;
}

bool Bd2JpsiKstDecayTreeProducer::isMuonFromJpsiID(pat::Muon m, reco::Vertex pVtx, pat::Muon trgMu) {
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

DEFINE_FWK_MODULE(Bd2JpsiKstDecayTreeProducer);
