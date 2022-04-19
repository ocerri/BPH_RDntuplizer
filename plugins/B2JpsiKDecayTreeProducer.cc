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
#define __dmB0_max__ 0.4 // loose cut

using namespace std;

class B2JpsiKDecayTreeProducer : public edm::EDProducer {

public:

    explicit B2JpsiKDecayTreeProducer(const edm::ParameterSet &iConfig);
    void AddTLVToOut(TLorentzVector, string, map<string, vector<float>>*);
    bool isMuonFromJpsiID(pat::Muon, reco::Vertex, pat::Muon);
    int trgMuIdx(pat::Muon, vector<pat::Muon>);

    ~B2JpsiKDecayTreeProducer() override {};

private:
    virtual void produce(edm::Event&, const edm::EventSetup&);

    // ----------member data ---------------------------
    edm::EDGetTokenT<vector<pat::PackedCandidate>> PFCandSrc_;
    edm::EDGetTokenT<vector<reco::Vertex>> vtxSrc_;
    edm::EDGetTokenT<vector<pat::Muon>> TrgMuonSrc_;
    edm::EDGetTokenT<vector<pat::Muon>> muonSrc_;
    edm::EDGetTokenT<reco::BeamSpot> beamSpotSrc_;

    double mass_Mu   = 0.10565;
    double mass_Pi   = 0.13957;
    double mass_K    = 0.49367;
    double mass_Jpsi = 3.09691;
    double mass_B0   = 5.27963;

    int verbose = 0;
};



B2JpsiKDecayTreeProducer::B2JpsiKDecayTreeProducer(const edm::ParameterSet &iConfig) :
  beamSpotSrc_( consumes<reco::BeamSpot> ( edm::InputTag("offlineBeamSpot") ) )
{
    verbose = iConfig.getParameter<int>( "verbose" );
    PFCandSrc_ = consumes<vector<pat::PackedCandidate>>(edm::InputTag("packedPFCandidates"));
    muonSrc_ = consumes<vector<pat::Muon>>(edm::InputTag("slimmedMuons"));
    vtxSrc_ = consumes<vector<reco::Vertex>>(edm::InputTag("offlineSlimmedPrimaryVertices"));
    TrgMuonSrc_ = consumes<vector<pat::Muon>> ( iConfig.getParameter<edm::InputTag>( "trgMuons" ) );

    produces<map<string, float>>("outputNtuplizer");
    produces<map<string, vector<float>>>("outputVecNtuplizer");
}


void B2JpsiKDecayTreeProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
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

    edm::Handle<reco::BeamSpot> beamSpotHandle;
    iEvent.getByToken(beamSpotSrc_, beamSpotHandle);

    vector<RefCountedKinematicTree> vecJpsiKinTree;
    vector<double> vecMass_mumu;
    vector<vtxu::kinFitResuts> vecRes_mumu;
    vector<pair<pat::Muon, pat::Muon>> vecJpsiMuons;
    vector<vector<reco::Vertex>> muonsSoftVertexes;
    for(uint i1 = 0; i1 < nMu -1; i1++){
      auto m1 = (*muonHandle)[i1];
      if (m1.innerTrack().isNull()) continue;

      vector<reco::Vertex> vtxMu1;
      for (auto vtx : (*vtxHandle)) {
        if (m1.isSoftMuon(vtx)) vtxMu1.push_back(vtx);
      }
      if (vtxMu1.size() == 0) continue;

      for(uint i2 = i1+1; i2 < nMu; i2++){
        auto m2 = (*muonHandle)[i2];
        if (m2.innerTrack().isNull()) continue;
        if(m1.charge() * m2.charge() != -1) continue;

        vector<reco::Vertex> possibleVtxs;
        for (auto vtx : vtxMu1) {
          if (m2.isSoftMuon(vtx)) possibleVtxs.push_back(vtx);
        }
        if (possibleVtxs.size() == 0) continue;

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

        vecJpsiKinTree.push_back(kinTree);
        vecMass_mumu.push_back(massMuPair);
        vecRes_mumu.push_back(res);
        vecJpsiMuons.push_back(make_pair(mup, mum));
        muonsSoftVertexes.push_back(possibleVtxs);
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

    int n_K = 0, n_B = 0;

    /*
    ############################################################################
                              Look for the K
    ############################################################################
    */
    for(uint i_K = 0; i_K < N_pfCand; ++i_K) {
      const pat::PackedCandidate & K = (*pfCandHandle)[i_K];
      if (!K.hasTrackDetails()) continue;

      if (abs(K.pdgId()) != 211) continue;
      if (K.isTrackerMuon() || K.isStandAloneMuon()) continue;

      //Require a minimum pt
      if(K.pt() < __pThad_min__) continue;
      // Require to be close to the trigger muon;
      auto K_tk = K.bestTrack();
      // if (fabs(K_tk->dz(primaryVtx.position()) - trgMu.muonBestTrack()->dz(primaryVtx.position())) > __dzMax__) continue;
      if (fabs( K_tk->dz(primaryVtx.position()) ) > __dzMax__) continue;
      // Require significant impact parameter
      auto dxy = K_tk->dxy(primaryVtx.position());
      double dxyErr = vtxu::dxyError(*K_tk,*beamSpotHandle);
      auto K_sigdxy_PV = fabs(dxy)/dxyErr;
      auto K_norm_chi2 = K_tk->normalizedChi2();
      auto K_N_valid_hits = K_tk->numberOfValidHits();
      if (K_sigdxy_PV < __sigIPpfCand_min__) continue;

      n_K++;
      if(verbose){cout << Form("K cand found i:%d, pt:%.2f, eta:%.2f, phi:%.2f ", i_K, K.pt(), K.eta(), K.phi()) << endl;}

      /*
      ############################################################################
                Make a B0 -> J/psi K using the candidates found above
      ############################################################################
      */
      for(uint i_J = 0; i_J < vecJpsiKinTree.size(); ++i_J) {
        auto JpsiKinTree = vecJpsiKinTree[i_J];
        JpsiKinTree->movePointerToTheTop();
        auto Jpsi = JpsiKinTree->currentParticle();
        auto vtxJpsi = JpsiKinTree->currentDecayVertex();

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
            auto cosT_B_PV = vtxu::computePointingCosTransverse(primaryVtx, vtxB, B);
            (*outv)["cosT_B_PV_"+tag].push_back(cosT_B_PV);
            auto d_vtxB_PV = vtxu::vtxsDistance(primaryVtx, vtxB);
            double sigd_vtxB_PV = fabs(d_vtxB_PV.first)/d_vtxB_PV.second;
            (*outv)["sigd_vtxB_PV_"+tag].push_back(sigd_vtxB_PV);
            auto dxy_vtxB_PV = vtxu::vtxsTransverseDistance(primaryVtx, vtxB);
            (*outv)["sigdxy_vtxB_PV_"+tag].push_back(fabs(dxy_vtxB_PV.first)/dxy_vtxB_PV.second);

            // B2JpsiKDecayTreeProducer::AddTLVToOut(vtxu::getTLVfromKinPart(B), string("B_"+tag), outv);
            auto pvec = B->currentState().globalMomentum();
            TLorentzVector out;
            out.SetXYZM(pvec.x(), pvec.y(), pvec.z(), m);
            (*outv)["B_"+tag+"_pt"].push_back(out.Pt());
            (*outv)["B_"+tag+"_eta"].push_back(out.Eta());
            (*outv)["B_"+tag+"_phi"].push_back(out.Phi());
          }
          else{
            (*outv)["cos_B_PV_"+tag].push_back(0);
            (*outv)["cosT_B_PV_"+tag].push_back(0);
            (*outv)["sigd_vtxB_PV_"+tag].push_back(0);
            (*outv)["sigdxy_vtxB_PV_"+tag].push_back(0);
            (*outv)["B_"+tag+"_pt"].push_back(0);
            (*outv)["B_"+tag+"_eta"].push_back(0);
            (*outv)["B_"+tag+"_phi"].push_back(0);
          }
        };


        RefCountedKinematicTree BKinTree;
        BKinTree = vtxu::FitB_mumuK(iSetup, vecJpsiMuons[i_J].first, vecJpsiMuons[i_J].second, K);
        auto res = vtxu::fitQuality(BKinTree, __PvalChi2Vtx_min__);
        if(!res.isValid) continue;
        auto mass = BKinTree->currentParticle()->currentState().mass();
        if(fabs(mass - mass_B0) > __dmB0_max__) continue;

        BKinTree->movePointerToTheTop();
        auto B = BKinTree->currentParticle();
        auto vtxB = BKinTree->currentDecayVertex();
        // Looking for the best vertex to associate it with
        uint i_best = 0;
        auto maxCos = vtxu::computePointingCos(muonsSoftVertexes[i_J][0], vtxB, B);
        for(uint i_vtx = 1; i_vtx < muonsSoftVertexes[i_J].size(); i_vtx++) {
          auto auxCos = vtxu::computePointingCos(muonsSoftVertexes[i_J][i_vtx], vtxB, B);
          if(auxCos > maxCos) {
            maxCos = auxCos;
            i_best = i_vtx;
          }
        }
        auto bestVtx = muonsSoftVertexes[i_J][i_best];

        addToOutBFit("mumuK", res, mass, BKinTree, bestVtx, &(*outputVecNtuplizer));

        BKinTree->movePointerToTheFirstChild();
        auto refit_Mu1 = BKinTree->currentParticle();
        auto p4_refit_Mu1 = vtxu::getTLVfromKinPart(refit_Mu1);
        BKinTree->movePointerToTheNextChild();
        auto refit_Mu2 = BKinTree->currentParticle();
        auto p4_refit_Mu2 = vtxu::getTLVfromKinPart(refit_Mu2);
        BKinTree->movePointerToTheNextChild();
        auto refit_K = BKinTree->currentParticle();
        auto p4_refit_K = vtxu::getTLVfromKinPart(refit_K);

        auto cos_Jpsi_PV = vtxu::computePointingCos(bestVtx, vtxJpsi, Jpsi);
        auto cosT_Jpsi_PV = vtxu::computePointingCosTransverse(bestVtx, vtxJpsi, Jpsi);
        auto d_vtxJpsi_PV = vtxu::vtxsDistance(bestVtx, vtxJpsi);
        auto sigd_vtxJpsi_PV = fabs(d_vtxJpsi_PV.first/d_vtxJpsi_PV.second);
        auto dxy_vtxJpsi_PV = vtxu::vtxsTransverseDistance(bestVtx, vtxJpsi);

        n_B++;

        /*
        ############################################################################
                              Compute analysis variables
        ############################################################################
        */
        AddTLVToOut(vtxu::getTLVfromCand(K, mass_K), string("K"), &(*outputVecNtuplizer));
        (*outputVecNtuplizer)["K_charge"].push_back(K.charge());
        (*outputVecNtuplizer)["K_sigdxy_PV"].push_back(K_sigdxy_PV);
        (*outputVecNtuplizer)["K_norm_chi2"].push_back(K_norm_chi2);
        (*outputVecNtuplizer)["K_N_valid_hits"].push_back(K_N_valid_hits);

        auto mup = vecJpsiMuons[i_J].first;
        AddTLVToOut(vtxu::getTLVfromTrack(*(mup.bestTrack()), mass_Mu), string("mup"), &(*outputVecNtuplizer));
        auto dxy_mup = mup.innerTrack()->dxy(primaryVtx.position());
        (*outputVecNtuplizer)["mup_dxy"].push_back(dxy_mup);
        (*outputVecNtuplizer)["mup_sigdxy_PV"].push_back(fabs(dxy_mup)/vtxu::dxyError(*mup.innerTrack(),primaryVtx));
        (*outputVecNtuplizer)["mup_isTrg"].push_back(trgMuIdx(mup, (*trgMuonsHandle) ));

        auto mum = vecJpsiMuons[i_J].second;
        AddTLVToOut(vtxu::getTLVfromTrack(*(mum.bestTrack()), mass_Mu), string("mum"), &(*outputVecNtuplizer));
        auto dxy_mum = mum.innerTrack()->dxy(primaryVtx.position());
        (*outputVecNtuplizer)["mum_dxy"].push_back(dxy_mum);
        (*outputVecNtuplizer)["mum_sigdxy_PV"].push_back(fabs(dxy_mum)/vtxu::dxyError(*mum.innerTrack(),primaryVtx));
        (*outputVecNtuplizer)["mum_isTrg"].push_back(trgMuIdx(mum, (*trgMuonsHandle) ));

        (*outputVecNtuplizer)["chi2_mumu"].push_back(vecRes_mumu[i_J].chi2);
        (*outputVecNtuplizer)["dof_mumu"].push_back(vecRes_mumu[i_J].dof);
        (*outputVecNtuplizer)["pval_mumu"].push_back(vecRes_mumu[i_J].pval);
        (*outputVecNtuplizer)["mass_mumu"].push_back(vecMass_mumu[i_J]);
        (*outputVecNtuplizer)["cos_Jpsi_PV"].push_back(cos_Jpsi_PV);
        (*outputVecNtuplizer)["cosT_Jpsi_PV"].push_back(cosT_Jpsi_PV);
        (*outputVecNtuplizer)["d_vtxJpsi_PV"].push_back(d_vtxJpsi_PV.first);
        (*outputVecNtuplizer)["sigd_vtxJpsi_PV"].push_back(sigd_vtxJpsi_PV);
        (*outputVecNtuplizer)["sigdxy_vtxJpsi_PV"].push_back(dxy_vtxJpsi_PV.first/dxy_vtxJpsi_PV.second);

        AddTLVToOut(p4_refit_Mu1, string("mupRefit"), &(*outputVecNtuplizer));
        AddTLVToOut(p4_refit_Mu2, string("mumRefit"), &(*outputVecNtuplizer));
        AddTLVToOut(p4_refit_Mu1+p4_refit_Mu1, string("JpsiRefit"), &(*outputVecNtuplizer));
        AddTLVToOut(p4_refit_K, string("KRefit"), &(*outputVecNtuplizer));

        /*
        ############################################################################
                Veto the presence of additional tracks in the D* mu vertex
        ############################################################################
        */
        uint N_compatible_tk = 0;
        vector<double> tksAdd_massVis = {};
        vector<double> tksAdd_pval = {};
        vector<double> tksAdd_charge = {};
        vector<double> tksAdd_pt = {};
        vector<double> tksAdd_eta = {};
        vector<double> tksAdd_phi = {};
        vector<double> tksAdd_dz = {};
        vector<double> tksAdd_sigdca_vtxB = {};
        vector<double> tksAdd_cos_PV = {};
        if(verbose) {cout << "Looking for additional tracks" << endl;}
        for(uint i_tk = 0; i_tk < N_pfCand; ++i_tk) {
          // Pion different from the pion from D0
          if( i_tk==i_K ) continue;

          const pat::PackedCandidate & ptk = (*pfCandHandle)[i_tk];
          if (!ptk.hasTrackDetails()) continue;
          //Require a charged hadron
          if (abs(ptk.pdgId()) != 211 ) continue;
          if (ptk.pt() < 0.3) continue;
          // Require to be close to the trigger muon;
          auto tk = ptk.bestTrack();
          if (fabs(tk->dz(primaryVtx.position()) - K.bestTrack()->dz(primaryVtx.position())) > __dzMax__) continue;

          GlobalPoint auxp(vtxB->position().x(), vtxB->position().y(), vtxB->position().z());
          auto dca = vtxu::computeDCA(iSetup, ptk, auxp);

          // Check if it makes a consistent vertex
          auto kinTree = vtxu::FitB_mumupiK(iSetup, vecJpsiMuons[i_J].first, vecJpsiMuons[i_J].second, ptk, K, false, false, false);
          auto res = vtxu::fitQuality(kinTree, __PvalChi2Vtx_min__);
          if(!res.isGood) continue;
          if (verbose) {cout << Form("Trk pars: pt=%.2f eta=%.2f phi=%.2f\n", ptk.pt(), ptk.eta(), ptk.phi());}

          kinTree->movePointerToTheTop();
          auto p_vis = kinTree->currentParticle();
          auto m_vis = p_vis->currentState().mass();
          if (m_vis > mass_B0 + __dmB0_max__) continue;

          kinTree->movePointerToTheFirstChild();
          kinTree->movePointerToTheNextChild();
          kinTree->movePointerToTheNextChild();
          auto refit_pi = kinTree->currentParticle();
          auto refit_pi_p4 = vtxu::getTLVfromKinPart(refit_pi);

          tksAdd_massVis.push_back(m_vis);
          tksAdd_pval.push_back(res.pval);
          tksAdd_charge.push_back(ptk.pdgId()/211);
          tksAdd_pt.push_back(refit_pi_p4.Pt());
          tksAdd_eta.push_back(refit_pi_p4.Eta());
          tksAdd_phi.push_back(refit_pi_p4.Phi());
          tksAdd_dz.push_back(tk->dz(bestVtx.position()));
          tksAdd_sigdca_vtxB.push_back(fabs(dca.first)/dca.second);
          tksAdd_cos_PV.push_back(vtxu::computePointingCos(bestVtx, vtxB, refit_pi));
          N_compatible_tk++;
        }
        if ( (*outputVecNtuplizer)["nTksAdd"].size() == 0 ) {
          (*outputVecNtuplizer)["tksAdd_massVis"] = {};
          (*outputVecNtuplizer)["tksAdd_pval"] = {};
          (*outputVecNtuplizer)["tksAdd_charge"] = {};
          (*outputVecNtuplizer)["tksAdd_pt"] = {};
          (*outputVecNtuplizer)["tksAdd_eta"] = {};
          (*outputVecNtuplizer)["tksAdd_phi"] = {};
          (*outputVecNtuplizer)["tksAdd_dz"] = {};
          (*outputVecNtuplizer)["tksAdd_sigdca_vtxB"] = {};
          (*outputVecNtuplizer)["tksAdd_cos_PV"] = {};
        }
        (*outputVecNtuplizer)["nTksAdd"].push_back(N_compatible_tk);
        if(verbose) {cout << "Compatible tracks: " << N_compatible_tk << endl;}
        for(uint i = 0; i < N_compatible_tk; i++){
          (*outputVecNtuplizer)["tksAdd_massVis"].push_back(tksAdd_massVis[i]);
          (*outputVecNtuplizer)["tksAdd_pval"].push_back(tksAdd_pval[i]);
          (*outputVecNtuplizer)["tksAdd_charge"].push_back(tksAdd_charge[i]);
          (*outputVecNtuplizer)["tksAdd_pt"].push_back(tksAdd_pt[i]);
          (*outputVecNtuplizer)["tksAdd_eta"].push_back(tksAdd_eta[i]);
          (*outputVecNtuplizer)["tksAdd_phi"].push_back(tksAdd_phi[i]);
          (*outputVecNtuplizer)["tksAdd_dz"].push_back(tksAdd_dz[i]);
          (*outputVecNtuplizer)["tksAdd_sigdca_vtxB"].push_back(tksAdd_sigdca_vtxB[i]);
          (*outputVecNtuplizer)["tksAdd_cos_PV"].push_back(tksAdd_cos_PV[i]);
        }

        int auxC = 0;
        for(auto auxN : (*outputVecNtuplizer)["nTksAdd"]) auxC += auxN;
        if( auxC != int ((*outputVecNtuplizer)["tksAdd_massVis"].size()) ){
          cout << "Number of tracks and tracks details lenght not matching" << endl;
          assert(false);
        }

        if(n_B >= 100) break;
      }

      if(n_B >= 100) break;
    }

    (*outputNtuplizer)["Run"] = iEvent.run();
    (*outputNtuplizer)["LumiBlock"] = iEvent.luminosityBlock();
    (*outputNtuplizer)["eventNumber"] = iEvent.id().event();

    (*outputNtuplizer)["n_K"] = n_K;
    (*outputNtuplizer)["n_Jpsi"] = vecJpsiKinTree.size();
    (*outputNtuplizer)["n_B"] = n_B;


    iEvent.put(move(outputNtuplizer), "outputNtuplizer");
    iEvent.put(move(outputVecNtuplizer), "outputVecNtuplizer");
    return;
}


void B2JpsiKDecayTreeProducer::AddTLVToOut(TLorentzVector v, string n, map<string, vector<float>>* outv) {
  (*outv)[n+"_pt"].push_back(v.Pt());
  // (*outv)[n+"_pz"].push_back(v.Pz());
  (*outv)[n+"_eta"].push_back(v.Eta());
  (*outv)[n+"_phi"].push_back(v.Phi());
  // (*outv)[n+"_P"].push_back(v.P());
  // (*outv)[n+"_E"].push_back(v.E());
  return;
}

int B2JpsiKDecayTreeProducer::trgMuIdx(pat::Muon m, vector<pat::Muon> trgMuons) {
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

bool B2JpsiKDecayTreeProducer::isMuonFromJpsiID(pat::Muon m, reco::Vertex pVtx, pat::Muon trgMu) {
  if(m.innerTrack().isNull()) return false;
  if (fabs(m.innerTrack()->dz(pVtx.position()) - trgMu.innerTrack()->dz(pVtx.position())) > __dzMax__) return false;

  if(m.innerTrack()->hitPattern().pixelLayersWithMeasurement() < 2) return false;
  if(!m.innerTrack()->quality(reco::TrackBase::highPurity)) return false;
  if(!m.isGood("TMOneStationTight")) return false;
  if(m.innerTrack()->normalizedChi2() > 1.8) return false;

  double dxy = m.innerTrack()->dxy(pVtx.position());
  float sigdxy = fabs(dxy)/vtxu::dxyError(*m.innerTrack(),pVtx);
  if (sigdxy < 2) return false;

  return true;
}

DEFINE_FWK_MODULE(B2JpsiKDecayTreeProducer);
