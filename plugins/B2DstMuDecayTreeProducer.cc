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
#define __dRMax__ 2.0
#define __cos_kpi_vtxMu_min__ 0.97 // Optimized in D0 fitting
#define __d_vtxkpi_vtxMu_min__ 0.02 // Optimized in D0 fitting
#define __dmD0_max__ 0.039 // 3*0.013 GeV(i.e. 3 sigma) form the D0 mass
#define __dmD0pis_max__ 0.0024 // 3*0.8 MeV (i.e 3 inflated sigma) from Dst mass
#define __mDstMu_max__ 7.0 // Some reasonable cut on the mass
#define __cos_DstMu_vtxMu_min__ 0.8 // Some loose cut tuned on MC
// #define __PvalChi2FakeVtx_min__ 0.90 // Very loose cut


using namespace std;

class B2DstMuDecayTreeProducer : public edm::EDProducer {

public:

    explicit B2DstMuDecayTreeProducer(const edm::ParameterSet &iConfig);
    void DeclareTLVToOut(string, map<string, vector<float>>*);
    void AddTLVToOut(TLorentzVector, string, map<string, vector<float>>*);

    ~B2DstMuDecayTreeProducer() override {};

private:
    virtual void produce(edm::Event&, const edm::EventSetup&);

    // ----------member data ---------------------------
    edm::EDGetTokenT<vector<pat::PackedCandidate>> PFCandSrc_;
    edm::EDGetTokenT<vector<reco::Vertex>> vtxSrc_;
    edm::EDGetTokenT<vector<pat::Muon>> TrgMuonSrc_;

    double mass_Mu  = 0.10565;
    double mass_Pi  = 0.13957;
    double mass_K   = 0.49367;
    double mass_D0  = 1.86483;
    double mass_Dst = 2.01026;
    double mass_B0  = 5.27961;

    int verbose = 0;
};



B2DstMuDecayTreeProducer::B2DstMuDecayTreeProducer(const edm::ParameterSet &iConfig)
{
    verbose = iConfig.getParameter<int>( "verbose" );
    PFCandSrc_ = consumes<vector<pat::PackedCandidate>>(edm::InputTag("packedPFCandidates"));
    vtxSrc_ = consumes<vector<reco::Vertex>>(edm::InputTag("offlineSlimmedPrimaryVertices"));
    TrgMuonSrc_ = consumes<vector<pat::Muon>> ( iConfig.getParameter<edm::InputTag>( "trgMuons" ) );

    produces<map<string, float>>("outputNtuplizer");
    produces<map<string, vector<float>>>("outputVecNtuplizer");
}


void B2DstMuDecayTreeProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    edm::Handle<vector<pat::PackedCandidate>> pfCandHandle;
    iEvent.getByToken(PFCandSrc_, pfCandHandle);
    unsigned int N_pfCand = pfCandHandle->size();

    edm::Handle<vector<reco::Vertex>> vtxHandle;
    iEvent.getByToken(vtxSrc_, vtxHandle);


    edm::Handle<vector<pat::Muon>> trgMuonsHandle;
    iEvent.getByToken(TrgMuonSrc_, trgMuonsHandle);

    // Output collection
    unique_ptr<map<string, float>> outputNtuplizer(new map<string, float>);
    unique_ptr<map<string, vector<float>>> outputVecNtuplizer(new map<string, vector<float>>);

    auto mu = (*trgMuonsHandle)[0];
    pat::PackedCandidate trgMu;
    uint i_trgMu = 0;
    for (i_trgMu = 0; i_trgMu < N_pfCand; ++i_trgMu) {
      const pat::PackedCandidate & p = (*pfCandHandle)[i_trgMu];
      if (fabs(mu.pt() - p.pt())/mu.pt() > 5e-3) continue;
      if (fabs(mu.eta() - p.eta()) > 5e-3) continue;
      if (fabs(vtxu::dPhi(mu.phi(), p.phi())) > 5e-3) continue;
      trgMu = p;
      break;
    }

    auto vtxMu = trgMu.vertexRef();
    GlobalPoint p_vtxMu(vtxMu->x(), vtxMu->y(), vtxMu->z());


    int n_K = 0, n_pi = 0, n_D0 = 0, n_pis = 0, n_Dst = 0, n_B = 0;

    (*outputVecNtuplizer)["chi2_kpi"] = {};
    (*outputVecNtuplizer)["mass_kpi"] = {};
    (*outputVecNtuplizer)["dca_kpi_vtxMu"] = {};
    (*outputVecNtuplizer)["sigdca_kpi_vtxMu"] = {};
    (*outputVecNtuplizer)["cos_kpi_vtxMu"] = {};
    (*outputVecNtuplizer)["d_vtxkpi_vtxMu"] = {};
    (*outputVecNtuplizer)["sigd_vtxkpi_vtxMu"] = {};

    (*outputVecNtuplizer)["chi2_D0pis"] = {};
    (*outputVecNtuplizer)["mass_D0pis"] = {};
    (*outputVecNtuplizer)["dca_D0pis_vtxMu"] = {};
    (*outputVecNtuplizer)["sigdca_D0pis_vtxMu"] = {};
    (*outputVecNtuplizer)["cos_D0pis_vtxMu"] = {};
    (*outputVecNtuplizer)["d_vtxD0pis_vtxMu"] = {};
    (*outputVecNtuplizer)["sigd_vtxD0pis_vtxMu"] = {};

    (*outputVecNtuplizer)["chi2_MuDst"] = {};
    (*outputVecNtuplizer)["mass_MuDst"] = {};
    (*outputVecNtuplizer)["cos_MuDst_vtxMu"] = {};

    DeclareTLVToOut("K", &(*outputVecNtuplizer));
    DeclareTLVToOut("pi", &(*outputVecNtuplizer));
    DeclareTLVToOut("pis", &(*outputVecNtuplizer));
    DeclareTLVToOut("D0", &(*outputVecNtuplizer));
    DeclareTLVToOut("Dst", &(*outputVecNtuplizer));
    DeclareTLVToOut("B", &(*outputVecNtuplizer));

    (*outputVecNtuplizer)["M2_miss"] = {};
    (*outputVecNtuplizer)["q2"] = {};
    (*outputVecNtuplizer)["Est_mu"] = {};

    if (verbose) {cout <<"-------------------- Evt -----------------------\n";}
    /*
    ############################################################################
                              Look for the K+ (321)
    ############################################################################
    */
    for(uint i_k = 0; i_k < N_pfCand; ++i_k) {
      if(i_k==i_trgMu) continue;

      const pat::PackedCandidate & k = (*pfCandHandle)[i_k];
      if (!k.hasTrackDetails()) continue;
      //Require a positive charged hadron
      if (k.pdgId() != 211 ) continue;
      // Require to be close to the trigger muon;
      if (fabs(k.dz() - trgMu.dz()) > __dzMax__) continue;
      if (vtxu::dR(k.phi(), trgMu.phi(), k.eta(), trgMu.eta()) > __dRMax__) continue;
      // Require significant displacement from PV
      // if (verbose) {cout << "###-### K number " << n_K << endl;}

      n_K++;

      /*
      ############################################################################
                               Look for the pi- (-211)
      ############################################################################
      */
      for(uint i_pi = 0; i_pi < N_pfCand; ++i_pi) {
        if(i_pi==i_trgMu || i_pi==i_k) continue;

        const pat::PackedCandidate & pi = (*pfCandHandle)[i_pi];
        if (!pi.hasTrackDetails()) continue;
        //Require a negative charged hadron
        if (pi.pdgId() != -211 ) continue;
        // Require to be close to the trigger muon;
        if (fabs(pi.dz() - trgMu.dz()) > __dzMax__) continue;
        if (vtxu::dR(pi.phi(), trgMu.phi(), pi.eta(), trgMu.eta()) > __dRMax__) continue;

        // Require to be close to the K;
        if (fabs(pi.dz() - k.dz()) > __dzMax__) continue;
        if (vtxu::dR(pi.phi(), k.phi(), pi.eta(), k.eta()) > __dRMax__) continue;

        //Fit the vertex
        auto D0KinTree = vtxu::FitD0(iSetup, pi, k, false, 0);
        bool accept_pi = false;
        double chi2_kpi;
        if(D0KinTree->isValid()) {
          D0KinTree->movePointerToTheTop();
          auto D0vtx = D0KinTree->currentDecayVertex();
          chi2_kpi = D0vtx->chiSquared();
          auto max_chi2 = TMath::ChisquareQuantile(__PvalChi2Vtx_max__, D0vtx->degreesOfFreedom());
          if (chi2_kpi > 0 && chi2_kpi < max_chi2) accept_pi = true;
        }
        if(!accept_pi) continue;
        n_pi++;

        auto D0 = D0KinTree->currentParticle();
        auto D0vtx = D0KinTree->currentDecayVertex();

        auto mass_kpi = D0->currentState().mass();

        auto dca = vtxu::computeDCA(D0->refittedTransientTrack(), p_vtxMu);
        auto dca_kpi_vtxMu = dca.first;
        auto sigdca_kpi_vtxMu = fabs(dca.first)/dca.second;

        TVector3 dvtx(D0vtx->position().x() - vtxMu->x(),
                      D0vtx->position().y() - vtxMu->y(),
                      D0vtx->position().z() - vtxMu->z()
                     );
        TVector3 pD0(D0->currentState().globalMomentum().x(),
                     D0->currentState().globalMomentum().y(),
                     D0->currentState().globalMomentum().z()
                     );
        double dalpha = dvtx.Angle(pD0);
        auto cos_kpi_vtxMu = cos(dalpha);

        auto dvtx_kpiPV = vtxu::vtxsDistance(vtxMu, D0vtx);
        auto d_vtxkpi_vtxMu = dvtx_kpiPV.first;
        auto sigd_vtxkpi_vtxMu = d_vtxkpi_vtxMu/dvtx_kpiPV.second;

        // Make selection cut to get good D0
        bool accept_D0 = cos_kpi_vtxMu > __cos_kpi_vtxMu_min__;
        accept_D0 &= d_vtxkpi_vtxMu > __d_vtxkpi_vtxMu_min__;
        accept_D0 &= fabs(mass_kpi - mass_D0) < __dmD0_max__;
        if(!accept_D0) continue;
        n_D0++;

        /*
        ############################################################################
                         Look for the soft pion to make a Dst-
        ############################################################################
        */
        for(uint i_pis = 0; i_pis < N_pfCand; ++i_pis) {
          // Pion different from the pion from D0
          if(i_pis==i_trgMu || i_pis==i_k || i_pis==i_pi) continue;

          const pat::PackedCandidate & pis = (*pfCandHandle)[i_pis];
          if (!pis.hasTrackDetails()) continue;
          //Require a negative charged hadron
          if (pis.pdgId() != -211 ) continue;

          // Require to be close to the trigger muon;
          if (fabs(pis.dz() - trgMu.dz()) > __dzMax__) continue;
          if (vtxu::dR(pis.phi(), trgMu.phi(), pis.eta(), trgMu.eta()) > __dRMax__) continue;

          // Fit the Dst vertex
          //  Refitting the D0 with mass constrint
          auto DstKinTree = vtxu::FitDst_fitD0wMassConstraint(iSetup, pis, pi, k, false, 0);
          //  Without refitting the D0 with mass constrint
          // auto DstKinTree = vtxu::FitDst(iSetup, pis, D0, false, 0);
          bool accept_pis = false;
          double chi2_D0pis;
          if(DstKinTree->isValid()) {
            DstKinTree->movePointerToTheTop();
            auto Dstvtx = DstKinTree->currentDecayVertex();
            chi2_D0pis = Dstvtx->chiSquared();
            auto max_chi2 = TMath::ChisquareQuantile(__PvalChi2Vtx_max__, Dstvtx->degreesOfFreedom());
            if (chi2_D0pis > 0 && chi2_D0pis < max_chi2) accept_pis = true;
          }
          if(!accept_pis) continue;
          n_pis++;

          auto Dst = DstKinTree->currentParticle();
          auto Dstvtx = DstKinTree->currentDecayVertex();

          auto mass_D0pis = Dst->currentState().mass();

          auto dca = vtxu::computeDCA(Dst->refittedTransientTrack(), p_vtxMu);
          auto dca_D0pis_vtxMu = dca.first;
          auto sigdca_D0pis_vtxMu = fabs(dca.first)/dca.second;

          TVector3 dvtx(Dstvtx->position().x() - vtxMu->x(),
                        Dstvtx->position().y() - vtxMu->y(),
                        Dstvtx->position().z() - vtxMu->z()
                       );
          TVector3 pDst(Dst->currentState().globalMomentum().x(),
                       Dst->currentState().globalMomentum().y(),
                       Dst->currentState().globalMomentum().z()
                       );
          double dalpha = dvtx.Angle(pDst);
          auto cos_D0pis_vtxMu = cos(dalpha);

          auto dvtx_D0pisPV = vtxu::vtxsDistance(vtxMu, Dstvtx);
          auto d_vtxD0pis_vtxMu = dvtx_D0pisPV.first;
          auto sigd_vtxD0pis_vtxMu = d_vtxD0pis_vtxMu/dvtx_D0pisPV.second;

          // Make selection cut to get good Dst
          bool accept_Dst = fabs(mass_D0pis - mass_Dst) < __dmD0pis_max__;
          if (!accept_Dst) continue;
          n_Dst++;
          if (verbose) {cout << "D* found\n";}

          //  Look for the Dst- and trigger muon to make a vertex
          auto BKinTree = vtxu::FitVtxMuDst(iSetup, Dst, trgMu, 0);
          bool accept_MuDst = false;
          double chi2_MuDst;
          if(BKinTree->isValid()) {
            BKinTree->movePointerToTheTop();
            auto Bvtx = BKinTree->currentDecayVertex();
            chi2_MuDst = Bvtx->chiSquared();
            auto max_chi2 = TMath::ChisquareQuantile(__PvalChi2Vtx_max__, Bvtx->degreesOfFreedom());
            if (chi2_MuDst > 0 && chi2_MuDst < max_chi2) accept_MuDst = true;
          }
          if(!accept_MuDst) continue;

          auto MuDst = BKinTree->currentParticle();
          auto Bvtx = BKinTree->currentDecayVertex();

          BKinTree->movePointerToTheFirstChild();
          auto refit_Mu = BKinTree->currentParticle();
          BKinTree->movePointerToTheNextChild();
          auto refit_Dst = BKinTree->currentParticle();

          auto mass_MuDst = MuDst->currentState().mass();

          TVector3 dvtxB(Bvtx->position().x() - vtxMu->x(),
                             Bvtx->position().y() - vtxMu->y(),
                             Bvtx->position().z() - vtxMu->z()
                            );
          TVector3 pMuDst(MuDst->currentState().globalMomentum().x(),
                          MuDst->currentState().globalMomentum().y(),
                          MuDst->currentState().globalMomentum().z()
                         );
          double dalphaB = dvtxB.Angle(pMuDst);
          auto cos_MuDst_vtxMu = cos(dalphaB);



          bool accept_B = mass_MuDst < __mDstMu_max__;
          accept_B &= cos_MuDst_vtxMu > __cos_DstMu_vtxMu_min__;
          if(!accept_B) continue;
          if (verbose) {cout << "B->D* mu found\n";}

          /*
          ############################################################################
                  Veto the presence of additional tracks in the D* mu vertex
          ############################################################################
          */
          // int N_compatible_tk = 0;
          // for(uint i_tk = 0; i_tk < N_pfCand; ++i_tk) {
          //   // Pion different from the pion from D0
          //   if( i_tk==i_trgMu || i_tk==i_k || i_tk==i_pi || i_tk==i_pis) continue;
          //
          //   const pat::PackedCandidate & ptk = (*pfCandHandle)[i_tk];
          //   //Require a charged candidate
          //   if ( ptk.charge() == 0 ) continue;
          //   if (!ptk.hasTrackDetails()) continue;
          //
          //   // Require to be close to the trigger muon;
          //   if (fabs(ptk.dz() - trgMu.dz()) > __dzMax__) continue;
          //   if (vtxu::dR(ptk.phi(), trgMu.phi(), ptk.eta(), trgMu.eta()) > __dRMax__) continue;
          //
          //   //  Look for the Dst- and trigger muon to make a vertex
          //   auto VtxKinTree = vtxu::FitVtxDstPi(iSetup, Dst, ptk, 0);
          //   // auto VtxKinTree = vtxu::FitVtxMuDstPi(iSetup, Dst, trgMu, ptk, 0);
          //   bool accept_Vtx = false;
          //   double chi2_vtx;
          //   if(VtxKinTree->isValid()) {
          //     VtxKinTree->movePointerToTheTop();
          //     auto vtx = VtxKinTree->currentDecayVertex();
          //     chi2_vtx = vtx->chiSquared();
          //     auto max_chi2 = TMath::ChisquareQuantile(__PvalChi2FakeVtx_min__, vtx->degreesOfFreedom());
          //     if (chi2_vtx > 0 && chi2_vtx < max_chi2) accept_Vtx = true;
          //   }
          //   if(!accept_Vtx) continue;
          //
          //   auto particle = VtxKinTree->currentParticle();
          //   auto mass_particle = particle->currentState().mass();
          //   if(verbose) {cout << "Vtx Mass: " << mass_particle << endl;}
          //
          //   // accept_Vtx &=  mass_particle >
          //   N_compatible_tk++;
          // }
          // if(verbose) {cout << "Compatible tracks: " << N_compatible_tk << endl;}



          if(!accept_B) continue;
          n_B++;


          /*
          ############################################################################
                                Compute analysis variables
          ############################################################################
          */
          auto p4_vis = vtxu::getTLVfromKinPart(MuDst);
          double pz_B_reco = p4_vis.Pz() * mass_B0/ p4_vis.M();
          auto B_vect = dvtxB * ( pz_B_reco / dvtxB.z() );
          TLorentzVector p4_B;
          p4_B.SetVectM(B_vect, mass_B0);
          auto M2_miss = (p4_B - p4_vis).M2();
          // DEGUB --> Use the refitted tracks instead of the orginal ones!
          auto p4_Dst = vtxu::getTLVfromKinPart(refit_Dst);
          auto p4_mu = vtxu::getTLVfromKinPart(refit_Mu);
          auto q2 = (p4_B - p4_Dst).M2();

          TLorentzVector p4st_mu(p4_mu);
          p4st_mu.Boost(-1*p4_B.BoostVector());
          auto Est_mu = p4st_mu.E();

          (*outputVecNtuplizer)["chi2_kpi"].push_back(chi2_kpi);
          (*outputVecNtuplizer)["mass_kpi"].push_back(mass_kpi);
          (*outputVecNtuplizer)["dca_kpi_vtxMu"].push_back(dca_kpi_vtxMu);
          (*outputVecNtuplizer)["sigdca_kpi_vtxMu"].push_back(sigdca_kpi_vtxMu);
          (*outputVecNtuplizer)["cos_kpi_vtxMu"].push_back(cos_kpi_vtxMu);
          (*outputVecNtuplizer)["d_vtxkpi_vtxMu"].push_back(d_vtxkpi_vtxMu);
          (*outputVecNtuplizer)["sigd_vtxkpi_vtxMu"].push_back(sigd_vtxkpi_vtxMu);

          (*outputVecNtuplizer)["chi2_D0pis"].push_back(chi2_D0pis);
          (*outputVecNtuplizer)["mass_D0pis"].push_back(mass_D0pis);
          (*outputVecNtuplizer)["dca_D0pis_vtxMu"].push_back(dca_D0pis_vtxMu);
          (*outputVecNtuplizer)["sigdca_D0pis_vtxMu"].push_back(sigdca_D0pis_vtxMu);
          (*outputVecNtuplizer)["cos_D0pis_vtxMu"].push_back(cos_D0pis_vtxMu);
          (*outputVecNtuplizer)["d_vtxD0pis_vtxMu"].push_back(d_vtxD0pis_vtxMu);
          (*outputVecNtuplizer)["sigd_vtxD0pis_vtxMu"].push_back(sigd_vtxD0pis_vtxMu);

          (*outputVecNtuplizer)["chi2_MuDst"].push_back(chi2_MuDst);
          (*outputVecNtuplizer)["mass_MuDst"].push_back(mass_MuDst);
          (*outputVecNtuplizer)["cos_MuDst_vtxMu"].push_back(cos_MuDst_vtxMu);

          AddTLVToOut(vtxu::getTLVfromCand(k, mass_K), string("K"), &(*outputVecNtuplizer));
          AddTLVToOut(vtxu::getTLVfromCand(pi, mass_Pi), string("pi"), &(*outputVecNtuplizer));
          AddTLVToOut(vtxu::getTLVfromCand(pis, mass_Pi), string("pis"), &(*outputVecNtuplizer));
          AddTLVToOut(vtxu::getTLVfromKinPart(D0), string("D0"), &(*outputVecNtuplizer));
          AddTLVToOut(p4_Dst, string("Dst"), &(*outputVecNtuplizer));
          AddTLVToOut(p4_B, string("B"), &(*outputVecNtuplizer));

          (*outputVecNtuplizer)["M2_miss"].push_back(M2_miss);
          (*outputVecNtuplizer)["q2"].push_back(q2);
          (*outputVecNtuplizer)["Est_mu"].push_back(Est_mu);
        }


      }

    }

    (*outputNtuplizer)["n_K"] = n_K;
    (*outputNtuplizer)["n_pi"] = n_pi;
    (*outputNtuplizer)["n_D0"] = n_D0;
    (*outputNtuplizer)["n_pis"] = n_pis;
    (*outputNtuplizer)["n_Dst"] = n_Dst;
    (*outputNtuplizer)["n_B"] = n_B;


    iEvent.put(move(outputNtuplizer), "outputNtuplizer");
    iEvent.put(move(outputVecNtuplizer), "outputVecNtuplizer");
    return;
}


void B2DstMuDecayTreeProducer::DeclareTLVToOut(string n, map<string, vector<float>>* outv) {
  (*outv)[n+"_pt"] = {};
  (*outv)[n+"_eta"] = {};
  (*outv)[n+"_phi"] = {};
  (*outv)[n+"_P"] = {};
  (*outv)[n+"_E"] = {};
  return;
}


void B2DstMuDecayTreeProducer::AddTLVToOut(TLorentzVector v, string n, map<string, vector<float>>* outv) {
  (*outv)[n+"_pt"].push_back(v.Pt());
  (*outv)[n+"_eta"].push_back(v.Eta());
  (*outv)[n+"_phi"].push_back(v.Phi());
  (*outv)[n+"_P"].push_back(v.P());
  (*outv)[n+"_E"].push_back(v.E());
  return;
}

DEFINE_FWK_MODULE(B2DstMuDecayTreeProducer);
