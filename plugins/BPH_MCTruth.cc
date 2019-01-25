#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include <TLorentzVector.h>
#include <TVector.h>
#include <TMatrix.h>


//
// class declaration
//

using namespace std;

class BPH_MCTruth : public edm::EDProducer {

public:

    explicit BPH_MCTruth(const edm::ParameterSet &iConfig);

    ~BPH_MCTruth() override {};

    vector<int> RecoMC_matching(edm::Handle<edm::View<reco::Candidate>>, reco::Candidate *);
    vector<int> RecoMC_matching(edm::Handle<edm::View<reco::Candidate>>, pat::PackedGenParticle*);

private:

    virtual void produce(edm::Event&, const edm::EventSetup&);

    // ----------member data ---------------------------

    edm::EDGetTokenT<vector<reco::GenParticle>> PrunedParticlesSrc_;
    edm::EDGetTokenT<vector<pat::PackedGenParticle>> PackedParticlesSrc_;
    edm::EDGetTokenT<edm::View<reco::Candidate>> PFCandSrc_;
    edm::EDGetTokenT<edm::View<reco::Candidate>> muonSrc_;

    int verbose = 1;
};



BPH_MCTruth::BPH_MCTruth(const edm::ParameterSet &iConfig):
  PrunedParticlesSrc_( consumes<std::vector<reco::GenParticle>> ( iConfig.getParameter<edm::InputTag>( "PrunedGenParticlesCollection" ) ) ),
  PackedParticlesSrc_( consumes<std::vector<pat::PackedGenParticle>> ( iConfig.getParameter<edm::InputTag>( "PackedGenParticlesCollection" ) ) ),
  PFCandSrc_( consumes<edm::View<reco::Candidate>> ( iConfig.getParameter<edm::InputTag>( "PFCandCollection" ) ) ),
  muonSrc_( consumes<edm::View<reco::Candidate>> ( iConfig.getParameter<edm::InputTag>( "muonCollection" ) ) )
{
    produces<pat::CompositeCandidateCollection>();
}


void BPH_MCTruth::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    // Get prunedGenParticles
    edm::Handle<std::vector<reco::GenParticle>> PrunedGenParticlesHandle;
    iEvent.getByToken(PrunedParticlesSrc_, PrunedGenParticlesHandle);
    unsigned int PrunedGenParticlesNumber = PrunedGenParticlesHandle->size();

    // Get packedGenParticles
    edm::Handle<std::vector<pat::PackedGenParticle>> PackedGenParticlesHandle;
    iEvent.getByToken(PackedParticlesSrc_, PackedGenParticlesHandle);
    unsigned int PackedGenParticlesNumber = PackedGenParticlesHandle->size();

    edm::Handle<edm::View<reco::Candidate>> pfCandHandle;
    iEvent.getByToken(PFCandSrc_, pfCandHandle);
    // unsigned int pfCandNumber = pfCandHandle->size();

    edm::Handle<edm::View<reco::Candidate>> muonHandle;
    iEvent.getByToken(muonSrc_, muonHandle);
    // unsigned int muonNumber = muonHandle->size();

    // Output collection
    std::unique_ptr<pat::CompositeCandidateCollection> result( new pat::CompositeCandidateCollection );

    if (verbose) {cout << "EVT" << endl;}

    for(unsigned int i = 0; i < PrunedGenParticlesNumber; i++) {
      const reco::GenParticle & p = (*PrunedGenParticlesHandle)[i];
      if ( abs(p.pdgId()) > 500 && abs(p.pdgId()) % 500 < 50 && verbose) {
        cout << "Found: " << p.pdgId() << " (-->" << flush;
        for(unsigned int k = 0; k < p.numberOfDaughters(); k++) {
            reco::Candidate * d = (reco::Candidate *) p.daughter(k);
            cout << " " << d->pdgId() << flush;
          }
        cout << ')' << endl;
      }

      if ( abs(p.pdgId()) == 521 ) {
        pat::CompositeCandidate Bmeson;
        // Bmeson.addDaughter( muon1, "muon1");

        reco::Candidate * D0 = 0;
        reco::Candidate * mu = 0;
        reco::Candidate * mu_reco = 0;
        reco::Candidate * nu = 0;
        pat::PackedGenParticle * K = 0;
        reco::Candidate * K_reco = 0;
        pat::PackedGenParticle * pi = 0;
        reco::Candidate * pi_reco = 0;


        int BtoD0munu = 0;
        if (p.numberOfDaughters() == 3) {
          for(unsigned int k = 0; k < p.numberOfDaughters(); k++) {
            reco::Candidate * d = (reco::Candidate *) p.daughter(k);
            if (verbose) {
              cout << "\t" << k << " " << d->pdgId() << "  " << d->pt() << "  " << d->eta() << "  " << d->phi() << endl;
            }
            if ( abs(d->pdgId()) == 13) {
              mu = (reco::Candidate *) p.daughter(k);
              auto i_match = RecoMC_matching(muonHandle, mu);
              if (i_match.size() >= 1) {
                if (verbose) {
                  const reco::Candidate & pr = (*muonHandle)[i_match[0]];
                  cout << "\t\t" << i_match[0] << " " << pr.pdgId() << "  " << pr.pt() << "  " << pr.eta() << "  " << pr.phi() << endl;
                  cout << "\t\t vz: " << pr.vz() << endl;
                }
                mu_reco = (reco::Candidate *) & ((*muonHandle)[i_match[0]]);
              }
              else if (verbose){
                cout << "\tNo matching muon found" << endl;
              }
            }
            else if ( abs(d->pdgId()) == 14) nu = (reco::Candidate *) p.daughter(k);
            else if ( abs(d->pdgId()) == 421) {
              D0 = (reco::Candidate *) p.daughter(k);
              if (verbose) {
                cout << "\t (-->" << flush;
                for(unsigned int k = 0; k < PackedGenParticlesNumber; k++) {
                  pat::PackedGenParticle* p = (pat::PackedGenParticle*) &((*PackedGenParticlesHandle)[k]);
                  if(p->mother(0) == D0) {
                    cout << " " << p->pdgId() << flush;
                  }
                }
                cout << ")" << endl;
              }
            }
          }
        }
        if ( D0!=0 && mu!=0 && nu!=0) {
          if (verbose) {cout << "\tIt's a BtoD0munu" << endl;}
          BtoD0munu = 1;
        }

        int D0toKpi = 0;
        if(BtoD0munu) {
          if (D0->numberOfDaughters() == 2) {
            for(unsigned int k = 0; k < PackedGenParticlesNumber; k++) {
              pat::PackedGenParticle* p = (pat::PackedGenParticle*) &((*PackedGenParticlesHandle)[k]);

              if(p->mother(0) == D0) {
                if (verbose) {
                  cout << "\t\t" << k << " " << p->pdgId() << "  " << p->pt() << "  " << p->eta() << "  " << p->phi() << endl;
                }

                if ( abs(p->pdgId()) == 211) pi = p;
                else if ( abs(p->pdgId()) == 321 ) K = p;

                auto i_match = RecoMC_matching(pfCandHandle, p);
                if (i_match.size() >= 1) {
                  if(verbose) {
                    const reco::Candidate & pr = (*pfCandHandle)[i_match[0]];
                    cout << "\t" << i_match[0] << " " << pr.pdgId() << "  " << pr.pt() << "  " << pr.eta() << "  " << pr.phi() << endl;
                    cout << "\t\t" << pr.vz() << endl;
                  }
                  if ( abs(p->pdgId()) == 211) pi_reco = (reco::Candidate *) & ((*pfCandHandle)[i_match[0]]);
                  if ( abs(p->pdgId()) == 321) K_reco = (reco::Candidate *) & ((*pfCandHandle)[i_match[0]]);
                }
              }
            }
          }
        }
        if (pi!=0 && K!=0) {
          if (verbose) {cout << "\tIt's a D0toKpi" << endl;}
          D0toKpi = 1;
        }

        Bmeson.addUserInt("B_pdgID", p.pdgId());
        Bmeson.addUserFloat("B_pt", p.pt());
        Bmeson.addUserFloat("B_eta", p.eta());
        Bmeson.addUserFloat("B_phi", p.phi());
        Bmeson.addUserInt("B_D0munu", BtoD0munu);

        if(BtoD0munu) {
          Bmeson.addUserInt("D0_pdgID", D0->pdgId());
          Bmeson.addUserFloat("D0_pt", D0->pt());
          Bmeson.addUserFloat("D0_eta", D0->eta());
          Bmeson.addUserFloat("D0_phi", D0->phi());
          Bmeson.addUserFloat("D0_dz", D0->vz());

          Bmeson.addUserInt("mu_pdgID", mu->pdgId());
          Bmeson.addUserFloat("mu_pt", mu->pt());
          Bmeson.addUserFloat("mu_eta", mu->eta());
          Bmeson.addUserFloat("mu_phi", mu->phi());
          Bmeson.addUserFloat("mu_dz", mu->vz());

          Bmeson.addUserInt("nu_pdgID", nu->pdgId());
          Bmeson.addUserFloat("nu_pt", nu->pt());
          Bmeson.addUserFloat("nu_eta", nu->eta());
          Bmeson.addUserFloat("nu_phi", nu->phi());
        }
        else {
          Bmeson.addUserInt("D0_pdgID", 0);
          Bmeson.addUserFloat("D0_pt", 0);
          Bmeson.addUserFloat("D0_eta", 0);
          Bmeson.addUserFloat("D0_phi", 0);
          Bmeson.addUserFloat("D0_dz", -999999);
          Bmeson.addUserInt("mu_pdgID", 0);
          Bmeson.addUserFloat("mu_pt", 0);
          Bmeson.addUserFloat("mu_eta", 0);
          Bmeson.addUserFloat("mu_phi", 0);
          Bmeson.addUserFloat("mu_dz", -999999);
          Bmeson.addUserInt("nu_pdgID", 0);
          Bmeson.addUserFloat("nu_pt", 0);
          Bmeson.addUserFloat("nu_eta", 0);
          Bmeson.addUserFloat("nu_phi", 0);
        }

        if(mu_reco!= 0) {
          Bmeson.addUserFloat("mu_reco_pt", mu_reco->pt());
          Bmeson.addUserFloat("mu_reco_eta", mu_reco->eta());
          Bmeson.addUserFloat("mu_reco_phi", mu_reco->phi());
          Bmeson.addUserFloat("mu_reco_dz", mu_reco->vz());
        }
        else {
          Bmeson.addUserFloat("mu_reco_pt", 0);
          Bmeson.addUserFloat("mu_reco_eta", 0);
          Bmeson.addUserFloat("mu_reco_phi", 0);
          Bmeson.addUserFloat("mu_reco_dz", -999999);
        }

        Bmeson.addUserInt("D0_Kpi", D0toKpi);
        if (D0toKpi) {
          Bmeson.addUserInt("pi_pdgID", pi->pdgId());
          Bmeson.addUserFloat("pi_pt", pi->pt());
          Bmeson.addUserFloat("pi_eta", pi->eta());
          Bmeson.addUserFloat("pi_phi", pi->phi());
          Bmeson.addUserFloat("pi_dz", pi->vz());

          Bmeson.addUserInt("K_pdgID", K->pdgId());
          Bmeson.addUserFloat("K_pt", K->pt());
          Bmeson.addUserFloat("K_eta", K->eta());
          Bmeson.addUserFloat("K_phi", K->phi());
          Bmeson.addUserFloat("K_dz", K->vz());
        }
        else {
          Bmeson.addUserInt("pi_pdgID", 0);
          Bmeson.addUserFloat("pi_pt", 0);
          Bmeson.addUserFloat("pi_eta", 0);
          Bmeson.addUserFloat("pi_phi", 0);
          Bmeson.addUserFloat("pi_dz", -999999);
          Bmeson.addUserInt("K_pdgID", 0);
          Bmeson.addUserFloat("K_pt", 0);
          Bmeson.addUserFloat("K_eta", 0);
          Bmeson.addUserFloat("K_phi", 0);
          Bmeson.addUserFloat("K_dz", -999999);
        }

        if(pi_reco!=0) {
          Bmeson.addUserFloat("pi_reco_pt", pi_reco->pt());
          Bmeson.addUserFloat("pi_reco_eta", pi_reco->eta());
          Bmeson.addUserFloat("pi_reco_phi", pi_reco->phi());
          Bmeson.addUserFloat("pi_reco_dz", pi_reco->vz());
        }
        else {
          Bmeson.addUserFloat("pi_reco_pt", 0);
          Bmeson.addUserFloat("pi_reco_eta", 0);
          Bmeson.addUserFloat("pi_reco_phi", 0);
          Bmeson.addUserFloat("pi_reco_dz", -999999);
        }

        if(K_reco!=0) {
          Bmeson.addUserFloat("K_reco_pt", K_reco->pt());
          Bmeson.addUserFloat("K_reco_eta", K_reco->eta());
          Bmeson.addUserFloat("K_reco_phi", K_reco->phi());
          Bmeson.addUserFloat("K_reco_dz", K_reco->vz());
        }
        else {
          Bmeson.addUserFloat("K_reco_pt", 0);
          Bmeson.addUserFloat("K_reco_eta", 0);
          Bmeson.addUserFloat("K_reco_phi", 0);
          Bmeson.addUserFloat("K_reco_dz", -999999);
        }

        result->push_back(Bmeson);
      }
    }

    iEvent.put(std::move(result));


    if (verbose) {cout << "--------------\n" << endl;}

}

vector<int> BPH_MCTruth::RecoMC_matching(edm::Handle<edm::View<reco::Candidate>> reco_handle, reco::Candidate * pg) {
  double max_DeltaR = 0.01;
  double max_Delta_pt_rel = 0.05;

  vector<int> out;

  unsigned int N = reco_handle->size();
  for (unsigned int k = 0; k < N; ++k) {
    const reco::Candidate & pr = (*reco_handle)[k];
    double dEta = pr.eta() - pg->eta();
    double dPhi = pr.phi() - pg->phi();
    double deltaR = sqrt(dEta*dEta + dPhi*dPhi);

    double dpt_rel = abs(pr.pt() - pg->pt())/pg->pt();

    if (dpt_rel < max_Delta_pt_rel && deltaR < max_DeltaR) {
      out.push_back(k);
      if(verbose) {cout << "\t\tdeltaR: " << deltaR << "   dpt_rel: " << dpt_rel << endl;}
    }
  }

  return out;
}

vector<int> BPH_MCTruth::RecoMC_matching(edm::Handle<edm::View<reco::Candidate>> reco_handle, pat::PackedGenParticle* pg) {
  double max_DeltaR = 0.01;
  double max_Delta_pt_rel = 0.05;

  vector<int> out;

  unsigned int N = reco_handle->size();
  for (unsigned int k = 0; k < N; ++k) {
    const reco::Candidate & pr = (*reco_handle)[k];
    double dEta = pr.eta() - pg->eta();
    double dPhi = pr.phi() - pg->phi();
    double deltaR = sqrt(dEta*dEta + dPhi*dPhi);

    double dpt_rel = abs(pr.pt() - pg->pt())/pg->pt();

    if (dpt_rel < max_Delta_pt_rel && deltaR < max_DeltaR) {
      out.push_back(k);
      if(verbose) {cout << "\t\tdeltaR: " << deltaR << "   dpt_rel: " << dpt_rel << endl;}
    }
  }

  return out;
}


DEFINE_FWK_MODULE(BPH_MCTruth);
