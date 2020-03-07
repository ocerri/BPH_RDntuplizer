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
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"

#include <iostream>
#include <string>
#include <map>

#include "Hammer/Hammer.hh"
#include "Hammer/Process.hh"
#include "Hammer/Particle.hh"

#include "TLorentzVector.h"

using namespace std;

class HammerWeightsProducer : public edm::EDProducer {

public:

    explicit HammerWeightsProducer(const edm::ParameterSet &iConfig);

    ~HammerWeightsProducer() {
      if(verbose) { cout << "[Hammer]: Computing rates" << endl;}
      vector<string> processes = {"B0D*-MuNu", "B0D*-TauNu"};
      // Getting overall rates
      for(auto proc : processes) {
        map<string, double> outRate;
        if(verbose) { cout << "Process: " << proc << endl;}
        outRate["den"] = hammer.getDenominatorRate(proc);
        if(outRate["den"] == 0) {
          if(verbose) { cout << "Not evaluated, skipping" << endl;}
          continue;
        }
        if(verbose) { cout << Form("Denominator rate: %1.3e", outRate["den"]) << endl;}

        map<string, double> settings;
        for(auto pars: paramsCLN) {
          string name = "delta_" + pars.first;
          settings[name] = 0;
        }
        hammer.setFFEigenvectors("BtoD*", "CLNVar", settings);
        outRate["Central"] = hammer.getRate(proc, "SchmeCLN");
        if(verbose) { cout << Form("Central rate: %1.3e (ratio = %.3f)", outRate["Central"], outRate["Central"]/outRate["den"]) << endl;}

        for(auto elem: paramsCLN) {
          for(int i=1; i <= 2; i++) {
            map<string, double> settings;
            for(auto pars: paramsCLN) {
              string name = "delta_" + pars.first;
              if(pars.first == elem.first) {
                settings[name] = paramsCLN[pars.first][i];
              }
              else settings[name] = 0;
            }
            hammer.setFFEigenvectors("BtoD*", "CLNVar", settings);
            auto rate = hammer.getRate(proc, "SchmeCLN");
            string var_name = elem.first;
            var_name += i==1? "Up" : "Down";
            outRate[var_name] = rate;
            if(verbose) {cout << var_name << Form(": %1.3e", rate) << endl;}
          }
        }

        edm::Service<TFileService> fs;
        TTree* tree = fs->make<TTree>( "Trate", Form("Rates from Hammer for %s", proc.c_str()));
        for(auto& kv : outRate) {
          auto k = kv.first;
          tree->Branch(k.c_str(), &(outRate[k]));
        }
        tree->Fill();
        break;
      }
    };

private:

    virtual void produce(edm::Event&, const edm::EventSetup&);

    // ----------member data ---------------------------

    edm::EDGetTokenT<vector<reco::GenParticle>> PrunedParticlesSrc_;
    edm::EDGetTokenT<int> indexBmcSrc_;

    Hammer::Hammer hammer;

    double mass_B0 = 5.27961;
    double mass_mu = 0.1056583745;
    double mass_K = 0.493677;
    double mass_pi = 0.13957018;

    int verbose = 0;

    map<string, vector<double>> paramsCLN {
      // from EX: https://hflav-eos.web.cern.ch/hflav-eos/semi/summer16/html/ExclusiveVcb/exclBtoDstar.html
      {"RhoSq", {1.205, +0.026, -0.026}},
      {"R1", {1.404, +0.032, -0.032}},
      {"R2", {0.854, +0.020, -0.020}},
      // from TH: https://journals.aps.org/prd/pdf/10.1103/PhysRevD.85.094025
      {"R0", {1.14, +0.11, -0.11}},
    };

    int N_evets_weights_produced = 0;
};



HammerWeightsProducer::HammerWeightsProducer(const edm::ParameterSet &iConfig)
{
    PrunedParticlesSrc_ = consumes<vector<reco::GenParticle>>(edm::InputTag("prunedGenParticles"));
    indexBmcSrc_ = consumes<int> (edm::InputTag("MCpart", "indexBmc"));

    verbose = iConfig.getParameter<int>( "verbose" );

    hammer.setUnits("GeV");
    auto decayOfInterest = iConfig.getParameter<vector<string>>( "decayOfInterest" );
    for(auto s : decayOfInterest) {
      if(verbose) {cout << "[Hammer]: Including decay " << s << endl;}
      hammer.includeDecay(s);
    }

    auto inputFFScheme_ = iConfig.getParameter<vector<string>>("inputFFScheme");
    if(verbose) {cout << "[Hammer]: Input scheme" << endl;}
    map<string, string> inputFFScheme;
    for(uint i = 0; i < inputFFScheme_.size(); i++) {
      if(i%2 == 1) continue;
      inputFFScheme[inputFFScheme_[i]] = inputFFScheme_[i+1];
      if(verbose){cout << "\t" << inputFFScheme_[i] << ": " << inputFFScheme_[i+1] << endl;}
    }
    hammer.setFFInputScheme(inputFFScheme);

    hammer.addFFScheme("SchmeCLN", {
                       // {"BD", "ISGW2"},
                       {"BD*", "CLNVar"}
                     });

    hammer.initRun();
    string centralValuesOpt = "BtoD*CLN: {";
    for(auto elem: paramsCLN) {
      centralValuesOpt += Form("%s: %f, ", elem.first.c_str(), elem.second[0]);
    }
    centralValuesOpt += "}";
    if (verbose) {cout << "[Hammer]: Central values\n\t" << centralValuesOpt << endl;}
    hammer.setOptions(centralValuesOpt);

    produces<map<string, float>>("outputNtuplizer");
}


void HammerWeightsProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    if (verbose) {cout << "-----------  Hammer weights ----------\n";}

    // Get prunedGenParticles
    edm::Handle<std::vector<reco::GenParticle>> PrunedGenParticlesHandle;
    iEvent.getByToken(PrunedParticlesSrc_, PrunedGenParticlesHandle);

    // Get the B index
    edm::Handle<int> indexBmcHandle;
    iEvent.getByToken(indexBmcSrc_, indexBmcHandle);
    int i_B = (*indexBmcHandle);
    if(i_B == -1){
      cout << "[ERROR]: Invalid B idx (i.e. no B MC set)" << endl;
      cerr << "[ERROR]: Invalid B idx (i.e. no B MC set)" << endl;
      if (N_evets_weights_produced == 0) return;
      else exit(1);
    }
    // cout << "i_B retieved: " << i_B << endl;

    // Output collection
    unique_ptr<map<string, float>> outputNtuplizer(new map<string, float>);

    // Initialize the Hammer event
    hammer.initEvent();
    Hammer::Process B2DstLNu_Dst2DPi;
    vector<size_t> Bvtx_idxs;
    int idxTau = -1;
    vector<size_t> Tauvtx_idxs;
    int idxDst = -1;
    vector<size_t> Dstvtx_idxs;

    TLorentzVector p4B;
    vector<TLorentzVector> pDauB;
    vector<TLorentzVector> pDauDst;
    bool DstFound = false;
    int idxDst_v = -1;

    auto p = (*PrunedGenParticlesHandle)[i_B];
    Hammer::Particle pB({p.energy(), p.px(), p.py(), p.pz()}, p.pdgId());
    if(verbose) {
      p4B.SetPxPyPzE(p.px(), p.py(), p.pz(), p.energy());
      cout << Form("B --> E:%.5f, Px:%.5f, Py:%.5f, Pz:%.5f PDG:%d", p.energy(), p.px(), p.py(), p.pz(), p.pdgId()) << endl;
    }
    auto idxB = B2DstLNu_Dst2DPi.addParticle(pB);
    for(auto d : p.daughterRefVector()) {
      Hammer::Particle B_dau({d->energy(), d->px(), d->py(), d->pz()}, d->pdgId());
      if(verbose) {
        TLorentzVector v;
        v.SetPxPyPzE(d->px(), d->py(), d->pz(), d->energy());
        pDauB.push_back(v);
        if(!DstFound){
          idxDst_v++;
          if (abs(d->pdgId()) == 413) DstFound = true;
        }
        cout << Form("%d --> E:%.5f, Px:%.5f, Py:%.5f, Pz:%.5f", d->pdgId(), d->energy(), d->px(), d->py(), d->pz()) << endl;
      }
      auto idx_d = B2DstLNu_Dst2DPi.addParticle(B_dau);
      Bvtx_idxs.push_back(idx_d);

      if(d->pdgId() == -15) {
        idxTau = idx_d;
        for (auto dd : d->daughterRefVector()) {
          if(verbose) { cout << Form("\t%d --> E:%.5f, Px:%.5f, Py:%.5f, Pz:%.5f", dd->pdgId(), dd->energy(), dd->px(), dd->py(), dd->pz()) << endl;}
          Hammer::Particle Tau_dau({dd->energy(), dd->px(), dd->py(), dd->pz()}, dd->pdgId());
          auto idx_dd = B2DstLNu_Dst2DPi.addParticle(Tau_dau);
          Tauvtx_idxs.push_back(idx_dd);
        }
      }
      else if(d->pdgId() == -413) {
        idxDst = idx_d;
        for (auto dd : d->daughterRefVector()) {
          if(verbose) {
            TLorentzVector v;
            v.SetPxPyPzE(dd->px(), dd->py(), dd->pz(), dd->energy());
            pDauDst.push_back(v);
            cout << Form("\t%d --> E:%.5f, Px:%.5f, Py:%.5f, Pz:%.5f", dd->pdgId(), dd->energy(), dd->px(), dd->py(), dd->pz()) << endl;
          }
          Hammer::Particle Dst_dau({dd->energy(), dd->px(), dd->py(), dd->pz()}, dd->pdgId());
          auto idx_dd = B2DstLNu_Dst2DPi.addParticle(Dst_dau);
          Dstvtx_idxs.push_back(idx_dd);
        }
      }
    }

    if(verbose) {
      TLorentzVector vtxB;
      vtxB = p4B;
      for(auto v : pDauB) {
        cout << "Mass: " << v.M() << endl;
        vtxB -= v;
      }
      cout << "B vertex coservation " << endl;
      vtxB.Print();

      TLorentzVector vtxDst;
      vtxDst = pDauB[idxDst_v];
      for(auto v : pDauDst) {
        cout << "Mass: " << v.M() << endl;
        vtxDst -= v;
      }
      cout << "Dst vertex coservation " << endl;
      vtxDst.Print();
    }

    B2DstLNu_Dst2DPi.addVertex(idxB, Bvtx_idxs);
    if(idxTau != -1) {
      B2DstLNu_Dst2DPi.addVertex(idxTau, Tauvtx_idxs);
    }
    if(idxDst != -1) {
      B2DstLNu_Dst2DPi.addVertex(idxDst, Dstvtx_idxs);
    }

    hammer.addProcess(B2DstLNu_Dst2DPi);
    hammer.processEvent();


    map<string, double> settings;
    for(auto pars: paramsCLN) {
      string name = "delta_" + pars.first;
      settings[name] = 0;
    }
    hammer.setFFEigenvectors("BtoD*", "CLNVar", settings);
    auto weights = hammer.getWeights("SchmeCLN");
    if(weights.empty()) return;
    else N_evets_weights_produced++;

    if(verbose) {cout << "CLNCentral: " << flush;}
    for(auto elem: weights) {
      if(isnan(elem.second)) {
        cout << "[ERROR]: CLNCentral nan weight: " << elem.second << endl;
        cerr << "[ERROR]: CLNCentral nan weight: " << elem.second << endl;
        assert(false);
      }
      (*outputNtuplizer)["wh_CLNCentral"] = elem.second;
      if(verbose) {cout << elem.second << endl;}
    }

    for(auto elem: paramsCLN) {
      for(int i=1; i <= 2; i++) {
        map<string, double> settings;
        for(auto pars: paramsCLN) {
          string name = "delta_" + pars.first;
          if(pars.first == elem.first) {
            settings[name] = paramsCLN[pars.first][i];
          }
          else settings[name] = 0;
        }
        hammer.setFFEigenvectors("BtoD*", "CLNVar", settings);
        auto weights = hammer.getWeights("SchmeCLN");
        string var_name = "CLN" + elem.first;
        var_name += i==1? "Up" : "Down";

        if(verbose) {cout << var_name << ": " << flush;}
        for(auto elem: weights) {
          (*outputNtuplizer)["wh_" + var_name] = elem.second;
          if(verbose) {cout << elem.second << endl;}
        }
      }
    }
    if(verbose) {cout << endl;}

    iEvent.put(move(outputNtuplizer), "outputNtuplizer");
    return;
}


DEFINE_FWK_MODULE(HammerWeightsProducer);
