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



        map<string, double> settingsCLN;
        for(auto n: parNameCLN) settingsCLN["delta_" + n] = 0;
        hammer.setFFEigenvectors("BtoD*", "CLNVar", settingsCLN);
        outRate["CLNCentral"] = hammer.getRate(proc, "SchmeCLN");
        if(verbose) { cout << Form("CLN central rate: %1.3e (ratio = %.3f)", outRate["CLNCentral"], outRate["CLNCentral"]/outRate["den"]) << endl;}
        // hammer.saveOptionCard("Opts.yml", false);

        for(int i=0; i<4; i++) { //Loop over eigenVar
          for (int j=0; j<2; j++) { //Loop over pos/neg direction
            map<string, double> settings;
            for (int k=0; k<4; k++) { //Loop over parameters
              settings["delta_" + parNameCLN[k]] = eigVarCLN[i][k][j];
            }

            hammer.setFFEigenvectors("BtoD*", "CLNVar", settings);
            auto rate = hammer.getRate(proc, "SchmeCLN");
            string var_name = "CLN" + varNameCLN[i];
            var_name += j==0? "Up" : "Down";
            outRate[var_name] = rate;
            if(verbose) {cout << var_name << Form(": %1.3e", rate) << endl;}
          }
        }



        map<string, double> settingsBLPR;
        for(auto n: parNameBLPR) settingsBLPR["delta_" + n] = 0;
        hammer.setFFEigenvectors("BtoD*", "BLPRVar", settingsBLPR);
        outRate["BLPRCentral"] = hammer.getRate(proc, "SchmeBLPR");
        if(verbose) { cout << Form("BLPR central rate: %1.3e (ratio = %.3f)", outRate["BLPRCentral"], outRate["BLPRCentral"]/outRate["den"]) << endl;}
        // hammer.saveOptionCard("Opts.yml", false);

        for(int i=0; i<7; i++) { //Loop over eigenVar
          for (int j=0; j<2; j++) { //Loop over pos/neg direction
            map<string, double> settings;
            for (int k=0; k<7; k++) { //Loop over parameters
              settings["delta_" + parNameBLPR[k]] = eigVarBLPR[i][k][j];
            }

            hammer.setFFEigenvectors("BtoD*", "BLPRVar", settings);
            auto rate = hammer.getRate(proc, "SchmeBLPR");
            string var_name = "BLPR" + varNameBLPR[i];
            var_name += j==0? "Up" : "Down";
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

    // ################ CLN parameters #############
    const vector<string> parNameCLN = {"RhoSq", "R1", "R2", "R0"};
    // rho2, R1 and R2 from Ex: https://arxiv.org/pdf/1909.12524.pdf
    // R0 from Th: https://journals.aps.org/prd/pdf/10.1103/PhysRevD.85.094025
    const double centralValCLN[4] = {1.122, 1.270, 0.852, 1.14};

    const vector<string> varNameCLN = {"eig1", "eig2", "eig3", "R0"};
    const double eigVarCLN[4][4][2] = {
      {{-0.02102485,  0.02102485}, {-0.02293032,  0.02293032}, { 0.01652843, -0.01652843}, {0, 0}},
      {{-0.01105955,  0.01105955}, { 0.01215074, -0.01215074}, { 0.00278883, -0.00278883}, {0, 0}},
      {{ 0.00341205, -0.00341205}, { 0.00159999, -0.00159999}, { 0.00655998, -0.00655998}, {0, 0}},
      {{0, 0}, {0, 0}, {0, 0}, {+0.11, -0.11}}
    };


    // ########## BLPR parameters ##############
    // Form arXiv:1703.05330v4, using scheme NoL+SR (No Lattice QCD, Yes QCD sum rules)
    const vector<string> parNameBLPR = {"RhoSq","chi21","chi2p","chi3p","eta1","etap","dV20"};
    const double centralValBLPR[7] = {1.19, -0.06, -0.00, 0.04, 0.35, -0.11, 0.};
    const vector<string> varNameBLPR = {"eig1", "eig2", "eig3", "eig4", "eig5", "eig6", "dV20"};
    const double eigVarBLPR[7][7][2] = {
      {{ 0.02322283, -0.02322283}, { 0.00022429, -0.00022429}, { 0.00023267, -0.00023267}, {-0.00079148,  0.00079148}, { 0.09664184, -0.09664184}, {-0.16885149,  0.16885149}, {0., 0.}},
      {{-0.04345272,  0.04345272}, { 0.00474643, -0.00474643}, { 0.00233232, -0.00233232}, {-0.00809314,  0.00809314}, {-0.09813503,  0.09813503}, {-0.06209618,  0.06209618}, {0., 0.}},
      {{ 0.06298508, -0.06298508}, {-0.00056467,  0.00056467}, {-0.00022681,  0.00022681}, { 0.00975046, -0.00975046}, {-0.02509061,  0.02509061}, {-0.00574474,  0.00574474}, {0., 0.}},
      {{ 0.002127,   -0.002127  }, { 0.00565397, -0.00565397}, { 0.00281448, -0.00281448}, {-0.01267561,  0.01267561}, { 0.0001571,  -0.0001571 }, { 0.00045325, -0.00045325}, {0., 0.}},
      {{-0.00073498,  0.00073498}, { 0.00643962, -0.00643962}, { 0.01769979, -0.01769979}, { 0.00668928, -0.00668928}, { 0.00041761, -0.00041761}, { 0.00013952, -0.00013952}, {0., 0.}},
      {{ 0.0006135,  -0.0006135 }, {-0.01742602,  0.01742602}, { 0.00855866, -0.00855866}, {-0.00577802,  0.00577802}, {-0.00036553,  0.00036553}, {-0.0001091,   0.0001091 }, {0., 0.}},
      {{0., 0.}, {0., 0.}, {0., 0.}, {0., 0.}, {0., 0.}, {0., 0.}, { 0.0001, -0.0001}}
    };



    int N_evets_weights_produced = 0;
    int N_evets_analyzed = 0;
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
    hammer.addFFScheme("SchmeCLN", {{"BD*", "CLNVar"}});
    hammer.addFFScheme("SchmeBLPR", {{"BD*", "BLPRVar"}});
    hammer.initRun();

    string centralValuesOpt = "BtoD*CLN: {";
    for(auto i=0; i<4; i++) {
      centralValuesOpt += Form("%s: %f, ", parNameCLN[i].c_str(), centralValCLN[i]);
    }
    centralValuesOpt += "}";
    if (verbose) {cout << "[Hammer]: CLN central values\n\t" << centralValuesOpt << endl;}
    hammer.setOptions(centralValuesOpt);

    centralValuesOpt = "BtoD*BLPR: {";
    for(auto i=0; i<7; i++) {
      centralValuesOpt += Form("%s: %f, ", parNameBLPR[i].c_str(), centralValBLPR[i]);
    }
    centralValuesOpt += "}";
    if (verbose) {cout << "[Hammer]: BLPR central values\n\t" << centralValuesOpt << endl;}
    hammer.setOptions(centralValuesOpt);

    produces<map<string, float>>("outputNtuplizer");
}


void HammerWeightsProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    if (verbose) {cout << "-----------  Hammer weights ----------\n";}
    N_evets_analyzed++;
    // Get prunedGenParticles
    edm::Handle<std::vector<reco::GenParticle>> PrunedGenParticlesHandle;
    iEvent.getByToken(PrunedParticlesSrc_, PrunedGenParticlesHandle);

    // Output collection
    unique_ptr<map<string, float>> outputNtuplizer(new map<string, float>);

    // Get the B index
    edm::Handle<int> indexBmcHandle;
    iEvent.getByToken(indexBmcSrc_, indexBmcHandle);
    int i_B = (*indexBmcHandle);
    if(i_B == -1){
      cout << "[ERROR]: Invalid B idx (i.e. no B MC set)" << endl;
      cerr << "[ERROR]: Invalid B idx (i.e. no B MC set)" << endl;
      (*outputNtuplizer)["wh_CLNCentral"] = 1e-9;
      for(int i=0; i<4; i++) { //Loop over eigenVar
        for (int j=0; j<2; j++) { //Loop over pos/neg direction
          string var_name = "CLN" + varNameCLN[i];
          var_name += j==0? "Up" : "Down";
          (*outputNtuplizer)["wh_" + var_name] = 0;
        }
      }

      (*outputNtuplizer)["wh_BLPRCentral"] = 1e-9;
      for(int i=0; i<7; i++) { //Loop over eigenVar
        for (int j=0; j<2; j++) { //Loop over pos/neg direction
          string var_name = "BLPR" + varNameBLPR[i];
          var_name += j==0? "Up" : "Down";
          (*outputNtuplizer)["wh_" + var_name] = 0;
        }
      }

      iEvent.put(move(outputNtuplizer), "outputNtuplizer");

      if (2*N_evets_weights_produced < N_evets_analyzed) return;
      else exit(1);
    }
    // cout << "i_B retieved: " << i_B << endl;

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


    map<string, double> settingsCLN;
    for(auto n: parNameCLN) settingsCLN["delta_" + n] = 0;
    hammer.setFFEigenvectors("BtoD*", "CLNVar", settingsCLN);
    auto weightsCLN = hammer.getWeights("SchmeCLN");
    if(!weightsCLN.empty()) {
      if(verbose) {cout << "CLNCentral: " << flush;}
      for(auto elem: weightsCLN) {
        if(isnan(elem.second)) {
          cout << "[ERROR]: CLNCentral nan weight: " << elem.second << endl;
          cerr << "[ERROR]: CLNCentral nan weight: " << elem.second << endl;
          assert(false);
        }
        (*outputNtuplizer)["wh_CLNCentral"] = elem.second;
        if(verbose) {cout << elem.second << endl;}
      }

      for(int i=0; i<4; i++) { //Loop over eigenVar
        for (int j=0; j<2; j++) { //Loop over pos/neg direction
          map<string, double> settings;
          for (int k=0; k<4; k++) { //Loop over parameters
            settings["delta_" + parNameCLN[k]] = eigVarCLN[i][k][j];
          }

          hammer.setFFEigenvectors("BtoD*", "CLNVar", settings);
          auto weights = hammer.getWeights("SchmeCLN");
          string var_name = "CLN" + varNameCLN[i];
          var_name += j==0? "Up" : "Down";

          if(verbose) {cout << var_name << ": " << flush;}
          for(auto elem: weights) {
            (*outputNtuplizer)["wh_" + var_name] = elem.second;
            if(verbose) {cout << elem.second << Form(" (Ratio to central = %.3f)", (*outputNtuplizer)["wh_" + var_name]/(*outputNtuplizer)["wh_CLNCentral"])<< endl;}
          }
        }
      }
      if(verbose) {cout << endl;}
    }




    map<string, double> settingsBLPR;
    for(auto n: parNameBLPR) settingsBLPR["delta_" + n] = 0;
    hammer.setFFEigenvectors("BtoD*", "BLPRVar", settingsBLPR);
    auto weightsBLPR = hammer.getWeights("SchmeBLPR");
    if(!weightsBLPR.empty()) {
      if(verbose) {cout << "BLPRCentral: " << flush;}
      for(auto elem: weightsBLPR) {
        if(isnan(elem.second)) {
          cout << "[ERROR]: BLPRCentral nan weight: " << elem.second << endl;
          cerr << "[ERROR]: BLPRCentral nan weight: " << elem.second << endl;
          assert(false);
        }
        (*outputNtuplizer)["wh_BLPRCentral"] = elem.second;
        if(verbose) {cout << elem.second << endl;}
      }

      for(int i=0; i<7; i++) { //Loop over eigenVar
        for (int j=0; j<2; j++) { //Loop over pos/neg direction
          map<string, double> settings;
          for (int k=0; k<7; k++) { //Loop over parameters
            settings["delta_" + parNameBLPR[k]] = eigVarBLPR[i][k][j];
          }

          hammer.setFFEigenvectors("BtoD*", "BLPRVar", settings);
          auto weights = hammer.getWeights("SchmeBLPR");
          string var_name = "BLPR" + varNameBLPR[i];
          var_name += j==0? "Up" : "Down";

          if(verbose) {cout << var_name << ": " << flush;}
          for(auto elem: weights) {
            (*outputNtuplizer)["wh_" + var_name] = elem.second;
            if(verbose) {cout << elem.second << endl;}
          }
        }
      }
      if(verbose) {cout << endl;}
    }


    if(!weightsCLN.empty() && !weightsBLPR.empty()) N_evets_weights_produced++;
    else if (weightsCLN.empty() && weightsBLPR.empty()) return;
    else {
      cout << "CLN is empty: " << weightsCLN.empty() << endl;
      cout << "BLPR is empty: " << weightsBLPR.empty() << endl;
      exit(666);
    }


    iEvent.put(move(outputNtuplizer), "outputNtuplizer");
    return;
}


DEFINE_FWK_MODULE(HammerWeightsProducer);
