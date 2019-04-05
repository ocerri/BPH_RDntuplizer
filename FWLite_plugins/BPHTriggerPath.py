import ROOT as rt
import numpy as np
import re
# load FWLite C++ libraries
rt.gSystem.Load("libFWCoreFWLite.so");
rt.gSystem.Load("libDataFormatsFWLite.so");
rt.FWLiteEnabler.enable()

# load FWlite python libraries
from DataFormats.FWLite import Handle

handle = {}
handle['TrgBits'] = [Handle("edm::TriggerResults"), ("TriggerResults","","HLT")]
handle['TrgObj'] = [Handle("std::vector<pat::TriggerObjectStandAlone>"), "slimmedPatTrigger"]
handle['Muons'] = [Handle("std::vector<pat::Muon>"), "slimmedMuons"]


def TriggerObj_matching(muon_list, obj, charge, verbose=False):
    max_DeltaR = 0.02
    max_Delta_pt_rel = 0.1
    bestM_DeltaR = max_DeltaR

    out = []
    for mu in muon_list:
        if charge!= 0 and mu.charge() != charge:
            continue

        dEta = mu.eta() - obj.eta()
        dPhi = mu.phi() - obj.phi()
        while np.abs(dPhi) > np.pi:
            dPhi -= np.sign(dPhi)*2*np.pi

        deltaR = np.sqrt(dEta*dEta + dPhi*dPhi)
        dpt_rel = abs(mu.pt() - obj.pt())/obj.pt();

        if dpt_rel < max_Delta_pt_rel and deltaR < max_DeltaR:
            if verbose:
                print '\t\tMuon matched with deltaR={:1.1e} and dpt_rel={:1.1e}'.format(deltaR, dpt_rel)
            if deltaR <= bestM_DeltaR:
                bestM_DeltaR = deltaR
                out = [mu] + out
            else: out.append(mu)
    return out

class BPHTriggerPath:
    def __init__(self,
                 mu_charge=0,
                 filter=True,
                 hltpath_template='HLT_Mu[0-9]+_IP[0-9]+.*',
                 mu_collection='hlt.*MuonCandidates::HLT',
                 produce_output=True
                 ):
        self.filter = filter
        self.mu_charge = mu_charge
        self.path_tmp = hltpath_template
        self.mu_collection = mu_collection
        self.produce_output = produce_output

    def process(self, event, output, verbose):
        out = output.evt_out

        prods = {}
        for k,v in handle.iteritems():
            event.getByLabel(v[1], v[0])
            prods[k] = v[0].product()

        pathNames = event.object().triggerNames(prods['TrgBits'])

        Triggered = False
        event.offline_BPH_trg = []
        for to in prods['TrgObj']:
            # if self.mu_charge and to.charge() != self.mu_charge:
            #     continue
            to.unpackPathNames(pathNames)
            BPH_paths = []
            for path in to.pathNames(True):
                if not re.match(self.path_tmp, path) is None:
                    BPH_paths.append(path)

            aux = re.match(self.mu_collection, to.collection())
            collection_ok = False if aux is None else True

            if len(BPH_paths)>0 and collection_ok:
                Triggered = True

                offline_matching_mu = TriggerObj_matching(prods['Muons'], to, self.mu_charge, verbose)

                dic = {}
                dic['BPH_paths'] = BPH_paths
                dic['matched_mu'] = offline_matching_mu
                dic['TrgObj'] = to
                event.offline_BPH_trg.append(dic)

                if verbose:
                    print "Trigger object pt %6.2f eta %+5.3f phi %+5.3f  " % (to.pt(),to.eta(),to.phi())
                    print 'Collection:', to.collection()
                    print 'Paths:', BPH_paths

        if self.produce_output:
            out['BPH_NTrgObj'] = len(event.offline_BPH_trg)
            out['BPH_TrgObJ_matched'] = 0
            out['trgMu_pt'] = -1
            out['trgMu_eta'] = -1
            out['trgMu_phi'] = -1

            for d in event.offline_BPH_trg:
                if len(d['matched_mu']) > 0:
                    out['BPH_TrgObJ_matched']+=1
                    if out['trgMu_pt'] == -1:
                        m = d['matched_mu'][0]
                        out['trgMu_pt'] = m.pt()
                        out['trgMu_eta'] = m.eta()
                        out['trgMu_phi'] = m.phi()




        if self.filter:
            return Triggered
        else:
            if self.produce_output:
                out['BPH_trg'] = 1 if Triggered else 0
            return True
