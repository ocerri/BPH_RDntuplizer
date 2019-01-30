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
handle['PFCand'] = [Handle('vector<pat::PackedCandidate>'), 'packedPFCandidates']
handle['Muons'] = [Handle('std::vector<pat::Muon>'), 'slimmedMuons']

def print_candidate(p):
    print 'Pt: {:.2f} Eta: {:.2f} Phi: {:.2f}'.format(p.pt(), p.eta(), p.phi())


class ReCo_Decay_B0_DstMu_D0piMu_KpipiMu:
    def __init__(self,
                 ):
        pass

    def process(self, event, output, verbose):
        out = output.evt_out

        prods = {}
        for k,v in handle.iteritems():
            event.getByLabel(v[1], v[0])
            prods[k] = v[0].product()

        print '---- TrgObj -----'
        for to in event.offline_BPH_trg:
            print 'Matched muons:', len(to['matched_mu'])
            for tm in to['matched_mu']:
                print_candidate(tm)

        print '-------- Muons ---------'
        for mu in prods['Muons']:
            print_candidate(mu)

        print '------- PF cand -------'
        N = 0
        for p in prods['PFCand']:
            if p.pdgId() != -211:
                continue
            N += 1

        print 'N: ', N

        print '\n\n\n'
        return True
