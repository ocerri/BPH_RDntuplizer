import argparse
import ROOT as rt
import numpy as np
import root_numpy as rtnp
from glob import glob

# load FWLite C++ libraries
rt.gSystem.Load("libFWCoreFWLite.so");
rt.gSystem.Load("libDataFormatsFWLite.so");
rt.FWLiteEnabler.enable()

# load FWlite python libraries
from DataFormats.FWLite import Handle, Events

#load specific moduels
from MC_Process import getTLorenzVector
from pdg_utils import getName


def parsing():
    parser = argparse.ArgumentParser()
    # parser.add_argument("-i", "--input_file", type=str, help="input root file, if -N is given XX is replaced with runNumber", nargs='+')

    parser.add_argument("-N", "--max_number", default=-1, help="max number of events", type=int)

    parser.add_argument("-B", "--batch", default=True, action='store_false', help="Root batch mode")
    parser.add_argument("-v", "--verbose", default=False, action='store_true', help="Activate verbose mode")

    args = parser.parse_args()
    return args

# files=['/afs/cern.ch/user/o/ocerri/cernbox/BPhysics/data/cmsMC_private/BPH_gammaTRIAL_Tag-Bm_D0kpmunu_Probe-Bp_D0kpmunu_GEN-SIM.root']
files=['/afs/cern.ch/user/o/ocerri/cernbox/BPhysics/data/cmsMC_private/BPH_gammaTRIAL_Tag-Bm_D0kpmunu_Probe-Bp_D0kpmunu_NoPU_10-2-3_v0/BPH_gammaTRIAL_Tag-Bm_D0kpmunu_Probe-Bp_D0kpmunu_GEN-SIM.root']

files = glob('/afs/cern.ch/user/o/ocerri/cernbox/BPhysics/data/cmsMC_private/BPH_Tag-Bm_D0kpmunu_Probe-Bp_D0kpmunu_NoPU_10-2-3_v0/jobs_out/*.root')

handleGenP  = Handle("std::vector<reco::GenParticle>")
# labelGenP = ('genParticles')
labelGenP = ('prunedGenParticles')


if __name__ == '__main__':
    rt.gErrorIgnoreLevel = 6000

    args = parsing()

    if args.batch:
        rt.gROOT.SetBatch()
    verbose = args.verbose

    try:
        n_Bp = np.array([0, 0, 0, 0, 0])
        n_Bm = np.array([0, 0, 0, 0, 0])

        n_Bp_cut = np.array([0, 0, 0, 0, 0])
        n_Bm_cut = np.array([0, 0, 0, 0, 0])

        out_tree = []

        for f in files:
            # open file (you can use 'edmFileUtil -d /store/whatever.root' to get the physical file name)
            # print "->Opening file",f.split()[0]
            events = Events(f.split()[0])


            for iev,event in enumerate(events):
                if args.max_number >= 0 and iev == args.max_number:
                    break

                # print '\nEvent {}: run {}, lumi {}, event {}'.format(iev,event.eventAuxiliary().run(), event.eventAuxiliary().luminosityBlock(),event.eventAuxiliary().event())

                event.getByLabel(labelGenP, handleGenP)
                genp = handleGenP.product()

                for p in genp:
                    if abs(p.pdgId()) == 521:
                    # if abs(p.pdgId()) > 500 and abs(p.pdgId()) % 500 < 50:
                        if verbose and 0 and iev==167:
                            out = 'Found: ' + getName(p.pdgId()) + '(-->'
                            for d in p.daughterRefVector():
                                out += ' '+getName(d.pdgId())
                            out += ')'
                            print out

                        # getTLorenzVector(p, verbose)
                        p4_sum = rt.TLorentzVector(0,0,0,0)
                        N_gamma = 0
                        sel = False
                        for d in p.daughterRefVector():
                            if d.pdgId() == 22:
                                if N_gamma == 0:
                                    gamma = d
                                elif gamma.pt() < d.pt():
                                    gamma = d

                                N_gamma += 1
                                continue

                            p4_sum += getTLorenzVector(d, 0)
                            if abs(d.pdgId()) == 13:
                                mu = d

                                d.this_sel=0
                                if d.pt()>6.5 and abs(d.eta()) < 2.5 :
                                    sel = True
                                    d.this_sel=1

                        e = 0
                        if N_gamma>0:
                            e = gamma.energy()

                        out = (mu.pt(), mu.eta(), mu.charge(), mu.this_sel, N_gamma, e, p4_sum.M() - p.mass())
                        out_tree.append(out)



                        # print 'Sum:', p4_sum.M()
                        # print '\n'

                        if p.pdgId()>0:
                            n_Bp[N_gamma] += 1
                        else:
                            n_Bm[N_gamma] += 1

                        if sel:
                            if p.pdgId()>0:
                                n_Bp_cut[N_gamma] += 1
                            else:
                                n_Bm_cut[N_gamma] += 1

                # if np.sum(n_Bp) != np.sum(n_Bm):
                #     print iev
                #     break

        print np.sum(n_Bp), n_Bp/float(np.sum(n_Bp))
        print np.sum(n_Bm), n_Bm/float(np.sum(n_Bm))
        print np.sum(n_Bp_cut), n_Bp_cut/float(np.sum(n_Bp_cut))
        print np.sum(n_Bm_cut), n_Bm_cut/float(np.sum(n_Bm_cut))

        fields = ['pt', 'eta', 'charge', 'trg', 'Ng', 'gEnergy', 'dM_B']
        dtypes = []
        for f in fields:
            dtypes.append((f, np.float))

        out_tree = np.array(out_tree, dtype=dtypes)
        rtnp.array2root(out_tree, 'gammatrial.root', mode='recreate')


    except KeyboardInterrupt:
        pass

    print "DONE"
