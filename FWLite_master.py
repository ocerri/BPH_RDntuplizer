import argparse
import ROOT as rt
import numpy as np
import root_numpy as rtnp
import imp
from time import time

# load FWLite C++ libraries
rt.gSystem.Load('libFWCoreFWLite.so');
rt.gSystem.Load('libDataFormatsFWLite.so');
rt.FWLiteEnabler.enable()

# load FWlite python libraries
from DataFormats.FWLite import Events

def parsing():
    parser = argparse.ArgumentParser()

    parser.add_argument('config', help='Config path')

    parser.add_argument('-N', '--max_accepted', default=-1, help='max number of events', type=int)
    parser.add_argument('-B', '--batch', default=True, action='store_false', help='Root batch mode')
    parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Activate verbose mode')

    args = parser.parse_args()
    return args

class OutputBauble:
    def __init__(self, outfile):
        self.global_out_list = []
        self.evt_out = {}
        self.out_names = None
        self.outfile = outfile

    def clear_evt_out(self):
        self.evt_out = {}

    def append_to_global_out(self):
        if self.out_names is None:
            self.out_names = self.evt_out.keys()

        aux = []
        for k in self.out_names:
            if k in self.evt_out.keys():
                aux.append(self.evt_out[k])
            else:
                aux.append(-99999)

        self.global_out_list.append(tuple(aux))

    def dump_to_tree(self):
        print 'Writing output form {} events'.format(len(self.global_out_list))
        print 'Writing output file:'
        print self.outfile

        if (not self.out_names is None) and len(self.out_names)>0 and len(self.global_out_list)>0:
            # idx_alphab_order = np.argsort(self.out_names)
            # self.out_names = self.out_names[idx_alphab_order]
            # for i in len(self.global_out_list):
            #     self.global_out_list[i] = self.global_out_list[i][idx_alphab_order]

            dtypes = [(f, np.float) for f in self.out_names]
            aux = np.array(self.global_out_list, dtype=dtypes)

            rtnp.array2root(aux, self.outfile, mode='recreate')
            print 'Done'
        else:
            if len(self.out_names) == 0:
                print 'No fields to be written'

            print 'No output produced'


if __name__ == '__main__':
    rt.gErrorIgnoreLevel = 6000

    args = parsing()

    if args.batch:
        rt.gROOT.SetBatch()
    verbose = args.verbose

    cfg = imp.load_source('cfg', args.config)

    output = OutputBauble(cfg.outfile)

    try:
        t_last_print = -99999999
        N_accepted = 0
        N_processed = 0
        stop = False

        for ifl, f in enumerate(cfg.files):
            # open file (you can use 'edmFileUtil -d /store/whatever.root' to get the physical file name)
            if verbose:
                print '->Opening file',f.split()[0]
            events = Events(f.split()[0])

            for iev,event in enumerate(events):
                dir(event)
                output.clear_evt_out()

                if verbose:
                    print '\nEvent {}: run {}, lumi {}, event {}'.format(iev,event.eventAuxiliary().run(), event.eventAuxiliary().luminosityBlock(),event.eventAuxiliary().event())
                elif time() - t_last_print > 5 and iev%100==0:
                    t_last_print = time()
                    print 'File {}/{} - Evt {}'.format(ifl+1, len(cfg.files), iev)

                for module in cfg.exe_seq:
                    accepted = module.process(event, output, verbose)
                    if not accepted:
                        break
                N_processed += 1
                if accepted:
                    N_accepted += 1
                    output.append_to_global_out()

                if args.max_accepted >= 0 and N_accepted == args.max_accepted:
                    stop = True

                if stop: break

            if stop: break


        print 'File {}/{} - Evt {}'.format(ifl+1, len(cfg.files), iev)

        print 'Processed events', N_processed
        print 'Accepted events', N_accepted
        output.dump_to_tree()

    except KeyboardInterrupt:
        output.dump_to_tree()
