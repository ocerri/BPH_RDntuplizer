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
# handle['Muons'] = [Handle('std::vector<pat::Muon>'), 'slimmedMuons']
# handle['Vtxs'] = [Handle('vector<reco::Vertex>'), 'offlineSlimmedPrimaryVertices']

def print_candidate(p, addon=''):
    print addon+'Pt: {:.2f} Eta: {:.2f} Phi: {:.2f}'.format(p.pt(), p.eta(), p.phi())

def deltaR(p1, p2):
    dphi = p1.phi()-p2.phi()
    while np.abs(dphi) > np.pi:
        dphi -= np.sign(dphi)*2*np.pi
    return np.hypot(dphi, p1.eta()-p2.eta())

def matchRECO2MC(p_MC, reco_collection):
    p_best = None
    dR_best = -1
    dpt_best = 0

    for p in reco_collection:
        if p.charge() == p_MC.charge():
            dR = deltaR(p, p_MC)
            dpt = (p.pt() - p_MC.pt())/ p_MC.pt()
            if (p_best is None) or dR < dR_best or (abs(dR-dR_best) < 1e-3 and dpt < dpt_best):
                p_best = p
                dR_best = dR
                dpt_best = dpt

    return [p_best, dR_best, dpt_best]

def addToOutput(p, tag, outmap):
    outmap[tag+'_pt'] = p.pt()
    outmap[tag+'_eta'] = p.eta()
    outmap[tag+'_phi'] = p.phi()
    outmap[tag+'_pdgId'] = p.pdgId()
    outmap[tag+'_dz'] = p.dz()
    outmap[tag+'_dxy'] = p.dxy()

    if p.hasTrackDetails():
        outmap[tag+'_dzError'] = p.dzError()
        outmap[tag+'_dxyError'] = p.dxyError()
        t = p.bestTrack()
        outmap[tag+'_Nchi2'] = t.normalizedChi2()
        outmap[tag+'_chi2'] = t.chi2()
        outmap[tag+'_ndof'] = t.ndof()
        outmap[tag+'_Nhits'] = t.numberOfValidHits()
    else:
        outmap[tag+'_dzError'] = -1
        outmap[tag+'_dxyError'] = -1
        outmap[tag+'_Nchi2'] = -1
        outmap[tag+'_chi2'] = -1
        outmap[tag+'_ndof'] = -1
        outmap[tag+'_Nhits'] = -1

def put_part(closest_particles, p, dR, dz, N_p_max):
    if len(closest_particles) == 0:
        closest_particles.append([p, dR, dz])
    else:
        idx = -1
        for i_s, cp in enumerate(closest_particles):
            if cp[1] > dR:
                idx = i_s
                break
        if idx == -1 and len(closest_particles) < N_p_max:
            closest_particles.append([p, dR, dz])
        else:
            closest_particles.insert(idx, [p, dR, dz])
            closest_particles = closest_particles[:-1]


class MuonNeighbourPUTracks:
    def __init__(self,
                 verbose=False,
                 N_part = 5
                ):
        self.verbose = verbose
        self.N_part = N_part

    def process(self, event, output, verbose):
        out = output.evt_out

        prods = {}
        for k,v in handle.iteritems():
            event.getByLabel(v[1], v[0])
            prods[k] = v[0].product()

        if verbose or self.verbose:
            print '------- Gathering the first particles around the muon ---------'

        pt_l = []
        eta_l = []
        phi_l = []
        for p in event.RECO_MCmatch.values():
            pt_l.append(p[0].pt())
            eta_l.append(p[0].eta())
            phi_l.append(p[0].phi())
        pt_l = np.array(pt_l)
        eta_l = np.array(eta_l)
        phi_l = np.array(phi_l)


        muon = event.RECO_MCmatch['mu'][0]

        closest_particles = []

        for p in prods['PFCand']:
            if np.min(pt_l - p.pt()) == 0 and np.min(eta_l - p.eta()) == 0 and np.min(phi_l - p.phi()) == 0:
                print 'Same found!'
                continue

            dz = np.abs(p.dz() - muon.dz())
            if dz < 2:
                dR = deltaR(muon, p)
                put_part(closest_particles, p, dR, dz, self.N_part)

        for ip, pc in enumerate(closest_particles):
            pn = 'pc'+str(ip)
            out[pn+'_pt'] = pc[0].pt()
            out[pn+'_eta'] = pc[0].eta()
            out[pn+'_phi'] = pc[0].phi()
            out[pn+'_dzmu'] = pc[2]
            out[pn+'_dRmu'] = pc[1]


        return True
