import ROOT as rt
import numpy as np
# load FWLite C++ libraries
rt.gSystem.Load("libFWCoreFWLite.so");
rt.gSystem.Load("libDataFormatsFWLite.so");
rt.FWLiteEnabler.enable()

# load FWlite python libraries
from DataFormats.FWLite import Handle

from pdg_utils import getName, getMeanLife

def isAncestor(a,p):
    if a == p:
        return True
    for i in xrange(0,p.numberOfMothers()):
        if isAncestor(a,p.mother(i)):
             return True
    return False

def getTLorenzVector(p, verbose=0):
    if verbose:
        print 'PtEtaPhiM:{:.3f}, {:.3f}, {:.3f}, {:.4f}'.format(p.pt(), p.eta(), p.phi(), p.mass())
    out =  rt.TLorentzVector()
    out.SetPtEtaPhiM(p.pt(), p.eta(), p.phi(), p.mass())
    return out

def getNuMomentum(d_flightB, p_vis, pt_B):
    d = d_flightB*(1/d_flightB.Mag())
    p_B = d * (pt_B/d.Perp())

    p_nu_vec = p_B - p_vis.Vect()
    p_nu = rt.TLorentzVector(p_nu_vec, p_nu_vec.Mag())
    return p_nu

handlePruned  = Handle("std::vector<reco::GenParticle>")
handlePacked  = Handle("std::vector<pat::PackedGenParticle>")
labelPruned = ("prunedGenParticles")
labelPacked = ("packedGenParticles")

c_light = 2.99792458e1 #cm/ns

class MC_Process:
    def __init__(self, probe_pdgId = 521, tag_pdgId = -521):
        self.probe_pdgId = probe_pdgId
        self.tag_pdgId = tag_pdgId

    def process(self, event, output, verbose):
        out = output.evt_out

        event.getByLabel(labelPacked, handlePacked)
        event.getByLabel(labelPruned, handlePruned)
        # get the product
        genp_packed = handlePacked.product()
        genp_pruned = handlePruned.product()

        Btag = None
        Bprobe = None

        for p in genp_pruned:
            if p.pdgId() in [self.probe_pdgId, self.tag_pdgId]:
            # if abs(p.pdgId()) > 500 and abs(p.pdgId()) % 500 < 50:
                if verbose:
                    msg = 'Found: ' + getName(p.pdgId()) + '(-->'
                    for d in p.daughterRefVector():
                        msg += ' '+getName(d.pdgId())
                    msg += ')'
                    print msg

                    getTLorenzVector(p, verbose)
                    for d in p.daughterRefVector():
                        print getName(d.pdgId())
                        getTLorenzVector(d, verbose)

                if p.pdgId() == self.probe_pdgId and p.numberOfDaughters()>1: #numberOfDaughters avoids B0 from oscillation
                    Bprobe = p
                else:
                    Btag = p


        if Btag is None or Bprobe is None:
            return False

        out['Btag_pt'] = Btag.pt()
        out['Bp_pt'] = Bprobe.pt()
        out['Bp_p'] = Bprobe.p()
        out['Bp_eta'] = Bprobe.eta()
        dphi = Bprobe.phi() - (Btag.phi()+np.pi)
        while dphi > np.pi:
            dphi -= 2*np.pi
        while dphi < -np.pi:
            dphi += 2*np.pi
        out['dphi_B'] = dphi

        vtx_prod = np.array([Bprobe.vx(), Bprobe.vy(), Bprobe.vz()])
        d = Bprobe.daughter(0)
        vtx_decay = np.array([d.vx(), d.vy(), d.vz()])
        d_flightB = vtx_decay - vtx_prod
        d_flightB = rt.TVector3(d_flightB[0], d_flightB[1], d_flightB[2])
        dl = np.linalg.norm(vtx_prod - vtx_decay)
        out['dl'] = dl
        out['Bp_tau'] = dl/c_light * Bprobe.mass()/Bprobe.p() #ns
        out['Bp_p_dl'] = Bprobe.mass() * dl / (getMeanLife(521) * np.log(2) * c_light)

        p4_vis = rt.TLorentzVector(0,0,0,0)
        N_neutral = 0

        def sumVisibleP4(p, vis_P4, N):
            # print getName(p.pdgId())
            if p.numberOfDaughters()>0:
                for d in p.daughterRefVector():
                    vis_P4, N = sumVisibleP4(d, vis_P4, N)
            elif p.charge() == 0:
                N += 1
                # print 'N+=1'
            else:
                vis_P4 += getTLorenzVector(p)
                # print 'Visible found, adding to p4'

            return vis_P4, N

        p4_vis, N_neutral = sumVisibleP4(Bprobe, p4_vis, N_neutral)

        out['N_neutral'] = N_neutral
        out['M_vis'] = p4_vis.M()
        out['M_miss'] = Bprobe.mass() - p4_vis.M()

        p4_Bprobe = getTLorenzVector(Bprobe)
        out['dTheta_visB'] = p4_vis.Angle(p4_Bprobe.Vect())

        #Find the D or D*
        for d in Bprobe.daughterRefVector():
            if abs(d.pdgId()) in [413, 421]:
                D = d
                break
            elif abs(d.pdgId()) == 423:
                for dd in d.daughterRefVector():
                    if abs(dd.pdgId()) == 421:
                        D = dd
                break

        #Find the muon
        for d in Bprobe.daughterRefVector():
            if abs(d.pdgId()) == 13:
                mu_probe = d
                break
            elif abs(d.pdgId()) == 15:
                for dd in d.daughterRefVector():
                    if abs(dd.pdgId()) == 13:
                        mu_probe = dd
                break
        out['mu_pt'] = mu_probe.pt()
        vtx_D = np.array([D.vx(), D.vy(), D.vz()])
        vtx_mu = np.array([mu_probe.vx(), mu_probe.vy(), mu_probe.vz()])

        dir_mu = np.array([mu_probe.px(), mu_probe.py(), mu_probe.pz()])
        dir_mu /= np.linalg.norm(dir_mu)

        out['mu_ip'] = np.linalg.norm(np.cross(dir_mu, vtx_mu-vtx_D))


        out['D_pt'] = D.pt()
        p4_D = getTLorenzVector(D)
        if d_flightB.Mag() > 0:
            out['D_pthat'] = p4_D.Pt(d_flightB*(1/d_flightB.Mag()))
            pnu_ext = getNuMomentum(d_flightB, p4_vis, Btag.pt())
            out['M_ext'] = (p4_vis + pnu_ext).M()
        else:
            out['D_pthat'] = -999
            out['M_ext'] = -999
            print '------------> Istantaneus decay'

        #Print muons info:
        print '----- MC info -----'
        for d in Btag.daughterRefVector():
            if abs(d.pdgId()) ==  13:
                print 'Tag Mu:'
                getTLorenzVector(d, True)
        print 'Probe Mu:'
        getTLorenzVector(mu_probe, True)

        return True
