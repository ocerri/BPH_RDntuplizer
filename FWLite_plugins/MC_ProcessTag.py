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

def addToOutput(p, tag, outmap):
    outmap[tag+'_pt'] = p.pt()
    outmap[tag+'_eta'] = p.eta()
    outmap[tag+'_phi'] = p.phi()
    outmap[tag+'_pdgId'] = p.pdgId()

handlePruned  = Handle("std::vector<reco::GenParticle>")
handlePacked  = Handle("std::vector<pat::PackedGenParticle>")
labelPruned = ("prunedGenParticles")
labelPacked = ("packedGenParticles")

c_light = 2.99792458e1 #cm/ns

class MC_ProcessTag:
    def __init__(self, tag_pdgId = -521):
        self.tag_pdgId = tag_pdgId

    def process(self, event, output, verbose):
        out = output.evt_out

        event.getByLabel(labelPacked, handlePacked)
        event.getByLabel(labelPruned, handlePruned)
        # get the product
        genp_packed = handlePacked.product()
        genp_pruned = handlePruned.product()

        Btag = None
        Btag = None

        for p in genp_pruned:
            if p.pdgId() == self.tag_pdgId:
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

                if p.numberOfDaughters()>1: #numberOfDaughters avoids B0 from oscillation
                    Btag = p


        if Btag is None:
            return False

        out['B_pt'] = Btag.pt()
        out['B_pz'] = Btag.pz()
        out['B_p'] = Btag.p()
        out['B_eta'] = Btag.eta()

        vtx_prod = np.array([Btag.vx(), Btag.vy(), Btag.vz()])
        d = Btag.daughter(0)
        vtx_decay = np.array([d.vx(), d.vy(), d.vz()])
        d_flightB = vtx_decay - vtx_prod
        if np.sum(d_flightB) == 0 and Btag.mother(0).pdgId() == -Btag.pdgId():
            Bbar = Btag.mother(0)
            vtx_prod = np.array([Bbar.vx(), Bbar.vy(), Bbar.vz()])
            d_flightB = vtx_decay - vtx_prod
        d_flightB = rt.TVector3(d_flightB[0], d_flightB[1], d_flightB[2])
        event.d_flightB_MC = d_flightB
        dl = np.linalg.norm(vtx_prod - vtx_decay)
        out['dl'] = dl
        out['B_tau'] = dl/c_light * Btag.mass()/Btag.p() #ns
        out['B_p_dl'] = Btag.mass() * dl / (getMeanLife(521) * np.log(2) * c_light)

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

        p4_vis, N_neutral = sumVisibleP4(Btag, p4_vis, N_neutral)

        out['N_neutral'] = N_neutral
        out['M_vis_MC'] = p4_vis.M()
        out['M_miss_MC'] = Btag.mass() - p4_vis.M()

        p4_Btag = getTLorenzVector(Btag)
        out['dTheta_visB_MC'] = p4_vis.Angle(p4_Btag.Vect())

        #Find the D or D*
        for d in Btag.daughterRefVector():
            if abs(d.pdgId()) in [413, 421]:
                D = d
                break
            elif abs(d.pdgId()) == 423:
                for dd in d.daughterRefVector():
                    if abs(dd.pdgId()) == 421:
                        D = dd
                break

        #Find the muon
        for d in Btag.daughterRefVector():
            if abs(d.pdgId()) == 13:
                mu_probe = d
                break
            elif abs(d.pdgId()) == 15:
                for dd in d.daughterRefVector():
                    if abs(dd.pdgId()) == 13:
                        mu_probe = dd
                break

        addToOutput(mu_probe, 'mu', out)
        vtx_D = np.array([D.vx(), D.vy(), D.vz()])
        vtx_mu = np.array([mu_probe.vx(), mu_probe.vy(), mu_probe.vz()])

        dir_mu = np.array([mu_probe.px(), mu_probe.py(), mu_probe.pz()])
        dir_mu /= np.linalg.norm(dir_mu)

        out['mu_ip_MC'] = np.linalg.norm(np.cross(dir_mu, vtx_mu-vtx_D))


        out['D_pt_MC'] = D.pt()
        p4_D = getTLorenzVector(D)
        out['q2_MC'] = (p4_Btag - p4_D).M2()
        if d_flightB.Mag() > 0:
            out['D_pthat_MC'] = p4_D.Pt(d_flightB*(1/d_flightB.Mag()))
            pnu_ext = getNuMomentum(d_flightB, p4_vis, Btag.pt())
            out['M_ext_MC'] = (p4_vis + pnu_ext).M()
        else:
            out['D_pthat_MC'] = -999
            out['M_ext_MC'] = -999
            print '------------> Istantaneus decay'
            print vtx_decay - vtx_prod
            for i in range(Btag.numberOfMothers()):
                print Btag.mother(i).pdgId()

        #Save in the event the interesting MC particles
        event.MC_part = {}
        event.MC_part['D'] = D
        event.MC_part['mu'] = mu_probe
        event.MC_part['B'] = Btag

        if abs(D.pdgId()) == 421:
            D0 = D
        else:
            for d in D.daughterRefVector():
                if abs(d.pdgId()) == 211:
                    event.MC_part['pisoft'] = d
                    addToOutput(d, 'pisoft', out)
                elif abs(d.pdgId()) == 421:
                    event.MC_part['D0'] = d
                    D0 = d
        for d in D0.daughterRefVector():
            if abs(d.pdgId()) == 211:
                event.MC_part['pi'] = d
                addToOutput(d, 'pi', out)
            elif abs(d.pdgId()) == 321:
                event.MC_part['K'] = d
                addToOutput(d, 'K', out)


        #Print muons info:
        # print '----- MC info -----'
        # for d in Btag.daughterRefVector():
        #     if abs(d.pdgId()) ==  13:
        #         print 'Tag Mu:'
        #         getTLorenzVector(d, True)
        # print 'Probe Mu:'
        # getTLorenzVector(mu_probe, True)


        return True
