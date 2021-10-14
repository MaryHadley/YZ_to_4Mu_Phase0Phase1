import ROOT
import sys
from DataFormats.FWLite import Events, Handle
from collections import OrderedDict
from array import array
from ROOT import TLorentzVector


# Make VarParsing object
# https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideAboutPythonConfigFile#VarParsing_Example
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.parseArguments()
print options

# Events takes either
# - single file name
# - list of file names
# - VarParsing options

# use Varparsing object
events = Events (options)

# Generated stuff
handlePruned  = Handle ("std::vector<reco::GenParticle>")
labelPruned=("genParticles")

# Create histograms, etc.
ROOT.gROOT.SetBatch()        # don't pop up canvases
ROOT.gROOT.SetStyle('Plain') # white background

# ntupla
branches = [
# Generated
        #    'g_mu_pt',
         #    'g_mu_eta',
             'numGenMu',
              'goodMuCount',
]


file_out = ROOT.TFile('nTuples_muonCheck.root', 'recreate')
file_out.cd()

ntuple   = ROOT.TNtuple('tree', 'tree', ':'.join(branches))
nTot=0
for event in events:
    goodMuCount=0
    print 'Processing event: %i...'%(nTot)

    tofill   = OrderedDict(zip(branches, [-99.]*len(branches)))
    # Generated stuff
    event.getByLabel (labelPruned, handlePruned)
    pruned_gen_particles = handlePruned.product()

    gen_mu = [pp for pp in pruned_gen_particles if abs(pp.pdgId()) == 13]
    g_mu_size = len(gen_mu)
    for gp in gen_mu:

        g_mu_lv = TLorentzVector(gp.px(), gp.py(), gp.pz(), gp.energy())
        #tofill['g_mu_pt'] = g_mu_lv.Pt()
        #tofill['g_mu_eta'] = g_mu_lv.Eta()
        #if g_mu_lv.Pt() < 2.5:
        #   print "YIKES"
        if g_mu_lv.Pt() >= 2.5 and abs(g_mu_lv.Eta()) <= 2.5:
          goodMuCount += 1
    tofill['numGenMu'] = g_mu_size
    tofill['goodMuCount'] = goodMuCount
    ntuple.Fill(array('f',tofill.values()))

    nTot += 1

# make a canvas, draw, and save it
file_out.cd()
ntuple.Write()

file_out.Close()
