from ROOT import *
#from tdrstyle import *

gROOT.Reset()

## ksyXmin = 0.
## ksyXmax = 2.3
## ksyNbins = 23
## ksyBinWidth = (ksyXmax - ksyXmin) / ksyNbins

## lamyXmin = 0.
## lamyXmax = 2.3
## lamyNbins = 23
## lamyBinWidth = (lamyXmax - lamyXmin) / lamyNbins


# Uncomment the line corresponding to the MC tune we're using
#mcalgo = 'tuneD6T'; mcweight = '11' ## 900 GeV
#mcalgo = 'tuneD6T'; mcweight = '1'  ## 7 TeV
mcalgo = 'pythia8'; mcweight = '2'
#mcalgo = 'tuneP0'; mcweight = '3'

# Same for energy
energy = '7TeV'
energy2 = '7000'
#energy = '900GeV'
#energy2 = '900'

ksdirr_data = 'ks_data_' + energy2 + '/plots/'
ksdirr_mc = 'ks_' + mcalgo + '_' + energy2 + '/plots/'
lamdirr_data = 'lam_data_' + energy2 + '/plots/'
lamdirr_mc = 'lam_' + mcalgo + '_' + energy2 + '/plots/'

ksy_data = ksdirr_data + 'ks_dNdy' + '.root'
ksy_mc = ksdirr_mc + 'ks_dNdy_rec_' + mcalgo + '_' + energy + '_weight' + mcweight + '.root'
lamy_data = lamdirr_data + 'lam_dNdy' + '.root'
lamy_mc = lamdirr_mc + 'lam_dNdy_rec_' + mcalgo + '_' + energy + '_weight' + mcweight + '.root'
lamRegy_data = lamdirr_data + 'lamReg_dNdy' + '.root'
lamRegy_mc = lamdirr_mc + 'lamReg_dNdy_rec_' + mcalgo + '_' + energy + '_weight' + mcweight + '.root'
lamBary_data = lamdirr_data + 'lamBar_dNdy' + '.root'
lamBary_mc = lamdirr_mc + 'lamBar_dNdy_rec_' + mcalgo + '_' + energy + '_weight' + mcweight + '.root'

fin = TFile(ksy_data, 'READ')
c1 = fin.Get('c1')
ksy_rec_data = c1.GetPrimitive('ks_dNdy')
ksdataRecEntries = ksy_rec_data.GetSumOfWeights()
fin.Close()

fin = TFile(ksy_mc, 'READ')
c1 = fin.Get('c1')
ksy_rec_mc = c1.GetPrimitive('ks_dNdy_rec')
ksmcRecEntries = ksy_rec_mc.GetSumOfWeights()
fin.Close()

fin = TFile(lamy_data, 'READ')
c1 = fin.Get('c1')
lamy_rec_data = c1.GetPrimitive('lam_dNdy')
lamdataRecEntries = lamy_rec_data.GetSumOfWeights()
fin.Close()

fin = TFile(lamy_mc, 'READ')
c1 = fin.Get('c1')
lamy_rec_mc = c1.GetPrimitive('lam_dNdy_rec')
lammcRecEntries = lamy_rec_mc.GetSumOfWeights()
fin.Close()

fin = TFile(lamRegy_data, 'READ')
c1 = fin.Get('c1')
lamRegy_rec_data = c1.GetPrimitive('lamReg_dNdy')
lamRegdataRecEntries = lamRegy_rec_data.GetSumOfWeights()
fin.Close()

fin = TFile(lamRegy_mc, 'READ')
c1 = fin.Get('c1')
lamRegy_rec_mc = c1.GetPrimitive('lamReg_dNdy_rec')
lamRegmcRecEntries = lamRegy_rec_mc.GetSumOfWeights()
fin.Close()

fin = TFile(lamBary_data, 'READ')
c1 = fin.Get('c1')
lamBary_rec_data = c1.GetPrimitive('lamBar_dNdy')
lamBardataRecEntries = lamBary_rec_data.GetSumOfWeights()
fin.Close()

fin = TFile(lamBary_mc, 'READ')
c1 = fin.Get('c1')
lamBary_rec_mc = c1.GetPrimitive('lamBar_dNdy_rec')
lamBarmcRecEntries = lamBary_rec_mc.GetSumOfWeights()
fin.Close()

ksScaleFactor = ksdataRecEntries / ksmcRecEntries
lamScaleFactor = lamdataRecEntries / lammcRecEntries
lamRegScaleFactor = lamRegdataRecEntries / lamRegmcRecEntries
lamBarScaleFactor = lamBardataRecEntries / lamBarmcRecEntries

ksy_weight = ksy_rec_data.Clone('ksy_weight')
lamy_weight = lamy_rec_data.Clone('lamy_weight')
lamRegy_weight = lamRegy_rec_data.Clone('lamRegy_weight')
lamBary_weight = lamBary_rec_data.Clone('lamBary_weight')

ksy_rec_mc.Scale( ksScaleFactor )
lamy_rec_mc.Scale( lamScaleFactor )
lamRegy_rec_mc.Scale( lamRegScaleFactor )
lamBary_rec_mc.Scale( lamBarScaleFactor )

ksy_weight.Divide( ksy_rec_data, ksy_rec_mc )
lamy_weight.Divide( lamy_rec_data, lamy_rec_mc )
lamRegy_weight.Divide( lamRegy_rec_data, lamBary_rec_mc )
lamBary_weight.Divide( lamBary_rec_data, lamBary_rec_mc )

can = TCanvas('c1', 'ks rapidity')
ksy_weight.Draw()
outFile = 'ksy_weight_' + energy + '_' + mcalgo + '.C'
can.SaveAs(outFile)
outFile = 'ksy_weight_' + energy + '_' + mcalgo + '.eps'
can.SaveAs(outFile)
outFile = 'ksy_weight_' + energy + '_' + mcalgo + '.root'
can.SaveAs(outFile)

can = TCanvas('c1', 'lam rapidity')
lamy_weight.Draw()
outFile = 'lamy_weight_' + energy + '_' + mcalgo + '.C'
can.SaveAs(outFile)
outFile = 'lamy_weight_' + energy + '_' + mcalgo + '.eps'
can.SaveAs(outFile)
outFile = 'lamy_weight_' + energy + '_' + mcalgo + '.root'
can.SaveAs(outFile)

can = TCanvas('c1', 'lam (part) rapidity')
lamRegy_weight.Draw()
outFile = 'lamRegy_weight_' + energy + '_' + mcalgo
can.SaveAs(outFile + '.C')
can.SaveAs(outFile + '.eps')
can.SaveAs(outFile + '.root')

can = TCanvas('c1', 'lam (anti) rapidity')
lamBary_weight.Draw()
outFile = 'lamBary_weight_' + energy + '_' + mcalgo
can.SaveAs(outFile + '.C')
can.SaveAs(outFile + '.eps')
can.SaveAs(outFile + '.root')
