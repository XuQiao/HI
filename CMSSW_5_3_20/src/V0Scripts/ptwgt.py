from ROOT import *
#from tdrstyle import *

gROOT.Reset()
## tdr = setTDRstyle()
## gStyle.SetStatFormat("5.5g");
## gStyle.SetStatY(0.94);
## gStyle.SetOptFit(111);
## gStyle.SetFitFormat("5.5g");
## gStyle.SetStatFontSize(0.03);
## gStyle.SetTitleSize(0.05, "XYZ");
## gStyle.SetLabelSize(0.04, "XYZ");
## gStyle.SetTitleYOffset(1.4);
## gStyle.SetStatH(0.16);
## gStyle.SetStatW(0.3);
## gStyle.SetTextSize(0.038);
## gStyle.SetTextFont(42);
## gStyle.SetOptTitle(0);
## gStyle.SetOptStat(0);


## ksptarray = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.2,3.4,3.6,3.8,4.0,4.4,4.8,5.2,5.6,6.0,7.0,8.0,10.0]
## ksptNbins = len(ksptarray) - 1
## ksptBinMinima = std.vector(double)()
## for pt in ksptarray:
##     ksptBinMinima.push_back( pt )
## ksptBinMinima.pop_back()
## kspt_default_binWidth = 0.1

## lamptarray = [0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.2,3.4,3.6,3.8,4.0,4.4,4.8,5.2,5.6,6.0,7.0,8.0,10.0]
## lamptNbins = len(lamptarray) - 1
## lamptBinMinima = std.vector(double)()
## for pt in lamptarray:
##     lamptBinMinima.push_back( pt )
## lamptBinMinima.pop_back()
## lampt_default_binWidth = 0.1

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

ksdirr_mc = 'ks_' + mcalgo + '_' + energy2 + '/plots/'
ksdirr_data = 'ks_data_' + energy2 + '/plots/'
lamdirr_mc = 'lam_' + mcalgo + '_' + energy2 + '/plots/'
lamdirr_data = 'lam_data_' + energy2 + '/plots/'

kspt_data = ksdirr_data + 'ks_dNdpT' + '.root'
kspt_mc = ksdirr_mc + 'ks_dNdpT_rec_' + mcalgo + '_' + energy + '_weight' + mcweight + '.root'
lampt_data = lamdirr_data + 'lam_dNdpT' + '.root'
lampt_mc = lamdirr_mc + 'lam_dNdpT_rec_' + mcalgo + '_' + energy + '_weight' + mcweight + '.root'
lamRegpt_data = lamdirr_data + 'lamReg_dNdpT' + '.root'
lamRegpt_mc = lamdirr_mc + 'lamReg_dNdpT_rec_' + mcalgo + '_' + energy + '_weight' + mcweight + '.root'
lamBarpt_data = lamdirr_data + 'lamBar_dNdpT' + '.root'
lamBarpt_mc = lamdirr_mc + 'lamBar_dNdpT_rec_' + mcalgo + '_' + energy + '_weight' + mcweight + '.root'

fin = TFile(kspt_data, 'READ')
c1 = fin.Get('c1')
kspt_rec_data = c1.GetPrimitive('ks_dNdpT_scale')
ksdatarecentries = kspt_rec_data.GetSumOfWeights()
fin.Close()

fin = TFile(kspt_mc, 'READ')
c1 = fin.Get('c1')
kspt_rec_mc = c1.GetPrimitive('ks_dNdpT_rec_scale')
ksmcrecentries = kspt_rec_mc.GetSumOfWeights()
fin.Close()

fin = TFile(lampt_data, 'READ')
c1 = fin.Get('c1')
lampt_rec_data = c1.GetPrimitive('lam_dNdpT_scale')
lamdatarecentries = lampt_rec_data.GetSumOfWeights()
fin.Close()

fin = TFile(lampt_mc, 'READ')
c1 = fin.Get('c1')
lampt_rec_mc = c1.GetPrimitive('lam_dNdpT_rec_scale')
lammcrecentries = lampt_rec_mc.GetSumOfWeights()
fin.Close()

fin = TFile(lamRegpt_data, 'READ')
c1 = fin.Get('c1')
lamRegpt_rec_data = c1.GetPrimitive('lamReg_dNdpT_scale')
lamRegdatarecentries = lamRegpt_rec_data.GetSumOfWeights()
fin.Close()

fin = TFile(lamRegpt_mc, 'READ')
c1 = fin.Get('c1')
lamRegpt_rec_mc = c1.GetPrimitive('lamReg_dNdpT_rec_scale')
lamRegmcrecentries = lamRegpt_rec_mc.GetSumOfWeights()
fin.Close()

fin = TFile(lamBarpt_data, 'READ')
c1 = fin.Get('c1')
lamBarpt_rec_data = c1.GetPrimitive('lamBar_dNdpT_scale')
lamBardatarecentries = lamBarpt_rec_data.GetSumOfWeights()
fin.Close()

fin = TFile(lamBarpt_mc, 'READ')
c1 = fin.Get('c1')
lamBarpt_rec_mc = c1.GetPrimitive('lamBar_dNdpT_rec_scale')
lamBarmcrecentries = lamBarpt_rec_mc.GetSumOfWeights()
fin.Close()

ksScaleFactor = ksdatarecentries/ksmcrecentries
lamScaleFactor = lamdatarecentries/lammcrecentries
lamRegScaleFactor = lamRegdatarecentries/lamRegmcrecentries
lamBarScaleFactor = lamBardatarecentries/lamBarmcrecentries

kspt_weight = kspt_rec_data.Clone('kspt_weight')
#kspt_scale = kspt_rec_mc.Clone('kspt_scale')
lampt_weight = lampt_rec_data.Clone('lampt_weight')
#lampt_scale = lampt_rec_mc.Clone('lampt_scale')
lamRegpt_weight = lamRegpt_rec_data.Clone('lamRegpt_weight')
lamBarpt_weight = lamBarpt_rec_data.Clone('lamBarpt_weight')

kspt_rec_mc.Scale( ksScaleFactor )
lampt_rec_mc.Scale( lamScaleFactor )
lamRegpt_rec_mc.Scale( lamRegScaleFactor )
lamBarpt_rec_mc.Scale( lamBarScaleFactor )

kspt_weight.Divide(kspt_rec_data, kspt_rec_mc)
lampt_weight.Divide(lampt_rec_data, lampt_rec_mc)
lamRegpt_weight.Divide(lamRegpt_rec_data, lamRegpt_rec_mc)
lamBarpt_weight.Divide(lamBarpt_rec_data, lamBarpt_rec_mc)

can = TCanvas('c1', 'ks pt')
kspt_weight.Draw()
outFile = 'kspt_weight_' + energy + '_' + mcalgo + '.C'
can.SaveAs(outFile)
outFile = 'kspt_weight_' + energy + '_' + mcalgo + '.eps'
can.SaveAs(outFile)
outFile = 'kspt_weight_' + energy + '_' + mcalgo + '.root'
can.SaveAs(outFile)

can = TCanvas('c1', 'lam pt')
lampt_weight.Draw()
outFile = 'lampt_weight_' + energy + '_' + mcalgo + '.C'
can.SaveAs(outFile)
outFile = 'lampt_weight_' + energy + '_' + mcalgo + '.eps'
can.SaveAs(outFile)
outFile = 'lampt_weight_' + energy + '_' + mcalgo + '.root'
can.SaveAs(outFile)

can = TCanvas('c1', 'lam (part) pt')
lamRegpt_weight.Draw()
outFile = 'lamRegpt_weight_' + energy + '_' + mcalgo
can.SaveAs(outFile + '.C')
can.SaveAs(outFile + '.eps')
can.SaveAs(outFile + '.root')

can = TCanvas('c1', 'lam (anti) pt')
lamBarpt_weight.Draw()
outFile = 'lamBarpt_weight_' + energy + '_' + mcalgo
can.SaveAs(outFile + '.C')
can.SaveAs(outFile + '.eps')
can.SaveAs(outFile + '.root')
