from ROOT import *
#from tdrstyle import *

gROOT.Reset()

# Uncomment the line corresponding to the MC tune we're using
mcalgo = 'tuneD6T'; mcweight = '1'
#mcalgo = 'pythia8'; mcweight = '2'
#mcalgo = 'tuneP0'; mcweight = '3'

# Same for energy
#energymc = '7'
#energydat = '7TeV'
energymc = '900'
energydat = '900'

dirr = '/nfs/data35/cms/keith/CMSSW_3_5_8/src/'

xi_data = dirr + 'Xi_' + energydat + '_data_histos.root'
xi_mc = dirr + 'efficienciesXi_' + energymc + '.root.real'

findat = TFile(xi_data, 'READ')
xipt_rec_data = numFitdata7TeVpT_copy.Clone('xipt_rec_data')
xiy_rec_data = numFitdata7TeVY_copy.Clone('xiy_rec_data')
xidataptentries = xipt_rec_data.GetSumOfWeights()
xidatayentries = xiy_rec_data.GetSumOfWeights()
finmc = TFile(xi_mc, 'READ')
xipt_rec_mc = numEffpT_copy.Clone('xipt_rec_mc')
xiy_rec_mc = numEffY_copy.Clone('xiy_rec_mc')
ximcptentries = xipt_rec_mc.GetSumOfWeights()
ximcyentries = xiy_rec_mc.GetSumOfWeights()

xipt_ScaleFactor = xidataptentries/ximcptentries
xipt_weight = xipt_rec_data.Clone('xipt_weight')
xipt_scale = xipt_rec_mc.Clone('xipt_scale')
xipt_rec_mc.Scale( xipt_ScaleFactor )
xipt_weight.Divide(xipt_rec_data, xipt_rec_mc)

xiy_ScaleFactor = xidatayentries/ximcyentries
xiy_weight = xiy_rec_data.Clone('xiy_weight')
xiy_scale = xiy_rec_mc.Clone('xiy_scale')
xiy_rec_mc.Scale( xiy_ScaleFactor )
xiy_weight.Divide(xiy_rec_data, xiy_rec_mc)

can = TCanvas('c1', 'Xi pt')
xipt_weight.Draw()
outFile = 'xipt_weight_' + energymc + '_' + mcalgo + '.C'
can.SaveAs(outFile)
outFile = 'xipt_weight_' + energymc + '_' + mcalgo + '.eps'
can.SaveAs(outFile)
outFile = 'xipt_weight_' + energymc + '_' + mcalgo + '.root'
can.SaveAs(outFile)

can = TCanvas('c1', 'Xi y')
xiy_weight.Draw()
outFile = 'xiy_weight_' + energymc + '_' + mcalgo + '.C'
can.SaveAs(outFile)
outFile = 'xiy_weight_' + energymc + '_' + mcalgo + '.eps'
can.SaveAs(outFile)
outFile = 'xiy_weight_' + energymc + '_' + mcalgo + '.root'
can.SaveAs(outFile)
