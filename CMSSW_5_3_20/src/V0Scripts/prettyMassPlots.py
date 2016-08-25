from ROOT import *
import math, time
import tdrstyle

tdrStyle = tdrstyle.setTDRStyle()
tdrStyle.SetOptStat(0)
tdrStyle.SetPadTopMargin(0.059)

TGaxis.SetMaxDigits(3)
gROOT.SetBatch(0)

isMC = True

dir_ = 'finalMassPlots/'

# Data
#ksfin_7000 = TFile('massHistos/data_dNdX_kshort_7000.root', 'READ')
#ksfin_900  = TFile('massHistos/data_dNdX_kshort_900.root',  'READ')
#lamfin_7000 = TFile('massHistos/data_dNdX_lambda.root', 'READ')
#lamfin_900  = TFile('massHistos/data_dNdX_lambda_900.root', 'READ')

# MC
ksfin_7000 = TFile('mcmass/sim_dNdX_kshort_pythia8.root', 'READ')
ksfin_900 = TFile('mcmass/sim_dNdX_kshort_tuneD6T.root', 'READ')
lamfin_7000 = TFile('mcmass/sim_dNdX_lambda_pythia8.root', 'READ')
lamfin_900 = TFile('mcmass/sim_dNdX_lambda_tuneD6T.root', 'READ')

ksmass7000 = ksfin_7000.Get('vee_withcuts')
lammass7000 = lamfin_7000.Get('vee_withcuts')
ksmass900 = ksfin_900.Get('vee_withcuts')
lammass900 = lamfin_900.Get('vee_withcuts')

# Coordinates for fit parameter TPave
pavex1 = 0.42
pavey1 = 0.35
pavex2 = 0.9
pavey2 = 0.7

## TPaves for labels
labltxt1_1 = TText()
labltxt1_2 = TText()
labltxt2_1 = TText()
labltxt2_2 = TText()

labl = TPaveText(0.55, 0.81, 0.85, 0.91, "brNDC")
labl.SetBorderSize(0)
labl.SetFillColor(0)
labl.SetFillStyle(0)
labl.SetTextSize(0.045)
labl2 = TPaveText(0.55, 0.81, 0.85, 0.91, "brNDC")
labl2.SetBorderSize(0)
labl2.SetFillColor(0)
labl2.SetFillStyle(0)
labl2.SetTextSize(0.045)
#labld6t = TPaveText(0.65, 0.81, 0.95, 0.91, "brNDC")
#labld6t.SetBorderSize(0)
#labld6t.SetFillColor(0)
#labld6t.SetFillStyle(0)
#labld6t.SetTextSize(0.045)

if isMC :
    labltxt1_1 = labl.AddText("CMS Simulation Pythia 8")
    labltxt1_2 = labl.AddText("#sqrt{s} = 7 TeV")
    labltxt2_1 = labl2.AddText("CMS Simulation Tune D6T")
    labltxt2_2 = labl2.AddText("#sqrt{s} = 900 GeV")
    #labltxt3_1 = labld6t.AddText("CMS Simulation tune D6T")
    #labltxt3_2 = labld6t.AddText("#sqrt{s} = 900GeV")
else :
    labltxt1_1 = labl.AddText("CMS Preliminary")
    labltxt1_2 = labl.AddText("#sqrt{s} = 7 TeV")
    labltxt2_1 = labl2.AddText("CMS Preliminary")
    labltxt2_2 = labl2.AddText("#sqrt{s} = 900 GeV")

    
ksMassXmin = 0.340
ksMassXmax = 0.740
ksMassNbins = 400
ksMassBinWidth = (ksMassXmax - ksMassXmin) / ksMassNbins
ksMassBinWidth_DG = ksMassBinWidth / sqrt(2*math.pi)

lamMassXmin = 1.09
lamMassXmax = 1.19
lamMassNbins = 200
lamMassBinWidth = (lamMassXmax - lamMassXmin) / lamMassNbins
lamMassBinWidth_DG = lamMassBinWidth / sqrt(2*math.pi)

tdrStyle.SetOptFit(0)


can = TCanvas()

ksmass7000.SetXTitle('#pi^{+}#pi^{-} invariant mass [GeV/c^{2}]')
ksmass7000.SetYTitle('Candidates / ' + str('%.3f' % ksMassBinWidth)
                     + ' GeV/c^{2}')
ksmass7000.GetXaxis().SetRangeUser(0.41, ksMassXmax)
ksmass7000.GetYaxis().SetTitleOffset(1.412)
#TGaxis.SetMaxDigits(1)

ksfit = ksmass7000.GetFunction('ksFit_DG_exp')
ksfit.SetNpx(50000)
ks7yield = ksfit.GetParameter(0)
ks7yielderr = ksfit.GetParError(0)
ks7mass = ksfit.GetParameter(1)*1000.
ks7masserr = ksfit.GetParError(1)*1000.
ks7sigma1 = ksfit.GetParameter(2)*1000.
ks7sigma1err = ksfit.GetParError(2)*1000.
ks7frac = ksfit.GetParameter(3)
ks7fracerr = ksfit.GetParError(3)
ks7sigma2 = ksfit.GetParameter(4)*1000.
ks7sigma2err = ksfit.GetParError(4)*1000.
ks7yieldstr = str('%.0f' % ks7yield)
ks7yielderrstr = str('%.0f' % ks7yielderr)
ks7massstr = str('%.3f' % ks7mass)
ks7masserrstr = str('%.3f' % ks7masserr)
ks7sigma1str = str('%.2f' % ks7sigma1)
ks7sigma1errstr = str('%.2f' % ks7sigma1err)
ks7sigma2str = str('%.2f' % ks7sigma2)
ks7sigma2errstr = str('%.2f' % ks7sigma2err)
ks7fracstr = str('%.3f' % ks7frac)
ks7fracerrstr = str('%.3f' % ks7fracerr)

print str(ks7mass) + '+-' + str(ks7masserr)

ks7pave = TPaveText(pavex1, pavey1, pavex2, pavey2, 'brNDC')
ks7pave.SetBorderSize(0)
ks7pave.SetFillStyle(0)
ks7pave.SetFillColor(0)
ks7pave.SetTextSize(0.04)
ks7pave.SetTextAlign(12)
ks7pave_t1 = ks7pave.AddText('Yield = ' + ks7yieldstr + ' #pm ' + ks7yielderrstr)
ks7pave_t15 = ks7pave.AddText('Mass = ' + ks7massstr + ' #pm ' + ks7masserrstr
                              + ' MeV/c^{2}')
ks7pave_t2 = ks7pave.AddText('#sigma_{1} = ' + ks7sigma1str
                               + ' #pm ' + ks7sigma1errstr + ' MeV/c^{2}')
ks7pave_t3 = ks7pave.AddText('Fraction with #sigma_{1}: ' + ks7fracstr)
                               #+ ' #pm ' + ks7fracerrstr)
ks7pave_t4 = ks7pave.AddText('#sigma_{2} = ' + ks7sigma2str
                               + ' #pm ' + ks7sigma2errstr + ' MeV/c^{2}')

## print 'K0S fit parameters:'
## print ks7yieldstr + '+-' + ks7yielderrstr
## print ks7sigma1str + '+-' + ks7sigma1errstr
## print ks7fracstr + '+- ' + ks7fracerrstr
## print ks7sigma2str + '+-' + ks7sigma2errstr
## print ''

ksmass7000.Draw('hist')
ksfit.Draw('same')
ks7pave.Draw('same')
labl.Draw('same')
if isMC:
    can.SaveAs(dir_ + "ksMass_7000_MC.png")
    can.SaveAs(dir_ + "ksMass_7000_MC.eps")
    can.SaveAs(dir_ + "ksMass_7000_MC.root")
    can.SaveAs(dir_ + "ksMass_7000_MC.pdf")
else:
    can.SaveAs(dir_ + "ksMass_7000.png")
    can.SaveAs(dir_ + "ksMass_7000.eps")
    can.SaveAs(dir_ + "ksMass_7000.root")
    can.SaveAs(dir_ + "ksMass_7000.pdf")

# Sleep for a minute so we can look at things
#time.sleep(60)

lammass7000.SetXTitle('p #pi^{-} (+ c.c.) invariant mass [GeV/c^{2}]')
lammass7000.SetYTitle('Candidates / ' + str('%.4f' % lamMassBinWidth)
                      + ' GeV/c^{2}')
lammass7000.GetXaxis().SetNdivisions(506)
#print "lam offset: " + str(lammass7000.GetYaxis().GetTitleOffset())
lammass7000.GetYaxis().SetTitleOffset(1.312)

lamfit = lammass7000.GetFunction('lamFit_DG_highbg')
lamfit.SetNpx(50000)
lam7mass = lamfit.GetParameter(1)*1000.
lam7masserr = lamfit.GetParError(1)*1000.
lam7sigma1 = lamfit.GetParameter(2)*1000.
lam7sigma1err = lamfit.GetParError(2)*1000.
lam7frac = lamfit.GetParameter(3)
lam7fracerr = lamfit.GetParError(3)
lam7sigma2 = lamfit.GetParameter(4)*1000.
lam7sigma2err = lamfit.GetParError(4)*1000.

lam7massstr = str('%.3f' % lam7mass)
lam7masserrstr = str('%.3f' % lam7masserr)
lam7sigma1str = str('%.2f' % lam7sigma1)
lam7sigma1errstr = str('%.2f' % lam7sigma1err)
lam7sigma2str = str('%.2f' % lam7sigma2)
lam7sigma2errstr = str('%.2f'% lam7sigma2err)
lam7fracstr = str('%.3f' % lam7frac)
lam7fracerrstr = str('%.3f' % lam7fracerr)
#print lam7sigma1str + ', ' + lam7fracstr + ', ' + lam7sigma2str
lam7yieldstr = str('%.0f' % lamfit.GetParameter(0))
lam7yielderrstr = str('%.0f' % lamfit.GetParError(0))

# Create TPave with fit quantities
lam7pave = TPaveText(pavex1, pavey1, pavex2, pavey2, 'brNDC')
lam7pave.SetBorderSize(0)
lam7pave.SetFillStyle(0)
lam7pave.SetFillColor(0)
lam7pave.SetTextSize(0.04)
lam7pave.SetTextAlign(12)
lam7pave_t1 = lam7pave.AddText('Yield = ' + lam7yieldstr + ' #pm ' + lam7yielderrstr)
lam7pave_t15 = lam7pave.AddText('Mass = ' + lam7massstr + ' #pm ' + lam7masserrstr
                                + ' MeV/c^{2}')
lam7pave_t2 = lam7pave.AddText('#sigma_{1} = ' + lam7sigma1str
                               + ' #pm ' + lam7sigma1errstr + ' MeV/c^{2}')
lam7pave_t3 = lam7pave.AddText('Fraction with #sigma_{1}: ' + lam7fracstr)
                               #+ ' #pm ' + lam7fracerrstr)
lam7pave_t4 = lam7pave.AddText('#sigma_{2} = ' + lam7sigma2str
                               + ' #pm ' + lam7sigma2errstr + ' MeV/c^{2}')

lammass7000.Draw('hist')
lamfit.Draw('same')
lam7pave.Draw('same')
labl.Draw('same')
if isMC:
    can.SaveAs(dir_ + "lambdaMass_7000MC.png")
    can.SaveAs(dir_ + "lambdaMass_7000MC.eps")
    can.SaveAs(dir_ + "lambdaMass_7000MC.root")
    can.SaveAs(dir_ + "lambdaMass_7000MC.pdf")
else:
    can.SaveAs(dir_ + "lambdaMass_7000.png")
    can.SaveAs(dir_ + "lambdaMass_7000.eps")
    can.SaveAs(dir_ + "lambdaMass_7000.root")
    can.SaveAs(dir_ + "lambdaMass_7000.pdf")

#time.sleep(100)

# Ks mass, 900 GeV
ksmass900.SetXTitle('#pi^{+}#pi^{-} invariant mass [GeV/c^{2}]')
ksmass900.SetYTitle('Candidates / ' + str('%.3f' % ksMassBinWidth)
                    + ' GeV/c^{2}')
ksmass900.GetXaxis().SetRangeUser(0.41, ksMassXmax)
ksmass900.GetYaxis().SetTitleOffset(1.312)
#TGaxis.SetMaxDigits(1)

#ksmass900.Draw()
#time.sleep(60)

ksfit900 = ksmass900.GetFunction('ksFit_DG_quad')
ksfit900.SetNpx(50000)
ks900yield = ksfit900.GetParameter(0)
ks900yielderr = ksfit900.GetParError(0)
ks900mass = ksfit900.GetParameter(1)*1000.
ks900masserr = ksfit900.GetParError(1)*1000.
ks900sigma1 = ksfit900.GetParameter(2)*1000.
ks900sigma1err = ksfit900.GetParError(2)*1000.
ks900frac = ksfit900.GetParameter(3)
ks900fracerr = ksfit900.GetParError(3)
ks900sigma2 = ksfit900.GetParameter(4)*1000.
ks900sigma2err = ksfit900.GetParError(4)*1000.
ks900yieldstr = str('%.0f' % ks900yield)
ks900yielderrstr = str('%.0f' % ks900yielderr)
ks900massstr = str('%.3f' % ks900mass)
ks900masserrstr = str('%.3f' % ks900masserr)
ks900sigma1str = str('%.2f' % ks900sigma1)
ks900sigma1errstr = str('%.2f' % ks900sigma1err)
ks900sigma2str = str('%.2f' % ks900sigma2)
ks900sigma2errstr = str('%.2f' % ks900sigma2err)
ks900fracstr = str('%.3f' % ks900frac)
ks900fracerrstr = str('%.3f' % ks900fracerr)

ks900pave = TPaveText(pavex1, pavey1, pavex2, pavey2, 'brNDC')
ks900pave.SetBorderSize(0)
ks900pave.SetFillStyle(0)
ks900pave.SetFillColor(0)
ks900pave.SetTextSize(0.04)
ks900pave.SetTextAlign(12)
ks900pave_t1 = ks900pave.AddText('Yield = ' + ks900yieldstr + ' #pm ' + ks900yielderrstr)
ks900pave_t15 = ks900pave.AddText('Mass = ' + ks900massstr + ' #pm ' + ks900masserrstr
                                  + ' MeV/c^{2}')
ks900pave_t2 = ks900pave.AddText('#sigma_{1} = ' + ks900sigma1str
                               + ' #pm ' + ks900sigma1errstr + ' MeV/c^{2}')
ks900pave_t3 = ks900pave.AddText('Fraction with #sigma_{1}: ' + ks900fracstr)
ks900pave_t4 = ks900pave.AddText('#sigma_{2} = ' + ks900sigma2str
                               + ' #pm ' + ks900sigma2errstr + ' MeV/c^{2}')

## print 'K0S fit parameters:'
## print ks900yieldstr + '+-' + ks900yielderrstr
## print ks900sigma1str + '+-' + ks900sigma1errstr
## print ks900fracstr + '+- ' + ks900fracerrstr
## print ks900sigma2str + '+-' + ks900sigma2errstr
## print ''

ksmass900.Draw('hist')
ksfit900.Draw('same')
ks900pave.Draw('same')
labl2.Draw('same')
if isMC:
    can.SaveAs(dir_ + "ksMass_900MC.png")
    can.SaveAs(dir_ + "ksMass_900MC.eps")
    can.SaveAs(dir_ + "ksMass_900MC.root")
    can.SaveAs(dir_ + "ksMass_900MC.pdf")
else:
    can.SaveAs(dir_ + "ksMass_900.png")
    can.SaveAs(dir_ + "ksMass_900.eps")
    can.SaveAs(dir_ + "ksMass_900.root")
    can.SaveAs(dir_ + "ksMass_900.pdf")

# Lam mass, 900 GeV
lammass900.SetXTitle('p #pi^{-} (+ c.c.) invariant mass [GeV/c^{2}]')
lammass900.SetYTitle('Candidates / ' + str('%.4f' % lamMassBinWidth)
                      + ' GeV/c^{2}')
lammass900.GetXaxis().SetNdivisions(506)
#print "lam offset: " + str(lammass900.GetYaxis().GetTitleOffset())
lammass900.GetYaxis().SetTitleOffset(1.312)

#lammass900.Draw()
#time.sleep(60)

lamfit900 = lammass900.GetFunction('lamFit_DG_highbg')
lamfit900.SetNpx(50000)
lam900mass = lamfit900.GetParameter(1)*1000.
lam900masserr = lamfit900.GetParError(1)*1000.
lam900sigma1 = lamfit900.GetParameter(2)*1000.
lam900sigma1err = lamfit900.GetParError(2)*1000.
lam900frac = lamfit900.GetParameter(3)
lam900fracerr = lamfit900.GetParError(3)
lam900sigma2 = lamfit900.GetParameter(4)*1000.
lam900sigma2err = lamfit900.GetParError(4)*1000.

lam900massstr = str('%.3f' % lam900mass)
lam900masserrstr = str('%.3f' % lam900masserr)
lam900sigma1str = str('%.2f' % lam900sigma1)
lam900sigma1errstr = str('%.2f' % lam900sigma1err)
lam900sigma2str = str('%.2f' % lam900sigma2)
lam900sigma2errstr = str('%.2f'% lam900sigma2err)
lam900fracstr = str('%.3f' % lam900frac)
lam900fracerrstr = str('%.3f' % lam900fracerr)
#print lam900sigma1str + ', ' + lam900fracstr + ', ' + lam900sigma2str
lam900yieldstr = str('%.0f' % lamfit900.GetParameter(0))
lam900yielderrstr = str('%.0f' % lamfit900.GetParError(0))

# Create TPave with fit quantities
lam900pave = TPaveText(pavex1, pavey1, pavex2, pavey2, 'brNDC')
lam900pave.SetBorderSize(0)
lam900pave.SetFillStyle(0)
lam900pave.SetFillColor(0)
lam900pave.SetTextSize(0.04)
lam900pave.SetTextAlign(12)
lam900pave_t1 = lam900pave.AddText('Yield = ' + lam900yieldstr + ' #pm '
                                   + lam900yielderrstr)
lam900pave_t15 = lam900pave.AddText('Mass = ' + lam900massstr + ' #pm '
                                    + lam900masserrstr + ' MeV/c^{2}')
lam900pave_t2 = lam900pave.AddText('#sigma_{1} = ' + lam900sigma1str
                               + ' #pm ' + lam900sigma1errstr + ' MeV/c^{2}')
lam900pave_t3 = lam900pave.AddText('Fraction with #sigma_{1}: ' + lam900fracstr)
                               #+ ' #pm ' + lam900fracerrstr)
lam900pave_t4 = lam900pave.AddText('#sigma_{2} = ' + lam900sigma2str
                               + ' #pm ' + lam900sigma2errstr + ' MeV/c^{2}')

lammass900.Draw('hist')
lamfit900.Draw('same')
lam900pave.Draw('same')
labl2.Draw('same')
if isMC:
    can.SaveAs(dir_ + "lambdaMass_900MC.png")
    can.SaveAs(dir_ + "lambdaMass_900MC.eps")
    can.SaveAs(dir_ + "lambdaMass_900MC.root")
    can.SaveAs(dir_ + "lambdaMass_900MC.pdf")
else:
    can.SaveAs(dir_ + "lambdaMass_900.png")
    can.SaveAs(dir_ + "lambdaMass_900.eps")
    can.SaveAs(dir_ + "lambdaMass_900.root")
    can.SaveAs(dir_ + "lambdaMass_900.pdf")


print 'Done.'
#time.sleep(60)
