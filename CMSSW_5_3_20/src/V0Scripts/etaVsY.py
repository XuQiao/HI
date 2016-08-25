from ROOT import *

gROOT.SetBatch(1)

fin = TFile('/nfs/data35/cms/v0ntuple/mc_7TeV_pythia8_may22fixes.root', 'READ')

fin.cd('analyzeKshort')
ksntpl = gDirectory.Get('ntuple')

ksMassConst = 0.497614
ksMassConst2 = ksMassConst * ksMassConst

ksEtaVsY = TH2F( 'ksEtaVsY', 'ksEtaVsY', 200, 0., 5., 200, 0., 5. )
ksEta = TH1F( 'ksEta', 'ksEta', 200, 0., 5. )
ksRap = TH1F( 'ksRapidity', 'ksRapidity', 200, 0., 5. )

#for ndx in xrange( ksntpl.GetEntriesFast() ):
for ndx in xrange( 500000 ):
    #print ndx
    ksTree = ksntpl.LoadTree(ndx)
    if ksTree < 0:
        break
    ksEntry = ksntpl.GetEntry(ndx)
    if ksEntry <= 0:
        break

    pZ = std.vector(float)( ksntpl.v0pZ )
    pT = std.vector(float)( ksntpl.v0pT )
    eta = std.vector(float)( ksntpl.v0Eta )

    for ndx2 in range( pZ.size() ):
        pz2 = pZ[ndx2] * pZ[ndx2]
        pt2 = pT[ndx2] * pT[ndx2]
        p2 = pz2 + pt2
        energy = sqrt( p2 + ksMassConst2 )

        y = 0.5*log( (energy + pZ[ndx2]) / (energy - pZ[ndx2]) )

        #print 'Rapidity: ' + str(y) + ', Pseudorapidity: ' + str(eta[ndx2])

        ksEtaVsY.Fill( fabs(y), fabs(eta[ndx2]) )
        ksEta.Fill( fabs(eta[ndx2]) )
        ksRap.Fill( fabs(y) )


can = TCanvas()
ksEtaVsY.SetXTitle('K0S rapidity y')
ksEtaVsY.SetYTitle('K0S pseudorapidity #eta')
ksEtaVsY.Draw()
can.SaveAs('etaVsY.png')

ksRap.Draw()
can.SaveAs('ksRapidity.png')

ksEta.Draw()
can.SaveAs('ksEta.png')

