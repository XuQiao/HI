#define seleXiData_cxx
#include "seleXiData.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <math.h>
#include <fstream>
#include <TAxis.h>
#include <TString.h>
#include <TPaveText.h>
#include <TText.h>
#include <TLatex.h>
#include <sys/stat.h> 
#include <sys/types.h>

using namespace std;

void seleXiData::Loop(){
   
  cout << " Cuts Used:\t " << "Xi3dIp " << Xi3dIpCut << "\t batPi3dIp " << batPi3dIpCut << "\t v0Pi3dIp" << v0Pi3dIpCut 
                        << "\t v0Proton3dIp "  << v0Proton3dIpCut << "\t XiSep " << XiSepCut << "\t VtxCL " << vtxCLCut << "\n";
 
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  if(nEvents < 0) nentries = fChain->GetEntriesFast();
  else
    nentries = nEvents;

  Long64_t nbytes = 0, nb = 0;
   
  // set pt bins and book histos
  SetBins();
  bookHistos();
   
  float mXi = 1.32171;
  float pi = TMath::Pi();
  float c = 0.0299793; // cm/ps
  int nCands = 0;
  int nSingleCands = 0;
  int nTotalCands = 0;
  int nEvtWCands = 0;
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);	     
    if (ientry < 0) break;			     
    nb = fChain->GetEntry(jentry);   nbytes += nb;  
      // if (Cut(ientry) < 0) continue;
    nCands = 0;
    int bestCand = -999;
    float current_batPiPt = 0.0;
    
    for ( unsigned int i = 0; i<nXI; i++ ) {

      //check NaNs first ...
      if (XiLongitudinalPCAPrimaryError->at(i) != XiLongitudinalPCAPrimaryError->at(i)) continue;
      if (XiTransversePCAPrimaryError->at(i) != XiTransversePCAPrimaryError->at(i)) continue;
      if (xiDecayVtxXE->at(i) != xiDecayVtxXE->at(i) )continue;
      if (xiDecayVtxYE->at(i) != xiDecayVtxYE->at(i) )continue;
      if (xiDecayVtxZE->at(i) != xiDecayVtxZE->at(i) )continue;
       
      //apply cuts
      if ( !trigHF || !trigTech34 || !trigTech36_39 || !trigScraping ) continue;
      if ( (*xiMass)[i] < 1.255 || (*xiMass)[i] > 1.455 ) continue;
      if ( fabs((*VMass)[i]-1.116)>0.008 ) continue;
      if ( (*xiVtxCL)[i] < vtxCLCut) continue;
      if ( (*VTrkPiQ)[i]!=(*batPiQ)[i]) continue;
       
      if(batPi3DIpSig->at(i) < batPi3dIpCut)continue;
      if(VTrkPi3DIpSig->at(i) < v0Pi3dIpCut)continue;
      if(VTrkP3DIpSig->at(i) < v0Proton3dIpCut)continue;
      if( XiFLsig3D->at(i) < XiSepCut ) continue;
      if( Xi3DIpSig->at(i) > Xi3dIpCut )continue;
																		    
      if ( priNTrk ==0 ) continue;
       
      nCands++;
       
      //select best cand based on highest batchelor track pT
      float batPiPt = sqrt( (batPiPx->at(i)*batPiPx->at(i)) +  (batPiPy->at(i)*batPiPy->at(i)) );
      if (batPiPt>current_batPiPt) {
        current_batPiPt = batPiPt;
        bestCand = i;
      }
    }

    if (bestCand>-1) {

      float xiP = sqrt( (*xiPx)[bestCand]*(*xiPx)[bestCand]+(*xiPy)[bestCand]*(*xiPy)[bestCand]+(*xiPz)[bestCand]*(*xiPz)[bestCand]);
      float xiPt = sqrt( (*xiPx)[bestCand]*(*xiPx)[bestCand]+(*xiPy)[bestCand]*(*xiPy)[bestCand]);
      float xiE = sqrt( mXi*mXi + xiP*xiP );
      float xiY = fabs( 0.5*log( (xiE+(*xiPz)[bestCand])/(xiE-(*xiPz)[bestCand]) ) );
       
      data->Fill( (*xiMass)[bestCand]*1000);
      dataOm->Fill( (*omMass)[bestCand]*1000);

      for(unsigned int it = 0; it < thePtBins.size() - 1; ++it){
        if(xiPt > thePtBins.at(it) && xiPt < thePtBins.at(it + 1))thePtHistoVector.at(it)->Fill( (*xiMass)[bestCand]*1000);
      }
       
      for(unsigned int it = 0; it < theYBins.size() - 1; ++it){
        if(xiY > theYBins.at(it) && xiY < theYBins.at(it + 1))theYHistoVector.at(it)->Fill( (*xiMass)[bestCand]*1000);
      }
       
      if(fabs((*xictauMPV)[bestCand]) < theCTauBins.at(1))theCTauHistoVector.at(0)->Fill( (*xiMass)[bestCand]*1000 );
      for(unsigned int it = 1; it < theCTauBins.size() - 1; ++it){
        if (fabs((*xictauMPV)[bestCand]) > theCTauBins.at(it) && fabs((*xictauMPV)[bestCand]) < theCTauBins.at(it + 1))theCTauHistoVector.at(it) ->Fill( (*xiMass)[bestCand]*1000 );
      }
       
      //do asymmetry here
      if((*batPiQ)[bestCand] > 0){
        for(unsigned int it = 0; it < thePtBins.size() - 1; ++it){
	  if(xiPt > thePtBins.at(it) && xiPt < thePtBins.at(it + 1))thePosPtHistoVector.at(it)->Fill( (*xiMass)[bestCand]*1000);
	}
        for(unsigned int it = 0; it < theYBins.size() - 1; ++it){
          if(xiY > theYBins.at(it) && xiY < theYBins.at(it + 1))thePosYHistoVector.at(it)->Fill( (*xiMass)[bestCand]*1000);
        }
      }
      if((*batPiQ)[bestCand] < 0){
        for(unsigned int it = 0; it < thePtBins.size() - 1; ++it){
	  if(xiPt > thePtBins.at(it) && xiPt < thePtBins.at(it + 1))theNegPtHistoVector.at(it)->Fill( (*xiMass)[bestCand]*1000);
	}
        for(unsigned int it = 0; it < theYBins.size() - 1; ++it){
          if(xiY > theYBins.at(it) && xiY < theYBins.at(it + 1))theNegYHistoVector.at(it)->Fill( (*xiMass)[bestCand]*1000);
        }
      }
    }
   if (nCands>0) nEvtWCands++;
   if (nCands==1) nSingleCands++;
   nTotalCands += nCands;
   nXiCands->Fill( nCands );
   }
   
   cout << "Total candidates = " << nTotalCands << " from " << nEvtWCands << " events and single cand events = " << nSingleCands << endl; 
   
   TAxis* theAxis = data->GetXaxis();
   TCanvas* selectedCanvas = new TCanvas("selectedCanvas", "selectedCanvas");
   string theSelectedOutFile = theOutDir + "/XiMassSelected.eps";
   string theSelectedOutFileOm = theOutDir + "/OmMassSelected.eps"; 
   data->GetXaxis()->SetNdivisions(505);
   data->Draw();  
   int bins = theAxis->GetNbins();;
   int x_min = (int)theAxis->GetXmin();
   int x_max = (int)theAxis->GetXmax();
   int x_range = x_max-x_min;
   int binwidth = x_range/bins;
   TString title = TString("#Lambda^{0} #pi^{-} invariant mass [MeV/c^{2}]");
    // build ARGUS + Gaussian
   TF1 *f1 = new TF1("f1", Form("%s%i%s%i%s%i%s%i%s","[0]*((x-", x_min, ")**(0.5)) +[1]*((x-", x_min, ")**(3/2)) + ((", x_range, ")/", bins, ")*gausn(2)"), x_min, x_max);

   f1->SetParameters(0.55, 0., data->GetEntries()*.5, 1320, 60);
   f1->SetParName(0,"bkg param 1");
   f1->SetParName(1,"bkg param 2");
   f1->SetParName(2,"yield");
   f1->SetParLimits(2,0,200000);
   f1->SetParName(3,"mean");
   f1->SetParName(4,"sigma");
   f1->SetParLimits(4,0,1000);
   
   data->GetXaxis()->SetTitle(title);
   data->GetYaxis()->SetTitle(Form("%s%i%s","Candidates / ", binwidth, " MeV/c^{2}"));
   data->SetStats(0);



   f1->SetParameters(0.6, 0., 50, 1320, 40);

   data->Fit("f1","LEM0");
   //f1->SetMaximum(data->GetMaximum()*1.05);
   f1->GetXaxis()->SetTitle(title);
   f1->GetYaxis()->SetTitle(Form("%s%i%s","Candidates / ", binwidth, " MeV/c^{2}"));
   f1->SetNpx(5000);
   //f1->Draw("C");
   data->Draw();
   f1->Draw("same");

   float nSig = data->GetFunction("f1")->GetParameter(2);
   float mean = data->GetFunction("f1")->GetParameter(3);
   float sigma = data->GetFunction("f1")->GetParameter(4);
   float nSigE = data->GetFunction("f1")->GetParError(2);
   float meanE = data->GetFunction("f1")->GetParError(3);
   float sigmaE = data->GetFunction("f1")->GetParError(4);

   TPaveText *stats = new TPaveText(0.47,0.65,0.89,0.85,"blNDC");
   stats->AddText(Form("%s%.1f%s%.1f", "Yield: ", nSig, " #pm ",nSigE));
   stats->AddText(Form("%s%.2f%s%.2f%s", "Mean: ", mean, " #pm ", meanE, " MeV/c^{2}"));
   stats->AddText(Form("%s%.2f%s%.2f%s", "Sigma: ", sigma, " #pm ", sigmaE, " MeV/c^{2}"));
   stats->AddText("Statistical uncertainties only");
   stats->SetTextSize(0.04);
   stats->SetTextAlign(11);
   stats->Draw("same");
    
   //get background function and draw it
   float bkgA = data->GetFunction("f1")->GetParameter(0);
   float bkgB = data->GetFunction("f1")->GetParameter(1);
   TF1 *f3 = new TF1("f3", Form("%s%i%s%i%s","[0]*((x-", x_min, ")**(0.5)) +[1]*((x-", x_min, ")**(3/2))"), x_min, x_max);

   f3->SetParameter(0, bkgA);
   f3->SetParameter(1, bkgB);
   f3->SetLineStyle(2);
   f3->Draw("same");

   TLatex *cms = new TLatex();
   cms->SetNDC();
   cms->SetTextSize(0.040);
   cms->DrawLatex(0.60,0.89,"CMS Preliminary");
   selectedCanvas->Print(theSelectedOutFile.c_str());

   dataOm->Draw();
   selectedCanvas->Print(theSelectedOutFileOm.c_str());

   string thePtBinsFileName = theOutDir + "/ptBinsFit.eps";
   DoNumFit(thePtHistoVector, numFitPt, thePtBinsFileName, true);
   
   string thePosPtBinsFileName = theOutDir + "/posPtBinsFit.eps";
   DoNumFit(thePosPtHistoVector, posNumFitPt, thePosPtBinsFileName, true);

   
   string theNegPtBinsFileName = theOutDir + "/negPtBinsFit.eps";
   DoNumFit(theNegPtHistoVector, negNumFitPt, theNegPtBinsFileName, true);

   string thePosYBinsFileName = theOutDir + "/PosYBinsFit.eps";
   DoNumFit(thePosYHistoVector, posNumFitY, thePosYBinsFileName);
   
   string theNegYBinsFileName = theOutDir + "/negYBinsFit.eps";
   DoNumFit(theNegYHistoVector, negNumFitY, theNegYBinsFileName);
   
   TCanvas* asymmetryCanvas = new TCanvas();
   asymmetryCanvas->Divide(2,1);
   string asymmFileName = theOutDir + "/asymmetryPt.eps";
   
   asymmetryCanvas->cd(1);
   TGraphAsymmErrors* theAsymmetryPt = new TGraphAsymmErrors();
   GetAsymm(posNumFitPt, negNumFitPt, theAsymmetryPt);
   theAsymmetryPt->GetXaxis()->SetTitle("#Xi p_{T}");
   theAsymmetryPt->GetYaxis()->SetTitle("#Xi^{-} - #Xi^{-} / #Xi^{+} + #Xi^{-}");
   theAsymmetryPt->Draw("ap");
   
   asymmetryCanvas->cd(2);
   TGraphAsymmErrors* theAsymmetryY = new TGraphAsymmErrors();
   GetAsymm(posNumFitY, negNumFitY, theAsymmetryY);
   theAsymmetryY->GetXaxis()->SetTitle("#Xi rapidity");
   theAsymmetryY->GetYaxis()->SetTitle("#Xi^{-} - #Xi^{-} / #Xi^{+} + #Xi^{-}");
   theAsymmetryY->Draw("ap");
   asymmetryCanvas->Print(asymmFileName.c_str());

   string theYBinsFileName = theOutDir + "/rapidityBinsFit.eps";
   DoNumFit(theYHistoVector, numFitY, theYBinsFileName);

   string theCTauBinsFileName = theOutDir + "/cTauBinsFit.eps";
   DoNumFit(theCTauHistoVector, numFitCTau, theCTauBinsFileName);
   
   TCanvas* XiEffCtauCanvas = new TCanvas("XiEffCtauCanvas", "XiEffCtauCanvas");
   string theCTauEffFileName = theOutDir + "/XiEffCTau.eps";
   XiEffCtauCanvas->Divide(2,2);
   XiEffCtauCanvas->cd(1);
   numFitCTau->SetMarkerColor(2);

   numFitCTau->GetXaxis()->SetTitle("#Xi c#tau (cm)");
   numFitCTau->GetYaxis()->SetTitle("Unscaled entries");
   numFitCTau->Draw();
   
   XiEffCtauCanvas->cd(2);
   effCtau->SetMarkerColor(2);
   effCtau->SetMarkerSize(0.8);
   effCtau->SetMarkerStyle(20);
   effCtau->GetXaxis()->SetTitle("#Xi c#tau (cm)");
   effCtau->GetYaxis()->SetTitle("Efficiency");
   effCtau->Draw();
   
   XiEffCtauCanvas->cd(3);
   cout << "nbins in effctau, numfitctua and numeffctauscaled = " << effCtau->GetXaxis()->GetNbins() << ", " <<
   numFitCTau->GetXaxis()->GetNbins() << ", " << numEffCtauScaled->GetXaxis()->GetNbins() << endl;
   GetNumEff(effCtau, numFitCTau, numEffCtauScaled);
   for (int i = 1; i<= numEffCtauScaled->GetXaxis()->GetNbins(); ++i){
     cout << "for bin " << i << " numeffctauscaled = " << numEffCtauScaled->GetBinContent(i) << endl;
   }
   numEffCtauScaled->GetXaxis()->SetTitle("#Xi c#tau (cm)");   
   numEffCtauScaled->GetYaxis()->SetTitle("Scaled entries/cm");
   TF1 *exp1 = new TF1("exp1", Form("%s%f%s","[0]*exp(-x/([1]*", c, "))"), theCTauBins.front(), theCTauBins.back());
   exp1->SetParameters(numEffCtauScaled->GetMaximum(), 164.0);
   exp1->SetParName(0,"Constant");
   exp1->SetParName(1,"tau [ps]");
   numEffCtauScaled->Fit(exp1,"EMI");
   numEffCtauScaled->Draw();
   
   XiEffCtauCanvas->cd(4)->SetLogy();
   TH1F *numEffCtauScaledCopy = (TH1F*)numEffCtauScaled->Clone("numEffCtauScaledCopy");
   numEffCtauScaledCopy->SetMinimum(numEffCtauScaledCopy->GetBinContent(numEffCtauScaledCopy->GetXaxis()->GetNbins())*0.5);
   numEffCtauScaledCopy->Draw();
   
   XiEffCtauCanvas->Print(theCTauEffFileName.c_str());

// %%%%%%%%% Rapidity plots %%%%%%%%
   TCanvas* XiEffYCanvas = new TCanvas("XiEffYCanvas", "XiEffYCanvas");
   string theYEffFileName = theOutDir + "/XiEffY.eps";
   
   XiEffYCanvas->Divide(2,2);
   XiEffYCanvas->cd(1);
   
   numFitY->SetMarkerColor(2);
   numFitY->GetXaxis()->SetTitle("#Xi rapidity");
   numFitY->GetYaxis()->SetTitle("Unscaled entries");
   numFitY->Draw();
   
   XiEffYCanvas->cd(2);
   effY->SetMarkerColor(2);
   effY->SetMarkerSize(0.8);
   effY->SetMarkerStyle(20);
   effY->GetXaxis()->SetTitle("#Xi rapidity");
   effY->GetYaxis()->SetTitle("Efficiency");
   effY->Draw();
   
   XiEffYCanvas->cd(3);
   cout << "size of effY before call = " << effY->GetXaxis()->GetNbins() << endl;
   GetNumEff(effY, numFitY, numEffYScaled, true);
   numEffYScaled->GetXaxis()->SetTitle("#Xi rapidity");   
   numEffYScaled->GetYaxis()->SetTitle("Scaled entries/cm");
   numEffYScaled->Draw();
   
   XiEffYCanvas->Print(theYEffFileName.c_str());

// %%%%%%%%% Xi^+ Rapidity plots %%%%%%%%
   string thePosYEffFileName = theOutDir + "/XiPosEffY.eps";
   TCanvas* thePosYEffCanvas = new TCanvas();
   thePosYEffCanvas->Divide(2,2);
   thePosYEffCanvas->cd(1);
   
   posNumFitY->SetMarkerColor(2);
   posNumFitY->GetXaxis()->SetTitle("#Xi rapiditPosXiY");
   posNumFitY->GetYaxis()->SetTitle("Unscaled entries");
   posNumFitY->Draw();
   
   thePosYEffCanvas->cd(2);
   effPosXiY->SetMarkerColor(2);
   effPosXiY->SetMarkerSize(0.8);
   effPosXiY->SetMarkerStyle(20);
   effPosXiY->GetXaxis()->SetTitle("#Xi rapiditPosXiY");
   effPosXiY->GetYaxis()->SetTitle("EfficiencPosXiY");
   effPosXiY->Draw();
   
   thePosYEffCanvas->cd(3);
   cout << "size of effPosXiY before call = " << effPosXiY->GetXaxis()->GetNbins() << endl;
   GetNumEff(effPosXiY, posNumFitY, numEffPosXiYScaled, true);
   numEffPosXiYScaled->GetXaxis()->SetTitle("#Xi rapiditPosXiY");   
   numEffPosXiYScaled->GetYaxis()->SetTitle("Scaled entries/cm");
   numEffPosXiYScaled->Draw();
   thePosYEffCanvas->Print(thePosYEffFileName.c_str());

// %%%%%%%%%% Xi^- rapidity plots    
   string theNegYEffFileName = theOutDir + "/XiNegEffY.eps";
   TCanvas* theNegYEffCanvas = new TCanvas();
   theNegYEffCanvas->Divide(2,2);
   theNegYEffCanvas->cd(1);
   
   negNumFitY->SetMarkerColor(2);
   negNumFitY->GetXaxis()->SetTitle("#Xi rapiditNegXiY");
   negNumFitY->GetYaxis()->SetTitle("Unscaled entries");
   negNumFitY->Draw();
   
   theNegYEffCanvas->cd(2);
   effNegXiY->SetMarkerColor(2);
   effNegXiY->SetMarkerSize(0.8);
   effNegXiY->SetMarkerStyle(20);
   effNegXiY->GetXaxis()->SetTitle("#Xi rapiditNegXiY");
   effNegXiY->GetYaxis()->SetTitle("EfficiencNegXiY");
   effNegXiY->Draw();
   
   theNegYEffCanvas->cd(3);
   cout << "size of effNegXiY before call = " << effNegXiY->GetXaxis()->GetNbins() << endl;
   GetNumEff(effNegXiY, negNumFitY, numEffNegXiYScaled, true);
   numEffNegXiYScaled->GetXaxis()->SetTitle("#Xi rapiditNegXiY");   
   numEffNegXiYScaled->GetYaxis()->SetTitle("Scaled entries/cm");
   numEffNegXiYScaled->Draw();
   theNegYEffCanvas->Print(theNegYEffFileName.c_str());
   
// %%%%%%%%%%% Pt Plots
   TCanvas* XiEffPtCanvas = new TCanvas("XiEffPtCanvas", "XiEffPtCanvas");
   string thePtEffFileName = theOutDir + "/XiEffPt.eps";
   XiEffPtCanvas->Divide(2,2);
   XiEffPtCanvas->cd(1);
   
   numFitPt->SetMarkerColor(2);
   numFitPt->GetXaxis()->SetTitle("#Xi p_{T}");
   numFitPt->GetYaxis()->SetTitle("Unscaled entries");
   numFitPt->Draw();
   
   XiEffPtCanvas->cd(2);
   effPt->SetMarkerColor(2);
   effPt->SetMarkerSize(0.8);
   effPt->SetMarkerStyle(20);
   effPt->GetXaxis()->SetTitle("#Xi p_{T}");
   effPt->GetYaxis()->SetTitle("Efficiency");
   effPt->Draw();
   
   XiEffPtCanvas->cd(3);
   GetNumEff(effPt, numFitPt, numEffpTScaled);
   numEffpTScaled->GetXaxis()->SetTitle("#Xi p_{T}");   
   numEffpTScaled->GetYaxis()->SetTitle("Scaled entries/cm");
   numEffpTScaled->Draw();
   
   XiEffPtCanvas->cd(4)->SetLogy();
   TH1F* numEffpTScaledCopy = (TH1F*)numEffpTScaled->Clone("numEffpTScaledCopy");
   numEffpTScaledCopy->SetMinimum(numEffpTScaledCopy->GetBinContent(numEffpTScaledCopy->GetXaxis()->GetNbins())*0.1);
   numEffpTScaledCopy->Draw();
   
   XiEffPtCanvas->Print(thePtEffFileName.c_str());

// %%%%%% Xi^+ Pt Plots
   
   string thePosPtEffFileName = theOutDir + "/XiPosEffPt.eps";
   TCanvas* thePosPtEffCanvas = new TCanvas();
   thePosPtEffCanvas->Divide(2,2);
   thePosPtEffCanvas->cd(1);
   
   posNumFitPt->SetMarkerColor(2);
   posNumFitPt->GetXaxis()->SetTitle("#Xi p_{T}");
   posNumFitPt->GetYaxis()->SetTitle("Unscaled entries");
   posNumFitPt->Draw();
   
   thePosPtEffCanvas->cd(2);
   effPosXiPt->SetMarkerColor(2);
   effPosXiPt->SetMarkerSize(0.8);
   effPosXiPt->SetMarkerStyle(20);
   effPosXiPt->GetXaxis()->SetTitle("#Xi p_{T}");
   effPosXiPt->GetYaxis()->SetTitle("Efficiency");
   effPosXiPt->Draw();
   
   thePosPtEffCanvas->cd(3);
   //numEffPosXipTScaled = GetNumEff(effPosXiPt, posNumFitPt);
   GetNumEff(effPosXiPt, posNumFitPt, numEffPosXipTScaled);
   numEffPosXipTScaled->GetXaxis()->SetTitle("#Xi p_{T}");   
   numEffPosXipTScaled->GetYaxis()->SetTitle("Scaled entries/cm");
   numEffPosXipTScaled->Draw();
   
   //numEffPosXipTScaledUncor = GetNumEff(effPtUncor, numFitPt);
   thePosPtEffCanvas->cd(4)->SetLogy();
   numEffPosXipTScaled->Draw();
   
   thePosPtEffCanvas->Print(thePosPtEffFileName.c_str());

// %%%%%%%%% X^- Pt Plots
   
   string theNegPtEffFileName = theOutDir + "/XiNegEffPt.eps";
   TCanvas* theNegPtEffCanvas = new TCanvas();
   theNegPtEffCanvas->Divide(2,2);
   theNegPtEffCanvas->cd(1);
   
   negNumFitPt->SetMarkerColor(2);
   negNumFitPt->GetXaxis()->SetTitle("#Xi p_{T}");
   negNumFitPt->GetYaxis()->SetTitle("Unscaled entries");
   negNumFitPt->Draw();
   
   theNegPtEffCanvas->cd(2);
   effNegXiPt->SetMarkerColor(2);
   effNegXiPt->SetMarkerSize(0.8);
   effNegXiPt->SetMarkerStyle(20);
   effNegXiPt->GetXaxis()->SetTitle("#Xi p_{T}");
   effNegXiPt->GetYaxis()->SetTitle("Efficiency");
   effNegXiPt->Draw();
   
   theNegPtEffCanvas->cd(3);
   //numEffNegXipTScaled = GetNumEff(effNegXiPt, negNumFitPt);
   GetNumEff(effNegXiPt, negNumFitPt, numEffNegXipTScaled);
   numEffNegXipTScaled->GetXaxis()->SetTitle("#Xi p_{T}");   
   numEffNegXipTScaled->GetYaxis()->SetTitle("Scaled entries/cm");
   numEffNegXipTScaled->Draw();
   
   //numEffNegXipTScaledUncor = GetNumEff(effPtUncor, numFitPt);
   theNegPtEffCanvas->cd(4)->SetLogy();
   numEffNegXipTScaled->Draw();
   
   theNegPtEffCanvas->Print(theNegPtEffFileName.c_str());
   
   
}


void seleXiData::DoNumFit(vector<TH1F*> theHistoVector, TH1F* theHisto, string plotName, bool restrictRange){
  TCanvas* theCanvas = new TCanvas();
  theCanvas->Divide(4,5);
  for(unsigned int vIt = 0; vIt < theHistoVector.size(); ++vIt){
    theCanvas->cd(vIt + 1);
    if(vIt<2&&restrictRange)
      FitXi(theHistoVector.at(vIt), true);
    else
      FitXi(theHistoVector.at(vIt));
    theHisto->SetBinContent(vIt + 1, theFitYield);
    theHisto->SetBinError(vIt + 1, theFitYieldError);
  }
  theCanvas->Print(plotName.c_str());
  delete theCanvas;
}

void seleXiData::GetNumEff(TH1F* effHisto, TH1F* numFitHisto, TH1F* theOutHisto, bool doFolded){
  cout << " from before fill: size of theOutHisto = " << theOutHisto->GetXaxis()->GetNbins() << " numFitHisto = " 
       << numFitHisto->GetXaxis()->GetNbins() << " effHisto = " << effHisto->GetXaxis()->GetNbins() << endl;

  TH1F* widthHisto = (TH1F*)effHisto->Clone("widthHisto");
  for (int i = 1; i <= effHisto->GetXaxis()->GetNbins(); i++) {
    float width = effHisto->GetXaxis()->GetBinWidth(i);
    if(doFolded) width = width*2;
    widthHisto->SetBinContent(i, width);
    widthHisto->SetBinError(i,0);
    cout << "for bin " << i << " width = " << width << " numfit = " << numFitHisto->GetBinContent(i) << " eff = " <<
    effHisto->GetBinContent(i) << endl;
  }
  theOutHisto->Divide(numFitHisto, effHisto);
  theOutHisto->Divide(theOutHisto, widthHisto);  
}

void seleXiData::GetAsymm(TH1F* thePosYieldsHisto, TH1F* theNegYieldsHisto, TGraphAsymmErrors* theAsymmGraph){
  TAxis* theXAxis = thePosYieldsHisto->GetXaxis();
  unsigned int theBins = theXAxis->GetNbins();
  float lowLimit = 0.;
  float hiLimit = 0.;
  for(unsigned int vIt = 0; vIt < theBins - 1; ++vIt){
    if(vIt == 0) lowLimit = theXAxis->GetBinLowEdge(vIt + 1);
    double asymmetryMinus = theNegYieldsHisto->GetBinContent(vIt + 1) - thePosYieldsHisto->GetBinContent(vIt + 1);							  
    double asymmetryPlus = thePosYieldsHisto->GetBinContent(vIt + 1) + theNegYieldsHisto->GetBinContent(vIt + 1);							  
    double asymmetry = asymmetryMinus / asymmetryPlus;													  
    theAsymmGraph->SetPoint(vIt, theXAxis->GetBinCenter(vIt + 1), asymmetry);							  
    theAsymmGraph->SetPointEXlow(vIt,theXAxis->GetBinCenter(vIt + 1)-theXAxis->GetBinLowEdge(vIt + 1));					  
    theAsymmGraph->SetPointEXhigh(vIt,theXAxis->GetBinLowEdge(vIt + 1)-theXAxis->GetBinCenter(vIt + 1)+theXAxis->GetBinWidth(vIt + 1)); 
    double plusminusErrorSqrd = thePosYieldsHisto->GetBinError(vIt + 1) * thePosYieldsHisto->GetBinError(vIt + 1) + 						  
    			       theNegYieldsHisto->GetBinError(vIt + 1) * theNegYieldsHisto->GetBinError(vIt + 1);							  
    double asymmError = (sqrt(plusminusErrorSqrd)/asymmetryPlus)*sqrt(1. + asymmetry*asymmetry);							  
    theAsymmGraph->SetPointEYhigh(vIt,asymmError);													  
    theAsymmGraph->SetPointEYlow(vIt,asymmError);
    hiLimit = theXAxis->GetBinUpEdge(vIt + 1);											  
  }
  theAsymmGraph->GetXaxis()->SetLimits(lowLimit, hiLimit); 
  theAsymmGraph->SetMaximum(1); 					     
  theAsymmGraph->SetMinimum(-1);					     
  theAsymmGraph->SetMarkerColor(2);					     
  theAsymmGraph->SetMarkerSize(0.8);					     
  theAsymmGraph->SetMarkerStyle(20);					     
}

void seleXiData::FitXi(TH1F* theXiHisto, bool restrictRange){
    //initialize results to 0
  theFitYield = 0.;
  theFitYieldError = 0.;
  
  TAxis* theAxis = theXiHisto->GetXaxis();
  int x_min = (int)theAxis->GetXmin();
  int x_max = (int)theAxis->GetXmax();
  int bins = theAxis->GetNbins();
  int x_range = x_max - x_min;
  TF1* xiFit = new TF1("xiFit", Form("%s%i%s%i%s%i%s","[0]*((x-", x_min, ")**([1])) + ((", x_range, ")/", bins, ")*gausn(2)"), x_min, x_max);
  xiFit->SetParameters(1.5, 0.1, theXiHisto->GetEntries()*.5, 1321.7, 3);
  xiFit->SetParName(0,"bkg param 1");
  xiFit->SetParName(1,"bkg param 2");
  xiFit->SetParName(2,"yield");
  xiFit->SetParLimits(2,0,200000);
  xiFit->SetParName(3,"mean");
  xiFit->SetParLimits(3,1318.0,1330.0);
  xiFit->SetParName(4,"sigma");
  xiFit->SetParLimits(4,1,6);

  if (restrictRange) {
    xiFit->SetParameters(20.0, -0.5, theXiHisto->GetEntries()*.5, 1321, 3);
    theXiHisto->Fit("xiFit","LEM","",1290,x_max);
  }
  else
    theXiHisto->Fit("xiFit","LEM");

  theFitYield = xiFit->GetParameter(2);
  theFitYieldError = xiFit->GetParError(2);

  if (theFitYieldError > theFitYield) {
    cout << "Fit Failed: CHECK IT. Setting error to have of value for now" << endl;
    theFitYieldError = theFitYield*0.5;
  }
  
}


void seleXiData::SetBins(){
   
   // Pt bins
   thePtBins.push_back(0.0);
   thePtBins.push_back(0.2);
   thePtBins.push_back(0.4);
   thePtBins.push_back(0.6);
   thePtBins.push_back(0.8);
   thePtBins.push_back(1.0);
   thePtBins.push_back(1.2);
   thePtBins.push_back(1.4);
   thePtBins.push_back(1.6);
   thePtBins.push_back(1.8);
   thePtBins.push_back(2.0);
   thePtBins.push_back(2.2);
   thePtBins.push_back(2.4);
   thePtBins.push_back(2.6);
   thePtBins.push_back(2.8);
   thePtBins.push_back(3.0);
   thePtBins.push_back(3.2);
   thePtBins.push_back(3.4);
   thePtBins.push_back(3.6);
   thePtBins.push_back(3.8);
   thePtBins.push_back(4.0);
   thePtBins.push_back(5.0);
   thePtBins.push_back(6.0);
   
   // rapidity
   theYBins.push_back(0.0);
   theYBins.push_back(0.2);
   theYBins.push_back(0.4);
   theYBins.push_back(0.6);
   theYBins.push_back(0.8);
   theYBins.push_back(1.0);
   theYBins.push_back(1.2);
   theYBins.push_back(1.4);
   theYBins.push_back(1.6);
   theYBins.push_back(1.8);
   theYBins.push_back(2.0);
   
   //ctau
   theCTauBins.push_back(0.0);
   theCTauBins.push_back(0.4);
   theCTauBins.push_back(0.8);
   theCTauBins.push_back(1.2);
   theCTauBins.push_back(1.6);
   theCTauBins.push_back(2.0);
   theCTauBins.push_back(2.5);
   theCTauBins.push_back(3.0);
   theCTauBins.push_back(4.0);
   theCTauBins.push_back(6.0);
   theCTauBins.push_back(10.0);
}

void seleXiData::bookHistos(){

  theEffFile = new TFile(theEffFileName.c_str());
   
  effCtau = (TH1F*)theEffFile->Get("effCtau");
  effCtauUncor = (TH1F*)theEffFile->Get("effCtauUncor");
   
  effY = (TH1F*)theEffFile->Get("effY");
  effYUncor = (TH1F*)theEffFile->Get("effYUncor");
   
  effPosXiY = (TH1F*)theEffFile->Get("effPosXiY");
  effPosXiYUncor = (TH1F*)theEffFile->Get("effPosXiYUncor");
   
  effNegXiY = (TH1F*)theEffFile->Get("effNegXiY");
  effNegXiYUncor = (TH1F*)theEffFile->Get("effNegXiYUncor");
    
  effPt = (TH1F*)theEffFile->Get("effpT");
  effPtUncor = (TH1F*)theEffFile->Get("effpTUncor");
   
  effPosXiPt = (TH1F*)theEffFile->Get("effPosXipT");
  effPosXiPtUncor = (TH1F*)theEffFile->Get("effPosXipTUncor");
   
  effNegXiPt = (TH1F*)theEffFile->Get("effNegXipT");
  effNegXiPtUncor = (TH1F*)theEffFile->Get("effNegXipTUncor");
   
  // open output root file
  std::string theOutRootFileName = theOutDir + "/xiOut.root";
  theOutRootFile = new TFile(theOutRootFileName.c_str(), "Recreate");
  theOutRootFile->cd();
  
  nXiCands = new TH1F("nXiCands","nXiCands", 20, 0, 20);
  
  for(unsigned int hIt = 0; hIt < thePtBins.size() - 1; ++hIt){
    int bins = 100;
    if (hIt==0) bins = 50;
    
    char* theHistoName = Form("ptBin%d",hIt+1);
    thePtHistoVector.push_back(new TH1F(theHistoName, theHistoName, bins, 1255, 1455));
    
    char* thePosHistoName = Form("posPtBin%d",hIt+1);
    thePosPtHistoVector.push_back(new TH1F(thePosHistoName, thePosHistoName, bins, 1255, 1455));
    
    char* theNegHistoName = Form("negPtBin%d",hIt+1);
    theNegPtHistoVector.push_back(new TH1F(theNegHistoName, theNegHistoName, bins, 1255, 1455));
  }
  
  for(unsigned int hIt = 0; hIt < theYBins.size() - 1; ++hIt){
    char* theHistoName = Form("yBin%d",hIt+1);
    theYHistoVector.push_back(new TH1F(theHistoName, theHistoName, 100, 1255, 1455));
    
    char* thePosHistoName = Form("posYBin%d",hIt+1);
    thePosYHistoVector.push_back(new TH1F(thePosHistoName, theHistoName, 100, 1255, 1455));
    
    char* theNegHistoName = Form("negYBin%d",hIt+1);
    theNegYHistoVector.push_back(new TH1F(theNegHistoName, theHistoName, 100, 1255, 1455));
  }
  
  for(unsigned int hIt = 0; hIt < theCTauBins.size() - 1; ++hIt){
    char* theHistoName = Form("CTauBin%d",hIt+1);
    theCTauHistoVector.push_back(new TH1F(theHistoName, theHistoName, 100, 1255, 1455));
  }
  
  data = new TH1F("data", "data", 100, 1255, 1455);
  dataOm = new TH1F("dataOm", "dataOm", 100, 1455, 1855);
  
  numEffpTScaled = new TH1F("numEffpTScaled","numEffpTScaled", thePtBins.size() - 1, &thePtBins[0]);
  numEffpTScaledUncor = new TH1F("numEffpTScaledUncor","numEffpTScaledUncor", thePtBins.size() - 1, &thePtBins[0]);
  numFitPt = new TH1F("numFitPt","numFitPt", thePtBins.size() - 1, &thePtBins[0]);
  
  numEffPosXipTScaled = new TH1F("numEffPosXipTScaled","numEffPosXipTScaled", thePtBins.size() - 1, &thePtBins[0]);
  numEffPosXipTScaledUncor = new TH1F("numEffPosXipTScaledUncor","numEffPosXipTScaledUncor", thePtBins.size() - 1, &thePtBins[0]);
  
  numEffNegXipTScaled = new TH1F("numEffNegXipTScaled","numEffNegXipTScaled", thePtBins.size() - 1, &thePtBins[0]);
  numEffNegXipTScaledUncor = new TH1F("numEffNegXipTScaledUncor","numEffNegXipTScaledUncor", thePtBins.size() - 1, &thePtBins[0]);
  
  posNumFitPt = new TH1F("posNumFitPt","posNumFitPt", thePtBins.size() - 1, &thePtBins[0]);
  negNumFitPt = new TH1F("negNumFitPt","negNumFitPt", thePtBins.size() - 1, &thePtBins[0]);
  asymmFitPt = new TH1F("asymmFitPt","asymmFitPt", thePtBins.size() - 1, &thePtBins[0]);
  
  numEffYScaled = new TH1F("numEffYScaled","numEffYScaled",theYBins.size() - 1, &theYBins[0]);
  numEffYScaledUncor = new TH1F("numEffYScaledUncor","numEffYScaledUncor", theYBins.size() - 1, &theYBins[0]);
  numFitY = new TH1F("numFitY","numFitY", theYBins.size() - 1, &theYBins[0]); 
  
  numEffPosXiYScaled = new TH1F("numEffPosXiYScaled","numEffPosXiYScaled", theYBins.size() - 1, &theYBins[0]);
  numEffPosXiYScaledUncor = new TH1F("numEffPosXiYScaledUncor","numEffPosXiYScaledUncor", theYBins.size() - 1, &theYBins[0]);
  
  numEffNegXiYScaled = new TH1F("numEffNegXiYScaled","numEffNegXiYScaled", theYBins.size() - 1, &theYBins[0]);
  numEffNegXiYScaledUncor = new TH1F("numEffNegXiYScaledUncor","numEffNegXiYScaledUncor", theYBins.size() - 1, &theYBins[0]);
  
  posNumFitY = new TH1F("posNumFitY","posNumFitY", theYBins.size() - 1, &theYBins[0]);
  negNumFitY = new TH1F("negNumFitY","negNumFitY", theYBins.size() - 1, &theYBins[0]);
  asymmFitY = new TH1F("asymmFitY","asymmFitY", theYBins.size() - 1, &theYBins[0]);
 
  numEffCtauScaled = new TH1F("numEffCtauScaled","numEffCtauScaled", theCTauBins.size() - 1, &theCTauBins[0]);
  numEffCtauScaledUncor = new TH1F("numEffCtauScaledUncor","numEffCtauScaledUncor", theCTauBins.size() - 1, &theCTauBins[0]);
  numFitCTau = new TH1F("numFitCTau","numFitCTau", theCTauBins.size() - 1, &theCTauBins[0]);   
  
}

void seleXiData::CleanVectors(){
  for(unsigned int it = 0; it < thePtHistoVector.size(); ++it){
    delete thePtHistoVector.at(it);
  }
  
  for(unsigned int it = 0; it < theYHistoVector.size(); ++it){
    delete theYHistoVector.at(it);
  }
  
  for(unsigned int it = 0; it < theCTauHistoVector.size(); ++it){
    delete theCTauHistoVector.at(it);
  }
  
  for(unsigned int it = 0; it < theFitVector.size(); ++it){
    delete theFitVector.at(it);
  }
  
  thePtBins.clear();
  theYBins.clear();
  theCTauBins.clear();
  
}
void seleXiData::SetEffFile(std::string theFile){
  theEffFileName = theFile;
}

void seleXiData::SetOutputDir(std::string theOutputDir){
  struct stat stFileInfo;
  int intStat;

  // Attempt to get the file attributes
  intStat = stat(theOutputDir.c_str(),&stFileInfo);
  if(intStat == 0) {
    theOutDir = theOutputDir;
  } else {
  // make dir
    intStat = mkdir(theOutputDir.c_str(), 0777);
    if(intStat == 0) theOutDir = theOutputDir;
    else std::cout << "Failed to created desired directory, using Default \n";
  }
}

