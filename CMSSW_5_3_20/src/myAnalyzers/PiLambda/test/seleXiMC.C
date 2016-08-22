#define seleXiMC_cxx
#include "seleXiMC.h"
#include "/nfs/data37/cms/v0ntuple_10_2010/ntrk_wgt.C"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <math.h>
#include <fstream>
#include <TF1.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <sys/stat.h> 
#include <sys/types.h>

using namespace std;

void seleXiMC::Loop(){
   
  cout << " Cuts Used:\t " << "Xi3dIp " << Xi3dIpCut << "\t batPi3dIp " << batPi3dIpCut << "\t v0Pi3dIp" << v0Pi3dIpCut 
                        << "\t v0Proton3dIp "  << v0Proton3dIpCut << "\t XiSep " << XiSepCut << "\t VtxCL " << vtxCLCut << "\n";
 
  
  SetBins();
  BookHistos();

  TStyle* gStyle = new TStyle();

  gStyle->SetStatFormat("5.5g");
  gStyle->SetStatY(0.94);
  gStyle->SetFitFormat("5.5g");
  gStyle->SetNdivisions(505, "XYZ");
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.03);
  gStyle->SetLabelOffset(0.010, "X");
  gStyle->SetTitleOffset(1.1,"Y");
  gStyle->SetFillColor(0);
  gStyle->SetTitleSize(0.055, "XYZ");
  gStyle->SetLabelSize(0.045, "XYZ");
  gStyle->SetTitleXOffset(1.0);
  gStyle->SetTitleYOffset(1.25);  
  gStyle->SetErrorX(0.5); 
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;   
  if(nEvents < 0) nentries = fChain->GetEntriesFast();
  else
     nentries = nEvents;   
   
  int nBins = 100;
  float min = 0.0;
  float max = 10.0;
  TH1F *numEff = new TH1F("numEff","numEff",nBins,min,max);
  TH1F *denEff = new TH1F("denEff","denEff",nBins,min,max);
  TH1F *effH = new TH1F("effH","effH",nBins,min,max);
  float mXi = 1.32171;
  float pi = TMath::Pi();
  float c = 0.0299793; // cm/ps
  int nCands = 0;
  int nSingleCands = 0;
  int nTotalCands = 0;
  int nEvtWCands = 0;

  int weight_choice = weightChoice;

  // first get efficiency vs ctau
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    double weightMC_ntrk = ntrk_wgt(nLooseTracks, weight_choice);
    
    for ( int i = 0; i<nGenXi; i++) {
      float E = sqrt( mXi*mXi + (*genXiP)[i]*(*genXiP)[i]);
      float theta = 2*atan(exp((*genXiEta)[i]));
      float pT = sin(theta)*(*genXiP)[i];
      float pL = cos(theta)*(*genXiP)[i];
      float Y = 0.5*log((E+pL)/(E-pL)); 
      
      if ( fabs(Y) > 2.0) continue;
      
      double weightMC_y = xipt_wgt( pT, weight_choice );
      double weightMC_pt = xiy_wgt( Y, weight_choice );

      if(weightType == "full"){
        theYWeight = weightMC_ntrk*weightMC_y;
	thePtWeight = weightMC_ntrk*weightMC_pt;
      }
      else{
        theYWeight = weightMC_ntrk;
	thePtWeight = weightMC_ntrk;
      }
      
            // for rapidity and pT we can consider all simulated Xi
      if ( !(processType==92||processType==93)) {
        denEffY->Fill( fabs(Y), theYWeight);
        denEffYUncor->Fill( fabs(Y) );
        allGenY->Fill( fabs(Y) );
        denEffpT->Fill(pT, thePtWeight);
        denEffpTUncor->Fill(pT);
        allGenpT->Fill(pT);
	
	//std::cout << " pion charge " << pionCharge << std::endl;      
        // for ctau it's only possible to consider the MC events with at least one simulated track
        if ( (*genXiL)[i] < -1 ) continue;
        float ctau = ((*genXiL)[i]*mXi)/((*genXiP)[i]);
        denEff->Fill(ctau);
        denEffCtau->Fill(ctau, theYWeight);
        denEffCtauUncor->Fill(ctau);
        
	//neg Xi
	if((*genXiQ)[i] < 0){
          denEffNegXipT->Fill(pT, thePtWeight);
          denEffNegXipTUncor->Fill(pT);
          allGenNegXipT->Fill(pT);
          denEffNegXiY->Fill( fabs(Y), theYWeight);
          denEffNegXiYUncor->Fill( fabs(Y) );
          allGenNegXiY->Fill( fabs(Y) );
	  
	}
	
	// pos Xi
	if((*genXiQ)[i] > 0){
          denEffPosXipT->Fill(pT, thePtWeight);
          denEffPosXipTUncor->Fill(pT);
          allGenPosXipT->Fill(pT);	   
          denEffPosXiY->Fill( fabs(Y), theYWeight);
          denEffPosXiYUncor->Fill( fabs(Y) );
          allGenPosXiY->Fill( fabs(Y) );
	}
      }
    }

 
        
    nCands = 0;	
    int bestCand = -999;
    float current_batPiPt = 0.0;
    
    for ( int i = 0; i<nXI; i++ ) {
    
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
      if ( (*xiVtxCL)[i] <  vtxCLCut) continue;
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

      data->Fill( (*xiMass)[bestCand]*1000 );
      dataOm->Fill( (*omMass)[bestCand]*1000 );
      //std::cout << "Filling with number of cands = " << nCands << "  xiMass = " << (*xiMass)[bestCand] << "  cand number " << bestCand << std::endl;
      float xiP = sqrt( (*xiPx)[bestCand]*(*xiPx)[bestCand]+(*xiPy)[bestCand]*(*xiPy)[bestCand]+(*xiPz)[bestCand]*(*xiPz)[bestCand]);
      float xiE = sqrt( mXi*mXi + xiP*xiP );
      float xiPt = sqrt( (*xiPx)[bestCand]*(*xiPx)[bestCand]+(*xiPy)[bestCand]*(*xiPy)[bestCand]);
      float xiY = fabs( 0.5*log( (xiE+(*xiPz)[bestCand])/(xiE-(*xiPz)[bestCand]) ) );
      
      // what to do about negative values?
      // for efficiency we take  the absolute value, so do that here
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
      
      float delRXiCut = 5.0;

      float batPiPt = sqrt( (batPiPx->at(bestCand)*batPiPx->at(bestCand)) +  (batPiPy->at(bestCand)*batPiPy->at(bestCand)) );
      float lamPiPt = sqrt( (VTrkPiPx->at(bestCand)*VTrkPiPx->at(bestCand)) +  (VTrkPiPy->at(bestCand)*VTrkPiPy->at(bestCand)) );
      float lamPPt = sqrt( (VTrkPPx->at(bestCand)*VTrkPPx->at(bestCand)) +  (VTrkPPy->at(bestCand)*VTrkPPy->at(bestCand)) );

      float delRXiPi = (0.01/batPiPt)+0.02;
      float delRLamPi = (0.01/lamPiPt)+0.02;
      float delRP = (0.01/lamPPt)+0.02;
      
      bool truthMatched = false;
      float genctau = -99.0;
      float genpT = -99.0;
      float genY = -99.0;
      for ( int j = 0; j<nGenXi; j++ ) {
        float delRbatPi = GetDeltaR((*batPiEta)[bestCand],(*batPiPhi)[bestCand], (*genBatPiEta)[j], (*genBatPiPhi)[j]);
        if (delRbatPi > pi) delRbatPi -= 2*pi;
        if (delRbatPi < -pi) delRbatPi += 2*pi; 	
    	
	float delRlamPi = GetDeltaR((*VTrkPiEta)[bestCand], (*VTrkPiPhi)[bestCand],(*genLamPiEta)[j], (*genLamPiPhi)[j]);
        if (delRlamPi > pi) delRlamPi -= 2*pi;
        if (delRlamPi < -pi) delRlamPi += 2*pi; 					 
        
	float delRlamP = GetDeltaR((*VTrkPEta)[bestCand], (*VTrkPPhi)[bestCand], (*genLamPEta)[j], (*genLamPPhi)[j]); 
        if (delRlamP > pi) delRlamP -= 2*pi;
        if (delRlamP < -pi) delRlamP += 2*pi;				       
        
	float delRXi = sqrt( ((*xiDecayVtxX)[bestCand]-(*genXiDecayVX)[j])*((*xiDecayVtxX)[bestCand]-(*genXiDecayVX)[j])+
                             ((*xiDecayVtxY)[bestCand]-(*genXiDecayVY)[j])*((*xiDecayVtxY)[bestCand]-(*genXiDecayVY)[j])+
			     ((*xiDecayVtxZ)[bestCand]-(*genXiDecayVZ)[j])*((*xiDecayVtxZ)[bestCand]-(*genXiDecayVZ)[j]));   
		   
        if ( (delRbatPi < delRXiPi) &&
             (delRlamPi < delRLamPi) && 
             (delRlamP < delRP) 
//	     && (delRXi < delRXiCut) 
	     ) {
          truthMatched = true;
          genctau = ((*genXiL)[j]*mXi)/((*genXiP)[j]);
          float E = sqrt( mXi*mXi + (*genXiP)[j]*(*genXiP)[j]);
          float theta = 2*atan(exp((*genXiEta)[j]));
          genpT = sin(theta)*(*genXiP)[j];
          float pL = cos(theta)*(*genXiP)[j];
          genY = fabs(0.5*log((E+pL)/(E-pL)));        	
	}
      }
      if (truthMatched) {  
        //std::cout << "Filling with ctau = " << (*xictauMPV)[bestCand] << std::endl;
        double weightMC_ntrk = ntrk_wgt(nLooseTracks, weight_choice); 

        // use true values of pT and Y to get weights (point is to get MC back to correct distribution)
        double weightMC_y = xipt_wgt( genpT, weight_choice );
        double weightMC_pt = xiy_wgt( genY, weight_choice );
	
	if(weightType == "full"){
      	  theYWeight = weightMC_ntrk*weightMC_y;
	  thePtWeight = weightMC_ntrk*weightMC_pt;
      	}
      	else{
      	  theYWeight = weightMC_ntrk;
	  thePtWeight = weightMC_ntrk;
      	}
	
	//this can be negative and my histo cuts off at 0.
        numEffCtau->Fill( fabs( (*xictauMPV)[bestCand] ), theYWeight);
        numEffCtauUncor->Fill( fabs( (*xictauMPV)[bestCand] ) );
	
	numEffY->Fill( xiY, theYWeight);
	numEffYUncor->Fill( xiY );	
	
	numEffpT->Fill( xiPt, thePtWeight);
        numEffpTUncor->Fill( xiPt );
	
	//pos Xi
	if((*batPiQ)[bestCand] > 0){
	  numEffPosXipT->Fill( xiPt, thePtWeight); 
          numEffPosXipTUncor->Fill( xiPt );			     
	  numEffPosXiY->Fill( xiY, theYWeight);
	  numEffPosXiYUncor->Fill( xiY );	  
	
	  
	}
	//neg Xi
	if((*batPiQ)[bestCand] < 0){
	  numEffNegXipT->Fill( xiPt, thePtWeight);
          numEffNegXipTUncor->Fill( xiPt );
	  numEffNegXiY->Fill( xiY, theYWeight); 
	  numEffNegXiYUncor->Fill( xiY );	  		  
	
	}
		
        numEff->Fill( fabs( (*xictauMPV)[bestCand] ) );
        dataTruth->Fill( (*xiMass)[bestCand]*1000 );
	if(xiPt < 0.2)dataTruthLowPt->Fill( (*xiMass)[bestCand]*1000 );

        // not related to efficiency, just for ctau comparison
        ctauMPV->Fill( (*xictauMPV)[bestCand] );
        ctauMPVRes->Fill( (*xictauMPV)[bestCand] - genctau );
        ctauMPVPull->Fill( ((*xictauMPV)[bestCand] - genctau)/(*xictauMPVE)[bestCand] );
        ctauMPVAbs->Fill( fabs( (*xictauMPV)[bestCand] ) );
        ctauMPVAbsRes->Fill( fabs( (*xictauMPV)[bestCand] ) - genctau );
        ctauMPVAbsPull->Fill( (fabs( (*xictauMPV)[bestCand] ) - genctau)/(*xictauMPVE)[bestCand] );
        ctau3D->Fill( (*xictau3D)[bestCand] );
        ctau3DRes->Fill( (*xictau3D)[bestCand] - genctau );
        ctau3DPull->Fill( ((*xictau3D)[bestCand] - genctau)/(*xictau3DE)[bestCand] );
        ctauB->Fill( (*xictauB)[bestCand] );
        ctauBRes->Fill( (*xictauB)[bestCand] - genctau );
        ctauBPull->Fill( ((*xictauB)[bestCand] - genctau)/(*xictauBE)[bestCand] );
      } else {
        dataNotTruth->Fill( (*xiMass)[bestCand]*1000 );
	if(xiPt < 0.2)dataNotTruthLowPt->Fill( (*xiMass)[bestCand]*1000 );
	
      }
    } 
  if (nCands>0) nEvtWCands++;
  if (nCands==1) nSingleCands++;
  nTotalCands += nCands;
  nXiCands->Fill( nCands );
  }
   
  cout << "Total candidates = " << nTotalCands << " from " << nEvtWCands << " events and single cand events = " << nSingleCands << endl; 
   
  TCanvas c1;

  data->GetXaxis()->SetNdivisions(505);
  data->Draw();  
  int bins = 100;
  int x_min = 1255;
  int x_max = 1455;
  int x_range = x_max-x_min;
  int binwidth = x_range/bins;
  TString title = TString("#Lambda^{0} #pi^{-} invariant mass [MeV/c^{2}]");
   // build ARGUS + Gaussian
  TF1 *f1 = new TF1("f1", Form("%s%i%s%i%s%i%s%i%s","[0]*((x-", x_min, ")**(0.5)) +[1]*((x-", x_min, ")**(3/2)) + ((", x_range, ")/", bins, ")*gausn(2)"), x_min, x_max);

  f1->SetParameters(0.55, 0., 120, 1320, 60);
  f1->SetParName(0,"bkg param 1");
  f1->SetParName(1,"bkg param 2");
  f1->SetParName(2,"yield");
  f1->SetParLimits(2,0,100000);
  f1->SetParName(3,"mean");
  f1->SetParName(4,"sigma");
  f1->SetParLimits(4,0,1000);
  

  
  data->GetXaxis()->SetTitle(title);
  data->GetYaxis()->SetTitle(Form("%s%i%s","Candidates / ", binwidth, " MeV/c^{2}"));
  data->SetStats(0);



  f1->SetParameters(0.6, 0., 50, 1320, 40);

  data->Fit("f1","LEM0");
  f1->GetXaxis()->SetTitle(title);
  f1->GetYaxis()->SetTitle(Form("%s%i%s","Candidates / ", binwidth, " MeV/c^{2}"));
  f1->SetNpx(5000);
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
  TF1 *f30 = new TF1("f3", Form("%s%i%s%i%s","[0]*((x-", x_min, ")**(0.5)) +[1]*((x-", x_min, ")**(3/2))"), x_min, x_max);

  f30->SetParameter(0, bkgA);
  f30->SetParameter(1, bkgB);
  f30->SetLineStyle(2);
  f30->Draw("same");

  TLatex *cms = new TLatex();
  cms->SetNDC();
  cms->SetTextSize(0.040);
  cms->DrawLatex(0.60,0.89,"CMS Preliminary");
  
  std::string xiMassSelectedOut = theOutDir + "/XiMass_selected_MC_" + theSqrtS + ".eps";
  std::string xiMassSelectedOutOm = theOutDir + "/OmMass_selected_MC_" + theSqrtS + ".eps";
  c1.SaveAs(xiMassSelectedOut.c_str());
  
  dataOm->Draw();
  c1.SaveAs(xiMassSelectedOutOm.c_str());
  
  TCanvas* c2 = new TCanvas("c2","c2",900,600);
  c2->Divide(3,2);
  c2->cd(1);
  denEff->GetXaxis()->SetTitle("Generated c#tau (cm)");
  TF1 *exp1 = new TF1("exp1", Form("%s%f%s","[0]*exp(-x/([1]*", c, "))"), theCTauBins.front(), theCTauBins.back());
  exp1->SetParameters(50000, 164.0);
  exp1->SetParName(0,"Constant");
  exp1->SetParName(1,"tau [ps]");
  denEff->Fit(exp1,"LEMI");
  denEff->Draw();

  c2->cd(4)->SetLogy();
  denEff->Draw();

  c2->cd(2);
  numEff->GetXaxis()->SetTitle("Truth-matched reco c#tau (cm)");
  numEff->Draw();
   
  c2->cd(5)->SetLogy();
  numEff->Draw();
   
  c2->cd(3);
  for (int i = 1; i< nBins+1; ++i) {
    float eff2;
    if ( denEff->GetBinContent(i) != 0 )
      eff2 = numEff->GetBinContent(i)/denEff->GetBinContent(i);
    else eff2 = 0;
    effH->SetBinContent(i,eff2);
  }
  effH->GetYaxis()->SetTitle("Efficiency");
  effH->GetXaxis()->SetTitle("c#tau (cm)");
  effH->Draw();
  
  c2->cd(6)->SetLogy();
  
  for (int i = 1; i< theCTauBins.size(); ++i) {
    float width = theCTauBins.at(i)-theCTauBins.at(i-1);
    allGenCtau->SetBinContent(i, denEffCtauUncor->GetBinContent(i)/width);
  }
  
  TH1F *allGenCtau_clone = (TH1F*)allGenCtau->Clone("allGenCtau_clone");
  TF1 *exp2 = new TF1("exp2", Form("%s%f%s","[0]*exp(-x/([1]*", c, "))"), theCTauBins.front(), theCTauBins.back());
  exp2->SetParameters(200000, 164.0);
  exp2->SetParName(0,"Constant");
  exp2->SetParName(1,"tau [ps]");
  allGenCtau->Fit(exp2,"LEMI");
  allGenCtau->GetXaxis()->SetTitle("Generated c#tau (cm)");
  allGenCtau->Draw();
  
  std::string XiEffMCOutput = theOutDir + "/XiEff_MC_" + theSqrtS + ".eps";
  c2->SaveAs(XiEffMCOutput.c_str());
  
  TCanvas* theCTauCanvas = new TCanvas("theCTauCanvas", "theCTauBins"); 														        
  string theCTauBinsFileName = theOutDir + "/cTauBinsFit.eps";																        
  theCTauCanvas->Divide(3,4);																				        
  for(unsigned int vIt = 0; vIt < theCTauHistoVector.size(); ++vIt){															        
    theCTauCanvas->cd(vIt + 1); 																			        
    FitXi(theCTauHistoVector.at(vIt));															        
    numFitMCCtau->SetBinContent(vIt + 1, theFitYield);								        
    numFitMCCtau->SetBinError(vIt + 1, theFitYieldError);								        
    theCTauHistoVector.at(vIt)->Draw(); 																		        
  }																							        
  theCTauCanvas->Print(theCTauBinsFileName.c_str());																	        

  TCanvas* c3 = new TCanvas("c3","c3", 600,600);
  c3->Divide(2,2);

  c3->cd(1);
  
  for (int i = 0; i < theCTauBins.size() - 1; i++) {
    std::cout << "for bin " << i+1;
    std::cout << " numEffCtau      = " <<  numEffCtau->GetBinContent(i+1) << " and denEffCtau      = " << denEffCtau->GetBinContent(i+1) << std::endl;
    std::cout << " numEffCtauUncordata = " <<  numEffCtauUncor->GetBinContent(i+1) << " and denEffCtauUncor = " << denEffCtauUncor->GetBinContent(i+1) << std::endl;
  }
  effCtau->Divide(numEffCtau,denEffCtau);
  effCtauUncor->Divide(numEffCtauUncor,denEffCtauUncor);
  numFitMCCtau->SetMarkerColor(2);
  numFitMCCtau->GetXaxis()->SetTitle("#Xi c#tau (cm)");
  numFitMCCtau->GetYaxis()->SetTitle("Unscaled entries");
  numFitMCCtau->Draw();
  c3->cd(2);
  effCtau->SetMarkerColor(2);
  effCtau->SetMarkerSize(0.8);
  effCtau->SetMarkerStyle(20);
  effCtau->GetXaxis()->SetTitle("#Xi c#tau (cm)");
  effCtau->GetYaxis()->SetTitle("Efficiency");
  effCtau->SetMinimum(0);
  effCtau->Draw("P");
  effCtauUncor->SetMarkerColor(4);
  effCtauUncor->SetMarkerSize(0.8);
  effCtauUncor->SetMarkerStyle(24);
  effCtauUncor->Draw("same,P");
  TLegend* leg = new TLegend(0.50,0.74,0.88,0.93);
  leg->AddEntry(effCtau,"Reweighted","p");
  leg->AddEntry(effCtauUncor,"Unreweighted","p");
  leg->SetTextSize(0.055);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetShadowColor(0);
  leg->Draw("same");
  c3->cd(3);
  for (int i = 0; i < theCTauBins.size() - 1; i++) {
    float width = theCTauBins.at(i+1)-theCTauBins.at(i);
    widthCtau->SetBinContent(i+1, width);
    widthCtau->SetBinError(i+1,0);
  }
  // want histograms for everything and then to do (yield/(efficiency*width))
  numEffCtauScaled->Divide(numFitMCCtau, effCtauUncor);
  numEffCtauScaled->Divide(numEffCtauScaled, widthCtau);
  numEffCtauScaled->GetXaxis()->SetTitle("#Xi c#tau (cm)");   
  numEffCtauScaled->GetYaxis()->SetTitle("Scaled entries/cm");
  TF1 *exp3 = new TF1("exp3", Form("%s%f%s","[0]*exp(-x/([1]*", c, "))"), theCTauBins.front(), theCTauBins.back());
  exp3->SetParameters(10, 164.0);
  exp3->SetParName(0,"Constant");
  exp3->SetParName(1,"tau [ps]");
  numEffCtauScaled->Fit(exp3,"EMI");
  numEffCtauScaled->Draw();
  c3->cd(4)->SetLogy();
  //numEffCtauScaled->SetMinimum(0.1);
  numEffCtauScaled->Draw();
  allGenCtau_clone->Draw("same");
  
  std::string XiEffCtauOutput = theOutDir + "/XiEffCtau_MC_" + theSqrtS + ".eps";
  c3->SaveAs(XiEffCtauOutput.c_str());
  TCanvas* c4 = new TCanvas("c4","c4", 900, 600);
  c4->Divide(2,1);
  c4->cd(1);
  dataTruth->GetXaxis()->SetTitle("#Xi mass from truth-matched events");
  dataTruth->Draw();
  c4->cd(2);
  dataNotTruth->GetXaxis()->SetTitle("#Xi mass from non-truth-matched events");
  dataNotTruth->Draw();
  
  std::string XiTruthCompOutput = theOutDir + "/XiTruthComp_MC_" + theSqrtS + ".eps";
  c4->SaveAs(XiTruthCompOutput.c_str());  
  
  c4->cd(1);
  dataTruthLowPt->GetXaxis()->SetTitle("#Xi mass from truth-matched events Pt < 0.2");
  dataTruthLowPt->Draw();
  c4->cd(2);
  dataNotTruthLowPt->GetXaxis()->SetTitle("#Xi mass from non-truth-matched events Pt < 0.2");
  dataNotTruthLowPt->Draw();
  
  XiTruthCompOutput = theOutDir + "/XiTruthComp_MC_bin1_" + theSqrtS + ".eps";
  c4->SaveAs(XiTruthCompOutput.c_str());  

  
  TCanvas* c5 = new TCanvas("c5","c5", 900,900);
  c5->Divide(3,4);
  c5->cd(1);
  ctauMPV->GetXaxis()->SetTitle("c#tau MPV (cm)");
  ctauMPV->Draw();
  c5->cd(2);
  ctauMPVRes->GetXaxis()->SetTitle("c#tau MPV residual (cm)");
  ctauMPVRes->Draw();
  c5->cd(3);
  ctauMPVPull->GetXaxis()->SetTitle("c#tau MPV pull");
  ctauMPVPull->Fit("gaus");
  ctauMPVPull->Draw();

  c5->cd(4);
  ctauMPVAbs->GetXaxis()->SetTitle("c#tau MPV fabs(c#tau) (cm)");
  ctauMPVAbs->Draw();
  c5->cd(5);
  ctauMPVAbsRes->GetXaxis()->SetTitle("c#tau residual fabs(c#tau) from reco (cm)");
  ctauMPVAbsRes->Draw();
  c5->cd(6);
  ctauMPVAbsPull->GetXaxis()->SetTitle("c#tau MPV fabs(c#tau) pull");
  ctauMPVAbsPull->Fit("gaus");
  ctauMPVAbsPull->Draw();

  c5->cd(7);
  ctau3D->GetXaxis()->SetTitle("c#tau 3D (cm)");
  ctau3D->Draw();
  c5->cd(8);
  ctau3DRes->GetXaxis()->SetTitle("c#tau 3D residual (cm)");
  ctau3DRes->Draw();
  c5->cd(9);
  ctau3DPull->GetXaxis()->SetTitle("c#tau 3D pull");
  ctau3DPull->Fit("gaus");
  ctau3DPull->Draw();

  c5->cd(10);
  ctauB->GetXaxis()->SetTitle("c#tau L (cm)");
  ctauB->Draw();
  c5->cd(11);
  ctauBRes->GetXaxis()->SetTitle("c#tau L residual (cm)");
  ctauBRes->Draw();
  c5->cd(12);
  ctauBPull->GetXaxis()->SetTitle("c#tau L pull");
  ctauBPull->Fit("gaus");
  ctauBPull->Draw();
  std::string XiCtauResOutput = theOutDir + "/XiCtauRes_MC_" + theSqrtS + ".eps";
  c5->SaveAs(XiCtauResOutput.c_str());

//////////////////////////////////////////////
// make rapidity plots
  
  TCanvas* theYCanvas = new TCanvas("theYCanvas", "theYBins"); 														        
  string theYBinsFileName = theOutDir + "/YBinsFit.eps";																        
  theYCanvas->Divide(3,4);																				        
  for(unsigned int vIt = 0; vIt < theYHistoVector.size(); ++vIt){															        
    theYCanvas->cd(vIt + 1); 																			        
    FitXi(theYHistoVector.at(vIt));															        
    numFitMCY->SetBinContent(vIt + 1, theFitYield);								        
    numFitMCY->SetBinError(vIt + 1, theFitYieldError);								        
    theYHistoVector.at(vIt)->Draw(); 																		        
  }																							        
  theYCanvas->Print(theYBinsFileName.c_str());																	        

  TCanvas* cYRes = new TCanvas("cYRes","cYRes", 1200,400);
  cYRes->Divide(3,1);

  cYRes->cd(1);
  
  for (int i = 0; i<theYBins.size() -1; i++) {
    std::cout << "for bin " << i+1;
    std::cout << " numEffY      = " <<  numEffY->GetBinContent(i+1) << " and denEffY      = " << denEffY->GetBinContent(i+1) << std::endl;
    std::cout << " numEffYUncor = " <<  numEffYUncor->GetBinContent(i+1) << " and denEffYUncor = " << denEffYUncor->GetBinContent(i+1) << std::endl;
  }

  effY->Divide(numEffY,denEffY);
  effYUncor->Divide(numEffYUncor,denEffYUncor);
  
  numFitMCY->SetMarkerColor(2);
  numFitMCY->SetMarkerSize(0.8);
  numFitMCY->SetMarkerStyle(20);
  numFitMCY->GetXaxis()->SetTitle("#Xi Y");
  numFitMCY->GetYaxis()->SetTitle("Unscaled entries");
  numFitMCY->SetMinimum(0);
  numFitMCY->Draw("P");

  cYRes->cd(2);
   
  effY->SetMarkerColor(2);
  effY->SetMarkerSize(0.8);
  effY->SetMarkerStyle(20);
  effY->GetXaxis()->SetTitle("#Xi Y");
  effY->GetYaxis()->SetTitle("Efficiency");
  effY->SetMinimum(0);
  effY->Draw("P");

  effYUncor->SetMarkerColor(4);
  effYUncor->SetMarkerSize(0.8);
  effYUncor->SetMarkerStyle(24);
  effYUncor->Draw("same,P");
  
  leg = new TLegend(0.20,0.14,0.58,0.33);
  leg->AddEntry(effY,"Reweighted","p");
  leg->AddEntry(effYUncor,"Unreweighted","p");
  leg->SetTextSize(0.055);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetShadowColor(0);
  leg->Draw("same");

  cYRes->cd(3);
  
  for (int i = 0; i<theYBins.size() - 1; i++) {
    float width = (theYBins.at(i+1)-theYBins.at(i))*2; //factor of 2 to accounted for folded distribution
    widthY->SetBinContent(i+1, width);
    widthY->SetBinError(i+1,0);
  }

  // want histograms for everything and then to do (yield/(efficiency*width))
  numEffYScaled->Divide(numFitMCY, effYUncor);
  numEffYScaled->Divide(numEffYScaled, widthY);
  
  numEffYScaled->GetXaxis()->SetTitle("#Xi Y");   
  numEffYScaled->GetYaxis()->SetTitle("Scaled entries/unit");
  numEffYScaled->SetMinimum(0);
  numEffYScaled->Draw();
  // must scale by 1/bin width to report results /unit rapidity
  allGenY->Scale(5);  //should be 5 to account for folded distribution 
  allGenY->Draw("same");
  
  std::string XiEffYOutput = theOutDir + "/XiEffY_MC_" + theSqrtS + ".eps";
  cYRes->SaveAs(XiEffYOutput.c_str()); 



  TCanvas* thePosXiYCanvas = new TCanvas("thePosXiYCanvas", "theYBins"); 														        
  string thePosXiYBinsFileName = theOutDir + "/PosXiYBinsFit.eps";																        
  thePosXiYCanvas->Divide(3,4);																				        
  for(unsigned int vIt = 0; vIt < thePosYHistoVector.size(); ++vIt){															        
    thePosXiYCanvas->cd(vIt + 1); 																			        
    FitXi(thePosYHistoVector.at(vIt));															        
    numFitPosXiMCY->SetBinContent(vIt + 1, theFitYield);								        
    numFitPosXiMCY->SetBinError(vIt + 1, theFitYieldError);								        
    thePosYHistoVector.at(vIt)->Draw(); 																		        
  }																							        
  thePosXiYCanvas->Print(thePosXiYBinsFileName.c_str());																	        

  TCanvas* cPosXiYRes = new TCanvas("cPosXiYRes","cPosXiYRes", 1200,400);
  cPosXiYRes->Divide(3,1);

  cPosXiYRes->cd(1);
  
  for (int i = 0; i<theYBins.size(); i++) {
    std::cout << "for bin " << i+1;
    std::cout << " numEffPosXiY      = " <<  numEffPosXiY->GetBinContent(i+1) << " and denEffPosXiY      = " << denEffPosXiY->GetBinContent(i+1) << std::endl;
    std::cout << " numEffPosXiYUncor = " <<  numEffPosXiYUncor->GetBinContent(i+1) << " and denEffPosXiYUncor = " << denEffPosXiYUncor->GetBinContent(i+1) << std::endl;
  }

  effPosXiY->Divide(numEffPosXiY,denEffPosXiY);
  effPosXiYUncor->Divide(numEffPosXiYUncor,denEffPosXiYUncor);
  
  numFitPosXiMCY->SetMarkerColor(2);
  numFitPosXiMCY->SetMarkerSize(0.8);
  numFitPosXiMCY->SetMarkerStyle(20);
  numFitPosXiMCY->GetXaxis()->SetTitle("#Xi PosXiY");
  numFitPosXiMCY->GetYaxis()->SetTitle("Unscaled entries");
  numFitPosXiMCY->SetMinimum(0);
  numFitPosXiMCY->Draw("P");

  cPosXiYRes->cd(2);
   
  effPosXiY->SetMarkerColor(2);
  effPosXiY->SetMarkerSize(0.8);
  effPosXiY->SetMarkerStyle(20);
  effPosXiY->GetXaxis()->SetTitle("#Xi PosXiY");
  effPosXiY->GetYaxis()->SetTitle("Efficiency");
  effPosXiY->SetMinimum(0);
  effPosXiY->Draw("P");

  effPosXiYUncor->SetMarkerColor(4);
  effPosXiYUncor->SetMarkerSize(0.8);
  effPosXiYUncor->SetMarkerStyle(24);
  effPosXiYUncor->Draw("same,P");
  
  leg = new TLegend(0.20,0.14,0.58,0.33);
  leg->AddEntry(effPosXiY,"Reweighted","p");
  leg->AddEntry(effPosXiYUncor,"Unreweighted","p");
  leg->SetTextSize(0.055);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetShadowColor(0);
  leg->Draw("same");

  cPosXiYRes->cd(3);
  
  for (int i = 0; i<theYBins.size() - 1; i++) {
    float width = (theYBins.at(i+1)-theYBins.at(i))*2;
    widthPosXiY->SetBinContent(i+1, width);
    widthPosXiY->SetBinError(i+1,0);
  }

  // want histograms for everything and then to do (yield/(efficiency*width))
  numEffPosXiYScaled->Divide(numFitPosXiMCY, effPosXiYUncor);
  numEffPosXiYScaled->Divide(numEffPosXiYScaled, widthPosXiY);
  
  numEffPosXiYScaled->GetXaxis()->SetTitle("#Xi PosXiY");   
  numEffPosXiYScaled->GetYaxis()->SetTitle("Scaled entries/unit");
  numEffPosXiYScaled->SetMinimum(0);
  numEffPosXiYScaled->Draw();
  // must scale by 1/bin width to report results /unit rapidity
  allGenPosXiY->Scale(5);
  allGenPosXiY->Draw("same");
  
  std::string XiEffPosXiYOutput = theOutDir + "/XiEffPosXiY_MC_" + theSqrtS + ".eps";
  
  cPosXiYRes->SaveAs(XiEffPosXiYOutput.c_str()); 


  TCanvas* theNegXiYCanvas = new TCanvas("theNegXiYCanvas", "theNegXiYBins"); 														        
  string theNegXiYBinsFileName = theOutDir + "/NegXiYBinsFit.eps";																        
  theNegXiYCanvas->Divide(3,4);																				        
  for(unsigned int vIt = 0; vIt < theNegYHistoVector.size(); ++vIt){															        
    theNegXiYCanvas->cd(vIt + 1); 																			        
    FitXi(theNegYHistoVector.at(vIt));															        
    numFitNegXiMCY->SetBinContent(vIt + 1, theFitYield);								        
    numFitNegXiMCY->SetBinError(vIt + 1, theFitYieldError);								        
    theNegYHistoVector.at(vIt)->Draw(); 																		        
  }																							        
  theNegXiYCanvas->Print(theNegXiYBinsFileName.c_str());																	        

  TCanvas* cNegXiYRes = new TCanvas("cNegXiYRes","cNegXiYRes", 1200,400);
  cNegXiYRes->Divide(3,1);

  cNegXiYRes->cd(1);
  
  for (int i = 0; i<theYBins.size(); i++) {
    std::cout << "for bin " << i+1;
    std::cout << " numEffNegXiY      = " <<  numEffNegXiY->GetBinContent(i+1) << " and denEffNegXiY      = " << denEffNegXiY->GetBinContent(i+1) << std::endl;
    std::cout << " numEffNegXiYUncor = " <<  numEffNegXiYUncor->GetBinContent(i+1) << " and denEffNegXiYUncor = " << denEffNegXiYUncor->GetBinContent(i+1) << std::endl;
  }

  effNegXiY->Divide(numEffNegXiY,denEffNegXiY);
  effNegXiYUncor->Divide(numEffNegXiYUncor,denEffNegXiYUncor);
  
  numFitNegXiMCY->SetMarkerColor(2);
  numFitNegXiMCY->SetMarkerSize(0.8);
  numFitNegXiMCY->SetMarkerStyle(20);
  numFitNegXiMCY->GetXaxis()->SetTitle("#Xi NegXiY");
  numFitNegXiMCY->GetYaxis()->SetTitle("Unscaled entries");
  numFitNegXiMCY->SetMinimum(0);
  numFitNegXiMCY->Draw("P");

  cNegXiYRes->cd(2);
   
  effNegXiY->SetMarkerColor(2);
  effNegXiY->SetMarkerSize(0.8);
  effNegXiY->SetMarkerStyle(20);
  effNegXiY->GetXaxis()->SetTitle("#Xi NegXiY");
  effNegXiY->GetYaxis()->SetTitle("Efficiency");
  effNegXiY->SetMinimum(0);
  effNegXiY->Draw("P");

  effNegXiYUncor->SetMarkerColor(4);
  effNegXiYUncor->SetMarkerSize(0.8);
  effNegXiYUncor->SetMarkerStyle(24);
  effNegXiYUncor->Draw("same,P");
  
  leg = new TLegend(0.20,0.14,0.58,0.33);
  leg->AddEntry(effNegXiY,"Reweighted","p");
  leg->AddEntry(effNegXiYUncor,"Unreweighted","p");
  leg->SetTextSize(0.055);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetShadowColor(0);
  leg->Draw("same");

  cNegXiYRes->cd(3);
  
  for (int i = 0; i<theYBins.size() - 1; i++) {
    float width = (theYBins.at(i+1)-theYBins.at(i))*2;
    widthNegXiY->SetBinContent(i+1, width);
    widthNegXiY->SetBinError(i+1,0);
  }

  // want histograms for everything and then to do (yield/(efficiency*width))
  numEffNegXiYScaled->Divide(numFitNegXiMCY, effNegXiYUncor);
  numEffNegXiYScaled->Divide(numEffNegXiYScaled, widthNegXiY);
  
  numEffNegXiYScaled->GetXaxis()->SetTitle("#Xi NegXiY");   
  numEffNegXiYScaled->GetYaxis()->SetTitle("Scaled entries/unit");
  numEffNegXiYScaled->SetMinimum(0);
  numEffNegXiYScaled->Draw();
  // must scale by 1/bin width to report results /unit rapidity
  allGenNegXiY->Scale(5);
  allGenNegXiY->Draw("same");
  
  std::string XiEffNegXiYOutput = theOutDir + "/XiEffNegXiY_MC_" + theSqrtS + ".eps";
  cNegXiYRes->SaveAs(XiEffNegXiYOutput.c_str()); 

//////////////////////////////////////////////
// make pT plots 

  
  TCanvas* thePtCanvas = new TCanvas("thePtCanvas", "thePtBins"); 														        
  string thePtBinsFileName = theOutDir + "/PtBinsFit.eps";																        
  thePtCanvas->Divide(4,5);																				        
  for(unsigned int vIt = 0; vIt < thePtHistoVector.size(); ++vIt){															        
    thePtCanvas->cd(vIt + 1); 																			        
    if(vIt<2)
      FitXi(thePtHistoVector.at(vIt), true);	
    else
      FitXi(thePtHistoVector.at(vIt));															        
    numFitMCpT->SetBinContent(vIt + 1, theFitYield);								        
    numFitMCpT->SetBinError(vIt + 1, theFitYieldError);								        
    thePtHistoVector.at(vIt)->Draw(); 																		        
  }																							        
  thePtCanvas->Print(thePtBinsFileName.c_str());																	        

  TCanvas* cpTRes = new TCanvas("cpTRes","cpTRes", 900,900);
  cpTRes->Divide(2,2);

  cpTRes->cd(1);
  
  for (int i = 0; i<thePtBins.size() -1; i++) {
    std::cout << "for bin " << i+1;
    std::cout << " numEffpT      = " <<  numEffpT->GetBinContent(i+1) << " and denEffpT      = " << denEffpT->GetBinContent(i+1) << std::endl;
    std::cout << " numEffpTUncor = " <<  numEffpTUncor->GetBinContent(i+1) << " and denEffpTUncor = " << denEffpTUncor->GetBinContent(i+1) << std::endl;
  }

  effpT->Divide(numEffpT,denEffpT);
  effpTUncor->Divide(numEffpTUncor,denEffpTUncor);
  
  numFitMCpT->SetMarkerColor(2);
  numFitMCpT->SetMarkerSize(0.8);
  numFitMCpT->SetMarkerStyle(20);
  numFitMCpT->GetXaxis()->SetTitle("#Xi pT (GeV)");
  numFitMCpT->GetYaxis()->SetTitle("Unscaled entries");
  numFitMCpT->SetMinimum(0);
  numFitMCpT->Draw("P");

  cpTRes->cd(2);
   
  effpT->SetMarkerColor(2);
  effpT->SetMarkerSize(0.8);
  effpT->SetMarkerStyle(20);
  effpT->GetXaxis()->SetTitle("#Xi pT (GeV)");
  effpT->GetYaxis()->SetTitle("Efficiency");
  effpT->SetMinimum(0);
  effpT->Draw("P");

  effpTUncor->SetMarkerColor(4);
  effpTUncor->SetMarkerSize(0.8);
  effpTUncor->SetMarkerStyle(24);
  effpTUncor->Draw("same,P");
  
  leg = new TLegend(0.50,0.14,0.88,0.33);
  leg->AddEntry(effpT,"Reweighted","p");
  leg->AddEntry(effpTUncor,"Unreweighted","p");
  leg->SetTextSize(0.055);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetShadowColor(0);
  leg->Draw("same");

  cpTRes->cd(3);

  for (int i = 0; i<thePtBins.size() - 1; i++) {
    float width = thePtBins.at(i+1)-thePtBins.at(i);
    widthpT->SetBinContent(i+1, width);
    widthpT->SetBinError(i+1,0);
  }

  // want histograms for everything and then to do (yield/(efficiency*width))
  numEffpTScaled->Divide(numFitMCpT, effpTUncor);
  numEffpTScaled->Divide(numEffpTScaled, widthpT);  
  
  numEffpTScaled->GetXaxis()->SetTitle("#Xi pT (GeV)");   
  numEffpTScaled->GetYaxis()->SetTitle("Scaled entries/GeV");
  numEffpTScaled->Draw();
  // must scale by 1/bin width to report results per unit GeV
  allGenpT->Scale(10);
  allGenpT->Draw("same");

  cpTRes->cd(4)->SetLogy();
  numEffpTScaled->Draw();
  allGenpT->Draw("same");
  
  std::string XiEffpTOutput = theOutDir + "/XiEffpT_MC_" + theSqrtS + ".eps";
  cpTRes->SaveAs(XiEffpTOutput.c_str());
  
  TCanvas* thePosXiPtCanvas = new TCanvas("thePosXiPtCanvas", "thePosXiPtBins"); 														        
  string thePosXiPtBinsFileName = theOutDir + "/PosXiPtBinsFit.eps";																        
  thePosXiPtCanvas->Divide(4,5);																				        
  for(unsigned int vIt = 0; vIt < thePosPtHistoVector.size(); ++vIt){															        
    thePosXiPtCanvas->cd(vIt + 1); 	
    if(vIt<2)
      FitXi(thePosPtHistoVector.at(vIt), true);	
    else
      FitXi(thePosPtHistoVector.at(vIt));																			        														        
    numFitPosXiMCpT->SetBinContent(vIt + 1, theFitYield);								        
    numFitPosXiMCpT->SetBinError(vIt + 1, theFitYieldError);								        
    thePosPtHistoVector.at(vIt)->Draw(); 																		        
  }																							        
  thePosXiPtCanvas->Print(thePosXiPtBinsFileName.c_str());																	        

  TCanvas* cPosXipTRes = new TCanvas("cpTRes","cpTRes", 900,900);
  cPosXipTRes->Divide(2,2);

  cPosXipTRes->cd(1);
  
  for (int i = 0; i<thePtBins.size() -1; i++) {
    std::cout << "for bin " << i+1;
    std::cout << " numEffPosXipT      = " <<  numEffPosXipT->GetBinContent(i+1) << " and denEffPosXipT      = " << denEffPosXipT->GetBinContent(i+1) << std::endl;
    std::cout << " numEffPosXipTUncor = " <<  numEffPosXipTUncor->GetBinContent(i+1) << " and denEffPosXipTUncor = " << denEffPosXipTUncor->GetBinContent(i+1) << std::endl;
  }

  effPosXipT->Divide(numEffPosXipT,denEffPosXipT);
  effPosXipTUncor->Divide(numEffPosXipTUncor,denEffPosXipTUncor);
  
  numFitPosXiMCpT->SetMarkerColor(2);
  numFitPosXiMCpT->SetMarkerSize(0.8);
  numFitPosXiMCpT->SetMarkerStyle(20);
  numFitPosXiMCpT->GetXaxis()->SetTitle("#Xi PosXipT (GeV)");
  numFitPosXiMCpT->GetYaxis()->SetTitle("Unscaled entries");
  numFitPosXiMCpT->SetMinimum(0);
  numFitPosXiMCpT->Draw("P");

  cPosXipTRes->cd(2);
   
  effPosXipT->SetMarkerColor(2);
  effPosXipT->SetMarkerSize(0.8);
  effPosXipT->SetMarkerStyle(20);
  effPosXipT->GetXaxis()->SetTitle("#Xi PosXipT (GeV)");
  effPosXipT->GetYaxis()->SetTitle("Efficiency");
  effPosXipT->SetMinimum(0);
  effPosXipT->Draw("P");

  effPosXipTUncor->SetMarkerColor(4);
  effPosXipTUncor->SetMarkerSize(0.8);
  effPosXipTUncor->SetMarkerStyle(24);
  effPosXipTUncor->Draw("same,P");
  
  leg = new TLegend(0.50,0.14,0.88,0.33);
  leg->AddEntry(effPosXipT,"Reweighted","p");
  leg->AddEntry(effPosXipTUncor,"Unreweighted","p");
  leg->SetTextSize(0.055);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetShadowColor(0);
  leg->Draw("same");

  cPosXipTRes->cd(3);

  for (int i = 0; i<thePtBins.size() - 1; i++) {
    float width = thePtBins.at(i+1)-thePtBins.at(i);
    widthPosXipT->SetBinContent(i+1, width);
    widthPosXipT->SetBinError(i+1,0);
  }

  // want histograms for everything and then to do (yield/(efficiency*width))
  numEffPosXipTScaled->Divide(numFitPosXiMCpT, effPosXipTUncor);
  numEffPosXipTScaled->Divide(numEffPosXipTScaled, widthPosXipT);  
  
  numEffPosXipTScaled->GetXaxis()->SetTitle("#Xi PosXipT (GeV)");   
  numEffPosXipTScaled->GetYaxis()->SetTitle("Scaled entries/GeV");
  numEffPosXipTScaled->Draw();
  // must scale by 1/bin width to report results per unit GeV
  allGenPosXipT->Scale(10);
  allGenPosXipT->Draw("same");

  cPosXipTRes->cd(4)->SetLogy();
  numEffPosXipTScaled->Draw();
  allGenPosXipT->Draw("same");
  
  std::string XiEffPosXipTOutput = theOutDir + "/XiEffPosXipT_MC_" + theSqrtS + ".eps";
  cPosXipTRes->SaveAs(XiEffPosXipTOutput.c_str());

  TCanvas* theNegXiPtCanvas = new TCanvas("theNegXiPtCanvas", "theNegXiPtBins"); 														        
  string theNegXiPtBinsFileName = theOutDir + "/NegXiPtBinsFit.eps";																        
  theNegXiPtCanvas->Divide(4,5);																				        
  for(unsigned int vIt = 0; vIt < theNegPtHistoVector.size(); ++vIt){															        
    theNegXiPtCanvas->cd(vIt + 1); 																			        
    if(vIt<2)
      FitXi(theNegPtHistoVector.at(vIt), true);	
    else
      FitXi(theNegPtHistoVector.at(vIt));	
    numFitNegXiMCpT->SetBinContent(vIt + 1, theFitYield);								        
    numFitNegXiMCpT->SetBinError(vIt + 1, theFitYieldError);								        
    theNegPtHistoVector.at(vIt)->Draw(); 																		        
  }																							        
  theNegXiPtCanvas->Print(theNegXiPtBinsFileName.c_str());																	        

  TCanvas* cNegXipTRes = new TCanvas("cpTRes","cpTRes", 900,900);
  cNegXipTRes->Divide(2,2);

  cNegXipTRes->cd(1);
  
  for (int i = 0; i<thePtBins.size() -1; i++) {
    std::cout << "for bin " << i+1;
    std::cout << " numEffNegXipT      = " <<  numEffNegXipT->GetBinContent(i+1) << " and denEffNegXipT      = " << denEffNegXipT->GetBinContent(i+1) << std::endl;
    std::cout << " numEffNegXipTUncor = " <<  numEffNegXipTUncor->GetBinContent(i+1) << " and denEffNegXipTUncor = " << denEffNegXipTUncor->GetBinContent(i+1) << std::endl;
  }

  effNegXipT->Divide(numEffNegXipT,denEffNegXipT);
  effNegXipTUncor->Divide(numEffNegXipTUncor,denEffNegXipTUncor);
  
  numFitNegXiMCpT->SetMarkerColor(2);
  numFitNegXiMCpT->SetMarkerSize(0.8);
  numFitNegXiMCpT->SetMarkerStyle(20);
  numFitNegXiMCpT->GetXaxis()->SetTitle("#Xi NegXipT (GeV)");
  numFitNegXiMCpT->GetYaxis()->SetTitle("Unscaled entries");
  numFitNegXiMCpT->SetMinimum(0);
  numFitNegXiMCpT->Draw("P");

  cNegXipTRes->cd(2);
   
  effNegXipT->SetMarkerColor(2);
  effNegXipT->SetMarkerSize(0.8);
  effNegXipT->SetMarkerStyle(20);
  effNegXipT->GetXaxis()->SetTitle("#Xi NegXipT (GeV)");
  effNegXipT->GetYaxis()->SetTitle("Efficiency");
  effNegXipT->SetMinimum(0);
  effNegXipT->Draw("P");

  effNegXipTUncor->SetMarkerColor(4);
  effNegXipTUncor->SetMarkerSize(0.8);
  effNegXipTUncor->SetMarkerStyle(24);
  effNegXipTUncor->Draw("same,P");
  
  leg = new TLegend(0.50,0.14,0.88,0.33);
  leg->AddEntry(effNegXipT,"Reweighted","p");
  leg->AddEntry(effNegXipTUncor,"Unreweighted","p");
  leg->SetTextSize(0.055);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetShadowColor(0);
  leg->Draw("same");

  cNegXipTRes->cd(3);

  for (int i = 0; i<thePtBins.size() - 1; i++) {
    float width = thePtBins.at(i+1)-thePtBins.at(i);
    widthNegXipT->SetBinContent(i+1, width);
    widthNegXipT->SetBinError(i+1,0);
  }

  // want histograms for everything and then to do (yield/(efficiency*width))
  numEffNegXipTScaled->Divide(numFitNegXiMCpT, effNegXipTUncor);
  numEffNegXipTScaled->Divide(numEffNegXipTScaled, widthNegXipT);  
  
  numEffNegXipTScaled->GetXaxis()->SetTitle("#Xi NegXipT (GeV)");   
  numEffNegXipTScaled->GetYaxis()->SetTitle("Scaled entries/GeV");
  numEffNegXipTScaled->Draw();
  // must scale by 1/bin width to report results per unit GeV
  allGenNegXipT->Scale(10);
  allGenNegXipT->Draw("same");

  cNegXipTRes->cd(4)->SetLogy();
  numEffNegXipTScaled->Draw();
  allGenNegXipT->Draw("same");
  
  std::string XiEffNegXipTOutput = theOutDir + "/XiEffNegXipT_MC_" + theSqrtS + ".eps";
  cNegXipTRes->SaveAs(XiEffNegXipTOutput.c_str());

}

void seleXiMC::SetBins(){
   
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

void seleXiMC::BookHistos(){
  std::string theOutRootFileName = theOutDir + "/efficienciesXi_" + theSqrtS + ".root";
  theOutRootFile = new TFile(theOutRootFileName.c_str(), "Recreate"); 
  theOutRootFile->cd(); 	

  data = new TH1F("data", "data", 100, 1255, 1455);
  dataOm = new TH1F("dataOm", "dataOm", 100, 1455, 1855);
  dataTruth = new TH1F("dataTruth", "dataTruth", 100, 1255, 1455);
  dataNotTruth = new TH1F("dataNotTruth", "dataNotTruth", 100, 1255, 1455);
  
  dataTruthLowPt = new TH1F("dataTruthLowPt", "dataTruthLowPt", 100, 1255, 1455);
  dataNotTruthLowPt = new TH1F("dataNotTruthLowPt", "dataNotTruthLowPt", 100, 1255, 1455);
  
  ctauMPV = new TH1F("ctauMPV", "ctauMPV", 108, -2.0, 25.0);
  ctauMPVRes = new TH1F("ctauMPVRes", "ctauMPVRes", 100, -5.0, 5.0);
  ctauMPVPull = new TH1F("ctauMPVPull", "ctauMPVPull", 100, -5.0, 5.0);
  ctauMPVAbs = new TH1F("ctauMPVAbs", "ctauMPVAbs", 108, -2.0, 25.0);
  ctauMPVAbsRes = new TH1F("ctauMPVAbsRes", "ctauMPVAbsRes", 100, -5.0, 5.0);
  ctauMPVAbsPull = new TH1F("ctauMPVAbsPull", "ctauMPVAbsPull", 100, -5.0, 5.0);
  ctau3D = new TH1F("ctau3D", "ctau3D", 108, -2.0, 25.0);
  ctau3DRes = new TH1F("ctau3DRes", "ctau3DRes", 100, -5.0, 5.0);
  ctau3DPull = new TH1F("ctau3DPull", "ctau3DPull", 100, -5.0, 5.0);
  ctauB = new TH1F("ctauB", "ctauB", 108, -2.0, 25.0);
  ctauBRes = new TH1F("ctauBRes", "ctauBRes", 100, -5.0, 5.0);
  ctauBPull = new TH1F("ctauBPull", "ctauBPull", 100, -5.0, 5.0);

  allGenCtau = new TH1F("allGenCtau","allGenCtau",theCTauBins.size() - 1, &theCTauBins[0]);
  numEffCtau = new TH1F("numEffCtau","numEffCtau",theCTauBins.size() - 1, &theCTauBins[0]);
  numEffCtau->Sumw2();
  denEffCtau = new TH1F("denEffCtau","denEffCtau",theCTauBins.size() - 1, &theCTauBins[0]);
  denEffCtau->Sumw2();
  effCtau = new TH1F("effCtau","effCtau",theCTauBins.size() - 1, &theCTauBins[0]);
  effCtauUncor = new TH1F("effCtauUncor","effCtauUncor",theCTauBins.size() - 1, &theCTauBins[0]);
  widthCtau = new TH1F("widthCtau", "widthCtau",theCTauBins.size() - 1, &theCTauBins[0]);
  numEffCtauUncor = new TH1F("numEffCtauUncor","numEffCtauUncor",theCTauBins.size() - 1, &theCTauBins[0]);
  numEffCtauUncor->Sumw2();
  denEffCtauUncor = new TH1F("denEffCtauUncor","denEffCtauUncor",theCTauBins.size() - 1, &theCTauBins[0]);
  denEffCtauUncor->Sumw2();
  numEffCtauScaled = new TH1F("numEffCtauScaled","numEffCtauScaled",theCTauBins.size() - 1, &theCTauBins[0]);
  numFitMCCtau = new TH1F("numFitMCCtau","numFitMCCtau",theCTauBins.size() - 1, &theCTauBins[0]);   
  
  allGenY = new TH1F("allGenY","allGenY",25,0,2.5);
  numEffY = new TH1F("numEffY","numEffY", theYBins.size() - 1, &theYBins[0] );
  numEffY->Sumw2();
  denEffY = new TH1F("denEffY","denEffY", theYBins.size() - 1, &theYBins[0] );
  denEffY->Sumw2();
  effY = new TH1F("effY","effY", theYBins.size() - 1, &theYBins[0] );
  effYUncor = new TH1F("effYUncor","effYUncor", theYBins.size() - 1, &theYBins[0] );
  widthY = new TH1F("widthY", "widthY", theYBins.size() - 1, &theYBins[0] );
  numEffYUncor = new TH1F("numEffYUncor","numEffYUncor", theYBins.size() - 1, &theYBins[0] );
  numEffYUncor->Sumw2();
  denEffYUncor = new TH1F("denEffYUncor","denEffYUncor", theYBins.size() - 1, &theYBins[0] );
  denEffYUncor->Sumw2();
  numEffYScaled = new TH1F("numEffYScaled","numEffYScaled", theYBins.size() - 1, &theYBins[0] );
  numFitMCY = new TH1F("numFitMCY","numFitMCY", theYBins.size() - 1, &theYBins[0] );
  
  allGenpT = new TH1F("allGenpT","allGenpT",100,0,10);
  numEffpT = new TH1F("numEffpT","numEffpT", thePtBins.size() - 1, &thePtBins[0] );
  numEffpT->Sumw2();
  effpT = new TH1F("effpT","effpT", thePtBins.size() - 1, &thePtBins[0] );
  effpTUncor = new TH1F("effpTUncor","effpTUncor", thePtBins.size() - 1, &thePtBins[0] );
  denEffpT = new TH1F("denEffpT","denEffpT", thePtBins.size() - 1, &thePtBins[0] );
  denEffpT->Sumw2();
  widthpT = new TH1F("widthpT", "widthpT", thePtBins.size() - 1, &thePtBins[0] );
  numEffpTUncor = new TH1F("numEffpTUncor","numEffpTUncor", thePtBins.size() - 1, &thePtBins[0] );
  numEffpTUncor->Sumw2();
  denEffpTUncor = new TH1F("denEffpTUncor","denEffpTUncor", thePtBins.size() - 1, &thePtBins[0] );
  denEffpTUncor->Sumw2();
  numEffpTScaled = new TH1F("numEffpTScaled","numEffpTScaled", thePtBins.size() - 1, &thePtBins[0] );
  numFitMCpT = new TH1F("numFitMCpT","numFitMCpT", thePtBins.size() - 1, &thePtBins[0] );      
  
  allGenPosXipT = new TH1F("allGenPosXipT","allGenPosXipT",100,0,10);
  numEffPosXipT = new TH1F("numEffPosXipT","numEffPosXipT", thePtBins.size() - 1, &thePtBins[0] );
  numEffPosXipT->Sumw2();
  effPosXipT = new TH1F("effPosXipT","effPosXipT", thePtBins.size() - 1, &thePtBins[0] );
  effPosXipTUncor = new TH1F("effPosXipTUncor","effPosXipTUncor", thePtBins.size() - 1, &thePtBins[0] );
  denEffPosXipT = new TH1F("denEffPosXipT","denEffPosXipT", thePtBins.size() - 1, &thePtBins[0] );
  denEffPosXipT->Sumw2();
  widthPosXipT = new TH1F("widthPosXipT", "widthPosXipT", thePtBins.size() - 1, &thePtBins[0] );
  numEffPosXipTUncor = new TH1F("numEffPosXipTUncor","numEffPosXipTUncor", thePtBins.size() - 1, &thePtBins[0] );
  numEffPosXipTUncor->Sumw2();
  denEffPosXipTUncor = new TH1F("denEffPosXipTUncor","denEffPosXipTUncor", thePtBins.size() - 1, &thePtBins[0] );
  denEffPosXipTUncor->Sumw2();
  numEffPosXipTScaled = new TH1F("numEffPosXipTScaled","numEffPosXipTScaled", thePtBins.size() - 1, &thePtBins[0] );
  numFitPosXiMCpT = new TH1F("numFitPosXiMCpT","numFitPosXiMCpT", thePtBins.size() - 1, &thePtBins[0] );      

  allGenNegXipT = new TH1F("allGenNegXipT","allGenNegXipT",100,0,10);
  numEffNegXipT = new TH1F("numEffNegXipT","numEffNegXipT", thePtBins.size() - 1, &thePtBins[0] );
  numEffNegXipT->Sumw2();
  effNegXipT = new TH1F("effNegXipT","effNegXipT", thePtBins.size() - 1, &thePtBins[0] );
  effNegXipTUncor = new TH1F("effNegXipTUncor","effNegXipTUncor", thePtBins.size() - 1, &thePtBins[0] );
  denEffNegXipT = new TH1F("denEffNegXipT","denEffNegXipT", thePtBins.size() - 1, &thePtBins[0] );
  denEffNegXipT->Sumw2();
  widthNegXipT = new TH1F("widthNegXipT", "widthNegXipT", thePtBins.size() - 1, &thePtBins[0] );
  numEffNegXipTUncor = new TH1F("numEffNegXipTUncor","numEffNegXipTUncor", thePtBins.size() - 1, &thePtBins[0] );
  numEffNegXipTUncor->Sumw2();
  denEffNegXipTUncor = new TH1F("denEffNegXipTUncor","denEffNegXipTUncor", thePtBins.size() - 1, &thePtBins[0] );
  denEffNegXipTUncor->Sumw2();
  numEffNegXipTScaled = new TH1F("numEffNegXipTScaled","numEffNegXipTScaled", thePtBins.size() - 1, &thePtBins[0] );
  numFitNegXiMCpT = new TH1F("numFitNegXiMCpT","numFitNegXiMCpT", thePtBins.size() - 1, &thePtBins[0] );      
  
  allGenPosXiY = new TH1F("allGenPosXiY","allGenPosXiY",100,0,2.5);
  numEffPosXiY = new TH1F("numEffPosXiY","numEffPosXiY", theYBins.size() - 1, &theYBins[0] );
  numEffPosXiY->Sumw2();
  effPosXiY = new TH1F("effPosXiY","effPosXiY", theYBins.size() - 1, &theYBins[0] );
  effPosXiYUncor = new TH1F("effPosXiYUncor","effPosXiYUncor", theYBins.size() - 1, &theYBins[0] );
  denEffPosXiY = new TH1F("denEffPosXiY","denEffPosXiY", theYBins.size() - 1, &theYBins[0] );
  denEffPosXiY->Sumw2();
  widthPosXiY = new TH1F("widthPosXiY", "widthPosXiY", theYBins.size() - 1, &theYBins[0] );
  numEffPosXiYUncor = new TH1F("numEffPosXiYUncor","numEffPosXiYUncor", theYBins.size() - 1, &theYBins[0] );
  numEffPosXiYUncor->Sumw2();
  denEffPosXiYUncor = new TH1F("denEffPosXiYUncor","denEffPosXiYUncor", theYBins.size() - 1, &theYBins[0] );
  denEffPosXiYUncor->Sumw2();
  numEffPosXiYScaled = new TH1F("numEffPosXiYScaled","numEffPosXiYScaled", theYBins.size() - 1, &theYBins[0] );
  numFitPosXiMCY = new TH1F("numFitPosXiMCY","numFitPosXiMCY", theYBins.size() - 1, &theYBins[0] );      

  allGenNegXiY = new TH1F("allGenNegXiY","allGenNegXiY",100,0,2.5);
  numEffNegXiY = new TH1F("numEffNegXiY","numEffNegXiY", theYBins.size() - 1, &theYBins[0] );
  numEffNegXiY->Sumw2();
  effNegXiY = new TH1F("effNegXiY","effNegXiY", theYBins.size() - 1, &theYBins[0] );
  effNegXiYUncor = new TH1F("effNegXiYUncor","effNegXiYUncor", theYBins.size() - 1, &theYBins[0] );
  denEffNegXiY = new TH1F("denEffNegXiY","denEffNegXiY", theYBins.size() - 1, &theYBins[0] );
  denEffNegXiY->Sumw2();
  widthNegXiY = new TH1F("widthNegXiY", "widthNegXiY", theYBins.size() - 1, &theYBins[0] );
  numEffNegXiYUncor = new TH1F("numEffNegXiYUncor","numEffNegXiYUncor", theYBins.size() - 1, &theYBins[0] );
  numEffNegXiYUncor->Sumw2();
  denEffNegXiYUncor = new TH1F("denEffNegXiYUncor","denEffNegXiYUncor", theYBins.size() - 1, &theYBins[0] );
  denEffNegXiYUncor->Sumw2();
  numEffNegXiYScaled = new TH1F("numEffNegXiYScaled","numEffNegXiYScaled", theYBins.size() - 1, &theYBins[0] );
  numFitNegXiMCY = new TH1F("numFitNegXiMCY","numFitNegXiMCY", theYBins.size() - 1, &theYBins[0] );      

  nXiCands = new TH1F("nXiCands","nXiCands", 20, 0, 20);

  
  for(unsigned int hIt = 0; hIt < thePtBins.size() - 1; ++hIt){
    char* theHistoName = Form("ptBin%d",hIt+1);
    thePtHistoVector.push_back(new TH1F(theHistoName, theHistoName, 100, 1255, 1455));
    
    char* thePosHistoName = Form("posPtBin%d",hIt+1);
    thePosPtHistoVector.push_back(new TH1F(thePosHistoName, thePosHistoName, 100, 1255, 1455));
    
    char* theNegHistoName = Form("negPtBin%d",hIt+1);
    theNegPtHistoVector.push_back(new TH1F(theNegHistoName, theNegHistoName, 100, 1255, 1455));
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
}

void seleXiMC::FitXi(TH1F* theXiHisto, bool restrictRange){
    //initialize results to 0
  theFitYield = 0.;
  theFitYieldError = 0.;

  
  TAxis* theAxis = theXiHisto->GetXaxis();
  int x_min = (int)theAxis->GetXmin();
  int x_max = (int)theAxis->GetXmax();
  int bins = theAxis->GetNbins();
  int x_range = x_max - x_min;
  TF1* xiFit = new TF1("xiFit", Form("%s%i%s%i%s%i%s","[0]*((x-", x_min, ")**([1])) + ((", x_range, ")/", bins, ")*gausn(2)"), x_min, x_max);
  xiFit->SetParameters(1.5, 0.1, theXiHisto->GetEntries()*.5, 1321, 3);
  xiFit->SetParName(0,"bkg param 1");
  xiFit->SetParName(1,"bkg param 2");
  xiFit->SetParName(2,"yield");
  xiFit->SetParLimits(2,0,20000);
  xiFit->SetParName(3,"mean");
  xiFit->SetParLimits(3,1310,1330);
  xiFit->SetParName(4,"sigma");
  xiFit->SetParLimits(4,1,6);
  if(restrictRange) {
    xiFit->SetParameters(20.0, -0.5, theXiHisto->GetEntries()*.5, 1321, 3);
    theXiHisto->Fit("xiFit", "LEM","",1290,x_max);
  } else
    theXiHisto->Fit("xiFit", "LEM");
  
  theFitYield = xiFit->GetParameter(2);
  theFitYieldError = xiFit->GetParError(2);
  
}


double seleXiMC::ntrk_wgt(int ntrk, int weight_choice)
{

  if( weight_choice == 1 ) { // Weights for tuneD6T MC
    const int nbins = 75;
    const Double_t xAxis[nbins+1] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 84, 88, 92, 96, 100, 110, 120, 140, 160, 200};
    const Double_t ntrk_weight[nbins] = {0.684823691845, 0.567588746548, 0.572231352329, 0.577732563019, 0.612788677216, 0.665794491768, 0.74260365963, 0.822943210602, 0.895672500134, 0.94077038765, 0.964128136635, 0.968468666077, 0.948569893837, 0.930590629578, 0.909639060497, 0.892322838306, 0.875799536705, 0.865711688995, 0.862775504589, 0.864894151688, 0.872865557671, 0.879333734512, 0.894657731056, 0.912622451782, 0.928457558155, 0.950745940208, 0.967661619186, 0.990352988243, 1.01151168346, 1.03049361706, 1.05319392681, 1.07638561726, 1.09800362587, 1.12640702724, 1.14030277729, 1.17644119263, 1.19989466667, 1.2310706377, 1.26455545425, 1.29289066792, 1.31875061989, 1.35096299648, 1.39227795601, 1.4200656414, 1.45285916328, 1.48266518116, 1.52858412266, 1.5482827425, 1.57303643227, 1.62603139877, 1.66228163242, 1.73159658909, 1.79345393181, 1.83892428875, 1.8832205534, 1.93724298477, 1.97469103336, 2.00871300697, 2.03477692604, 2.06109476089, 2.08581829071, 2.08820271492, 2.09942793846, 2.10594010353, 2.13290381432, 2.12786817551, 2.15150284767, 2.18847417831, 2.23482131958, 2.29027318954, 2.51525402069, 3.07639288902, 4.86200761795, 12.9713077545, 64.1660614014};
    int ndx;
    for( ndx = 1; ndx < nbins+1; ndx++ ) {
      if( ntrk < xAxis[ndx] ) break;
    }
    ndx += -1;
    if( ndx >= nbins ) ndx = nbins-1;
    return ntrk_weight[ndx];
  }
  else if( weight_choice == 11 ) {// Weights for tuneD6T_900 MC
    const int nbins = 75;
    const Double_t xAxis[nbins+1] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 84, 88, 92, 96, 100, 110, 120, 140, 160, 200};
    const Double_t ntrk_weight[nbins] = {0.807127356529, 0.659041583538, 0.698155522346, 0.713301122189, 0.755888938904, 0.81724011898, 0.890749454498, 0.955834448338, 1.01427793503, 1.03914952278, 1.02860558033, 1.00652229786, 0.97242128849, 0.938611805439, 0.904110908508, 0.877858459949, 0.855360090733, 0.845526218414, 0.851646363735, 0.852675378323, 0.870900571346, 0.895979106426, 0.912231385708, 0.943045377731, 0.978449225426, 0.998540341854, 1.0323805809, 1.06876373291, 1.09394979477, 1.13242733479, 1.16942465305, 1.20131146908, 1.22183001041, 1.26455259323, 1.29725158215, 1.33382105827, 1.37354123592, 1.4025156498, 1.43412911892, 1.49138796329, 1.52684915066, 1.56775450706, 1.5998250246, 1.64751839638, 1.70411026478, 1.72870481014, 1.79020178318, 1.86345207691, 1.90034878254, 2.02073407173, 2.0822160244, 2.26331019402, 2.48342299461, 2.76777458191, 2.95888948441, 3.3442966938, 3.88835382462, 4.46009063721, 5.05966234207, 5.81398391724, 6.88521909714, 8.14224147797, 10.2224674225, 12.6465578079, 14.1106805801, 18.7609806061, 32.9513435364, 46.8789825439, 75.7474136353, 119.637054443, 200.0, 200.0, 200.0, 200.0, 200.0};
    int ndx;
    for( ndx = 1; ndx < nbins+1; ndx++ ) {
      if( ntrk < xAxis[ndx] ) break;
    }
    ndx += -1;
    if( ndx >= nbins ) ndx = nbins-1;
    return ntrk_weight[ndx];
  }
  return 1.;
}


double seleXiMC::xipt_wgt( double pt, int weight_choice ) {

  if( !weight_choice ) return 1.;

  else if( weight_choice == 1 ) {//D6T - updated for 362
    const int nbins = 22;
    const Double_t xAxis1[nbins+1] = {0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 5, 6};
    const Double_t pt_weight[nbins] = {1.644238, 1.157256, 1.104384, 1.275661, 1.285989, 1.429563, 1.303506, 1.265767, 1.03562, 0.9872718, 0.8572736, 0.8666133, 0.7739058, 0.6379386, 0.5731549, 0.5783848, 0.5446806, 0.5363548, 0.4679532, 0.4986373, 0.4317229, 0.4418376};
    int ndx;
    for( ndx = 1; ndx < nbins+1; ndx++ ) {
      if( pt < xAxis1[ndx] ) break;
    }
    ndx += -1;
    if ( ndx >= nbins ) ndx = nbins-1;
    return pt_weight[ndx];
  }

  else if( weight_choice == 11 ) {//900 GeV D6T not updated with new bins!
    const int nbins = 22;
    const Double_t xAxis1[nbins+1] = {0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 5, 6};
    const Double_t pt_weight[nbins] = {1.167261, 0.7844442, 0.8747076, 0.9220581, 1.05339, 1.077754, 1.231445, 1.005943, 1.076825, 1.050165, 1.207763, 0.9013405, 0.937835, 1.061079, 0.9664335, 0.8140177, 0.7338285, 0.8254001, 0.6236658, 1.150751, 0.7896331, 0.6822104} ;
    int ndx;
    for( ndx = 1; ndx < nbins+1; ndx++ ) {
      if( pt < xAxis1[ndx] ) break;
    }
    ndx += -1;
    if ( ndx >= nbins ) ndx = nbins-1;
    return pt_weight[ndx];
  }

  else return 1.0;

}


double seleXiMC::xiy_wgt( double y, int weight_choice ) {

  if( !weight_choice ) return 1.;

  else if( weight_choice == 1 ) {//D6T - updated for 362
    const int nbins = 10;
    const Double_t xAxis1[nbins+1] = {0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2};
    const Double_t y_weight[nbins] = {1.069187, 1.074479, 1.054569, 1.040567, 1.094412, 1.104698, 1.026344, 0.7999424, 0.8214608, 0.7656174};
    int ndx;
    for( ndx = 1; ndx < nbins+1; ndx++ ) {
      if( y < xAxis1[ndx] ) break;
    }
    ndx += -1;
    if ( ndx >= nbins ) ndx = nbins-1;
    return y_weight[ndx];
  }

  else if( weight_choice == 11 ) {//900 GeV D6T
    const int nbins = 11;
    const Double_t xAxis1[nbins+1] = {0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.3};
    const Double_t y_weight[nbins] = {0.9951, 1.009328, 1.109952, 1.114444, 0.9711025, 0.9514874, 0.8525454, 0.9190069, 0.9860223, 0.9657202};
    int ndx;
    for( ndx = 1; ndx < nbins+1; ndx++ ) {
      if( y < xAxis1[ndx] ) break;
    }
    ndx += -1;
    if ( ndx >= nbins ) ndx = nbins-1;
    return y_weight[ndx];
  }
  
 else return 1.0;  

}



void seleXiMC::SetOutputDir(std::string theOutputDir){
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

float seleXiMC::GetDeltaR(double eta1, double phi1, double eta2, double phi2){
  float dEta = eta1 - eta2;
  float dPhi = phi1 - phi2;
  float dR = sqrt(dEta*dEta + dPhi*dPhi);
  return dR;
}

