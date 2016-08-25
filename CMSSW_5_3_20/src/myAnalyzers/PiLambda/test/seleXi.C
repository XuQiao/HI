#define seleXi_cxx
#include "seleXi.h"
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

void seleXi::Loop(){
   
  cout << " Cuts Used:\t " << "Xi3dIp " << Xi3dIpCut << "\t batPi3dIp " << batPi3dIpCut << "\t v0Pi3dIp" << v0Pi3dIpCut 
                        << "\t v0Proton3dIp "  << v0Proton3dIpCut << "\t XiSep " << XiSepCut << "\t VtxCL " << vtxCLCut << "\n";
 
  
  SetBins();
  BookHistos();
  
/*
  gROOT->Reset();
  setTDRStyle();
*/
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
      
      /*
      //match genXi to a recoXi
      for(unsigned int nRecoXi = 0; nRecoXi < nXI; ++nRecoXi){
        //check DR between reco and gen tracks
      }
      */
      if ( fabs(Y) > 2.0) continue;
      
      double weightMC_y = xipt_wgt( pT, weight_choice );
      double weightMC_pt = xiy_wgt( Y, weight_choice );
      // for rapidity and pT we can consider all simulated Xi
      if ( !(processType==92||processType==93)) {
        denEffY->Fill( fabs(Y), weightMC_ntrk*weightMC_y);
        //denEffY->Fill( fabs(Y), weightMC_ntrk);
        denEffYUncor->Fill( fabs(Y) );
        allGenY->Fill( fabs(Y) );
        denEffpT->Fill(pT, weightMC_ntrk*weightMC_pt);
        //denEffpT->Fill(pT, weightMC_ntrk);
        denEffpTUncor->Fill(pT);
        allGenpT->Fill(pT);
	
	//std::cout << " pion charge " << pionCharge << std::endl;      
        // for ctau it's only possible to consider the MC events with at least one simulated track
        if ( (*genXiL)[i] < -1 ) continue;
        float ctau = ((*genXiL)[i]*mXi)/((*genXiP)[i]);
        denEff->Fill(ctau);
        denEffCtau->Fill(ctau, weightMC_ntrk*weightMC_y);
        //denEffCtau->Fill(ctau, weightMC_ntrk);
        denEffCtauUncor->Fill(ctau);
        
	//neg Xi
	if((*genXiQ)[i] < 0){
          denEffNegXipT->Fill(pT, weightMC_ntrk*weightMC_pt);
          //denEffNegXipT->Fill(pT, weightMC_ntrk);
          denEffNegXipTUncor->Fill(pT);
          allGenNegXipT->Fill(pT);
          denEffNegXiY->Fill( fabs(Y), weightMC_ntrk*weightMC_y);
          //denEffNegXiY->Fill( fabs(Y), weightMC_ntrk);
          denEffNegXiYUncor->Fill( fabs(Y) );
          allGenNegXiY->Fill( fabs(Y) );
	  
	}
	
	// pos Xi
	if((*genXiQ)[i] > 0){
          denEffPosXipT->Fill(pT, weightMC_ntrk*weightMC_pt);
          //denEffPosXipT->Fill(pT, weightMC_ntrk);
          denEffPosXipTUncor->Fill(pT);
          allGenPosXipT->Fill(pT);	   
          denEffPosXiY->Fill( fabs(Y), weightMC_ntrk*weightMC_y);
          //denEffPosXiY->Fill( fabs(Y), weightMC_ntrk);
          denEffPosXiYUncor->Fill( fabs(Y) );
          allGenPosXiY->Fill( fabs(Y) );
	}
      }
    }
    



        
    int nCands = 0;
    int bestCand = -999;
    int worstCand = -999;
    float bestSep = -999.0;
    float worstSep = 999.0;
    float bestV = 999.0;
    float worstV = -999.0;
        
    for ( int i = 0; i<nXI; i++ ) {
      //apply cuts
      //check NaNs first ...
      
      if (XiLongitudinalPCAPrimaryError->at(i) != XiLongitudinalPCAPrimaryError->at(i)) continue;
      if (XiTransversePCAPrimaryError->at(i) != XiTransversePCAPrimaryError->at(i)) continue;
      if (xiDecayVtxXE->at(i) != xiDecayVtxXE->at(i) )continue;
      if (xiDecayVtxYE->at(i) != xiDecayVtxYE->at(i) )continue;
      if (xiDecayVtxZE->at(i) != xiDecayVtxZE->at(i) )continue;
      
      if ( !trigHF || !trigTech34 || !trigTech36_39 || !trigScraping ) continue;
      if ( (*xiMass)[i] < 1.255 || (*xiMass)[i] > 1.455 ) continue;
      if ( fabs((*VMass)[i]-1.116)>0.008 ) continue;
      if ( (*xiVtxCL)[i] <  vtxCLCut) continue;
      if ( (*VTrkPiQ)[i]!=(*batPiQ)[i]) continue;
      
      if(batPi3DIpSig->at(i) < batPi3dIpCut)continue;
      if(VTrkPi3DIpSig->at(i) < v0Pi3dIpCut)continue;
      if(VTrkP3DIpSig->at(i) < v0Proton3dIpCut)continue;
      
      float sepX = (*xiDecayVtxX)[i]-(*priRefitVtxX)[i];
      float sepY = (*xiDecayVtxY)[i]-(*priRefitVtxY)[i];
      float sepZ = (*xiDecayVtxZ)[i]-(*priRefitVtxZ)[i];
      float sepXE2 = (*xiDecayVtxXE)[i]*(*xiDecayVtxXE)[i]+(*priRefitVtxXE)[i]*(*priRefitVtxXE)[i];
      float sepYE2 = (*xiDecayVtxYE)[i]*(*xiDecayVtxYE)[i]+(*priRefitVtxYE)[i]*(*priRefitVtxYE)[i];
      float sepZE2 = (*xiDecayVtxZE)[i]*(*xiDecayVtxZE)[i]+(*priRefitVtxZE)[i]*(*priRefitVtxZE)[i];
      float sepXi = sqrt(sepX*sepX+sepY*sepY+sepZ*sepZ);      
      float sepXiE = sqrt( sepXE2 + sepYE2 + sepZE2 );
      if ((sepXi/sepXiE) < XiSepCut) continue;
      
      float xi3DIp = sqrt(XiLongitudinalPCAPrimary->at(i) * XiLongitudinalPCAPrimary->at(i) + XiTransversePCAPrimary->at(i) * XiTransversePCAPrimary->at(i)); 
      float xi3DIpError = sqrt(((XiLongitudinalPCAPrimary->at(i) * XiLongitudinalPCAPrimaryError->at(i) * XiLongitudinalPCAPrimaryError->at(i)) + 	      
        		       (XiTransversePCAPrimary->at(i) * XiTransversePCAPrimaryError->at(i) * XiTransversePCAPrimaryError->at(i)))/xi3DIp);	      
      float xi3DIpSig = xi3DIp/xi3DIpError;														      
      																			      
      if ( xi3DIpSig > Xi3dIpCut ) continue;
      if ( priNTrk ==0 ) continue;

      //std::cout << "passed cand" << std::endl;
      if ( i>-1 ) {
        //pick best cand based on V decay radius - Xi decay radius
        float vR = sqrt( ((*VDecayVtxX)[i]-priVtxX)*((*VDecayVtxX)[i]-priVtxX) + ((*VDecayVtxY)[i]-priVtxY)*((*VDecayVtxY)[i]-priVtxY) );
        //std::cout << (*VDecayVtxX)[i] << "," << priVtxX << "," << (*VDecayVtxX)[i] << "," << priVtxX << std::endl;
        float xiR = sqrt( ((*xiDecayVtxX)[i]-priVtxX)*((*xiDecayVtxX)[i]-priVtxX) + ((*xiDecayVtxY)[i]-priVtxY)*((*xiDecayVtxY)[i]-priVtxY) );
        //std::cout << xiDecayVtxX << "," << priVtxX << "," << xiDecayVtxX << "," << priVtxX << std::endl;	
        float sep = vR-xiR;

        //if ( sep > bestSep ) {
        //  bestCand = i;
        //  bestSep = sep;
        //}
        //if ( sep < worstSep ) {
        //  worstCand = i;
        //  worstSep = sep;
        //}
        if ( fabs( (*VMass)[i]-1.11568) < bestV ) {
          bestCand = i;
          bestV = fabs( (*VMass)[i]-1.11568);
        }
        if ( fabs( (*VMass)[i]-1.11568) > worstV ) {
          worstCand = i;
          worstV = fabs( (*VMass)[i]-1.11568);
        }
      }
      nCands++;
      nTotalCands++;

    }
    if (nCands > 0) nEvtWCands++;
    if ( nCands==1 ) {
      nSingleCands++;


//bestCand = i;

      data->Fill( (*xiMass)[bestCand]*1000 );
      dataOm->Fill( (*omMass)[bestCand]*1000 );
      //std::cout << "Filling with number of cands = " << nCands << "  xiMass = " << (*xiMass)[bestCand] << "  cand number " << bestCand << std::endl;
      dataBad->Fill( (*xiMass)[worstCand]*1000 );
      //std::cout << "Filling with number of cands = " << nCands << "  xiMass = " << (*xiMass)[worstCand] << "  cand number " << worstCand << std::endl;
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

      float delRPi = 0.05;
      float delRP = 0.05;
      float delRXiCut = 5.0;
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
		   
        if ( (delRbatPi < delRPi) &&
             (delRlamPi < delRPi) && 
             (delRlamP < delRP) &&
	     (delRXi < delRXiCut) ) {
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
	
	//std::cout << "true pT = " << genpT << " and reco pT = " << xiPt << std::endl;
        //std::cout << "true Y = " << genY << " and reco Y = " << xiY << std::endl;
	
        //this can be negative and my histo cuts off at 0.
        numEffCtau->Fill( fabs( (*xictauMPV)[bestCand] ), weightMC_ntrk*weightMC_y );
        //numEffCtau->Fill( fabs( (*xictauMPV)[bestCand] ), weightMC_ntrk );
        numEffCtauUncor->Fill( fabs( (*xictauMPV)[bestCand] ) );
	
	numEffY->Fill( xiY, weightMC_ntrk*weightMC_y );
	//numEffY->Fill( xiY, weightMC_ntrk );
	numEffYUncor->Fill( xiY );	
	
	numEffpT->Fill( xiPt, weightMC_ntrk*weightMC_pt );
	//numEffpT->Fill( xiPt, weightMC_ntrk );
        numEffpTUncor->Fill( xiPt );
	
	//pos Xi
	if((*batPiQ)[bestCand] > 0){
	  numEffPosXipT->Fill( xiPt, weightMC_ntrk*weightMC_pt ); 
	  //numEffPosXipT->Fill( xiPt, weightMC_ntrk ); 
          numEffPosXipTUncor->Fill( xiPt );			     
	  numEffPosXiY->Fill( xiY, weightMC_ntrk*weightMC_y );
	  //numEffPosXiY->Fill( xiY, weightMC_ntrk );
	  numEffPosXiYUncor->Fill( xiY );	  
	
	  
	}
	//neg Xi
	if((*batPiQ)[bestCand] < 0){
	  numEffNegXipT->Fill( xiPt, weightMC_ntrk*weightMC_pt );
	  //numEffNegXipT->Fill( xiPt, weightMC_ntrk );
          numEffNegXipTUncor->Fill( xiPt );
	  numEffNegXiY->Fill( xiY, weightMC_ntrk*weightMC_y ); 
	  //numEffNegXiY->Fill( xiY, weightMC_ntrk ); 
	  numEffNegXiYUncor->Fill( xiY );	  		  
	
	}
		
        numEff->Fill( fabs( (*xictauMPV)[bestCand] ) );
        dataTruth->Fill( (*xiMass)[bestCand]*1000 );
	if(xiPt < 0.6 && xiPt > 0.4)dataTruthLowPt->Fill( (*xiMass)[bestCand]*1000 );

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
	if(xiPt < 0.6 && xiPt > 0.4)dataNotTruthLowPt->Fill( (*xiMass)[bestCand]*1000 );
	
      }
    }  
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
  f1->SetParLimits(2,0,30000);
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
  TText *text1 = stats->AddText(Form("%s%.1f%s%.1f", "Yield: ", nSig, " #pm ",nSigE));
  TText *text2 = stats->AddText(Form("%s%.2f%s%.2f%s", "Mean: ", mean, " #pm ", meanE, " MeV/c^{2}"));
  TText *text3 = stats->AddText(Form("%s%.2f%s%.2f%s", "Sigma: ", sigma, " #pm ", sigmaE, " MeV/c^{2}"));
  TText *text4 = stats->AddText("Statistical uncertainties only");
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
  dataTruthLowPt->GetXaxis()->SetTitle("#Xi mass from truth-matched events 0.4 < Pt < 0.6");
  dataTruthLowPt->Draw();
  c4->cd(2);
  dataNotTruthLowPt->GetXaxis()->SetTitle("#Xi mass from non-truth-matched events 0.4 < Pt < 0.6");
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

  
  TFile *fout;
  std::string theRootOutFile = theOutDir + "/efficienciesXi_" + theSqrtS + ".root";
  fout = TFile::Open(theRootOutFile.c_str(), "RECREATE");

  fout->cd();
  TH1F* data_copy = (TH1F*)data->Clone("data_copy");
  TH1F* dataOm_copy = (TH1F*)dataOm->Clone("dataOm_copy");
  TH1F* denEff_copy = (TH1F*)denEff->Clone("denEff_copy");
  TH1F* numEff_copy = (TH1F*)numEff->Clone("numEff_copy");
  TH1F* effH_copy = (TH1F*)effH->Clone("effH_copy");
  
  TH1F* dataTruth_copy = (TH1F*)dataTruth->Clone("dataTruth_copy");
  TH1F* dataNotTruth_copy = (TH1F*)dataNotTruth->Clone("dataNotTruth_copy");
  
  TH1F* effCtau_copy = (TH1F*)effCtau->Clone("effCtau_copy");
  TH1F* allGenCtau_copy = (TH1F*)allGenCtau->Clone("allGenCtau_copy");
  TH1F* numEffCtau_copy = (TH1F*)numEffCtau->Clone("numEffCtau_copy");
  TH1F* denEffCtau_copy = (TH1F*)denEffCtau->Clone("denEffCtau_copy");
  TH1F* denEffCtauUncor_copy = (TH1F*)denEffCtauUncor->Clone("denEffCtauUncor_copy");
  TH1F* numEffCtauScaled_copy = (TH1F*)numEffCtauScaled->Clone("numEffCtauScaled_copy");  
  TH1F* effCtauUncor_copy = (TH1F*)effCtauUncor->Clone("effCtauUncor_copy");
  
  
  TH1F* allGenpT_copy = (TH1F*)allGenpT->Clone("allGenpT_copy");
  TH1F* effpT_copy = (TH1F*)effpT->Clone("effpT_copy");
  TH1F* numEffpT_copy = (TH1F*)numEffpT->Clone("numEffpT_copy");
  TH1F* denEffpT_copy = (TH1F*)denEffpT->Clone("denEffpT_copy");
  TH1F* denEffpTUncor_copy = (TH1F*)denEffpTUncor->Clone("denEffpTUncor_copy"); 
  TH1F* numEffpTScaled_copy = (TH1F*)numEffpTScaled->Clone("numEffpTScaled_copy");
  TH1F* effpTUncor_copy = (TH1F*)effpTUncor->Clone("effpTUncor_copy");

  TH1F* allGenPosXipT_copy = (TH1F*)allGenPosXipT->Clone("allGenPosXipT_copy");
  TH1F* effPosXipT_copy = (TH1F*)effPosXipT->Clone("effPosXipT_copy");
  TH1F* numEffPosXipT_copy = (TH1F*)numEffPosXipT->Clone("numEffPosXipT_copy");
  TH1F* denEffPosXipT_copy = (TH1F*)denEffPosXipT->Clone("denEffPosXipT_copy");
  TH1F* numEffPosXipTScaled_copy = (TH1F*)numEffPosXipTScaled->Clone("numEffPosXipTScaled_copy");
  TH1F* effPosXipTUncor_copy = (TH1F*)effPosXipTUncor->Clone("effPosXipTUncor_copy");

  TH1F* allGenNegXipT_copy = (TH1F*)allGenNegXipT->Clone("allGenNegXipT_copy");
  TH1F* effNegXipT_copy = (TH1F*)effNegXipT->Clone("effNegXipT_copy");
  TH1F* numEffNegXipT_copy = (TH1F*)numEffNegXipT->Clone("numEffNegXipT_copy");
  TH1F* denEffNegXipT_copy = (TH1F*)denEffNegXipT->Clone("denEffNegXipT_copy");
  TH1F* numEffNegXipTScaled_copy = (TH1F*)numEffNegXipTScaled->Clone("numEffNegXipTScaled_copy");
  TH1F* effNegXipTUncor_copy = (TH1F*)effNegXipTUncor->Clone("effNegXipTUncor_copy");
  
  TH1F* allGenY_copy = (TH1F*)allGenY->Clone("allGenY_copy");
  TH1F* effY_copy = (TH1F*)effY->Clone("effY_copy");
  TH1F* numEffY_copy = (TH1F*)numEffY->Clone("numEffY_copy");
  TH1F* denEffY_copy = (TH1F*)denEffY->Clone("denEffY_copy");
  TH1F* denEffYUncor_copy = (TH1F*)denEffYUncor->Clone("denEffYUncor_copy");
  TH1F* numEffYScaled_copy = (TH1F*)numEffYScaled->Clone("numEffYScaled_copy");
  TH1F* effYUncor_copy = (TH1F*)effYUncor->Clone("effYUncor_copy");  

  TH1F* allGenPosXiY_copy = (TH1F*)allGenPosXiY->Clone("allGenPosXiY_copy");
  TH1F* effPosXiY_copy = (TH1F*)effPosXiY->Clone("effPosXiY_copy");
  TH1F* numEffPosXiY_copy = (TH1F*)numEffPosXiY->Clone("numEffPosXiY_copy");
  TH1F* denEffPosXiY_copy = (TH1F*)denEffPosXiY->Clone("denEffPosXiY_copy");
  TH1F* numEffPosXiYScaled_copy = (TH1F*)numEffPosXiYScaled->Clone("numEffPosXiYScaled_copy");
  TH1F* effPosXiYUncor_copy = (TH1F*)effPosXiYUncor->Clone("effPosXiYUncor_copy");  

  TH1F* allGenNegXiY_copy = (TH1F*)allGenNegXiY->Clone("allGenNegXiY_copy");
  TH1F* effNegXiY_copy = (TH1F*)effNegXiY->Clone("effNegXiY_copy");
  TH1F* numEffNegXiY_copy = (TH1F*)numEffNegXiY->Clone("numEffNegXiY_copy");
  TH1F* denEffNegXiY_copy = (TH1F*)denEffNegXiY->Clone("denEffNegXiY_copy");
  TH1F* numEffNegXiYScaled_copy = (TH1F*)numEffNegXiYScaled->Clone("numEffNegXiYScaled_copy");
  TH1F* effNegXiYUncor_copy = (TH1F*)effNegXiYUncor->Clone("effNegXiYUncor_copy");  

  fout->Write();
}

void seleXi::SetBins(){
   
   // Pt bins
   //thePtBins.push_back(0.4);
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

void seleXi::BookHistos(){
  std::string theOutRootFileName = theOutDir + "/xiOut.root";	      
  theOutRootFile = new TFile(theOutRootFileName.c_str(), "Recreate"); 
  theOutRootFile->cd(); 	
  data = new TH1F("data", "data", 100, 1255, 1455);
  dataOm = new TH1F("dataOm", "dataOm", 100, 1455, 1855);
  dataBad = new TH1F("dataBad", "dataBad", 100, 1255, 1455);
  dataTruth = new TH1F("dataTruth", "dataTruth", 100, 1255, 1455);
  dataNotTruth = new TH1F("dataNotTruth", "dataNotTruth", 100, 1255, 1455);
  
  dataTruthLowPt = new TH1F("dataTruth", "dataTruth", 100, 1255, 1455);
  dataNotTruthLowPt = new TH1F("dataNotTruth", "dataNotTruth", 100, 1255, 1455);
  
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

void seleXi::FitXi(TH1F* theXiHisto){
    //initialize results to 0
  theFitYield = 0.;
  theFitYieldError = 0.;

  
  TAxis* theAxis = theXiHisto->GetXaxis();
  int x_min = (int)theAxis->GetXmin();
  int x_max = (int)theAxis->GetXmax();
  int bins = theAxis->GetNbins();
  int x_range = x_max - x_min;
  TF1* xiFit = new TF1("xiFit", Form("%s%i%s%i%s%i%s%i%s","[0]*((x-", x_min, ")**(0.5)) +[1]*((x-", x_min, ")**(3/2)) + ((", x_range, ")/", bins, ")*gausn(2)"), x_min, x_max);
  xiFit->SetParameters(0.15, 0., theXiHisto->GetEntries()*.5, 1321, 3);
  xiFit->SetParName(0,"bkg param 1");
  xiFit->SetParName(1,"bkg param 2");
  xiFit->SetParName(2,"yield");
  xiFit->SetParLimits(2,0,20000);
  xiFit->SetParName(3,"mean");
  xiFit->SetParName(4,"sigma");
  xiFit->SetParLimits(4,1,6);
  theXiHisto->Fit("xiFit", "LEM");
  
  theFitYield = xiFit->GetParameter(2);
  theFitYieldError = xiFit->GetParError(2);
  
}

double seleXi::ntrk_wgt(int ntrk, int weight_choice) {


  if ( weight_choice == 1) { // tuneD6T, 7 TeV, no primary vertex
    const int nbins = 75;
    const Double_t xAxis1[nbins+1] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 84, 88, 92, 96, 100, 110, 120, 140, 160, 200}; 
    const Double_t ntrk_weight[nbins] = {0.6511886,0.5558178,0.5587921,0.5788214,0.6277799,0.7018213,0.7895163,0.868456,0.9275547,0.9485932,0.947093,0.9298214,0.9037366,0.8802626,0.8631144,0.8523428,0.8489775,0.8536014,0.8612283,0.8771418,0.8990253,0.914389,0.9375328,0.9600454,0.9862425,1.010019,1.032768,1.05441,1.086065,1.114157,1.139019,1.165673,1.196096,1.240339,1.264464,1.297802,1.339239,1.369484,1.40494,1.442174,1.481435,1.51453,1.545663,1.574016,1.623789,1.657342,1.686259,1.716283,1.73799,1.774361,1.80698,1.861587,1.905641,1.933265,1.963547,1.986199,1.999039,2.009495,2.012197,2.014793,2.01983,2.010925,2.016263,1.99989,2.016479,2.001244,2.041022,2.083046,2.157776,2.285894,2.603289,3.525266,6.353901,21.28761,132.5661};
    int ndx;
    for ( ndx = 1; ndx <  nbins+1; ndx++ ) {
      if ( ntrk < xAxis1[ndx] ) break;
    }
    ndx += -1;
    if ( ndx >= nbins ) ndx = nbins-1;
    return ntrk_weight[ndx];


  } else if ( weight_choice == 2) { // pythia 8, 7 TeV, no primary vertex
    const int nbins = 75;
    Double_t xAxis1[nbins+1] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 84, 88, 92, 96, 100, 110, 120, 140, 160, 200}; 
    const Double_t ntrk_weight[nbins] = {0.7583848,1.364867,1.932691,2.132744,1.986891,1.710771,1.436661,1.215889,1.076008,0.964021,0.8933408,0.8513227,0.8186476,0.8019238,0.7942611,0.7953302,0.8001245,0.8126709,0.8273689,0.8487046,0.8663292,0.8793734,0.9006234,0.9189676,0.9422619,0.9514298,0.9713438,0.9724804,0.9914346,1.000026,1.007341,1.003133,1.020771,1.031123,1.037387,1.040028,1.043577,1.044883,1.052531,1.056744,1.063543,1.063997,1.069665,1.063039,1.067783,1.081269,1.073623,1.076372,1.082199,1.080386,1.088779,1.090386,1.097871,1.100716,1.113368,1.115761,1.118757,1.127538,1.155251,1.159525,1.18015,1.190233,1.21353,1.226885,1.246372,1.291205,1.352633,1.420856,1.481615,1.584355,1.766239,2.111093,2.865317,4.697073,8.972905};
    int ndx;
    for( ndx = 1; ndx < nbins+1; ndx++ ) {
      if( ntrk < xAxis1[ndx] ) break;
    }
    ndx += -1;
    if ( ndx >= nbins ) ndx = nbins-1;
    return ntrk_weight[ndx];


  } else if ( weight_choice == 3) { // tuneP0, 7 TeV, no primary vertex
    const int nbins = 75;
    Double_t xAxis1[nbins+1] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 84, 88, 92, 96, 100, 110, 120, 140, 160, 200}; 
    // 358 results:   const Double_t ntrk_weight[nbins] = {0.6113107, 0.6113107, 0.6113107, 0.6053922, 0.6343482, 0.6575175, 0.6798984, 0.6848243, 0.6929444, 0.7154517, 0.7475836, 0.7882146, 0.8322517, 0.8857094, 0.9359018, 0.9730012, 1.014162, 1.043932, 1.068149, 1.089636, 1.104945, 1.121409, 1.140087, 1.142924, 1.163503, 1.18115, 1.189013, 1.20008, 1.211411, 1.21069, 1.223679, 1.233176, 1.240755, 1.254538, 1.267727, 1.258991, 1.281829, 1.2781, 1.29623, 1.289593, 1.295538, 1.309365, 1.318696, 1.312257, 1.319332, 1.322189, 1.343504, 1.330006, 1.358613, 1.359512, 1.376825, 1.393088, 1.41033, 1.438139, 1.450585, 1.477302, 1.507699, 1.542305, 1.573349, 1.608912, 1.665798, 1.701407, 1.764046, 1.81462, 1.85989, 1.992459, 2.148923, 2.303676, 2.520458, 2.791737, 3.343563, 4.507767, 6.967413, 14.67809, 40.87044};
    const Double_t ntrk_weight[nbins] = {0.6600672,0.5262316,0.5291787,0.5576365,0.6042867,0.648715,0.6794875,0.6958973,0.7119868,0.7347366,0.770914,0.8134168,0.8596728,0.9161927,0.9645801,1.007189,1.048433,1.087396,1.113726,1.138973,1.155527,1.172295,1.186989,1.21233,1.223847,1.235815,1.254067,1.259997,1.267299,1.280645,1.294102,1.301874,1.308915,1.326859,1.334744,1.322762,1.346465,1.355852,1.364503,1.373898,1.376463,1.36924,1.388751,1.383131,1.409655,1.419957,1.409392,1.416712,1.427729,1.436744,1.446881,1.470225,1.499659,1.526046,1.535751,1.558639,1.584844,1.628592,1.664813,1.703127,1.756106,1.785631,1.857253,1.903782,1.952839,2.081257,2.245752,2.433308,2.629362,2.908451,3.459108,4.647993,7.277652,15.88272,42.22966};
    int ndx;
    for( ndx = 1; ndx < nbins+1; ndx++ ) {
      if( ntrk < xAxis1[ndx] ) break;
    }
    ndx += -1;
    if ( ndx >= nbins ) ndx = nbins-1;
    return ntrk_weight[ndx];

  } else if ( weight_choice == 11) { // D6T 900 GeV
    const int nbins = 75;
    const Double_t xAxis1[nbins+1] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 84, 88, 92, 96, 100, 110, 120, 140, 160, 200}; 
    // 358 results:    const Double_t ntrk_weight[nbins] = {0.7120804, 0.7120804, 0.7120804, 0.7560743, 0.8187093, 0.9024728, 0.9692311, 1.021291, 1.037656, 1.022344, 0.9795186, 0.9338584, 0.8938458, 0.858408, 0.8363768, 0.8210503, 0.8186527, 0.8290459, 0.8487917, 0.8681399, 0.8971176, 0.9226106, 0.9572074, 0.981635, 1.016574, 1.040278, 1.075404, 1.09967, 1.138376, 1.162668, 1.189298, 1.210598, 1.250772, 1.290503, 1.32034, 1.356436, 1.389858, 1.443522, 1.466343, 1.517786, 1.556459, 1.614047, 1.675907, 1.735429, 1.816951, 1.89242, 2.009592, 2.075383, 2.109613, 2.256741, 2.475007, 2.815125, 3.172367, 3.698223, 4.235982, 4.996351, 6.223412, 7.4198, 9.075213, 10.99792, 14.65506, 16.56022, 22.13489, 28.80532, 33.0667, 56.87799, 113.7677, 178.7183, 200.0, 200.0, 200.0, 200.0, 200.0, 200.0, 200.0};
    // Frankenstein
    //const Double_t ntrk_weight[nbins] = {0.8121125,0.726478,0.7873486,0.830931,0.8797048,0.9200454,0.9332961,0.924447,0.9142907,0.8979387,0.8837417,0.8802965,0.882481,0.8913741,0.9080257,0.9301792,0.952966,0.9873891,1.018473,1.046005,1.077875,1.096652,1.126845,1.149197,1.166814,1.177694,1.203604,1.21139,1.225665,1.229918,1.24111,1.256053,1.267938,1.27574,1.270057,1.286223,1.299729,1.295669,1.318001,1.316429,1.327181,1.329701,1.327778,1.372345,1.379321,1.379946,1.387677,1.414408,1.438566,1.443804,1.458181,1.508292,1.551891,1.586478,1.666062,1.729609,1.780928,1.854244,1.96099,2.126867,2.182382,2.305395,2.394886,2.585145,2.878618,3.129298,3.77357,4.551805,5.211093,7.68188,9.144209,16.42705,37.79441, 40.0, 50.0};
    const Double_t ntrk_weight[nbins] = {0.8008507, 0.6800452, 0.7053241, 0.729585, 0.7904063, 0.868402, 0.9467011, 1.007388, 1.038141, 1.027983, 0.9986413, 0.9547896, 0.9153704, 0.8825721, 0.8606647, 0.8488371, 0.8442871, 0.8612404, 0.8795851, 0.8989031, 0.9376439, 0.9658297, 0.999901, 1.0371, 1.070309, 1.10116, 1.159267, 1.182447, 1.225471, 1.245174, 1.285313, 1.316964, 1.355525, 1.39009, 1.434817, 1.472543, 1.504442, 1.538746, 1.594955, 1.622222, 1.665173, 1.703172, 1.769843, 1.806658, 1.920535, 1.953628, 2.023154, 2.081792, 2.169408, 2.304216, 2.430937, 2.7744, 3.072114, 3.438184, 3.916324, 4.685513, 5.367588, 6.116444, 8.080085, 10.04357, 11.66427, 13.59402, 19.28118, 23.5173, 28.68605, 48.488, 60.63322, 156.9263, 138.4146, 527.6962};
    int ndx;
    for ( ndx = 1; ndx <  nbins+1; ndx++ ) {
      if ( ntrk < xAxis1[ndx] ) break;
    }
    ndx += -1;
    if ( ndx >= nbins ) ndx = nbins-1;
    return ntrk_weight[ndx];

  } else if ( weight_choice == 12) { // pythia8 900 GeV
    const int nbins = 75;
    const Double_t xAxis1[nbins+1] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 84, 88, 92, 96, 100, 110, 120, 140, 160, 200}; 
    // 358 results:    const Double_t ntrk_weight[nbins] = {0.7120804, 0.7120804, 0.7120804, 0.7560743, 0.8187093, 0.9024728, 0.9692311, 1.021291, 1.037656, 1.022344, 0.9795186, 0.9338584, 0.8938458, 0.858408, 0.8363768, 0.8210503, 0.8186527, 0.8290459, 0.8487917, 0.8681399, 0.8971176, 0.9226106, 0.9572074, 0.981635, 1.016574, 1.040278, 1.075404, 1.09967, 1.138376, 1.162668, 1.189298, 1.210598, 1.250772, 1.290503, 1.32034, 1.356436, 1.389858, 1.443522, 1.466343, 1.517786, 1.556459, 1.614047, 1.675907, 1.735429, 1.816951, 1.89242, 2.009592, 2.075383, 2.109613, 2.256741, 2.475007, 2.815125, 3.172367, 3.698223, 4.235982, 4.996351, 6.223412, 7.4198, 9.075213, 10.99792, 14.65506, 16.56022, 22.13489, 28.80532, 33.0667, 56.87799, 113.7677, 178.7183, 200.0, 200.0, 200.0, 200.0, 200.0, 200.0, 200.0};
    // Frankenstein
    //const Double_t ntrk_weight[nbins] = {0.8121125,0.726478,0.7873486,0.830931,0.8797048,0.9200454,0.9332961,0.924447,0.9142907,0.8979387,0.8837417,0.8802965,0.882481,0.8913741,0.9080257,0.9301792,0.952966,0.9873891,1.018473,1.046005,1.077875,1.096652,1.126845,1.149197,1.166814,1.177694,1.203604,1.21139,1.225665,1.229918,1.24111,1.256053,1.267938,1.27574,1.270057,1.286223,1.299729,1.295669,1.318001,1.316429,1.327181,1.329701,1.327778,1.372345,1.379321,1.379946,1.387677,1.414408,1.438566,1.443804,1.458181,1.508292,1.551891,1.586478,1.666062,1.729609,1.780928,1.854244,1.96099,2.126867,2.182382,2.305395,2.394886,2.585145,2.878618,3.129298,3.77357,4.551805,5.211093,7.68188,9.144209,16.42705,37.79441, 40.0, 50.0};
    // Pythia 8
    const Double_t ntrk_weight[nbins] = {0.8698299, 1.133929, 1.505682, 1.568508, 1.420571, 1.217741, 1.03395, 0.9052682, 0.8148503, 0.7643023, 0.740001, 0.7338403, 0.7476363, 0.7723145, 0.8071364, 0.8464819, 0.8839312, 0.9316061, 0.9664785, 1.004444, 1.03491, 1.06311, 1.085548, 1.110256, 1.128619, 1.153715, 1.171373, 1.176932, 1.196115, 1.20767, 1.230101, 1.250178, 1.27388, 1.284152, 1.293236, 1.310968, 1.336388, 1.357607, 1.398819, 1.412168, 1.433548, 1.461427, 1.480794, 1.534181, 1.546544, 1.574622, 1.586213, 1.663781, 1.687251, 1.722918, 1.769625, 1.856937, 1.924528, 2.051649, 2.182902, 2.305244, 2.396302, 2.599706, 2.813186, 3.080441, 3.193056, 3.361284, 3.633754, 3.643556, 4.37814, 4.373274, 5.211234, 6.923011, 7.01898, 8.638114, 12.14617, 14.47733, 20.18703};
    int ndx;
    for ( ndx = 1; ndx <  nbins+1; ndx++ ) {
      if ( ntrk < xAxis1[ndx] ) break;
    }
    ndx += -1;
    if ( ndx >= nbins ) ndx = nbins-1;
    return ntrk_weight[ndx];

  } else if ( weight_choice == 13) { // P0 900 GeV
    const int nbins = 75;
    const Double_t xAxis1[nbins+1] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 84, 88, 92, 96, 100, 110, 120, 140, 160, 200}; 
    // 358 results:    const Double_t ntrk_weight[nbins] = {0.7120804, 0.7120804, 0.7120804, 0.7560743, 0.8187093, 0.9024728, 0.9692311, 1.021291, 1.037656, 1.022344, 0.9795186, 0.9338584, 0.8938458, 0.858408, 0.8363768, 0.8210503, 0.8186527, 0.8290459, 0.8487917, 0.8681399, 0.8971176, 0.9226106, 0.9572074, 0.981635, 1.016574, 1.040278, 1.075404, 1.09967, 1.138376, 1.162668, 1.189298, 1.210598, 1.250772, 1.290503, 1.32034, 1.356436, 1.389858, 1.443522, 1.466343, 1.517786, 1.556459, 1.614047, 1.675907, 1.735429, 1.816951, 1.89242, 2.009592, 2.075383, 2.109613, 2.256741, 2.475007, 2.815125, 3.172367, 3.698223, 4.235982, 4.996351, 6.223412, 7.4198, 9.075213, 10.99792, 14.65506, 16.56022, 22.13489, 28.80532, 33.0667, 56.87799, 113.7677, 178.7183, 200.0, 200.0, 200.0, 200.0, 200.0, 200.0, 200.0};
    // Frankenstein
    //const Double_t ntrk_weight[nbins] = {0.8121125,0.726478,0.7873486,0.830931,0.8797048,0.9200454,0.9332961,0.924447,0.9142907,0.8979387,0.8837417,0.8802965,0.882481,0.8913741,0.9080257,0.9301792,0.952966,0.9873891,1.018473,1.046005,1.077875,1.096652,1.126845,1.149197,1.166814,1.177694,1.203604,1.21139,1.225665,1.229918,1.24111,1.256053,1.267938,1.27574,1.270057,1.286223,1.299729,1.295669,1.318001,1.316429,1.327181,1.329701,1.327778,1.372345,1.379321,1.379946,1.387677,1.414408,1.438566,1.443804,1.458181,1.508292,1.551891,1.586478,1.666062,1.729609,1.780928,1.854244,1.96099,2.126867,2.182382,2.305395,2.394886,2.585145,2.878618,3.129298,3.77357,4.551805,5.211093,7.68188,9.144209,16.42705,37.79441, 40.0, 50.0};
    const Double_t ntrk_weight[nbins] = {0.8558056, 0.6407431, 0.6697215, 0.7030422, 0.7507221, 0.7818462, 0.7956233, 0.7966417, 0.800523, 0.8108643, 0.8315852, 0.8687197, 0.9109246, 0.9578448, 1.003101, 1.048422, 1.079929, 1.116, 1.143368, 1.169906, 1.198146, 1.204486, 1.233517, 1.246131, 1.257286, 1.271608, 1.294502, 1.301131, 1.319941, 1.311432, 1.334922, 1.343651, 1.354945, 1.378523, 1.380522, 1.393882, 1.407326, 1.413659, 1.449031, 1.457864, 1.47807, 1.490332, 1.503246, 1.524135, 1.549498, 1.583332, 1.586316, 1.646677, 1.681149, 1.726774, 1.730233, 1.819326, 1.922571, 2.02323, 2.159004, 2.253778, 2.333775, 2.60474, 2.75771, 3.069827, 3.395023, 3.56286, 3.920214, 4.371885, 4.660389, 5.794777, 6.776816, 8.194007, 12.1028, 15.38036, 19.6759, 78.68833, 43.88887};
    int ndx;
    for ( ndx = 1; ndx <  nbins+1; ndx++ ) {
      if ( ntrk < xAxis1[ndx] ) break;
    }
    ndx += -1;
    if ( ndx >= nbins ) ndx = nbins-1;
    return ntrk_weight[ndx];

  } else if ( weight_choice == 14) { // Frankenstein 2, 900 GeV
    const int nbins = 75;
    const Double_t xAxis1[nbins+1] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 84, 88, 92, 96, 100, 110, 120, 140, 160, 200}; 
    // 358 results:    const Double_t ntrk_weight[nbins] = {0.7120804, 0.7120804, 0.7120804, 0.7560743, 0.8187093, 0.9024728, 0.9692311, 1.021291, 1.037656, 1.022344, 0.9795186, 0.9338584, 0.8938458, 0.858408, 0.8363768, 0.8210503, 0.8186527, 0.8290459, 0.8487917, 0.8681399, 0.8971176, 0.9226106, 0.9572074, 0.981635, 1.016574, 1.040278, 1.075404, 1.09967, 1.138376, 1.162668, 1.189298, 1.210598, 1.250772, 1.290503, 1.32034, 1.356436, 1.389858, 1.443522, 1.466343, 1.517786, 1.556459, 1.614047, 1.675907, 1.735429, 1.816951, 1.89242, 2.009592, 2.075383, 2.109613, 2.256741, 2.475007, 2.815125, 3.172367, 3.698223, 4.235982, 4.996351, 6.223412, 7.4198, 9.075213, 10.99792, 14.65506, 16.56022, 22.13489, 28.80532, 33.0667, 56.87799, 113.7677, 178.7183, 200.0, 200.0, 200.0, 200.0, 200.0, 200.0, 200.0};
    // Frankenstein
    //const Double_t ntrk_weight[nbins] = {0.8121125,0.726478,0.7873486,0.830931,0.8797048,0.9200454,0.9332961,0.924447,0.9142907,0.8979387,0.8837417,0.8802965,0.882481,0.8913741,0.9080257,0.9301792,0.952966,0.9873891,1.018473,1.046005,1.077875,1.096652,1.126845,1.149197,1.166814,1.177694,1.203604,1.21139,1.225665,1.229918,1.24111,1.256053,1.267938,1.27574,1.270057,1.286223,1.299729,1.295669,1.318001,1.316429,1.327181,1.329701,1.327778,1.372345,1.379321,1.379946,1.387677,1.414408,1.438566,1.443804,1.458181,1.508292,1.551891,1.586478,1.666062,1.729609,1.780928,1.854244,1.96099,2.126867,2.182382,2.305395,2.394886,2.585145,2.878618,3.129298,3.77357,4.551805,5.211093,7.68188,9.144209,16.42705,37.79441, 40.0, 50.0};
    const Double_t ntrk_weight[nbins] = {0.8471166, 0.7747086, 0.8534236, 0.8906906, 0.920615, 0.9248136, 0.9061951, 0.878546, 0.8519395, 0.8337365, 0.8266824, 0.8304441, 0.8443437, 0.8644182, 0.8893438, 0.9174435, 0.9415596, 0.9770847, 1.005003, 1.034081, 1.066713, 1.08815, 1.116346, 1.140832, 1.16103, 1.184116, 1.214378, 1.22483, 1.249853, 1.256156, 1.283204, 1.301786, 1.324212, 1.344955, 1.359037, 1.379283, 1.401386, 1.419371, 1.461308, 1.475588, 1.500716, 1.524153, 1.54962, 1.586371, 1.62218, 1.653815, 1.671273, 1.738827, 1.777745, 1.833496, 1.872641, 1.995974, 2.109257, 2.253059, 2.425031, 2.591801, 2.722196, 3.009602, 3.292291, 3.675596, 3.965288, 4.200694, 4.66012, 4.959784, 5.639791, 6.356228, 7.522874, 9.708465, 11.52256, 14.54791, 19.84394, 32.93715, 36.74176};
    int ndx;
    for ( ndx = 1; ndx <  nbins+1; ndx++ ) {
      if ( ntrk < xAxis1[ndx] ) break;
    }
    ndx += -1;
    if ( ndx >= nbins ) ndx = nbins-1;
    return ntrk_weight[ndx];

  } else if ( weight_choice == 21) { // tuneD6T, 7 TeV, no primary vertex, no HF
    const int nbins = 75;
    Double_t xAxis1[nbins+1] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 84, 88, 92, 96, 100, 110, 120, 140, 160, 200}; 

    const Double_t ntrk_weight[nbins] = {1.311728,0.6880677,0.6397273,0.6173317,0.641782,0.6975554,0.7766888,0.8539821,0.9190512,0.9518477,0.9597228,0.9496967,0.9283509,0.9063073,0.8897493,0.8771079,0.8715684,0.8724656,0.8770419,0.8890017,0.9071344,0.9184986,0.9391029,0.9580962,0.9813579,1.001935,1.021457,1.04107,1.069826,1.095968,1.118448,1.142968,1.171203,1.212617,1.234977,1.26599,1.305322,1.333993,1.367479,1.402984,1.440134,1.471124,1.500586,1.527247,1.575015,1.607348,1.634649,1.662788,1.68336,1.718395,1.748973,1.801136,1.843314,1.868926,1.898041,1.919752,1.931917,1.941794,1.943892,1.946087,1.950732,1.941953,1.947075,1.931155,1.947152,1.932253,1.970545,2.010956,2.083081,2.206729,2.513102,3.402988,6.133433,20.54897,127.9663};
    int ndx;
    for( ndx = 1; ndx < nbins+1; ndx++ ) {
      if( ntrk < xAxis1[ndx] ) break;
    }
    ndx += -1;
    if ( ndx >= nbins ) ndx = nbins-1;
    return ntrk_weight[ndx];

  } else if ( weight_choice == 22) { // pythia 8, 7 TeV, no primary vertex, no HF
    const int nbins = 75;
    Double_t xAxis1[nbins+1] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 84, 88, 92, 96, 100, 110, 120, 140, 160, 200}; 

    const Double_t ntrk_weight[nbins] = {1.32619,1.21876,1.490618,1.591567,1.56972,1.437679,1.280853,1.122554,1.013539,0.9228153,0.8616974,0.8238635,0.7940758,0.7788984,0.7731082,0.7731055,0.7781804,0.7908646,0.8053773,0.8256868,0.8420982,0.8550634,0.8760663,0.8939818,0.9164134,0.9260742,0.9437479,0.945897,0.9638105,0.9733848,0.9799216,0.9762368,0.9936696,1.002548,1.00889,1.011781,1.014923,1.016529,1.023796,1.027992,1.034122,1.034637,1.039988,1.033122,1.037984,1.051276,1.043341,1.046684,1.051882,1.050065,1.05786,1.059385,1.066966,1.069138,1.081687,1.084029,1.087105,1.095639,1.122232,1.126333,1.146405,1.156017,1.178718,1.19162,1.210628,1.253956,1.313625,1.37973,1.438995,1.538787,1.715444,2.050254,2.782814,4.561827,8.714542};
    int ndx;
    for( ndx = 1; ndx < nbins+1; ndx++ ) {
      if( ntrk < xAxis1[ndx] ) break;
    }
    ndx += -1;
    if ( ndx >= nbins ) ndx = nbins-1;
    return ntrk_weight[ndx];

  } else if ( weight_choice == 41) { // tuneD6T, 7 TeV, with primary vertex
    const int nbins = 75;
    Double_t xAxis1[nbins+1] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 84, 88, 92, 96, 100, 110, 120, 140, 160, 200}; 
    // 358 results   const Double_t ntrk_weight[nbins] = {0.6312897, 0.6312897, 0.6312897, 0.6253497, 0.6658684, 0.7294919, 0.8098004, 0.8764595, 0.9188986, 0.9345537, 0.9242891, 0.9012419, 0.8726596, 0.8504198, 0.8333755, 0.8186717, 0.8168494, 0.8206679, 0.828048, 0.8415845, 0.8565698, 0.8766649, 0.8987519, 0.9157242, 0.939483, 0.9614064, 0.9815759, 1.004863, 1.031907, 1.053937, 1.078441, 1.113513, 1.136552, 1.168183, 1.203721, 1.229559, 1.268941, 1.302564, 1.334454, 1.364148, 1.395204, 1.43777, 1.457255, 1.491248, 1.528079, 1.560139, 1.582776, 1.616776, 1.654698, 1.673776, 1.707374, 1.763208, 1.798846, 1.824169, 1.856153, 1.884028, 1.890027, 1.915008, 1.905394, 1.912814, 1.92458, 1.921458, 1.918613, 1.932631, 1.932532, 1.94503, 1.978028, 2.038563, 2.129746, 2.289332, 2.617837, 3.585295, 6.558398, 21.86197, 128.8731};
    const Double_t ntrk_weight[nbins] = {0.5346584,0.5346584,0.5346584,0.5568761,0.6095921,0.6816129,0.7660205,0.842346,0.8989263,0.9193846,0.9182254,0.9016897,0.8765742,0.8538699,0.837409,0.8270569,0.8238369,0.8283648,0.8358256,0.8513013,0.8726104,0.8874779,0.9100027,0.9318814,0.9573484,0.9804721,1.00253,1.023539,1.054297,1.081563,1.105746,1.131595,1.161169,1.204135,1.227576,1.259902,1.300146,1.329525,1.363973,1.400137,1.438247,1.470389,1.500601,1.528155,1.576439,1.609052,1.637118,1.666232,1.68735,1.722651,1.754339,1.807355,1.850101,1.876932,1.906338,1.928337,1.940794,1.950954,1.953578,1.956098,1.960988,1.952343,1.957513,1.941629,1.957735,1.942944,1.981562,2.022363,2.094916,2.219302,2.52744,3.422568,6.168799,20.66746,128.7042};

    int ndx;
    for( ndx = 1; ndx < nbins+1; ndx++ ) {
      if( ntrk < xAxis1[ndx] ) break;
    }
    ndx += -1;
    if ( ndx >= nbins ) ndx = nbins-1;
    return ntrk_weight[ndx];

  } else if ( weight_choice == 42) { // pythia 8, 7 TeV, with primary vertex
    const int nbins = 75;
    Double_t xAxis1[nbins+1] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 84, 88, 92, 96, 100, 110, 120, 140, 160, 200}; 
    // 358 results:   const Double_t ntrk_weight[nbins] = {2.880105, 2.880105, 2.880105, 2.764258, 2.340879, 1.875745, 1.520436, 1.264521, 1.097916, 0.9896734, 0.9144784, 0.8653475, 0.8321328, 0.8135421, 0.808488, 0.8013878, 0.8105248, 0.8233051, 0.8336208, 0.8513218, 0.8631456, 0.8780952, 0.9041756, 0.9101756, 0.9304845, 0.9438393, 0.9505529, 0.9695289, 0.976975, 0.9879417, 0.9956397, 0.9993244, 1.004792, 1.011292, 1.028647, 1.018952, 1.03085, 1.027771, 1.034154, 1.034367, 1.04629, 1.046516, 1.046808, 1.043099, 1.049885, 1.047694, 1.053975, 1.055165, 1.067647, 1.063897, 1.066069, 1.075896, 1.078819, 1.07949, 1.082706, 1.09668, 1.100897, 1.120641, 1.126834, 1.142747, 1.174013, 1.170221, 1.196914, 1.219088, 1.246631, 1.279601, 1.339855, 1.398506, 1.495295, 1.587536, 1.775251, 2.169005, 2.898582, 4.783298, 9.599502};
    const Double_t ntrk_weight[nbins] = {1.97617,1.97617,1.97617,2.137671,1.978618,1.706341,1.437364,1.219161,1.081008,0.9694974,0.8990903,0.8572648,0.8245798,0.8078692,0.8003042,0.801486,0.8063687,0.8190495,0.8339258,0.8554677,0.8733004,0.8864059,0.9078891,0.9264084,0.9499294,0.9592131,0.9792666,0.9804129,0.9995499,1.008208,1.015627,1.011362,1.029179,1.039629,1.045962,1.048591,1.052184,1.053516,1.061246,1.065507,1.072357,1.072823,1.078528,1.071866,1.076624,1.090248,1.082533,1.085281,1.091185,1.089351,1.097825,1.099446,1.106978,1.109854,1.122615,1.125031,1.128048,1.136907,1.16485,1.169159,1.189955,1.200122,1.223606,1.237079,1.256728,1.301933,1.363872,1.432661,1.493925,1.597519,1.780907,2.128634,2.889124,4.7361,9.04746};
    int ndx;
    for( ndx = 1; ndx < nbins+1; ndx++ ) {
      if( ntrk < xAxis1[ndx] ) break;
    }
    ndx += -1;
    if ( ndx >= nbins ) ndx = nbins-1;
    return ntrk_weight[ndx];

  }
  
  return 1.0;


}

double seleXi::xipt_wgt( double pt, int weight_choice ) {

  if( !weight_choice ) return 1.;

  else if( weight_choice == 1 ) {//D6T - updated for 362
    const int nbins = 19;
    const Double_t xAxis1[nbins+1] = {0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 6};
    //const Double_t pt_weight[nbins] = {2.670758,1.685691, 1.731396, 1.767882, 1.622466, 1.380356, 1.241386, 1.073554, 1.009175, 0.9572005, 0.8854983, 0.7768946, 0.637305, 0.6001311, 0.665923, 0.6931358, 0.5890431, 0.592141, 0.5389072};
    const Double_t pt_weight[nbins] = {1.685691, 1.685691, 1.731396, 1.767882, 1.622466, 1.380356, 1.241386, 1.073554, 1.009175, 0.9572005, 0.8854983, 0.7768946, 0.637305, 0.6001311, 0.665923, 0.6931358, 0.5890431, 0.592141, 0.5389072};
    int ndx;
    for( ndx = 1; ndx < nbins+1; ndx++ ) {
      if( pt < xAxis1[ndx] ) break;
    }
    ndx += -1;
    if ( ndx >= nbins ) ndx = nbins-1;
    return pt_weight[ndx];
  }

  else if( weight_choice == 2 ) {//pythia8 - updated for 362
    const int nbins = 19;
    Double_t xAxis1[nbins+1] = {0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 6}; 
    //const Double_t pt_weight[nbins] = {1.270543,0.9315437,0.9560547,1.024832,1.230498,1.116464,0.9981943,1.14744,1.068846,1.074101,1.09059,0.9112481,0.8775955,0.9163617,0.7983331,0.6481371,0.8093392,1.276358,0.6881005};
    const Double_t pt_weight[nbins] = {0.9315437,0.9315437,0.9560547,1.024832,1.230498,1.116464,0.9981943,1.14744,1.068846,1.074101,1.09059,0.9112481,0.8775955,0.9163617,0.7983331,0.6481371,0.8093392,1.276358,0.6881005};
    int ndx;
    for( ndx = 1; ndx < nbins+1; ndx++ ) {
      if( pt < xAxis1[ndx] ) break;
    }
    ndx += -1;
    if ( ndx >= nbins ) ndx = nbins-1;
    return pt_weight[ndx];
  }

  else if( weight_choice == 3 ) {//P0 - updated for 362
    const int nbins = 19;
    Double_t xAxis1[nbins+1] = {0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 6}; 
    //const Double_t pt_weight[nbins] = {1.519527, 0.8600993, 0.8628748, 0.9376627, 0.9961686, 0.9873122, 1.069006, 1.231733, 1.124576, 1.07279, 1.532467, 0.9103652, 0.9224715, 0.8550856, 1.080019, 0.8233612, 0.8688515, 0.7712246, 0.6567522};
    const Double_t pt_weight[nbins] = {0.8600993, 0.8600993, 0.8628748, 0.9376627, 0.9961686, 0.9873122, 1.069006, 1.231733, 1.124576, 1.07279, 1.532467, 0.9103652, 0.9224715, 0.8550856, 1.080019, 0.8233612, 0.8688515, 0.7712246, 0.6567522};
    int ndx;
    for( ndx = 1; ndx < nbins+1; ndx++ ) {
      if( pt < xAxis1[ndx] ) break;
    }
    ndx += -1;
    if ( ndx >= nbins ) ndx = nbins-1;
    return pt_weight[ndx];
  }
  else if( weight_choice == 11 ) {//900 GeV D6T not updated with new bins!
    const int nbins = 10;
    const Double_t xAxis1[nbins+1] = {1, 1.4, 1.6, 1.8, 2, 2.2, 2.6, 3, 3.4, 4, 6};
    const Double_t pt_weight[nbins] = {1.246124, 0.9699563, 0.9253977, 1.141979, 0.734954, 1.162812, 0.9290018, 2.069333, 0.8191425, 0.2775709};
    int ndx;
    for( ndx = 1; ndx < nbins+1; ndx++ ) {
      if( pt < xAxis1[ndx] ) break;
    }
    ndx += -1;
    if ( ndx >= nbins ) ndx = nbins-1;
    return pt_weight[ndx];
  }
  else if( weight_choice == 14 ) {//900 mixed
    const int nbins = 18;
    const Double_t xAxis1[nbins+1] = {0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 6};
    const Double_t pt_weight[nbins] = {0.6583757, 0.862677, 0.8540292, 1.213436, 0.9447124, 1.083224, 1.134622, 0.9059398, 1.417759, 0.8883265, 1.614526, 1.127046, 0.9391207, 1.021772, 2.095197, 0.5425083, 1.259516, 0.7086425};

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

double seleXi::xiy_wgt( double y, int weight_choice ) {

  if( !weight_choice ) return 1.;

  else if( weight_choice == 1 ) {//D6T - updated for 362
    const int nbins = 10;
    const Double_t xAxis1[nbins+1] = {0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2};
    const Double_t y_weight[nbins] = {1.05773, 0.9023893, 0.98598, 0.969235, 1.051303, 1.10063, 1.090076, 1.032312, 0.9118068, 0.9854679};
    int ndx;
    for( ndx = 1; ndx < nbins+1; ndx++ ) {
      if( y < xAxis1[ndx] ) break;
    }
    ndx += -1;
    if ( ndx >= nbins ) ndx = nbins-1;
    return y_weight[ndx];
  }

  else if( weight_choice == 2 ) {//pythia8 - updated for 362
    const int nbins = 10;
    Double_t xAxis1[nbins+1] = {0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2}; 
    const Double_t y_weight[nbins] = {1.180628,0.9691568,0.7950849,0.9720812,0.8712323,1.192058,1.01109,1.072142,0.9948435,1.13731};
    int ndx;
    for( ndx = 1; ndx < nbins+1; ndx++ ) {
      if( y < xAxis1[ndx] ) break;
    }
    ndx += -1;
    if ( ndx >= nbins ) ndx = nbins-1;
    return y_weight[ndx];
  
  }

  else if( weight_choice == 3 ) {//P0 - updated for 362
    const int nbins = 10;
    Double_t xAxis1[nbins+1] = {0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2}; 
    const Double_t y_weight[nbins] = {1.04669, 1.138241, 0.9094625, 1.010044, 0.9789823, 1.109949, 0.8777753, 0.9280546, 1.101257, 0.9248295};
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
    const Double_t y_weight[nbins] = {0.8902479, 1.27547, 1.162203, 0.947529, 0.9979188, 1.008283, 0.9632742, 0.7555642, 1.002395, 1.05369};
    int ndx;
    for( ndx = 1; ndx < nbins+1; ndx++ ) {
      if( y < xAxis1[ndx] ) break;
    }
    ndx += -1;
    if ( ndx >= nbins ) ndx = nbins-1;
    return y_weight[ndx];
  }
  
  else if( weight_choice == 14 ) {//900 GeV mixed - updated for 362 
    const int nbins = 11;
    const Double_t xAxis1[nbins+1] = {0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.3};
    const Double_t y_weight[nbins] = {0.9902304, 1.051384, 0.978337, 1.039913, 0.9531046, 1.024396, 0.9074398, 0.7681192, 1.121877, 1.283838};
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

void seleXi::SetOutputDir(std::string theOutputDir){
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

float seleXi::GetDeltaR(double eta1, double phi1, double eta2, double phi2){
  float dEta = eta1 - eta2;
  float dPhi = phi1 - phi2;
  float dR = sqrt(dEta*dEta + dPhi*dPhi);
  return dR;
}

