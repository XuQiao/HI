#include "hiForest.h"
#include "utilities.h"
#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TH1D.h>
#include <TNtuple.h>
#include <iostream>
#include <TLegend.h>
#include <TLatex.h>
#include "TrackCorrector.C"

TH1D* divideHistosDiffBins(TH1D* h_Num, TH1D* h_Den);
void normalizeByBinWidth(TH1D *histo);
void fillHistograms(HiForest * c_minbias, JetData data, int & numev_minbias, float & numev_minbias_trigCorr, int & numev_0_14_minbias, float & numev_0_14_minbias_trigCorr, int & numev_14_22_minbias, float & numev_14_22_minbias_trigCorr, int & numev_22_32_minbias, float & numev_22_32_minbias_trigCorr, int & numev_32_X_minbias, TH1D *hPartPt_minbias, TH1D *hPartPt_minbias_trkCorr, TH1D *hPartPt_minbias_trkCorr_trigCorr, TH1D *hPartPt_minbias_trkCorr_trigCorr_smallerBins, TH1D *hPartPt_0_14_minbias, TH1D *hPartPt_0_14_minbias_trkCorr, TH1D *hPartPt_0_14_minbias_trkCorr_trigCorr, TH1D*hRecMultiplicity);
//------------------------

void CombineSpectra_minbias_EtaCM_P13_P18_NoOLDAlignmentRuns_NewCorrv4()
{

   gROOT->Reset();
//   gROOT->ProcessLine(".x rootlogonChristof.C");
   gROOT->ForceStyle();
   gStyle->SetPalette(1);
   gStyle->SetTitleYOffset(1.30);
   gStyle->SetPadRightMargin(0.16);

   bool doSave = true;
   TFile * f_output = new TFile("CombineSpectra_minbias_EtaCM_P13_P18_NoOLDAlignmentRuns_NewCorrv4.root","recreate");

   // Define the input file and HiForest
  // char *infName_1 = "root://eoscms//eos/cms/store/caf/user/yjlee/pPb2013/promptReco/PA2013_HiForest_PromptReco_KrisztianMB_JSonPPb_forestv84.root";
   char *infName_1 = "/mnt/hadoop/cms/store/user/ymao/hiForest/PA2013/PA2013_HiForest_PromptReco_KrisztianMB_JSonPPb_forestv84.root";
//   char *infName_2 = "root://eoscms//eos/cms/store/group/phys_heavyions/kjung/MinBiasUPCForest_v71/MergedForest_withCones_MinBiasUPC_v71_2.root";
   HiForest *c_minbias_1 = new HiForest(infName_1,"",cPPb);
//   HiForest *c_minbias_2 = new HiForest(infName_2,"",cPPb);
//   c_minbias_1->doTrackCorrections=1;
   c_minbias_1->InitTree();
   cout <<"Input file" << infName_1<<endl ;

   std::vector<Double_t> ptBins;
   const Double_t small = 1e-3;
   Double_t pt;
   for(pt =   0.  ; pt <   48.-small; pt +=  2./2.) ptBins.push_back(pt);
   for(pt =   48.; pt <   60.-small; pt +=  4./2. ) ptBins.push_back(pt); 
   for(pt =   60.; pt <   76.-small; pt +=  8./2. ) ptBins.push_back(pt); 
   for(pt =   76.; pt <   140.-small; pt +=  16./2. ) ptBins.push_back(pt); 
   for(pt =   140.; pt <   364.-small; pt +=  32./2. ) ptBins.push_back(pt); 
   ptBins.push_back(364.0);

   std::vector<Double_t> ptBins_smaller;
   Double_t pt_smaller;
   for(pt_smaller = 0.2; pt_smaller < 20.0-small; pt_smaller += 0.1 ) ptBins_smaller.push_back(pt_smaller);
   for(pt_smaller = 20.0; pt_smaller < 40.0-small; pt_smaller += 0.2 ) ptBins_smaller.push_back(pt_smaller);
   for(pt_smaller = 40.0; pt_smaller < 103.6-small; pt_smaller += 0.4 ) ptBins_smaller.push_back(pt_smaller);
   ptBins_smaller.push_back(103.6);

/*
   std::vector<Double_t> ptBins_part;
   Double_t pt_part;
   for(pt_part =   0.2  ; pt_part <   1.2-small; pt_part +=  0.05) ptBins_part.push_back(pt_part);
   for(pt_part =   1.2; pt_part <   2.4-small; pt_part +=  0.1 ) ptBins_part.push_back(pt_part); //12 bins
   for(pt_part =   2.4; pt_part <   7.2-small; pt_part +=  0.4 ) ptBins_part.push_back(pt_part); //12 bins
   for(pt_part =   7.2; pt_part <  16.8-small; pt_part +=  1.6 ) ptBins_part.push_back(pt_part);//it was 3.2
   for(pt_part =  16.8; pt_part <  30.0-small; pt_part +=  6.6 ) ptBins_part.push_back(pt_part);
   for(pt_part =  30.0; pt_part <  100.0-small; pt_part +=  17.5 ) ptBins_part.push_back(pt_part);
   for(pt_part =  100.0; pt_part <  300.0-small; pt_part += 25. ) ptBins_part.push_back(pt_part);
   ptBins_part.push_back(300.0);
*/
   //orig
//   Int_t NumOfPtBins_part = 34;
//   double ptBins_part[] = {0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 3.2, 4.0, 4.8, 5.6, 6.4, 7.2, 9.6, 12.0, 14.4, 19.2, 24.0, 28.8, 35.2, 41.6, 48.0, 60.8, 73.6, 86.4, 103.6};//, 120.8, 138., 155.2, 172.4};
   Int_t NumOfPtBins_part = 40;
   double ptBins_part[] = {0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 3.2, 4.0, 4.8, 5.6, 6.4, 7.2, 9.6, 12.0, 14.4, 19.2, 24.0, 28.8, 35.2, 41.6, 48.0, 60.8, 73.6, 86.4, 103.6, 120.8, 138., 155.2, 172.4, 189.6, 206.8};

   TTree * t = new TTree("t","whatever");
   JetData data(t,1);

   TH1D * hPartPt_minbias = new TH1D("hPartPt_minbias","hPartPt_minbias",NumOfPtBins_part-1,&ptBins_part[0]);
   TH1D * hPartPt_minbias_trkCorr = new TH1D("hPartPt_minbias_trkCorr","hPartPt_minbias_trkCorr",NumOfPtBins_part-1,&ptBins_part[0]);
   TH1D * hPartPt_minbias_trkCorr_trigCorr = new TH1D("hPartPt_minbias_trkCorr_trigCorr","hPartPt_minbias_trkCorr_trigCorr",NumOfPtBins_part-1,&ptBins_part[0]);
   TH1D * hPartPt_minbias_trkCorr_trigCorr_smallerBins = new TH1D("hPartPt_minbias_trkCorr_trigCorr_smallerBins","hPartPt_minbias_trkCorr_trigCorr_smallerBins",ptBins_smaller.size()-1,&ptBins_smaller[0]);
   TH1D * hPartPt_0_14_minbias = new TH1D("hPartPt_0_14_minbias","hPartPt_0_14_minbias",NumOfPtBins_part-1,&ptBins_part[0]);
   TH1D * hPartPt_0_14_minbias_trkCorr = new TH1D("hPartPt_0_14_minbias_trkCorr","hPartPt_0_14_minbias_trkCorr",NumOfPtBins_part-1,&ptBins_part[0]);
   TH1D * hPartPt_0_14_minbias_trkCorr_trigCorr = new TH1D("hPartPt_0_14_minbias_trkCorr_trigCorr","hPartPt_0_14_minbias_trkCorr_trigCorr",NumOfPtBins_part-1,&ptBins_part[0]);
   TH1D * hNumEv_minbias = new TH1D("hNumEv_minbias","hNumEv_minbias",1,0.,1.);
   TH1D * hNumEv_minbias_trigCorr = new TH1D("hNumEv_minbias_trigCorr","hNumEv_minbias_trigCorr",1,0.,1.);
   TH1D * hNumEv_0_14_minbias = new TH1D("hNumEv_0_14_minbias","hNumEv_0_14_minbias",1,0.,1.);
   TH1D * hNumEv_0_14_minbias_trigCorr = new TH1D("hNumEv_0_14_minbias_trigCorr","hNumEv_0_14_minbias_trigCorr",1,0.,1.);
   TH1D * hNumEv_14_22_minbias = new TH1D("hNumEv_14_22_minbias","hNumEv_14_22_minbias",1,0.,1.);
   TH1D * hNumEv_14_22_minbias_trigCorr = new TH1D("hNumEv_14_22_minbias_trigCorr","hNumEv_14_22_minbias_trigCorr",1,0.,1.);
   TH1D * hNumEv_22_32_minbias = new TH1D("hNumEv_22_32_minbias","hNumEv_22_32_minbias",1,0.,1.);
   TH1D * hNumEv_22_32_minbias_trigCorr = new TH1D("hNumEv_22_32_minbias_trigCorr","hNumEv_22_32_minbias_trigCorr",1,0.,1.);
   TH1D * hNumEv_32_X_minbias = new TH1D("hNumEv_32_X_minbias","hNumEv_32_X_minbias",1,0.,1.);

   TH1D * hRecMultiplicity = new TH1D("hRecMultiplicity","hRecMultiplicity",300.,0.,300.);

   hPartPt_minbias->Sumw2();
   hPartPt_minbias_trkCorr->Sumw2();
   hPartPt_minbias_trkCorr_trigCorr->Sumw2();
   hPartPt_minbias_trkCorr_trigCorr_smallerBins->Sumw2();
   hPartPt_0_14_minbias->Sumw2();
   hPartPt_0_14_minbias_trkCorr->Sumw2();
   hPartPt_0_14_minbias_trkCorr_trigCorr->Sumw2();

   hRecMultiplicity->Sumw2();

   int numev_minbias = 0;
   float numev_minbias_trigCorr = 0;
   int numev_0_14_minbias = 0;
   float numev_0_14_minbias_trigCorr = 0;
   int numev_14_22_minbias = 0;
   float numev_14_22_minbias_trigCorr = 0;
   int numev_22_32_minbias = 0;
   float numev_22_32_minbias_trigCorr = 0;
   int numev_32_X_minbias = 0;

   fillHistograms(c_minbias_1, data, numev_minbias, numev_minbias_trigCorr, numev_0_14_minbias, numev_0_14_minbias_trigCorr, numev_14_22_minbias, numev_14_22_minbias_trigCorr, numev_22_32_minbias, numev_22_32_minbias_trigCorr, numev_32_X_minbias, hPartPt_minbias, hPartPt_minbias_trkCorr, hPartPt_minbias_trkCorr_trigCorr, hPartPt_minbias_trkCorr_trigCorr_smallerBins, hPartPt_0_14_minbias, hPartPt_0_14_minbias_trkCorr, hPartPt_0_14_minbias_trkCorr_trigCorr,hRecMultiplicity);
//   fillHistograms(c_minbias_2, data, numev_minbias, numev_minbias_trigCorr, numev_0_14_minbias, numev_0_14_minbias_trigCorr, numev_14_22_minbias, numev_14_22_minbias_trigCorr, numev_22_32_minbias, numev_22_32_minbias_trigCorr, numev_32_X_minbias, hPartPt_minbias, hPartPt_minbias_trkCorr, hPartPt_minbias_trkCorr_trigCorr, hPartPt_minbias_trkCorr_trigCorr_smallerBins, hPartPt_0_14_minbias, hPartPt_0_14_minbias_trkCorr, hPartPt_0_14_minbias_trkCorr_trigCorr,hRecMultiplicity);

   hNumEv_minbias->SetBinContent(1,numev_minbias);
   hNumEv_minbias_trigCorr->SetBinContent(1,numev_minbias_trigCorr);
   hNumEv_0_14_minbias->SetBinContent(1,numev_0_14_minbias);
   hNumEv_0_14_minbias_trigCorr->SetBinContent(1,numev_0_14_minbias_trigCorr);
   hNumEv_14_22_minbias->SetBinContent(1,numev_14_22_minbias);
   hNumEv_14_22_minbias_trigCorr->SetBinContent(1,numev_14_22_minbias_trigCorr);
   hNumEv_22_32_minbias->SetBinContent(1,numev_22_32_minbias);
   hNumEv_22_32_minbias_trigCorr->SetBinContent(1,numev_22_32_minbias_trigCorr);
   hNumEv_32_X_minbias->SetBinContent(1,numev_32_X_minbias);
   normalizeByBinWidth(hPartPt_minbias);
   normalizeByBinWidth(hPartPt_minbias_trkCorr);
   normalizeByBinWidth(hPartPt_minbias_trkCorr_trigCorr);
   normalizeByBinWidth(hPartPt_minbias_trkCorr_trigCorr_smallerBins);
   normalizeByBinWidth(hPartPt_0_14_minbias);
   normalizeByBinWidth(hPartPt_0_14_minbias_trkCorr);
   normalizeByBinWidth(hPartPt_0_14_minbias_trkCorr_trigCorr);

   f_output->cd();
   hPartPt_minbias->Write();
   hPartPt_minbias_trkCorr->Write();
   hPartPt_minbias_trkCorr_trigCorr->Write();
   hPartPt_minbias_trkCorr_trigCorr_smallerBins->Write();
   hPartPt_0_14_minbias->Write();
   hPartPt_0_14_minbias_trkCorr->Write();
   hPartPt_0_14_minbias_trkCorr_trigCorr->Write();
   hNumEv_minbias->Write();
   hNumEv_minbias_trigCorr->Write();
   hNumEv_0_14_minbias->Write();
   hNumEv_0_14_minbias_trigCorr->Write();
   hNumEv_14_22_minbias->Write();
   hNumEv_14_22_minbias_trigCorr->Write();
   hNumEv_22_32_minbias->Write();
   hNumEv_22_32_minbias_trigCorr->Write();
   hNumEv_32_X_minbias->Write();
   hRecMultiplicity->Write();
   f_output->Close();
}

TH1D* divideHistosDiffBins(TH1D* h_Num, TH1D* h_Den) {
   TH1D *h_ratio = (TH1D*)h_Num->Clone("h_ratio");
   h_ratio->Reset();
   for(int i = 1; i <= h_Num->GetNbinsX(); i++) {
      float content = h_Num->GetBinContent(i);
      float error = h_Num->GetBinError(i);
      float center = h_Num->GetBinCenter(i);
      int which_bin_in_h_Den = h_Den->FindBin(center);
      float content_h_Den = h_Den->GetBinContent(which_bin_in_h_Den);

      if(content_h_Den==0)
         continue;

      h_ratio->SetBinContent(i,content/content_h_Den);      
      h_ratio->SetBinError(i,error/content_h_Den);      
   }
   return h_ratio;
}


void normalizeByBinWidth(TH1D *histo) {
   for(int i = 1; i <= histo->GetNbinsX(); i++) {
      float content = histo->GetBinContent(i);
      float error = histo->GetBinError(i);
      histo->SetBinContent(i,content/histo->GetBinWidth(i));
      histo->SetBinError(i,error/histo->GetBinWidth(i));
   }
}

void fillHistograms(HiForest * c_minbias, JetData data, int & numev_minbias, float & numev_minbias_trigCorr, int & numev_0_14_minbias, float & numev_0_14_minbias_trigCorr, int & numev_14_22_minbias, float & numev_14_22_minbias_trigCorr, int & numev_22_32_minbias, float & numev_22_32_minbias_trigCorr, int & numev_32_X_minbias, TH1D *hPartPt_minbias, TH1D *hPartPt_minbias_trkCorr, TH1D *hPartPt_minbias_trkCorr_trigCorr, TH1D *hPartPt_minbias_trkCorr_trigCorr_smallerBins, TH1D *hPartPt_0_14_minbias, TH1D *hPartPt_0_14_minbias_trkCorr, TH1D *hPartPt_0_14_minbias_trkCorr_trigCorr, TH1D* hRecMultiplicity) {

   TrackCorrector corr("trackCorrections_HIN12017v4_HijingCombined.root");
   corr.load("trkCorr_HIN12017");
   corr.setOption1(true);
   corr.setOption2(true);

   for (int i=0;i<c_minbias->GetEntries();i++) {
//   for (int i=0;i<500000;i++) {

      c_minbias->GetEntry(i);

      data.leadingJetPt = -1;
      data.leadingJetIt = -1;
      data.leadingJetTrackMax = -1;

      if (i % 2000 == 0) cout <<i<<" / "<<c_minbias->GetEntries()<<endl;
    //  if (i % 100 == 0) cout <<i<<" / "<<c_minbias->GetEntries()<<endl;

      if(c_minbias->evt.run<210676 || c_minbias->evt.run>211256) //211256: last pPb run (Pb goes to +eta)
         continue;

      //event selection
      bool event_accepted = true;
      if(!(c_minbias->skim.phfPosFilter1 && c_minbias->skim.phfNegFilter1
          && c_minbias->skim.pBeamScrapingFilter
          && c_minbias->skim.pprimaryvertexFilter
          && c_minbias->skim.pVertexFilterCutGplus
          && TMath::Abs(c_minbias->evt.vz)<15.
          )
        ) event_accepted = false;

      if(event_accepted == false)
         continue;

      //Jets
      for (int j=0;j<c_minbias->akPu3PF.nref;j++) {
         if (fabs(c_minbias->akPu3PF.jteta[j])>2.5) continue;
//         if (fabs(c_minbias->akPu3PF.jtpt[j])<30) continue;
//         if (c_minbias->akPu3PF.trackMax[j]<4) continue; //Do we need this?
         if (c_minbias->akPu3PF.rawpt[j]<15) continue;
         if (c_minbias->akPu3PF.jtpt[j]>data.leadingJetPt) {
            data.leadingJetPt = c_minbias->akPu3PF.jtpt[j];
            data.leadingJetEta = c_minbias->akPu3PF.jteta[j];
            data.leadingJetPhi = c_minbias->akPu3PF.jtphi[j];
            data.leadingJetTrackMax = c_minbias->akPu3PF.trackMax[j];
            data.leadingJetIt = j;
         }
      }

      //Tracks for event weights; eta-pt cut removed
      int trackMult = 0;
      for(int j=0;j<c_minbias->track.nTrk;j++) {
         if(!((c_minbias->track.highPurity[j])
             && (fabs(c_minbias->track.trkDz1[j]/c_minbias->track.trkDzError1[j])<3)
             && (fabs(c_minbias->track.trkDxy1[j]/c_minbias->track.trkDxyError1[j])<3)
             && (c_minbias->track.trkPtError[j]/c_minbias->track.trkPt[j]<0.1)
            ))
            continue;            
         trackMult++;
      }
      // Don't analyze 0 multiplicity events; correction added later
      if(trackMult==0)
         continue;
      double evtWeight = 1.;
      evtWeight = corr.getEventWeight(trackMult);

//      hRecMultiplicity->Fill(trackMult);

      //Tracks for event classification: same as for analysis
      bool trackAbove14 = false;
      bool trackAbove22 = false;
      bool trackAbove32 = false;
      double rapidityShift = 0.465; //eta_lab = eta_cm - y_NN
      for(int j=0;j<c_minbias->track.nTrk;j++) {
         if(!((c_minbias->track.trkPt[j]>0.2) 
             && (c_minbias->track.trkEta[j]<(1.8-rapidityShift))
             && (c_minbias->track.trkEta[j]>(1.3-rapidityShift))
             && (c_minbias->track.highPurity[j])
             && (fabs(c_minbias->track.trkDz1[j]/c_minbias->track.trkDzError1[j])<3)
             && (fabs(c_minbias->track.trkDxy1[j]/c_minbias->track.trkDxyError1[j])<3)
             && (c_minbias->track.trkPtError[j]/c_minbias->track.trkPt[j]<0.1)
           ))
            continue;            
           if(c_minbias->track.trkPt[j]>14.)
              trackAbove14 = true;
           if(c_minbias->track.trkPt[j]>22.)
              trackAbove22 = true;
           if(c_minbias->track.trkPt[j]>32.)
              trackAbove32 = true;
      }

      //Tracks for analysis
      for(int j=0;j<c_minbias->track.nTrk;j++) {
         if(!((c_minbias->track.trkPt[j]>0.2) 
             && (c_minbias->track.trkEta[j]<(1.8-rapidityShift))
             && (c_minbias->track.trkEta[j]>(1.3-rapidityShift))
             && (c_minbias->track.highPurity[j])
             && (fabs(c_minbias->track.trkDz1[j]/c_minbias->track.trkDzError1[j])<3)
             && (fabs(c_minbias->track.trkDxy1[j]/c_minbias->track.trkDxyError1[j])<3)
             && (c_minbias->track.trkPtError[j]/c_minbias->track.trkPt[j]<0.1)
           ))
            continue;            

         if(data.leadingJetPt==-1) //for corrections
            data.leadingJetPt=10.;
         double trkWeight = 1.;
         trkWeight = corr.getWeight(c_minbias->track.trkPt[j],c_minbias->track.trkEta[j],data.leadingJetPt);

         hPartPt_minbias->Fill(c_minbias->track.trkPt[j]);
         hPartPt_minbias_trkCorr->Fill(c_minbias->track.trkPt[j],trkWeight);
         hPartPt_minbias_trkCorr_trigCorr->Fill(c_minbias->track.trkPt[j],trkWeight*evtWeight);
         hPartPt_minbias_trkCorr_trigCorr_smallerBins->Fill(c_minbias->track.trkPt[j],trkWeight*evtWeight);
         if(!trackAbove14) {
            hPartPt_0_14_minbias->Fill(c_minbias->track.trkPt[j]);
            hPartPt_0_14_minbias_trkCorr->Fill(c_minbias->track.trkPt[j],trkWeight);
            hPartPt_0_14_minbias_trkCorr_trigCorr->Fill(c_minbias->track.trkPt[j],trkWeight*evtWeight);
         }
      } //Tracks end

      // event numbers
      numev_minbias++;
      numev_minbias_trigCorr += evtWeight;
      if(!trackAbove14) {
         numev_0_14_minbias++;
         numev_0_14_minbias_trigCorr += evtWeight;
      }
      if(trackAbove14 && !trackAbove22) {
         numev_14_22_minbias++;
         numev_14_22_minbias_trigCorr += evtWeight;
      }
      if(trackAbove22 && !trackAbove32) {
         numev_22_32_minbias++;
         numev_22_32_minbias_trigCorr += evtWeight;
      }
      if(trackAbove32)
         numev_32_X_minbias++;

   }//event
   double zeroMultFraction = 1.;
   zeroMultFraction = corr.getZeroMultFrac();
   numev_minbias_trigCorr /= (1.-zeroMultFraction);

   //approximation:
   numev_0_14_minbias_trigCorr /= (1.-zeroMultFraction);
}

