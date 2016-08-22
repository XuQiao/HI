#include "/home/xuq/Documents/HI/RpA/Quality/root_setting.h"
#include "/home/xuq/Documents/HI/RpA/TreeAna/produceandcheck/file.h"

const double binbound_pt[]={ 3, 4, 5, 7, 9, 12, 15, 18,21,24,28,32,37,43,49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,1000};
int Nbin_pt=sizeof(binbound_pt)/sizeof(double)-1;
int ilist=0;//0,6,12,14
int ilist1=6, ilist2=14;
const double cut[15]={4,200,200,200,200,200,0.05,0.6,0.08,0.15,0.50,0.3,1.01,1.13,0.975};
bool isMC=0;
double JetIDcut[2];
double xrange_JetIDcut[2];
TH1F* makehisto(int ilist, int ieta, double cut){
TString JetID;
TString JetIDName = JetIDNameList[ilist];
int ilistw;
if(ilist== 0)  {    JetID = "Max p_{T}^{charged}"; JetIDcut[0]=cut; JetIDcut[1]=300; ilistw=0;}
if(ilist== 6)  {    JetID = "Max p_{T}^{charged}/p_{T}^{jet}"; JetIDcut[0]=cut; JetIDcut[1]=1; ilistw=1;}
else if(ilist== 7)  {    JetID = "#Sum p_{T}^{charged}/p_{T}^{jet}"; JetIDcut[0] = 0; JetIDcut[1]=cut;}
else if(ilist== 8)  {   JetID = "Max p_{T}^{neutral}/p_{T}^{jet}"; JetIDcut[0] = 0;  JetIDcut[1]=cut;}
else if(ilist== 9)  {    JetID = "#Sum p_{T}^{neutral}/p_{T}^{jet}"; JetIDcut[0] = 0; JetIDcut[1]=cut;}
else if(ilist== 10)  {   JetID = "Max p_{T}^{photons}/p_{T}^{jet}";     JetIDcut[0] = 0; JetIDcut[1]=cut;}
else if(ilist== 11) {    JetID = "#Sum p_{T}^{photons}/p_{T}^{jet}";     JetIDcut[0] = 0; JetIDcut[1]=cut;}
else if(ilist== 12) {    JetID = "#frac{#Sigma p_{T}^{charged}+#Sigma p_{T}^{photons}+#Sigma p_{T}^{neutral}+#Sigma p_{T}^{#mu}+#Sigma p_{T}^{e}}{p_{T}^{jet}}";JetIDcut[0]=0; JetIDcut[1]=cut; ilistw=3;}
else if(ilist== 13) {    JetID = "#Sigma p_{T}^{charged}+#Sigma p_{T}^{photons}+#Sigma p_{T}^{neutral}+#Sigma p_{T}^{#mu}+#Sigma p_{T}^{e}/p_{T}^{raw}"; JetIDcut[0]=0;JetIDcut[1]=cut;}
else if(ilist== 14)   {  JetID = "Max p_{T}^{neutral} /Max(#Sigma p_{T}^{charged},#Sigma p_{T}^{neutral})";JetIDcut[0]=0;JetIDcut[1]=cut;ilistw=2;}

if(ilist==0){
double binbound_JetID[]={0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,300};}

else if(ilist==12 || ilist==13 || ilist==14 ){
double binbound_JetID[]={0,0.2,0.4,0.6,0.7,0.8,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1.0,1.02,1.04,1.06,1.1,1.15,1.2,1.3,1.4,1.6,1.8,2.};}
else{
double binbound_JetID[]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.};}
//double binbound_JetID[]={0,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,0.4,0.425,0.45,0.475,0.5};
//int Nbin_JetID=sizeof(binbound_JetID)/sizeof(double)-1;
TString histonameID=Form("jetpt%sEtaBin%s",JetIDName.Data(),etabinnameswap[ieta].Data());
if(isMC) 
TH2F* histo2F=(TH2F*)fMCPPbJetID->Get(histonameID);
else
TH2F* histo2F=(TH2F*)fDataPPbJetID->Get(histonameID);

xrange_JetIDcut[0]=JetIDcut[0]+1e-5;
xrange_JetIDcut[1]=JetIDcut[1]-1e-5;

TH1F* histo=(TH1F*)histo2F->ProjectionX("histo",histo2F->GetYaxis()->FindBin(xrange_JetIDcut[0]),histo2F->GetYaxis()->FindBin(xrange_JetIDcut[1]));
histo=(TH1F*)histo->Rebin(Nbin_pt,"histo",binbound_pt);
//normalizeByBinWidth(histo);
histo->SetMarkerSize(1.4);
histo->SetTitle(Form("%s",JetID.Data()));

return histo;
}

void DrawJetIDsystematic()
{ 
  bool DoSave = kFALSE;
  double EtaMin[8] = {-2.0,-1.5,-1.0,-0.5,0.5,1.0,1.5,-1.0}; 
  double EtaMax[8] = {-1.5,-1.0,-0.5, 0.5,1.0,1.5,2.0, 1.0};
 
  double EtaMinRev[8];
  double EtaMaxRev[8];

  for(int j=0; j<8; j++){
     EtaMinRev[j] = -EtaMax[j];
     EtaMaxRev[j] = -EtaMin[j];
     }//swap eta convention

  int MarkerColor[8] = {1,2,4,6,7,8,40,46};

//  double ptbins[22] = {32,37,43,49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430}; 
  
//  TString TxtDir = "/afs/cern.ch/user/h/honi/CMSSW_5_3_8_HI_patch2/src/JetRpA0418/NewJetID/txtfile";
TLatex *T2=new TLatex();
T2->SetTextAlign(12);
T2->SetTextSize(0.03);
T2->SetTextColor(1);
T2->SetTextFont(42);
T2->SetNDC();

  double bincenter[8][21];
  double sys[8][21];

  TH1F *hJetID0[8];
  TH1F *hJetIDUp[8];
  TH1F *hJetIDDown[8];

  TH1F *hMax;
  TH1F *hMin;
  double MaxValueUp[100];
  double MinValueDown[100];
 
  for(int ieta = 0; ieta<8; ieta++){
if(ieta!=0 && ieta!=4) continue;
//    JetID[ieta].open(Form("%s/JetIDsystematic%.f_%.f.txt",TxtDir.Data(),10*EtaMin[ieta],10*EtaMax[ieta]));
//    hJetIDUp[ieta]=new TH1F(Form("hJetIDup%.f_%.f",10*EtaMin[ieta],10*EtaMax[ieta]),Form("hJetIDup%.f_%.f",10*EtaMin[ieta],10*EtaMax[ieta]),21,ptbins);
//    hJetIDDown[ieta]=new TH1F(Form("hJetIDdown%.f_%.f",10*EtaMin[ieta],10*EtaMax[ieta]),Form("hJetIDdown%.f_%.f",10*EtaMin[ieta],10*EtaMax[ieta]),21,ptbins);

hJetID0[ieta] = makehisto(ilist,ieta,cut[ilist]);
if(ilist==6){
double cutUp=cut[ilist]*0.8;
double cutDown=cut[ilist]*1.2;
}
else{
double cutUp=cut[ilist]*1.02;
double cutDown=cut[ilist]*0.98;
}
hJetIDDown[ieta] = makehisto(ilist2,ieta,cut[ilist2]);
hJetIDUp[ieta] = makehisto(ilist1,ieta,cut[ilist1]);

hJetIDDown[ieta]->Divide(hJetID0[ieta]);
hJetIDUp[ieta]->Divide(hJetID0[ieta]);

    for(int ibin = 1; ibin<=hJetID0[ieta]->GetNbinsX(); ibin++){
       hJetIDUp[ieta]->SetBinContent(ibin,TMath::Abs((hJetIDUp[ieta]->GetBinContent(ibin)-1))*100);
       hJetIDDown[ieta]->SetBinContent(ibin,-TMath::Abs((hJetIDDown[ieta]->GetBinContent(ibin)-1))*100);
       }//end loop for ibin
  }//end loop for ieta
 
  hMax = (TH1F *)hJetIDUp[0]->Clone(hJetIDUp[0]->GetName());
  hMin = (TH1F *)hJetIDDown[0]->Clone(hJetIDDown[0]->GetName());

  for(int ibin=0; ibin<=hJetID0[0]->GetNbinsX(); ibin++){
      MaxValueUp[ibin]=0;
      MinValueDown[ibin]=0;
      for(int ieta=0; ieta<8; ieta++){
if(ieta!=0 && ieta!=4) continue;
         if(MaxValueUp[ibin]<hJetIDUp[ieta]->GetBinContent(ibin)) MaxValueUp[ibin] = hJetIDUp[ieta]->GetBinContent(ibin);
         if(MinValueDown[ibin]>hJetIDDown[ieta]->GetBinContent(ibin)) MinValueDown[ibin] = hJetIDDown[ieta]->GetBinContent(ibin);
       }//find the max value for each bin
      hMax->SetBinContent(ibin,MaxValueUp[ibin]);
      hMin->SetBinContent(ibin,MinValueDown[ibin]);

     } 
    

  TCanvas *c1 = new TCanvas("c1","c1",0,0,650,650);
  makeMultiPanelCanvas(c1,1,1,0.0,0.0,0.12,0.01,0.06);
  gStyle->SetOptStat(0);
  c1->cd();
  TH1F *hDum = new TH1F("hDum","",1000,0,1000);
  hDum->SetMaximum(15.0);
  hDum->SetMinimum(-15.0);
  hDum->GetXaxis()->SetTitle("p_{T}[GeV/c]");
  hDum->GetXaxis()->SetRangeUser(47.5,500);
  hDum->GetXaxis()->SetTitleSize(0.045);
  hDum->GetXaxis()->SetLabelSize(0.035);
  hDum->GetXaxis()->CenterTitle();
  hDum->GetYaxis()->SetTitle("Relative Sys" );
  hDum->GetYaxis()->SetTitleSize(0.045);
  hDum->GetYaxis()->SetLabelSize(0.035);
  hDum->GetYaxis()->CenterTitle();
  hDum->DrawCopy();
  
  hMax->SetFillStyle(3001);
  hMax->SetFillColor(5);
  hMax->DrawCopy("sameHIST");

  hMin->SetFillStyle(3001);
  hMin->SetFillColor(5);
  hMin->DrawCopy("sameHIST");
  
  TLegend *t1 = new TLegend(0.25,0.72,0.55,0.90);
  t1->SetFillColor(0);
  t1->SetBorderSize(0);
  t1->SetFillStyle(0);
  t1->SetTextFont(63);
  t1->SetTextSize(17);  

  TLegend *t2 = new TLegend(0.55,0.72,0.85,0.90);
  t2->SetFillColor(0);
  t2->SetBorderSize(0);
  t2->SetFillStyle(0);
  t2->SetTextFont(63);
  t2->SetTextSize(17);


  for(int ieta=0; ieta<8; ieta++){
if(ieta!=0 && ieta!=4) continue;
     hJetIDUp[ieta]->SetMarkerSize(1.2);
     hJetIDUp[ieta]->SetLineWidth(2);
     hJetIDUp[ieta]->SetMarkerColor(MarkerColor[ieta]);
     hJetIDUp[ieta]->SetLineColor(MarkerColor[ieta]);
     hJetIDUp[ieta]->DrawCopy("sameHIST");

     hJetIDDown[ieta]->SetMarkerSize(1.2);
     hJetIDDown[ieta]->SetLineWidth(2);
     hJetIDDown[ieta]->SetMarkerColor(MarkerColor[ieta]);
     hJetIDDown[ieta]->SetLineColor(MarkerColor[ieta]);
     hJetIDDown[ieta]->DrawCopy("sameHIST");
    if(ieta<4) t1->AddEntry(hJetIDUp[ieta],Form("%.1f<#eta_{CM}<%.1f",EtaMinRev[ieta],EtaMaxRev[ieta]),"l");
    else  t2->AddEntry(hJetIDUp[ieta],Form("%.1f<#eta_{CM}<%.1f",EtaMinRev[ieta],EtaMaxRev[ieta]),"l");
   }
/*	if(ilist==6)
     T2->DrawLatex(0.2,0.2,Form("%.2f#pm%.2f<%s<1",cut[ilist],TMath::Abs(cutUp-cutDown)/2,hJetID0[0]->GetTitle()));
	else
     T2->DrawLatex(0.2,0.2,Form("0<%s<%.2f#pm%.2f",hJetID0[0]->GetTitle(),cut[ilist],TMath::Abs(cutUp-cutDown)/2));
*/
	
     T2->DrawLatex(0.2,0.15,Form("%.2f<%s<1",cut[ilist1],hJetID0[0]->GetTitle()));
     T2->DrawLatex(0.2,0.2,Form("0<%s<%.2f",hJetIDUp[0]->GetTitle(),cut[ilist2]));
     t1->Draw("same");  
     t2->Draw("same");
     gPad->RedrawAxis();

   if(DoSave){
if(isMC)  c1->SaveAs(Form("MC%sAllEtaBinSystematic.pdf",JetIDNameList[ilist].Data()));
else {c1->SaveAs(Form("Data%sAllEtaBinSystematic.pdf",JetIDNameList[ilist].Data()));
 c1->SaveAs(Form("Data%sAllEtaBinSystematic.png",JetIDNameList[ilist].Data()));
}
}
}
void makeMultiPanelCanvas(TCanvas*& canv,
                          const Int_t columns,
                          const Int_t rows,
                          const Float_t leftOffset,
                          const Float_t bottomOffset,
                          const Float_t leftMargin,
                          const Float_t bottomMargin,
                          const Float_t edge)
{
  if(canv == 0){
    Error("makeMultiPanelCanvas","GotNullCanvas.");
    return;
  }
  canv->Clear();

  TPad* pad[columns][rows];
  Float_t Xlow[columns];
   Float_t Xlow[columns];
  Float_t Xup[columns];
  Float_t Ylow[rows];
  Float_t Yup[rows];
 Float_t PadWidth =
  (1.0-leftOffset)/((1.0/(1.0-leftMargin)) +
  (1.0/(1.0-edge))+(Float_t)columns-2.0);
  Float_t PadHeight =
  (1.0-bottomOffset)/((1.0/(1.0-bottomMargin)) +
  (1.0/(1.0-edge))+(Float_t)rows-2.0);
  Xlow[0] = leftOffset;
  Xup[0] = leftOffset + PadWidth/(1.0-leftMargin);
  Xup[columns-1] = 1;
  Xlow[columns-1] = 1.0-PadWidth/(1.0-edge);
  Yup[0] = 1;
  Ylow[0] = 1.0-PadHeight/(1.0-edge);
  Ylow[rows-1] = bottomOffset;
  Yup[rows-1] = bottomOffset + PadHeight/(1.0-bottomMargin);

  for(Int_t i=1;i<columns-1;i++) {
    Xlow[i] = Xup[0] + (i-1)*PadWidth;
    Xup[i] = Xup[0] + (i)*PadWidth;
    }
  Int_t ct = 0;
  for(Int_t i=rows-2;i>0;i--) {
    Ylow[i] = Yup[rows-1] + ct*PadHeight;
    Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
    ct++;
  }
  TString padName;
  for(Int_t i=0;i<columns;i++) {
    for(Int_t j=0;j<rows;j++) {
   canv->cd();
      padName = Form("p_%d_%d",i,j);
      pad[i][j] = new TPad(padName.Data(),padName.Data(),
      Xlow[i],Ylow[j],Xup[i],Yup[j]);
      if(i==0) pad[i][j]->SetLeftMargin(leftMargin);
        else pad[i][j]->SetLeftMargin(0);
      if(i==(columns-1)) pad[i][j]->SetRightMargin(edge);
        else pad[i][j]->SetRightMargin(0);
      if(j==0) pad[i][j]->SetTopMargin(edge);
        else pad[i][j]->SetTopMargin(0);
      if(j==(rows-1)) pad[i][j]->SetBottomMargin(bottomMargin+0.3);
        else pad[i][j]->SetBottomMargin(0);
      pad[i][j]->Draw();
      pad[i][j]->cd();
      pad[i][j]->SetNumber(columns*j+i+1);
  }
 }
}

