#include "/home/xuq/Documents/HI/RpA/Quality/root_setting.h"
#include "/home/xuq/Documents/HI/RpA/TreeAna/produceandcheck/file.h"

const int Npoint=1000;
int ilist0=0, ilist1=6, ilist2=14;
double JetIDcut[2];
double xrange_JetIDcut[2];
int ieta=0;
TH1D* makehisto(int ilist){
TString JetID;
TString JetIDName = JetIDNameList[ilist];
double binbound_JetID[100];
int ilistw;
if(ilist== 0)  {    JetID = "Max p_{T}^{charged}"; JetIDcut[0]=4; JetIDcut[1]=300; ilistw=0;}
if(ilist== 6)  {    JetID = "Max p_{T}^{charged}/p_{T}^{jet}"; JetIDcut[0]=0.05; JetIDcut[1]=1; ilistw=1;}
else if(ilist== 7)  {    JetID = "#Sum p_{T}^{charged}/p_{T}^{jet}"; JetIDcut[0] = 0; JetIDcut[1]=0.6;}
else if(ilist== 8)  {   JetID = "Max p_{T}^{neutral}/p_{T}^{jet}"; JetIDcut[0] = 0;  JetIDcut[1]=0.08;}
else if(ilist== 9)  {    JetID = "#Sum p_{T}^{neutral}/p_{T}^{jet}"; JetIDcut[0] = 0; JetIDcut[1]=0.15;}
else if(ilist== 10)  {   JetID = "Max p_{T}^{photons}/p_{T}^{jet}";     JetIDcut[0] = 0; JetIDcut[1]=0.50;}
else if(ilist== 11) {    JetID = "#Sum p_{T}^{photons}/p_{T}^{jet}";     JetIDcut[0] = 0; JetIDcut[1]=0.3;}
else if(ilist== 12) {    JetID = "#frac{#Sigma p_{T}^{charged}+#Sigma p_{T}^{photons}+#Sigma p_{T}^{neutral}+#Sigma p_{T}^{#mu}+#Sigma p_{T}^{e}}{p_{T}^{jet}}";JetIDcut[0]=0; JetIDcut[1]=1.01; ilistw=3;}
else if(ilist== 13) {    JetID = "#frac{#Sigma p_{T}^{charged}+#Sigma p_{T}^{photons}+#Sigma p_{T}^{neutral}+#Sigma p_{T}^{#mu}+#Sigma p_{T}^{e}}{p_{T}^{raw}}"; JetIDcut[0]=0;JetIDcut[1]=1.13;}
else if(ilist== 14)   {  JetID = "Max p_{T}^{neutral} /Max(#Sigma p_{T}^{charged},#Sigma p_{T}^{neutral})";JetIDcut[0]=0;JetIDcut[1]=0.975;ilistw=2;}

if(ilist==0){
double binbound_JetID[]={0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,300};}
if(ilist==12 || ilist==13 || ilist==14 ){
double binbound_JetID[]={0,0.2,0.4,0.6,0.7,0.8,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1.0,1.02,1.04,1.06,1.1,1.15,1.2,1.3,1.4,1.6,1.8,2.};}
else{
double binbound_JetID[]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.};}
//double binbound_JetID[]={0,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,0.4,0.425,0.45,0.475,0.5};
int Nbin_JetID=sizeof(binbound_JetID)/sizeof(double)-1;
TString histonameID=Form("jetpt%sEtaBin%s",JetIDName.Data(),etabinnameswap[ieta].Data());
TH2F* hdata2F=(TH2F*)fDataPPbJetID->Get(histonameID);

xrange_JetIDcut[0]=JetIDcut[0]+1e-4;
xrange_JetIDcut[1]=JetIDcut[1]-1e-4;
TH1D* histo=(TH1D*)hdata2F->ProjectionX("histo",hdata2F->GetYaxis()->FindBin(xrange_JetIDcut[0]),hdata2F->GetYaxis()->FindBin(xrange_JetIDcut[1]));
histo=(TH1D*)histo->Rebin(Nbin_pt,"histo",binbound_pt);
normalizeByBinWidth(histo);
histo->SetMarkerSize(1.4);
if(ilist!=14)
histo->SetTitle(Form("%.2f < %s < %.2f",JetIDcut[0],JetID.Data(),JetIDcut[1]));
else
histo->SetTitle(Form("%.2f < %s < %.3f",JetIDcut[0],JetID.Data(),JetIDcut[1]));

return histo;
}

void DrawDiffJetID(){
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);
TH1F* hFrame=new TH1F("","",1000,0,1000);
c1 = new TCanvas("c1","",600,600);
makeMultiPanelCanvas(c1,1,1,-0.12,0.03,0.13,0.1,0.03);
hFrame->GetXaxis()->SetTitle("p_{T} [GeV/c]");
hFrame->GetXaxis()->SetNdivisions(505);
if(ieta!=0 && ieta!=6)
hFrame->GetXaxis()->SetLimits(47.5,445);
else
hFrame->GetXaxis()->SetLimits(47.5,330);

c1->cd(1);
fixedFontHist(hFrame,1.2,1.7);
hFrame->SetTitle("");
hFrame->GetXaxis()->SetTitle("p_{T}^{jet} [GeV/c]");
hFrame->GetYaxis()->SetTitle("Number of Events");
hFrame->GetYaxis()->SetRangeUser(2e-1,5e7);

hFrame->DrawCopy();
TH1D* histo0 = makehisto(ilist0);
TH1D* histo1 = makehisto(ilist1);
TH1D* histo2 = makehisto(ilist2);
TString JetID0 = histo0->GetTitle();
TString JetID1 = histo1->GetTitle();
//JetID0.Prepend("Cut0: ");
//JetID1.Prepend("Cut1: ");
TString JetID2 = histo2->GetTitle();
//JetID2.Prepend("Cut2: ");
histo0->SetMarkerStyle(20);
histo0->SetMarkerColor(1);
histo0->SetLineColor(1);
histo1->SetMarkerStyle(25);
histo1->SetMarkerColor(4);
histo1->SetLineColor(4);
histo2->SetMarkerStyle(33);
histo2->SetMarkerColor(2);
histo2->SetLineColor(2);
//histo0->Draw("E1same");
//histo1->Draw("E1same");
//histo2->Draw("E1same");
TLegend *leg1=new TLegend(0.25,0.70,0.90,0.90);
TLegend *leg2=new TLegend(0.18,0.70,0.85,0.95);
leg1->SetBorderSize(0);
leg2->SetBorderSize(0);
leg1->SetFillColor(0);
leg2->SetFillColor(0);
leg1->SetTextSize(0.028);
leg2->SetTextSize(0.028);
leg1->AddEntry(histo0,JetID0,"lp");
leg1->AddEntry(histo1,JetID1,"lp");
leg1->AddEntry(histo2,JetID2,"lp");
//leg1->Draw("same");
TLatex *T=new TLatex();
T->SetNDC();
T->SetTextAlign(12);
T->SetTextSize(0.04);
T->SetTextColor(1);
T->SetTextFont(42);
TH1D* ratio1=(TH1D*)histo1->Clone("ratio1");
TH1D* ratio2=(TH1D*)histo2->Clone("ratio2");
hFrame->GetYaxis()->SetTitle(Form("Ratio between different cut"));
hFrame->GetYaxis()->SetRangeUser(0.8,1.29);
hFrame->DrawCopy();

ratio1->Divide(histo0);
ratio2->Divide(histo0);
leg2->AddEntry(ratio1,Form("#frac{%s}{%s}",JetID1.Data(),JetID0.Data()),"lp");
leg2->AddEntry(ratio2,Form("#frac{%s}{%s}",JetID2.Data(),JetID0.Data()),"lp");
leg2->Draw("same");
 for(int ibin = 1; ibin<=ratio1->GetNbinsX(); ibin++){
	ratio1->SetBinError(ibin,ratio1->GetBinContent(ibin)*histo0->GetBinError(ibin)/histo0->GetBinContent(ibin));
	ratio2->SetBinError(ibin,ratio2->GetBinContent(ibin)*histo0->GetBinError(ibin)/histo0->GetBinContent(ibin));
       if(ibin>=ratio1->GetXaxis()->FindBin(45) && ibin<=ratio1->GetXaxis()->FindBin(500));
       }//end loop for ibin

ratio1->Draw("E1same");
ratio2->Draw("E1same");
TLine *l =new TLine(47.5,1,500,1);
l->SetLineStyle(2);
l->SetLineColor(1);
l->Draw("same");
T->SetTextSize(0.05);
T->DrawLatex(0.25,0.20,etastring[ieta]);
c1->Print(Form("/home/xuq/Documents/HI/RpA/TreeAna/JetID/DrawDiffJetID%s.png",etabinnameswap[ieta].Data()));
c1->Print(Form("/home/xuq/Documents/HI/RpA/TreeAna/JetID/DrawDiffJetID%s.pdf",etabinnameswap[ieta].Data()));
}
