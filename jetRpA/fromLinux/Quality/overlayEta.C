#include "root_setting.h"
#include "file.h"

void overlayEta(){
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);

if(ilist<6 || ilist>13){
cout<<"Wrong JetID variable!"<<endl;
exit();
}

if(ilist==12){
double binbound_JetID[]={0,0.2,0.4,0.6,0.7,0.8,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1.0,1.02,1.04,1.06,1.1,1.15,1.2,1.3,1.4,1.6,1.8,2.};}
else if(ilist==13){
double binbound_JetID[]={0,0.5,0.9,0.92,0.94,0.96,0.98,1.0,1.01,1.02,1.03,1.04,1.06,1.08,1.1,1.2,1.4,1.6,1.8,2.};}
else if(ilist==6){
double binbound_JetID[]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,1.};}
else{
double binbound_JetID[]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.};}
int Nbin_JetID=sizeof(binbound_JetID)/sizeof(double)-1;

const int netabin=8;
const double etabin[netabin]={-2.2,-1.2,-0.7,-0.3,0.3,0.7,1.2,2.2};

TFile *file=fdata;
TString histoname[netabin];
TH2F* histo2F[netabin]; TH1D* histo[netabin]; TH1D* rehisto[netabin];	TH1D* ratio[netabin-1];
c1 = new TCanvas("c1","",800,800);
makeMultiPanelCanvas(c1,1,2,0.12,0.1,0.12,0.1,0.03);

for(int i=0;i<netabin-1;i++){
histoname[i]=Form("jetpt%sEtaBin%.f_%.f_Cen0-100%%",JetIDName.Data(),etabin[i]*10,etabin[i+1]*10);

histo2F[i]=(TH2F*)file->Get(histoname[i]);
//TH1D* histo1=(TH1D*)histo12F->ProjectionX("histo1",histo12F->GetXaxis()->FindBin(binbound_pt[5]),histo12F->GetXaxis()->FindBin(binbound_pt[Nbin_pt]),"e");
histo[i]=(TH1D*)histo2F[i]->ProjectionY(Form("histo_%d",i),histo2F[i]->GetXaxis()->FindBin(binbound_pt[5]),histo2F[i]->GetXaxis()->FindBin(binbound_pt[Nbin_pt]),"e");
//histo[i]=histo2F[i]->ProjectionY(Form("histo_%d",i));
rehisto[i]=(TH1D*)histo[i]->Clone(Form("rehisto_%d",i));
rehisto[i]=(TH1D*)histo[i]->Rebin(Nbin_JetID,Form("rehisto_%d",i),binbound_JetID);
normalizeByBinWidth(rehisto[i]);
rehisto[i]->Scale(1/rehisto[i]->Integral());
rehisto[i]->GetYaxis()->SetTitle("Event Fraction");
fixedFontHist(rehisto[i],1.5,2);
}

c1->cd(1)->SetLogy();
int j=1,k=3,l=5,m=4;
rehisto[j]->GetXaxis()->SetTitle(JetID.Data());
rehisto[j]->GetXaxis()->SetRangeUser(0.85,1.2);
rehisto[j]->SetMarkerStyle(20);
rehisto[j]->SetMarkerSize(1.2);
rehisto[k]->SetMarkerSize(1.2);
rehisto[l]->SetMarkerSize(1.2);
rehisto[m]->SetMarkerSize(1.2);
rehisto[k]->SetMarkerStyle(24);
rehisto[l]->SetMarkerStyle(29);
rehisto[m]->SetMarkerStyle(21);
rehisto[j]->SetMarkerColor(1);
rehisto[k]->SetMarkerColor(2);
rehisto[m]->SetMarkerColor(3);
rehisto[l]->SetMarkerColor(4);
rehisto[j]->SetMinimum(1e-8);
rehisto[j]->SetMaximum(1);
rehisto[j]->SetTitle("");
rehisto[j]->Draw("E1");
rehisto[k]->Draw("E1same");
rehisto[l]->Draw("E1same");
rehisto[m]->Draw("E1same");

TLegend *leg1=new TLegend(0.75,0.80,0.9,0.95);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetTextSize(0.04);
leg1->AddEntry(rehisto[j],Form("%.1f<#eta_{CM}<%.1f",etabin[j],etabin[j+1]),"lp");
leg1->AddEntry(rehisto[k],Form("%.1f<#eta_{CM}<%.1f",etabin[k],etabin[k+1]),"lp");
leg1->AddEntry(rehisto[l],Form("%.1f<#eta_{CM}<%.1f",etabin[l],etabin[l+1]),"lp");
leg1->AddEntry(rehisto[m],Form("%.1f<#eta_{CM}<%.1f",etabin[m],etabin[m+1]),"lp");
leg1->Draw("same");
TLatex *T=new TLatex(0.50,0.90,Form("p_{T}^{jet}(GeV/c): %.f-%.f",binbound_pt[5],binbound_pt[Nbin_pt]));
//TLatex *T=new TLatex(0.50,0.90,Form(""));
T->SetNDC();
T->SetTextAlign(12);
T->SetTextSize(0.05);
T->SetTextColor(1);
T->SetTextFont(42);
T->Draw("same");

c1->cd(2);
TH1D* ratio1=(TH1D*)rehisto[j]->Clone("ratio1");
TH1D* ratio2=(TH1D*)rehisto[l]->Clone("ratio2");
ratio1->Divide(rehisto[k]);
ratio2->Divide(rehisto[k]);
ratio1->SetMarkerStyle(20);
ratio1->SetMarkerSize(1.2);
ratio1->SetMarkerColor(1);
ratio1->SetLineColor(1);
ratio2->SetMarkerStyle(24);
ratio2->SetMarkerSize(1.2);
ratio2->SetMarkerColor(2);
ratio2->SetLineColor(2);
ratio1->SetMinimum(0);
ratio1->SetMaximum(2);
ratio1->SetTitle("");
ratio1->GetYaxis()->SetTitle("ratio");
ratio1->Draw("E1");
ratio2->Draw("E1same");
TLegend *leg2=new TLegend(0.70,0.80,0.85,0.95);
leg2->SetBorderSize(0);
leg2->SetFillColor(0);
leg2->SetTextSize(0.04);
leg2->AddEntry(ratio1,Form("%.1f<#eta_{CM}<%.1f/%.1f<#eta_{CM}<%.1f",etabin[j],etabin[j+1],etabin[k],etabin[k+1]),"lp");
leg2->AddEntry(ratio2,Form("%.1f<#eta_{CM}<%.1f/%.1f<#eta_{CM}<%.1f",etabin[l],etabin[l+1],etabin[k],etabin[k+1]),"lp");
leg2->Draw("same");
TLine *line =new TLine(binbound_JetID[0],1,binbound_JetID[Nbin_JetID],1);
//TLine *l =new TLine(binbound_pt[5],1,binbound_pt[Nbin_pt],1);
//TLine *l =new TLine(-3,1,3,1);
line->SetLineStyle(2);
line->SetLineColor(1);
line->Draw("same");

//c1->Print(Form("pic/JetIDcut/%s/overlay_Eta.png",JetIDName.Data()));


}
