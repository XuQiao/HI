#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"

const int Nbin_pt=219;
double binbound_pt[Nbin_pt];

//const double binbound_pt[]={3, 4, 5,6, 7, 9, 12, 15, 18, 22, 27, 33, 39, 47, 55, 64,74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 362, 429, 692, 1000};
//int Nbin_pt=sizeof(binbound_pt)/sizeof(double)-1;

TH1* maniTH(TH1* histo){
histo->SetMarkerStyle(20);
histo->SetMarkerSize(1);
histo->SetMarkerColor(1);
histo->SetLineColor(1);
normalizeByBinWidth(histo);
return histo;
}

void pthat(){
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);

//double binbound_pt[]={30,40,50,60,70,80,90,100,110,120,140,160,180,200,220,260,300,350,500};
//const double binbound_pt[]={0,5,10,15,20,30,45,60,75,90,105,120,140,160,180,200,220,260,300,400,600,1000};
for(int i=0;i<=Nbin_pt;i++){
if(i<100)
binbound_pt[i]=i;
else if(i<150)
binbound_pt[i]=100+(i-100)*2;
else if(i<180)
binbound_pt[i]=200+(i-150)*5;
else if(i<205)
binbound_pt[i]=350+(i-180)*10;
else if(i<215)
binbound_pt[i]=600+(i-205)*20;
else
binbound_pt[i]=800+(i-215)*50;
cout<<binbound_pt[i]<<" ,";
}


TString filename1="MCPPbNoGplus_ak3PFDiJetMassJetNoResidualPt30GenChargeTrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_merged.root";
TString filename2="MCPPNoGplus_ak3PFDiJetMassJetNoResidualPt30GenChargeTrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_merged.root";

TString dir1="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/ak3PF/NoResidual/NoGplusNewpthat";
TString dir2="~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/ak3PF/NoResidual/NoGplusNewpthat";

TString histoname="Pthat";
TString histonameAf="PthatAfXw";
//TString histoname1=Form("jetpt%s_0-100%%",JetIDName.Data());

TFile *file1=TFile::Open(Form("%s/%s",dir1.Data(),filename1.Data()));
TFile *file2=TFile::Open(Form("%s/%s",dir2.Data(),filename2.Data()));

TH1D* histo1=(TH1D*)file1->Get(histoname);
TH1D* histo1Af=(TH1D*)file1->Get(histonameAf);
TH1D* histo2=(TH1D*)file2->Get(histoname);
TH1D* histo2Af=(TH1D*)file2->Get(histonameAf);

histo1=(TH1D*)histo1->Rebin(Nbin_pt,histo1->GetName(),binbound_pt);	maniTH(histo1); 
histo1Af=(TH1D*)histo1Af->Rebin(Nbin_pt,histo1Af->GetName(),binbound_pt);	maniTH(histo1Af);       
histo2=(TH1D*)histo2->Rebin(Nbin_pt,histo2->GetName(),binbound_pt);	maniTH(histo2); 
histo2Af=(TH1D*)histo2Af->Rebin(Nbin_pt,histo2Af->GetName(),binbound_pt);	maniTH(histo2Af);

c1 = new TCanvas("c1"," ",600,600);
makeMultiPanelCanvas(c1,2,2,0,0,0.18,0.18,0.03);

TH1F* hFrame=new TH1F("","",1000,0,1000);
fixedFontHist(hFrame,2,2.2);
hFrame->SetTitle("");
hFrame->GetXaxis()->SetTitle("");
hFrame->GetXaxis()->SetNdivisions(505);
hFrame->GetYaxis()->SetTitle("# of Events");
hFrame->GetXaxis()->SetLimits(0,999);
hFrame->GetYaxis()->SetRangeUser(1.01e-2,1e6);
c1->cd(1)->SetLogy();
TLatex *T1=new TLatex(0.35,0.92,"Before cross section weighting");
T1->SetNDC();
T1->SetTextAlign(12);
T1->SetTextSize(0.05);
T1->SetTextColor(1);
T1->SetTextFont(42);
hFrame->GetXaxis()->SetTitle("");
hFrame->DrawCopy();
histo1->Draw("E1same");
T1->Draw("same");
TLegend *leg1=new TLegend(0.60,0.75,0.85,0.90);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetTextSize(0.05);
leg1->AddEntry(histo1,"PYTHIA+HIJING","lp");
leg1->Draw("same");
c1->cd(3)->SetLogy();
TLatex *T2=new TLatex(0.35,0.92,"After cross section weighting");
T2->SetNDC();
T2->SetTextAlign(12);
T2->SetTextSize(0.05);
T2->SetTextColor(1);
T2->SetTextFont(42);
hFrame->GetXaxis()->SetTitle("p_{T}^{hat} (GeV/c)");
hFrame->GetYaxis()->SetTitle("Event Probability");
hFrame->GetYaxis()->SetRangeUser(1e-13,0.99);
hFrame->DrawCopy();
histo1Af->Draw("E1same");
T2->Draw("same");
leg1->Draw("same");

c1->cd(2)->SetLogy();
hFrame->GetYaxis()->SetRangeUser(1.01e-2,1e6);
hFrame->GetXaxis()->SetTitle("");
hFrame->DrawCopy();
histo2->Draw("E1same");
TLegend *leg2=new TLegend(0.60,0.75,0.85,0.90);
leg2->SetBorderSize(0);
leg2->SetFillColor(0);
leg2->SetTextSize(0.05);
leg2->AddEntry(histo1,"PYTHIA","lp");
leg2->Draw("same");
T1->Draw("same");

c1->cd(4)->SetLogy();
hFrame->GetXaxis()->SetTitle("p_{T}^{hat} (GeV/c)");
hFrame->GetYaxis()->SetTitle("Event Probability");
hFrame->GetYaxis()->SetRangeUser(1e-13,0.99);
hFrame->DrawCopy();
histo2Af->Draw("E1same");
T2->Draw("same");
leg2->Draw("same");
c1->Print("pic/pthat.png");
c1->Print("pic/pthat.pdf");


}
