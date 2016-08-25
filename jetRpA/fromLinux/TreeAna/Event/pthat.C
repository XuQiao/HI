#include "/home/xuq/Documents/HI/RpA/Quality/root_setting.h"
#include "/home/xuq/Documents/HI/RpA/TreeAna/produceandcheck/file.h"


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

const int Nbin_pt=219;
double binbound_pt[220];
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
//cout<<binbound_pt[i]<<" ,";
}


TString histoname="Pthat";
TString histonameAf="PthatAfXw";
TH1F* histo1=(TH1F*)fMCxSecPPb->Get(histoname);
TH1F* histo1Af=(TH1F*)fMCxSecPPb->Get(histonameAf);
TH1F* histo2=(TH1F*)fMCxSecPP->Get(histoname);
TH1F* histo2Af=(TH1F*)fMCxSecPP->Get(histonameAf);

histo1=(TH1F*)histo1->Rebin(Nbin_pt,histo1->GetName(),binbound_pt);	maniTH(histo1); 
histo1Af=(TH1F*)histo1Af->Rebin(Nbin_pt,histo1Af->GetName(),binbound_pt);	maniTH(histo1Af);       
histo2=(TH1F*)histo2->Rebin(Nbin_pt,histo2->GetName(),binbound_pt);	maniTH(histo2); 
histo2Af=(TH1F*)histo2Af->Rebin(Nbin_pt,histo2Af->GetName(),binbound_pt);	maniTH(histo2Af);

TCanvas* c1 = new TCanvas("c1"," ",600,600);
makeMultiPanelCanvas(c1,2,2,0,0,0.18,0.18,0.03);

TH1F* hFrame=new TH1F("","",1000,0,1000);
fixedFontHist(hFrame,2,2.2);
hFrame->SetTitle("");
hFrame->GetXaxis()->SetTitle("");
hFrame->GetXaxis()->SetNdivisions(505);
hFrame->GetYaxis()->SetTitle("# of Events");
hFrame->GetXaxis()->SetLimits(0,999);
hFrame->GetYaxis()->SetRangeUser(1.01e-2,1e7);
hFrame->GetXaxis()->SetRangeUser(0,500);
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
hFrame->GetYaxis()->SetRangeUser(1.01e-2,1e7);
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
//c1->Print("pthat.png");
c1->Print("pthat.pdf");
}
