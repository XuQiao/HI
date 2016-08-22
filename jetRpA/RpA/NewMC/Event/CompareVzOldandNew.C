#include "/home/xuq7/HI/jetRpA/RpA/Quality/root_setting.h"
#include "/home/xuq7/HI/jetRpA/RpA/NewMC/produceandcheck/file.h"

void CompareVzOldandNew(){
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);
TString coll="PPb";
bool Save=kFALSE;
 TF1 * fVz = new TF1("fVx","[0]+[1]*x+[2]*TMath::Power(x,2)+[3]*TMath::Power(x,3)+[4]*TMath::Power(x,4)", -15., 15.);
if(coll=="PPb"){
fVz->SetParameters(1.47442e+00, -2.83100e-03, -1.19295e-02, 1.10312e-05, 2.64814e-05); //! new official MC >QM14
//fVz->SetParameters(1.60182e+00,1.08425e-03,-1.29156e-02,-7.24899e-06,2.80750e-05);
//fVz->SetParameters(1.66731e+00,-2.43367e-03,-1.42488e-02,7.40147e-06,3.22477e-05);
TH1F* histodata=(TH1F*)fDataPPb->Get("Vz");
TH1F* histo1=(TH1F*)fMCPPb->Get("Vz");
TH1F* histo2=(TH1F*)fMCPPbOld->Get("Vz");
TString data="Pb going positive side";
}
else if(coll=="PbP"){
//fVz->SetParameters(1.54398e+00, -8.56155e-03, -1.40026e-02, 4.01020e-05, 3.47683e-05); //latest parameterization
fVz->SetParameters(1.49736e+00, -6.93060e-03, -1.26864e-02, 2.98693e-05, 2.89538e-05); //! new official MC after QM14 
TH1F* histodata=(TH1F*)fDataPbP->Get("Vz");
TH1F* histo1=(TH1F*)fMCPbP->Get("Vz");
TH1F* histo2=(TH1F*)fMCPbPOld->Get("Vz");
TString data="Proton going positive side";
}
else
fVz->SetParameters(1.,0,0,0,0);
histo1->SetName(Form("%sNewMC_unweighted",coll.Data()));
histo2->SetName(Form("%sOldMC_unweighted",coll.Data()));
histodata->SetName(Form("%sData",coll.Data()));
histo1->Rebin(5);
histo2->Rebin(5);
histodata->Rebin(5);
//histo1=(TH1F*)histo1->Rebin(Nbin_vz,"histo1",binbound_vz);
//histo2=(TH1F*)histo2->Rebin(Nbin_vz,"histo2",binbound_vz);
histo1->Scale(1/histo1->Integral());
histo2->Scale(1/histo2->Integral());
histodata->Scale(1/histodata->Integral());
histo1->SetMarkerStyle(24);
histo1->SetMarkerSize(1.2);
histo1->SetMarkerColor(2);
histo1->SetLineColor(2);
histo2->SetMarkerStyle(30);
histo2->SetMarkerColor(4);
histo2->SetMarkerSize(1.2);
//histo2->SetFillStyle(3004);
histo2->SetFillColor(2);
histo2->SetLineColor(2);
//histo2->SetMarkerColor(0);
//histo2->SetMarkerColor(4);
histodata->SetMarkerStyle(20);
histodata->SetMarkerSize(1.2);
histodata->SetMarkerColor(1);
histodata->SetLineColor(1);

TCanvas* c1 = new TCanvas("c1"," ",500,500);
//TCanvas* c2 = new TCanvas("c2"," ",500,500);
makeMultiPanelCanvas(c1,1,1,-0.16,0,0.16,0.14,0.03);
//makeMultiPanelCanvas(c2,1,1,-0.16,0,0.16,0.14,0.03);

TH1F* hFrame=new TH1F("","",400,-20,20);
fixedFontHist(hFrame,1.1,1.7);
hFrame->SetTitle("");
hFrame->GetXaxis()->SetTitleSize(0.05);
hFrame->GetYaxis()->SetTitleSize(0.05);
hFrame->GetXaxis()->SetTitle("Vz (cm)");
hFrame->GetYaxis()->SetTitle("Event Fraction");
hFrame->GetXaxis()->SetLimits(-15,15);
//hFrame->GetXaxis()->SetLimits(-3,3);
hFrame->GetYaxis()->SetRangeUser(0,5e-2);
c1->cd(1);
hFrame->DrawCopy();
histodata->Draw("same");
histo1->Draw("same");
histo2->Draw("same");
TLegend *leg1=new TLegend(0.30,0.80,0.85,0.92);
TLegend *leg2=new TLegend(0.30,0.85,0.85,0.95);
leg1->SetBorderSize(0);
leg2->SetBorderSize(0);
leg1->SetFillColor(0);
leg2->SetFillColor(0);
leg1->SetTextSize(0.04);
leg2->SetTextSize(0.04);
leg1->AddEntry(histo1,"Before Vz weighting New MC","lp");
leg1->AddEntry(histo2,"Before Vz weighting Old MC","lp");
//leg1->AddEntry(histo2,"New MC Before Vz weighting","lp");
leg1->AddEntry(histodata,data,"lp");
leg1->Draw("same");
TLatex *T1=new TLatex(0.25,0.92,"");
T1->SetNDC();
T1->SetTextAlign(12);
T1->SetTextSize(0.05);
T1->SetTextColor(1);
T1->SetTextFont(42);
T1->Draw("same");
/*c2->cd(1);
hFrame->GetYaxis()->SetTitle("Data/MC");
hFrame->GetXaxis()->SetLimits(-15,15);
hFrame->GetYaxis()->SetRangeUser(0,2);
hFrame->DrawCopy();
TH1F* ratio = (TH1F*)histodata->Clone(Form("%sratio",coll.Data()));*/
//ratio->SetName("ratio");
/*ratio->Divide(histo1);
ratio->Fit(fCen);
cout<<fCen->GetNDF()<<endl;
cout<<fCen->GetChisquare()<<endl;
for(int icent=0;icent<=6;icent++)
cout<<fCen->GetParameter(icent)<<",";
cout<<endl;*/
/*ratio->DrawCopy("same");
fVz->SetLineColor(2);
fVz->Draw("same");
leg2->AddEntry(ratio,data,"lp");
leg2->AddEntry(fVz,"reweighting function","lp");
leg2->Draw("same");
*/
if(Save){
TFile *fout = new TFile("VzandcentCompare.root","Update");
fout->cd();
ratio->Write("",TObject::kOverwrite);
histo1->Write("",TObject::kOverwrite);
histo2->Write("",TObject::kOverwrite);
histodata->Write("",TObject::kOverwrite);
fout->Close();
c2->Print(Form("ratio_Vz_%s.png",coll.Data()));
c1->Print(Form("VzCompare_Data%s.png",coll.Data()));
c1->Print(Form("VzCompare_Data%s.pdf",coll.Data()));
}
c1->Print(Form("VzCompareOldvsNew_Data%s.gif",coll.Data()));

}
