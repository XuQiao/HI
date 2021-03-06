#include "/home/xuq/Documents/HI/RpA/Quality/root_setting.h"
#include "/home/xuq/Documents/HI/RpA/TreeAna/produceandcheck/file.h"

void CompareCent(){
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);
TString coll="PbP";
bool Save=kTRUE;
TF1 * fCen = new TF1("fCen","[0]*exp([1]+[2]*x+[3]*x*x+[4]*x*x*x+[5]*x*x*x*x+[6]*x*x*x*x*x)", 0., 100.);
/*fCen->FixParameter(0,1.20916e-02);
fCen->FixParameter(1,5.02157e+00);
fCen->FixParameter(2, -3.38300e-02);
fCen->FixParameter(3, 1.87647e-03);
fCen->FixParameter(4, -6.76442e-05);
fCen->FixParameter(5, 9.08602e-07);
fCen->FixParameter(6,-4.01536e-09); //YX's */
if(coll=="PPb"){
fCen->SetParameters(1.20916e-02, 5.02157e+00, -3.38300e-02, 1.87647e-03, -6.76442e-05, 9.08602e-07, -4.01536e-09);//parameterize on 05.03 after approval
TH1F* histodata=(TH1F*)fDataPPb->Get("Cent");
TH1F* histo1=(TH1F*)fMCPPb->Get("Cent");
TH1F* histo2=(TH1F*)fMCPPb->Get("CentW");
TString data="Pb going positive side";
}
//fCen->SetParameters(8.68073e-03, 5.09356e+00, -1.33053e-02, 1.46904e-03, -6.99681e-05, 1.06721e-06, -5.21398e-09); //original
else if(coll=="PbP"){
//fCen->SetParameters(1.05408e-02, 5.27477e+00, -8.03382e-02, 3.51669e-03, -8.85332e-05, 1.08917e-06, -4.90091e-09);
fCen->SetParameters(1.14851e-02, 5.31172e+00, -8.52366e-02, 3.00268e-03, -6.04667e-05, 6.24105e-07, -2.43580e-09);
TH1F* histodata=(TH1F*)fDataPbP->Get("Cent");
TH1F* histo1=(TH1F*)fMCPbP->Get("Cent");
TH1F* histo2=(TH1F*)fMCPbP->Get("CentW");
TString data="Proton going positive side";
}
histo1->SetName("MC_unweighted");
histo2->SetName("MC_weighted");
histodata->SetName("Data");
//histo1=(TH1F*)histo1->Rebin(Nbin_vz,"histo1",binbound_vz);
//histo2=(TH1F*)histo2->Rebin(Nbin_vz,"histo2",binbound_vz);
histo1->Scale(1/histo1->Integral());
histo2->Scale(1/histo2->Integral());
histodata->Scale(1/histodata->Integral());
histo1->SetMarkerStyle(24);
histo1->SetMarkerSize(1.2);
histo1->SetMarkerColor(2);
histo1->SetLineColor(2);
histo2->SetMarkerStyle(0);
histo2->SetMarkerSize(0);
histo2->SetFillStyle(3004);
histo2->SetFillColor(2);
histo2->SetLineColor(2);
histo2->SetMarkerColor(0);
//histo2->SetMarkerColor(2);
histodata->SetMarkerStyle(20);
histodata->SetMarkerSize(1.2);
histodata->SetMarkerColor(1);
histodata->SetLineColor(1);

TCanvas* c1 = new TCanvas("c1"," ",500,500);
TCanvas* c2 = new TCanvas("c2"," ",500,500);
makeMultiPanelCanvas(c1,1,1,-0.1,0,0.16,0.14,0.03);
makeMultiPanelCanvas(c2,1,1,-0.1,0,0.16,0.14,0.03);

TH1F* hFrame=new TH1F("","",20000,-1000,1000);
fixedFontHist(hFrame,1.2,1.7);
hFrame->SetTitle("");
hFrame->GetXaxis()->SetTitleSize(0.05);
hFrame->GetYaxis()->SetTitleSize(0.05);
hFrame->GetXaxis()->SetTitle("Centrality");
hFrame->GetYaxis()->SetTitle("Event Fraction");
hFrame->GetXaxis()->SetLimits(0,100);
//hFrame->GetXaxis()->SetLimits(-3,3);
hFrame->GetYaxis()->SetRangeUser(0,5.2e-2);
c1->cd(1);
hFrame->DrawCopy();
histodata->Draw("same");
histo1->Draw("same");
histo2->Draw("HIST same");
TLegend *leg1=new TLegend(0.40,0.80,0.85,0.92);
TLegend *leg2=new TLegend(0.40,0.80,0.85,0.92);
leg1->SetBorderSize(0);
leg2->SetBorderSize(0);
leg1->SetFillColor(0);
leg2->SetFillColor(0);
leg1->SetTextSize(0.04);
leg2->SetTextSize(0.04);
leg1->AddEntry(histo1,"Before Centrality weighting","lp");
leg1->AddEntry(histo2,"After Centrality weighting","lfp");
leg1->AddEntry(histodata,data,"lp");
leg1->Draw("same");
TLatex *T1=new TLatex(0.25,0.92,"");
T1->SetNDC();
T1->SetTextAlign(12);
T1->SetTextSize(0.05);
T1->SetTextColor(1);
T1->SetTextFont(42);
T1->Draw("same");
c2->cd(1);
hFrame->GetYaxis()->SetTitle("Data/MC");
hFrame->GetXaxis()->SetLimits(0,100);
hFrame->GetYaxis()->SetRangeUser(0,2);
hFrame->DrawCopy();
TH1F* ratio = (TH1F*)histodata->Clone("ratio");
//ratio->SetName("ratio");
ratio->Divide(histo1);
/*ratio->Fit(fCen);
cout<<fCen->GetNDF()<<endl;
cout<<fCen->GetChisquare()<<endl;
for(int icent=0;icent<=6;icent++)
cout<<fCen->GetParameter(icent)<<",";
cout<<endl;*/
ratio->DrawCopy("same");
fCen->SetLineColor(2);
fCen->Draw("same");
leg2->AddEntry(ratio,data,"lp");
leg2->AddEntry(fCen,"reweighting function","lp");
leg2->Draw("same");
if(Save){
TFile *fout = new TFile("centCompare.root","RECREATE");
fout->cd();
ratio->Write();
histo1->Write();
histodata->Write();
fout->Close();
c2->Print(Form("ratio_Cent_%s.png",coll.Data()));
c1->Print(Form("CentCompare_Data%s.png",coll.Data()));
c1->Print(Form("CentCompare_Data%s.pdf",coll.Data()));
}

}
