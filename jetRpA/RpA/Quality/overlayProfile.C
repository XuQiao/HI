#include "file.h"
#include "root_setting.h"

void overlayProfile(){
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);

if(ilist== 0)   JetID = "chargedMax";
else if(ilist== 1)      JetID = "chargedSum";
else if(ilist== 2)      JetID = "neutralMax";
else if(ilist== 3) JetID = "neutralSum";
else if(ilist== 4) JetID = "photonMax";
else if(ilist== 5) JetID = "photonSum";
else if(ilist== 6)      JetID = "chargedMax/p_{T}";
else if(ilist== 7)      JetID = "chargedSum/p_{T}";
else if(ilist== 8)      JetID = "neutralMax/p_{T}";
else if(ilist== 9)      JetID = "neutralSum/p_{T}";
else if(ilist== 10)     JetID = "photonMax/p_{T}";
else if(ilist== 11)     JetID = "photonSum/p_{T}";
else if(ilist== 12)     JetID = "(#Sigma h^{#pm}+#Sigma #gamma +#Sigma h^{0}+#Sigma #mu+#Sigma e)/p_{T}^{jet}";
else if(ilist== 13)     JetID = "(#Sigma h^{#pm}+#Sigma #gamma +#Sigma h^{0}+#Sigma #mu+#Sigma e)/p_{T}^{raw}";
else if(ilist== 14)     JetID = "neutralMax/(neutralMax+chargedMax+photonMax)";
else if(ilist== 15)     JetID = "charged Multiplicity";
else if(ilist== 16)     JetID = "neutral Multiplicity";
else if(ilist== 17)     JetID = "photon Multiplicity";
else{   exit();}

TProfile* hdata=(TProfile*)hdata2F->ProfileX("hdata");
TProfile* hPPb=(TProfile*)hPPb2F->ProfileX("hPPb");
rehdata=(TProfile*)hdata->Clone("rehdata");
rehPPb=(TProfile*)hPPb->Clone("rehPPb");
rehdata=(TProfile*)hdata->Rebin(Nbin_pt,"rehdata",binbound_pt);
rehPPb=(TProfile*)hPPb->Rebin(Nbin_pt,"rehPPb",binbound_pt);
rehPPb->SetMarkerStyle(0);
rehPPb->SetMarkerColor(0);
rehPPb->SetLineColor(2);
rehPPb->SetFillStyle(3004);
rehPPb->SetFillColor(2);
rehdata->SetMarkerStyle(20);
rehdata->SetMarkerSize(1.2);
rehdata->SetMarkerColor(1);
for(int i=0;i<rehdata->GetNbinsX();i++){
cout<<rehdata->GetBinContent(i)<<endl;
}
TH1F* hFrame=new TH1F("","",1000,0,1000);
hFrame->GetXaxis()->SetRangeUser(26.9,692);
c1 = new TCanvas("c1","",600,600);
makeMultiPanelCanvas(c1,1,1,0,0,0.12,0.12,0.03);
hFrame->SetMaximum(TMath::Max(rehdata->GetMaximum(),rehPPb->GetMaximum())*1.2);
hFrame->SetMinimum(TMath::Min(rehdata->GetMinimum(),rehPPb->GetMinimum())/1.2);
hFrame->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
hFrame->GetYaxis()->SetTitle(Form("%s",JetID.Data()));
hFrame->SetTitle("");
fixedFontHist(hFrame,1.5,1.3);
hFrame->DrawCopy();
TLatex *T1=new TLatex(0.50,0.85,"");
T1->SetNDC();
T1->SetTextAlign(12);
T1->SetTextSize(0.05);
T1->SetTextColor(1);
T1->SetTextFont(42);
T1->Draw("same");
rehdata->Draw("E1same");
rehPPb->Draw("HIST E1same");
TLegend *leg1=new TLegend(0.60,0.2,0.80,0.35);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetTextSize(0.04);
leg1->AddEntry(rehdata,"Data","lp");
leg1->AddEntry(rehPPb,"PYTHIA+HIJING","lpf");
leg1->Draw("same");
c1->Print(Form("pic/JetIDcut/%s/%s/overlay_Profile.png",algo.Data(),JetIDName.Data()));

}
