#include "/home/xuq/Documents/HI/RpA/Quality/root_setting.h"
#include "/home/xuq/Documents/HI/RpA/TreeAna/produceandcheck/file.h"

void P3overlayProfile(){
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);
TString svar="other";
TString coll="PbP";
int ieta=7;
if(svar=="Max"){
int listsvar[]={0,2,4,6,8,10};
double ymin1=0.1,ymax1=150;
double ymin2=0,ymax2=0.37;
TString svar1="Max (GeV/c)";
TString svar2="Max / p_{T}";}

else if(svar=="Sum"){
int listsvar[]={1,3,5,7,9,11};
double ymin1=0.1,ymax1=260;
double ymin2=0,ymax2=0.72;
TString svar1="Sum (GeV/c)";
TString svar2="Sum / p_{T}";}

else if(svar=="Mult"){
int listsvar[]={15,16,17};
double ymin1=0,ymax1=25;
TString svar1="Mult (GeV/c)";}

else if(svar=="other"){
int listsvar[]={12,13,14};
double ymin2=0,ymax2=1.25;
TString svar2="Jet ID Variables";}

else exit();

TH1F* hFrame=new TH1F("","",1000,0,1000);
int N=sizeof(listsvar)/sizeof(int);
if(N==6){
c1 = new TCanvas("c1"," ",1000,650);
makeMultiPanelCanvas(c1,3,2,0,0,0.22,0.16,0.03);
fixedFontHist(hFrame,1.8,2.2);
}
else if(N==3){
c1 = new TCanvas("c1"," ",1000,400);
makeMultiPanelCanvas(c1,3,1,0,0,0.16,0.14,0.03);
fixedFontHist(hFrame,1.2,1.5);
}
else exit();

for(int k=0;k<sizeof(listsvar)/sizeof(int);k++){
int ilist=listsvar[k];
c1->cd(k+1);
TString JetID,Unit;     TString JetIDName=JetIDNameList[ilist];

if(coll=="PPb"){
TString histoname=Form("jetpt%sEtaBin%s",JetIDName.Data(),etabinnameswap[ieta].Data());
TH2F* hdata2F=(TH2F*)fDataPPbJetID->Get(histoname);
TH2F* hPPb2F=(TH2F*)fMCPPbJetID->Get(histoname);
}
if(coll=="PbP"){
TString histoname=Form("jetpt%sEtaBin%s",JetIDName.Data(),etabinname[ieta].Data());
TH2F* hdata2F=(TH2F*)fDataPbPJetIDOld->Get(histoname);
TH2F* hPPb2F=(TH2F*)fMCPbPJetID->Get(histoname);
}


if(ilist== 0)   JetID = "Max p_{T}^{charged}";
else if(ilist== 1)      JetID = "#Sigma p_{T}^{charged}";
else if(ilist== 2)      JetID = "Max p_{T}^{neutral}";
else if(ilist== 3) JetID = "#Sigma p_{T}^{neutral}";
else if(ilist== 4) JetID = "Max p_{T}^{photon}";
else if(ilist== 5) JetID = "#Sigma p_{T}^{photon}";
else if(ilist== 6)      JetID = "Max p_{T}^{charged}/p_{T}^{jet}";
else if(ilist== 7)      JetID = "#Sigma p_{T}^{charged}/p_{T}^{jet}";
else if(ilist== 8)      JetID = "Max p_{T}^{neutral}/p_{T}^{jet}";
else if(ilist== 9)      JetID = "#Sigma p_{T}^{neutral}/p_{T}^{jet}";
else if(ilist== 10)     JetID = "Max p_{T}^{photon}/p_{T}^{jet}";
else if(ilist== 11)     JetID = "#Sigma p_{T}^{photon}/p_{T}^{jet}";
else if(ilist== 12)     JetID = "(#Sigma p_{T}^{charged}+#Sigma p_{T}^{photon} +#Sigma p_{T}^{neutral}+#Sigma p_{T}^{#mu}+#Sigma p_{T}^{e})/p_{T}^{jet}";
else if(ilist== 13)     JetID = "(#Sigma p_{T}^{charged}+#Sigma p_{T}^{photon} +#Sigma p_{T}^{neutral}+#Sigma p_{T}^{#mu}+#Sigma p_{T}^{e})/p_{T}^{raw}";
else if(ilist== 14)     JetID = "Max p_{T}^{neutral}/Max(#Sigma p_{T}^{charged}, #Sigma p_{T}^{neutral})";
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
//for(int i=0;i<rehdata->GetNbinsX();i++){
//cout<<rehdata->GetBinContent(i)<<endl;
//}
hFrame->GetXaxis()->SetTitleSize(0.05);
hFrame->GetXaxis()->SetLabelSize(0.05);
hFrame->GetYaxis()->SetTitleSize(0.05);
hFrame->GetYaxis()->SetLabelSize(0.05);
hFrame->GetXaxis()->SetLimits(xrange_pt[0],xrange_pt[1]);
if(JetID.Contains("p_{T}")&&svar!="Mult"){
hFrame->GetYaxis()->SetRangeUser(ymin2,ymax2);
hFrame->GetYaxis()->SetTitle(Form("%s",svar2.Data()));}
else{
hFrame->GetYaxis()->SetRangeUser(ymin1,ymax1);
hFrame->GetYaxis()->SetTitle(Form("%s",svar1.Data()));}
hFrame->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
hFrame->SetTitle("");
hFrame->DrawCopy();
TLatex *T1=new TLatex();
T1->SetNDC();
T1->SetTextAlign(12);
if(svar=="other"){
T1->SetTextSize(0.042);
T1->DrawLatex(0.18,0.90,JetID.Data());
}
else{
T1->SetTextSize(0.05);
T1->DrawLatex(0.20,0.88,JetID.Data());
}
T1->SetTextColor(1);
T1->SetTextFont(42);
rehdata->Draw("E1same");
rehPPb->Draw("HIST E1same");
TLegend *leg1=new TLegend(0.32,0.68,0.62,0.83);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetTextSize(0.05);
if(coll=="PPb"){
leg1->AddEntry(rehdata,"Pb going positive Data","lpf");
leg1->AddEntry(rehPPb,"PYTHIA+HIJING","lpf");
}
else if(coll=="PbP"){
leg1->AddEntry(rehdata,"Proton going positive Data","lpf");
leg1->AddEntry(rehPPb,"HIJING+PYTHIA","lpf");
}

if(svar=="other"){
if(k==2) leg1->Draw("same");}
else {
if(k==0) leg1->Draw("same");}
}
c1->Print(Form("pic/%soverlay_%s_Profile.png",coll.Data(),svar.Data()));
c1->Print(Form("pic/%soverlay_%s_Profile.pdf",coll.Data(),svar.Data()));
}
