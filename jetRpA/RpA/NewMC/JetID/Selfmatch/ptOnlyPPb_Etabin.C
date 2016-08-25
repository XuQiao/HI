#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TMath.h"
#include "TLine.h"
#include "TObject.h"
#include "iostream"
#include "fstream"
#include "../../../Quality/root_setting.h"
#include "../../produceandcheck/file.h"

TH1* myRebin(TH1* histo,int Nbin, double binbound[]){
//rehisto->Scale(1/rehisto->Integral());

TH1* rehisto=histo->Clone("rehisto");
TH1* rehisto=(TH1*)histo->Rebin(Nbin,histo->GetName(),binbound);
TString histoname=Form("%s",histo->GetName());
if(!histoname.Contains("Prof")){
normalizeByBinWidth(rehisto);
}

return rehisto;
}

void ptOnlyPPb_Etabin(){

using namespace std;
gStyle->SetOptStat(kFALSE);
//gStyle->SetOptStat(100);
gStyle->SetErrorX(0);
gStyle->SetLabelFont(70);

//------------------------------------------------------------Get Histograms---------------------------------------------

c6 = new TCanvas("c6"," ",1200,600);
c7 = new TCanvas("c7"," ",1200,600);
makeMultiPanelCanvas(c6,4,2,0,0,0.25,0.2,0.03);
makeMultiPanelCanvas(c7,4,2,0,0,0.25,0.2,0.03);
TH1F* hFrame1=new TH1F("","",1000,0,1000);
TH1F* hFrame2=new TH1F("","",500,0,5);
fixedFontHist(hFrame1,2.0,3.0);
fixedFontHist(hFrame2,2.0,3.0);
hFrame1->SetTitle("");
hFrame2->SetTitle("");
TLegend *leg1=new TLegend(0.36,0.82,0.60,0.95);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetLineWidth(0);
leg1->SetTextSize(0.065);
TLegend *leg2=new TLegend(0.62,0.23,0.75,0.53);
leg2->SetBorderSize(0);
leg2->SetFillColor(0);
leg2->SetLineWidth(0);
leg2->SetTextSize(0.055);
TLatex T1;
T1.SetNDC();
T1.SetTextSize(0.065);
T1.SetTextFont(42);


for(int i=0;i<Neta;i++){
	if(canvas[i]==0){
		double ybase=0.15;	double xbase=0.28;
	}
	else if(canvas[i]>4){
		double ybase=0.30;	double xbase=0.05;
	}
	else{
		double ybase=0.15;	double xbase=0.05;
	}
    if(i==7){
	TH1F *hPPb_pt= (TH1F*)fMCPPbSm->Get(Form("jetpt"));
	TH1F *hPPb_pt_real= (TH1F*)fMCPPbSm->Get(Form("jetpt_real"));
	TH1F *hPPb_pt_sm_real= (TH1F*)fMCPPbSm->Get(Form("jetpt_selfmatch_real"));
	TH1F *hPPb_pt_sm1_real= (TH1F*)fMCPPbSm->Get(Form("jetpt_selfmatch1_real"));
	TH1F *hPPb_pt_fake= (TH1F*)fMCPPbSm->Get(Form("jetpt_fake"));
	TH1F *hPPb_pt_sm_fake= (TH1F*)fMCPPbSm->Get(Form("jetpt_selfmatch_fake"));
	TH1F *hPPb_pt_sm1_fake= (TH1F*)fMCPPbSm->Get(Form("jetpt_selfmatch1_fake"));
    }
    else{
	TH1F *hPPb_pt= (TH1F*)fMCPPbSm->Get(Form("jetptEtaBin%s",etabinnameswap[i].Data()));
	TH1F *hPPb_pt_real= (TH1F*)fMCPPbSm->Get(Form("jetptEtaBin%s_real",etabinnameswap[i].Data()));
//	TH1F *hPPb_pt_sm_real= (TH1F*)fMCPPbSm->Get(Form("jetptEtaBin_selfmatch_real%s",etabinnameswap[i].Data()));
//	TH1F *hPPb_pt_sm1_real= (TH1F*)fMCPPbSm->Get(Form("jetptEtaBin_selfmatch1_real%s",etabinnameswap[i].Data()));
	TH1F *hPPb_pt_fake= (TH1F*)fMCPPbSm->Get(Form("jetptEtaBin%s_fake",etabinnameswap[i].Data()));
	TH1F *hPPb_pt_sm_fake= (TH1F*)fMCPPbSm->Get(Form("jetptEtaBin%s_selfmatch_fake",etabinnameswap[i].Data()));
	TH1F *hPPb_pt_sm1_fake= (TH1F*)fMCPPbSm->Get(Form("jetptEtaBin%s_selfmatch1_fake",etabinnameswap[i].Data()));
    }


//---------------------------------------------------Rebin and Format histograms-----------------------------------------
TH1F* rehisto_hPPb_pt=(TH1F*)myRebin(hPPb_pt,Nbin_pt,binbound_pt);
TH1F* rehisto_hPPb_pt_real=(TH1F*)myRebin(hPPb_pt_real,Nbin_pt,binbound_pt);
TH1F* rehisto_hPPb_pt_fake=(TH1F*)myRebin(hPPb_pt_fake,Nbin_pt,binbound_pt);

TH1F* rehisto_hPPb_pt_sm_fake=(TH1F*)myRebin(hPPb_pt_sm_fake,Nbin_pt,binbound_pt);
TH1F* rehisto_hPPb_pt_sm1_fake=(TH1F*)myRebin(hPPb_pt_sm1_fake,Nbin_pt,binbound_pt);

rehisto_hPPb_pt->SetMarkerStyle(20);
rehisto_hPPb_pt->SetMarkerSize(1.5);
rehisto_hPPb_pt->SetMarkerColor(1);
rehisto_hPPb_pt->SetLineColor(1);
rehisto_hPPb_pt->SetLineColor(1);

//---------------------------------------PPb pT f/i&r/i ratio before and after cut-----------------------------------

c6->cd(canvas[i]+1)->SetGridx();
c6->cd(canvas[i]+1)->SetLogy();
    if(canvas[i]==0  || canvas[i]==4){
	hFrame1->GetYaxis()->SetTitle(Form("Ratio"));
        hFrame1->SetLabelSize(18,"Y");
    }
    else{
        hFrame1->SetLabelSize(0,"Y");
        hFrame1->GetYaxis()->SetTitle("");
    }
    if(canvas[i]>=4)
        hFrame1->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
hFrame1->GetXaxis()->SetLimits(28,599);
hFrame1->SetMaximum(5);
hFrame1->SetMinimum(5e-8);
hFrame1->DrawCopy();
TH1F* ratio_hPPb_pt_fakevsInc=(TH1F*)rehisto_hPPb_pt_fake->Clone("ratio_hPPb_pt_fakevsInc");
ratio_hPPb_pt_fakevsInc->Divide(rehisto_hPPb_pt);
ratio_hPPb_pt_fakevsInc->SetTitle("");
ratio_hPPb_pt_fakevsInc->SetMarkerStyle(20);
ratio_hPPb_pt_fakevsInc->SetMarkerSize(1.5);
ratio_hPPb_pt_fakevsInc->SetMarkerColor(2);
ratio_hPPb_pt_fakevsInc->SetLineColor(2);

TH1F* ratio_hPPb_pt_sm_fakevsInc=(TH1F*)rehisto_hPPb_pt_sm_fake->Clone("ratio_hPPb_pt_sm_fakevsInc");
ratio_hPPb_pt_sm_fakevsInc->Divide(rehisto_hPPb_pt);
ratio_hPPb_pt_sm_fakevsInc->SetTitle("");
ratio_hPPb_pt_sm_fakevsInc->SetMarkerStyle(28);
ratio_hPPb_pt_sm_fakevsInc->SetMarkerSize(1.5);
ratio_hPPb_pt_sm_fakevsInc->SetMarkerColor(4);
ratio_hPPb_pt_sm_fakevsInc->SetLineColor(4);

TH1F* ratio_hPPb_pt_sm1_fakevsInc=(TH1F*)rehisto_hPPb_pt_sm1_fake->Clone("ratio_hPPb_pt_sm1_fakevsInc");
ratio_hPPb_pt_sm1_fakevsInc->Divide(rehisto_hPPb_pt);
ratio_hPPb_pt_sm1_fakevsInc->SetTitle("");
ratio_hPPb_pt_sm1_fakevsInc->SetMarkerStyle(33);
ratio_hPPb_pt_sm1_fakevsInc->SetMarkerSize(1.5);
ratio_hPPb_pt_sm1_fakevsInc->SetMarkerColor(7);
ratio_hPPb_pt_sm1_fakevsInc->SetLineColor(7);

TH1F* ratio_hPPb_pt_realvsInc=(TH1F*)rehisto_hPPb_pt_real->Clone("ratio_hPPb_pt_realvsInc");
ratio_hPPb_pt_realvsInc->Divide(rehisto_hPPb_pt);
ratio_hPPb_pt_realvsInc->SetMarkerStyle(29);
ratio_hPPb_pt_realvsInc->SetMarkerSize(1.5);
ratio_hPPb_pt_realvsInc->SetMarkerColor(4);
ratio_hPPb_pt_realvsInc->SetLineColor(4);
if(canvas[i]!=4){
T1.DrawLatex(0.4,ybase-0.06,etastring[i]);
ratio_hPPb_pt_fakevsInc->Draw("same");
ratio_hPPb_pt_realvsInc->Draw("same");
ratio_hPPb_pt_sm_fakevsInc->Draw("same");
ratio_hPPb_pt_sm1_fakevsInc->Draw("same");
TLine *l = new TLine(hFrame1->GetXaxis()->GetXmin(),1,hFrame1->GetXaxis()->GetXmax(),1);
l->SetLineStyle(2);
l->Draw("same");
}
if(canvas[i]==4){
T1.DrawLatex(0.28,0.73,Form("No JetID Cut"));
leg1->Draw("same");
TLegend* leg3 = new TLegend(*leg2);
leg3->Clear();
leg3->SetX1(0.45);
leg3->SetX2(0.9);
leg3->AddEntry(ratio_hPPb_pt_fakevsInc,"fake Jet/Inc Jet","lp");
leg3->AddEntry(ratio_hPPb_pt_realvsInc,"real Jet/Inc Jet","lp");
leg3->AddEntry(ratio_hPPb_pt_sm_fakevsInc,"fake match/Inc","lp");
leg3->AddEntry(ratio_hPPb_pt_sm1_fakevsInc,"fake match+#Delta p_{T}/Inc","lp");
leg3->Draw("same");
}

//-------------------------------PPb pT distribution before and after cut-------------------------------------------
c7->cd(canvas[i]+1)->SetGridx();
c7->cd(canvas[i]+1)->SetLogy();
    if(canvas[i]==0  || canvas[i]==4){
	hFrame1->GetYaxis()->SetTitle(Form("Jet Yield"));
        hFrame1->SetLabelSize(18,"Y");
    }
    else{
        hFrame1->SetLabelSize(0,"Y");
	hFrame1->GetYaxis()->SetTitle("");
    }
    if(canvas[i]>=4)
        hFrame1->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
	
   hFrame1->GetXaxis()->SetNdivisions(510);
   hFrame1->GetXaxis()->SetLimits(28,599);
   hFrame1->SetMinimum(2e-15);
   hFrame1->SetMaximum(5e-3);
   hFrame1->DrawCopy();

rehisto_hPPb_pt->SetMarkerStyle(24);
rehisto_hPPb_pt->SetMarkerSize(1.5);
rehisto_hPPb_pt_real->SetMarkerStyle(29);
rehisto_hPPb_pt_real->SetMarkerColor(4);
rehisto_hPPb_pt_real->SetLineColor(4);
rehisto_hPPb_pt_real->SetMarkerSize(1.5);
rehisto_hPPb_pt_fake->SetMarkerStyle(20);
rehisto_hPPb_pt_fake->SetMarkerColor(2);
rehisto_hPPb_pt_fake->SetLineColor(2);
rehisto_hPPb_pt_fake->SetMarkerSize(1.5);

rehisto_hPPb_pt_sm_fake->SetMarkerStyle(28);
rehisto_hPPb_pt_sm_fake->SetMarkerColor(4);
rehisto_hPPb_pt_sm_fake->SetLineColor(4);
rehisto_hPPb_pt_sm_fake->SetMarkerSize(1.5);

rehisto_hPPb_pt_sm1_fake->SetMarkerStyle(33);
rehisto_hPPb_pt_sm1_fake->SetMarkerColor(7);
rehisto_hPPb_pt_sm1_fake->SetLineColor(7);
rehisto_hPPb_pt_sm1_fake->SetMarkerSize(1.5);

if(canvas[i]!=4){
T1.DrawLatex(0.4,ybase-0.06,etastring[i]);
rehisto_hPPb_pt->Draw("E1same");
rehisto_hPPb_pt_real->Draw("E1same");
rehisto_hPPb_pt_fake->Draw("E1same");
rehisto_hPPb_pt_sm_fake->Draw("E1same");
rehisto_hPPb_pt_sm1_fake->Draw("E1same");
}
if(canvas[i]==4){
T1.DrawLatex(0.28,0.73,Form("No JetID Cut"));
leg1->Draw("same");
leg2->AddEntry(rehisto_hPPb_pt,"Inclusive","lp");
leg2->AddEntry(rehisto_hPPb_pt_fake,"fake","lp");
leg2->AddEntry(rehisto_hPPb_pt_real,"real","lp");
leg2->AddEntry(rehisto_hPPb_pt_sm_fake,"fake match","lp");
leg2->AddEntry(rehisto_hPPb_pt_sm1_fake,"fake match+#Delta p_{T}","lp");
leg2->Draw("same");
}
}
c6->Print(Form("ratio_jetpt_befcut_Etabin.png"));
c7->Print(Form("jetpt_befcut_Etabin.png"));
}
