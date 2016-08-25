#include "../../../Quality/root_setting.h"
#include "../../produceandcheck/file.h"
TH1* myRebin(TH1* histo,int Nbin, double binbound[]){
TH1* rehisto=histo->Clone("rehisto");
TH1* rehisto=(TH1*)histo->Rebin(Nbin,histo->GetName(),binbound);
normalizeByBinWidth(rehisto);
return rehisto;
}

void Compareids(){
    gStyle->SetOptStat(kFALSE);
  //  gStyle->SetOptStat(100);
    gStyle->SetErrorX(0);
    gStyle->SetLabelFont(70);
    TString coll = "MCPPb";
    c1 = new TCanvas("c1"," ",1200,600);
    makeMultiPanelCanvas(c1,4,2,0,0,0.25,0.2,0.03);
    TH1F* hFrame=new TH1F("","",1000,0,1000);
    fixedFontHist(hFrame,2.0,3.0);
    TLegend *leg1=new TLegend(0.62,0.23,0.75,0.53);
    leg1->SetTextSize(0.055);
    leg1->SetBorderSize(0);
    leg1->SetFillColor(0);
    leg1->SetLineWidth(0);
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
	TH1F *hPPb_pt= (TH1F*)fMCPPb->Get(Form("jetpt"));
	TH1F *hPPb_pt_real= (TH1F*)fMCPPb->Get(Form("jetpt_real1"));
	TH1F *hPPb_pt_fake= (TH1F*)fMCPPb->Get(Form("jetpt_fake1"));
	TH1F *hPPb_pt_ids= (TH1F*)fMCPPbids->Get(Form("jetpt"));
	TH1F *hPPb_pt_ids_real= (TH1F*)fMCPPbids->Get(Form("jetpt_real1"));
	TH1F *hPPb_pt_ids_fake= (TH1F*)fMCPPbids->Get(Form("jetpt_fake1"));
    }
    else{
	TH1F *hPPb_pt= (TH1F*)fMCPPb->Get(Form("jetptEtaBin%s",etabinnameswap[i].Data()));
	TH1F *hPPb_pt_real= (TH1F*)fMCPPb->Get(Form("jetptEtaBin%s_real1",etabinnameswap[i].Data()));
	TH1F *hPPb_pt_fake= (TH1F*)fMCPPb->Get(Form("jetptEtaBin%s_fake1",etabinnameswap[i].Data()));
	TH1F *hPPb_pt_ids= (TH1F*)fMCPPbids->Get(Form("jetptEtaBin%s",etabinnameswap[i].Data()));
	TH1F *hPPb_pt_ids_real= (TH1F*)fMCPPbids->Get(Form("jetptEtaBin%s_real1",etabinnameswap[i].Data()));
	TH1F *hPPb_pt_ids_fake= (TH1F*)fMCPPbids->Get(Form("jetptEtaBin%s_fake1",etabinnameswap[i].Data()));
    }
TH1F* rehisto_hPPb_pt=(TH1F*)myRebin(hPPb_pt,Nbin_pt,binbound_pt);
TH1F* rehisto_hPPb_pt_real=(TH1F*)myRebin(hPPb_pt_real,Nbin_pt,binbound_pt);
TH1F* rehisto_hPPb_pt_fake=(TH1F*)myRebin(hPPb_pt_fake,Nbin_pt,binbound_pt);
TH1F* rehisto_hPPb_pt_ids=(TH1F*)myRebin(hPPb_pt_ids,Nbin_pt,binbound_pt);
TH1F* rehisto_hPPb_pt_ids_real=(TH1F*)myRebin(hPPb_pt_ids_real,Nbin_pt,binbound_pt);
TH1F* rehisto_hPPb_pt_ids_fake=(TH1F*)myRebin(hPPb_pt_ids_fake,Nbin_pt,binbound_pt);

c1->cd(canvas[i]+1);
c1->cd(canvas[i]+1)->SetGridx();
    if(canvas[i]==0  || canvas[i]==4){
        hFrame->GetYaxis()->SetTitle("Selection efficiency");
        hFrame->SetLabelSize(18,"Y");
    }
    else{
        hFrame->SetLabelSize(0,"Y");
        hFrame->GetYaxis()->SetTitle("");
    }
    if(canvas[i]>=4)
        hFrame->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
    hFrame->GetXaxis()->SetRangeUser(28,599);
    hFrame->GetXaxis()->SetNdivisions(510);
    hFrame->SetMinimum(0.98);
    hFrame->SetMaximum(1.01);
    hFrame->DrawCopy();
    if(canvas[i]!=4){
        TH1D* ratio_hPPb_pt_ids = (TH1D*)rehisto_hPPb_pt_ids->Clone("ratio_hPPb_pt_ids");
        TH1D* ratio_hPPb_pt_ids_real = (TH1D*)rehisto_hPPb_pt_ids_real->Clone("ratio_hPPb_pt_ids_real");
        TH1D* ratio_hPPb_pt_ids_fake = (TH1D*)rehisto_hPPb_pt_ids_fake->Clone("ratio_hPPb_pt_ids_fake");
        ratio_hPPb_pt_ids->Divide(rehisto_hPPb_pt);
        ratio_hPPb_pt_ids_real->Divide(rehisto_hPPb_pt_real);
        ratio_hPPb_pt_ids_fake->Divide(rehisto_hPPb_pt_fake);
        ratio_hPPb_pt_ids->SetMarkerSize(1.4);
        ratio_hPPb_pt_ids->SetMarkerStyle(24);
        ratio_hPPb_pt_ids->SetMarkerColor(1);
        ratio_hPPb_pt_ids->SetLineColor(1);
        ratio_hPPb_pt_ids_real->SetMarkerSize(1.4);
        ratio_hPPb_pt_ids_real->SetMarkerStyle(29);
        ratio_hPPb_pt_ids_real->SetMarkerColor(4);
        ratio_hPPb_pt_ids_real->SetLineColor(4);
        ratio_hPPb_pt_ids_fake->SetMarkerSize(1.4);
        ratio_hPPb_pt_ids_fake->SetMarkerStyle(20);
        ratio_hPPb_pt_ids_fake->SetMarkerColor(2);
        ratio_hPPb_pt_ids_fake->SetLineColor(2);
        getRidYError(ratio_hPPb_pt_ids);
        getRidYError(ratio_hPPb_pt_ids_real);
        getRidYError(ratio_hPPb_pt_ids_fake);
        ratio_hPPb_pt_ids_real->Draw("Psame");
        //ratio_hPPb_pt_ids_fake->Draw("Psame");
        ratio_hPPb_pt_ids->Draw("Psame");
        T1.SetTextSize(0.065);
        T1.DrawLatex(0.4,ybase,etastring[i]);
        TLine *l =new TLine(28,1,600,1);
        l->SetLineStyle(2);
        l->Draw("same");
    }
    if(canvas[i]==4){
        leg1->AddEntry(ratio_hPPb_pt_ids,"Inclusive","lp");
        leg1->AddEntry(ratio_hPPb_pt_ids_fake,"Fake","lp");
        leg1->AddEntry(ratio_hPPb_pt_ids_real,"Real","lp");
        leg1->Draw("same");
    }
}
}

