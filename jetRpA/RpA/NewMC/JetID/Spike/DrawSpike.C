#include "../../../Quality/root_setting.h"
#include "../../produceandcheck/file.h"
void DrawSpike(){
    TString coll = "Data";
    gStyle->SetOptStat(0);
    TCanvas *c1  = new TCanvas("c1","c1",650,600);
    c1->SetRightMargin(0.14);
    c1->SetLeftMargin(0.08);
    TCanvas *c2  = new TCanvas("c2","c2",650,600);
    c2->SetRightMargin(0.14);
    c2->SetLeftMargin(0.13);
    if(coll.Contains("Data")){
    TH2F* Spikeetaphi = (TH2F*)fDataPPbJetID->Get("Spikeetaphi");
    TH2F* Spikeetapt = (TH2F*)fDataPPbJetID->Get("Spikeetapt");
    }
    else{
	if(coll.Contains("real")){
	    TH2F* Spikeetaphi = (TH2F*)fMCPPb->Get("Spikeetaphi_real1");
	    TH2F* Spikeetapt = (TH2F*)fMCPPb->Get("Spikeetapt_real1");
	}
	else if(coll.Contains("fake")){
	    TH2F* Spikeetaphi = (TH2F*)fMCPPb->Get("Spikeetaphi_fake1");
	    TH2F* Spikeetapt = (TH2F*)fMCPPb->Get("Spikeetapt_fake1");
	}
	else{
	    TH2F* Spikeetaphi = (TH2F*)fMCPPb->Get("Spikeetaphi");
	    TH2F* Spikeetapt = (TH2F*)fMCPPb->Get("Spikeetapt");
    	}
    }
    c1->cd();
    c1->SetLogz();
    Spikeetaphi->Rebin2D(10,5);
    Spikeetaphi->SetTitle("");
    Spikeetaphi->GetXaxis()->SetTitle("jet_#eta");
    Spikeetaphi->GetYaxis()->SetTitle("jet_#phi");
    if(coll.Contains("MC")){
    double max = Spikeetaphi->GetMaximum();
    double min = 1e-11;
    }
    else{
    double max = 5e4;
    double min = 1;
    }
    Spikeetaphi->SetMaximum(max);
    Spikeetaphi->SetMinimum(min);
    Spikeetaphi->GetXaxis()->SetRangeUser(-3.5,3);
    Spikeetaphi->GetYaxis()->SetRangeUser(-4,4);
    Spikeetaphi->Draw("colz");
    c2->cd();
    c2->SetLogz();
    c2->SetLogy();
    Spikeetapt->SetTitle("");
    Spikeetapt->GetYaxis()->SetTitleOffset(1.5);
    Spikeetapt->Rebin2D(10,10);
    Spikeetapt->GetXaxis()->SetTitle("jet_#eta");
    Spikeetapt->GetYaxis()->SetTitle("jet_p_{T}");
    Spikeetapt->SetMaximum(max);
    Spikeetapt->SetMinimum(min);
    Spikeetapt->GetXaxis()->SetRangeUser(-3.5,3);
    Spikeetapt->GetYaxis()->SetRangeUser(28,599);
    Spikeetapt->GetYaxis()->SetMoreLogLabels();
    Spikeetapt->GetYaxis()->SetNoExponent();
    Spikeetapt->Draw("colz");
    bool save=1;
    if(save){
        c1->Print(Form("Spikeetaphi_%s.png",coll.Data()));
        c2->Print(Form("Spikeetapt_%s.png",coll.Data()));
    }
}

