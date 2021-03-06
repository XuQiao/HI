#include "../../Quality/root_setting.h"
#include "../produceandcheck/file.h"
void Drawselfmatch(){
	gStyle->SetOptStat(kFALSE);
	gStyle->SetOptFit(kFALSE);
	gStyle->SetErrorX(0);
	TH1D* hNmatched = (TH1D*)fMCPPb->Get("hNmatched");
	TH1D* hdeltaR = (TH1D*)fMCPPb->Get("hdeltaR");
	TH2F* hptdeltaptfrc = (TH2F*)fMCPPb->Get("hptdeltaptfrc");
	TH2F* hptdeltaptfrcmatch1 = (TH2F*)fMCPPb->Get("hptdeltaptfrcmatch1");
	TH2F* hptdeltaptfrcmatch20 = (TH2F*)fMCPPb->Get("hptdeltaptfrcmatch20");
	TH2F* hptdeltaptfrcmatch21 = (TH2F*)fMCPPb->Get("hptdeltaptfrcmatch21");
	TH2F* hptgenpt = (TH2F*)fMCPPb->Get("hptgenpt");
	TH2F* hptgenptmatch1 = (TH2F*)fMCPPb->Get("hptgenptmatch1");
	TH2F* hptgenptmatch20 = (TH2F*)fMCPPb->Get("hptgenptmatch20");
	TH2F* hptgenptmatch21 = (TH2F*)fMCPPb->Get("hptgenptmatch21");
	hdeltaR->Rebin(10);
	TCanvas *c1 = new TCanvas("c1","c1",600,600);
	TCanvas *c2 = new TCanvas("c2","c2",600,600);
	TCanvas *c3 = new TCanvas("c3","c3",600,600);
	TCanvas *c4 = new TCanvas("c4","c4",600,600);
	TCanvas *c5 = new TCanvas("c5","c5",1200,600);
	TCanvas *c6 = new TCanvas("c6","c6",600,600);
	TCanvas *c7 = new TCanvas("c7","c7",600,600);
	TCanvas *c8 = new TCanvas("c8","c8",1200,600);
	TLatex tl;
	tl.SetNDC();
	tl.SetTextSize(0.032);
	c1->cd()->SetLogy();
	hFrame = new TH1D("hFrame","hFrame",1000,0,1000);
	hFrame->SetTitle("");
	hFrame->GetXaxis()->SetTitle("matched gen");
	hFrame->GetYaxis()->SetTitle("Number of jets");
	hFrame->GetXaxis()->SetRangeUser(0,4);
	hFrame->GetYaxis()->SetRangeUser(1,1e9);
	hFrame->DrawCopy();
	fixedFontHist(hFrame,1.2,1.4);
	hNmatched->SetMarkerSize(1.2);
	hNmatched->SetMarkerStyle(20);
	hNmatched->Draw("Psame");
	c2->cd()->SetLogy();
	hFrame1 = new TH1D("hFrame1","hFrame1",100,0,10);
	fixedFontHist(hFrame1,1.2,1.4);
	hFrame1->SetTitle("");
	hFrame1->GetXaxis()->SetTitle("#Delta R");
	hFrame1->GetYaxis()->SetTitle("Number of jets");
	hFrame1->GetXaxis()->SetRangeUser(0,10);
	hFrame1->GetYaxis()->SetRangeUser(1,1e9);
	hFrame1->DrawCopy();
	hdeltaR->SetMarkerSize(1.2);
	hdeltaR->SetMarkerStyle(20);
	hdeltaR->Draw("Psame");
	TLine *l = new TLine(0.3,0,0.3,1e9);
	l->SetLineStyle(2);
	l->Draw("same");	
	c3->cd()->SetLogz();
	c3->cd()->SetLogy();
	c3->SetRightMargin(0.12);
	hFrame2 = new TH2F("hFrame2","hFrame2",1000,0,1000,100,0.,10.);
	fixedFontHist(hFrame2,1.2,1.3);
	hFrame2->SetTitle("");
	hFrame2->GetXaxis()->SetTitle("p^{reco}_{T}");
	hFrame2->GetYaxis()->SetTitle("|p^{reco}_{T}-p^{gen}_{T}|/p^{reco}_{T}");
	hFrame2->GetXaxis()->SetRangeUser(10,300);
	hFrame2->GetYaxis()->SetRangeUser(0,10);
	hFrame2->DrawCopy();
	tl.DrawLatex(0.4,0.6,"One jet matched case");

	hptdeltaptfrcmatch1->Draw("colz same");
	
	c4->cd()->SetLogz();
	c4->SetRightMargin(0.12);
	hFrame3 = new TH2F("hFrame3","hFrame3",1000,0,1000,1000,0.,1000.);
	fixedFontHist(hFrame3,1.2,1.3);
	hFrame3->SetTitle("");
	hFrame3->GetXaxis()->SetTitle("p^{reco}_{T}");
	hFrame3->GetYaxis()->SetTitle("p^{gen}_{T}");
	hFrame3->GetXaxis()->SetRangeUser(10,300);
	hFrame3->GetYaxis()->SetRangeUser(10,300);
	hFrame3->DrawCopy();
	hptgenptmatch1->Draw("colz same");
	tl.DrawLatex(0.5,0.3,"One jet matched case");

	c5->Divide(2,1);
	c5->cd(1)->SetLogz();
	c5->cd(1)->SetRightMargin(0.12);
	hFrame3->DrawCopy();
	tl.DrawLatex(0.45,0.2,"Two jet matched, first");
	hptgenptmatch20->Draw("colz same");
	c5->cd(2)->SetLogz();
	c5->cd(2)->SetRightMargin(0.12);
	hFrame3->DrawCopy();
	hptgenptmatch21->Draw("colz same");
	tl.DrawLatex(0.25,0.7,"Two jet matched, second");

	c6->cd()->SetLogz();
	c6->cd()->SetLogy();
	c6->cd()->SetRightMargin(0.12);
	hFrame2->DrawCopy();
	hptdeltaptfrc->Draw("colz same");
	tl.DrawLatex(0.5,0.3,"All jets");

	c7->cd()->SetLogz();
	c7->SetRightMargin(0.12);
	hFrame3->DrawCopy();
	hptgenpt->Draw("colz same");
	tl.DrawLatex(0.5,0.3,"All jets");

	c8->Divide(2,1);
	c8->cd(1)->SetLogz();
	c8->cd(1)->SetLogy();
	c8->cd(1)->SetRightMargin(0.12);
	hFrame2->DrawCopy();
	tl.DrawLatex(0.45,0.2,"Two jet matched, first");
	hptdeltaptfrcmatch20->Draw("colz same");
	c8->cd(2)->SetLogz();
	c8->cd(2)->SetLogy();
	c8->cd(2)->SetRightMargin(0.12);
	hFrame2->DrawCopy();
	hptdeltaptfrcmatch21->Draw("colz same");
	tl.DrawLatex(0.25,0.7,"Two jet matched, second");

	c1->Print("pic/hNmatched.png");
	c2->Print("pic/hdeltaR.png");
	c3->Print("pic/ptvsdeltaptfraction1.png");
	c4->Print("pic/ptvsgenpt1.png");
	c5->Print("pic/ptvsgenpt2.png");
	c6->Print("pic/ptvsdeltaptfraction.png");
	c7->Print("pic/ptvsgenpt.png");
	c8->Print("pic/ptvsdeltaptfraction2.png");
}

