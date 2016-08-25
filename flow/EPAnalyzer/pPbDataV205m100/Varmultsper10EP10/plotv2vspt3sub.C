#include "par.h"
void plotv2vspt3sub(){
	TCanvas *c1 = new TCanvas();
	V2vsPt->SetLineColor(2);
	V2vsPt->SetTitle("v_{2} vs momentum");
	V2vsPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	V2vsPt->GetYaxis()->SetTitle("v_{2}");
	V2vsPt->GetYaxis()->SetRangeUser(0,0.15);
	V2vsPt->Draw();
	TGraphErrors* gr[neta];
        TLegend *tl = new TLegend(0.6,0.6,0.8,0.95);
	int color[neta]={1,2,3,4,6};
	int marker[neta]={20,24,29,34,25};
	for(int ieta=0;ieta<2;ieta++){
		gr[ieta] = plot(ieta,color[ieta],marker[ieta]);
		gr[ieta]->Draw("Psame");
		tl->AddEntry(gr[ieta],Form("EP |eta|>%.1f gap",etap[ieta]),"lp");
	}
        TLatex t;
        t.SetNDC();
        t.SetTextSize(0.04);
        t.DrawLatex(0.2,0.7,"3-sub event method");
	tl->SetFillColor(0);
	tl->SetBorderSize(0);
	tl->AddEntry(V2vsPt,"input v_{2}","lp");
	tl->Draw("same");
	c1->Print("v2vspt3sub.png");
}

TGraphErrors* plot(int ieta, int color,int marker){
	int ibin=0;
	TFile *f = TFile::Open("mergedVobs.root");
	TVectorD *vecDv2 = (TVectorD*)f->Get(Form("D_%d/E_%d/v2",ibin,ieta));
	TVectorD *vecDv2obs = (TVectorD*)f->Get(Form("D_%d/E_%d/v2obs",ibin,ieta));
	TVectorD *vecDv2sub3 = (TVectorD*)f->Get(Form("D_%d/E_%d/v23sub",ibin,ieta));
	TVectorD *vecDv2sub3err = (TVectorD*)f->Get(Form("D_%d/E_%d/v23suberr",ibin,ieta));
	TVectorD *vecDavgpt = (TVectorD*)f->Get(Form("D_%d/avgpt",ibin));
	double *avgpt = vecDavgpt->GetMatrixArray();
	double *v2 = vecDv2->GetMatrixArray();
	double *v2sub3 = vecDv2sub3->GetMatrixArray();
	double *v2sub3err = vecDv2sub3err->GetMatrixArray();
	//double *v2err = vecDv2err->GetMatrixArray();
	TGraphErrors *gr=new TGraphErrors(npt,avgpt,v2sub3,0,v2sub3err);
	gr->SetTitle("v_{2} vs momentum");
	gr->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	gr->GetYaxis()->SetTitle("v_{2}");
	gr->SetMarkerSize(1);
	gr->SetMarkerColor(color);
	gr->SetLineColor(color);
	gr->SetMarkerStyle(marker);
	f->Close();
	return gr;
}

