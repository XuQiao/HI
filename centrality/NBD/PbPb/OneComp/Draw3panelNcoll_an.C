#include "/home/xuq7/HI/centrality/NBD/parameter.h"
void Draw3panelNcoll_an(){
        gStyle->SetOptStat(kFALSE);
	TString outG="G1.root";
	TString dirname,name;
	TCanvas *c1 = new TCanvas("","",1000,400);
	c1->Divide(3);
	TString str="Ncoll";
	TGraphErrors *ds_graph;
	TGraphErrors *ss_graph;
	TGraphErrors *Ntrk_graph;
	fds = TFile::Open(Form("%s/%s",cenvardir[0].Data(),outG.Data()));
	fss = TFile::Open(Form("%s/%s",cenvardir[1].Data(),outG.Data()));
	fNtrk = TFile::Open(Form("%s/%s",cenvardir[2].Data(),outG.Data()));
        for(int sth=0;sth<3;sth++){
	if(sth==0)	dirname="std";
	else if(sth==1)	dirname="Gri055";
	else	dirname="Gri101";
        c1->cd(sth+1)->SetFillColor(10);
        c1->cd(sth+1)->SetFrameFillColor(0);
        c1->cd(sth+1)->SetFrameBorderSize(0);
        c1->cd(sth+1)->SetFrameBorderMode(0);
        c1->cd(sth+1)->SetLeftMargin(0.15);
        c1->cd(sth+1)->SetBottomMargin(0.15);
        c1->cd(sth+1)->SetTopMargin(0.02);
        c1->cd(sth+1)->SetRightMargin(0.01);
        c1->cd(sth+1)->SetTicks(-1);

	int Gth=0;
	if(Gth==0)
        	name = "G0";
        else if(Gth<nGlau)
                name = Form("Glau_%d",Gth);
        else
                name = Form("bin_%d",Gth-nGlau+1);

	ds_graph = (TGraphErrors*)fds->Get(Form("%s/%s_graph",dirname.Data(),str.Data()));
	ss_graph = (TGraphErrors*)fss->Get(Form("%s/%s_graph",dirname.Data(),str.Data()));
	Ntrk_graph = (TGraphErrors*)fNtrk->Get(Form("%s/%s_graph",dirname.Data(),str.Data()));
        TVectorD *centbin = (TVectorD*)fds->Get(Form("%s/%s/centbin",dirname.Data(),name.Data()));
        TVectorD *kpoint = (TVectorD*)fds->Get(Form("%s/%s/kpoint",dirname.Data(),name.Data()));
        TVectorD *method = (TVectorD*)fds->Get(Form("%s/%s/method",dirname.Data(),name.Data()));
	int N = centbin->GetUpb();
	TH1D* hist = new TH1D("","",N,0,N);
        hist->GetXaxis()->SetNdivisions(502);
        if(sth==0)	hist->SetYTitle(Form("<%s> and systematic errors",str.Data()));
        else	hist->SetYTitle("");
        hist->SetMinimum(-0.001);
        hist->SetMaximum(29.99);
        hist->GetXaxis()->CenterTitle(0);
        hist->GetYaxis()->CenterTitle(1);
        hist->GetYaxis()->SetTitleOffset(1.1);
        hist->GetXaxis()->SetTitleOffset(1.4);
        hist->GetXaxis()->SetTitleSize(0.056);
        hist->GetYaxis()->SetTitleSize(0.056);
        hist->GetXaxis()->SetLabelSize(0.05);
        hist->GetYaxis()->SetLabelSize(0.05);
        hist->GetXaxis()->SetLabelOffset(99);
	if((*method)[0]==0)
        hist->SetXTitle("Centrality");
        else{
	cout<<"Cannot make graphs with different x axis"<<endl,exit();
	}
        hist->Draw();

	ds_graph->SetTitle("g1");
	ds_graph->SetMarkerStyle(20);
	ds_graph->SetMarkerColor(1);
	ds_graph->SetLineColor(1);
	ds_graph->SetLineWidth(2);
	ds_graph->SetMarkerSize(1.2);
	ds_graph->Draw("Psameez");

	ss_graph->SetTitle("g2");
	ss_graph->SetMarkerStyle(33);
	ss_graph->SetMarkerColor(2);
	ss_graph->SetLineColor(2);
	ss_graph->SetLineWidth(2);
	ss_graph->SetMarkerSize(1.2);
	ss_graph->Draw("Psameez");

	Ntrk_graph->SetTitle("g3");
	Ntrk_graph->SetMarkerStyle(34);
	Ntrk_graph->SetMarkerColor(4);
	Ntrk_graph->SetLineColor(4);
	Ntrk_graph->SetLineWidth(2);
	Ntrk_graph->SetMarkerSize(1.2);
	Ntrk_graph->Draw("Psameez");
	
	TString label[N];
	for(int i=0;i<N;i++)
	if((*method)[0]==0)label[i] = Form("%.2f-%.2f%%",(*centbin)[i]*100,(*centbin)[i+1]*100);
	else label[i] = Form("%.2f-%.2f",(*kpoint)[i],(*kpoint)[i+1]);

    	TLatex *tex1= new TLatex(0.1,26.05,"CMS Preliminary pPb \ #sqrt{s_{NN}} = 5.02 TeV");
    	TLatex *tex2= new TLatex(2,5,Form("%s",Grilist[sth].Data()));
    	tex1->SetTextColor(1);
    	tex2->SetTextColor(1);
    	tex1->SetTextFont(42);
    	tex2->SetTextFont(42);
   	tex1->SetTextSize(0.05);
   	tex2->SetTextSize(0.045);
    	tex1->Draw();
    	tex2->Draw();

	double y = gPad->GetUymin();
// - 0.2*h->GetYaxis()->GetBinWidth(1);
	TText t;
	t.SetTextAngle(45);
	t.SetTextSize(0.03);
	t.SetTextAlign(33);
	for (int i=0;i<N;i++) {
		double x = hist->GetXaxis()->GetBinCenter(i+1);
		t.DrawText(x,y,label[i]);
	}
	TLegend *leg0 = new TLegend(0.28,0.63,0.50,0.85);
	leg0->SetFillColor(10);
	leg0->SetBorderSize(0.035);
	leg0->SetTextFont(42);
	leg0->SetTextSize(0.047);
	leg0->AddEntry(ds_graph,cenvariable[0],"p");
	leg0->AddEntry(ss_graph,cenvariable[1],"p");
	leg0->AddEntry(Ntrk_graph,cenvariable[2],"p");
	if(sth==2)	leg0->Draw();	
	}
	c1->SaveAs(Form("3panelNcoll_an.png"));	
	c1->SaveAs(Form("3panelNcoll_an.pdf"));	
}
