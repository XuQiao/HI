#include "/home/xuq7/CMSSW_6_2_3_patch1/src/Centrality/NBD/parameter.h"
void Draw3panelNcoll(){
        gStyle->SetOptStat(kFALSE);
	TString cenvarname;
	TString outG="G1.root";
	TString dirname,name;
	TCanvas *c1 = new TCanvas("","",1000,400);
	c1->Divide(3);
	TString str="Ncoll";
	TGraphErrors *graph;// *Gri055_graph, *Gri101_graph;
        for(int cenvar=0;cenvar<3;cenvar++){
        cenvarname = cenvarlist[cenvar];
	TFile *f = TFile::Open(Form("%s/%s",cenvardir[cenvar].Data(),outG.Data()));
        c1->cd(cenvar+1)->SetFillColor(10);
        c1->cd(cenvar+1)->SetFrameFillColor(0);
        c1->cd(cenvar+1)->SetFrameBorderSize(0);
        c1->cd(cenvar+1)->SetFrameBorderMode(0);
        c1->cd(cenvar+1)->SetLeftMargin(0.15);
        c1->cd(cenvar+1)->SetBottomMargin(0.15);
        c1->cd(cenvar+1)->SetTopMargin(0.02);
        c1->cd(cenvar+1)->SetRightMargin(0.02);
        c1->cd(cenvar+1)->SetTicks(-1);

	int Gth=0;
	if(Gth==0)
        	name = "G0";
        else if(Gth<nGlau)
                name = Form("Glau_%d",Gth);
        else
                name = Form("bin_%d",Gth-nGlau+1);

	graph = (TGraphErrors*)f->Get(Form("std/%s_graph",str.Data()));
	Gri055_graph = (TGraphErrors*)f->Get(Form("Gri055/%s_graph",str.Data()));
	Gri101_graph = (TGraphErrors*)f->Get(Form("Gri101/%s_graph",str.Data()));
        TVectorD *centbin = (TVectorD*)f->Get(Form("std/%s/centbin",name.Data()));
        TVectorD *kpoint = (TVectorD*)f->Get(Form("std/%s/kpoint",name,Data()));
        TVectorD *method = (TVectorD*)f->Get(Form("std/%s/method",name.Data()));
	int N = centbin->GetUpb();
	TH1D* hist = new TH1D("","",N,0,N);
        hist->GetXaxis()->SetNdivisions(502);
        if(cenvar==0)	hist->SetYTitle(Form("<%s> and systematic errors",str.Data()));
        else	hist->SetYTitle("");
        hist->SetMinimum(-0.001);
        hist->SetMaximum(29.99);
        hist->GetXaxis()->CenterTitle(0);
        hist->GetYaxis()->CenterTitle(1);
        hist->GetYaxis()->SetTitleOffset(1.1);
        hist->GetXaxis()->SetTitleOffset(1.1);
        hist->GetXaxis()->SetTitleSize(0.056);
        hist->GetYaxis()->SetTitleSize(0.056);
        hist->GetXaxis()->SetLabelSize(0.05);
        hist->GetYaxis()->SetLabelSize(0.05);
        hist->GetXaxis()->SetLabelOffset(99);
	if((*method)[0]==0)
        hist->SetXTitle("Centrality");
        else
        hist->SetXTitle(cenvariable[cenvar]);
        hist->Draw();

	graph->SetTitle("g1");
	graph->SetMarkerStyle(20);
	graph->SetMarkerColor(1);
	graph->SetLineColor(1);
	graph->SetLineWidth(2);
	graph->SetMarkerSize(1.2);
	graph->Draw("Psameez");

	Gri055_graph->SetTitle("g2");
	Gri055_graph->SetMarkerStyle(33);
	Gri055_graph->SetMarkerColor(2);
	Gri055_graph->SetLineColor(2);
	Gri055_graph->SetLineWidth(2);
	Gri055_graph->SetMarkerSize(1.2);
	Gri055_graph->Draw("Psameez");

	Gri101_graph->SetTitle("g3");
	Gri101_graph->SetMarkerStyle(34);
	Gri101_graph->SetMarkerColor(4);
	Gri101_graph->SetLineColor(4);
	Gri101_graph->SetLineWidth(2);
	Gri101_graph->SetMarkerSize(1.2);
	Gri101_graph->Draw("Psameez");
	
	TString label[N];
	for(int i=0;i<N;i++)
	if((*method)[0]==0)label[i] = Form("%.2f-%.2f%%",(*centbin)[i]*100,(*centbin)[i+1]*100);
	else label[i] = Form("%.2f-%.2f",(*kpoint)[i],(*kpoint)[i+1]);

    	TLatex *tex1= new TLatex(0.6,26.05,"CMS Preliminary pPb \ #sqrt{s_{NN}} = 5.02 TeV");
    	TLatex *tex2= new TLatex(3,5,Form("%s",cenvariable[cenvar].Data()));
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
	leg0->AddEntry(graph,"standard","p");
	leg0->AddEntry(Gri055_graph,"Gribov #Omega=0.55","p");
	leg0->AddEntry(Gri101_graph,"Gribov #Omega=1.01","p");
	if(cenvar==2)	leg0->Draw();	
	}
	c1->SaveAs(Form("3panelNcoll.png"));	
	c1->SaveAs(Form("3panelNcoll.pdf"));	
}
