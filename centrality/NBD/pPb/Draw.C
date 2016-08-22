#include "/afs/cern.ch/user/q/qixu/CMSSW_6_2_5/src/Centrality/NBD_Glauber_fit/NBD/parameter.h"
void Draw(){ 
 TCanvas *c1 = new TCanvas("c1","c1",1,1,550,460);
  c1->SetFillColor(10);
  c1->SetFrameFillColor(0);
  c1->SetFrameBorderSize(0);
  c1->SetFrameBorderMode(0);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.15);
  c1->SetTopMargin(0.02);
  c1->SetRightMargin(0.02);
  gStyle->SetOptStat(0);
  c1->SetTicks(-1);

	const int N=5;
	const int method=1;
 TString str="Ncoll";
 TString Glaustr="Gri101";
 TH1D* hist = new TH1D("","",N,0,N);
 hist->GetXaxis()->SetNdivisions(502);
if(method==0)
 hist->SetXTitle("Centrality");
else
 hist->SetXTitle("Multiplicity");
 hist->SetYTitle(Form("<%s> and systematic errors",str.Data()));
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
 hist->Draw();
	
	TFile *fds = TFile::Open(Form("double_side/G.root"));
	TFile *fplus = TFile::Open(Form("PlusEta4/G.root"));
	TFile *fntrk = TFile::Open(Form("Ntrk/G.root"));
	TGraphErrors* gr_ds = (TGraphErrors*)fds->Get(Form("%s/%s_graph",Glaustr.Data(),str.Data()));
	TGraphErrors* gr_plus = (TGraphErrors*)fplus->Get(Form("%s/%s_graph",Glaustr.Data(),str.Data()));
	TGraphErrors* gr_ntrk = (TGraphErrors*)fntrk->Get(Form("%s/%s_graph",Glaustr.Data(),str.Data()));
       // TVectorD *centbin = (TVectorD*)f->Get(Form("std/G0/centbin"));
       // TVectorD *kpoint = (TVectorD*)f->Get(Form("std/G0/kpoint"));

gr_ds->SetTitle("g1");
gr_ds->SetMarkerStyle(20);
gr_ds->SetMarkerColor(1);
gr_ds->SetLineColor(1);
gr_ds->SetLineWidth(2);
gr_ds->SetMarkerSize(1.2);
gr_ds->Draw("Psameez");

gr_plus->SetTitle("g2");
gr_plus->SetMarkerStyle(33);
gr_plus->SetMarkerColor(2);
gr_plus->SetLineColor(2);
gr_plus->SetLineWidth(2);
gr_plus->SetMarkerSize(1.2);
gr_plus->Draw("Psameez");

gr_ntrk->SetTitle("g3");
gr_ntrk->SetMarkerStyle(34);
gr_ntrk->SetMarkerColor(4);
gr_ntrk->SetLineColor(4);
gr_ntrk->SetLineWidth(2);
gr_ntrk->SetMarkerSize(1.2);
gr_ntrk->Draw("Psameez");
TString label[N];
for(int i=0;i<N;i++)
//	if(method==0)label[i] = Form("%.2f-%.2f%%",(*centbin)[i]*100,(*centbin)[i+1]*100);
//	else label[i] = Form("%.2f-%.2f",(*kpoint)[i],(*kpoint)[i+1]);
	label[i]=Form("%d",i);

    TLatex *tex1= new TLatex(0.6,26.05,"CMS Preliminary pPb \ #sqrt{s_{NN}} = 5.02 TeV");
    TLatex *tex2= new TLatex(4,9.05,Glaustr.Data());
    tex1->SetTextColor(1);
    tex2->SetTextColor(1);
    tex1->SetTextFont(42);
    tex2->SetTextFont(42);
    tex1->SetTextSize(0.05);
    tex2->SetTextSize(0.05);
    tex1->Draw();
    tex2->Draw();

double y = gPad->GetUymin() - 0.2*hist->GetYaxis()->GetBinWidth(1);
   TText t;
   t.SetTextAngle(0);
   t.SetTextSize(0.05);
   t.SetTextAlign(11);
   for (int i=0;i<N;i++) {
      double x = hist->GetXaxis()->GetBinCenter(i+1);
      t.DrawText(x,y,label[i]);
   }
TLegend *leg0 = new TLegend(0.28,0.63,0.50,0.85);
    leg0->SetFillColor(10);
    leg0->SetBorderSize(0.035);
    leg0->SetTextFont(42);
    leg0->SetTextSize(0.047);
    leg0->AddEntry(gr_ds,"HF both side","p");
    leg0->AddEntry(gr_plus,"HF pb going side","p");
    leg0->AddEntry(gr_ntrk,"Ntrkoffline","p");
	leg0->Draw();	
c1->SaveAs(Form("%s_%s.png",str.Data(),Glaustr.Data()));

}
