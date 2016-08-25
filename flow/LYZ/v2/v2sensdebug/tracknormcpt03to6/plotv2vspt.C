#include "/home/xuq7/HI/jetRpA/RpA/Quality/root_setting.h"
#include "M150120/par.h"
const int ntotbin=5;
const int trkpointmin[ntotbin] = {120,150,185,220,260};
const int trkpointmax[ntotbin] = {150,185,220,260,300};
const int ndir = 6;
const int marker[] = {20,24,27,30,31,29};
const int color[] = {2,3,4,6,7,8};
int ibin=0;
void plotv2vspt(){

const int npt24=9;
const double pt[npt24]={0.4,0.7,1.2,1.7,2.2,2.7,3.4,4.4,5.4};
const double v24[5][npt24]={
{0.02809,0.04599,0.06869,0.08440,0.09715,0.09357,0.07924,0.04701,0.01620},
{0.02682,0.04447,0.06771,0.08586,0.09932,0.09782,0.08714,0.06139,0.05068},
{0.02582,0.04324,0.06644,0.08396,0.09725,0.1032,0.09603,0.07668,0.0605},
{0.02507,0.04214,0.06471,0.08294,0.09241,0.09917,0.09663,0.08479,0.06271},
};
const double v24err[5][npt24]={
{0.0007,0.00057,0.00085,0.0014,0.0022,0.0034,0.0041,0.0072,0.011},
{0.00058,0.00043,0.00063,0.00099,0.0016,0.0025,0.0031,0.0054,0.0086},
{0.00053,0.00038,0.00056,0.00088,0.0014,0.0022,0.0027,0.0049,0.0076},
{0.00079,0.00057,0.00084,0.0013,0.002,0.0032,0.004,0.0075,0.012},
};

c1 = new TCanvas("c1"," ",1200,700);
makeMultiPanelCanvas(c1,3,2,0,0,0.25,0.2,0.03);
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetErrorX(0);
    TH1D *hFrame = new TH1D("","",80,-1,7);
    hFrame->SetTitle("v_{2} vs momentum");
    hFrame->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hFrame->GetYaxis()->SetTitle("v_{2}");	
    hFrame->GetXaxis()->SetTitleSize(0.04);
    hFrame->GetYaxis()->SetTitleSize(0.04);
    hFrame->GetXaxis()->SetRangeUser(-0.3,6.3);
    hFrame->SetMaximum(0.2);
    TGraphErrors* grProd_loop[ndir];
    TLegend *tl = new TLegend(0.2,0.4,0.6,0.80);
for(int i=0;i<ntotbin;i++){
	c1->cd(i+1);
	if(i!=ntotbin-1)
	TGraphErrors *gr24=new TGraphErrors(npt24,pt,v24[i],0,v24err[i]);
        for(int iloop=0;iloop<ndir;iloop++)
	grProd_loop[iloop]=plotTG(i,iloop,marker[iloop],color[iloop]);
//	TGraphErrors *grProd=plotTG(i,0,-1,24,4);
	gr24->SetMarkerSize(1.3);
	gr24->SetMarkerColor(1);
	gr24->SetMarkerStyle(20);
	hFrame->Draw();
	gr24->Draw("Psame");
        for(int iloop=0;iloop<ndir;iloop++)
	grProd_loop[iloop]->Draw("Psame");
//	grProd->Draw("Psame");
	tl->SetFillColor(0);
	tl->SetBorderSize(0);
	tl->SetTextSize(0.05);
            if(i==1)
	tl->AddEntry(gr24,"v2 4-particle cumu","lp");
        for(int iloop=0;iloop<ndir;iloop++){
            TFile *fV2 = TFile::Open(Form("M%d%d/parin_%d.root",trkpointmax[i],trkpointmin[i],iloop));
            TVectorD *Vmean = (TVectorD*)fV2->Get(Form("D_%d/Vmean",ibin));
            TVectorD *Vcalcmean = (TVectorD*)fV2->Get(Form("D_%d/Vcalcmean",ibin));
            TVectorD *r0 = (TVectorD*)fV2->Get(Form("D_%d/D_0/r0",ibin));
            double r0mean = r0->Sum()/r0->GetNrows();
            if(i==1)
            tl->AddEntry(grProd_loop[iloop],Form("LYZ V^{mean}_{2} = %.3f, r^{mean}_{0}=%.3f",(*Vcalcmean)[0],r0mean),"lp");
        }
	if(i==0 || i==3) 
		TLatex *tlx2 = new TLatex(0.3,0.8,Form("%d<Ntrkoffline<%d",trkpointmin[i],trkpointmax[i]));
	else
		TLatex *tlx2 = new TLatex(0.1,0.8,Form("%d<Ntrkoffline<%d",trkpointmin[i],trkpointmax[i]));
	tlx2->SetNDC();
	tlx2->Draw("same");
	//tl->Draw("same");
	}
	c1->cd(ntotbin+1);
	TLatex *tlx0 = new TLatex(0.12,0.3,Form("%s","track"));
        TLatex *tlx1 = new TLatex(0.12,0.25,Form("%.1f<p_{T}<%.1f (GeV/c)",0.3,6.0));
	tlx0->SetNDC();
	tlx1->SetNDC();
	tlx0->SetTextSize(0.045);
	tlx1->SetTextSize(0.045);
	tlx2->SetTextSize(0.045);
	hFrame->Draw();
	tlx0->Draw("same");
	tlx1->Draw("same");
	tl->Draw("same");
	c1->Print(Form("v2vspt.png"));

}

TGraphErrors* plotTG(int i, int iloop, int marker, int color){
        TFile *f = TFile::Open(Form("M%d%d/mergedv_Prod2_sub_loop%d.root",trkpointmax[i],trkpointmin[i],iloop));
        TVectorD *vecDv2 = (TVectorD*)f->Get(Form("D_%d/vmeanmean",ibin));
        TVectorD *vecDv2err = (TVectorD*)f->Get(Form("D_%d/sigmavmeanmean",ibin));
        TVectorD *vecDavgpt = (TVectorD*)f->Get(Form("D_%d/avgavgpt",ibin));

        double *avgpt = vecDavgpt->GetMatrixArray();
        double *v2 = vecDv2->GetMatrixArray();
        double *v2err = vecDv2err->GetMatrixArray();
        int npt = vecDavgpt->GetNrows();
	TGraphErrors *gr=new TGraphErrors(npt,avgpt,v2,0,v2err);
	gr->SetMarkerSize(1.3);
        gr->SetMarkerStyle(marker);
        gr->SetMarkerColor(color);
        gr->SetLineColor(color);
	f->Close();
	return gr;
}


