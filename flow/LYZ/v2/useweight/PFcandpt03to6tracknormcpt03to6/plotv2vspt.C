#include "/home/xuq7/HI/jetRpA/RpA/Quality/root_setting.h"
void plotv2vspt(){

const int ntotbin=5;
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

const int trkpointmin[ntotbin] = {120,150,185,220,260};
const int trkpointmax[ntotbin] = {150,185,220,260,300};
int ibin=0;
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

for(int i=0;i<ntotbin;i++){
	TFile *fSum = TFile::Open(Form("M%d%d/mergedv_Prod2.root",trkpointmax[i],trkpointmin[i]));
	TFile *fProd = TFile::Open(Form("M%d%d/mergedv_Prod2.root",trkpointmax[i],trkpointmin[i]));
	TVectorD *vecDv2_Sum = (TVectorD*)fSum->Get(Form("D_%d/vmean",ibin));
	TVectorD *vecDv2err_Sum = (TVectorD*)fSum->Get(Form("D_%d/deltavmean",ibin));
	TVectorD *vecDavgpt_Sum = (TVectorD*)fSum->Get(Form("D_%d/avgpt",ibin));
	
	TVectorD *vecDv2_Prod = (TVectorD*)fProd->Get(Form("D_%d/vmean",ibin));
	TVectorD *vecDv2err_Prod = (TVectorD*)fProd->Get(Form("D_%d/deltavmean",ibin));
	TVectorD *vecDavgpt_Prod = (TVectorD*)fProd->Get(Form("D_%d/avgpt",ibin));

	double *avgpt_Sum = vecDavgpt_Sum->GetMatrixArray();
	double *v2_Sum = vecDv2_Sum->GetMatrixArray();
	double *v2err_Sum = vecDv2err_Sum->GetMatrixArray();

	double *avgpt_Prod = vecDavgpt_Prod->GetMatrixArray();
	double *v2_Prod = vecDv2_Prod->GetMatrixArray();
	double *v2err_Prod = vecDv2err_Prod->GetMatrixArray();
	int npt = vecDavgpt_Sum->GetNrows();
	
	c1->cd(i+1);
	if(i!=ntotbin-1)
	TGraphErrors *gr24=new TGraphErrors(npt24,pt,v24[i],0,v24err[i]);
	TGraphErrors *grSum=new TGraphErrors(npt,avgpt_Sum,v2_Sum,0,v2err_Sum);
	TGraphErrors *grProd=new TGraphErrors(npt,avgpt_Prod,v2_Prod,0,v2err_Prod);
	gr24->SetMarkerSize(1.3);
	gr24->SetMarkerColor(1);
	gr24->SetMarkerStyle(20);
	grSum->SetMarkerStyle(24);
	grProd->SetMarkerStyle(29);
	grSum->SetMarkerColor(2);
	grProd->SetMarkerColor(4);
	grSum->SetLineColor(2);
	grProd->SetLineColor(4);
	hFrame->Draw();
	gr24->Draw("Psame");
	//grSum->Draw("Psame");
	grProd->Draw("Psame");
	TLegend *tl = new TLegend(0.4,0.5,0.7,0.65);
	tl->SetFillColor(0);
	tl->SetBorderSize(0);
	tl->SetTextSize(0.05);
	tl->AddEntry(gr24,"v2 4-particle cum","lp");
	//tl->AddEntry(grSum,"LYZ Sum method","lp");
	tl->AddEntry(grProd,"LYZ Prod method","lp");
	if(i==0 || i==3) 
		TLatex *tlx2 = new TLatex(0.3,0.8,Form("%d<Ntrkoffline<%d",trkpointmin[i],trkpointmax[i]));
	else
		TLatex *tlx2 = new TLatex(0.1,0.8,Form("%d<Ntrkoffline<%d",trkpointmin[i],trkpointmax[i]));
	tlx2->SetNDC();
	tlx2->Draw("same");
	fSum->Close();
	fProd->Close();
	}
	c1->cd(ntotbin+1);
	TLatex *tlx0 = new TLatex(0.12,0.3,Form("track normal cut"));
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
	c1->Print("v2vspt.png");

}

