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
const double v26[5][npt24] = {
    {0.028648,0.0407009,0.0618194,0.0765241,0.0847754,0.107657,0.0916163,0.0584143,0.0498903},
    {0.0259366,0.0390773,0.0608349,0.0813255,0.0909296,0.0952114,0.0895709,0.0483789,0.0592095},
    {0.0237041,0.03821,0.0617625,0.080257,0.0951352,0.103015,0.104138,0.085058,0.0647899},
    {0.019040,0.039077,0.056074,0.075613,0.086659,0.113795,0.114627,0.086099,0.142567}
};
const double v26err[5][npt24] = {
    {0.00480585,0.00447563,0.00638442,0.00597334,0.00861569,0.0164925,0.0144436,0.0393584},
    {0.00339086,0.00295584,0.00468087,0.00858749,0.00599718,0.00754908,0.0075991,0.0169601,0.0216677},
    {0.00300052,0.00300538,0.00372043,0.00441082,0.00602353,0.00701487,0.00885424,0.00739591,0.0147496},
    {0.00339856,0.00219286,0.00523808,0.00539561,0.0063706,0.0102556,0.0148949,0.0266261,0.0296179}
};
const double v28[5][npt24] = {
    {0.0345912,0.0367697,0.0667635,0.065211,0.0741344,0.0898253,-0.0555855,0.151465,0.279243},
    {0.021038,0.0325744,0.0787748,0.0870162,0.0792368,0.111926,0.0580171,-0.00399227,-0.0692344},
    {0.0201575,0.0333076,0.0658784,0.0729223,0.0792012,0.0841368,0.0846554,0.0600356,6.81017e-05},
    {-0.008018,0.030227,0.021554,0.133730,0.138579,0.241780,0.171231,-0.370900,-0.776637}
};
const double v28err[5][npt24] = {
    {0.00946981,0.0117817,0.0133057,0.0129061,0.0254731,0.0297223,0.0296153,0.0769674,0.0810091},
    {0.00903331,0.00737717,0.00841947,0.0116376,0.0157374,0.016964,0.0269157,0.045247,0.0541035},
    {0.00863054,0.00988199,0.0122853,0.0161234,0.0234761,0.0275355,0.0310692,0.0574622,0.0818537},
    {0.0332846,0.0324435,0.0522983,0.059827,0.0521528,0.0913221,0.083228,0.228437,0.334497}
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
    hFrame->GetXaxis()->SetTitleSize(0.06);
    hFrame->GetYaxis()->SetTitleSize(0.06);	
    hFrame->GetXaxis()->SetLabelSize(0.06);
    hFrame->GetYaxis()->SetLabelSize(0.06);	
    hFrame->GetXaxis()->CenterTitle();
    hFrame->GetYaxis()->CenterTitle();	
    hFrame->GetXaxis()->SetRangeUser(-0.3,6.3);
    hFrame->SetMaximum(0.35);

for(int i=0;i<ntotbin;i++){
	TFile *fProd = TFile::Open(Form("M%d%d/mergedv_Prod2_sub.root",trkpointmax[i],trkpointmin[i]));
	//TFile *fProd = TFile::Open(Form("M%d%d/mergedv_Prod2.root",trkpointmax[i],trkpointmin[i]));
	//TVectorD *vecDv2_Proderr = (TVectorD*)fProderr->Get(Form("D_%d/vmeanmean",ibin));
	//TVectorD *vecDv2err_Proderr = (TVectorD*)fProderr->Get(Form("D_%d/sigmavmeanmean",ibin));
	//TVectorD *vecDavgpt_Proderr = (TVectorD*)fProderr->Get(Form("D_%d/avgavgpt",ibin));
	/*
	TVectorD *vecDv2_Prod = (TVectorD*)fProd->Get(Form("D_%d/vmean",ibin));
	TVectorD *vecDv2err_Prod = (TVectorD*)fProd->Get(Form("D_%d/deltavmean",ibin));
	TVectorD *vecDavgpt_Prod = (TVectorD*)fProd->Get(Form("D_%d/avgpt",ibin));
	double *avgpt_Prod = vecDavgpt_Prod->GetMatrixArray();
	double *v2_Prod = vecDv2_Prod->GetMatrixArray();
	double *v2err_Prod = vecDv2err_Prod->GetMatrixArray();
        */

	TVectorD *vecDv2_Prod = (TVectorD*)fProd->Get(Form("D_%d/vmeanmean",ibin));
	//TVectorD *vecDv2err_Prod = (TVectorD*)fProd->Get(Form("D_%d/deltavmeanmean",ibin));
	TVectorD *vecDv2err_Prod = (TVectorD*)fProd->Get(Form("D_%d/sigmavmeanmean",ibin));
	TVectorD *vecDavgpt_Prod = (TVectorD*)fProd->Get(Form("D_%d/avgavgpt",ibin));
	
        double *avgpt_Prod = vecDavgpt_Prod->GetMatrixArray();
	double *v2_Prod = vecDv2_Prod->GetMatrixArray();
	double *v2err_Prod = vecDv2err_Prod->GetMatrixArray();

	int npt = vecDavgpt_Prod->GetNrows();

	c1->cd(i+1);
	if(i!=ntotbin-1){
	TGraphErrors *gr24=new TGraphErrors(npt24,pt,v24[i],0,v24err[i]);
	TGraphErrors *gr26=new TGraphErrors(npt24,pt,v26[i],0,v26err[i]);
	TGraphErrors *gr28=new TGraphErrors(npt24,pt,v28[i],0,v28err[i]);
        }
	TGraphErrors *grProd=new TGraphErrors(npt,avgpt_Prod,v2_Prod,0,v2err_Prod);
	gr24->SetMarkerSize(1.4);
	gr24->SetMarkerColor(1);
	gr24->SetMarkerStyle(20);
	gr24->SetLineColor(1);
	gr26->SetMarkerSize(1.4);
	gr26->SetMarkerColor(4);
	gr26->SetMarkerStyle(34);
	gr26->SetLineColor(4);
	gr28->SetMarkerSize(1.4);
	gr28->SetMarkerColor(2);
	gr28->SetMarkerStyle(33);
	gr28->SetLineColor(2);
	grProd->SetMarkerStyle(20);
	grProd->SetMarkerColor(6);
	grProd->SetLineColor(6);
	hFrame->Draw();
	if(i!=ntotbin-1){
	gr24->Draw("Psame");
	gr26->Draw("Psame");
	gr28->Draw("Psame");
        }
	grProd->Draw("Psame");
	TLegend *tl = new TLegend(0.4,0.4,0.7,0.65);
	tl->SetFillColor(0);
	tl->SetBorderSize(0);
	tl->SetTextSize(0.05);
	tl->AddEntry(gr24,"v_{2}{4} cumulant","lp");
	tl->AddEntry(gr26,"v_{2}{6} cumulant","lp");
	tl->AddEntry(gr28,"v_{2}{8} cumulant","lp");
	tl->AddEntry(grProd,"v_{2}{LYZ}","lp");
	if(i==0 || i==3) 
		TLatex *tlx2 = new TLatex(0.3,0.8,Form("%d<Ntrkoffline<%d",trkpointmin[i],trkpointmax[i]));
	else
		TLatex *tlx2 = new TLatex(0.1,0.8,Form("%d<Ntrkoffline<%d",trkpointmin[i],trkpointmax[i]));
	tlx2->SetNDC();
	tlx2->Draw("same");
	fProd->Close();
	}
	c1->cd(ntotbin+1);
	TLatex *tlx0 = new TLatex(0.12,0.3,Form("track"));
        TLatex *tlx1 = new TLatex(0.12,0.25,Form("%.1f<#eta<%.1f",-2.4,2.4));
	tlx0->SetNDC();
	tlx1->SetNDC();
	tlx0->SetTextSize(0.045);
	tlx1->SetTextSize(0.045);
	tlx2->SetTextSize(0.045);
        hFrame->Draw();
	//tlx0->Draw("same");
	tlx1->Draw("same");
	tl->Draw("same");
	c1->Print("v2vspt.png");
	c1->Print("v2vspt.pdf");

}

