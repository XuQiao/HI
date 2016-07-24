#include "/home/xuq7/HI/jetRpA/RpA/Quality/root_setting.h"
const int ntotbin=5;
const int trkpointmin[ntotbin] = {120,150,185,220,260};
const int trkpointmax[ntotbin] = {150,185,220,260,300};
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

c1 = new TCanvas("c1"," ",1200,800);
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
    hFrame->GetXaxis()->SetRangeUser(-0.3,10.3);
    hFrame->SetMaximum(0.24);

for(int i=0;i<ntotbin;i++){
	TGraphErrors* gr24=new TGraphErrors(npt24,pt,v24[i],0,v24err[i]);
//    	TGraphErrors* grPFfull=plot("PFcandpt01to10",0,4,33,i);
//    	TGraphErrors* grPFon=plot("PFcandpt01to10/OnlyCharged",0,4,29,i);
//    	TGraphErrors* grtracklc01to10=plot("tracklcpt01to10",0,2,34,i);
    	TGraphErrors* grtracknormcpt03to6=plot("tracknormcpt03to6",0,2,21,i);
//  	TGraphErrors* grtracknormc03to6Sum=plot("tracknormcpt03to6",1,4,29,i);
    	TGraphErrors* grtracknormcpt03to3tracknormcpt03to6=plot("tracknormcpt03to3tracknormcpt03to6",0,2,29,i);
   // 	TGraphErrors* grtracknormc03to3=plot("tracknormcpt03to3",0,4,34,i);
    	TGraphErrors* grPFcandpt01to10tracknormcpt03to6=plot("PFcandpt01to10tracknormcpt03to6",0,4,34,i);
    	TGraphErrors* grPFcandpt03to6tracknormcpt03to6=plot("PFcandpt03to6tracknormcpt03to6",0,4,33,i);
//    	TGraphErrors* grPFcandpt03to3tracknormcpt03to6=plot("PFcandpt03to3tracknormcpt03to6",0,3,33,i);
        gr24->SetMarkerSize(1.6);
        gr24->SetMarkerColor(1);
        gr24->SetMarkerStyle(20);
        gr24->SetLineColor(1);
	c1->cd(i+1);
	hFrame->Draw();
	//grtracklc01to10->Draw("Psame");
	//grtracknormc03to6->Draw("Psame");
	//grPFfull->Draw("Psame");
	//grPFon->Draw("Psame");
   	grtracknormcpt03to6->Draw("Psame");
   	//grtracknormc03to6Sum->Draw("Psame");
   	grtracknormcpt03to3tracknormcpt03to6->Draw("Psame");
   	//grtracknormc03to3->Draw("Psame");
   	grPFcandpt01to10tracknormcpt03to6->Draw("Psame");
   	grPFcandpt03to6tracknormcpt03to6->Draw("Psame");
   	//grPFcandpt03to3tracknormcpt03to6->Draw("Psame");
	if(i!=ntotbin-1) gr24->Draw("Psame");
        if(i==0 || i==3)
                TLatex *tlx2 = new TLatex(0.3,0.8,Form("%d<Ntrkoffline<%d",trkpointmin[i],trkpointmax[i]));
        else
                TLatex *tlx2 = new TLatex(0.1,0.8,Form("%d<Ntrkoffline<%d",trkpointmin[i],trkpointmax[i]));
        tlx2->SetNDC();
        tlx2->Draw("same");
	}
        c1->cd(ntotbin+1);
        //TLatex *tlx1 = new TLatex(0.12,0.60,Form("%.1f<p_{T}<%.1f (GeV/c)",0.1,10.0));
        TLegend *tl = new TLegend(0.2,0.40,0.5,0.85);
        tl->SetFillColor(0);
        tl->SetBorderSize(0);
        tl->SetTextSize(0.065);
        tl->AddEntry(gr24,"v2 4-particle cum","lp");
      //  tl->AddEntry(grPFfull,"LYZ PF candidate full","lp");
      //  tl->AddEntry(grPFon,"LYZ PF candidate charged","lp");
      //  tl->AddEntry(grtracklc01to10,"LYZ track loose cut","lp");
        // tl->AddEntry(hFrame,Form("%.1f<p_{T}<%.1f (GeV/c)",0.1,10.0),"");
        tl->AddEntry(grtracknormcpt03to6,Form("track %.1f < p_{T} < %.1f",0.3,6.0),"lp");
        tl->AddEntry(grtracknormcpt03to3tracknormcpt03to6,Form("track %.1f < p_{T} < %.1f",0.3,3.0),"lp");
       //tl->AddEntry(grtracknormc03to6Sum,"LYZ general track Sum","lp");
       tl->AddEntry(grPFcandpt01to10tracknormcpt03to6,Form("PF cand %.1f < p_{T} < %.1f",0.1,10.0),"lp");
       tl->AddEntry(grPFcandpt03to6tracknormcpt03to6,Form("PF cand %.1f < p_{T} < %.1f",0.3,6.0),"lp");
       //tl->AddEntry(grPFcandpt03to3tracknormcpt03to6,Form("PF cand %.1f < p_{T} < %.1f",0.3,3.0),"lp");
        tl->AddEntry(hFrame,Form("%.1f<p_{T}<%.1f (GeV/c)",0.3,6.0),"");
        //tl->AddEntry(grtracknormc03to3tracknormcpt03to6,"LYZ general track 0.3 to 6","lp");
        //tl->AddEntry(grtracknormc03to3,"LYZ general track ","lp");
        //tl->AddEntry(hFrame,Form("%.1f<p_{T}<%.1f (GeV/c)",0.3,3.0),"");
        //tlx1->SetNDC();
        //tlx1->SetTextSize(0.065);
        tlx2->SetTextSize(0.065);
        hFrame->Draw();
        //tlx1->Draw("same");
        tl->Draw("same");
        c1->Print("v2vspt.png");
}

TGraphErrors* plot(TString dir, int isSum, int color, int marker, int i){

	if(isSum==0) TFile *f = TFile::Open(Form("%s/M%d%d/mergedv_Prod2.root",dir.Data(),trkpointmax[i],trkpointmin[i]));
	else TFile *f = TFile::Open(Form("%s/M%d%d/mergedv_Prod.root",dir.Data(),trkpointmax[i],trkpointmin[i]));
	int ibin=0;
	
	TVectorD *vecDv2 = (TVectorD*)f->Get(Form("D_%d/vmean",ibin));
	TVectorD *vecDv2err = (TVectorD*)f->Get(Form("D_%d/deltavmean",ibin));
	TVectorD *vecDavgpt = (TVectorD*)f->Get(Form("D_%d/avgpt",ibin));

	double *avgpt = vecDavgpt->GetMatrixArray();
	double *v2 = vecDv2->GetMatrixArray();
	double *v2err = vecDv2err->GetMatrixArray();
	int npt = vecDavgpt->GetNrows();
	
	TGraphErrors *gr=new TGraphErrors(npt,avgpt,v2,0,v2err);
	gr->SetMarkerStyle(marker);
	gr->SetMarkerSize(1.6);
	gr->SetMarkerColor(color);
	gr->SetLineColor(color);
	f->Close();
	return gr;
}

