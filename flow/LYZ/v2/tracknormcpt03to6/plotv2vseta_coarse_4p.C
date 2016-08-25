#include "/home/xuq7/HI/jetRpA/RpA/Quality/root_setting.h"
void plotv2vseta_coarse_4p(){

const int ntotbin=4;
float sys = 0.115;
const int trkpointmin[ntotbin] = {120,150,185,220};
const int trkpointmax[ntotbin] = {150,185,220,260};
c1 = new TCanvas("c1"," ",1200,350);
makeMultiPanelCanvas(c1,4,1,0,0,0.15,0.15,0.01);
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetErrorX(0);
    TH1D *hFrame = new TH1D("","",80,-3,3);
    hFrame->SetTitle("");
    hFrame->SetTitle("v_{2} vs #eta");
    hFrame->GetXaxis()->SetTitle("#eta");
    hFrame->GetYaxis()->SetTitle("v_{2}");	
    hFrame->GetXaxis()->SetTitleSize(0.06);
    hFrame->GetYaxis()->SetTitleSize(0.06);
    hFrame->GetXaxis()->SetLabelSize(0.06);
    hFrame->GetYaxis()->SetLabelSize(0.06);
    hFrame->GetXaxis()->CenterTitle();
    hFrame->GetYaxis()->CenterTitle();
    hFrame->GetYaxis()->SetNdivisions(505);
    hFrame->GetXaxis()->SetRangeUser(-2.5,2.5);
    hFrame->SetMaximum(0.18);

int xbin=0;
for(int i=0;i<ntotbin;i++){
	//TFile *fProderr = TFile::Open(Form("M%d%d/mergedv_Prod2_eta_sub.root",trkpointmax[i],trkpointmin[i]));
	TFile *fProderr = TFile::Open(Form("M%d%d/mergedv_Prod2_coarse_eta_sub.root",trkpointmax[i],trkpointmin[i]));
	//TFile *fProd = TFile::Open(Form("M%d%d/mergedv_Prod2_eta.root",trkpointmax[t-1],trkpointmin[t-1]));
	TVectorD *vecDv2_Proderr = (TVectorD*)fProderr->Get(Form("D_%d/vmeanmean",xbin));
	TVectorD *vecDv2err_Proderr = (TVectorD*)fProderr->Get(Form("D_%d/sigmavmeanmean",xbin));
	TVectorD *vecDavgeta_Proderr = (TVectorD*)fProderr->Get(Form("D_%d/avgavgeta",xbin));
	
	//TVectorD *vecDv2_Prod = (TVectorD*)fProd->Get(Form("D_%d/vmean",xbin));
	//TVectorD *vecDv2err_Prod = (TVectorD*)fProd->Get(Form("D_%d/deltavmean",xbin));
	//TVectorD *vecDavgeta_Prod = (TVectorD*)fProd->Get(Form("D_%d/avgeta",xbin));

	double *avgeta_Proderr = vecDavgeta_Proderr->GetMatrixArray();
	double *v2_Proderr = vecDv2_Proderr->GetMatrixArray();
	double *v2err_Proderr = vecDv2err_Proderr->GetMatrixArray();

	//double *avgeta_Prod = vecDavgeta_Prod->GetMatrixArray();
	//double *v2_Prod = vecDv2_Prod->GetMatrixArray();
	//double *v2err_Prod = vecDv2err_Prod->GetMatrixArray();
	int neta = vecDavgeta_Proderr->GetNrows();
	
	c1->cd(i+1);
	//if(i!=ntotbin-1)
	//TGraphErrors *grProd=new TGraphErrors(neta,avgeta_Prod,v2_Prod,0,v2err_Prod);
	TGraphErrors *grProd=new TGraphErrors(neta,avgeta_Proderr,v2_Proderr,0,v2err_Proderr);
	//TGraphErrors *gr24=new TGraphErrors(neta24,eta,v24[i],0,v24err[i]);
	//TGraphErrors *gr26=new TGraphErrors(neta24,eta,v26[i],0,v26err[i]);
	//TGraphErrors *gr28=new TGraphErrors(neta24,eta,v28[i],0,v28err[i]);
	//gr24->SetMarkerSize(1.4);
	//gr24->SetMarkerColor(1);
	//gr24->SetMarkerStyle(20);
	//gr24->SetLineColor(1);
	//gr26->SetMarkerSize(1.4);
	//gr26->SetMarkerColor(4);
	//gr26->SetMarkerStyle(34);
	//gr26->SetLineColor(4);
	//gr28->SetMarkerSize(1.4);
	//gr28->SetMarkerColor(2);
	//gr28->SetMarkerStyle(33);
	//gr28->SetLineColor(2);
	grProd->SetMarkerSize(1.5);
	grProd->SetMarkerStyle(20);
	grProd->SetMarkerColor(2);
	grProd->SetLineColor(2);
	hFrame->Draw();
	//gr24->Draw("Psame");
	//gr26->Draw("Psame");
	//gr28->Draw("Psame");
	grProd->Draw("Psame");
	TLegend *tl = new TLegend(0.2,0.70,0.4,0.95);
	tl->SetFillColor(0);
	tl->SetBorderSize(0);
	tl->SetTextSize(0.055);
	tl->AddEntry(grProd,"v_{2}{LYZ}","lp");
	//tl->AddEntry(gr24,"v_{2}{4} cumulant","lp");
	//tl->AddEntry(gr26,"v_{2}{6} cumulant","lp");
	//tl->AddEntry(gr28,"v_{2}{8} cumulant","lp");
            TLatex *tlx0 = new TLatex(0.2,0.85,Form("CMS pPb #sqrt{s_{NN}} = 5.02TeV"));
	tlx0->SetNDC();
	tlx0->SetTextSize(0.065);
	if(i==0)
	tlx0->Draw("same");
	if(i==0) 
		TLatex *tlx2 = new TLatex(0.6,0.2,Form("%d<N_{trk}^{offline}<%d",trkpointmin[i],trkpointmax[i]));
	else
		TLatex *tlx2 = new TLatex(0.5,0.2,Form("%d<N_{trk}^{offline}<%d",trkpointmin[i],trkpointmax[i]));
	if(i==1) tl->Draw("same");
	tlx2->SetNDC();
	tlx2->SetTextSize(0.055);
	tlx2->Draw("same");
        for(int ieta=0; ieta<neta; ieta++){
        double px1 = avgeta_Proderr[ieta]-0.18;
        double py1 = v2_Proderr[ieta]*(1-sys);
        double px2 = avgeta_Proderr[ieta]+0.18;
        double py2 = v2_Proderr[ieta]*(1+sys);
        TBox *boxv2 = new TBox(px1,py1,px2,py2);
        boxv2->SetFillColor(16);
        boxv2->SetFillStyle(1001);
        boxv2->SetLineColor(16);
        boxv2->SetLineWidth(1);
        boxv2->Draw("lsame");
      }
        grProd->Draw("Psame");
	fProderr->Close();
        //TLatex *tlx1 = new TLatex(0.12,0.25,Form("%.1f<p_{T}<%.1f (GeV/c)",0.3,6.0));
	//tlx1->SetNDC();
	//tlx1->SetTextSize(0.045);
        }
	//tlx1->Draw("same");
	c1->Print("v2vseta_coarse_4p.png");
	c1->Print("v2vseta_coarse_4p.pdf");
}

