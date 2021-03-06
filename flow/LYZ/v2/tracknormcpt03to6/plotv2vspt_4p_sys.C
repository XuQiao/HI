#include "par.h"
#include "/home/xuq7/HI/jetRpA/RpA/Quality/root_setting.h"
void plotv2vspt_4p_sys(){

const int ntotbin=5;
const int trkpointmin[ntotbin] = {120,150,185,220,260};
const int trkpointmax[ntotbin] = {150,185,220,260,300};
const int xbin=0;
c1 = new TCanvas("c1"," ",1200,350);
makeMultiPanelCanvas(c1,4,1,0,0,0.15,0.15,0.01);
c2 = new TCanvas("c2"," ",1200,350);
makeMultiPanelCanvas(c2,4,1,0,0,0.15,0.15,0.01);
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetErrorX(0);
    TH1D *hFrame = new TH1D("","",80,-1,7);
    hFrame->SetTitle("");
    hFrame->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hFrame->GetXaxis()->SetTitleSize(0.06);
    hFrame->GetYaxis()->SetTitleSize(0.06);
    hFrame->GetXaxis()->SetLabelSize(0.06);
    hFrame->GetYaxis()->SetLabelSize(0.06);
    hFrame->GetYaxis()->SetNdivisions(505);
    hFrame->GetXaxis()->CenterTitle();
    hFrame->GetYaxis()->CenterTitle();
    hFrame->GetXaxis()->SetRangeUser(-0.3,6.3);
for(int i=0;i<ntotbin;i++){
	TFile *fProd = TFile::Open(Form("M%d%d/mergedv_Prod2_sub.root",trkpointmax[i],trkpointmin[i],trkpointmax[i],trkpointmin[i]));
	TFile *fProdtc = TFile::Open(Form("../tracktc/tracknormcpt03to6/M%d%d/mergedv_Prod2_sub.root",trkpointmax[i],trkpointmin[i]));
	TFile *fProdlc = TFile::Open(Form("../tracklc/tracknormcpt03to6/M%d%d/mergedv_Prod2_sub.root",trkpointmax[i],trkpointmin[i]));
	TFile *fProdzl = TFile::Open(Form("../trackzvtxl/tracknormcpt03to6/M%d%d/mergedv_Prod2_sub.root",trkpointmax[i],trkpointmin[i]));
	TFile *fProdzs = TFile::Open(Form("../trackzvtxs/tracknormcpt03to6/M%d%d/mergedv_Prod2_sub.root",trkpointmax[i],trkpointmin[i]));
	//TVectorD *vecDv2_Proderr = (TVectorD*)fProderr->Get(Form("D_%d/vmeanmean",ibin));
	//TVectorD *vecDv2err_Proderr = (TVectorD*)fProderr->Get(Form("D_%d/sigmavmeanmean",ibin));
	//TVectorD *vecDavgpt_Proderr = (TVectorD*)fProderr->Get(Form("D_%d/avgavgavgpt",ibin));
	
	TVectorD *vecDv2_Prod = (TVectorD*)fProd->Get(Form("D_%d/vmeanmean",xbin));
	TVectorD *vecDv2err_Prod = (TVectorD*)fProd->Get(Form("D_%d/sigmavmeanmean",xbin));
	TVectorD *vecDavgpt_Prod = (TVectorD*)fProd->Get(Form("D_%d/avgavgpt",xbin));

	TVectorD *vecDv2_Prodtc = (TVectorD*)fProdtc->Get(Form("D_%d/vmeanmean",xbin));
	TVectorD *vecDv2err_Prodtc = (TVectorD*)fProdtc->Get(Form("D_%d/sigmavmeanmean",xbin));
	TVectorD *vecDavgpt_Prodtc = (TVectorD*)fProdtc->Get(Form("D_%d/avgavgpt",xbin));

	TVectorD *vecDv2_Prodlc = (TVectorD*)fProdlc->Get(Form("D_%d/vmeanmean",xbin));
	TVectorD *vecDv2err_Prodlc = (TVectorD*)fProdlc->Get(Form("D_%d/sigmavmeanmean",xbin));
	TVectorD *vecDavgpt_Prodlc = (TVectorD*)fProdlc->Get(Form("D_%d/avgavgpt",xbin));
	
        TVectorD *vecDv2_Prodzl = (TVectorD*)fProdzl->Get(Form("D_%d/vmeanmean",xbin));
	TVectorD *vecDv2err_Prodzl = (TVectorD*)fProdzl->Get(Form("D_%d/sigmavmeanmean",xbin));
	TVectorD *vecDavgpt_Prodzl = (TVectorD*)fProdzl->Get(Form("D_%d/avgavgpt",xbin));

	TVectorD *vecDv2_Prodzs = (TVectorD*)fProdzs->Get(Form("D_%d/vmeanmean",xbin));
	TVectorD *vecDv2err_Prodzs = (TVectorD*)fProdzs->Get(Form("D_%d/sigmavmeanmean",xbin));
	TVectorD *vecDavgpt_Prodzs = (TVectorD*)fProdzs->Get(Form("D_%d/avgavgpt",xbin));

	//double *avgpt_Proderr = vecDavgpt_Proderr->GetMatrixArray();
	//double *v2_Proderr = vecDv2_Proderr->GetMatrixArray();
	//double *v2err_Proderr = vecDv2err_Proderr->GetMatrixArray();

	double *avgpt_Prod = vecDavgpt_Prod->GetMatrixArray();
	double *v2_Prod = vecDv2_Prod->GetMatrixArray();
	double *v2err_Prod = vecDv2err_Prod->GetMatrixArray();
	double *v2_Prodtc = vecDv2_Prodtc->GetMatrixArray();
	double *v2err_Prodtc = vecDv2err_Prodtc->GetMatrixArray();
	double *v2_Prodlc = vecDv2_Prodlc->GetMatrixArray();
	double *v2err_Prodlc = vecDv2err_Prodlc->GetMatrixArray();
	double *v2_Prodzl = vecDv2_Prodzl->GetMatrixArray();
	double *v2err_Prodzl = vecDv2err_Prodzl->GetMatrixArray();
	double *v2_Prodzs = vecDv2_Prodzs->GetMatrixArray();
	double *v2err_Prodzs = vecDv2err_Prodzs->GetMatrixArray();
	int npt = vecDavgpt_Prod->GetNrows();
        
	//if(i!=ntotbin-1)
	TGraphErrors *grProd=new TGraphErrors(npt,avgpt_Prod,v2_Prod,0,v2err_Prod);
	TGraphErrors *grProdtc=new TGraphErrors(npt,avgpt_Prod,v2_Prodtc,0,v2err_Prodtc);
	TGraphErrors *grProdlc=new TGraphErrors(npt,avgpt_Prod,v2_Prodlc,0,v2err_Prodlc);
	TGraphErrors *grProdzl=new TGraphErrors(npt,avgpt_Prod,v2_Prodzl,0,v2err_Prodzl);
	TGraphErrors *grProdzs=new TGraphErrors(npt,avgpt_Prod,v2_Prodzs,0,v2err_Prodzs);
	grProd->SetMarkerSize(1.5);
	grProd->SetMarkerStyle(20);
	grProd->SetMarkerColor(1);
	grProd->SetLineColor(1);
	grProdtc->SetMarkerSize(1.5);
	grProdtc->SetMarkerStyle(20);
	grProdtc->SetMarkerColor(2);
	grProdtc->SetLineColor(2);
	grProdlc->SetMarkerSize(1.5);
	grProdlc->SetMarkerStyle(24);
	grProdlc->SetMarkerColor(2);
	grProdlc->SetLineColor(2);
	grProdzl->SetMarkerSize(1.5);
	grProdzl->SetMarkerStyle(34);
	grProdzl->SetMarkerColor(4);
	grProdzl->SetLineColor(4);
	grProdzs->SetMarkerSize(1.5);
	grProdzs->SetMarkerStyle(33);
	grProdzs->SetMarkerColor(4);
	grProdzs->SetLineColor(4);
	TGraphErrors *grProdr=(TGraphErrors*)grProd->Clone("grProdr");
	TGraphErrors *grProdtcr=(TGraphErrors*)grProdtc->Clone("grProdtcr");
	TGraphErrors *grProdlcr=(TGraphErrors*)grProdlc->Clone("grProdlcr");
	TGraphErrors *grProdzlr=(TGraphErrors*)grProdzl->Clone("grProdzlr");
	TGraphErrors *grProdzsr=(TGraphErrors*)grProdzs->Clone("grProdzsr");
        for(int j = 0;j<grProd->GetN();j++){
            double x = grProd->GetX()[j];
            double y = grProd->GetY()[j];
            double y_sys = grProdtc->GetY()[j];
            grProdtcr->SetPoint(j,x,y_sys/y);
            grProdtcr->SetPointError(j,0,0);
            double y_sys = grProdlc->GetY()[j];
            grProdlcr->SetPoint(j,x,y_sys/y);
            grProdlcr->SetPointError(j,0,0);
            double y_sys = grProdzl->GetY()[j];
            grProdzlr->SetPoint(j,x,y_sys/y);
            grProdzlr->SetPointError(j,0,0);
            double y_sys = grProdzs->GetY()[j];
            grProdzsr->SetPoint(j,x,y_sys/y);
            grProdzsr->SetPointError(j,0,0);
            grProdr->SetPoint(j,x,1.);
            grProdr->SetPointError(j,0,0);
        }
	c1->cd(i+1);
        hFrame->SetMaximum(0.35);
        hFrame->SetMinimum(0);
        hFrame->GetYaxis()->SetTitle("v_{2}");	
	hFrame->DrawCopy();
	grProd->Draw("Psame");
	//grProdtc->Draw("Psame");
	//grProdlc->Draw("Psame");
	grProdzl->Draw("Psame");
	grProdzs->Draw("Psame");
	TLegend *tl = new TLegend(0.2,0.70,0.4,0.95);
	tl->SetFillColor(0);
	tl->SetBorderSize(0);
	tl->SetTextSize(0.055);
	tl->AddEntry(grProd,"v_{2}{LYZ}","lp");
	//tl->AddEntry(grProdtc,"v_{2}{LYZ} loose cut ","lp");
	//tl->AddEntry(grProdlc,"v_{2}{LYZ} tight cut","lp");
	tl->AddEntry(grProdzl,"v_{2}{LYZ} 3 < zvertex < 15 ","lp");
	//tl->AddEntry(grProdzs,"v_{2}{LYZ} zvertex < 3","lp");
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
        c2->cd(i+1);
        hFrame->SetMaximum(1.2);
        hFrame->SetMinimum(0.8);
        hFrame->GetYaxis()->SetTitle("v_{2} systematics");	
	hFrame->DrawCopy();
	grProdr->Draw("Psame");
	grProdtcr->Draw("Psame");
	grProdlcr->Draw("Psame");
	grProdzlr->Draw("Psame");
	grProdzsr->Draw("Psame");
	if(i==0)
	tlx0->Draw("same");
	if(i==1) tl->Draw("same");
	tlx2->Draw("same");
	fProd->Close();
        }
	c1->Print("v2vspt_4p_sys.png");
	c1->Print("v2vspt_4p_sys.pdf");
	c2->Print("v2vspt_4p_sysr.png");
	c2->Print("v2vspt_4p_sysr.pdf");
}

