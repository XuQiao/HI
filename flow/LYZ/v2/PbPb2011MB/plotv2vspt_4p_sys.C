#include "par.h"
#include "/home/xuq7/HI/jetRpA/RpA/Quality/root_setting.h"
void plotv2vspt_4p_sys(){

const int ntotbin=4;
const int npt24=9;
const double pt[npt24]={0.4,0.7,1.2,1.7,2.2,2.7,3.4,4.4,5.4};
const double v24[5][npt24]={
    {0.0436326,0.0795412,0.117455,0.135615,0.138952,0.137704,0.12716,0.104461,0.102546},
    {0.0441134,0.081454,0.120992,0.142609,0.148646,0.144631,0.13748,0.115206,0.102368},
    {0.0444633,0.0827859,0.124137,0.145269,0.149668,0.150845,0.144292,0.115638,0.127838},
    {0.045803,0.083701,0.127181,0.148425,0.151695,0.152144,0.135538,0.115144,0.0939622}
};
const double v24err[5][npt24]={
    {0.000639539,0.000826055,0.00117583,0.00173971,0.00147411,0.00425427,0.00290752,0.0075635,0.0180777},
    {0.000489939,0.000748004,0.00106544,0.00130725,0.0021476,0.00300211,0.00460689,0.00430641,0.00958347},
    {0.000539929,0.000592858,0.00159546,0.001461,0.00161371,0.00332517,0.00381725,0.00692403,0.0137116},
    {0.000685138,0.0011504,0.00182689,0.00135556,0.00211262,0.00191136,0.00361452,0.00609861,0.0112248}
};
const double v26[5][npt24]={
    {0.0436647,0.0766372,0.115253,0.134962,0.138116,0.132017,0.129223,0.0959764,0.159984},
    {0.0427129,0.07916,0.119421,0.142728,0.147721,0.142104,0.138136,0.111479,0.0930553},
    {0.0439986,0.0812754,0.121722,0.141086,0.144929,0.14576,0.141362,0.103954,0.118998},
    {0.045338,0.081817,0.124592,0.145770,0.148010,0.151483,0.135710,0.115941,0.088239}
};
const double v26err[5][npt24]={
    {0.00117633,0.00172456,0.00243925,0.00331087,0.00340786,0.00552076,0.0053438,0.00634,0.023902},
    {0.000714298,0.00115609,0.00178786,0.00201899,0.00247538,0.0042914,0.00441189,0.00375012,0.00991865},
    {0.0007454,0.00108873,0.00211174,0.00242773,0.002964,0.00406841,0.00325586,0.0105999,0.0164428},
    {0.000950693,0.00149145,0.00207434,0.00228235,0.0025587,0.00320457,0.00459367,0.00899185,0.00916074}
};
const double v28[5][npt24]={
    {0.0413359,0.0745479,0.111287,0.1315,0.13296,0.120331,0.132252,0.105274,0.158618},
    {0.0424299,0.0761756,0.115879,0.137643,0.138677,0.13397,0.125922,0.103297,0.115589},
    {0.0424748,0.079193,0.118646,0.143186,0.147939,0.139611,0.137666,0.118945,0.0998064},
    {0.044857,0.081403,0.122831,0.143983,0.145980,0.150662,0.140130,0.104239,0.101058}
};
const double v28err[5][npt24]={
    {0.00174636,0.00232066,0.00378455,0.00539576,0.00686651,0.0070001,0.00972887,0.0107575,0.0177152},
    {0.0014467,0.00177893,0.00310776,0.00351237,0.00391883,0.00406891,0.00434652,0.00735751,0.0165999},
    {0.000964082,0.00149995,0.00232247,0.00286737,0.00303328,0.00467016,0.00455719,0.00512068,0.00835018},
    {0.000613315,0.00112984,0.0017138,0.00182123,0.00204996,0.0027923,0.00233936,0.00455531,0.00993196}
};


c1 = new TCanvas("c1"," ",1200,350);
makeMultiPanelCanvas(c1,4,1,0,0,0.15,0.15,0.01);
c2 = new TCanvas("c2"," ",1200,350);
makeMultiPanelCanvas(c2,4,1,0,0,0.15,0.15,0.01);
c3 = new TCanvas("c3"," ",1200,350);
makeMultiPanelCanvas(c3,4,1,0,0,0.15,0.15,0.01);
c1_ = new TCanvas("c1_"," ",1200,350);
makeMultiPanelCanvas(c1_,4,1,0,0,0.15,0.15,0.01);
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
int shift=2;
for(int i=shift;i<shift+ntotbin;i++){
        int t = ntotbin-(i-shift);
	//TFile *fProderr = TFile::Open(Form("M%d%d/mergedv_Prod2_sub.root",trkpointmax[i],trkpointmin[i]));
	TFile *fProd = TFile::Open(Form("mergedv_Prod2.root"));
	TFile *fProdtc = TFile::Open(Form("../tracktc/PbPb2011MB/mergedv_Prod2.root"));
	TFile *fProdlc = TFile::Open(Form("../tracklc/PbPb2011MB/mergedv_Prod2.root"));
	TFile *fProdzl = TFile::Open(Form("../trackzvtxl/PbPb2011MB/mergedv_Prod2.root"));
	TFile *fProdzs = TFile::Open(Form("../trackzvtxs/PbPb2011MB/mergedv_Prod2.root"));
	TFile *fProdpu = TFile::Open(Form("../PUrej/PbPb2011MB/mergedv_Prod2.root"));
	//TVectorD *vecDv2_Proderr = (TVectorD*)fProderr->Get(Form("D_%d/vmeanmean",ibin));
	//TVectorD *vecDv2err_Proderr = (TVectorD*)fProderr->Get(Form("D_%d/sigmavmeanmean",ibin));
	//TVectorD *vecDavgpt_Proderr = (TVectorD*)fProderr->Get(Form("D_%d/avgavgpt",ibin));
	
	TVectorD *vecDv2_Prod = (TVectorD*)fProd->Get(Form("D_%d/vmean",i));
	TVectorD *vecDv2err_Prod = (TVectorD*)fProd->Get(Form("D_%d/deltavmean",i));
	TVectorD *vecDavgpt_Prod = (TVectorD*)fProd->Get(Form("D_%d/avgpt",i));

	TVectorD *vecDv2_Prodtc = (TVectorD*)fProdtc->Get(Form("D_%d/vmean",i));
	TVectorD *vecDv2err_Prodtc = (TVectorD*)fProdtc->Get(Form("D_%d/deltavmean",i));
	TVectorD *vecDavgpt_Prodtc = (TVectorD*)fProdtc->Get(Form("D_%d/avgpt",i));

	TVectorD *vecDv2_Prodlc = (TVectorD*)fProdlc->Get(Form("D_%d/vmean",i));
	TVectorD *vecDv2err_Prodlc = (TVectorD*)fProdlc->Get(Form("D_%d/deltavmean",i));
	TVectorD *vecDavgpt_Prodlc = (TVectorD*)fProdlc->Get(Form("D_%d/avgpt",i));
	
        TVectorD *vecDv2_Prodzl = (TVectorD*)fProdzl->Get(Form("D_%d/vmean",i));
	TVectorD *vecDv2err_Prodzl = (TVectorD*)fProdzl->Get(Form("D_%d/deltavmean",i));
	TVectorD *vecDavgpt_Prodzl = (TVectorD*)fProdzl->Get(Form("D_%d/avgpt",i));

	TVectorD *vecDv2_Prodzs = (TVectorD*)fProdzs->Get(Form("D_%d/vmean",i));
	TVectorD *vecDv2err_Prodzs = (TVectorD*)fProdzs->Get(Form("D_%d/deltavmean",i));
	TVectorD *vecDavgpt_Prodzs = (TVectorD*)fProdzs->Get(Form("D_%d/avgpt",i));
	
        TVectorD *vecDv2_Prodpu = (TVectorD*)fProdpu->Get(Form("D_%d/vmean",i));
	TVectorD *vecDv2err_Prodpu = (TVectorD*)fProdpu->Get(Form("D_%d/deltavmean",i));
	TVectorD *vecDavgpt_Prodpu = (TVectorD*)fProdpu->Get(Form("D_%d/avgpt",i));

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
	double *v2_Prodpu = vecDv2_Prodpu->GetMatrixArray();
	double *v2err_Prodpu = vecDv2err_Prodpu->GetMatrixArray();
	int npt = vecDavgpt_Prod->GetNrows();
        
	//if(i!=ntotbin-1)
	TGraphErrors *grProd=new TGraphErrors(npt,avgpt_Prod,v2_Prod,0,v2err_Prod);
	TGraphErrors *grProdtc=new TGraphErrors(npt,avgpt_Prod,v2_Prodtc,0,v2err_Prodtc);
	TGraphErrors *grProdlc=new TGraphErrors(npt,avgpt_Prod,v2_Prodlc,0,v2err_Prodlc);
	TGraphErrors *grProdzl=new TGraphErrors(npt,avgpt_Prod,v2_Prodzl,0,v2err_Prodzl);
	TGraphErrors *grProdzs=new TGraphErrors(npt,avgpt_Prod,v2_Prodzs,0,v2err_Prodzs);
	TGraphErrors *grProdpu=new TGraphErrors(npt,avgpt_Prod,v2_Prodpu,0,v2err_Prodpu);
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
	grProdpu->SetMarkerSize(1.5);
	grProdpu->SetMarkerStyle(33);
	grProdpu->SetMarkerColor(4);
	grProdpu->SetLineColor(4);

	TGraphErrors *grProdr=(TGraphErrors*)grProd->Clone("grProdr");
	TGraphErrors *grProdtcr=(TGraphErrors*)grProdtc->Clone("grProdtcr");
	TGraphErrors *grProdlcr=(TGraphErrors*)grProdlc->Clone("grProdlcr");
	TGraphErrors *grProdzlr=(TGraphErrors*)grProdzl->Clone("grProdzlr");
	TGraphErrors *grProdzsr=(TGraphErrors*)grProdzs->Clone("grProdzsr");
	TGraphErrors *grProdpur=(TGraphErrors*)grProdpu->Clone("grProdpur");
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
            double y_sys = grProdpu->GetY()[j];
            grProdpur->SetPoint(j,x,y_sys/y);
            grProdpur->SetPointError(j,0,0);
            grProdr->SetPoint(j,x,1.);
            grProdr->SetPointError(j,0,0);
        }
	c1->cd(t);
        hFrame->SetMaximum(0.35);
        hFrame->SetMinimum(0);
        hFrame->GetYaxis()->SetTitle("v_{2}");	
	hFrame->DrawCopy();
	grProd->Draw("Psame");
	grProdtc->Draw("Psame");
	grProdlc->Draw("Psame");
	TLegend *tl = new TLegend(0.2,0.70,0.4,0.95);
	tl->SetFillColor(0);
	tl->SetBorderSize(0);
	tl->SetTextSize(0.055);
	tl->AddEntry(grProd,"v_{2}{LYZ}","lp");
	tl->AddEntry(grProdtc,"v_{2}{LYZ} loose cut ","lp");
	tl->AddEntry(grProdlc,"v_{2}{LYZ} tight cut","lp");
        TLatex *tlx0 = new TLatex(0.2,0.85,Form("CMS PbPb #sqrt{s_{NN}} = 2.76TeV"));
	tlx0->SetNDC();
	tlx0->SetTextSize(0.065);
	if(t==1)
	tlx0->Draw("same");
	if(t==1) 
		TLatex *tlx2 = new TLatex(0.6,0.2,Form("%d<N_{trk}^{offline}<%d",trkbin[i+1],trkbin[i]));
	else
		TLatex *tlx2 = new TLatex(0.5,0.2,Form("%d<N_{trk}^{offline}<%d",trkbin[i+1],trkbin[i]));
	if(t==2) tl->Draw("same");
	tlx2->SetNDC();
	tlx2->SetTextSize(0.055);
	tlx2->Draw("same");

        c2->cd(t);
        hFrame->SetMaximum(0.35);
        hFrame->SetMinimum(0);
        hFrame->GetYaxis()->SetTitle("v_{2}");	
	hFrame->DrawCopy();
	grProd->Draw("Psame");
	grProdzl->Draw("Psame");
	grProdzs->Draw("Psame");
	TLegend *tl = new TLegend(0.2,0.70,0.4,0.95);
	tl->SetFillColor(0);
	tl->SetBorderSize(0);
	tl->SetTextSize(0.055);
	tl->AddEntry(grProd,"v_{2}{LYZ}","lp");
	tl->AddEntry(grProdzl,"v_{2}{LYZ} 3 < zvertex < 15 ","lp");
	tl->AddEntry(grProdzs,"v_{2}{LYZ} zvertex < 3","lp");
        TLatex *tlx0 = new TLatex(0.2,0.85,Form("CMS PbPb #sqrt{s_{NN}} = 2.76TeV"));
	tlx0->SetNDC();
	tlx0->SetTextSize(0.065);
	if(t==1)
	tlx0->Draw("same");
	if(t==1) 
		TLatex *tlx2 = new TLatex(0.6,0.2,Form("%d<N_{trk}^{offline}<%d",trkbin[i+1],trkbin[i]));
	else
		TLatex *tlx2 = new TLatex(0.5,0.2,Form("%d<N_{trk}^{offline}<%d",trkbin[i+1],trkbin[i]));
	if(t==2) tl->Draw("same");
	tlx2->SetNDC();
	tlx2->SetTextSize(0.055);
	tlx2->Draw("same");

        c3->cd(t);
        hFrame->SetMaximum(0.35);
        hFrame->SetMinimum(0);
        hFrame->GetYaxis()->SetTitle("v_{2}");	
	hFrame->DrawCopy();
	grProd->Draw("Psame");
	grProdpu->Draw("Psame");
	tl->AddEntry(grProdzs,"v_{2}{LYZ} zvertex < 3","lp");
        TLatex *tlx0 = new TLatex(0.2,0.85,Form("CMS PbPb #sqrt{s_{NN}} = 2.76TeV"));
	tlx0->SetNDC();
	tlx0->SetTextSize(0.065);
	if(t==1)
	tlx0->Draw("same");
	if(t==1) 
		TLatex *tlx2 = new TLatex(0.6,0.2,Form("%d<N_{trk}^{offline}<%d",trkbin[i+1],trkbin[i]));
	else
		TLatex *tlx2 = new TLatex(0.5,0.2,Form("%d<N_{trk}^{offline}<%d",trkbin[i+1],trkbin[i]));
	if(t==2) tl->Draw("same");
	tlx2->SetNDC();
	tlx2->SetTextSize(0.055);
	tlx2->Draw("same");

        c1_->cd(t);
        hFrame->SetMaximum(1.2);
        hFrame->SetMinimum(0.8);
        hFrame->GetYaxis()->SetTitle("v_{2} systematics");	
	hFrame->DrawCopy();
	grProdr->Draw("Psame");
	grProdtcr->Draw("Psame");
	grProdlcr->Draw("Psame");
	grProdzlr->Draw("Psame");
	grProdzsr->Draw("Psame");
	if(t==1)
	tlx0->Draw("same");
	if(t==2) tl->Draw("same");
	tlx2->Draw("same");
	fProd->Close();
        }
	c1->Print("v2vspt_4p_sys.png");
	c1->Print("v2vspt_4p_sys.pdf");
	c2->Print("v2vspt_4p_sys.png");
	c2->Print("v2vspt_4p_sys.pdf");
	c3->Print("v2vspt_4p_sys.png");
	c3->Print("v2vspt_4p_sys.pdf");
	c1_->Print("v2vspt_4p_sysr.png");
	c1_->Print("v2vspt_4p_sysr.pdf");
}

