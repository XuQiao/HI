#include "/home/xuq7/HI/jetRpA/RpA/Quality/root_setting.h"
const int ntotbin=5;
const int ndir = 6;
const int trkpointmin[ntotbin] = {120,150,185,220,260};
const int trkpointmax[ntotbin] = {150,185,220,260,300};
const int marker[ndir] = {20,24,27,30,31,29};
const int color[ndir] = {2,3,4,6,7,8};
const int ntheta=5;
void plotr0vstheta_Ntrk(){

const int nbin24 = 12;
const double avgtrkbin[nbin24]={44.36,54.33,68.63,88.39,108.2,131.3,162.1,196.6,227.5,247.2,269.2,301.2};
const double V24[nbin24]={0.02965,0.03913,0.04832,0.04941,0.04822,0.04955,0.049,0.04805,0.04709,0.04665,0.04772,0.04797};
const double V24err[nbin24]={0.005967,0.004262,0.001496,0.001351,0.001599,0.0003922,0.0003065,0.0002939,0.0004568,0.0008684,0.001471,0.004329};
const int startpoint = 5;


c1 = new TCanvas("c1"," ",1200,800);
makeMultiPanelCanvas(c1,3,2,0,0,0.25,0.2,0.03);
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetErrorX(0);
    TH1D *hFrame = new TH1D("","",300,-1,2);
    hFrame->SetTitle("");
    hFrame->GetXaxis()->SetTitle("#theta");
    hFrame->GetYaxis()->SetTitle("r_{0}");	
    hFrame->GetXaxis()->SetTitleSize(0.04);
    hFrame->GetYaxis()->SetTitleSize(0.04);
    hFrame->GetXaxis()->SetRangeUser(-0.1,1.5);
    hFrame->SetMaximum(0.45);
    hFrame->SetMinimum(0.01);
    TGraphErrors* gr[ndir];
for(int i=0;i<ntotbin;i++){
	c1->cd(i+1);
	hFrame->Draw();
        for(int iloop = 0; iloop<ndir; iloop++){
    	gr[iloop]=plot(iloop,color[iloop],marker[iloop],i);
   	gr[iloop]->Draw("Psame");
        }
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
        for(int iloop = 0; iloop<ndir; iloop++)
            tl->AddEntry(gr[iloop],Form("iloop = %d", iloop),"lp");
        tl->SetFillColor(0);
        tl->SetBorderSize(0);
        tl->SetTextSize(0.065);

        tlx2->SetTextSize(0.065);
        hFrame->Draw();
        //tlx1->Draw("same");
        tl->Draw("same");
        c1->Print("r0vstheta_Ntrk.png");
}

TGraphErrors* plot(int iloop, int color, int marker, int i){

    TFile *f = TFile::Open(Form("M%d%d/parin_%d.root",trkpointmax[i],trkpointmin[i],iloop));
    int xbin=0;
    int xpt=0;
    const int nn = 2;
        //	TVectorD *vecNtrk = (TVectorD*)f->Get("avgtrk");
       // TVectorD *vecavgmult = (TVectorD*)f->Get("avgmultall");
       // TVectorD *vecr0mean_0=(TVectorD*)f->Get(Form("D_%d/Vmean",xbin));
        TVectorD *vecr0_0=(TVectorD*)f->Get(Form("D_%d/D_%d/r0",xbin,xpt));
       // TVectorD *vecr0err_0=(TVectorD*)f->Get(Form("D_%d/D_%d/deltaV",xbin,xpt));
        double *r0_0=vecr0_0->GetMatrixArray();
        //double *r0err_0=vecr0err_0->GetMatrixArray();
        double theta[ntheta];
        for(int itheta=0;itheta<ntheta;itheta++){
            theta[itheta]=itheta*TMath::Pi()/ntheta/nn;
        }

        gr=new TGraphErrors(ntheta,theta,r0_0,0,0);
	
	gr->SetMarkerStyle(marker);
	gr->SetMarkerSize(1.6);
	gr->SetMarkerColor(color);
	gr->SetLineColor(color);
	f->Close();
	return gr;
}

