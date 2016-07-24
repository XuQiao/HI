#include "M150120/par.h"
#include "/home/xuq7/HI/jetRpA/RpA/Quality/root_setting.h"

void plotr0vstheta(){
	TFile *f;
        int isSum=0;
        const int ntotbin=5;
        const int trkpointmin[ntotbin] = {120,150,185,220,260};
        const int trkpointmax[ntotbin] = {150,185,220,260,300};
        int trkbin=1;
	int xbin=0;
        c1 = new TCanvas("c1"," ",1200,700);
        makeMultiPanelCanvas(c1,3,2,0,0,0.25,0.2,0.03);
        gStyle->SetOptFit(1);
        gStyle->SetOptStat(0);
        gStyle->SetOptTitle(0);
        gStyle->SetErrorX(0);
        TH1D *hFrame = new TH1D("","",100,0,2);
        hFrame->SetTitle("");
        hFrame->GetXaxis()->SetTitle("#theta");
        hFrame->GetYaxis()->SetTitle("r_{0}");
        hFrame->GetYaxis()->SetTitleOffset(1.1);
        hFrame->GetXaxis()->SetTitleSize(0.04);
        hFrame->GetYaxis()->SetTitleSize(0.04);
        hFrame->GetXaxis()->SetRangeUser(0,1.5);
        hFrame->SetMinimum(0.05);
        hFrame->SetMaximum(0.30);
        for(int trkbin=0;trkbin<ntotbin; trkbin++){
	if(isSum==0){
	    f = TFile::Open(Form("M%d%d/mergedV_Prod.root",trkpointmax[trkbin],trkpointmin[trkbin]));
	}
	else{
	    f = TFile::Open(Form("M%d%d/mergedV_Sum.root",trkpointmax[trkbin],trkpointmin[trkbin]));
	}
	TVectorD* vecr0 = (TVectorD*)f->Get(Form("D_%d/D_0/r0",xbin));

        double *r0 = vecr0->GetMatrixArray();
        double r0mean = getmean(r0,ntheta);
        double theta[ntheta];
        for(int itheta=0;itheta<ntheta;itheta++){
            theta[itheta]=itheta*TMath::Pi()/ntheta/nn;
        }
        int maxper10 = findmaxper(r0,ntheta,r0mean);
        double maxper = (double)(maxper10+1)/10;
        c1->cd(trkbin+1);
        hFrame->Draw();
        TGraph *gr0theta = new TGraph(ntheta,theta,r0);
        gr0theta->SetMarkerStyle(20);
        gr0theta->SetMarkerSize(1.3);
        gr0theta->SetMarkerColor(1);
        gr0theta->SetLineColor(1);
	gr0theta->Draw("Psame");
        TLine *lup = new TLine(gr0theta->GetXaxis()->GetXmin(),r0mean*(1+maxper), gr0theta->GetXaxis()->GetXmax(),r0mean*(1+maxper));
        TLine *ldown = new TLine(gr0theta->GetXaxis()->GetXmin(),r0mean*(1-maxper), gr0theta->GetXaxis()->GetXmax(),r0mean*(1-maxper));
        TLine *l = new TLine(hFrame->GetXaxis()->GetXmin(),r0mean, hFrame->GetXaxis()->GetXmax(),r0mean);
        l->SetLineStyle(2);
        lup->SetLineStyle(2);
        ldown->SetLineStyle(2);
        l->SetLineWidth(1.2);
        lup->SetLineWidth(1.2);
        ldown->SetLineWidth(1.2);
        TLatex *tl = new TLatex();
 //     tl->SetNDC();
        tl->SetTextFont(42);
        tl->SetTextSize(0.04);
 //     tl->SetBorderStyle(0);
        tl->DrawLatex(0,r0mean*(1+maxper),Form("mean up %.f%%",maxper*100));
        tl->DrawLatex(0,r0mean*(1-maxper),Form("mean down %.f%%",maxper*100));
        tl->SetNDC();
        tl->DrawLatex(0.7,0.85,Form("Multiplicity %d to %d",trkpointmin[trkbin],trkpointmax[trkbin]));
        l->Draw("same");
        lup->Draw("same");
        ldown->Draw("same");
        }
        c1->cd(ntotbin+1);
        TLatex *tlx0 = new TLatex(0.12,0.3,Form("PF candi"));
        TLatex *tlx1 = new TLatex(0.12,0.25,Form("%.1f<p_{T}<%.1f (GeV/c)",0.3,6.0));
        tlx0->SetNDC();
        tlx1->SetNDC();
        tlx0->SetTextSize(0.045);
        tlx1->SetTextSize(0.045);
        hFrame->Draw();
        tlx0->Draw("same");
        tlx1->Draw("same");
	if(isSum==0)c1->SaveAs("hr0theta_Prod.png");
	else c1->SaveAs("hr0theta_Sum.png");
}

int findmaxper(double* x, int n, double xmean){
    double temp=0;
    for(int i=0;i<n;i++){
        if(TMath::Abs(x[i]-xmean)>temp)
            temp=TMath::Abs(x[i]-xmean);
    }
    return (int)(temp/xmean*10);
}

double getmean(double*x, int n){
    double temp=0;
    for(int i=0;i<n;i++){
        temp+=x[i];
    }
    return temp/n;
}
