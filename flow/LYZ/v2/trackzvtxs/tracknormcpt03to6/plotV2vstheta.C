#include "M150120/par.h"
#include "/home/xuq7/HI/jetRpA/RpA/Quality/root_setting.h"

void plotV2vstheta(){
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
        TH1D *hFrame = new TH1D("","",220,-0.2,2);
        hFrame->SetTitle("");
        hFrame->GetXaxis()->SetTitle("#theta");
        hFrame->GetYaxis()->SetTitle("reference V_{2}");
        hFrame->GetYaxis()->SetTitleOffset(1.1);
        hFrame->GetXaxis()->SetTitleSize(0.04);
        hFrame->GetYaxis()->SetTitleSize(0.04);
        hFrame->GetXaxis()->SetRangeUser(-0.2,1.4);
        hFrame->SetMinimum(0.033);
        hFrame->SetMaximum(0.076);
        for(int trkbin=0;trkbin<ntotbin; trkbin++){
	if(isSum==0){
	    f = TFile::Open(Form("M%d%d/mergedV_Prod.root",trkpointmax[trkbin],trkpointmin[trkbin]));
	}
	else{
	    f = TFile::Open(Form("M%d%d/mergedV_Sum.root",trkpointmax[trkbin],trkpointmin[trkbin]));
	}
	TVectorD* vecVmean = (TVectorD*)f->Get(Form("Vmean",xbin));
        TVectorD* vecV = (TVectorD*)f->Get(Form("D_%d/V",xbin));
        TVectorD* vecdeltaV = (TVectorD*)f->Get(Form("D_%d/deltaV",xbin));

        double Vmean = (*vecVmean)[0];
        double *V = vecV->GetMatrixArray();
        double *deltaV = vecdeltaV->GetMatrixArray();
        double theta[ntheta];
        for(int itheta=0;itheta<ntheta;itheta++){
            theta[itheta]=itheta*TMath::Pi()/ntheta/nn;
        }
        int maxper10 = findmaxper(V,ntheta,Vmean);
        double maxper = (double)(maxper10+1)/10;
        c1->cd(trkbin+1);
        hFrame->Draw();
        TGraphErrors *gV2theta = new TGraphErrors(ntheta,theta,V,0,deltaV);
        gV2theta->SetMarkerStyle(20);
        gV2theta->SetMarkerSize(1.3);
        gV2theta->SetMarkerColor(1);
        gV2theta->SetLineColor(1);
	gV2theta->Draw("Psame");
        TLine *lup = new TLine(gV2theta->GetXaxis()->GetXmin(),Vmean*(1+maxper), gV2theta->GetXaxis()->GetXmax(),Vmean*(1+maxper));
        TLine *ldown = new TLine(gV2theta->GetXaxis()->GetXmin(),Vmean*(1-maxper), gV2theta->GetXaxis()->GetXmax(),Vmean*(1-maxper));
        TLine *l = new TLine(hFrame->GetXaxis()->GetXmin(),Vmean, hFrame->GetXaxis()->GetXmax(),Vmean);
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
        tl->DrawLatex(0,Vmean*(1+maxper),Form("mean up %.f%%",maxper*100));
        tl->DrawLatex(0,Vmean*(1-maxper),Form("mean down %.f%%",maxper*100));
        tl->SetNDC();
        tl->DrawLatex(0.5,0.85,Form("Multiplicity %d to %d",trkpointmin[trkbin],trkpointmax[trkbin]));
        l->Draw("same");
        lup->Draw("same");
        ldown->Draw("same");
        }
        c1->cd(ntotbin+1);
        TLatex *tlx0 = new TLatex(0.12,0.3,Form("track"));
        TLatex *tlx1 = new TLatex(0.12,0.25,Form("%.1f<p_{T}<%.1f (GeV/c)",0.3,6.0));
        tlx0->SetNDC();
        tlx1->SetNDC();
        tlx0->SetTextSize(0.045);
        tlx1->SetTextSize(0.045);
        hFrame->Draw();
        tlx0->Draw("same");
        tlx1->Draw("same");
	if(isSum==0)c1->SaveAs("hV2theta_Prod.png");
	else c1->SaveAs("hV2theta_Sum.png");
}

int findmaxper(double* x, int n, double xmean){
    double temp=0;
    for(int i=0;i<n;i++){
        if(TMath::Abs(x[i]-xmean)>temp)
            temp=TMath::Abs(x[i]-xmean);
    }
    return (int)(temp/xmean*10);
}
