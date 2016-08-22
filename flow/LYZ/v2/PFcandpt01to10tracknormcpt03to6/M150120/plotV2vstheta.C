#include "par.h"
void plotV2vstheta(){
	TFile *f;
        int isSum=0;
	if(isSum==0){
	    f = TFile::Open("mergedV_Prod.root");
	}
	else{
	    f = TFile::Open("mergedV_Sum.root");
	}
	int xbin=0;

	TVectorD* vecVmean = (TVectorD*)f->Get(Form("D_%d/Vmean",xbin));
        TVectorD* vecV = (TVectorD*)f->Get(Form("D_%d/D_0/V",xbin));

        double Vmean = (*vecVmean)[0];
        double *V = vecV->GetMatrixArray();
        double theta[ntheta];
        for(int itheta=0;itheta<ntheta;itheta++){
            theta[itheta]=itheta*TMath::Pi()/ntheta/nn;

        }
        double maxper = findmaxper(V,ntheta,Vmean);
        cout<<maxper<<endl;
        TGraph *gV2theta = new TGraph(ntheta,theta,V);
        gV2theta->SetMarkerStyle(20);
        gV2theta->SetMarkerSize(1.3);
        gV2theta->SetMarkerColor(1);
        gV2theta->SetLineColor(1);
        gV2theta->SetMinimum(Vmean*(1-maxper*2));
        gV2theta->SetMaximum(Vmean*(1+maxper*2));
	gV2theta->Draw("AP");
        TLine *lup = new TLine(gV2theta->GetXaxis()->GetXmin(),Vmean*(1+maxper*1.2), gV2theta->GetXaxis()->GetXmax(),Vmean*(1+maxper*1.2));
        TLine *ldown = new TLine(gV2theta->GetXaxis()->GetXmin(),Vmean*(1-maxper*1.2), gV2theta->GetXaxis()->GetXmax(),Vmean*(1-maxper*1.2));
        TLine *l = new TLine(gV2theta->GetXaxis()->GetXmin(),Vmean, gV2theta->GetXaxis()->GetXmax(),Vmean);
        l->SetLineStyle(2);
        lup->SetLineStyle(2);
        ldown->SetLineStyle(2);
        l->SetLineWidth(1.2);
        lup->SetLineWidth(1.2);
        ldown->SetLineWidth(1.2);
        l->Draw("same");
        lup->Draw("same");
        ldown->Draw("same");
	if(isSum==0)c1->SaveAs("hV2theta_Prod.png");
	else c1->SaveAs("hV2theta_Sum.png");
}

double findmaxper(double* x, int n, double xmean){
    double temp=0;
    for(int i=0;i<n;i++){
        if(TMath::Abs(x[i]-xmean)>temp)
            temp=TMath::Abs(x[i]-xmean);
    }
    return temp/xmean;
}
