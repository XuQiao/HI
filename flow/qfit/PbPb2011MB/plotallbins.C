#include "par.h"
void plotallbins(){
    const int nnonf = 2;
    int fixv2=0;
    int fixg2=1;
    double avgmult[nnonf][nbin],avgtrk[nnonf][nbin], v2[nnonf][nbin],v2err[nnonf][nbin],g2[nnonf][nbin],g2err[nnonf][nbin],v2calc[nnonf][nbin],g2calc[nnonf][nbin];
    for(int ibin=0;ibin<nbin;ibin++){
        for(int inonf=0;inonf<nnonf;inonf++){
        TFile *f = TFile::Open(Form("qfitV.root"));
        TVectorD *r;
        if(!inonf)
        r = (TVectorD*)f->Get(Form("r_%d_%d_%d",ibin,fixv2,fixg2));
        else
        r = (TVectorD*)f->Get(Form("rnonf_%d_%d_%d",ibin,fixv2,fixg2));
        g2[inonf][ibin]=(*r)[0];
        g2err[inonf][ibin]=(*r)[1];
        v2[inonf][ibin]=TMath::Abs((*r)[2]);
        v2err[inonf][ibin]=(*r)[3];
        avgmult[inonf][ibin]=(*r)[4];
        avgtrk[inonf][ibin]=(*r)[5];
        if(!inonf){
        v2calc[inonf][ibin]=(*r)[6];
        g2calc[inonf][ibin]=(*r)[7];
        }
        }
    }
    TGraphErrors *grv = new TGraphErrors(nbin,avgtrk[0],v2[0],0,v2err[0]);
    TGraphErrors *grvnonf = new TGraphErrors(nbin,avgtrk[0],v2[1],0,v2err[1]);
    TGraphErrors *grvcalc = new TGraphErrors(nbin,avgtrk[0],v2calc[0],0,v2err[0]);
    TGraphErrors *grg = new TGraphErrors(nbin,avgtrk[0],g2[0],0,g2err[0]);
    TGraphErrors *grgnonf = new TGraphErrors(nbin,avgtrk[0],g2[1],0,g2err[1]);
    TGraphErrors *grgcalc = new TGraphErrors(nbin,avgtrk[0],g2calc[0],0,g2err[0]);
    TLegend *tl = new TLegend(0.40,0.62,0.80,0.88);
    tl->SetTextSize(0.035);
    tl->SetBorderSize(0);
    tl->SetFillColor(0);

    TCanvas *c1 = new TCanvas("c1","c1",600,600);
    grv->SetTitle("");
    grv->GetXaxis()->SetTitle("Ntrkoffline");
    grv->GetYaxis()->SetTitle("v_{2}");
    grv->GetYaxis()->SetRangeUser(0.0,0.2);
    grv->SetMarkerStyle(20);
    grv->SetMarkerSize(1.4);
    grvnonf->SetMarkerStyle(24);
    grvnonf->SetMarkerColor(2);
    grvnonf->SetLineColor(2);
    grvcalc->SetMarkerStyle(29);
    grvcalc->SetMarkerColor(4);
    grvcalc->SetLineColor(4);
    grvcalc->SetMarkerSize(1.4);
    tl->AddEntry(grv,"result by fitting","lp");
    //tl->AddEntry(grvcalc,"result by calculating","lp");
    grv->Draw("AP");
    //grvnonf->Draw("Psame");
    //grvcalc->Draw("Psame");
    tl->Draw("same");

    TCanvas *c2 = new TCanvas("c2","c2",600,600);
    grg->SetTitle("");
    grg->GetXaxis()->SetTitle("Ntrkoffline");
    grg->GetYaxis()->SetTitle("g_{2}");
    grg->GetYaxis()->SetRangeUser(-1,5);
    grg->SetMarkerStyle(20);
    grg->SetMarkerSize(1.4);
    grgnonf->SetMarkerStyle(24);
    grgnonf->SetMarkerColor(2);
    grgnonf->SetLineColor(2);
    grgcalc->SetMarkerStyle(29);
    grgcalc->SetMarkerColor(4);
    grgcalc->SetLineColor(4);
    grgcalc->SetMarkerSize(1.4);
    grg->Draw("AP");
    //grgnonf->Draw("Psame");
    //grgcalc->Draw("Psame");
    tl->Draw("same");

    c1->Print("v2vsNtrk.png");
    c2->Print("g2vsNtrk.png");
}

