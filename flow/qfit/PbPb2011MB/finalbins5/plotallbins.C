#include "par.h"
void plotallbins(){
    
    const int nbin24 = 15;
    const double avgtrkbin[nbin24]={34.34,44.37,54.39,69.11,89.21,109.3,134.1,166.5,201.6,229.4,249.4,269.4,289.4,309.4,334.3};
    const double V24[nbin24]={0.0398,0.05425,0.05399,0.05996,0.0641,0.06684,0.07102,0.07426,0.07692,0.07842,0.08031,0.08057,0.08251,0.08261,0.08355};
    const double V24err[nbin24]={0.01727,0.005501,0.004521,0.001912,0.001923,0.002005,0.002131,0.002228,0.002308,0.002353,0.002409,0.002417,0.002475,0.002478,0.002507};
    const int nbin22 = 18;
    const double avgtrkbin22[nbin22]={6.967,14.22,24.3,34.34,44.37,54.39,69.1,89.21,109.3,134.1,166.5,201.6,229.4,249.4,269.4,289.4,309.4,334.3};
    const double V22[nbin22]={0.04862,0.05233,0.05476,0.05762,0.06139,0.06476,0.07128,0.07666,0.08179,0.08609,0.09026,0.09364,0.09594,0.09672,0.09802,0.09909,0.0999,0.1005};
    const double V22err[nbin22]={0.004403,0.001615,0.001643,0.001729,0.001842,0.001943,0.002138,0.0023,0.002454,0.002583,0.002708,0.002809,0.002978,0.002902,0.002941,0.002973,0.002997,0.003015};

    const int nnonf = 2;
    int fixv2=0;
    int fixg2=1;
    double avgmult[nnonf][nbin],avgtrk[nnonf][nbin], v2[nnonf][nbin],v2err[nnonf][nbin],g2[nnonf][nbin],g2err[nnonf][nbin];
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
        }
    }
    TGraphErrors *grvfinal = new TGraphErrors(nbin,avgtrk[0],v2[0],0,v2err[0]);
    TGraphErrors *grvnonffinal = new TGraphErrors(nbin,avgtrk[0],v2[1],0,v2err[1]);
    TGraphErrors *grgfinal = new TGraphErrors(nbin,avgtrk[0],g2[0],0,g2err[0]);
    TGraphErrors *grgnonfinal = new TGraphErrors(nbin,avgtrk[0],g2[1],0,g2err[1]);

    double avgmult_[nnonf][nbin],avgtrk_[nnonf][nbin], v2_[nnonf][nbin],v2err_[nnonf][nbin],g2_[nnonf][nbin],g2err_[nnonf][nbin];
    const int nbinsf = 17;
    for(int ibin=0;ibin<nbinsf;ibin++){
        for(int inonf=0;inonf<nnonf;inonf++){
        TFile *f = TFile::Open(Form("../qfitV.root"));
        TVectorD *r;
        if(!inonf)
        r = (TVectorD*)f->Get(Form("r_%d_%d_%d",ibin,fixv2,fixg2));
        else
        r = (TVectorD*)f->Get(Form("rnonf_%d_%d_%d",ibin,fixv2,fixg2));
        g2_[inonf][ibin]=(*r)[0];
        g2err_[inonf][ibin]=(*r)[1];
        v2_[inonf][ibin]=TMath::Abs((*r)[2]);
        v2err_[inonf][ibin]=(*r)[3];
        avgmult_[inonf][ibin]=(*r)[4];
        avgtrk_[inonf][ibin]=(*r)[5];
        }
    }
    TGraphErrors *grv24=new TGraphErrors(nbin24,avgtrkbin,V24,0,V24err);
    TGraphErrors *grv22=new TGraphErrors(nbin22,avgtrkbin22,V22,0,V22err);
    TGraphErrors *grv = new TGraphErrors(nbinsf,avgtrk_[0],v2_[0],0,v2err_[0]);
    TGraphErrors *grvnonf = new TGraphErrors(nbinsf,avgtrk_[0],v2_[1],0,v2err_[1]);
    TGraphErrors *grg = new TGraphErrors(nbinsf,avgtrk_[0],g2_[0],0,g2err_[0]);
    TGraphErrors *grgnonf = new TGraphErrors(nbinsf,avgtrk_[0],g2_[1],0,g2err_[1]);

    TLegend *tl = new TLegend(0.15,0.75,0.45,0.88);
    tl->SetTextSize(0.035);
    tl->SetBorderSize(0);
    tl->SetFillColor(0);
    TLegend *tl1 = new TLegend(*tl);

    TCanvas *c1 = new TCanvas("c1","c1",600,600);
    grvfinal->SetTitle("");
    grvfinal->GetXaxis()->SetTitle("Ntrkoffline");
    grvfinal->GetYaxis()->SetTitle("v_{2}");
    grvfinal->GetYaxis()->SetRangeUser(0,0.15);
    grvfinal->SetMarkerStyle(20);
    grvfinal->SetMarkerSize(1);
//    grvnonf->SetMarkerStyle(24);
//    grvnonf->SetMarkerColor(2);
//    grvnonf->SetLineColor(2);
    grv->SetMarkerStyle(27);
    grv->SetMarkerColor(2);
    grv->SetLineColor(2);
    grv24->SetMarkerStyle(24);
    grv24->SetMarkerColor(2);
    grv24->SetLineColor(2);
    grv22->SetMarkerStyle(29);
    grv22->SetMarkerColor(8);
    grv22->SetLineColor(8);
    tl->AddEntry(grv,"bin width coarser","p");
    tl->AddEntry(grvfinal,"bin width ~ 5","p");
    tl->AddEntry(grv24,"v_{2}{4}","p");
    tl->AddEntry(grv22,"v_{2}{2}","p");
    grvfinal->Draw("AP");
    grv->Draw("Psame");
    grv24->Draw("Psame");
    grv22->Draw("Psame");
    tl->Draw("same");

    TCanvas *c2 = new TCanvas("c2","c2",600,600);
    grgfinal->SetTitle("");
    grgfinal->GetXaxis()->SetTitle("Ntrkoffline");
    grgfinal->GetYaxis()->SetTitle("g_{2}");
    grgfinal->GetYaxis()->SetRangeUser(-1,5);
    grgfinal->SetMarkerStyle(20);
    grgfinal->SetMarkerSize(1);
//    grgnonf->SetMarkerStyle(24);
//    grgnonf->SetMarkerColor(2);
//    grgnonf->SetLineColor(2);
    grg->SetMarkerStyle(27);
    grg->SetMarkerColor(2);
    grg->SetLineColor(2);
    grgfinal->Draw("AP");
    grg->Draw("Psame");
    tl1->AddEntry(grv,"bin width coarser","p");
    tl1->AddEntry(grvfinal,"bin width ~ 5","p");
    tl1->Draw("same");

    c1->Print("v2vsNtrk.png");
    c2->Print("g2vsNtrk.png");
}

