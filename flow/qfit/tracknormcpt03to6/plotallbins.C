void plotallbins(){
    const int ntotbins = 5;
    const int nnonf = 2;
    int xbin=0;
    TString dir[ntotbins]={"M150120","M185150","M220185","M260220","M300260"};
    double avgmult[nnonf][ntotbins],avgtrk[nnonf][ntotbins], v2[nnonf][ntotbins],v2err[nnonf][ntotbins],g2[nnonf][ntotbins],g2err[nnonf][ntotbins];
    for(int ibin=0;ibin<ntotbins;ibin++){
        for(int inonf=0;inonf<nnonf;inonf++){
        TFile *f = TFile::Open(Form("%s/qfitV.root",dir[ibin].Data()));
        TVectorD *r;
        if(!inonf)
        r = (TVectorD*)f->Get(Form("r_%d",xbin));
        else
        r = (TVectorD*)f->Get(Form("rnonf_%d",xbin));
        g2[inonf][ibin]=(*r)[0];
        g2err[inonf][ibin]=(*r)[1];
        v2[inonf][ibin]=TMath::Abs((*r)[2]);
        v2err[inonf][ibin]=(*r)[3];
        avgmult[inonf][ibin]=(*r)[4];
        avgtrk[inonf][ibin]=(*r)[5];
        }
    }
    TGraphErrors *grv = new TGraphErrors(ntotbins,avgtrk[0],v2[0],0,v2err[0]);
    TGraphErrors *grvnonf = new TGraphErrors(ntotbins,avgtrk[0],v2[1],0,v2err[1]);
    TGraphErrors *grg = new TGraphErrors(ntotbins,avgtrk[0],g2[0],0,g2err[0]);
    TGraphErrors *grgnonf = new TGraphErrors(ntotbins,avgtrk[0],g2[1],0,g2err[1]);

    TCanvas *c1 = new TCanvas("c1","c1",600,600);
    grv->SetTitle("");
    grv->GetXaxis()->SetTitle("Ntrkoffline");
    grv->GetYaxis()->SetTitle("v_{2}");
    grv->GetYaxis()->SetRangeUser(0.03,0.08);
    grv->SetMarkerStyle(20);
    grv->SetMarkerSize(1);
    grvnonf->SetMarkerStyle(24);
    grvnonf->SetMarkerColor(2);
    grvnonf->SetLineColor(2);
    grv->Draw("AP");
    grvnonf->Draw("Psame");

    TCanvas *c2 = new TCanvas("c2","c2",600,600);
    grg->SetTitle("");
    grg->GetXaxis()->SetTitle("Ntrkoffline");
    grg->GetYaxis()->SetTitle("g_{2}");
    grg->GetYaxis()->SetRangeUser(0,5);
    grg->SetMarkerStyle(20);
    grg->SetMarkerSize(1);
    grgnonf->SetMarkerStyle(24);
    grgnonf->SetMarkerColor(2);
    grgnonf->SetLineColor(2);
    grg->Draw("AP");
    grgnonf->Draw("Psame");

    c1->Print("v2vsNtrk.png");
    c2->Print("g2vsNtrk.png");
}

