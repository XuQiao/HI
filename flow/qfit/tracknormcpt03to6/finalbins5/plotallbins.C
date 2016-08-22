void plotallbins(){
    const int nbin24 = 12;
    const double avgtrkbin[nbin24]={44.36,54.33,68.63,88.39,108.2,131.3,162.1,196.6,227.5,247.2,269.2,301.2};
    const double V24[nbin24]={0.02965,0.03913,0.04832,0.04941,0.04822,0.04955,0.049,0.04805,0.04709,0.04665,0.04772,0.04797};
    const double V24err[nbin24]={0.005967,0.004262,0.001496,0.001351,0.001599,0.0003922,0.0003065,0.0002939,0.0004568,0.0008684,0.001471,0.004329};
    const int nbin22 = 18;
    const double avgtrkbin22[nbin22]={6.917,14.37,24.41,34.39,44.36,54.33,68.63,88.39,108.2,131.3,162.1,196.6,227.5,247.2,267.1,287,306.8,328.6};
    const double V22[nbin22]={0.05367,0.05428,0.05379,0.05501,0.05651,0.05814,0.06038,0.0629,0.06421,0.06633,0.06777,0.06887,0.06927,0.0697,0.06977,0.06896,0.06706,0.06621};
 
   const double V22err[nbin22]={0.002093,0.002117,0.002098,0.002145,0.002204,0.002267,0.002355,0.002453,0.002504,0.002587,0.002643,0.002686,0.00381,0.003833,0.005163,0.005103,0.004962,0.0049};
    const int ntot = 5;
    const int nnonf = 2;
    TString dir[ntot]={"M150120","M185150","M220185","M260220","M300260"};
    const int nbin[ntot] = {6,7,7,8,8};
    const int nbin_[ntot] = {0,6,13,20,28};
    const int ntotbinsf = 36;
    int fixv2=0;
    double avgmult_[nnonf][ntotbinsf], avgtrk_[nnonf][ntotbinsf],v2_[nnonf][ntotbinsf],v2err_[nnonf][ntotbinsf],g2_[nnonf][ntotbinsf],g2err_[nnonf][ntotbinsf];
    for(int ibin=0;ibin<ntot;ibin++){
        for(int i=0;i<nbin[ibin];i++){
        for(int inonf=0;inonf<nnonf;inonf++){
        TFile *f = TFile::Open(Form("%s/qfitV.root",dir[ibin].Data()));
        if(!inonf)
        TVectorD *r = (TVectorD*)f->Get(Form("r_%d_%d",i,fixv2));
        else
        TVectorD *r = (TVectorD*)f->Get(Form("rnonf_%d_%d",i,fixv2));
        g2_[inonf][nbin_[ibin]+i]=(*r)[0];
        g2err_[inonf][nbin_[ibin]+i]=(*r)[1];
        v2_[inonf][nbin_[ibin]+i]=TMath::Abs((*r)[2]);
        v2err_[inonf][nbin_[ibin]+i]=(*r)[3];
        avgmult_[inonf][nbin_[ibin]+i]=(*r)[4];
        avgtrk_[inonf][nbin_[ibin]+i]=(*r)[5];
        if(ibin==ntot-1 && i==0){
            v2_[inonf][nbin_[ibin]+i]=-999;
            v2err_[inonf][nbin_[ibin]+i]=0;
        }
        }
        }
    }

    TGraphErrors *grv24=new TGraphErrors(nbin24,avgtrkbin,V24,0,V24err);
    TGraphErrors *grv22=new TGraphErrors(nbin22,avgtrkbin22,V22,0,V22err);
    TGraphErrors *grvfinal = new TGraphErrors(ntotbinsf,avgtrk_[0],v2_[0],0,v2err_[0]);
    TGraphErrors *grvnonffinal = new TGraphErrors(ntotbinsf,avgtrk_[0],v2_[1],0,v2err_[1]);
    TGraphErrors *grgfinal = new TGraphErrors(ntotbinsf,avgtrk_[0],g2_[0],0,g2err_[0]);
    TGraphErrors *grgnonffinal = new TGraphErrors(ntotbinsf,avgtrk_[0],g2_[1],0,g2err_[1]);

    const int ntotbins = 5;
    int xbin=0;
    double avgmult[nnonf][ntotbinsf], avgtrk[nnonf][ntotbinsf],v2[nnonf][ntotbinsf],v2err[nnonf][ntotbinsf],g2[nnonf][ntotbinsf],g2err[nnonf][ntotbinsf];
    for(int ibin=0;ibin<ntotbins;ibin++){
        for(int inonf=0;inonf<nnonf;inonf++){
        TFile *f = TFile::Open(Form("../%s/qfitV.root",dir[ibin].Data()));
        if(!inonf)
        TVectorD *r = (TVectorD*)f->Get(Form("r_%d",xbin));
        else
        TVectorD *r = (TVectorD*)f->Get(Form("rnonf_%d",xbin));
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
    
    TLegend *tl = new TLegend(0.15,0.75,0.45,0.88);
    tl->SetTextSize(0.035);
    tl->SetBorderSize(0);
    tl->SetFillColor(0);

    TCanvas *c1 = new TCanvas("c1","c1",600,600);
    grv->SetTitle("");
    grv->GetXaxis()->SetTitle("Ntrkoffline");
    grv->GetYaxis()->SetTitle("v_{2}");
    grv->GetYaxis()->SetRangeUser(0.,0.12);
    grv->GetXaxis()->SetLimits(30,350);
    grv->SetMarkerStyle(20);
    grv->SetMarkerSize(1);
    grv24->SetMarkerStyle(24);
    grv24->SetMarkerColor(2);
    grv24->SetLineColor(2);
    grvfinal->SetMarkerColor(4);
    grvfinal->SetMarkerStyle(27);
    grvfinal->SetLineColor(4);
    grv22->SetMarkerStyle(29);
    grv22->SetMarkerColor(8);
    grv22->SetLineColor(8);
    tl->AddEntry(grv,"bin width ~ 40","p");
   // tl->AddEntry(grvnonf,"bin width = 20 with non-flow","p");
    tl->AddEntry(grvfinal,"bin width = 5","p");
    tl->AddEntry(grv24,"v_{2}{4}","p");
    tl->AddEntry(grv22,"v_{2}{2}","p");
    grv->Draw("AP");
    grvfinal->Draw("Psame");
    grv24->Draw("Psame");
    grv22->Draw("Psame");
    //grvnonf->Draw("Psame");
    //grvnonffinal->Draw("Psame");
    tl->Draw("same");

    TCanvas *c2 = new TCanvas("c2","c2",600,600);
    grg->SetTitle("");
    grg->GetXaxis()->SetTitle("Ntrkoffline");
    grg->GetYaxis()->SetTitle("g_{2}");
    grg->GetYaxis()->SetRangeUser(0,3);
    grg->SetMarkerStyle(20);
    grg->SetMarkerSize(1);
    grgnonf->SetMarkerStyle(24);
    grgnonf->SetMarkerColor(2);
    grgnonf->SetLineColor(2);
    grgfinal->SetMarkerColor(4);
    grgfinal->SetMarkerStyle(27);
    grgfinal->SetLineColor(4);
    grgnonffinal->SetMarkerStyle(29);
    grgnonffinal->SetMarkerColor(8);
    grgnonffinal->SetLineColor(8);
    grg->Draw("AP");
    grgfinal->Draw("Psame");
    //grgnonf->Draw("Psame");
    //grgnonffinal->Draw("Psame");
    tl->Draw("same");

    c1->Print(Form("v2vsNtrk_fixv2%d.png",fixv2));
    c2->Print(Form("g2vsNtrk_fixv2%d.png",fixv2));
}

