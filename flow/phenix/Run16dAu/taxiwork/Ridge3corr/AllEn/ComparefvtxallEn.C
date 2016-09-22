void ComparefvtxallEn(){
    gStyle->SetOptStat(kFALSE);
    TString en[4] = {"200GeV","62GeV","39GeV","20GeV"};
    TFile *f[4];
    TH1D* hphi[4];
    f[0] = TFile::Open("../../../work/200GeV/output_EPfvtxtrk.root","ReadOnly");
    f[1] = TFile::Open("../../../work/62GeV/output_EPfvtxtrk.root","ReadOnly");
    f[2] = TFile::Open("../../../work/39GeV/output_EPfvtxtrk.root","ReadOnly");
    f[3] = TFile::Open("../../../work/20GeV/output_EPfvtxtrk.root","ReadOnly");

    TCanvas *c1 = new TCanvas();
    c1->cd();
    TH1F* h = new TH1F("","",400,-4,4);
    h->GetYaxis()->SetRangeUser(0,0.5);
    h->GetXaxis()->SetTitle("#phi");
    //h->GetYaxis()->SetTitle("dN_{trk}/d#phi");
    h->GetYaxis()->SetTitle("Norlimalized per event");
    h->Draw();
    TLegend *leg = new TLegend(0.6,0.6,0.8,0.8);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
for(int i=0;i<4;i++){
    TH2F* h1 = (TH2F*)f[i]->Get("hfvtxtrk");
    TH1F* h2 = (TH1F*)f[i]->Get("hnfvtxtrk");
    int nevent = h2->Integral();
    hphi[i] = h1->ProjectionY(Form("hphi_%d",i),50,100);
    cout<<h1->Integral()<<" "<<h2->Integral()<<endl;
    hphi[i]->Scale(1./nevent*hphi[i]->GetBinWidth(1));
    hphi[i]->SetMarkerSize(1.3);
    hphi[i]->SetMarkerStyle(20+i);
    hphi[i]->SetMarkerColor(i+1);
    hphi[i]->SetLineColor(i+1);
    hphi[i]->DrawCopy("Psame");
    leg->AddEntry(hphi[i],en[i]);
}
    leg->Draw("same");
    c1->Print("FVTXtrknphi.png");
}

