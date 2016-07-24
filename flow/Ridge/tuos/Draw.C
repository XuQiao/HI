void Draw(){
    gStyle->SetOptStat(kFALSE);
    TFile *f = TFile::Open(Form("correlationSTEG.root"));
    TH2F* hr = (TH2F*)f->Get(Form("correlation"));
    hr->SetTitle("");
    hr->GetXaxis()->SetRangeUser(-4,4);
    hr->GetYaxis()->SetRangeUser(-1.5,4);
 //   hr->GetZaxis()->SetRangeUser(0,3e6);
    hr->GetXaxis()->SetTitle("#Delta#eta");
    hr->GetXaxis()->CenterTitle();
    hr->GetYaxis()->SetTitle("#Delta#phi");
    hr->GetYaxis()->CenterTitle();
    hr->GetZaxis()->SetTitle("#frac{1}{N_{trig}}#frac{d^{2}N^{pair}}{d#Delta#eta d#Delta#phi}");
    TCanvas *c1 = new TCanvas();
    hr->Draw("surf1");
    TLatex l;
    l.SetNDC();
    l.SetTextSize(0.04);
    c1->Print("ridge.png");
}

