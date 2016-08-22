#include "par.h"
void Draw(){
    gStyle->SetOptStat(kFALSE);
    TFile *f = TFile::Open(Form("Ana_merged.root"));
    int xbin=0;
    TH2F* hr = (TH2F*)f->Get(Form("D_%d/hr",xbin));
    TVectorD* vavgtrk = (TVectorD*)f->Get(Form("avgtrk"));
    hr->SetTitle("");
    hr->GetXaxis()->SetRangeUser(-4,4);
    hr->GetYaxis()->SetRangeUser(-1,4.5);
//    hr->GetZaxis()->SetRangeUser();
    hr->GetXaxis()->SetTitle("#Delta#eta");
    hr->GetXaxis()->CenterTitle();
    hr->GetYaxis()->SetTitle("#Delta#phi");
    hr->GetYaxis()->CenterTitle();
    hr->GetZaxis()->SetTitle("#frac{1}{N_{trig}}#frac{d^{2}N^{pair}}{d#Delta#eta d#Delta#phi}");
    hr->Draw("surf1");
    TLatex l;
    l.SetNDC();
    l.SetTextSize(0.04);
    l.DrawLatex(0.1,0.9,Form("%d < N_{trk}^{offline} <%d",trkbin[xbin],trkbin[xbin+1]));
    l.DrawLatex(0.5,0.9,Form("Average N_{trk}^{offline} = %.1f",(*vavgtrk)[xbin]));
    l.DrawLatex(0.1,0.85,Form("%.1f < p_{T}^{trig} < %.1f",ptmin,ptmax));
    l.DrawLatex(0.5,0.85,Form("%.1f < p_{T}^{ass} < %.1f",ptmin,ptmax));
}

