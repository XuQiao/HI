#include "Dauwave.h"
void plotDaufunc(){
    const int N=8;
    Dauwave dau(N);
    const int np=1000;
    int m=100;
    int n = (int)np/m;
    double x[np],yphi[np],ypsi[np];
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            x[j+i*n] = i*TMath::Power(2,-j);
            yphi[j+i*n] = dau.phi(i,j);
            ypsi[j+i*n] = dau.psi(i,j);
        }
    }
    TGraph *g1 = new TGraph(np,x,yphi);
    TGraph *g2 = new TGraph(np,x,ypsi);
    TLatex T;
    T.SetTextSize(0.05);
    T.SetNDC();
    TCanvas *c1 = new TCanvas();
    c1->cd();
    g1->SetTitle("");
    g1->GetXaxis()->SetRangeUser(-0.1,2*N-1);
    g1->GetYaxis()->SetRangeUser(TMath::Min(TMath::MinElement(np,g1->GetY()),TMath::MinElement(np,g2->GetY()))*1.1,TMath::Max(TMath::MaxElement(np,g1->GetY()),TMath::MaxElement(np,g2->GetY()))*1.1);
    g1->GetXaxis()->SetTitle("x");
//    g1->GetYaxis()->SetTitle("#varphi(x)");
    g1->SetMarkerSize(0.8);
    g1->SetMarkerStyle(20);
    g2->SetMarkerSize(0.8);
    g2->SetMarkerStyle(24);
    g2->SetMarkerColor(2);
    g1->Draw("AP");
    g2->Draw("Psame");
    TLine *l = new TLine(-0.1,0,2*N-1,0);
    l->SetLineStyle(2);
    l->Draw("same");
    T.DrawLatex(0.7,0.2,Form("order of %d",N));
    TLegend *leg = new TLegend(0.6,0.7,0.7,0.8);
    leg->SetTextSize(0.045);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(g1,"#varphi(x)","lp");
    leg->AddEntry(g2,"#psi(x)","lp");
    leg->Draw("same");
/*    TCanvas *c2 = new TCanvas();
    c2->cd();
    g2->GetXaxis()->SetRangeUser(-0.1,2*N-1);
    g2->GetXaxis()->SetTitle("x");
    g2->GetYaxis()->SetTitle("#psi(x)");
    g2->SetMarkerSize(0.8);
    g2->SetMarkerStyle(20);
    g2->Draw("AP");
    T.DrawLatex(0.7,0.2,Form("order of %d",N));*/
    c1->Print(Form("waveletfunc_N%d.png",N));
//    c2->Print(Form("wavelet_N%d.png",N));
}
