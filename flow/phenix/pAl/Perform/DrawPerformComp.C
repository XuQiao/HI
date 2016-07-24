#include "/phenix/u/xuq/util/SimplifyLife.C"
void DrawPerformComp(){
    gStyle->SetErrorX(0);
    gStyle->SetOptStat(0);
    TFile *fmb = new TFile("mergedFull.root","ReadOnly");
    TFile *fhm = new TFile("merged_FVtx.root","ReadOnly");
    const int ncav = 2;
    TCanvas *c1[ncav];
    for(int i=0;i<ncav;i++){
        c1[i] = new TCanvas();
    }
    TH2F* hpc1hitsbbc_mb = (TH2F*)fmb->Get("hbbcnvtx_0");
    TH2F* hpc1hitsbbc_hm = (TH2F*)fhm->Get("hbbcnvtx_0");
    TH1F* hbbc_mb = (TH1F*)hpc1hitsbbc_mb->ProjectionX("hbbc_mb",0,-1);
    TH1F* hbbc_hm = (TH1F*)hpc1hitsbbc_hm->ProjectionX("hbbc_hm",0,-1);
    hbbc_mb->GetXaxis()->SetLimits(0,200);
    //hbbc_mb->Rebin(5);
    //hbbc_hm->Rebin(5);
    //hbbc_mb->Scale(1./hbbc_mb->Integral());
    //hbbc_hm->Scale(1./hbbc_hm->Integral());
    c1[0]->cd();
    c1[0]->SetLogy();
    hbbc_mb->Draw();
    SetTitle(*hbbc_mb,"bbc charge sum","# of events","");
    //SetTitle(*hbbc_mb,"nvtx layer 1","# of events","");
    SetXRange(*hbbc_mb,0,100);
    SetYRange(*hbbc_mb,1,1e9);
    SetStyle(*hbbc_mb,0.6,1,20,0,0);
    SetStyle(*hbbc_hm,0.6,2,20,0,0);
    TLegend *leg = new TLegend(0.6,0.64,0.8,0.8);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetTextSize(0.045);
    leg->AddEntry(hbbc_mb,"pp minbias","pl");
   // leg->AddEntry(hbbc_mb_scale,"pp minbias * 50","pl");
    leg->AddEntry(hbbc_hm,"pp high-mult","pl");
    leg->Draw("same");
    hbbc_hm->Draw("same");
    c1[0]->Print("fig/hbbc_comp.png");

    c1[1]->cd();
    TH1F* hbbc_ratio = (TH1F*)hbbc_hm->Clone("hbbc_ratio");
    hbbc_ratio->Divide(hbbc_mb);
    SetTitle(*hbbc_ratio,"bbc charge sum","ratio hm/mb","");
    //SetTitle(*hbbc_ratio,"nvtx layer 1","# of events","");
    SetXRange(*hbbc_ratio,0,100);
    SetYRange(*hbbc_ratio,0,3);
    hbbc_ratio->Draw();
    c1[1]->Print("fig/hbbc_ratio.png");
/*
TH1F* hnvtx = (TH1F*)hbbcnvtx[0]->ProjectionY("hnvtx",0,-1);
TH1F* hnvtx_hm = (TH1F*)hbbcnvtx_hm[0]->ProjectionY("hnvtx_hm",0,-1);
TH1F* hbbc = (TH1F*)hbbcnvtx[0]->ProjectionX("hbbc",0,-1);
TH1F* hbbc_hm = (TH1F*)hbbcnvtx_hm[0]->ProjectionX("hbbc_hm",0,-1);
c1[22]->cd();
c1[22]->SetLogy();
hnvtx->Rebin(4);
hnvtx_hm->Rebin(4);
hnvtx->Scale(1./hnvtx->Integral());
hnvtx_hm->Scale(1./hnvtx_hm->Integral());
SetRange(*hnvtx,0,1e-11,200,10);
SetTitle(*hnvtx,"#hits in cluster layer 1","normalized","");
SetStyle(*hnvtx,1.2,1,20,0,0);
SetStyle(*hnvtx_hm,1.2,2,24,0,0);
hnvtx->Draw("P");
hnvtx_hm->Draw("Psame");
c1[22]->Print("fig/hnvtx.png");
c1[23]->cd();
c1[23]->SetLogy();
SetTitle(*hbbc,"bbc charge sum","normalized","");
hbbc->Rebin(5);
hbbc_hm->Rebin(5);
hbbc->Scale(1./hbbc->Integral());
hbbc_hm->Scale(1./hbbc_hm->Integral());
SetRange(*hbbc,0,1e-11,200,10);
SetStyle(*hbbc,1.2,1,20,0,0);
SetStyle(*hbbc_hm,1.2,2,24,0,0);
hbbc->Draw("P");
hbbc_hm->Draw("Psame");
c1[23]->Print("fig/hbbc.png");
*/


}
