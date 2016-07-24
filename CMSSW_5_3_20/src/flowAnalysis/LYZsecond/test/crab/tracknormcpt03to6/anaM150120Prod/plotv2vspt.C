#include "par.h"
void plotv2vspt(){

int ibin=0;
const int npt24 = 9;
const double pt24[npt24] = {0.4,0.7,1.2,1.7,2.2,2.7,3.4,4.4,5.4};
const double v24[nbin][npt24]={0.02809,0.04599,0.06869,0.0844,0.09715,0.09357,0.07924,0.04701,0.0162};
const double v24err[nbin][npt24]={0.00077,0.00057,0.00085,0.0014,0.0022,0.0034,0.0041,0.0072,0.011};
TFile *f = TFile::Open("mergedv_Prod.root");
TVectorD *vecDpt = (TVectorD*)f->Get(Form("D_%d/avgpt",ibin));
TVectorD *vecDv2 = (TVectorD*)f->Get(Form("D_%d/vmean",ibin));
TVectorD *vecDv2err = (TVectorD*)f->Get(Form("D_%d/deltavmean",ibin));
double *pt = vecDpt->GetMatrixArray();
double *v2 = vecDv2->GetMatrixArray();
double *v2err = vecDv2err->GetMatrixArray();
double v2_[npt24], v2err_[npt24];
for(int i=0;i<npt24;i++){
	v2_[i]=v2[i+1];
	v2err_[i]=v2err[i+1];
}
TCanvas *c1 = new TCanvas;
TGraphErrors *gr=new TGraphErrors(npt24,pt24,v2_,0,v2err_);
TGraphErrors *grv24=new TGraphErrors(npt24,pt24,v24[ibin],0,v24err[ibin]);
TLatex *tl = new TLatex(0.2,0.2,Form("%d < Ntrkoffline < %d",trkbin[ibin],trkbin[ibin+1]));
tl->SetNDC();
gr->SetTitle("v_{2} vs momentum");
gr->GetXaxis()->SetTitle("p_{T} (GeV/c)");
gr->GetYaxis()->SetTitle("v_{2}");
gr->SetMarkerSize(1);
grv24->SetMarkerSize(1);
gr->SetMarkerColor(1);
grv24->SetMarkerColor(2);
gr->SetMarkerStyle(20);
grv24->SetMarkerStyle(20);
TLegend *leg = new TLegend(0.2,0.7,0.4,0.8);
leg->SetTextSize(0.05);
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->AddEntry(gr,"LYZ method","P");
leg->AddEntry(grv24,"4-particle cumulant","P");
gr->Draw("AP");
grv24->Draw("Psame");
tl->Draw("same");
leg->Draw("same");
c1->Print("v2vspt_M185150.png");
tl->Draw();
f->Close();

}

