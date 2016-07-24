#include <iomanip>
#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"
void checkRun(){
 TCanvas *c1 = new TCanvas("c1a", "c1",0,0,600,600);
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetErrorX(0);
    c1->Range(0,0,1,1);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
    c1->SetTickx(1);
    c1->SetTicky(1);
    c1->SetLeftMargin(0.13);
    c1->SetRightMargin(0.13);
    c1->SetTopMargin(0.05);
    c1->SetBottomMargin(0.1);
    c1->SetFrameFillStyle(0);
    c1->SetFrameBorderMode(0);

    gStyle->SetOptStat(0);
    gStyle->SetPadBottomMargin(0.04);
    gStyle->SetPadTopMargin   (0.025);
    gStyle->SetPadLeftMargin  (0.04);
    gStyle->SetPadRightMargin (0.025);
    gStyle->SetPadTickX       (1);
    gStyle->SetPadTickY       (1);
ofstream fstr("RunEvent.txt");
const int N=1;
TString TrigName[N]={"Jet20"};
TH1I* hRun[N];
TGraph* g[N];
int Nrun[N];
int irun=0;
int vdm=211561, JSonstart=210676, JSonend=211631;
int runnumber[N][100];
int Nevent[N][100];
for(int i=0;i<N;i++){
f[i]=TFile::Open(Form("/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/GlobalEvent/%sEvSel_v2.root",TrigName[i].Data()));
fstr<<TrigName[i]<<":"<<endl;
hRun[i]=(TH1I*)f[i]->Get("hRun");
for(int ibin=1;ibin<hRun[i]->GetNbinsX();ibin++){
if(hRun[i]->GetBinContent(ibin)!=0){
Nrun[i]++;
fstr<<hRun[i]->GetXaxis()->GetBinLowEdge(ibin);
fstr<<'\t'<<hRun[i]->GetBinContent(ibin)<<endl;
runnumber[i][irun]=hRun[i]->GetXaxis()->GetBinLowEdge(ibin);
Nevent[i][irun]=hRun[i]->GetBinContent(ibin);
//fstr<<"to"<<hRun[i]->GetXaxis()->GetBinUpEdge(ibin);
irun++;
}
}
cout<<(double)hRun[i]->GetBinContent(hRun[i]->GetXaxis()->FindBin(vdm))/hRun[i]->Integral(hRun[i]->GetXaxis()->FindBin(JSonstart),hRun[i]->GetXaxis()->FindBin(JSonend))<<endl;
g[i]=new TGraph(Nrun[i],runnumber[i],Nevent[i]);
g[i]->GetXaxis()->SetTitle("Run number");
g[i]->GetYaxis()->SetTitle("Number of Events");
g[i]->SetMarkerColor(1);
g[i]->SetMarkerSize(1.2);
g[i]->SetMarkerStyle(20);
TLine *l1=new TLine(JSonstart,0,JSonstart,2e5);
TLine *l2=new TLine(JSonend,0,JSonend,2e5);
TLine *l3=new TLine(vdm,0,vdm,2e5);
TLine *l4=new TLine(211631,0,211631,2e5);
l1->SetLineStyle(2);
l2->SetLineStyle(2);
l3->SetLineStyle(2);
l1->SetLineColor(2);
l2->SetLineColor(2);
l3->SetLineColor(4);
g[i]->Draw("AP");
l1->Draw("same");
l2->Draw("same");
l3->Draw("same");
}

}
