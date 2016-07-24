void fitandcorr(char origfile[], char ratiofile[], char fitfunction[], char orighistoname[],char ratiohistoname[], char genname[]){

char buf[1024],gentrname[1024];
sprintf(gentrname,"%s_tr",orighistoname);
TFile *fRatio = new TFile(ratiofile,"update");
TFile *forig = new TFile(origfile,"update");
TH1F *histogentr = (TH1F*)fRatio->Get(gentrname);
TH1F *histoRatio = (TH1F*)fRatio->Get(ratiohistoname);
TH1F *histogentrtest = (TH1F*)histogentr->Clone();
TH1F *histoRatiofit = (TH1F*)histoRatio->Clone();
TH1F *histoorig = (TH1F*)forig->Get(orighistoname);

TF1 *fitEff = new TF1("fitEff",fitfunction,0,histoRatio->GetXaxis()->GetXmax());
histoRatio->Fit("fitEff","R");

TCanvas *c1=new TCanvas("c1");
c1->cd();
gStyle->SetOptFit(100);
gStyle->SetOptStat(kFALSE);
histoRatio->Draw();
c1->Update();
TPaveStats *st = (TPaveStats*)c1->GetPrimitive("stats");
st->SetX1NDC(0.68);
st->SetX2NDC(0.88);
st->SetY1NDC(0.62);
st->SetY2NDC(0.72);
st->SetLineColor(0);
c1->Update();
c1->Modified();
Double_t binsize=histoRatiofit->GetBinWidth(1);

for(int i=0;i<histoRatiofit->GetNbinsX();i++)
histoRatiofit->SetBinContent(i,fitEff->Integral(i*binsize,(i+1)*binsize)/binsize);

fRatio->cd();
sprintf(buf,"%s_fit",ratiohistoname);
histoRatiofit->Write(buf,TObject::kOverwrite);
sprintf(buf,"%s_corr",orighistoname);
histogentrtest->Divide(histoRatiofit);
histogentrtest->Write(buf,TObject::kOverwrite);

forig->cd();
sprintf(buf,"%s_%s_corr",orighistoname,genname);
histoorig->Divide(histoRatiofit);
histoorig->Write(buf,TObject::kOverwrite);

delete c1;

fRatio->Close();
forig->Close();

}

void Eff_Cor(){
fitandcorr("pPbHistPlusEta4.root","pPbHist_Hijing_TrandEs.root","[0]-[1]*exp([2]*x)-[3]*exp(-[4]*x)","hHFEnergy4","hEff4","Hijing");
fitandcorr("pPbHistPlusEta4.root","pPbHist_Hijing_TrandEs.root","[0]-[1]*exp([2]*x)-[3]*exp(-[4]*x)","hHFEnergyPlus4","hEffPlus4","Hijing");
fitandcorr("pPbHistPlusEta4.root","pPbHist_Hijing_TrandEs.root","[0]-[1]*exp([2]*x)-[3]*x","hNtrack","hEffNtrack","Hijing");
fitandcorr("pPbHistPlusEta4.root","pPbHist_Epos_TrandEs.root","[0]-[1]*exp([2]*x)-[3]*exp(-[4]*x)","hHFEnergy4","hEff4","Epos");
fitandcorr("pPbHistPlusEta4.root","pPbHist_Epos_TrandEs.root","[0]-[1]*exp([2]*x)-[3]*exp(-[4]*x)","hHFEnergyPlus4","hEffPlus4","Epos");
fitandcorr("pPbHistPlusEta4.root","pPbHist_Epos_TrandEs.root","[0]-[1]*exp([2]*x)-[3]*x","hNtrack","hEffNtrack","Epos");
}
