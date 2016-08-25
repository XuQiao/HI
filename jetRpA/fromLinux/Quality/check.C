void check(){
//cout<<fixed<<setprecision(3);
double binbound_pt[]={0,5,10,15,20,30,45,60,75,90,105,120,140,160,180,200,220,260,300,400,600,1000};
int nbinspt=sizeof(binbound_pt)/sizeof(double)-1;
TFile *f=TFile::Open("~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/mergedCSdiff_MCPP_AllAkPu3PFIncJetNoResidual2013HistIterTrkCorrtest0EtaCut10_pthatcutnewbin_HP04_hiforest77_hiSignal.root");
//TFile *f=TFile::Open("~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/mergedCSdiff_MCPPb_AllAkPu3PFIncJetNoResidual2013HistIterTrkCorrtest0EtaCut10_pthatcutnewbin_HP04_prod16_v77_merged_forest_0.root");
TH2F* histo2F=(TH2F*)f->Get("jetptSumSumrawpt_fake_0-100%");
TH2F *h =(TH2F*)histo2F->Clone("h"); 
h->GetXaxis()->SetRangeUser(binbound_pt[5],binbound_pt[nbinspt]);
/*TH2F *h=new TH2F("oldrebin",histo2F->GetTitle(),nbinspt,binbound_pt,500,0.,5.);
	h->Sumw2();
   TAxis *xaxis = histo2F->GetXaxis();
   TAxis *yaxis = histo2F->GetYaxis();
   for(int j=1;j<=yaxis->GetNbins();j++) {
      for(int i=1;i<=xaxis->GetNbins();i++) {
         h->Fill(xaxis->GetBinCenter(i),yaxis->GetBinCenter(j),histo2F->GetBinContent(i,j));
      }
   }
*/
h->Draw("colz");
cout<<h->Integral(0,-1,0,-1)<<endl;
TProfile* hnewProf=(TProfile*)h->ProfileX("hnewProf");
hnewProf->SetMarkerStyle(20);
hnewProf->SetMarkerColor(1);
//hnewProf->Draw();
//const int nbinspt=histo2F->GetNbinsX();
TProfile* hProf=histo2F->ProfileX("hProf");
TProfile* rehProf=hProf->Rebin(nbinspt,"rehProf",binbound_pt);
//rehProf->Draw("same");

TH1D* histo1D;
//histo1D=h->ProjectionY("histo1D",0,0,"e");
//cout<<histo1D->Integral(0,-1)<<'\t';

//for(int i=0;i<nbinspt;i++){
i=7;
histo1D=h->ProjectionY("histo1D",h->GetXaxis()->FindBin(binbound_pt[i]+1e-4),h->GetXaxis()->FindBin(binbound_pt[i+1]-1e-4),"e");
//histo1D=h->ProjectionY("histo1D",-1,-1);
histo1D->Draw();
if(histo1D->GetEntries()!=0){
cout<<binbound_pt[i]<<"-"<<binbound_pt[i+1]<<'\t';
cout<<histo1D->Integral(0,500)<<'\t';
//cout<<hnewProf->GetBinContent(i)<<'\t';
//cout<<hnewProf->GetBinError(i)<<'\t';
//cout<<hnewProf->GetBinError(i)/hnewProf->GetBinContent(i)<<'\t'<<"| ";
/*cout<<histo1D->GetMean()<<'\t';
cout<<histo1D->GetRMS()<<'\t';
cout<<histo1D->GetRMS()/TMath::Sqrt(histo1D->GetEntries())<<'\t';
cout<<histo1D->GetRMS()/TMath::Sqrt(histo1D->GetEntries())/histo1D->GetMean()<<endl;*/
//cout<<endl;
//}
}

}
