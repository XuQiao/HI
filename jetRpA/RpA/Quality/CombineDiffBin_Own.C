TH1 * normalizeByBinWidth(TH1 *histo) {
    for(int i = 1; i <= histo->GetNbinsX(); i++) {
        float content = histo->GetBinContent(i);
        float error = histo->GetBinError(i);
        histo->SetBinContent(i,content/histo->GetBinWidth(i));
        histo->SetBinError(i,error/histo->GetBinWidth(i));
    }
    return histo ;
}


void CombineDiffBin_Own(){
TString effTab = "2013HistIterTrkCorrtest";
double trackcut=0;
TString coll="PP";

const int Nfile =8;
        const int pthat[]={15, 30, 50, 80,120, 170,220, 280};
        const double cs[]={5.335E-01, 3.378E-02, 3.778E-03, 4.412E-04, 6.147E-05,1.018E-05,2.477E-06,6.160E-07};
        //2013 pp @ 5.02 TeV MC sample with pp tracking
        const double Nevt[]={197175, 199375, 182175, 199375, 164900, 200975, 180200,194375};
 int N;       
    TFile * f[Nfile] ;
    TString name[Nfile] ;
    TList * list[Nfile] ;
    double weight[Nfile] ;
int filestart=0;
    const char * kDir = "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID";
    for(int ifile=filestart; ifile<Nfile; ifile++){

	name[ifile]=Form("MC%s_AkPu3PFJetGenCharge%s%.fEtaCut10_pt%d_HP04_hiforest77_hiSignal.root",coll.Data(),effTab.Data(), trackcut ,pthat[ifile]);
	  f[ifile] = TFile::Open(Form("%s/%s", kDir, name[ifile].Data()), "readonly");

      list[ifile] =(TList*)f[ifile]->GetListOfKeys();
      if(ifile==filestart){
        N = list[ifile]->GetEntries();
cout<<"N = "<<N<<endl;
      }
        if(ifile==8)weight[ifile]=cs[ifile]/Nevt[ifile];
        else weight[ifile]=(cs[ifile]-cs[ifile+1])/Nevt[ifile];

        cout <<"coll = " << coll << " ifile = " << ifile << " weight = " << weight[ifile] <<endl ;
    }


      const int Nhist = N ;
      TH1 * h[Nfile][Nhist];
      TH1 * htot[Nhist];
    for(int ifile=filestart; ifile<Nfile; ifile++){
      for(int ihist = 0 ; ihist < Nhist ; ihist++){
         h[ifile][ihist] = (TH1*)f[ifile]->Get(list[ifile]->At(ihist)->GetName()) ;
        //h[ifile][ihist]->Sumw2();
                h[ifile][ihist]->Scale(weight[ifile]);
                      }
                          }
        

 TString outname=Form("merged2CSdiff_MC%s_AkPu3PFIncJetGenCharge%s%.fEtaCut10_HP04_hiforest77_hiSignal.root",coll.Data(),effTab.Data(), trackcut);

 TFile * fmerge = new TFile(Form("%s/%s", kDir,outname.Data()), "RECREATE");   
 TCanvas  * c1 = new TCanvas("c1", "c1", 800, 800);
    c1->SetFillColor(0);
    c1->SetBorderSize(0);
    c1->SetFrameBorderMode(0);
    gStyle->SetOptStat(0);

    TLegend *t1=new TLegend(0.65,0.58,0.85,0.88);
    t1->SetFillColor(0);
    t1->SetBorderSize(0);
    t1->SetFillStyle(0);
    t1->SetTextFont(63);
    t1->SetTextSize(17);

    TH1D * hFrame = new TH1D("hFrame","frame",500,0.,500.);
  hFrame-> SetAxisRange(30,500.,"X");
    hFrame->SetAxisRange(1.e-10, 10., "Y") ;

hFrame->GetXaxis()->SetTitle("p_{T} (GeV/c)");

    const Int_t dis = 9 ;

    for(int ihist = 0 ; ihist < Nhist ; ihist++){

            for(int ifile=filestart; ifile<Nfile; ifile++){
                if(ifile==filestart)
                    htot[ihist] = (TH1*)h[ifile][ihist]->Clone(Form("%s",list[ifile]->At(ihist)->GetName()));
                else
                    htot[ihist]->Add(h[ifile][ihist]);

            }

        if(htot[ihist]->GetEntries()>0){
    htot[ihist]->Write();
        }

 //pick only one histogram for display
      if(ihist==dis){
              c1->cd() ;
              c1->SetLogy() ;
              cout <<list[filestart]->At(dis)->GetName()<<endl ;
              hFrame->DrawCopy();

	    normalizeByBinWidth(htot[dis]);
            htot[dis]->SetLineColor(1);
            htot[dis]->SetMarkerColor(1);
            htot[dis]->SetMarkerStyle(20);
            htot[dis]->SetMarkerSize(0.2);
            htot[dis]->SetTitle(Form("%s",list[filestart]->At(dis)->GetName()));
            htot[dis]->DrawCopy("same PE");
	for(int i=1; i<htot[dis]->GetNbinsX()+1;i++)
	cout<<htot[dis]->GetBinError(i)/htot[dis]->GetBinContent(i)<<endl;
            t1->AddEntry(htot[dis], "merged", "P");
            for(int ifile=filestart; ifile<Nfile; ifile++){
 
              normalizeByBinWidth(h[ifile][dis]);
                h[ifile][dis]->SetLineColor(ifile+2);
                h[ifile][dis]->SetMarkerColor(ifile+2);
                h[ifile][dis]->SetMarkerStyle(21+ifile);
           	h[ifile][dis]->SetMarkerSize(0.2);
                h[ifile][dis]->DrawCopy("same PE");
                t1->AddEntry(h[ifile][dis], Form("\#hat{p}_{T}=%dGeV/c",pthat[ifile]), "P");
            }

            t1->Draw("same");
            c1->Update();

        }
    }
   fmerge->Close() ;

}
