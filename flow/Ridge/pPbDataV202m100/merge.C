#include "par.h"
void merge(){
	TVectorD Nevent;	Nevent.ResizeTo(nbin);  Nevent.Zero();
     //   TVectorD totmultall;	totmultall.ResizeTo(nbin);      totmultall.Zero();
     //   TVectorD avgmultall;	avgmultall.ResizeTo(nbin);      avgmultall.Zero();
        TVectorD tottrk;	tottrk.ResizeTo(nbin);      tottrk.Zero();
//	TVectorD totptall;      totptall.ResizeTo(nbin);    totptall.Zero();
//	TVectorD totetaall;      totetaall.ResizeTo(nbin);    totetaall.Zero();
        TVectorD avgtrk;	avgtrk.ResizeTo(nbin);      avgtrk.Zero();
//	TVectorD avgmult;       avgmult.ResizeTo(nbin);
            
        TH2F* s[nbin];
        TH2F* b[nbin];
        for(int ibin=0;ibin<nbin;ibin++){
            s[ibin] = new TH2F(Form("s_%d",ibin),Form("signal",ibin),detastep,detamin,detamax,dphistep,dphimin,dphimax);
            s[ibin]->Sumw2();
            b[ibin]  = new TH2F(Form("b_%d",ibin), "background",detastep,detamin,detamax,dphistep,dphimin,dphimax);
            b[ibin]->Sumw2();
        }
        TFile *fout = new TFile(Form("Anav3_merged.root"),"Recreate");
        TFile *f[nFileAll];
        for(int ifile=; ifile<63; ifile++){
                f[ifile] = TFile::Open(Form("%s/Anav3_%d.root",outdir.Data(),ifile));
                if(!f[ifile]) continue;
		TVectorD* Nevent_t =  (TVectorD*)f[ifile]->Get(Form("Nevent"));
	//	TVectorD* totmultall_t =  (TVectorD*)f[ifile]->Get(Form("totmultall"));
		TVectorD* tottrk_t =  (TVectorD*)f[ifile]->Get(Form("tottrk"));
	//	TVectorD* totptall_t =  (TVectorD*)f[ifile]->Get(Form("totptall"));
//		TVectorD* totetaall_t =  (TVectorD*)f[ifile]->Get(Form("totetaall"));
		for(int ibin=0;ibin<nbin;ibin++){
			//	totptall[ibin] += (*totptall_t)[ibin];
			//	totetaall[ibin] += (*totetaall_t)[ibin];
			        Nevent[ibin] += (*Nevent_t)[ibin];
			  //      totmultall[ibin] += (*totmultall_t)[ibin];	
			        tottrk[ibin] += (*tottrk_t)[ibin];	
                                TH1D* s_t = (TH1D*)f[ifile]->Get(Form("signal_%d",ibin));
                                TH1D* b_t = (TH1D*)f[ifile]->Get(Form("background_%d",ibin));
                                s[ibin]->Add(s_t);
                                b[ibin]->Add(b_t);
		}
		f[ifile]->Close();
        }
    for(int ibin=0;ibin<nbin;ibin++){
    avgtrk[ibin] = tottrk[ibin]/Nevent[ibin];
    TH1D* hr = (TH1D*)s[ibin]->Clone("hr");
    hr->Divide(b[ibin]);
    hr->Scale(b[ibin]->GetBinContent(b[ibin]->FindBin(0,0)));

    fout->cd();
    Nevent.Write("Nevent");
    avgtrk.Write("avgtrk");
    TDirectory *dir0 = (TDirectory*)fout->mkdir(Form("D_%d",ibin));
    dir0->cd();
    s[ibin]->Write("s");
    b[ibin]->Write("b");
    hr->Write();
    }
}
