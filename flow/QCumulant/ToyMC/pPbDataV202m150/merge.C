#include <complex>
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
        TVectorD qc22, qc24, qc26, qc28;
        qc22.ResizeTo(nbin); qc22.Zero();
        qc24.ResizeTo(nbin); qc24.Zero();
        qc26.ResizeTo(nbin); qc26.Zero();
        qc28.ResizeTo(nbin); qc28.Zero();
        TVectorD c22, c24, c26, c28;
        c22.ResizeTo(nbin); c22.Zero();
        c24.ResizeTo(nbin); c24.Zero();
        c26.ResizeTo(nbin); c26.Zero();
        c28.ResizeTo(nbin); c28.Zero();
        TVectorD v22, v24, v26, v28;
        v22.ResizeTo(nbin); v22.Zero();
        v24.ResizeTo(nbin); v24.Zero();
        v26.ResizeTo(nbin); v26.Zero();
        v28.ResizeTo(nbin); v28.Zero();
            
        TH1* sumreals[nbin];
        TH1* sumimags[nbin];
        TH1* weights[nbin];
        for(int ibin=0;ibin<nbin;ibin++){
            sumreals[ibin]  = new TH1D(Form("sumreals",ibin), "Re(C{n})", maxN-2+1, 2+.5, maxN+1+.5);
            sumimags[ibin]  = static_cast<TH1*>(sumreals[ibin]->Clone(Form("sumimags",ibin)));
            weights[ibin]  = static_cast<TH1*>(sumreals[ibin]->Clone(Form("weights",ibin)));
        }
        TFile *fout = new TFile(Form("%s_merged.root",method.Data()),"Recreate");
        TFile *f[nFileAll];
        for(int ifile=0; ifile<nFileAll; ifile++){
                f[ifile] = TFile::Open(Form("%s/%s_%d.root",outdir.Data(),method.Data(),ifile));
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
                                TH1D* sumreals_t = (TH1D*)f[ifile]->Get(Form("sumreals_%d",ibin));
                                TH1D* sumimags_t = (TH1D*)f[ifile]->Get(Form("sumimags_%d",ibin));
                                TH1D* weights_t = (TH1D*)f[ifile]->Get(Form("weights_%d",ibin));
                                sumreals[ibin]->Add(sumreals_t);
                                sumimags[ibin]->Add(sumimags_t);
                                weights[ibin]->Add(weights_t);
		}
		f[ifile]->Close();
        }
    cout<<"v2{gen} = "<<inV2<<endl;
    for(int ibin=0;ibin<nbin;ibin++){
    avgtrk[ibin] = tottrk[ibin]/Nevent[ibin];
    TH1D* hQCr = (TH1D*)sumreals[ibin]->Clone("hQCr");
    TH1D* hQCi = (TH1D*)sumimags[ibin]->Clone("hQCi");
    hQCr->Divide(weights[ibin]);
    hQCi->Divide(weights[ibin]);

    complex<double> qc22com(hQCr->GetBinContent(1),hQCi->GetBinContent(1));
    complex<double> qc24com(hQCr->GetBinContent(3),hQCi->GetBinContent(3));
    complex<double> qc26com(hQCr->GetBinContent(5),hQCi->GetBinContent(5));
    complex<double> qc28com(hQCr->GetBinContent(7),hQCi->GetBinContent(7));
    qc22[ibin] = abs(qc22com);
    qc24[ibin] = abs(qc24com);
    qc26[ibin] = abs(qc26com);
    qc28[ibin] = abs(qc28com);
    c22[ibin] = qc22[ibin];
    c24[ibin] = qc24[ibin] - 2*qc22[ibin]*qc22[ibin];
    c26[ibin] = qc26[ibin] - 9*qc22[ibin]*qc24[ibin] + 12*qc22[ibin]*qc22[ibin]*qc22[ibin];
    c28[ibin] = qc28[ibin] - 16*qc26[ibin]*qc22[ibin] - 18*qc24[ibin]*qc24[ibin] + 144*qc24[ibin]*qc22[ibin]*qc22[ibin] - 144*qc22[ibin]*qc22[ibin]*qc22[ibin]*qc22[ibin];
    v22[ibin] = sqrt(c22[ibin]);
    v24[ibin] = TMath::Power(-c24[ibin],1./4);
    v26[ibin] = TMath::Power(1./4*c26[ibin],1./6);
    v28[ibin] = TMath::Power(-1./33*c28[ibin],1./8);
    //v24[ibin] = TMath::Power(qc24[ibin],1./4);
    //v26[ibin] = TMath::Power(qc26[ibin],1./6);

    cout<<"avgtrk = "<<avgtrk[ibin]<<endl;
    cout<<"v2{2} = "<<v22[ibin]<<endl;
    cout<<"v2{4} = "<<v24[ibin]<<endl;
    cout<<"v2{6} = "<<v26[ibin]<<endl;
    cout<<"v2{8} = "<<v28[ibin]<<endl;
    
    fout->cd();
    Nevent.Write("Nevent");
    avgtrk.Write("avgtrk");
    qc22.Write("qc22"); qc24.Write("qc24");   qc26.Write("qc26");   qc28.Write("qc28");
    c22.Write("c22"); c24.Write("c24");   c26.Write("c26");   c28.Write("c28");
    v22.Write("v22"); v24.Write("v24");   v26.Write("v26");   v28.Write("v28");
    TDirectory *dir0 = (TDirectory*)fout->mkdir(Form("D_%d",ibin));
    dir0->cd();
    sumreals[ibin]->Write();
    sumimags[ibin]->Write();
    weights[ibin]->Write();
    }
}
