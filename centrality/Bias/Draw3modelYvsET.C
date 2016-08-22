#include "parameter.h"
#include "par.h"
void Draw3modelYvsET(){
	gStyle->SetOptStat(kFALSE);	
	int sth=0, Gth=0;
        TFile *f = TFile::Open(outG);
        if(sth==0){TString dirname = "std";}
        else if(sth==1){TString dirname ="Gri055";}
        else {TString dirname ="Gri101";}
        TObjString* dataname = (TObjString*)f->Get(Form("dataname"));
        TObjString* histoname = (TObjString*)f->Get(Form("histoname"));
                TString name;
        if(Gth==0)
                name = "G0";
        else if(Gth<nGlau)
                name = Form("Glau_%d",Gth);
        else
                name = Form("bin_%d",Gth-nGlau+1);
        TObjString* Glaubername = (TObjString*)f->Get(Form("%s/%s/Glaubername",dirname.Data(),name.Data()));
        TVectorD* k0 = (TVectorD*)f->Get(Form("%s/%s/k0",dirname.Data(),name.Data()));
        TVectorD* theta0 = (TVectorD*)f->Get(Form("%s/%s/theta0",dirname.Data(),name.Data()));
        TVectorD* xmin = (TVectorD*)f->Get(Form("%s/%s/xmin",dirname.Data(),name.Data()));
        TVectorD* xmax = (TVectorD*)f->Get(Form("%s/%s/xmax",dirname.Data(),name.Data()));
        TVectorD* thetabest = (TVectorD*)f->Get(Form("%s/%s/thetabest",dirname.Data(),name.Data()));
        TVectorD* kbest = (TVectorD*)f->Get(Form("%s/%s/kbest",dirname.Data(),name.Data()));
	
	TFile *fGlauber = TFile::Open(Glaubername->GetName());
	//(*k0)[0]=1.39;	(*kbest)[0]=0.425;
	//(*theta0)[0]=3.41;	(*thetabest)[0]=1.30;
	int xNcoll = 5;
	int xNpart = xNcoll + 1; 
        double k=(*k0)[0]+(*kbest)[0]*(xNpart-2);
        double theta=(*theta0)[0]+(*thetabest)[0]*TMath::Log(xNpart-1);
        gammafun = new TF1("gammafun","TMath::GammaDist(x,[0],0,[1])",0,100);
	gammafunevt = new TF1("gammafunevt","TMath::GammaDist(x,[0],0,[1])",0,100);
	gammafunnucl = new TF1("gammafunnucl","TMath::GammaDist(x,[0],0,[1])",0,100);
	gammafunnuclNcoll = new TF1("gammafunnuclNcoll","TMath::GammaDist(x,[0],0,[1])",0,100);
        gammafun->SetParameter(0,k);   //[1]: k value
        gammafun->SetParameter(1,theta);       //[2]: theta value
	double kevt = (*k0)[0]-(*kbest)[0];
        gammafunevt->SetParameter(0,kevt);   //[1]: k value
        gammafunevt->SetParameter(1,theta);       //[2]: theta value
        gammafunnucl->SetParameter(0,(*kbest)[0]);   //[1]: k value
        gammafunnucl->SetParameter(1,theta);   //[2]: theta value
        gammafunnuclNcoll->SetParameter(0,(*kbest)[0]*xNcoll);   //[1]: k value
        gammafunnuclNcoll->SetParameter(1,theta);   //[2]: theta value
	cout<<k<<"\t"<<theta<<endl;
	cout<<kevt<<"\t"<<(*kbest)[0]<<endl;
	cout<<"mean="<<gammafunnucl->Mean(0,100)<<endl;
	TCanvas *c1 = new TCanvas();
	gammafun->SetLineColor(1);
	gammafunevt->SetLineColor(2);
	gammafunnucl->SetLineColor(4);
	TLegend *leg = new TLegend(0.4,0.55,0.7,0.85);
	leg->AddEntry(gammafun,"gamma function for total","lp");
	leg->AddEntry(gammafunevt,"gamma function for proton","lp");
	leg->AddEntry(gammafunnucl,"gamma function for nucleus","lp");
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
	leg->SetTextSize(0.04);
	gammafun->GetYaxis()->SetRangeUser(0,0.5);
	gammafun->SetTitle("gamma function");
	gammafun->GetXaxis()->SetTitle("x");
	gammafun->GetYaxis()->SetTitle("Probability");
	gammafun->Draw();	
	gammafunevt->Draw("same");	
	gammafunnucl->Draw("same");	
	leg->Draw("same");


	TH2D *UCM = new TH2D("UCM","UCM",1000,0,100,2000,0,40);
	TH2D *PCM = new TH2D("PCM","PCM",1000,0,100,2000,0,40);
	TH2D *VCM = new TH2D("VCM","VCM",1000,0,100,2000,0,40);
	TH1D *h = new TH1D("","",10000,0,100);
	TTree *t = (TTree*)fGlauber->Get("nt_p_Pb");
        Float_t Ncoll, Npart, B;        Long_t Nevent;
	t->SetBranchAddress("Ncoll",&Ncoll);
        t->SetBranchAddress("Npart",&Npart);
        t->SetBranchAddress("B",&B);

        Nevent = (Long_t) t->GetEntries();
        Long_t Ev;      Int_t Bino;     Double_t Para_nucl, Para_p, Para_evt, Bi_Para_nucl, Bi_Para_evt;
        for (Ev=0; Ev<1e7; Ev++){
     	  	if(Ev%1000000==0)       cout<<"\t"<<"Have run "<<Ev<<" events"<<endl;
 //       	t->GetEntry(Ev);
//		if(Npart!=xNpart) continue;
		//Para_nucl = 0; //make sure that Para doesn't accuthetalate through loops
		Para_evt = 0; //make sure that Para doesn't accuthetalate through loops
               	Para_nucl = gammafunnuclNcoll->GetRandom();
               	Para_p = gammafunnuclNcoll->GetRandom();
//                	Para_nucl += Bi_Para_nucl;
//                }
                for(Bino=0; Bino<xNcoll; Bino++){
		Bi_Para_evt = gammafunevt->GetRandom();
		Para_evt += Bi_Para_evt;
		}
		h->Fill(Para_nucl);
		UCM->Fill(Para_nucl,xNcoll);
		PCM->Fill(Para_nucl,1.0/2/(*kbest)[0]/theta*(Para_nucl+Para_p));
		VCM->Fill(Para_nucl,1.0/2*(Para_nucl/(*kbest)[0]/theta+xNcoll));
	}
	TCanvas *c2 = new TCanvas();
	c2->SetLogz();
	PCM->GetYaxis()->SetRangeUser(0,16.5);
	PCM->GetYaxis()->SetTitle("Y/C");
	PCM->GetXaxis()->SetRangeUser(0,75);
	UCM->SetLineWidth(2);
	PCM->GetXaxis()->SetTitle("HF #Sigma E_{T}^{j}");
	PCM->SetTitle("");
	UCM->SetLineStyle(4);
	UCM->SetLineWidth(1.5);
	PCM->Draw("colz");
	UCM->Draw("Csame");
	VCM->SetLineStyle(2);
	VCM->SetLineWidth(1.5);
	VCM->Draw("Csame");

        TLegend *leg = new TLegend(0.56, 0.15, 0.8, 0.35);
        leg->SetFillColor(10);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0.035);
        leg->SetTextFont(42);
        leg->SetTextSize(0.04);
        leg->AddEntry(UCM,"Uncorrelated Model","lp");
        leg->AddEntry(PCM,"Part. Corr. Model","lp");
        leg->AddEntry(VCM,"Var. Corr. Model","lp");
        leg->Draw("same");
        TLatex *tex= new TLatex(0.7,0.45,Form("Ncoll = %d",xNcoll));
        tex->SetNDC();
        tex->SetTextColor(1);
        tex->SetTextSize(0.045);
        tex->SetTextFont(42);
        tex->Draw("same");


	c1->Print("gammafundist_CMS.png");
	c2->Print("paperfig2_CMS.png");
}
