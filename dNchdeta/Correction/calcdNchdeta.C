#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TVectorD.h"
#include "TObject.h"
#include "TObjArray.h"
#include <iostream>
#include <iomanip>
#include "SetupBinning.h"
#include "SetupTree.h"
using namespace std;
TObjArray* calcNtracks(int);
TH2F* calcwtrack(int);
TH1D* calcwevent(int);
TH1D* calcNevtsel(int);
void calcdNchdeta(int type=0, bool dotrack=1, bool doevent=1){ //-1 for genMC, 0 for recoMC, 1 for data
	TH1::SetDefaultSumw2();
        TString MCtype = "Hijing";
	TString stype;
	if(type==-1) stype=Form("%s_gen",MCtype.Data());	else if(type==0) stype=Form("%s_recoMC",MCtype.Data()); else stype="data";
	TH1D* dNchdeta = new TH1D("dNchdeta","dNchdeta",neta,etabin);
	TH1D* hwevent = calcwevent(type);
	TH1D* hdnom = new TH1D("hdnom","hdnom",neta,etabin);
	TH1D* hNevtsel = calcNevtsel(type);
	TH1D* hnom = (TH1D*)hNevtsel->Clone("hnom");
        if(doevent)
            hnom->Multiply(hwevent);
	TVectorD nom;	nom.ResizeTo(1);
	nom[0] = hnom->Integral();
	cout<<"nom = "<<nom[0]<<endl;
	double fMV = 0.;
	TObjArray* SKNtracks = calcNtracks(type);
	TH2F* hwtrack = calcwtrack(type);
	for(int ieta=0;ieta<neta;ieta++){
            hdnom->SetBinContent(ieta,0);
            hdnom->SetBinError(ieta,0);
        }
	for(int iMult=0;iMult<nMult;iMult++){
	    TH2F* htemp1 = (TH2F*)SKNtracks->At(iMult);
            if(dotrack)
                htemp1->Multiply(hwtrack);
	    TH1D* htemp2 = htemp1->ProjectionY("htemp2",htemp1->GetXaxis()->FindBin(ptmin),-1,"e");
            if(doevent)
	        htemp2->Scale(hwevent->GetBinContent(iMult+1));
	    hdnom->Add(htemp2);
	}
	for(int ieta=0;ieta<neta;ieta++){
	    dNchdeta->SetBinContent(ieta,hdnom->GetBinContent(ieta)/hdnom->GetBinWidth(ieta));
	    dNchdeta->SetBinError(ieta,hdnom->GetBinError(ieta)/hdnom->GetBinWidth(ieta));
	}
	dNchdeta->Scale(1.0/((1.+fMV)*nom[0]));
	TH2F* hNtracks = (TH2F*)SKNtracks->At(1);
        if(type>=0)
	    TFile *fout = new TFile(Form("output/Corr%s_track%d_event%d_pt%.1f.root",stype.Data(),dotrack,doevent,ptmin),"Recreate");
        else
	    TFile *fout = new TFile(Form("output/Corr%s_pt%.1f.root",stype.Data(),ptmin),"Recreate");
	fout->cd();
	dNchdeta->Write();
	hNtracks->Write();
	hwevent->Write();
	nom.Write("nom",TObject::kOverwrite);
	hwtrack->Write();
	hNevtsel->Write();
	SKNtracks->Write("SKNtracks",TObject::kOverwrite | TObject::kSingleKey);
	fout->Close();
}	

TObjArray* calcNtracks(int type){
	TH1::SetDefaultSumw2();
	treeInt *t = new treeInt(type);
	t->Setup();
	TH2F* hNtracks[nMult];
	for(int iMult=0;iMult<nMult;iMult++){
		hNtracks[iMult]  = new TH2F(Form("hNtracks_%d",iMult),Form("hNtracks_%d",iMult),npt,ptbin,neta,etabin);
	}
	TObjArray* SKNtracks = new TObjArray(0);
	int xMult=-1;
	for(int Ev=0; Ev<t->GetEntries();Ev++){
		t->GetEntry(Ev);
		double M0=0;
		if(type>=0){
			if(!(t->HLT_PAZeroBiasPixel_SingleTrack_v1 && t->pVertexFilterCutGplus && t->pBeamScrapingFilter && t->phfPosFilter1 && t->phfNegFilter1 && t->pprimaryvertexFilter && TMath::Abs(t->vz)<15)) continue;
			for(int itrk=0;itrk<t->nTrk;itrk++){
				if(!(t->highPurity[itrk] && t->trkPtError[itrk]/t->trkPt[itrk]<0.1 && fabs(t->trkDxy1[itrk]/t->trkDxyError1[itrk])<3.0 &&  fabs(t->trkDz1[itrk]/t->trkDzError1[itrk])<3.0)) continue;
				if(t->trkPt[itrk]<=0.1 || TMath::Abs(t->trkEta[itrk])>2.4 ) continue;
					M0++;
			}
		}
		else{
			for(int imult=0;imult<t->mult;imult++){
				if(t->pt[imult]<=0.1 || TMath::Abs(t->eta[imult])>2.4 ) continue;
					M0++;
			}
		}
		for(int iMult=0;iMult<nMult;iMult++){
			if(M0>=Multbin[iMult] && M0<Multbin[iMult+1]){
				xMult = iMult;
				continue;
			}
		}
		if(xMult<0 || xMult==nMult) continue;
		if(type>=0){
			for(int itrk=0;itrk<t->nTrk;itrk++){
				if(!(t->highPurity[itrk] && t->trkPtError[itrk]/t->trkPt[itrk]<0.1 && fabs(t->trkDxy1[itrk]/t->trkDxyError1[itrk])<3.0 &&  fabs(t->trkDz1[itrk]/t->trkDzError1[itrk])<3.0)) continue;
				if(t->trkPt[itrk]<=0.1 || TMath::Abs(t->trkEta[itrk])>2.4 ) continue;
				hNtracks[xMult]->Fill(t->trkPt[itrk],t->trkEta[itrk]);
			}
		}
		else{
			for(int imult=0;imult<t->mult;imult++){
				if(t->pt[imult]<=0.1 || TMath::Abs(t->eta[imult])>2.4 || t->chg[imult]==0 ) continue;
                                hNtracks[xMult]->Fill(t->pt[imult],t->eta[imult]);
                        }
		}
	}
	for(int iMult=0;iMult<nMult;iMult++){
		SKNtracks->Add(hNtracks[iMult]);
	}
	return SKNtracks;
}

TH2F* calcwtrack(int type){
	TH1::SetDefaultSumw2();
	TString corrname1 = "/home/xuq7/HI/dNchdeta/Correction/trkEff_pp_all_42X_origin.root";
	TString corrname2 = "/home/xuq7/HI/dNchdeta/Correction/TrackCorrections_HIJING_538_OFFICIAL_Mar24.root";
	TFile *fcorr1 = TFile::Open(corrname1);
//	TFile *fcorr2 = TFile::Open(corrname2);
        TFile *fout = new TFile("trackCorr.root","Recreate");
        fcorr1->cd();
	TH2F* rTotalEff3D = (TH2F*)fcorr1->Get("rTotalEff3D");
	TGraphAsymmErrors* gFakPt = (TGraphAsymmErrors*)fcorr1->Get("gFakPt");

	TH2F* hwtrack  = new TH2F("hwtrack","hwtrack",npt,ptbin,neta,etabin);
	TH2F* hF  = new TH2F("hF","hF",npt,ptbin,neta,etabin);
	TH2F* hR  = new TH2F("hR","hR",npt,ptbin,neta,etabin);
	TH2F* heptrack  = new TH2F("heptrack","heptrack",npt,ptbin,neta,etabin);
	for(int binx=0;binx<npt;binx++){
                    for(int biny=0;biny<neta;biny++){
			if(type>=0){
//				if(rTotalEff3D->GetBinContent(binx,biny)==0)	continue;
				hF->SetBinContent(binx,biny,gFakPt->GetY()[binx]);
				hR->SetBinContent(binx,biny,1-gFakPt->GetY()[binx]);
				hF->SetBinError(hF->GetBin(binx,biny),gFakPt->GetErrorYlow(binx));
				hR->SetBinError(hR->GetBin(binx,biny),gFakPt->GetErrorYlow(binx));
				heptrack->SetBinContent(binx,biny,rTotalEff3D->GetBinContent(biny,binx));
				heptrack->SetBinError(heptrack->GetBin(binx,biny), rTotalEff3D->GetBinError(biny,binx));

			}
			else{
				hwtrack->SetBinContent(binx,biny,1.);
				hwtrack->SetBinError(binx,biny,0.);
			}
		}
	}
	if(type>=0){
		hwtrack = (TH2F*)hR->Clone("hwtrack");
                hwtrack->SetName("hwtrack");
		hwtrack->Divide(heptrack);
	}
/*	for(int binx=0;binx<npt;binx++){
                    for(int biny=0;biny<neta;biny++){
                        if(binx<10 && biny<10)
                        cout<<binx<<"\t"<<biny<<"\t"<<hF->GetBinContent(binx,biny)<<"\t"<<hR->GetBinContent(binx,biny)<<"\t"<<heptrack->GetBinContent(binx,biny)<<"\t"<<hwtrack->GetBinContent(binx,biny)<<"\t"<<rTotalEff3D->GetBinContent(binx,biny)<<endl;
                    }
        }
        heptrack->Draw("colz");*/
        fout->cd();
        hF->Write();
        hR->Write();
        heptrack->Write();
        hwtrack->Write();
	return hwtrack;
        fout->Close();

}

TH1D* calcwevent(int type){
	TH1::SetDefaultSumw2();
        treeInt *t = new treeInt(type);
        t->Setup();
	TH1D* hNevent_bef = new TH1D("hNevent_bef","hNevent_bef",nMult,Multbin);
	TH1D* hNevent_aft = new TH1D("hNevent_aft","hNevent_aft",nMult,Multbin);
	TH1D* hepPV = new TH1D("hepPV","hepPV",nMult,Multbin);
	TH1D* hwevent = new TH1D("hwevent","hwevent",nMult,Multbin);
	if(type>=0){
		for(int iMult=1;iMult<=nMult;iMult++){
			hepPV->SetBinContent(iMult,1.);
			hepPV->SetBinError(iMult,0.);
		}
		for(int Ev=0; Ev<t->GetEntries();Ev++){
        		t->GetEntry(Ev);
			double M0=0;
			for(int itrk=0;itrk<t->nTrk;itrk++){
				if(!(t->highPurity[itrk] && t->trkPtError[itrk]/t->trkPt[itrk]<0.1 && fabs(t->trkDxy1[itrk]/t->trkDxyError1[itrk])<3.0 &&  fabs(t->trkDz1[itrk]/t->trkDzError1[itrk])<3.0)) continue;
                        	if(t->trkPt[itrk]<=0.1 || TMath::Abs(t->trkEta[itrk])>2.4 ) continue;
	                        	M0++;
			}
        	        if(!(TMath::Abs(t->vz)<15)) continue;
			hNevent_bef->Fill(M0);
                	if(!(t->HLT_PAZeroBiasPixel_SingleTrack_v1 && t->pVertexFilterCutGplus && t->pBeamScrapingFilter && t->phfPosFilter1 && t->phfNegFilter1 && t->pprimaryvertexFilter && TMath::Abs(t->vz)<15)) continue;
			hNevent_aft->Fill(M0);
		}
		TH1D* heptrig = (TH1D*)hNevent_aft->Clone("heptrig");
		heptrig->Divide(hNevent_bef);	//notice heptrig should be less than or equal to 1, BayesDivide used
		hwevent = (TH1D*)hNevent_bef->Clone("hwevent");
                hwevent->SetName("hwevent");
		TH1D* htemp = (TH1D*)hNevent_aft->Clone("htemp");
		htemp->Multiply(hepPV);
		hwevent->Divide(htemp);	
	}
	else{
	for(int iMult=1;iMult<=nMult;iMult++){
        	hwevent->SetBinContent(iMult,1.);
                hwevent->SetBinError(iMult,0.);
                }
	}
	return hwevent;


}

TH1D* calcNevtsel(int type){
	TH1::SetDefaultSumw2();
	TH1D* hNevtsel = new TH1D("hNevtsel","hNevtsel",nMult,Multbin);
        treeInt *t = new treeInt(type);
        t->Setup();
	for(int Ev=0; Ev<t->GetEntries();Ev++){
                t->GetEntry(Ev);
                double M0=0;
		if(type>=0){
        	        if(!(t->HLT_PAZeroBiasPixel_SingleTrack_v1 && t->pVertexFilterCutGplus && t->pBeamScrapingFilter && t->phfPosFilter1 && t->phfNegFilter1 && t->pprimaryvertexFilter && TMath::Abs(t->vz)<15)) continue;
	                for(int itrk=0;itrk<t->nTrk;itrk++){
				if(!(t->highPurity[itrk] && t->trkPtError[itrk]/t->trkPt[itrk]<0.1 && fabs(t->trkDxy1[itrk]/t->trkDxyError1[itrk])<3.0 &&  fabs(t->trkDz1[itrk]/t->trkDzError1[itrk])<3.0)) continue;
        	                if(t->trkPt[itrk]<=0.1 || TMath::Abs(t->trkEta[itrk])>2.4 ) continue;
                	                M0++;
	                }
               		hNevtsel->Fill(M0);
        	}
		else{
	                for(int imult=0;imult<t->mult;imult++){
				 if(t->pt[imult]<=0.1 || TMath::Abs(t->eta[imult])>2.4 || t->chg[imult]==0 ) continue;
                                        M0++;
                        }
                        hNevtsel->Fill(M0);
		}
	}
	return hNevtsel;
}
