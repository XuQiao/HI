#include "/afs/cern.ch/user/q/qixu/CMSSW_6_2_5/src/Centrality/NBD_Glauber_fit/NBD/parameter.h"
#include "par.h"
void getpar(){ 
	   if (!TClass::GetDict("NBD")) {
                gROOT->ProcessLine(".L /afs/cern.ch/user/q/qixu/CMSSW_6_2_5/src/Centrality/NBD_Glauber_fit/NBD/NBDclass.C+");
        }
/*
	TFile *f = new TFile("G.root","Recreate");
	TString str;
	for(int sth=0;sth<3;sth++){
		if(sth==0)	str="std";
		else if(sth==1)	str="Gri055";
		else	str="Gri101";
        	TDirectory *dir = f->mkdir(str);
		for(int Gth=0; Gth<nGlau; Gth++){
			TFile *fs = TFile::Open(Form("G%d%d.root",sth,Gth),"Read");
			TList *tl = fs->GetListOfKeys();
			TKey *key = (TKey*)tl->FindObject("Glau");
			NBD *temp = (NBD*)key->ReadObj();
			temp->SetName(Form("Glau%d",Gth));
			dir->cd();
			temp->Write();
			fs->Close();
		}
	}
		f->Close();


        TFile *f0 = TFile::Open("G.root");
        for(int sth=0;sth<3;sth++){
                if(sth==0)      str="std";
                else if(sth==1) str="Gri055";
                else    str="Gri101";
	TDirectory *dir = f0->GetDirectory(str);
        TList *def = new TList();
	TList *tl = dir->GetListOfKeys();
        for(int iGlau=0; iGlau<nGlau; iGlau++){
                key = (TKey*) tl->FindObject(Form("Glau%d",iGlau));
                NBD* temp = (NBD*)key->ReadObj();
                def->Add(temp);
        }
	cout<<dir->GetName()<<endl;
        for(int iGlau=0; iGlau<nGlau; iGlau++){
        	NBD *n = (NBD*)def->At(iGlau);
		double mubest = n->mubest;
		double kbest = n->kbest;
		cout<<mubest<<","<<kbest<<"},{";
	}
	}
*/
	int sth=atoi(getenv("STH"));
	int Gth=atoi(getenv("GTH"));
	TFile *f0 = TFile::Open(Form("G%d%d.root",sth,Gth));
        TList *tl = f0->GetListOfKeys();
                key = (TKey*) tl->At(0);
                NBD* n = (NBD*)key->ReadObj();
                cout<<n->mubest<<","<<n->kbest<<endl;
                cout<<n->chis<<","<<n->Ndf<<endl;
	for(int i=0;i<N;i++){
	cout<<i<<" kpoint \t"<<n->kpoint[i]<<endl;
	cout<<i<<" centbin \t"<<n->centbin[i]<<endl;
	}
	for(int i=0;i<N-1;i++){
	cout<<i<<"\t"<<n->NcollAver[i]<<endl;
	}
}

