#include <vector>
//#include "/afs/cern.ch/user/q/qixu/CMSSW_6_2_5/src/Centrality/NBD_Glauber_fit/NBD/NBDclass.C"
#include "/afs/cern.ch/user/q/qixu/CMSSW_6_2_5/src/Centrality/NBD_Glauber_fit/NBD/parameter.h"
#include "par.h"
void merge(){
	if (!TClass::GetDict("NBD")) {
		gROOT->ProcessLine(".L /afs/cern.ch/user/q/qixu/CMSSW_6_2_5/src/Centrality/NBD_Glauber_fit/NBD/NBDclass.C++");
	}

	TFile *f;
	f=TFile::Open(Form("G.root"),"ReadOnly");
	TFile *fg = new TFile("graph.root","Update");
	TDirectory *dir;
	TList *def, *tl;
	TKey *key;
	TString dirname;
	for(int sth=0; sth<3; sth++){
	if(sth==0) dirname = "std";
        else if(sth==1) dirname ="Gri055";
        else  dirname ="Gri101";
	dir = f->GetDirectory(dirname);
	def = new TList();
	int iGlau=0;
	//int iGlau=atoi(getenv("GTH"));
	tl = dir->GetListOfKeys();
	for(iGlau=0; iGlau<nGlau; iGlau++){
		key = ((TKey*) tl->At(iGlau));
		NBD* temp = (NBD*)key->ReadObj();
		def->Add(temp);
	}
	
	NBD *n0 = (NBD*)def->At(0);
	n0->calcsys(def);
	TString op;
	for(int option=0;option<3;option++){
	if(option==0) op=Form("%s_Ncoll",dirname.Data());
	else if(option==1) op=Form("%s_Npart",dirname.Data());
	else op=Form("%s_B",dirname.Data());
	TGraphErrors* gr = n0->DrawGraph(option);
	gr->Write(op);
	}
	}
	fg->Close();
	f->Close();
}
