#include "TString.h"
#include "TTree.h"
#include "TFile.h"

class treeInt {
public:
	treeInt(TString);
	~treeInt();
	void Setup();
	void SetupExtra();
	int GetEntries();
	void GetEntry(int);
	double pt[10000],eta[10000],phi[10000];
	Int_t mult,ntrk;

private:
	TFile *f;
	TTree* tree;
};

treeInt::treeInt(TString _filename){
        f = TFile::Open(_filename);
}

treeInt::~treeInt(){
	delete tree;
	f->Close();
//	delete f;
}

void 
treeInt::Setup(){
        tree=(TTree*)f->Get("demo/TrackTree");
        tree->SetBranchAddress("mult",&mult);
        tree->SetBranchAddress("ntrk",&ntrk);
        tree->SetBranchAddress("pt",pt);
        tree->SetBranchAddress("eta",eta);
        tree->SetBranchAddress("phi",phi);
        //f->Close();
}

void treeInt::SetupExtra(){
        tree->LoadBaskets();
        tree->SetMaxVirtualSize(4e10);
}

int treeInt::GetEntries(){
	return tree->GetEntries();
}

void treeInt::GetEntry(int _Ev){
	if(_Ev>=this->GetEntries()) exit(0);
//	if(_Ev%50000==0) cout<<"Processed "<<_Ev<<" events"<<endl;
	tree->GetEntry(_Ev);
}
