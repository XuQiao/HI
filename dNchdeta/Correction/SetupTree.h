#include "TString.h"
#include "TTree.h"
#include "TFile.h"

class treeInt {
public:
	treeInt(int);
	~treeInt();
	void Setup();
	int GetEntries();
	void GetEntry(int);
	Float_t hiHF, vz;
	Float_t hiHFminus4, hiHFplus4;
	Float_t pt[20000],eta[20000];
        Int_t chg[20000];
	Int_t HLT_PAZeroBiasPixel_SingleTrack_v1, phltPixelClusterShapeFilter, phfPosFilter1, phfNegFilter1, pprimaryvertexFilter, pBeamScrapingFilter, pVertexFilterCutGplus;
	Int_t nTrk, mult;
	Bool_t highPurity[10000];
	Float_t trkPtError[10000], trkDxy1[10000], trkDz1[10000], trkDxyError1[10000], trkDzError1[10000];
	Float_t trkPt[10000], trkEta[10000];

private:
	int type;
	TFile *f;
	TTree* hitree;
	TTree* pptree;
	TTree* skimtree;
	TTree* hlttree;
	TTree* hftree;
	TTree* gentree;
};

treeInt::treeInt(int _type){
	type = _type;
	TString filename;
	if(type<=0)
        	filename = "/store/user/tuos/pPb_MC_MinBiasTree_v4_Hijing_180k.root";
        //        filename = "/store/user/tuos/pPb_MC_MinBiasTree_v4_Epos_189k.root";
	else	
	        filename = "/store/user/tuos/pPb_MinBiasTree_v5_211256_json.root";
        f = TFile::Open(filename);
}

treeInt::~treeInt(){
	f->Close();
//	delete f;
	delete hitree;
	delete pptree;
	delete skimtree;
	delete hlttree;
	delete hftree;
        if(type<0)
	    delete gentree;
}

void 
treeInt::Setup(){
        hitree=(TTree*)f->Get("hiEvtAnalyzer/HiTree");
        pptree=(TTree*)f->Get("pptracks/trackTree");
        skimtree=(TTree*)f->Get("skimanalysis/HltTree");
        hlttree=(TTree*)f->Get("hltanalysis/HltTree");
        hftree=(TTree*)f->Get("hfrechits/hfTree");
        if(type<0)
                gentree=(TTree*)f->Get("genparticles/hi");
        //tree->SetBranchAddress("hiNtracks",&Ntrack);//number of tracks
        //tree->SetBranchAddress("hiHF",&hiHF);//HF energy
        //tree->SetBranchAddress("hiHFplusEta4",&hiHFplus4);//HF energy positive 4 to 5
        //tree->SetBranchAddress("hiHFminusEta4",&hiHFminus4);//HF energy negative -5 to -4
        hitree->SetBranchAddress("vz",&vz);
        //tree->SetBranchAddress("n",&n);//HF hits
        hlttree->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1",&HLT_PAZeroBiasPixel_SingleTrack_v1);
        skimtree->SetBranchAddress("phltPixelClusterShapeFilter",&phltPixelClusterShapeFilter);
        skimtree->SetBranchAddress("phfPosFilter1",&phfPosFilter1);
        skimtree->SetBranchAddress("phfNegFilter1",&phfNegFilter1);
        skimtree->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter);
        skimtree->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter);
        skimtree->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus);
        pptree->SetBranchAddress("nTrk",&nTrk);
        pptree->SetBranchAddress("highPurity",highPurity);
        pptree->SetBranchAddress("trkPt",trkPt);
        pptree->SetBranchAddress("trkEta",trkEta);
        pptree->SetBranchAddress("trkPtError",trkPtError);
        pptree->SetBranchAddress("trkDxy1",trkDxy1);
        pptree->SetBranchAddress("trkDxyError1",trkDxyError1);
        pptree->SetBranchAddress("trkDz1",trkDz1);
        pptree->SetBranchAddress("trkDzError1",trkDzError1);
        if(type<0){
                gentree->SetBranchAddress("mult",&mult);
                gentree->SetBranchAddress("pt",pt);
                gentree->SetBranchAddress("eta",eta);
                gentree->SetBranchAddress("chg",chg);
        }
        //f->Close();
}

int treeInt::GetEntries(){
	return hitree->GetEntries();
}

void treeInt::GetEntry(int _Ev){
	if(_Ev>=this->GetEntries()) exit(0);
	if(_Ev%50000==0) cout<<"Processed "<<_Ev<<" events"<<endl;
	hitree->GetEntry(_Ev);
	pptree->GetEntry(_Ev);
	skimtree->GetEntry(_Ev);
	hlttree->GetEntry(_Ev);
	hftree->GetEntry(_Ev);
	if(type<0)
	gentree->GetEntry(_Ev);
}
