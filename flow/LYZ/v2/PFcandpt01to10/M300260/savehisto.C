#include "par.h"
void savehisto(){
	string tbin = getenv("DIR");
	TFile *fout = new TFile(Form("histomerged.root",tbin.c_str()),"Recreate");
	TString dir= Form("/lio/lfs/cms/store/user/qixu/flow/PACorrHM/skim/PFcandi/crab/%s",tbin.c_str());
	//system(Form("hadd %s/mergedTree.root %s/FlowHM*/*",dir.Data(),dir.Data()));
	TFile *mergefile = TFile::Open(Form("%s/mergedTree.root",dir.Data()));
	TTree *tree = (TTree*)mergefile->Get("demo/PFTree");
	long nevent = tree->GetEntries();
	int mult, ntrk;
        double phi[10000], eta[10000], pt[10000];
	int id[10000];
	tree->SetBranchAddress("nPFpart",&mult);
        tree->SetBranchAddress("ntrk",&ntrk);
	tree->SetBranchAddress("pfpt",pt);
        tree->SetBranchAddress("pfphi",phi);
	tree->SetBranchAddress("pfeta",eta);
	tree->SetBranchAddress("pfid",id);
	
	TH1D* hpt[nbin];
	TH1D* heta[nbin];
	TH1D* hphi[nbin];	
	TH1D* hid[nbin];	
	for(int ibin=0;ibin<nbin;ibin++){
		heta[ibin] = new TH1D("heta","heta",1000,-5,5);
		hphi[ibin] = new TH1D("hphi","hphi",800,-4,4);
		hpt[ibin] = new TH1D("hpt","hpt",1500,0,15);
		hid[ibin] = new TH1D("hid","hid",10,0,10);
	}
	for(long i=0;i<nevent;i++){
		tree->GetEntry(i);
		if(i%10000==0)		cout<<"has processed "<<i<<" events"<<endl;
		int xbin=-1;
		for(int j=0;j<nbin;j++)
                	if(ntrk<trkbin[j]&&ntrk>=trkbin[j+1])
                        	xbin=j;
                if(xbin<0 || xbin==nbin) continue;
		for(int imult=0;imult<mult;imult++){
			hid[xbin]->Fill(id[imult]);
			if(id[imult]>5 || id[imult]<1) continue;
			hpt[xbin]->Fill(pt[imult]);
			heta[xbin]->Fill(eta[imult]);
			hphi[xbin]->Fill(phi[imult]);
		}
	}
	
	for(int ibin=0;ibin<nbin;ibin++){
		fout->cd();
		TDirectory *dirp = fout->mkdir(Form("D_%d",ibin));
		dirp->cd();
		hpt[ibin]->Write("hpfpt",TObject::kOverwrite);
		heta[ibin]->Write("hpfeta",TObject::kOverwrite);
		hphi[ibin]->Write("hpfphi",TObject::kOverwrite);
		hid[ibin]->Write("hpfid",TObject::kOverwrite);
	}
	fout->Close();
}
