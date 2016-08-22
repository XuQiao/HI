#include "par.h"
void savehisto(){
	string tbin = getenv("DIR");
	TFile *fout = new TFile(Form("histomerged.root",tbin.c_str()),"Recreate");
	TString dir= Form("/lio/lfs/cms/store/user/qixu/flow/PACorrHM/skim/tracknormcpt03to6/multi%s",tbin.c_str());
//	system(Form("hadd %s/mergedTree.root %s/FlowLYZHM*/*",dir.Data(),dir.Data()));
	TFile *mergefile = TFile::Open(Form("%s/mergedTree.root",dir.Data()));
	TTree *tree = (TTree*)mergefile->Get("demo/TrackTree");
	long nevent = tree->GetEntries();
	int mult, ntrk;
        double phi[10000], eta[10000], pt[10000];
	tree->SetBranchAddress("mult",&mult);
        tree->SetBranchAddress("ntrk",&ntrk);
	tree->SetBranchAddress("pt",pt);
        tree->SetBranchAddress("phi",phi);
	tree->SetBranchAddress("eta",eta);
	
	TH1D* hpt[nbin];
	TH1D* heta[nbin];
	TH1D* hphi[nbin];	
	for(int ibin=0;ibin<nbin;ibin++){
		heta[ibin] = new TH1D("heta","heta",600,-3,3);
		hphi[ibin] = new TH1D("hphi","hphi",800,-4,4);
		hpt[ibin] = new TH1D("hpt","hpt",1500,0,15);
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
			hpt[xbin]->Fill(pt[imult]);
			heta[xbin]->Fill(eta[imult]);
			hphi[xbin]->Fill(phi[imult]);
		}
	}
/*	TCanvas *c1 = new TCanvas();
	TCanvas *c2 = new TCanvas();	
	TCanvas *c3 = new TCanvas();
	c1->cd();
	c1->SetLogy();
	hpt->Draw();
	c2->cd();	
	heta->Draw();
	c3->cd();	
	hphi->Draw();
*/	
	for(int ibin=0;ibin<nbin;ibin++){
		fout->cd();
		TDirectory *dirp = fout->mkdir(Form("D_%d",ibin));
		dirp->cd();
		hpt[ibin]->Write("hpt",TObject::kOverwrite);
		heta[ibin]->Write("heta",TObject::kOverwrite);
		hphi[ibin]->Write("hphi",TObject::kOverwrite);
	}
	fout->Close();
}
