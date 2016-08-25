#include "par.h"
void savehisto(){
	string tbin = getenv("DIR");
	TFile *fout = new TFile(Form("histomerged.root",tbin.c_str()),"Recreate");
	TString dir= Form("/lio/lfs/cms/store/user/qixu/flow/PACorrHM/skim/PFcandi/crab/%s",tbin.c_str());
//	system(Form("hadd %s/mergedTree.root %s/FlowLYZHM*/*",dir.Data(),dir.Data()));
	TFile *mergefile = TFile::Open(Form("%s/mergedTree.root",dir.Data()));
	TTree *tree = (TTree*)mergefile->Get("demo/PFTree");
	long nevent = tree->GetEntries();
	int nPFpart, ntrk;
        const int nID=6;
        int pfid[10000];
        double pfphi[10000], pfeta[10000], pfpt[10000];
	tree->SetBranchAddress("nPFpart",&nPFpart);
        tree->SetBranchAddress("ntrk",&ntrk);
	tree->SetBranchAddress("pfpt",pfpt);
        tree->SetBranchAddress("pfphi",pfphi);
	tree->SetBranchAddress("pfeta",pfeta);
	tree->SetBranchAddress("pfid",pfid);
	
	TH1D* hpt[nbin][nID];
	TH1D* heta[nbin][nID];
	TH1D* hphi[nbin][nID];	
	for(int ibin=0;ibin<nbin;ibin++){
            for(int iID=0;iID<nID;iID++){
		heta[ibin][iID] = new TH1D(Form("heta_%d",iID),Form("heta_%d",iID),1000,-5,5);
		hphi[ibin][iID] = new TH1D(Form("hphi_%d",iID),Form("hphi_%d",iID),800,-4,4);
		hpt[ibin][iID] = new TH1D(Form("hpt_%d",iID),Form("hpt_%d",iID),1500,0,15);
	    }
        }
	for(long i=0;i<nevent;i++){
		tree->GetEntry(i);
		if(i%10000==0)		cout<<"has processed "<<i<<" events"<<endl;
		int xbin=-1;
		for(int j=0;j<nbin;j++)
                	if(ntrk<trkbin[j]&&ntrk>=trkbin[j+1])
                        	xbin=j;
                if(xbin<0 || xbin==nbin) continue;
		for(int imult=0;imult<nPFpart;imult++){
                    if(pfid[imult]<1 || pfid[imult]>5) continue;
			hpt[xbin][0]->Fill(pfpt[imult]);
			heta[xbin][0]->Fill(pfeta[imult]);
			hphi[xbin][0]->Fill(pfphi[imult]);
                        int xid=pfid[imult];
                        hpt[xbin][xid]->Fill(pfpt[imult]);
                        heta[xbin][xid]->Fill(pfeta[imult]);
                        hphi[xbin][xid]->Fill(pfphi[imult]);
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
            for(int iID=0;iID<nID;iID++){
		hpt[ibin][iID]->Write();
		heta[ibin][iID]->Write();
		hphi[ibin][iID]->Write();
            }
	}
	fout->Close();
}
