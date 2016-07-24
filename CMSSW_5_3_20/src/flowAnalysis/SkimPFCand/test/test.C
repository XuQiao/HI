#include "iostream"
#include "fstream"
void test(){
	TH1::SetDefaultSumw2();
	TFile *f = TFile::Open("pPbskimAnaTree.root");
	TTree *t = (TTree*)f->Get("demo/PFTree");
	TFile *fout = new TFile("histo.root","Recreate");
	ofstream fstr("sample.txt");
	int ntrk, mult, nPFpart;
	int pfid[10000];
	int countpfid[8];
	int count=0;
	double pfpt[10000], pfeta[10000], pfphi[10000], pt[10000], eta[10000], phi[10000];
	int neve = t->GetEntries();
	t->SetBranchAddress("ntrk",&ntrk);
	t->SetBranchAddress("mult",&mult);
	t->SetBranchAddress("nPFpart",&nPFpart);
	t->SetBranchAddress("pfid",pfid);
	t->SetBranchAddress("pfpt",pfpt);
	t->SetBranchAddress("pfeta",pfeta);
	t->SetBranchAddress("pfphi",pfphi);
	t->SetBranchAddress("pt",pt);
	t->SetBranchAddress("eta",eta);
	t->SetBranchAddress("phi",phi);
	TH1D* hpt = new TH1D("hpt","",200,0,20);
	TH1D* hpfpt = new TH1D("hpfpt","",200,0,20);
	TH1D* hpfpt_id1 = new TH1D("hpfpt_id1","",200,0,20);
	TH1D* heta = new TH1D("heta","",100,-5,5);
	TH1D* hpfeta = new TH1D("hpfeta","",100,-5,5);
	TH1D* hpfeta_id1 = new TH1D("hpfeta_id1","",100,-5,5);
	TH1D* hphi = new TH1D("hphi","",80,-4,4);
	TH1D* hpfphi = new TH1D("hpfphi","",80,-4,4);
	TH1D* hpfphi_id1 = new TH1D("hpfphi_id1","",80,-4,4);
	fstr<<"ntrk"<<"\t"<<"mult"<<"\t"<<"nPFpart"<<"\t"<<"ptid"<<"\t"<<"pfpt"<<"\t"<<"pfeta"<<"\t"<<"pfphi"<<"\t"<<"pt"<<"\t"<<"eta"<<"\t"<<"phi"<<endl;
	for(int i=0;i<8;i++)
		countpfid[i]=0;
	for(int i=0;i<neve;i++){
		t->GetEntry(i);
		for(int j=0;j<mult;j++){
			fstr<<ntrk<<"\t"<<mult<<"\t"<<nPFpart<<"\t"<<pfid[j]<<"\t"<<pfpt[j]<<"\t"<<pfeta[j]<<"\t"<<pfphi[j]<<"\t"<<pt[j]<<"\t"<<eta[j]<<"\t"<<phi[j]<<endl;
			hpt->Fill(pt[j]);	hpfpt->Fill(pfpt[j]);	
			heta->Fill(eta[j]);	hpfeta->Fill(pfeta[j]); 
			hphi->Fill(phi[j]);	hpfphi->Fill(pfphi[j]); 
			if(pfid[j]==1 || pfid[j]==5){
				hpfpt_id1->Fill(pfpt[j]);
				hpfeta_id1->Fill(pfeta[j]);
				hpfphi_id1->Fill(pfphi[j]);
			}
			countpfid[pfid[j]]++;
			count++;
		}
		for(int j=mult;j<nPFpart;j++){
			fstr<<ntrk<<"\t"<<mult<<"\t"<<nPFpart<<"\t"<<pfid[j]<<"\t"<<pfpt[j]<<"\t"<<pfeta[j]<<"\t"<<pfphi[j]<<"\t"<<"\t"<<"\t"<<"\t"<<"\t"<<"\t"<<endl;
			hpfpt->Fill(pfpt[j]);   
                        hpfeta->Fill(pfeta[j]);
                        hpfphi->Fill(pfphi[j]); 
		        if(pfid[j]==1 || pfid[j]==5){
                		hpfpt_id1->Fill(pfpt[j]);
	                  	hpfeta_id1->Fill(pfeta[j]);
	                       	hpfphi_id1->Fill(pfphi[j]);
                	}
			countpfid[pfid[j]]++;
		}
	}
	for(int i=0;i<8;i++)
	cout<<"particle type "<<i<<", "<<countpfid[i]<<endl;
	cout<<count<<endl;
	fout->cd();
	hpt->Write();	heta->Write();	hphi->Write();
	hpfpt->Write();	hpfeta->Write();	hpfphi->Write();
	hpfpt_id1->Write();	hpfeta_id1->Write();	hpfphi_id1->Write();
}
