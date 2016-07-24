void Comparehistos(){
TFile *fold=new TFile("/scratch/xuq7/RpA/TreeAna/PPbdataNoEvtSel_ppReco_akPu3PF_AlljetTrigKurtCombinationWeight_JetPt0noIPupperCut.root");
TFile *fnew=new TFile("/scratch/xuq7/RpA/TreeAna/skimfile.root");
const int Ntrig=5;
TString trigName[Ntrig]={"Jet20","Jet40","Jet60","Jet80","Jet100"};

TTree *ntold=fold->Get("nt");
TTree *nt=fnew->Get("nt");
TFile *fout=new TFile("codeCompare.root","RECREATE");

TH1F *jetptold	= new TH1F("jetptold","jetptold",1000,0,1000); 
TH1F *jetpt	= new TH1F("jetpt","jetpt",1000,0,1000); 
TH1F* jetetaold = new TH1F("jetetaold","jetetaold",100,-5,5);
TH1F* jeteta    = new TH1F("jeteta","jeteta",100,-5,5);
TH1F* vz = new TH1F("vz","vz",600,-30,30);
TH1I* cent = new TH1I("cent","cent",100,0,100);
TH1I *nrefold	= new TH1I("nrefold","nrefold",20,0,20); 
TH1I *nref	= new TH1I("nref","nref",20,0,20); 

double tjetptold, tjetetaold;
double tjetpt[100], tjeteta[100];
int tnrefold, tnref, tcent;
double tvz;
  ntold->SetBranchAddress("jtpt",&tjetptold);
  ntold->SetBranchAddress("jteta",&tjetetaold);
  ntold->SetBranchAddress("nref",&tnrefold);
  nt->SetBranchAddress("jtpt",tjetpt);
  nt->SetBranchAddress("jteta",tjeteta);
  nt->SetBranchAddress("nref",&tnref);
  nt->SetBranchAddress("vz",&tvz);
  nt->SetBranchAddress("hiBin",&tcent);

for(Long_t i; i<nt->GetEntries(); i++){
if(i%100000==0)
cout<<"has processed "<<i<<" th Events"<<endl; 
nt->GetEntry(i);
	for(int j4i=0;j4i<tnref;j4i++){
	jetpt->Fill(tjetpt[j4i]);
	jeteta->Fill(tjeteta[j4i]);
	}
		
	nref->Fill(tnref);	
	vz->Fill(tvz);
	cent->Fill(tcent);
}

for(Long_t iold; iold<ntold->GetEntries(); iold++){
if(iold%100000==0)
cout<<"has processed "<<iold<<" th Events"<<endl; 
ntold->GetEntry(iold);
        jetptold->Fill(tjetptold);
        jetetaold->Fill(tjetetaold);
	nrefold->Fill(tnrefold);
}

fout->cd();
jetpt->Write();
jeteta->Write();
nref->Write();
vz->Write();
cent->Write();
jetptold->Write();
jetetaold->Write();
nrefold->Write();
}
