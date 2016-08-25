void ComparehistosEvS(){
TFile *fold=new TFile("/scratch/xuq7/RpA/TreeAna/PPbdata_ppReco_akPu3PF_AlljetTrigKurtCombinationWeight_JetPt0noIPupperCut.root");
TFile *fnew=new TFile("/scratch/xuq7/RpA/TreeAna/skimfile.root");
const int Ntrig=5;
TString trigName[Ntrig]={"Jet20","Jet40","Jet60","Jet80","Jet100"};

TTree *ntold=fold->Get("nt");
TTree *nt=fnew->Get("nt");
TFile *fout=new TFile("codeCompareEvS.root","RECREATE");

TH1F *jetptold	= new TH1F("jetptold","jetptold",1000,0,1000); 
TH1F *jetpt	= new TH1F("jetpt","jetpt",1000,0,1000); 
TH1F* jetetaold = new TH1F("jetetaold","jetetaold",100,-5,5);
TH1F* jeteta    = new TH1F("jeteta","jeteta",100,-5,5);
TH1F* jetphi = new TH1F("jetphi","jetphi",100,-TMath::Pi(),TMath::Pi());
TH1F* jetphiold = new TH1F("jetphiold","jetphiold",100,-TMath::Pi(),TMath::Pi());
TH1F* vz = new TH1F("vz","vz",600,-30,30);
TH1I* cent = new TH1I("cent","cent",100,0,100);
TH1I *nrefold	= new TH1I("nrefold","nrefold",20,0,20); 
TH1I *nref	= new TH1I("nref","nref",20,0,20); 

double tjetptold, tjetetaold, tjetphiold, trawptold;
double tjetpt[100], tjeteta[100], tjetphi[100], trawpt[100];
int tnrefold, tnref, tcent;
double tvz;	int run;
int pHBHENoiseFilter, pprimaryvertexFilter, pPAcollisionEventSelectionPA;
int HLT_PAJet20_NoJetID_v1,HLT_PAJet40_NoJetID_v1, HLT_PAJet60_NoJetID_v1, HLT_PAJet80_NoJetID_v1, HLT_PAJet100_NoJetID_v1;

  ntold->SetBranchAddress("jtpt",&tjetptold);
  ntold->SetBranchAddress("jteta",&tjetetaold);
  ntold->SetBranchAddress("jtphi",&tjetphiold);
  ntold->SetBranchAddress("rawpt",&trawptold);
  ntold->SetBranchAddress("nref",&tnrefold);
  nt->SetBranchAddress("jtpt",tjetpt);
  nt->SetBranchAddress("jteta",tjeteta);
  nt->SetBranchAddress("jtphi",tjetphi);
  nt->SetBranchAddress("rawpt",trawpt);
  nt->SetBranchAddress("nref",&tnref);
  nt->SetBranchAddress("vz",&tvz);
  nt->SetBranchAddress("hiBin",&tcent);
  nt->SetBranchAddress("run",&run);
  nt->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);
  nt->SetBranchAddress("pPAcollisionEventSelectionPA",&pPAcollisionEventSelectionPA);
  nt->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter);
  nt->SetBranchAddress("HLT_PAJet20_noJetID_v1",&HLT_PAJet20_NoJetID_v1);
  nt->SetBranchAddress("HLT_PAJet40_noJetID_v1",&HLT_PAJet40_NoJetID_v1);
  nt->SetBranchAddress("HLT_PAJet60_noJetID_v1",&HLT_PAJet60_NoJetID_v1);
  nt->SetBranchAddress("HLT_PAJet80_noJetID_v1",&HLT_PAJet80_NoJetID_v1);
  nt->SetBranchAddress("HLT_PAJet100_noJetID_v1",&HLT_PAJet100_NoJetID_v1);



for(Long_t i; i<nt->GetEntries(); i++){
if(i%100000==0)
cout<<"has processed "<<i<<" th Events"<<endl; 
nt->GetEntry(i);
	if(fabs(vz)>15.) continue;
	if(!pHBHENoiseFilter || !pprimaryvertexFilter || !pPAcollisionEventSelectionPA ) continue;
            if(!HLT_PAJet20_NoJetID_v1 && !HLT_PAJet40_NoJetID_v1 && !HLT_PAJet60_NoJetID_v1 && !HLT_PAJet80_NoJetID_v1 && !HLT_PAJet100_NoJetID_v1 ) continue;
	for(int j4i=0;j4i<tnref;j4i++){
	if(trawpt[j4i]<22.||fabs(tjeteta[ij])>0.3) continue ;
	jetpt->Fill(tjetpt[j4i]);
	jeteta->Fill(tjeteta[j4i]);
	jetphi->Fill(tjetphi[j4i]);
	}
		
	nref->Fill(tnref);	
	vz->Fill(tvz);
	cent->Fill(tcent);
}

for(Long_t iold; iold<ntold->GetEntries(); iold++){
if(iold%100000==0)
cout<<"has processed "<<iold<<" th Events"<<endl; 
ntold->GetEntry(iold);
	if(trawptold<22.||fabs(tjetetaold)>0.3) continue ;
        jetptold->Fill(tjetptold);
        jetetaold->Fill(tjetetaold);
        jetphiold->Fill(tjetphiold);
	nrefold->Fill(tnrefold);
}

fout->cd();
jetpt->Write();
jeteta->Write();
jetphi->Write();
nref->Write();
vz->Write();
cent->Write();
jetptold->Write();
jetetaold->Write();
jetphiold->Write();
nrefold->Write();
}
