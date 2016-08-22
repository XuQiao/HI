void Calc(Char_t *Generatorfile,Char_t *outfile, Char_t *condition, ofstream &fouttxt){

fouttxt<<condition<<" efficiency:"<<endl;

TFile *fin = new TFile(Generatorfile);
TFile *fout =  new TFile(outfile,"update"); 

TH1F *hHFEnergy4 = new TH1F("hHFEnergy","HFEnergy Distribution before Trigger selection",1000,0,100);
hHFEnergy4->Sumw2();
TH1F *hHFEnergy4_tr = new TH1F("hHFEnergy_tr","HFEnergy Distribution after Trigger selection",1000,0,100); 
hHFEnergy4_tr->Sumw2();
TH1F *hHFEnergyPlus4 = new TH1F("hHFEnergyPlus4","HFEnergy Distribution before Trigger selection",1000,0,100);
hHFEnergyPlus4->Sumw2();
TH1F *hHFEnergyPlus4_tr = new TH1F("hHFEnergyPlus4_tr","HFEnergy Distribution after Trigger selection",1000,0,100); 
hHFEnergyPlus4_tr->Sumw2();
TH1F *hNtrack = new TH1F("hNtrack","Ntrack offline Distribution before Trigger selection",300,0,300);
hNtrack->Sumw2();
TH1F *hNtrack_tr = new TH1F("hNtrack_tr","Ntrack offline Distribution after Trigger selection",300,0,300); 
hNtrack_tr->Sumw2();

TTree *tree;
fin->cd();
tree=(TTree*)hiEvtAnalyzer->Get("HiTree");
tree->AddFriend("pptrack=pptracks/trackTree");
tree->AddFriend("skimTree=skimanalysis/HltTree");
tree->AddFriend("hltTree=hltanalysis/HltTree");
tree->AddFriend("hfTree=hfrechits/hfTree");

Float_t hiHF, vz;
Float_t hiHFminus4, hiHFplus4;
Long_t Ev,Nevent;
Int_t Ntrack, n, HLT_PAZeroBiasPixel_SingleTrack_v1, phltPixelClusterShapeFilter, phfPosFilter1, phfNegFilter1, pprimaryvertexFilter, pBeamScrapingFilter, pVertexFilterCutGplus;
Int_t Nskim=0;

Nevent=tree->GetEntries();
cout<<Nevent<<endl;
tree->SetBranchAddress("hiNtracks",&Ntrack);//number of tracks
tree->SetBranchAddress("hiHF",&hiHF);//HF energy
tree->SetBranchAddress("hiHFplusEta4",&hiHFplus4);//HF energy positive 4 to 5
tree->SetBranchAddress("hiHFminusEta4",&hiHFminus4);//HF energy negative -5 to -4
tree->SetBranchAddress("vz",&vz);
tree->SetBranchAddress("n",&n);//HF hits
tree->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1",&HLT_PAZeroBiasPixel_SingleTrack_v1);
tree->SetBranchAddress("phltPixelClusterShapeFilter",&phltPixelClusterShapeFilter);
tree->SetBranchAddress("phfPosFilter1",&phfPosFilter1);
tree->SetBranchAddress("phfNegFilter1",&phfNegFilter1);
tree->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter);
tree->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter);
tree->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus);

for(Ev=0;Ev<Nevent;Ev++){
if(Ev%5000==0)
cout<<"Have run "<<Ev<<" Events"<<endl;

tree->GetEntry(Ev);

hHFEnergy4->Fill(hiHFplus4+hiHFminus4);
hHFEnergyPlus4->Fill(hiHFplus4);
hNtrack->Fill(Ntrack);

if(condition=="TR"||condition=="ES"||condition=="TRandES"){

if(condition=="TR" && !(HLT_PAZeroBiasPixel_SingleTrack_v1)){
	Nskim++;
	continue;
	}

else if(condition=="ES" && !(pVertexFilterCutGplus && pBeamScrapingFilter && phfPosFilter1 && phfNegFilter1 && pprimaryvertexFilter && TMath::Abs(vz)<15)){
	Nskim++;
	continue;
	}	

else if(condition=="TRandES" && !(HLT_PAZeroBiasPixel_SingleTrack_v1 && pVertexFilterCutGplus && pBeamScrapingFilter && phfPosFilter1 && phfNegFilter1 && pprimaryvertexFilter && TMath::Abs(vz)<15)){
	Nskim++;
	continue;
	}

else{
hHFEnergy4_tr->Fill(hiHFplus4+hiHFminus4);
hHFEnergyPlus4_tr->Fill(hiHFplus4);
hNtrack_tr->Fill(Ntrack);
}

}

else{
cout<<"Error!"<<endl;
break;
}

}

fouttxt<<outfile<<endl;
fouttxt<<"before "<<condition<<" efficiency correction:"<<hHFEnergy4->GetEntries()<<endl;
fouttxt<<"after "<<condition<<" efficiency correction:"<<hHFEnergy4_tr->GetEntries()<<endl;
fouttxt<<"ratio:"<<(double)hHFEnergy4_tr->GetEntries()/hHFEnergy4->GetEntries()<<endl;
fouttxt<<"-----------------------------------------------"<<endl;

fout->cd();
hHFEnergy4->Write("hHFEnergy4",TObject::kOverwrite);
hHFEnergy4_tr->Write("hHFEnergy4_tr",TObject::kOverwrite);

hHFEnergyPlus4->Write("hHFEnergyPlus4",TObject::kOverwrite);
hHFEnergyPlus4_tr->Write("hHFEnergyPlus4_tr",TObject::kOverwrite);

hNtrack->Write("hNtrack",TObject::kOverwrite);
hNtrack_tr->Write("hNtrack_tr",TObject::kOverwrite);

//---------------------------------------------------------------------

TH1F *hEff4 = (TH1F*)hHFEnergy4_tr->Clone();
TGraphAsymmErrors *gEff4 = new TGraphAsymmErrors();
hEff4->SetTitle("Double side |#eta|>4");
hEff4->GetXaxis()->SetTitle("HFEnergy");
hEff4->GetYaxis()->SetTitle("Trigger Efficiency");
gEff4->BayesDivide(hHFEnergy4_tr,hHFEnergy4);
TH1F *hEffPlus4 = (TH1F*)hHFEnergyPlus4_tr->Clone();
hEffPlus4->SetTitle("#eta>4");
hEffPlus4->GetXaxis()->SetTitle("HFEnergy");
hEffPlus4->GetYaxis()->SetTitle("Trigger Efficiency");
TH1F *hEffNtrack = (TH1F*)hNtrack_tr->Clone();
hEffNtrack->SetTitle("");
hEffNtrack->GetXaxis()->SetTitle("Ntrack offline");
hEffNtrack->GetYaxis()->SetTitle("Trigger Efficiency");

for(int ibin=1;ibin<=hHFEnergy4->GetNbinsX();ibin++){
gEff4->SetPointEXlow(ibin,0);
gEff4->SetPointEXhigh(ibin,0);}

hEffPlus4->Divide(hHFEnergyPlus4);
hEffNtrack->Divide(hNtrack);

hEff4->Write("hEff4",TObject::kOverwrite);
hEffPlus4->Write("hEffPlus4",TObject::kOverwrite);
hEffNtrack->Write("hEffNtrack",TObject::kOverwrite);
gEff4->Write("gEff4",TObject::kOverwrite);

fin->Close();
fout->Close();

gROOT->Reset();

}

void Trig_eff(){

Char_t *Hijing_file = "/store/user/tuos/pPb_MC_MinBiasTree_v4_Hijing_180k.root";
Char_t *Hijing_outfile1 = "pPbHist_Hijing_Es.root";
Char_t *Hijing_outfile2 = "pPbHist_Hijing_Tr.root";
Char_t *Hijing_outfile3 = "pPbHist_Hijing_TrandEs.root";

Char_t *Epos_file = "/store/user/tuos/pPb_MC_MinBiasTree_v4_Epos_189k.root";
Char_t *Epos_outfile1 = "pPbHist_Epos_Es.root";
Char_t *Epos_outfile2 = "pPbHist_Epos_Tr.root";
Char_t *Epos_outfile3 = "pPbHist_Epos_TrandEs.root";

Char_t *condition1 = "ES";
Char_t *condition2 = "TR";
Char_t *condition3 = "TRandES";

ofstream fouttxt;
Char_t *outtxt="Eff_Result.txt";
fouttxt.open(outtxt);

//Calc(Hijing_file,Hijing_outfile1,condition1,fouttxt);
//Calc(Epos_file,Epos_outfile1,condition1,fouttxt);
//Calc(Hijing_file,Hijing_outfile2,condition2,fouttxt);
//Calc(Epos_file,Epos_outfile2,condition2,fouttxt);
Calc(Hijing_file,Hijing_outfile3,condition3,fouttxt);
Calc(Epos_file,Epos_outfile3,condition3,fouttxt);

fouttxt.close();
}
