#include <string>
#include "ReadFromForest.h"

IntVar::IntVar():
tree(NULL),
HFp3(-9999),
HFp4(-9999),
HFm3(-9999),
HFm4(-9999),
ET(-9999),
eta(-9999),
vtz(-9999),
Bin(-9999),
Ntrack(-9999),
vz(-9999),
HFhit(-9999),
HLT_PAZeroBiasPixel_SingleTrack_v1(-9999),
phltPixelClusterShapeFilter(-9999),
phfPosFilter1(-9999),
phfNegFilter1(-9999),
pprimaryvertexFilter(-9999),
pBeamScrapingFilter(-9999),
pVertexFilterCutGplus(-9999)
{
}

IntVar::~IntVar(){
    delete tree;
}

void IntVar::Init(vector<TString> sin){
    string tr1 = "hiEvtAnalyzer/HiTree";
    string tr2 = "ppTrack/trackTree";
    string tr3 = "skimanalysis/HltTree";
    string tr4 = "hltanalysis/HltTree";
//    string tr5 = "hfrechits/hfTree";

    tree = new TChain(tr1.c_str());
    TChain *t2 = new TChain(tr2.c_str());
    TChain *t3 = new TChain(tr3.c_str());
    TChain *t4 = new TChain(tr4.c_str());
//    TChain *t5 = new TChain(tr5.c_str());

    for(unsigned int ifile = 0;ifile < sin.size(); ifile++){
        tree -> AddFile(sin[ifile]);
        t2 -> AddFile(sin[ifile]);
        t3 -> AddFile(sin[ifile]);
        t4 -> AddFile(sin[ifile]);
//        t5 -> AddFile(sin[ifile]);
    }

    tree -> AddFriend(t2,tr2.c_str());
    tree -> AddFriend(t3,tr3.c_str());
    tree -> AddFriend(t4,tr4.c_str());
//    tree -> AddFriend(t5,tr5.c_str());
    cout << tree -> GetEntries() << endl;
    tree -> SetBranchAddress("hiNtracks",&Ntrack);
    tree -> SetBranchAddress("hiHFplusEta4",&HFp4);
    tree -> SetBranchAddress("hiHFplus",&HFp3);
    tree -> SetBranchAddress("hiHFminusEta4",&HFm4);
    tree -> SetBranchAddress("hiHFminus",&HFm3);
    tree -> SetBranchAddress("hiET",&ET);
    tree -> SetBranchAddress("hiBin",&Bin);
    tree->SetBranchAddress("vz",&vz);
    tree->SetBranchAddress("hiHFhit",&HFhit);//HF hits
  //  tree->SetBranchAddress("HLT_PAZeroBias_SinglePixelTrack_v1",&HLT_PAZeroBiasPixel_SingleTrack_v1);
    tree->SetBranchAddress("L1_MinimumBiasHF0_OR_BptxAND_Final",&HLT_PAZeroBiasPixel_SingleTrack_v1);
  //  tree->SetBranchAddress("phltPixelClusterShapeFilter",&phltPixelClusterShapeFilter);
  //  tree->SetBranchAddress("phfPosFilter1",&phfPosFilter1);
  //  tree->SetBranchAddress("phfNegFilter1",&phfNegFilter1);
    tree->SetBranchAddress("pPAprimaryVertexFilter",&pprimaryvertexFilter);
    tree->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter);
    tree->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus);
}

void IntVar::Fill(){
    TString name = "";
    int nevent = tree -> GetEntries();
    name = Form("hNtrack");
    TH1F* hNtrack = new TH1F(name,name,100,0,100);
    name = Form("hHFp3");
    TH1F* hHFp3 = new TH1F(name,name,1000,0,1000);
    name = Form("hHFp4");
    TH1F* hHFp4 = new TH1F(name,name,1000,0,1000);
    name = Form("hHFm3");
    TH1F* hHFm3 = new TH1F(name,name,1000,0,1000);
    name = Form("hHFm4");
    TH1F* hHFm4 = new TH1F(name,name,1000,0,1000);
    name = Form("hET");
    TH1F* hET = new TH1F(name,name,100,0,100);
    name = Form("hBin");
    TH1F* hBin = new TH1F(name,name,200,0,200);
    name = Form("hBinpV");
    TH2F* hBinpV = new TH2F(name,name,200,0,200,2,0,2);
    name = Form("hBinpB");
    TH2F* hBinpB = new TH2F(name,name,200,0,200,2,0,2);
    name = Form("hBinpp");
    TH2F* hBinpp = new TH2F(name,name,200,0,200,2,0,2);
    name = Form("hBinvz");
    TH2F* hBinvz = new TH2F(name,name,200,0,200,40,-20,20);
    name = Form("hvz");
    TH1F* hvz = new TH1F(name,name,200,-20,20);
    h1.push_back(hNtrack);
    h1.push_back(hHFp3);
    h1.push_back(hHFp4);
    h1.push_back(hHFm3);
    h1.push_back(hHFm4);
    h1.push_back(hET);
    h1.push_back(hBin);
    h1.push_back(hvz);
    h2.push_back(hBinpV);
    h2.push_back(hBinpB);
    h2.push_back(hBinpp);
    h2.push_back(hBinvz);
    for(int ievt = 0; ievt < nevent; ievt ++){
    tree -> GetEntry(ievt);
    if(ievt % 5000==0) cout << "ievt = " << ievt << endl;
    //if(!(HLT_PAZeroBiasPixel_SingleTrack_v1 && pVertexFilterCutGplus && pBeamScrapingFilter && phfPosFilter1 && phfNegFilter1 && pprimaryvertexFilter && TMath::Abs(vz)<15)) continue;
   // cout << HLT_PAZeroBiasPixel_SingleTrack_v1 << pVertexFilterCutGplus << pBeamScrapingFilter << pprimaryvertexFilter << vz << endl;
    if(!(pVertexFilterCutGplus && pBeamScrapingFilter && pprimaryvertexFilter && TMath::Abs(vz)<30)) continue;
        hNtrack -> Fill(Ntrack);
        hHFp3 -> Fill(HFp3);
        hHFp4 -> Fill(HFp4);
        hHFm3 -> Fill(HFm3);
        hHFm4 -> Fill(HFm4);
        hET -> Fill(ET);
        hBin -> Fill(Bin);
        hvz -> Fill(vz);
        hBinpV -> Fill(Bin, pVertexFilterCutGplus);
        hBinpB -> Fill(Bin, pBeamScrapingFilter);
        hBinpp -> Fill(Bin, pprimaryvertexFilter);
        hBinvz -> Fill(Bin, vz);
    }
}

void IntVar::Write(TString sout){
    TFile *fout = new TFile(sout.Data(),"recreate");
    fout->cd();
    for(unsigned int ih1 = 0; ih1 < h1.size(); ih1++){
        h1[ih1]->Write();
    }
    for(unsigned int ih2 = 0; ih2 < h2.size(); ih2++){
        h2[ih2]->Write();
    }
    fout->Close();
}
