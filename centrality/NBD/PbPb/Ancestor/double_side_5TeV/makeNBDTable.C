#include <memory>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <map>
#include <TF1.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TNtuple.h>
#include <TChain.h>
#include <TFile.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TTree.h>
#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"

using namespace std;

void makeNBDTable(const int nbins = 200, const string label = "HFtowersDoubleTrunc", const char * tag = "CentralityTable_HFtowersDoubleTrunc_NBDGlauber_sigma74_eff0_v0", int eff = 0) {

 TH1::SetDefaultSumw2();
 const char* outfilename = "tables_Glauber2015_PbPb.root";
 TFile * outf = new TFile(outfilename,"recreate");
 outf->cd();
 TDirectory* dir = outf->mkdir(tag);
 dir->cd();
 TNtuple* nt = new TNtuple("nt","","b:Npart:Ncoll");
 CentralityBins * outputTable = new CentralityBins(Form("run%d",1), tag, nbins);
 outputTable->table_.reserve(nbins);

 ofstream txtfile("out/output.txt");
 txtfile << "Input Glauber tree: " << infilename << endl << "Input MC HiForest HIJING" << endl;

 //Setting up variables and branches
 double binboundaries[nbins+1];
 vector<float> values;

 float b, npart, ncoll, nhard, parameter;

 t->SetBranchAddress("B",&b);
 t->SetBranchAddress("Npart",&npart);
 t->SetBranchAddress("Ncoll",&ncoll);
 nhard = 0;

 //Event loop 1
 unsigned int Nevents = t->GetEntries();
 txtfile << "Number of events = " << Nevents << endl << endl;
 for(unsigned int iev = 0; iev < Nevents; iev++) {
   if(iev%20000 == 0) cout<<"Processing event: " << iev << endl;
   t->GetEntry(iev);

   if (binHFplusTrunc) parameter = getHFplusByPt(npart);
   if (binHFminusTrunc) parameter = getHFminusByPt(npart);
   if (binTracks) parameter = getTracksByGen(npart);

   values.push_back(parameter);
 }

 if(label.compare("b") == 0) sort(values.begin(),values.end(),descend);
 else sort(values.begin(),values.end());

 //Finding the bin boundaries
 txtfile << "-------------------------------------" << endl;
 txtfile << label.data() << " based cuts are: " << endl;
 txtfile << "(";

 int size = values.size();
 binboundaries[nbins] = values[size-1];

 for(int i = 0; i < nbins; i++) {
   int entry = (int)(i*(size/nbins));
   if(entry < 0 || i == 0) binboundaries[i] = 0;
   else binboundaries[i] = values[entry];
   txtfile << binboundaries[i] << ", ";
 }
 txtfile << binboundaries[nbins] << ")" << endl << "-------------------------------------" << endl;

 // Determining Glauber results in various bins
 TH2D* hNpart = new TH2D("hNpart","",nbins,binboundaries,40,0,40);
 TH2D* hNcoll = new TH2D("hNcoll","",nbins,binboundaries,40,0,40);
 TH2D* hNhard = new TH2D("hNhard","",nbins,binboundaries,50,0,50);
 TH2D* hb = new TH2D("hb","",nbins,binboundaries,600,0,30);

 for(unsigned int iev = 0; iev < Nevents; iev++) {
   if( iev % 20000 == 0 ) cout<<"Processing event : " << iev << endl;
   t->GetEntry(iev);

   if (binHFplusTrunc) parameter = getHFplusByPt(npart);
   if (binHFminusTrunc) parameter = getHFminusByPt(npart);
   if (binTracks) parameter = getTracksByGen(npart);

   hNpart->Fill(parameter,npart);
   hNcoll->Fill(parameter,ncoll);
   hNhard->Fill(parameter,nhard);
   hb->Fill(parameter,b);

   int bin = hNpart->GetXaxis()->FindBin(parameter) - 1;
   if(bin < 0) bin = 0;
   if(bin >= nbins) bin = nbins - 1;
   nt->Fill(parameter, et, bin, b, npart, ncoll, nhard);
 }

 TCanvas *cf = new TCanvas();
 TF1* fGaus = new TF1("fb","gaus(0)",0,2);
 fitSlices(hNpart,fGaus);
 fitSlices(hNcoll,fGaus);
 fitSlices(hNhard,fGaus);
 fitSlices(hb,fGaus);

 TH1D* hNpartMean = (TH1D*)gDirectory->Get("hNpart_1");
 TH1D* hNpartSigma = (TH1D*)gDirectory->Get("hNpart_2");
 TH1D* hNcollMean = (TH1D*)gDirectory->Get("hNcoll_1");
 TH1D* hNcollSigma = (TH1D*)gDirectory->Get("hNcoll_2");
 TH1D* hNhardMean = (TH1D*)gDirectory->Get("hNhard_1");
 TH1D* hNhardSigma = (TH1D*)gDirectory->Get("hNhard_2");
 TH1D* hbMean = (TH1D*)gDirectory->Get("hb_1");
 TH1D* hbSigma = (TH1D*)gDirectory->Get("hb_2");

 txtfile<<"-------------------------------------"<<endl;
 txtfile<<"# Bin NpartMean NpartSigma NcollMean NcollSigma bMean bSigma BinEdge"<<endl;
 for(int i = 0; i < nbins; i++){
   int ii = nbins-i;
   outputTable->table_[i].n_part_mean = hNpartMean->GetBinContent(ii);
   outputTable->table_[i].n_part_var = hNpartSigma->GetBinContent(ii);
   outputTable->table_[i].n_coll_mean = hNcollMean->GetBinContent(ii);
   outputTable->table_[i].n_coll_var = hNcollSigma->GetBinContent(ii);
   outputTable->table_[i].b_mean = hbMean->GetBinContent(ii);
   outputTable->table_[i].b_var = hbSigma->GetBinContent(ii);
   outputTable->table_[i].n_hard_mean = hNhardMean->GetBinContent(ii);
   outputTable->table_[i].n_hard_var = hNhardSigma->GetBinContent(ii);
   outputTable->table_[i].bin_edge = binboundaries[ii-1];

   txtfile << i << " " << hNpartMean->GetBinContent(ii) << " " << hNpartSigma->GetBinContent(ii) << " " << hNcollMean->GetBinContent(ii) << " " << hNcollSigma->GetBinContent(ii) << " " << hbMean->GetBinContent(ii) << " " <<hbSigma->GetBinContent(ii) << " " << binboundaries[ii-1] << " " <<endl;
 }
 txtfile<<"-------------------------------------"<<endl;

 outf->cd();
 dir->cd();
 outputTable->Write();
 nt->Write();
 for(int ih = 0; ih < nhist; ih++) {
   Npart_vs_PtGenPlusEta4[ih]->Write();
   PtGenPlusEta4_vs_HFplusEta4[ih]->Write();
   Npart_vs_PtGenMinusEta4[ih]->Write();
   PtGenMinusEta4_vs_HFminusEta4[ih]->Write();
   Npart_vs_Ngentrk[ih]->Write();
   Ngentrk_vs_Ntracks[ih]->Write();
 }
 outf->Write();
 txtfile.close();

}
