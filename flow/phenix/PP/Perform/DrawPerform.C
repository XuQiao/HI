#include "/phenix/u/xuq/util/SimplifyLife.C"
void DrawPerform(){
gStyle->SetErrorX(0);
 TFile *fmb = new TFile("merged.root","ReadOnly");
 TFile *fhm = new TFile("merged_FVtx.root","ReadOnly");
const int ncav = 26;
TCanvas *c1[ncav];
for(int i=0;i<ncav;i++){
c1[i] = new TCanvas();
}

//Tracks
  TH2F *hcntetaphi = (TH2F*)fmb->Get("hcntetaphi");
c1[0]->cd();
hcntetaphi->Draw("colz");
c1[0]->Print("fig/hcntetaphi.png");

  TH1F *hcntpt = (TH1F*)fmb->Get("hcntpt");
  TH1F *hcntpt_hm = (TH1F*)fhm->Get("hcntpt");
c1[1]->cd();
c1[1]->SetLogy();
hcntpt->Rebin(4);
hcntpt_hm->Rebin(4);
SetStyle(*hcntpt,1.2,1,20,0,0);
SetStyle(*hcntpt_hm,1.2,2,24,0,0);
hcntpt->Draw("P");
hcntpt_hm->Draw("Psame");
c1[1]->Print("fig/hcntpt.png");

  TH2F* pc3dphidz = (TH2F*)fmb->Get("pc3dphidz");
TH1F* pc3dphi = (TH1F*)pc3dphidz->ProjectionX("pc3dphi",0,-1); 
TH1F* pc3dz = (TH1F*)pc3dphidz->ProjectionY("pc3dz",0,-1); 
c1[2]->cd();
SetStyle(*pc3dphi,1.2,1,20,0,0);
pc3dphi->Draw("P");
c1[2]->Print("fig/pc3dphi.png");

c1[3]->cd();
SetStyle(*pc3dz,1.2,1,20,0,0);
pc3dz->Rebin(4);
pc3dz->Draw("P");
c1[3]->Print("fig/pc3dz.png");
  
TH2F* hntracknmpc = (TH2F*)fmb->Get("hntracknmpc");
TH2F* hntracknmpc_hm = (TH2F*)fhm->Get("hntracknmpc");
TH1F* hntrack = (TH1F*)hntracknmpc->ProjectionX("hntrack",0,-1);
TH1F* hntrack_hm = (TH1F*)hntracknmpc_hm->ProjectionX("hntrack_hm",0,-1);
c1[4]->cd();
c1[4]->SetLogy();
hntrack->Scale(1./hntrack->Integral());
hntrack_hm->Scale(1./hntrack_hm->Integral());
SetStyle(*hntrack,1.2,1,20,0,0);
SetStyle(*hntrack_hm,1.2,2,24,0,0);
SetRange(*hntrack,0,1e-11,20,10);
hntrack->Draw();
hntrack_hm->Draw("Psame");
c1[4]->Print("fig/ntrack.png");
/*
//tof
  TH2F* tofdphidz = (TH2F*)f->Get("tofdphidz");
TH1F* tofdphi = (TH1F*)tofdphidz->ProjectionX("tofdphi",0,-1); 
TH1F* tofdz = (TH1F*)tofdphidz->ProjectionY("tofdz",0,-1); 
c1[4]->cd();
tofdphi->Draw();
c1[4]->Print("fig/tofdphi.png");

c1[5]->cd();
tofdz->Draw();
c1[5]->Print("fig/tofdz.png");

  TH2F* tofwdphidz = (TH2F*)f->Get("tofwdphidz");
  TH2F* ttofqpratio = (TH2F*)f->Get("ttofqpratio");
c1[6]->cd();
ttofqpratio->Draw("colz");
c1[6]->Print("fig/ttofqpratio.png");
  TH2F* m2qpratio = (TH2F*)f->Get("m2qpratio");


//vtx
  TH2F* hcluetaphi[4];
   c1[7]->Divide(2,2);
for(int i=0;i<4;i++){
   hcluetaphi[i] = (TH2F*)f->Get(Form("hcluetaphi_%d",i));
   c1[7]->cd(i+1);
   hcluetaphi[i]->Draw("colz");
   TLatex t;
   t.SetNDC();
   t.DrawLatex(0.1,0.8,Form("layer %d",i));
}
   c1[7]->Print("fig/hcluetaphi.png");
//bbc
  TH1F* bbcet = (TH1F*)f->Get("bbcet");
  c1[8]->cd();
  c1[8]->SetLogy();
  bbcet->Draw();
  c1[8]->Print("fig/bbcet.png");
//fvtx
  TH2D *DCAxydis[2];
  TH2D *DCAxy2dis[2];
  TH2D *DCAcentdis[2];
  c1[9]->Divide(1,2);
  c1[10]->Divide(1,2);
  c1[11]->Divide(1,2);
for(int iarm=0;iarm<2;iarm++){
  DCAxydis[iarm] = (TH2D*)f->Get(Form("DCAxydis_%d",iarm));
  DCAxy2dis[iarm] = (TH2D*)f->Get(Form("DCAxy2dis_%d",iarm));
  DCAcentdis[iarm] = (TH2D*)f->Get(Form("DCAcentdis_%d",iarm));
  c1[9]->cd(iarm+1);
  c1[9]->SetLogz();
  DCAxydis[iarm]->Draw("colz");
  c1[10]->cd(iarm+1);
  c1[10]->SetLogz();
  DCAxy2dis[iarm]->Draw("colz");
  c1[11]->cd(iarm+1);
  c1[11]->SetLogz();
  DCAcentdis[iarm]->Draw("colz");
}
  c1[9]->Print("fig/DCAxydis.png");
  c1[10]->Print("fig/DCAxy2dis.png");
  c1[11]->Print("fig/DCAcentdis.png");
 // TH1F *fvtxdphidis;
 // TH1F *fvtxdphidis2;

c1[12]->Divide(1,2);
  TH2D *hvtx0etaz = (TH2D*) f->Get("hvtx0etaz");
  TH2D *hvtx1etaz = (TH2D*) f->Get("hvtx1etaz");
c1[12]->cd(1);
hvtx0etaz->Draw("colz");
c1[12]->cd(2);
hvtx1etaz->Draw("colz");
c1[12]->Print("fig/hvtxetaz.png");

c1[13]->Divide(1,2);
  TH2D *hvtx0etaphi = (TH2D*)f->Get("hvtx0etaphi");;
  TH2D *hvtx1etaphi = (TH2D*)f->Get("hvtx1etaphi");;
c1[13]->cd(1);
hvtx0etaphi->Draw("colz");
c1[13]->cd(2);
hvtx1etaphi->Draw("colz");
c1[13]->Print("fig/hvtxetaphi.png");

//mpc
 
TH2F* mpcetdis;
TH2F* mpcetetasouth[1];
TH2F* mpcetetanorth[1];
TH2F*  mpc_south_cent;
TH2F*  mpc_north_cent;
TH2F*  mpc_south_north;
TH2F *south_mpc_north_mpc;


//correlate
  TH2F* hvtxzfvtxz = (TH2F*)f->Get("hvtxzfvtxz");
c1[14]->cd();
c1[14]->SetLogz();
hvtxzfvtxz->Draw("colz");
c1[14]->Print("fig/hvtxzfvtxz.png");
  TH2F* hpc1hitsbbc = (TH2F*)f->Get("hpc1hitsbbc");
   TH1F* hpc1hits = (TH1F*)hpc1hitsbbc->ProjectionX("hpc1hits",0,-1);
c1[15]->cd();
c1[15]->SetLogy();
hpc1hits->Draw();
c1[15]->Print("fig/hpc1hits.png");
   TH1F* hbbc = (TH1F*)hpc1hitsbbc->ProjectionY("hbbc",0,-1);
c1[16]->cd();
c1[16]->SetLogy();
hbbc->Draw();
c1[16]->Print("fig/hbbc.png");
  TH2F* hnpc3hitsntof = (TH2F*)f->Get("hnpc3hitsntof");
   TH1F* hnpc3hits = (TH1F*)hnpc3hitsntof->ProjectionX("hnpc3hits",0,-1);
c1[17]->cd();
c1[17]->SetLogy();
hnpc3hits->Draw();
c1[17]->Print("fig/hnpc3hits.png");
   TH1F* hntof = (TH1F*)hnpc3hitsntof->ProjectionY("hntof",0,-1);
c1[18]->cd();
c1[18]->SetLogy();
hntof->Draw();
c1[18]->Print("fig/hntof.png");
  TH2F* hbbcnbbc = (TH2F*)f->Get("hbbcnbbc");
c1[19]->cd();
c1[19]->SetLogz();
hbbcnbbc->Draw("colz");
c1[19]->Print("fig/hbbcnbbc.png");

  TH2F* hbbcsbbcn = (TH2F*)f->Get("hbbcsbbcn");
c1[20]->cd();
c1[20]->SetLogz();
hbbcsbbcn->Draw("colz");
c1[20]->Print("fig/hbbcsbbcn.png");
*/
  TH2F* hnvtxnfvtxtrk[4];
  TH2F* hbbcsnvtx[4];
  TH2F* hbbcnnvtx[4];
  TH2F* hbbcnvtx[4];
  TH2F* hbbcnvtx_hm[4];
c1[21]->Divide(2,2);
for(int i=0;i<4;i++){
  hnvtxnfvtxtrk[i] = (TH2F*)fmb->Get(Form("hnvtxnfvtxtrk_%d",i));
  hbbcsnvtx[i] = (TH2F*)fmb->Get(Form("hbbcsnvtx_%d",i));
  hbbcnnvtx[i] = (TH2F*)fmb->Get(Form("hbbcnvtx_%d",i));
  hbbcnvtx[i] = (TH2F*)fmb->Get(Form("hbbcnvtx_%d",i));
  hbbcnvtx_hm[i] = (TH2F*)fhm->Get(Form("hbbcnvtx_%d",i));
c1[21]->cd(i+1);
c1[21]->SetLogz();
hnvtxnfvtxtrk[i]->Draw("colz");
TLatex t;
t.SetNDC();
t.DrawLatex(0.1,0.8,Form("layer %d",i));
}
c1[21]->Print("fig/hnvtxnfvtxtrk.png");

TH1F* hnvtx = (TH1F*)hbbcnvtx[0]->ProjectionY("hnvtx",0,-1);
TH1F* hnvtx_hm = (TH1F*)hbbcnvtx_hm[0]->ProjectionY("hnvtx_hm",0,-1);
TH1F* hbbc = (TH1F*)hbbcnvtx[0]->ProjectionX("hbbc",0,-1);
TH1F* hbbc_hm = (TH1F*)hbbcnvtx_hm[0]->ProjectionX("hbbc_hm",0,-1);
c1[22]->cd();
c1[22]->SetLogy();
hnvtx->Rebin(4);
hnvtx_hm->Rebin(4);
hnvtx->Scale(1./hnvtx->Integral());
hnvtx_hm->Scale(1./hnvtx_hm->Integral());
SetRange(*hnvtx,0,1e-11,200,10);
SetTitle(*hnvtx,"#hits in cluster layer 1","normalized","");
SetStyle(*hnvtx,1.2,1,20,0,0);
SetStyle(*hnvtx_hm,1.2,2,24,0,0);
hnvtx->Draw("P");
hnvtx_hm->Draw("Psame");
c1[22]->Print("fig/hnvtx.png");
c1[23]->cd();
c1[23]->SetLogy();
SetTitle(*hbbc,"bbc charge sum","normalized","");
hbbc->Rebin(5);
hbbc_hm->Rebin(5);
hbbc->Scale(1./hbbc->Integral());
hbbc_hm->Scale(1./hbbc_hm->Integral());
SetRange(*hbbc,0,1e-11,200,10);
SetStyle(*hbbc,1.2,1,20,0,0);
SetStyle(*hbbc_hm,1.2,2,24,0,0);
hbbc->Draw("P");
hbbc_hm->Draw("Psame");
c1[23]->Print("fig/hbbc.png");
c1[24]->cd();
c1[24]->SetLogz();
SetTitle(*hbbcnvtx_hm[0],"bbc charge sum","#hits in cluster layer 1","");
hbbcnvtx_hm[0]->Draw("colz");
c1[24]->Print("fig/hbbcnvtx_hm.png");


/*
  TH2F* hnbbcnclu = (TH2F*)f->Get("hnbbcnclu");
TH1F* hnbbc = (TH1F*)hnbbcnclu->ProjectionX("hnbbc",0,-1);
TH1F* hnclu = (TH1F*)hnbbcnclu->ProjectionY("hnclu",0,-1);
c1[23]->cd();
c1[23]->SetLogy();
hnbbc->Draw();
c1[23]->Print("fig/hnbbc.png");
c1[24]->cd();
c1[24]->SetLogy();
hnclu->Draw();
c1[24]->Print("fig/hnclu.png");

  TH2F *hnfvtxtrkbbc = (TH2F*)f->Get("hnfvtxtrkbbc");
c1[25]->cd();
c1[25]->SetLogz();
hnfvtxtrkbbc->Draw("colz");
c1[25]->Print("fig/hnfvtxtrkbbc.png");

  TH2F *hnfvtxtrksnmpcs = (TH2F*)f->Get("hnfvtxtrksnmpcs");
  TH2F *hnfvtxtrknnmpcn = (TH2F*)f->Get("hnfvtxtrknnmpcn");
  TH2F *south_mpc_south_bbc;
  TH2F *north_mpc_north_bbc;
*/  


}
