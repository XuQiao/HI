//Extremely simple plotting of some variables from the PHOBOS v1.5 Glauber Ntuple File
{
  // Setup
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);

  //Flag to save GIFs ?
  Int_t nFlagSaveGIFs = 1 ; // 0 = don't save,  1 = save

  // open ntuple from glauber
  //TFile f("Phob_Glau_PbPb_test_v15_10k.root");
  TFile f("Phob_Glau_PbPb_test_v15_1M.root");


  //-----------------------Npart----------------------
/*
  TCanvas *c1 = new TCanvas("c1","c1",0,0,500,500);
  c1->SetLogy(1);
  c1->cd();
  TH1D *hNpart = new TH1D("hNpart","Npart;Npart;Counts",500,0,500);
  nt_Pb_Pb->Project("hNpart","Npart");
  hNpart->Draw();
  if(nFlagSaveGIFs==1)  c1->SaveAs("c1_v15.gif");

  //-----------------------Ncoll----------------------

  TCanvas *c2 = new TCanvas("c2","c2",100,0,500,500);
  c2->SetLogy(1);
  c2->cd();
  TH1D *hNcoll = new TH1D("hNcoll","Ncoll;Ncoll;Counts",2500,0,2500);
  nt_Pb_Pb->Project("hNcoll","Ncoll");
  hNcoll->Draw();
  if(nFlagSaveGIFs==1)  c2->SaveAs("c2_v15.gif");

  //-----------------------b----------------------

  TCanvas *c3 = new TCanvas("c3","c3",200,0,500,500);
  //  c3->SetLogy(1);
  c3->cd();
  TH1D *hb = new TH1D("hb","Impact Parameter;b;Counts",200,0,20);
  nt_Pb_Pb->Project("hb","B");
  hb->Draw();
  if(nFlagSaveGIFs==1)  c3->SaveAs("c3_v15.gif");
*/
  
  //------------------Ncoll vs Npart (scatter)----------------------

  TCanvas *c4 = new TCanvas("c4","c4",0,50,500,500);
  c4->SetLogz(1);
  c4->cd();
  TH2D *hNpartNcoll = new TH2D("hNpartNcoll","Npart vs Ncoll;Npart;Ncoll",400,0,400,2500,0,2500);
  nt_Pb_Pb->Project("hNpartNcoll","Ncoll:Npart");
  hNpartNcoll->Draw("colz");
  if(nFlagSaveGIFs==1)  c4->SaveAs("c4_v15.gif");

  //------------------Ncoll vs Npart (profile)----------------------

  TCanvas *c5 = new TCanvas("c5","c5",100,50,500,500);
  c5->cd();
  TProfile *pNpartNcoll = new TProfile("pNpartNcoll","Profile Npart vs Ncoll;Npart;Ncoll",400,0,400);
  nt_Pb_Pb->Project("pNpartNcoll","Ncoll:Npart");
  pNpartNcoll->Draw();
  if(nFlagSaveGIFs==1)  c5->SaveAs("c5_v15.gif");


  TCanvas *c5x = new TCanvas("c5x","c5x",100,50,500,500);
  c5x->cd();
  TH1D *pNpartNcollx = new TH1D("pNpartNcollx","Profile Npart vs Ncollx;Npart;Ncollx",400,0,400);
  double np;
  double nc;
  for(int i=0;i<400;i++){
     //cout<<"i ="<<i<<"  ncoll: "<<pNpartNcoll->GetBinContent(i)<<endl;
  pNpartNcollx->SetBinContent(i,sqrt(pNpartNcoll->GetBinContent(i)));
  }
  pNpartNcollx->Draw();
  if(nFlagSaveGIFs==1)  c5x->SaveAs("c5x_v15.gif");



  //------------------EccPART vs Npart (scatter)----------------------

  TCanvas *c6 = new TCanvas("c6","c6",0,100,500,500);
  c6->SetLogz(1);
  c6->cd();
  TH2D *hNpartEccPART = new TH2D("hNpartEccPART","Npart vs EccPART;Npart;EccPART",400,0,400,1000,0,1);
  nt_Pb_Pb->Project("hNpartEccPART","EccPART:Npart");
  hNpartEccPART->Draw("colz");
  if(nFlagSaveGIFs==1)  c6->SaveAs("c6_v15.gif");

  //------------------EccPART vs Npart (profile)----------------------

  TCanvas *c7 = new TCanvas("c7","c7",100,100,500,500);
  c7->cd();
  TProfile *pNpartEccPART = new TProfile("pNpartEccPART","Profile Npart vs EccPART;Npart;EccPART",400,0,400);
  nt_Pb_Pb->Project("pNpartEccPART","EccPART:Npart");
  pNpartEccPART->Draw();
  if(nFlagSaveGIFs==1)  c7->SaveAs("c7_v15.gif");

  //------------------SPART vs Npart (scatter)----------------------
/*
  TCanvas *c8 = new TCanvas("c8","c8",0,150,500,500);
  c8->SetLogz(1);
  c8->cd();
  TH2D *hNpartSPART = new TH2D("hNpartSPART","Npart vs SPART;Npart;SPART",400,0,400,500,0,50);
  nt_Pb_Pb->Project("hNpartSPART","SPART:Npart");
  hNpartSPART->Draw("colz");
  if(nFlagSaveGIFs==1)  c8->SaveAs("c8_v15.gif");

  //------------------SPART vs Npart (profile)----------------------

  TCanvas *c9 = new TCanvas("c9","c9",100,150,500,500);
  c9->cd();
  TProfile *pNpartSPART = new TProfile("pNpartSPART","Profile Npart vs SPART;Npart;SPART",400,0,400);
  nt_Pb_Pb->Project("pNpartSPART","SPART:Npart");
  pNpartSPART->Draw();
  if(nFlagSaveGIFs==1)  c9->SaveAs("c9_v15.gif");


  //------------------S12RP_Npart vs Npart (scatter)----------------------

  TCanvas *c10 = new TCanvas("c10","c10",0,200,500,500);
  c10->SetLogz(1);
  c10->cd();
  TH2D *hNpartS12RP_Npart = new TH2D("hNpartS12RP_Npart","Npart vs S12RP_Npart;Npart;S12RP_Npart",400,0,400,2000,0,200);
  nt_Pb_Pb->Project("hNpartS12RP_Npart","S12RP_Npart:Npart");
  hNpartS12RP_Npart->Draw("colz");
  if(nFlagSaveGIFs==1)  c10->SaveAs("c10_v15.gif");

  //------------------SPART vs Npart (profile)----------------------

  TCanvas *c11 = new TCanvas("c11","c11",100,200,500,500);
  c11->cd();
  TProfile *pNpartS12RP_Npart = new TProfile("pNpartS12RP_Npart","Profile Npart vs S12RP_Npart;Npart;S12RP_Npart",400,0,400);
  nt_Pb_Pb->Project("pNpartS12RP_Npart","S12RP_Npart:Npart");
  pNpartS12RP_Npart->Draw();
  if(nFlagSaveGIFs==1)  c11->SaveAs("c11_v15.gif");

  //------------------Straight Line Fit------------------

  TF1 *fStraightLineReference = new TF1("fStraightLineReference","x",0,200);
  fStraightLineReference->SetLineColor(2);

  //------------------S12RP_Npart vs SPART4 (scatter)----------------------

  TCanvas *c12 = new TCanvas("c12","c12",0,250,500,500);
  c12->SetLogz(1);
  c12->cd();
  TH2D *hSPART4S12RP_Npart = new TH2D("hSPART4S12RP_Npart","SPART4 vs S12RP_Npart;SPART4;S12RP_Npart",2000,0,200,2000,0,200);
  nt_Pb_Pb->Project("hSPART4S12RP_Npart","S12RP_Npart:SPART4");
  hSPART4S12RP_Npart->Draw("colz");
  fStraightLineReference->Draw("same");
  if(nFlagSaveGIFs==1)  c12->SaveAs("c12_v15.gif");

  //------------------S12RP_Npart vs SPART4 (profile)----------------------

  TCanvas *c13 = new TCanvas("c13","c13",100,250,500,500);
  c13->cd();
  TProfile *pSPART4S12RP_Npart = new TProfile("pSPART4S12RP_Npart","Profile SPART4 vs S12RP_Npart;SPART4;S12RP_Npart",200,0,200);
  nt_Pb_Pb->Project("pSPART4S12RP_Npart","S12RP_Npart:SPART4");
  pSPART4S12RP_Npart->Draw();
  fStraightLineReference->Draw("same");
  if(nFlagSaveGIFs==1)  c13->SaveAs("c13_v15.gif");
  */

/*
  //------------------Ecc1PART (1D)----------------------

  TCanvas *c14a = new TCanvas("c14a","c14a",0,50,500,500);
  c14a->cd();
  TH1D *hEcc1PART = new TH1D("hEcc1PART","Ecc1PART;Ecc1PART;Counts",100,0,1);
  nt_Pb_Pb->Project("hEcc1PART","Ecc1PART");
  hEcc1PART->Draw();
  c14a->Update();
  if(nFlagSaveGIFs==1)  c14a->SaveAs("c14a_v15.gif");


  //------------------Ecc2PART (1D)----------------------

  TCanvas *c14b = new TCanvas("c14b","c14b",100,50,500,500);
  c14b->cd();
  TH1D *hEcc2PART = new TH1D("hEcc2PART","Ecc2PART;Ecc2PART;Counts",100,0,1);
  nt_Pb_Pb->Project("hEcc2PART","Ecc2PART");
  hEcc2PART->Draw();
  c14b->Update();
  if(nFlagSaveGIFs==1)  c14b->SaveAs("c14b_v15.gif");


  //------------------Ecc3PART (1D)----------------------

  TCanvas *c14c = new TCanvas("c14c","c14c",200,50,500,500);
  c14c->cd();
  TH1D *hEcc3PART = new TH1D("hEcc3PART","Ecc3PART;Ecc3PART;Counts",100,0,1);
  nt_Pb_Pb->Project("hEcc3PART","Ecc3PART");
  hEcc3PART->Draw();
  c14c->Update();
  if(nFlagSaveGIFs==1)  c14c->SaveAs("c14c_v15.gif");


  //------------------Ecc4PART (1D)----------------------

  TCanvas *c14d = new TCanvas("c14d","c14d",300,50,500,500);
  c14d->cd();
  TH1D *hEcc4PART = new TH1D("hEcc4PART","Ecc4PART;Ecc4PART;Counts",100,0,1);
  nt_Pb_Pb->Project("hEcc4PART","Ecc4PART");
  hEcc4PART->Draw();
  c14d->Update();
  if(nFlagSaveGIFs==1)  c14d->SaveAs("c14d_v15.gif");


  //------------------Ecc5PART (1D)----------------------

  TCanvas *c14e = new TCanvas("c14e","c14e",400,50,500,500);
  c14e->cd();
  TH1D *hEcc5PART = new TH1D("hEcc5PART","Ecc5PART;Ecc5PART;Counts",100,0,1);
  nt_Pb_Pb->Project("hEcc5PART","Ecc5PART");
  hEcc5PART->Draw();
  c14e->Update();
  if(nFlagSaveGIFs==1)  c14e->SaveAs("c14e_v15.gif");


  //------------------Ecc6PART (1D)----------------------

  TCanvas *c14f = new TCanvas("c14f","c14f",500,50,500,500);
  c14f->cd();
  TH1D *hEcc6PART = new TH1D("hEcc6PART","Ecc6PART;Ecc6PART;Counts",100,0,1);
  nt_Pb_Pb->Project("hEcc6PART","Ecc6PART");
  hEcc6PART->Draw();
  c14f->Update();
  if(nFlagSaveGIFs==1)  c14f->SaveAs("c14f_v15.gif");




  //------------------Psi1PART (1D)----------------------

  TCanvas *c15a = new TCanvas("c15a","c15a",0,150,500,500);
  c15a->cd();
  TH1D *hPsi1PART = new TH1D("hPsi1PART","Psi1PART;Psi1PART;Counts",100.*TMath::Pi(),0,2.0*TMath::Pi());
  nt_Pb_Pb->Project("hPsi1PART","Psi1PART");
  hPsi1PART->Draw();
  c15a->Update();
  if(nFlagSaveGIFs==1)  c15a->SaveAs("c15a_v15.gif");


  //------------------Psi2PART (1D)----------------------

  TCanvas *c15b = new TCanvas("c15b","c15b",100,150,500,500);
  c15b->cd();
  TH1D *hPsi2PART = new TH1D("hPsi2PART","Psi2PART;Psi2PART;Counts",100.*TMath::Pi(),0,2.0*TMath::Pi());
  nt_Pb_Pb->Project("hPsi2PART","Psi2PART");
  hPsi2PART->Draw();
  c15b->Update();
  if(nFlagSaveGIFs==1)  c15b->SaveAs("c15b_v15.gif");


  //------------------Psi3PART (1D)----------------------

  TCanvas *c15c = new TCanvas("c15c","c15c",200,150,500,500);
  c15c->cd();
  TH1D *hPsi3PART = new TH1D("hPsi3PART","Psi3PART;Psi3PART;Counts",100.*TMath::Pi(),0,2.0*TMath::Pi());
  nt_Pb_Pb->Project("hPsi3PART","Psi3PART");
  hPsi3PART->Draw();
  c15c->Update();
  if(nFlagSaveGIFs==1)  c15c->SaveAs("c15c_v15.gif");


  //------------------Psi4PART (1D)----------------------

  TCanvas *c15d = new TCanvas("c15d","c15d",300,150,500,500);
  c15d->cd();
  TH1D *hPsi4PART = new TH1D("hPsi4PART","Psi4PART;Psi4PART;Counts",100.*TMath::Pi(),0,2.0*TMath::Pi());
  nt_Pb_Pb->Project("hPsi4PART","Psi4PART");
  hPsi4PART->Draw();
  c15d->Update();
  if(nFlagSaveGIFs==1)  c15d->SaveAs("c15d_v15.gif");


  //------------------Psi5PART (1D)----------------------

  TCanvas *c15e = new TCanvas("c15e","c15e",400,150,500,500);
  c15e->cd();
  TH1D *hPsi5PART = new TH1D("hPsi5PART","Psi5PART;Psi5PART;Counts",100.*TMath::Pi(),0,2.0*TMath::Pi());
  nt_Pb_Pb->Project("hPsi5PART","Psi5PART");
  hPsi5PART->Draw();
  c15e->Update();
  if(nFlagSaveGIFs==1)  c15e->SaveAs("c15e_v15.gif");


  //------------------Psi6PART (1D)----------------------

  TCanvas *c15f = new TCanvas("c15f","c15f",500,150,500,500);
  c15f->cd();
  TH1D *hPsi6PART = new TH1D("hPsi6PART","Psi6PART;Psi6PART;Counts",100.*TMath::Pi(),0,2.0*TMath::Pi());
  nt_Pb_Pb->Project("hPsi6PART","Psi6PART");
  hPsi6PART->Draw();
  c15f->Update();
  if(nFlagSaveGIFs==1)  c15f->SaveAs("c15f_v15.gif");

*/

  /*

  //------------------Ecc1PART vs Npart (2D and profile)----------------------

  TCanvas *c16a = new TCanvas("c16a","c16a",0,300,500,500);
  c16a->SetLogz(1);
  c16a->cd();
  TH2D *hNpart_Ecc1PART = new TH2D("hNpart_Ecc1PART","Ecc1PART vs Npart;Npart;Ecc1PART",500,0,500,100,0,1);
  nt_Pb_Pb->Project("hNpart_Ecc1PART","Ecc1PART:Npart","Npart>2");
  hNpart_Ecc1PART->Draw("colz");
  TProfile *pNpartEcc1PART = new TProfile("pNpartEcc1PART","Profile Npart vs Ecc1PART;Npart;Ecc1PART",400,0,400);
  nt_Pb_Pb->Project("pNpartEcc1PART","Ecc1PART:Npart");
  pNpartEcc1PART->Draw("same");
  c16a->Update();
  if(nFlagSaveGIFs==1)  c16a->SaveAs("c16a_v15.gif");


  //------------------Ecc2PART vs Npart (2D and profile)----------------------

  TCanvas *c16b = new TCanvas("c16b","c16b",100,300,500,500);
  c16b->SetLogz(1);
  c16b->cd();
  TH2D *hNpart_Ecc2PART = new TH2D("hNpart_Ecc2PART","Ecc2PART vs Npart;Npart;Ecc2PART",500,0,500,100,0,1);
  nt_Pb_Pb->Project("hNpart_Ecc2PART","Ecc2PART:Npart","Npart>2");
  hNpart_Ecc2PART->Draw("colz");
  TProfile *pNpartEcc2PART = new TProfile("pNpartEcc2PART","Profile Npart vs Ecc2PART;Npart;Ecc2PART",400,0,400);
  nt_Pb_Pb->Project("pNpartEcc2PART","Ecc2PART:Npart");
  pNpartEcc2PART->Draw("same");
  c16b->Update();
  if(nFlagSaveGIFs==1)  c16b->SaveAs("c16b_v15.gif");


  //------------------Ecc3PART vs Npart (2D and profile)----------------------

  TCanvas *c16c = new TCanvas("c16c","c16c",200,300,500,500);
  c16c->SetLogz(1);
  c16c->cd();
  TH2D *hNpart_Ecc3PART = new TH2D("hNpart_Ecc3PART","Ecc3PART vs Npart;Npart;Ecc3PART",500,0,500,100,0,1);
  nt_Pb_Pb->Project("hNpart_Ecc3PART","Ecc3PART:Npart","Npart>2");
  hNpart_Ecc3PART->Draw("colz");
  TProfile *pNpartEcc3PART = new TProfile("pNpartEcc3PART","Profile Npart vs Ecc3PART;Npart;Ecc3PART",400,0,400);
  nt_Pb_Pb->Project("pNpartEcc3PART","Ecc3PART:Npart");
  pNpartEcc3PART->Draw("same");
  c16c->Update();
  if(nFlagSaveGIFs==1)  c16c->SaveAs("c16c_v15.gif");


  //------------------Ecc4PART vs Npart (2D and profile)----------------------

  TCanvas *c16d = new TCanvas("c16d","c16d",300,300,500,500);
  c16d->SetLogz(1);
  c16d->cd();
  TH2D *hNpart_Ecc4PART = new TH2D("hNpart_Ecc4PART","Ecc4PART vs Npart;Npart;Ecc4PART",500,0,500,100,0,1);
  nt_Pb_Pb->Project("hNpart_Ecc4PART","Ecc4PART:Npart","Npart>2");
  hNpart_Ecc4PART->Draw("colz");
  TProfile *pNpartEcc4PART = new TProfile("pNpartEcc4PART","Profile Npart vs Ecc4PART;Npart;Ecc4PART",400,0,400);
  nt_Pb_Pb->Project("pNpartEcc4PART","Ecc4PART:Npart");
  pNpartEcc4PART->Draw("same");
  c16d->Update();
  if(nFlagSaveGIFs==1)  c16d->SaveAs("c16d_v15.gif");


  //------------------Ecc5PART vs Npart (2D and profile)----------------------

  TCanvas *c16e = new TCanvas("c16e","c16e",400,300,500,500);
  c16e->SetLogz(1);
  c16e->cd();
  TH2D *hNpart_Ecc5PART = new TH2D("hNpart_Ecc5PART","Ecc5PART vs Npart;Npart;Ecc5PART",500,0,500,100,0,1);
  nt_Pb_Pb->Project("hNpart_Ecc5PART","Ecc5PART:Npart","Npart>2");
  hNpart_Ecc5PART->Draw("colz");
  TProfile *pNpartEcc5PART = new TProfile("pNpartEcc5PART","Profile Npart vs Ecc5PART;Npart;Ecc5PART",400,0,400);
  nt_Pb_Pb->Project("pNpartEcc5PART","Ecc5PART:Npart");
  pNpartEcc5PART->Draw("same");
  c16e->Update();
  if(nFlagSaveGIFs==1)  c16e->SaveAs("c16e_v15.gif");


  //------------------Ecc6PART vs Npart (2D and profile)----------------------

  TCanvas *c16f = new TCanvas("c16f","c16f",500,300,500,500);
  c16f->SetLogz(1);
  c16f->cd();
  TH2D *hNpart_Ecc6PART = new TH2D("hNpart_Ecc6PART","Ecc6PART vs Npart;Npart;Ecc6PART",500,0,500,100,0,1);
  nt_Pb_Pb->Project("hNpart_Ecc6PART","Ecc6PART:Npart","Npart>2");
  hNpart_Ecc6PART->Draw("colz");
  TProfile *pNpartEcc6PART = new TProfile("pNpartEcc6PART","Profile Npart vs Ecc6PART;Npart;Ecc6PART",400,0,400);
  nt_Pb_Pb->Project("pNpartEcc6PART","Ecc6PART:Npart");
  pNpartEcc6PART->Draw("same");
  c16f->Update();
  if(nFlagSaveGIFs==1)  c16f->SaveAs("c16f_v15.gif");



  //------------------Psi1PART (2D)----------------------

  TCanvas *c17a = new TCanvas("c17a","c17a",0,450,500,500);
  c17a->cd();
  c17a->SetLogz(1);
  TH2D *hNpart_Psi1PART = new TH2D("hNpart_Psi1PART","Psi1PART vs Npart;Npart;Psi1PART",500,0,500,100.*TMath::Pi(),0,2.0*TMath::Pi());
  nt_Pb_Pb->Project("hNpart_Psi1PART","Psi1PART:Npart","Npart>2");
  hNpart_Psi1PART->Draw("colz");
  c17a->Update();
  if(nFlagSaveGIFs==1)  c17a->SaveAs("c17a_v15.gif");


  //------------------Psi2PART (2D)----------------------

  TCanvas *c17b = new TCanvas("c17b","c17b",100,450,500,500);
  c17b->cd();
  c17b->SetLogz(1);
  TH2D *hNpart_Psi2PART = new TH2D("hNpart_Psi2PART","Psi2PART vs Npart;Npart;Psi2PART",500,0,500,100.*TMath::Pi(),0,2.0*TMath::Pi());
  nt_Pb_Pb->Project("hNpart_Psi2PART","Psi2PART:Npart","Npart>2");
  hNpart_Psi2PART->Draw("colz");
  c17b->Update();
  if(nFlagSaveGIFs==1)  c17b->SaveAs("c17b_v15.gif");


  //------------------Psi3PART (2D)----------------------

  TCanvas *c17c = new TCanvas("c17c","c17c",200,450,500,500);
  c17c->SetLogz(1);
  c17c->cd();
  TH2D *hNpart_Psi3PART = new TH2D("hNpart_Psi3PART","Psi3PART vs Npart;Npart;Psi3PART",500,0,500,100.*TMath::Pi(),0,2.0*TMath::Pi());
  nt_Pb_Pb->Project("hNpart_Psi3PART","Psi3PART:Npart","Npart>2");
  hNpart_Psi3PART->Draw("colz");
  c17c->Update();
  if(nFlagSaveGIFs==1)  c17c->SaveAs("c17c_v15.gif");


  //------------------Psi4PART (2D)----------------------

  TCanvas *c17d = new TCanvas("c17d","c17d",300,450,500,500);
  c17d->SetLogz(1);
  c17d->cd();
  TH2D *hNpart_Psi4PART = new TH2D("hNpart_Psi4PART","Psi4PART vs Npart;Npart;Psi4PART",500,0,500,100.*TMath::Pi(),0,2.0*TMath::Pi());
  nt_Pb_Pb->Project("hNpart_Psi4PART","Psi4PART:Npart","Npart>2");
  hNpart_Psi4PART->Draw("colz");
  c17d->Update();
  if(nFlagSaveGIFs==1)  c17d->SaveAs("c17d_v15.gif");

  //------------------Psi5PART (2D)----------------------

  TCanvas *c17e = new TCanvas("c17e","c17e",400,450,500,500);
  c17e->SetLogz(1);
  c17e->cd();
  TH2D *hNpart_Psi5PART = new TH2D("hNpart_Psi5PART","Psi5PART vs Npart;Npart;Psi5PART",500,0,500,100.*TMath::Pi(),0,2.0*TMath::Pi());
  nt_Pb_Pb->Project("hNpart_Psi5PART","Psi5PART:Npart","Npart>2");
  hNpart_Psi5PART->Draw("colz");
  c17e->Update();
  if(nFlagSaveGIFs==1)  c17e->SaveAs("c17e_v15.gif");

  //------------------Psi6PART (2D)----------------------

  TCanvas *c17f = new TCanvas("c17f","c17f",500,450,500,500);
  c17f->SetLogz(1);
  c17f->cd();
  TH2D *hNpart_Psi6PART = new TH2D("hNpart_Psi6PART","Psi6PART vs Npart;Npart;Psi6PART",500,0,500,100.*TMath::Pi(),0,2.0*TMath::Pi());
  nt_Pb_Pb->Project("hNpart_Psi6PART","Psi6PART:Npart","Npart>2");
  hNpart_Psi6PART->Draw("colz");
  c17f->Update();
  if(nFlagSaveGIFs==1)  c17f->SaveAs("c17f_v15.gif");
  */


  //------------------EccCum2 vs Npart (scatter)----------------------

  TCanvas *c34 = new TCanvas("c34","c34",0,100,500,500);
  c34->SetLogz(1);
  c34->cd();
  TH2D *hNpartEccCum2 = new TH2D("hNpartEccCum2","Npart vs EccCum2;Npart;EccCum2",400,0,400,1000,0,1);
  nt_Pb_Pb->Project("hNpartEccCum2","EccCum2:Npart");
  hNpartEccCum2->Draw("colz");
  if(nFlagSaveGIFs==1)  c34->SaveAs("c34_v15.gif");

  //------------------EccCum4 vs Npart (scatter)----------------------

  TCanvas *c35 = new TCanvas("c35","c35",0,100,500,500);
  c35->SetLogz(1);
  c35->cd();
  TH2D *hNpartEccCum4 = new TH2D("hNpartEccCum4","Npart vs EccCum4;Npart;EccCum4",400,0,400,1000,0,1);
  nt_Pb_Pb->Project("hNpartEccCum4","EccCum4:Npart");
  hNpartEccCum4->Draw("colz");
  if(nFlagSaveGIFs==1)  c34->SaveAs("c34_v15.gif");

  //------------------EccCum vs Npart (profile)----------------------
  TCanvas *c34b = new TCanvas("c34b","c34b",100,100,610,500);
  c34b->cd();
  TProfile *pNpartEccCum2 = new TProfile("pNpartEccCum2","Profile Npart vs EccCum2;Npart;EccCum2",400,0,400);
  nt_Pb_Pb->Project("pNpartEccCum2","EccCum2:Npart");
  //pNpartEccCum2->Draw();
  TH1D *pNpartEccCum2b = new TH1D("pNpartEcc3","Profile N_{part} vs #epsilon;N_{part};#epsilon",400,0,400);
const int ReBins = 20;
  for(int i=0;i<400;i++){
      //if(i<200)cout<<pNpartEccCum2->GetBinContent(i)<<endl;
      //pNpartEccCum2b->SetBinContent(i,sqrt(pNpartEccCum2->GetBinContent(i))/ReBins);
      pNpartEccCum2b->SetBinContent(i,sqrt(pNpartEccCum2->GetBinContent(i)));
      //pNpartEccCum2b->SetBinError(i,sqrt(pNpartEccCum2->GetBinError(i)));
    }

  TProfile *pNpartEccCum4 = new TProfile("pNpartEccCum4","Profile Npart vs EccCum4;Npart;EccCum4",400,0,400);
  nt_Pb_Pb->Project("pNpartEccCum4","EccCum4:Npart");
  pNpartEccCum4->Draw("same");
  TH1D *pNpartEccCum4b = new TH1D("pNpartEccCum4b","Profile Npart vs EccCum4;Npart;EccCum4",400,0,400);
  double temp;
  for(int i=0;i<400;i++){
      temp=pNpartEccCum2->GetBinContent(i);
      if(2*temp*temp-pNpartEccCum4->GetBinContent(i)>0)
      //pNpartEccCum4b->SetBinContent(i,pow(2*temp*temp-pNpartEccCum4->GetBinContent(i),1.0/4)/ReBins);
      pNpartEccCum4b->SetBinContent(i,pow(2*temp*temp-pNpartEccCum4->GetBinContent(i),1.0/4));
      //pNpartEccCum2b->SetBinError(i,sqrt(pNpartEccCum2->GetBinError(i)));
      //cout<<pow(2*temp*temp-pNpartEccCum4->GetBinContent(i),1.0/4)<<endl;
    }
//////////
  TH1D *pNpartEccPARTb = new TH1D("pNpartEccPARTb","Profile Npart vs EccPART;Npart;EccPART",400,0,400);
  for(int i=0;i<400;i++){
      //pNpartEccPARTb->SetBinContent(i,pNpartEccPART->GetBinContent(i)/ReBins);
      pNpartEccPARTb->SetBinContent(i,pNpartEccPART->GetBinContent(i));
    }
  pNpartEccCum2b->SetStats(0);
  //pNpartEccCum2b->GetYaxis()>CenterTitle();
  //pNpartEccCum2b->GetXaxis()>CenterTitle();
  //pNpartEccCum2b->Rebin(ReBins);
  //pNpartEccPARTb->Rebin(ReBins);
  //pNpartEccCum4b->Rebin(ReBins);
  pNpartEccCum2b->Draw("L");
  pNpartEccPARTb->SetLineColor(2);
  pNpartEccPARTb->Draw("Lsame"); 
  pNpartEccCum4b->SetLineColor(4);
  pNpartEccCum4b->Draw("Lsame");

    TLegend *leg = new TLegend(0.6,0.6,0.88,0.86);
    leg->SetFillColor(10);
    leg->SetBorderSize(0.035);
    leg->SetTextFont(42);
    //leg->SetTextColor(2);
    leg->SetTextSize(0.055);
    //leg->SetLineStyle(0.06);
    leg->AddEntry(pNpartEccCum2b,"#epsilon_{part}{2}","l");
    leg->AddEntry(pNpartEccPARTb,"#epsilon_{part}","l");
    leg->AddEntry(pNpartEccCum4b,"#epsilon_{part}{4}","l");
    leg->Draw();

  if(nFlagSaveGIFs==1)  c34b->SaveAs("c34b_v15.gif");
}
