{
#include </afs/cern.ch/user/k/kaulmer/tdrstyle.C>

/*
Jobs to run to do systematics study:

Qsub -l lnxfarm_quad4 -e -N doEff9_nom -o doEff_systVar_900GeVMC_nominal.log doEff systVar_900GeVMC_nominal 11 900GeV 3 3 3 3 4 0.01

Qsub -l lnxfarm_quad4 -e -N doEff9_ip2 -o doEff_systVar_900GeVMC_ip2.log doEff systVar_900GeVMC_ip2 11 900GeV 3 2 2 2 4 0.01
Qsub -l lnxfarm_quad4 -e -N doEff9_ip4 -o doEff_systVar_900GeVMC_ip4.log doEff systVar_900GeVMC_ip4 11 900GeV 3 4 4 4 4 0.01

Qsub -l lnxfarm_quad4 -e -N doEff9_CL001 -o doEff_systVar_900GeVMC_CL001.log doEff systVar_900GeVMC_CL001 11 900GeV 3 3 3 3 4 0.001
Qsub -l lnxfarm_quad4 -e -N doEff9_CL02 -o doEff_systVar_900GeVMC_CL02.log doEff systVar_900GeVMC_CL02 11 900GeV 3 3 3 3 4 0.02

Qsub -l lnxfarm_quad4 -e -N doEff9_sepXi3 -o doEff_systVar_900GeVMC_sepXi3.log doEff systVar_900GeVMC_sepXi3 11 900GeV 3 3 3 3 3 0.01
Qsub -l lnxfarm_quad4 -e -N doEff9_sepXi5 -o doEff_systVar_900GeVMC_sepXi5.log doEff systVar_900GeVMC_sepXi5 11 900GeV 3 3 3 3 5 0.01

Qsub -l lnxfarm_quad4 -e -N doEff9_ipXi2 -o doEff_systVar_900GeVMC_ipXi2.log doEff systVar_900GeVMC_ipXi2 11 900GeV 2 3 3 3 4 0.01
Qsub -l lnxfarm_quad4 -e -N doEff9_ipXi4 -o doEff_systVar_900GeVMC_ipXi4.log doEff systVar_900GeVMC_ipXi4 11 900GeV 4 3 3 3 4 0.01


Qsub -l lnxfarm_quad4 -e -N doEff7_nom -o doEff_systVar_7TeVMC_nominal.log doEff systVar_7TeVMC_nominal 1 7TeV 3 3 3 3 4 0.01

Qsub -l lnxfarm_quad4 -e -N doEff7_ip2 -o doEff_systVar_7TeVMC_ip2.log doEff systVar_7TeVMC_ip2 1 7TeV 3 2 2 2 4 0.01
Qsub -l lnxfarm_quad4 -e -N doEff7_ip4 -o doEff_systVar_7TeVMC_ip4.log doEff systVar_7TeVMC_ip4 1 7TeV 3 4 4 4 4 0.01

Qsub -l lnxfarm_quad4 -e -N doEff7_CL001 -o doEff_systVar_7TeVMC_CL001.log doEff systVar_7TeVMC_CL001 1 7TeV 3 3 3 3 4 0.001
Qsub -l lnxfarm_quad4 -e -N doEff7_CL02 -o doEff_systVar_7TeVMC_CL02.log doEff systVar_7TeVMC_CL02 1 7TeV 3 3 3 3 4 0.02

Qsub -l lnxfarm_quad4 -e -N doEff7_sepXi3 -o doEff_systVar_7TeVMC_sepXi3.log doEff systVar_7TeVMC_sepXi3 1 7TeV 3 3 3 3 3 0.01
Qsub -l lnxfarm_quad4 -e -N doEff7_sepXi5 -o doEff_systVar_7TeVMC_sepXi5.log doEff systVar_7TeVMC_sepXi5 1 7TeV 3 3 3 3 5 0.01

Qsub -l lnxfarm_quad4 -e -N doEff7_ipXi2 -o doEff_systVar_7TeVMC_ipXi2.log doEff systVar_7TeVMC_ipXi2 1 7TeV 2 3 3 3 4 0.01
Qsub -l lnxfarm_quad4 -e -N doEff7_ipXi4 -o doEff_systVar_7TeVMC_ipXi4.log doEff systVar_7TeVMC_ipXi4 1 7TeV 4 3 3 3 4 0.01





Qsub -l lnxfarm_quad4 -e -N doAn9_nom -o doAn_systVar_900GeVData_nominal.log doAn systVar_900GeVData_nominal systVar_900GeVMC_nominal 900GeV 3 3 3 3 4 0.01

Qsub -l lnxfarm_quad4 -e -N doAn9_ip2 -o doAn_systVar_900GeVData_ip2.log doAn systVar_900GeVData_ip2 systVar_900GeVMC_ip2 900GeV 3 2 2 2 4 0.01
Qsub -l lnxfarm_quad4 -e -N doAn9_ip4 -o doAn_systVar_900GeVData_ip4.log doAn systVar_900GeVData_ip4 systVar_900GeVMC_ip4 900GeV 3 4 4 4 4 0.01

Qsub -l lnxfarm_quad4 -e -N doAn9_CL001 -o doAn_systVar_900GeVData_CL001.log doAn systVar_900GeVData_CL001 systVar_900GeVMC_CL001 900GeV 3 3 3 3 4 0.001
Qsub -l lnxfarm_quad4 -e -N doAn9_CL02 -o doAn_systVar_900GeVData_CL02.log doAn systVar_900GeVData_CL02 systVar_900GeVMC_CL02 900GeV 3 3 3 3 4 0.02

Qsub -l lnxfarm_quad4 -e -N doAn9_sepXi3 -o doAn_systVar_900GeVData_sepXi3.log doAn systVar_900GeVData_sepXi3 systVar_900GeVMC_sepXi3 900GeV 3 3 3 3 3 0.01
Qsub -l lnxfarm_quad4 -e -N doAn9_sepXi5 -o doAn_systVar_900GeVData_sepXi5.log doAn systVar_900GeVData_sepXi5 systVar_900GeVMC_sepXi5 900GeV 3 3 3 3 5 0.01

Qsub -l lnxfarm_quad4 -e -N doAn9_ipXi2 -o doAn_systVar_900GeVData_ipXi2.log doAn systVar_900GeVData_ipXi2 systVar_900GeVMC_ipXi2 900GeV 2 3 3 3 4 0.01
Qsub -l lnxfarm_quad4 -e -N doAn9_ipXi4 -o doAn_systVar_900GeVData_ipXi4.log doAn systVar_900GeVData_ipXi4 systVar_900GeVMC_ipXi4 900GeV 4 3 3 3 4 0.01


Qsub -l lnxfarm_quad4 -e -N doAn7_nom -o doAn_systVar_7TeVData_nominal.log doAn systVar_7TeVData_nominal systVar_7TeVMC_nominal 7TeV 3 3 3 3 4 0.01

Qsub -l lnxfarm_quad4 -e -N doAn7_ip2 -o doAn_systVar_7TeVData_ip2.log doAn systVar_7TeVData_ip2 systVar_7TeVMC_ip2 7TeV 3 2 2 2 4 0.01
Qsub -l lnxfarm_quad4 -e -N doAn7_ip4 -o doAn_systVar_7TeVData_ip4.log doAn systVar_7TeVData_ip4 systVar_7TeVMC_ip4 7TeV 3 4 4 4 4 0.01

Qsub -l lnxfarm_quad4 -e -N doAn7_CL001 -o doAn_systVar_7TeVData_CL001.log doAn systVar_7TeVData_CL001 systVar_7TeVMC_CL001 7TeV 3 3 3 3 4 0.001
Qsub -l lnxfarm_quad4 -e -N doAn7_CL02 -o doAn_systVar_7TeVData_CL02.log doAn systVar_7TeVData_CL02 systVar_7TeVMC_CL02 7TeV 3 3 3 3 4 0.02

Qsub -l lnxfarm_quad4 -e -N doAn7_sepXi3 -o doAn_systVar_7TeVData_sepXi3.log doAn systVar_7TeVData_sepXi3 systVar_7TeVMC_sepXi3 7TeV 3 3 3 3 3 0.01
Qsub -l lnxfarm_quad4 -e -N doAn7_sepXi5 -o doAn_systVar_7TeVData_sepXi5.log doAn systVar_7TeVData_sepXi5 systVar_7TeVMC_sepXi5 7TeV 3 3 3 3 5 0.01

Qsub -l lnxfarm_quad4 -e -N doAn7_ipXi2 -o doAn_systVar_7TeVData_ipXi2.log doAn systVar_7TeVData_ipXi2 systVar_7TeVMC_ipXi2 7TeV 2 3 3 3 4 0.01
Qsub -l lnxfarm_quad4 -e -N doAn7_ipXi4 -o doAn_systVar_7TeVData_ipXi4.log doAn systVar_7TeVData_ipXi4 systVar_7TeVMC_ipXi4 7TeV 4 3 3 3 4 0.01



*/

  // Style and other ROOT commands
  gROOT->Reset();
  setTDRStyle();
  gStyle->SetErrorX(0.5); 
  gStyle->SetTitleSize(0.05, "XYZ");
  gStyle->SetLabelSize(0.045, "XYZ");
  gStyle->SetTitleXOffset(1.1);
  gStyle->SetTitleYOffset(1.5);

  gStyle->SetPadLeftMargin(0.20);
  gStyle->SetTitleYOffset(1.5);
  TGaxis::SetMaxDigits(3);

// get 7TeV plots
  TFile* fin1 = TFile::Open("systVar_7TeVData_nominal/xiOut.root", "READ");
  TH1F* pt7 = numEffpTScaled; 
  TH1F* y7 = numEffYScaled;     
  
  TFile* f1 = TFile::Open("systVar_7TeVData_ip2/xiOut.root");
  TH1F* pt7_trackIP2 = numEffpTScaled;
  TH1F* y7_trackIP2 = numEffYScaled;
  TFile* f1 = TFile::Open("systVar_7TeVData_ip4/xiOut.root");
  TH1F* pt7_trackIP4 = numEffpTScaled;
  TH1F* y7_trackIP4 = numEffYScaled;

  TFile* f2 = TFile::Open("systVar_7TeVData_CL001/xiOut.root");
  TH1F* pt7_vtxCL001 = numEffpTScaled;
  TH1F* y7_vtxCL001 = numEffYScaled;
  TFile* f2 = TFile::Open("systVar_7TeVData_CL02/xiOut.root");
  TH1F* pt7_vtxCL02 = numEffpTScaled;
  TH1F* y7_vtxCL02 = numEffYScaled;  
  
  TFile* f2 = TFile::Open("systVar_7TeVData_ipXi2/xiOut.root");
  TH1F* pt7_ipXi2 = numEffpTScaled;
  TH1F* y7_ipXi2 = numEffYScaled;
  TFile* f2 = TFile::Open("systVar_7TeVData_ipXi4/xiOut.root");
  TH1F* pt7_ipXi4 = numEffpTScaled;
  TH1F* y7_ipXi4 = numEffYScaled;

  TFile* f2 = TFile::Open("systVar_7TeVData_sepXi3/xiOut.root");
  TH1F* pt7_sepXi3 = numEffpTScaled;
  TH1F* y7_sepXi3 = numEffYScaled;
  TFile* f2 = TFile::Open("systVar_7TeVData_sepXi5/xiOut.root");
  TH1F* pt7_sepXi5 = numEffpTScaled;
  TH1F* y7_sepXi5 = numEffYScaled;
  
// get  900 GeV plots
  TFile* fin2 = TFile::Open("systVar_900GeVData_nominal/xiOut.root", "READ");
  TH1F* pt900 = numEffpTScaled;
  TH1F* y900 = numEffYScaled;

  TFile* f1 = TFile::Open("systVar_900GeVData_ip2/xiOut.root");
  TH1F* pt900_trackIP2 = numEffpTScaled;
  TH1F* y900_trackIP2 = numEffYScaled;
  TFile* f1 = TFile::Open("systVar_900GeVData_ip4/xiOut.root");
  TH1F* pt900_trackIP4 = numEffpTScaled;
  TH1F* y900_trackIP4 = numEffYScaled;

  TFile* f2 = TFile::Open("systVar_900GeVData_CL001/xiOut.root");
  TH1F* pt900_vtxCL001 = numEffpTScaled;
  TH1F* y900_vtxCL001 = numEffYScaled;
  TFile* f2 = TFile::Open("systVar_900GeVData_CL02/xiOut.root");
  TH1F* pt900_vtxCL02 = numEffpTScaled;
  TH1F* y900_vtxCL02 = numEffYScaled;
  
  TFile* f2 = TFile::Open("systVar_900GeVData_ipXi2/xiOut.root");
  TH1F* pt900_ipXi2 = numEffpTScaled;
  TH1F* y900_ipXi2 = numEffYScaled;
  TFile* f2 = TFile::Open("systVar_900GeVData_ipXi4/xiOut.root");
  TH1F* pt900_ipXi4 = numEffpTScaled;
  TH1F* y900_ipXi4 = numEffYScaled;

  TFile* f2 = TFile::Open("systVar_900GeVData_sepXi3/xiOut.root");
  TH1F* pt900_sepXi3 = numEffpTScaled;
  TH1F* y900_sepXi3 = numEffYScaled;
  TFile* f2 = TFile::Open("systVar_900GeVData_sepXi5/xiOut.root");
  TH1F* pt900_sepXi5 = numEffpTScaled;
  TH1F* y900_sepXi5 = numEffYScaled;
  
  TLatex* cmsPre = new TLatex();
  cmsPre->SetNDC();
  cmsPre->SetTextSize(0.045);
  cmsPre->SetTextAlign(22);
  cmsPre->SetText(0.78,0.90,"CMS preliminary");

  TLatex* cmsPreUL = new TLatex();
  cmsPreUL->SetNDC();
  cmsPreUL->SetTextSize(0.045);
  cmsPreUL->SetTextAlign(22);
  cmsPreUL->SetText(0.45,0.90,"CMS preliminary");



  TCanvas* c1 = new TCanvas("c1", "Xi distributions", 900, 1200); 
  c1->Divide(2,4);
  c1->cd(1);
  
  TLegend* leg = new TLegend(0.65, 0.7, 0.9, 0.8);
  leg->SetTextSize(0.035);
  leg->SetFillColor(0);
  pt900_trackIP2->Divide(pt900_trackIP2,pt900);
  pt900_trackIP2->GetYaxis()->SetTitle("Ratio track IP > 2#sigma/nominal");
  pt900_trackIP2->GetXaxis()->SetTitle("#Xi p_{T} (GeV)");
  pt900_trackIP2->SetLineColor(2);
  pt900_trackIP2->SetMarkerColor(2);
  pt900_trackIP2->SetMarkerStyle(21);
  pt900_trackIP2->Draw();  
  pt7_trackIP2->Divide(pt7_trackIP2,pt7);
  pt7_trackIP2->SetLineColor(4);
  pt7_trackIP2->SetMarkerColor(4);
  pt7_trackIP2->SetMarkerStyle(20);
  pt7_trackIP2->Draw("same");   
  leg->AddEntry(pt7_trackIP2, "#sqrt{s} = 7 TeV", "pl");
  leg->AddEntry(pt900_trackIP2, "#sqrt{s} = 900 GeV", "pl");
  leg->Draw();

  c1->cd(2);
  pt900_trackIP4->Divide(pt900_trackIP4,pt900);
  pt900_trackIP4->GetYaxis()->SetTitle("Ratio track IP > 4#sigma/nominal");
  pt900_trackIP4->GetXaxis()->SetTitle("#Xi p_{T} (GeV)");
  pt900_trackIP4->SetLineColor(2);
  pt900_trackIP4->SetMarkerColor(2);
  pt900_trackIP4->SetMarkerStyle(21);
  pt900_trackIP4->Draw();  
  
  pt7_trackIP4->Divide(pt7_trackIP4,pt7);
  pt7_trackIP4->SetLineColor(4);
  pt7_trackIP4->SetMarkerColor(4);
  pt7_trackIP4->SetMarkerStyle(20);
  pt7_trackIP4->Draw("same");   

  c1->cd(3);  
  pt900_vtxCL001->Divide(pt900_vtxCL001,pt900);
  pt900_vtxCL001->GetYaxis()->SetTitle("Ratio #Xi^{-} vtxCL>0.001/nominal");
  pt900_vtxCL001->GetXaxis()->SetTitle("#Xi p_{T} (GeV)");
  pt900_vtxCL001->SetLineColor(2);
  pt900_vtxCL001->SetMarkerColor(2);
  pt900_vtxCL001->SetMarkerStyle(21);
  pt900_vtxCL001->Draw();   

  pt7_vtxCL001->Divide(pt7_vtxCL001,pt7);
  pt7_vtxCL001->SetLineColor(4);
  pt7_vtxCL001->SetMarkerColor(4);
  pt7_vtxCL001->SetMarkerStyle(20);
  pt7_vtxCL001->Draw("same");  
  
  c1->cd(4);  
  pt900_vtxCL02->Divide(pt900_vtxCL02,pt900);
  pt900_vtxCL02->GetYaxis()->SetTitle("Ratio #Xi^{-} vtxCL>0.02/nominal");
  pt900_vtxCL02->GetXaxis()->SetTitle("#Xi p_{T} (GeV)");
  pt900_vtxCL02->SetLineColor(2);
  pt900_vtxCL02->SetMarkerColor(2);
  pt900_vtxCL02->SetMarkerStyle(21);
  pt900_vtxCL02->Draw();   

  pt7_vtxCL02->Divide(pt7_vtxCL02,pt7);
  pt7_vtxCL02->SetLineColor(4);
  pt7_vtxCL02->SetMarkerColor(4);
  pt7_vtxCL02->SetMarkerStyle(20);
  pt7_vtxCL02->Draw("same");  
  
  c1->cd(5);  
  pt900_ipXi2->Divide(pt900_ipXi2,pt900);
  pt900_ipXi2->GetYaxis()->SetTitle("Ratio #Xi^{-} IP > 2#sigma/nominal");
  pt900_ipXi2->GetXaxis()->SetTitle("#Xi p_{T} (GeV)");
  pt900_ipXi2->SetLineColor(2);
  pt900_ipXi2->SetMarkerColor(2);
  pt900_ipXi2->SetMarkerStyle(21);
  pt900_ipXi2->Draw();   

  pt7_ipXi2->Divide(pt7_ipXi2,pt7);
  pt7_ipXi2->SetLineColor(4);
  pt7_ipXi2->SetMarkerColor(4);
  pt7_ipXi2->SetMarkerStyle(20);
  pt7_ipXi2->Draw("same");  
  
  c1->cd(6);  
  pt900_ipXi4->Divide(pt900_ipXi4,pt900);
  pt900_ipXi4->GetYaxis()->SetTitle("Ratio #Xi^{-} IP > 4#sigma/nominal");
  pt900_ipXi4->GetXaxis()->SetTitle("#Xi p_{T} (GeV)");
  pt900_ipXi4->SetLineColor(2);
  pt900_ipXi4->SetMarkerColor(2);
  pt900_ipXi4->SetMarkerStyle(21);
  pt900_ipXi4->Draw();   

  pt7_ipXi4->Divide(pt7_ipXi4,pt7);
  pt7_ipXi4->SetLineColor(4);
  pt7_ipXi4->SetMarkerColor(4);
  pt7_ipXi4->SetMarkerStyle(20);
  pt7_ipXi4->Draw("same");  

  c1->cd(7);  
  pt900_sepXi3->Divide(pt900_sepXi3,pt900);
  pt900_sepXi3->GetYaxis()->SetTitle("Ratio #Xi^{-} vtx sep > 3#sigma/nominal");
  pt900_sepXi3->GetXaxis()->SetTitle("#Xi p_{T} (GeV)");
  pt900_sepXi3->SetLineColor(2);
  pt900_sepXi3->SetMarkerColor(2);
  pt900_sepXi3->SetMarkerStyle(21);
  pt900_sepXi3->Draw();   

  pt7_sepXi3->Divide(pt7_sepXi3,pt7);
  pt7_sepXi3->SetLineColor(4);
  pt7_sepXi3->SetMarkerColor(4);
  pt7_sepXi3->SetMarkerStyle(20);
  pt7_sepXi3->Draw("same");  
  
  c1->cd(8);  
  pt900_sepXi5->Divide(pt900_sepXi5,pt900);
  pt900_sepXi5->GetYaxis()->SetTitle("Ratio #Xi^{-} vtx sep > 5#sigma/nominal");
  pt900_sepXi5->GetXaxis()->SetTitle("#Xi p_{T} (GeV)");
  pt900_sepXi5->SetLineColor(2);
  pt900_sepXi5->SetMarkerColor(2);
  pt900_sepXi5->SetMarkerStyle(21);
  pt900_sepXi5->Draw();   

  pt7_sepXi5->Divide(pt7_sepXi5,pt7);
  pt7_sepXi5->SetLineColor(4);
  pt7_sepXi5->SetMarkerColor(4);
  pt7_sepXi5->SetMarkerStyle(20);
  pt7_sepXi5->Draw("same");  
  
  c1->SaveAs("XiCuts_pt.eps");

  TCanvas* c2 = new TCanvas("c2", "Xi distributions", 900, 1200); 
  c2->Divide(2,4);
  c2->cd(1);
  
  TLegend* leg = new TLegend(0.65, 0.7, 0.9, 0.8);
  leg->SetTextSize(0.035);
  leg->SetFillColor(0);
  
  y900_trackIP2->Divide(y900_trackIP2,y900);
  y900_trackIP2->GetYaxis()->SetTitle("Ratio track IP > 2#sigma/nominal");
  y900_trackIP2->GetXaxis()->SetTitle("#Xi rapidity");
  y900_trackIP2->SetLineColor(2);
  y900_trackIP2->SetMarkerColor(2);
  y900_trackIP2->SetMarkerStyle(21);
  y900_trackIP2->Draw();  
  
  y7_trackIP2->Divide(y7_trackIP2,y7);
  y7_trackIP2->SetLineColor(4);
  y7_trackIP2->SetMarkerColor(4);
  y7_trackIP2->SetMarkerStyle(20);
  y7_trackIP2->Draw("same");   
  
  leg->AddEntry(y7_trackIP2, "#sqrt{s} = 7 TeV", "pl");
  leg->AddEntry(y900_trackIP2, "#sqrt{s} = 900 GeV", "pl");
  leg->Draw();

  c2->cd(2);
  y900_trackIP4->Divide(y900_trackIP4,y900);
  y900_trackIP4->GetYaxis()->SetTitle("Ratio track IP > 4#sigma/nominal");
  y900_trackIP4->GetXaxis()->SetTitle("#Xi rapidity");
  y900_trackIP4->SetLineColor(2);
  y900_trackIP4->SetMarkerColor(2);
  y900_trackIP4->SetMarkerStyle(21);
  y900_trackIP4->Draw();  
  
  y7_trackIP4->Divide(y7_trackIP4,y7);
  y7_trackIP4->SetLineColor(4);
  y7_trackIP4->SetMarkerColor(4);
  y7_trackIP4->SetMarkerStyle(20);
  y7_trackIP4->Draw("same");   

  c2->cd(3);  
  y900_vtxCL001->Divide(y900_vtxCL001,y900);
  y900_vtxCL001->GetYaxis()->SetTitle("Ratio #Xi^{-} vtxCL>0.001/nominal");
  y900_vtxCL001->GetXaxis()->SetTitle("#Xi rapidity");
  y900_vtxCL001->SetLineColor(2);
  y900_vtxCL001->SetMarkerColor(2);
  y900_vtxCL001->SetMarkerStyle(21);
  y900_vtxCL001->Draw();   

  y7_vtxCL001->Divide(y7_vtxCL001,y7);
  y7_vtxCL001->SetLineColor(4);
  y7_vtxCL001->SetMarkerColor(4);
  y7_vtxCL001->SetMarkerStyle(20);
  y7_vtxCL001->Draw("same");  
  
  c2->cd(4);  
  y900_vtxCL02->Divide(y900_vtxCL02,y900);
  y900_vtxCL02->GetYaxis()->SetTitle("Ratio #Xi^{-} vtxCL>0.02/nominal");
  y900_vtxCL02->GetXaxis()->SetTitle("#Xi rapidity");
  y900_vtxCL02->SetLineColor(2);
  y900_vtxCL02->SetMarkerColor(2);
  y900_vtxCL02->SetMarkerStyle(21);
  y900_vtxCL02->Draw();   

  y7_vtxCL02->Divide(y7_vtxCL02,y7);
  y7_vtxCL02->SetLineColor(4);
  y7_vtxCL02->SetMarkerColor(4);
  y7_vtxCL02->SetMarkerStyle(20);
  y7_vtxCL02->Draw("same");  
  
  c2->cd(5);  
  y900_ipXi2->Divide(y900_ipXi2,y900);
  y900_ipXi2->GetYaxis()->SetTitle("Ratio #Xi^{-} IP > 2#sigma/nominal");
  y900_ipXi2->GetXaxis()->SetTitle("#Xi rapidity");
  y900_ipXi2->SetLineColor(2);
  y900_ipXi2->SetMarkerColor(2);
  y900_ipXi2->SetMarkerStyle(21);
  y900_ipXi2->Draw();   

  y7_ipXi2->Divide(y7_ipXi2,y7);
  y7_ipXi2->SetLineColor(4);
  y7_ipXi2->SetMarkerColor(4);
  y7_ipXi2->SetMarkerStyle(20);
  y7_ipXi2->Draw("same");  
  
  c2->cd(6);  
  y900_ipXi4->Divide(y900_ipXi4,y900);
  y900_ipXi4->GetYaxis()->SetTitle("Ratio #Xi^{-} IP > 4#sigma/nominal");
  y900_ipXi4->GetXaxis()->SetTitle("#Xi rapidity");
  y900_ipXi4->SetLineColor(2);
  y900_ipXi4->SetMarkerColor(2);
  y900_ipXi4->SetMarkerStyle(21);
  y900_ipXi4->Draw();   

  y7_ipXi4->Divide(y7_ipXi4,y7);
  y7_ipXi4->SetLineColor(4);
  y7_ipXi4->SetMarkerColor(4);
  y7_ipXi4->SetMarkerStyle(20);
  y7_ipXi4->Draw("same");  

  c2->cd(7);  
  y900_sepXi3->Divide(y900_sepXi3,y900);
  y900_sepXi3->GetYaxis()->SetTitle("Ratio #Xi^{-} vtx sep > 3#sigma/nominal");
  y900_sepXi3->GetXaxis()->SetTitle("#Xi rapidity");
  y900_sepXi3->SetLineColor(2);
  y900_sepXi3->SetMarkerColor(2);
  y900_sepXi3->SetMarkerStyle(21);
  y900_sepXi3->Draw();   

  y7_sepXi3->Divide(y7_sepXi3,y7);
  y7_sepXi3->SetLineColor(4);
  y7_sepXi3->SetMarkerColor(4);
  y7_sepXi3->SetMarkerStyle(20);
  y7_sepXi3->Draw("same");  
  
  c2->cd(8);  
  y900_sepXi5->Divide(y900_sepXi5,y900);
  y900_sepXi5->GetYaxis()->SetTitle("Ratio #Xi^{-} vtx sep > 5#sigma/nominal");
  y900_sepXi5->GetXaxis()->SetTitle("#Xi rapidity");
  y900_sepXi5->SetLineColor(2);
  y900_sepXi5->SetMarkerColor(2);
  y900_sepXi5->SetMarkerStyle(21);
  y900_sepXi5->Draw();   

  y7_sepXi5->Divide(y7_sepXi5,y7);
  y7_sepXi5->SetLineColor(4);
  y7_sepXi5->SetMarkerColor(4);
  y7_sepXi5->SetMarkerStyle(20);
  y7_sepXi5->Draw("same");  
  
  c2->SaveAs("XiCuts_y.eps");


}
