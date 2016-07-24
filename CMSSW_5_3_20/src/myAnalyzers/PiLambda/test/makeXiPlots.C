#include "tsallis.C"
//#include </afs/cern.ch/user/k/kaulmer/tdrstyle.C>
#include </usr/users/ulmerk/stuff/tdrstyle_keith.C>
void makeXiPlots()
{

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
//  TGaxis::SetMaxDigits(3);

  ofstream outTextFile;
  outTextFile.open ("fitinfo.txt");

  TString outputDir = TString("fitResults_900GeVMC_nom_p01p02_singleCand/");
  gSystem->mkdir(outputDir);
  
  float c = 0.0299793; // cm/ps
  TF1 *exp1 = new TF1("exp1", "[0]*exp(x*[1])", 0.0, 10.0);
  exp1->SetParameters(700000,-0.16);

  TFile* fin1 = TFile::Open("7TeVData_nom_delRp06_test/xiOut.root", "READ");
  TH1F* pT7TeV = numEffpTScaled;
  TH1F* pT7TeVFit = numFitPt;
  TH1F* y7TeV = numEffYScaled;     
  TH1F* y7TeVFit = numFitY;
  TH1F* ct7TeV = numEffCtauScaled;
  TH1F* ct7TeVFit = numFitCTau;
  TH1F* thePosYieldsPt7TeVHisto = posNumFitPt;
  thePosYieldsPt7TeVHisto->Rebin(2);
  TH1F* theNegYieldsPt7TeVHisto = negNumFitPt;
  theNegYieldsPt7TeVHisto->Rebin(2);
  TH1F* thePosCorrYieldsPt7TeVHisto = numEffPosXipTScaled;
  thePosCorrYieldsPt7TeVHisto->Rebin(2);
  TH1F* theNegCorrYieldsPt7TeVHisto = numEffNegXipTScaled;
  theNegCorrYieldsPt7TeVHisto->Rebin(2);
  TH1F* thePosYieldsY7TeVHisto = posNumFitY;
  thePosYieldsY7TeVHisto->Rebin(2);
  TH1F* theNegYieldsY7TeVHisto = negNumFitY;
  theNegYieldsY7TeVHisto->Rebin(2);
  TH1F* thePosCorrYieldsY7TeVHisto = numEffPosXiYScaled;
  thePosCorrYieldsY7TeVHisto->Rebin(2);
  TH1F* theNegCorrYieldsY7TeVHisto = numEffNegXiYScaled;
  theNegCorrYieldsY7TeVHisto->Rebin(2);
  TH1F* data7TeV = data;

  TFile* fin2 = TFile::Open("900GeVData_nom_p01p02_singleCand/xiOut.root", "READ");
  TH1F* pT900 = numEffpTScaled;
  TH1F* pT900Fit = numFitPt;
  TH1F* y900 = numEffYScaled;
  TH1F* y900Fit = numFitY;
  TH1F* ct900 = numEffCtauScaled;
  TH1F* ct900Fit = numFitCTau;
  TH1F* thePosYieldsPt900GeVHisto = posNumFitPt;
  thePosYieldsPt900GeVHisto->Rebin(2);
  TH1F* theNegYieldsPt900GeVHisto = negNumFitPt;
  theNegYieldsPt900GeVHisto->Rebin(2);
  TH1F* thePosCorrYieldsPt900GeVHisto = numEffPosXipTScaled;
  thePosCorrYieldsPt900GeVHisto->Rebin(2);
  TH1F* theNegCorrYieldsPt900GeVHisto = numEffNegXipTScaled;
  theNegCorrYieldsPt900GeVHisto->Rebin(2);
  TH1F* thePosYieldsY900GeVHisto = posNumFitY;
  thePosYieldsY900GeVHisto->Rebin(2);
  TH1F* theNegYieldsY900GeVHisto = negNumFitY;
  theNegYieldsY900GeVHisto->Rebin(2);
  TH1F* thePosCorrYieldsY900GeVHisto = numEffPosXiYScaled;
  thePosCorrYieldsY900GeVHisto->Rebin(2);
  TH1F* theNegCorrYieldsY900GeVHisto = numEffNegXiYScaled;
  theNegCorrYieldsY900GeVHisto->Rebin(2);
  TH1F* data900 = data;

  TFile* fin3 = TFile::Open("7TeVMC_nom_p01p02_singleCand/efficienciesXi_7TeV.root", "READ");
  TH1F* pT7TeVMCGen = allGenpT;
  TH1F* pT7TeVMC = numEffpTScaled; 
  TH1F* pt7TeVMCforRatio = denEffpTUncor;  
  TH1F* pt7TeVMCforFitVal = denEffpTUncor->Clone("pt7TeVMCforFitVal");
  TH1F* effpT7TeV = effpT;
  TH1F* y7TeVMCGen = allGenY;
  TH1F* y7TeVMC = numEffYScaled;
  TH1F* y7TeVMCforRatio = denEffYUncor;  
  TH1F* y7TeVMCforFitVal = denEffYUncor->Clone("y7TeVMCforFitVal");
  TH1F* effY7TeV = effY;
  TH1F* ct7TeVMC = numEffCtauScaled;
  TH1F* ct7TeVMCforFitVal = denEffCtauUncor->Clone("ct7TeVMCforFitVal");
  TH1F* effCtau7TeV = effCtau;

  TFile* fin4 = TFile::Open("900GeVMC_nom_p01p02_singleCand/efficienciesXi_900GeV.root", "READ");
  TH1F* pT900MCGen = allGenpT;
  TH1F* pT900MC = numEffpTScaled;
  TH1F* pt900MCforRatio = denEffpTUncor;  
  TH1F* pt900MCforFitVal = denEffpTUncor->Clone("pt900MCforFitVal");  
  TH1F* effpT900 = effpT;
  TH1F* y900MCGen = allGenY;
  TH1F* y900MC = numEffYScaled;
  TH1F* y900MCforRatio = denEffYUncor;  
  TH1F* y900MCforFitVal = denEffYUncor->Clone("y900MCforFitVal");
  TH1F* effY900 = effY;
  TH1F* ct900MC = numEffCtauScaled;
  TH1F* ct900MCforFitVal = denEffCtauUncor->Clone("ct900forFitVal");
  TH1F* effCtau900 = effCtau;

  TFile* fin5 = TFile::Open("900GeVMC_nom_delRp06_noEnhancedMC_test/efficienciesXi_900GeV.root", "READ");
  TH1F* pT900MCGenNoEn = allGenpT;
  TH1F* y900MCGenNoEn = allGenY;

  //TFile* fin6 = TFile::Open("systVar_900GeVMC_looseDelR_D6T/efficienciesXi_900GeV.root", "READ");
  TH1F* pT900MCGenD6T = allGenpT;
  TH1F* y900MCGenD6T = allGenY;

  //TFile* fin7 = TFile::Open("systVar_900GeVMC_looseDelR_P0/efficienciesXi_900GeV.root", "READ");
  TH1F* pT900MCGenp0 = allGenpT;
  TH1F* y900MCGenp0 = allGenY;
  
  //TFile* fin8 = TFile::Open("systVar_900GeVMC_looseDelR_pythia8/efficienciesXi_900GeV.root", "READ");
  TH1F* pT900MCGenp8 = allGenpT;
  TH1F* y900MCGenp8 = allGenY;

  //TFile* fin9 = TFile::Open("systVar_7TeVMC_looseDelR_P0/efficienciesXi_7TeV.root", "READ");
  TH1F* pT7TeVMCGenp0 = allGenpT;
  TH1F* y7TeVMCGenp0 = allGenY;
  
  //TFile* fin10 = TFile::Open("systVar_7TeVMC_looseDelR_pythia8/efficienciesXi_7TeV.root", "READ");
  TH1F* pT7TeVMCGenp8 = allGenpT;
  TH1F* y7TeVMCGenp8 = allGenY;

  Double_t mXi = 1.3213;

  gStyle->SetOptStat(0);

  float nEvt7 = 37100402;
  float nEvt900 = 9953161;
  float lumi7 = 460895;  //milibarns  //still old numbers
  float lumi900 = 264427;  //milibarns  //still old numbers
  

  float systUncory = 0.095;
  float systCory7 = 0.105;
  float systCory9 = 0.177;
  float systRaty7 = 0.121;  //corr+uncorr
  float systRaty9 = 0.188;  //corr+uncorr
  float systUncorpt = 0.095;
  float systCorpt = 0.076;
  float systRatpt = 0.098; //corr+uncorr
   
  float nEvt7MC = 46788056;
  float nEvt900MC = 110266412;
  float nEvt900MCnoEn = 8365068;
  float lumi7MC = 0;
  float lumi900MC = 0;

  float nEvt900MCD6T = 118631480;
  float nEvt7MCp0 = 8819948;
  float nEvt900MCp0 = 7812712;
  float nEvt7MCp8 = 8090336;
  float nEvt900MCp8 = 7057211;

  TLatex* cmsPre = new TLatex();
  cmsPre->SetNDC();
  cmsPre->SetTextSize(0.045);
  cmsPre->SetTextAlign(22);
  cmsPre->SetText(0.78,0.90,"CMS Preliminary");

  TLatex* cmsPreUL = new TLatex();
  cmsPreUL->SetNDC();
  cmsPreUL->SetTextSize(0.045);
  cmsPreUL->SetTextAlign(22);
  cmsPreUL->SetText(0.45,0.90,"CMS Preliminary");

  TLatex* cmsPreUL2 = new TLatex();
  cmsPreUL2->SetNDC();
  cmsPreUL2->SetTextSize(0.045);
  cmsPreUL2->SetTextAlign(22);
  cmsPreUL2->SetText(0.385,0.90,"CMS Preliminary");



// first make the Cascade full mass plots
  TCanvas* c0 = new TCanvas("c0", "Xi mass", 600, 600); 
  c0->SetLeftMargin(0.17);
  c0->SetRightMargin(0.03);

  data7TeV->GetXaxis()->SetNdivisions(505);
  data7TeV->GetYaxis()->SetTitleOffset(1.65);
  data7TeV->Draw();  
  int bins = 100;
  int x_min = 1255;
  int x_max = 1455;
  int x_range = x_max-x_min;
  int binwidth = x_range/bins;
 TString title = TString("#Lambda^{0} #pi^{-} invariant mass [MeV/c^{2}]");
   // build ARGUS + Gaussian
  TF1 *f1 = new TF1("f1", Form("%s%i%s%i%s%i%s%i%s","[0]*((x-", x_min, ")**(0.5)) +[1]*((x-", x_min, ")**(3/2)) + ((", x_range, ")/", bins, ")*gausn(2)"), x_min, x_max);

//  f1->SetParameters(0.55, 0., 120, 1320, 60);
  f1->SetParName(0,"bkg param 1");
  f1->SetParName(1,"bkg param 2");
  f1->SetParName(2,"yield");
  f1->SetParLimits(2,0,100000);
  f1->SetParName(3,"mean");
  f1->SetParName(4,"sigma");
  f1->SetParLimits(4,0,1000);
  

  
  data7TeV->GetXaxis()->SetTitle(title);
  data7TeV->GetYaxis()->SetTitle(Form("%s%i%s","Candidates / ", binwidth, " MeV/c^{2}"));
  data7TeV->SetStats(0);



//  f1->SetParameters(0.6, 0., 18062.8, 1322, 3.25);
  f1->SetParameters(68.2, -0.01, 18062, 1322, 3.25);

  data7TeV->Fit("f1","LEM0");
  f1->GetXaxis()->SetTitle(title);
  f1->GetYaxis()->SetTitle(Form("%s%i%s","Candidates / ", binwidth, " MeV/c^{2}"));
  f1->SetNpx(5000);
  data7TeV->Draw();
  f1->Draw("same");

  float nSig = data7TeV->GetFunction("f1")->GetParameter(2);
  float mean = data7TeV->GetFunction("f1")->GetParameter(3);
  float sigma = data7TeV->GetFunction("f1")->GetParameter(4);
  float nSigE = data7TeV->GetFunction("f1")->GetParError(2);
  float meanE = data7TeV->GetFunction("f1")->GetParError(3);
  float sigmaE = data7TeV->GetFunction("f1")->GetParError(4);

  TPaveText *stats = new TPaveText(0.47,0.60,0.89,0.80,"blNDC");
  TText *text1 = stats->AddText(Form("%s%.1f%s%.1f", "Yield: ", nSig, " #pm ",nSigE));
  TText *text2 = stats->AddText(Form("%s%.2f%s%.2f%s", "Mean: ", mean, " #pm ", meanE, " MeV/c^{2}"));
  TText *text3 = stats->AddText(Form("%s%.2f%s%.2f%s", "Sigma: ", sigma, " #pm ", sigmaE, " MeV/c^{2}"));
  TText *text4 = stats->AddText("Statistical uncertainties only");
  stats->SetFillColor(0);
  stats->SetTextSize(0.038);
  stats->SetTextAlign(11);
  stats->Draw("same");



  //get background function and draw it
  float bkgA = data7TeV->GetFunction("f1")->GetParameter(0);
  float bkgB = data7TeV->GetFunction("f1")->GetParameter(1);
  TF1 *f30 = new TF1("f3", Form("%s%i%s%i%s","[0]*((x-", x_min, ")**(0.5)) +[1]*((x-", x_min, ")**(3/2))"), x_min, x_max);

  f30->SetParameter(0, bkgA);
  f30->SetParameter(1, bkgB);
  f30->SetLineStyle(2);
  f30->Draw("same");

  TLatex *cms = new TLatex();
  cms->SetNDC();
  cms->SetTextSize(0.040);
  cms->DrawLatex(0.60,0.89,"CMS Preliminary");
  TLatex *cms2 = new TLatex();
  cms2->SetNDC();
  cms2->SetTextSize(0.040);
  cms2->DrawLatex(0.653,0.85,"#sqrt{s} = 7 TeV");

  c0->SaveAs(outputDir+"XiMass_selected_data_7.eps");
  
  
//  f1->SetParameters(0.55, 0., 120, 1320, 60);
  data900->GetXaxis()->SetNdivisions(505);
  data900->GetYaxis()->SetTitleOffset(1.65);
  data900->Draw();  
  
  data900->GetXaxis()->SetTitle(title);
  data900->GetYaxis()->SetTitle(Form("%s%i%s","Candidates / ", binwidth, " MeV/c^{2}"));
  data900->SetStats(0);

//  f1->SetParameters(0.6, 0., 50, 1320, 40);
  f1->SetParameters(8.8, 0.0015, 26943, 1322, 3.36);

  data900->Fit("f1","LEM0");
  f1->GetXaxis()->SetTitle(title);
  f1->GetYaxis()->SetTitle(Form("%s%i%s","Candidates / ", binwidth, " MeV/c^{2}"));
  f1->SetNpx(5000);
  data900->Draw();
  f1->Draw("same");

  float nSig = data900->GetFunction("f1")->GetParameter(2);
  float mean = data900->GetFunction("f1")->GetParameter(3);
  float sigma = data900->GetFunction("f1")->GetParameter(4);
  float nSigE = data900->GetFunction("f1")->GetParError(2);
  float meanE = data900->GetFunction("f1")->GetParError(3);
  float sigmaE = data900->GetFunction("f1")->GetParError(4);

  TPaveText *stats = new TPaveText(0.47,0.60,0.89,0.80,"blNDC");
  TText *text1 = stats->AddText(Form("%s%.1f%s%.1f", "Yield: ", nSig, " #pm ",nSigE));
  TText *text2 = stats->AddText(Form("%s%.2f%s%.2f%s", "Mean: ", mean, " #pm ", meanE, " MeV/c^{2}"));
  TText *text3 = stats->AddText(Form("%s%.2f%s%.2f%s", "Sigma: ", sigma, " #pm ", sigmaE, " MeV/c^{2}"));
  TText *text4 = stats->AddText("Statistical uncertainties only");
  stats->SetFillColor(0);
  stats->SetTextSize(0.038);
  stats->SetTextAlign(11);
  stats->Draw("same");



  //get background function and draw it
  float bkgA = data900->GetFunction("f1")->GetParameter(0);
  float bkgB = data900->GetFunction("f1")->GetParameter(1);
  TF1 *f30 = new TF1("f3", Form("%s%i%s%i%s","[0]*((x-", x_min, ")**(0.5)) +[1]*((x-", x_min, ")**(3/2))"), x_min, x_max);

  f30->SetParameter(0, bkgA);
  f30->SetParameter(1, bkgB);
  f30->SetLineStyle(2);
  f30->Draw("same");

  TLatex *cms = new TLatex();
  cms->SetNDC();
  cms->SetTextSize(0.040);
  cms->DrawLatex(0.60,0.89,"CMS Preliminary");
  TLatex *cms2 = new TLatex();
  cms2->SetNDC();
  cms2->SetTextSize(0.040);
  cms2->DrawLatex(0.631,0.85,"#sqrt{s} = 900 GeV");

  c0->SaveAs(outputDir+"XiMass_selected_data_900.eps");



// second make the uncorrected, efficiency and corrected yield distribution plots

  TCanvas* c2 = new TCanvas("c2", "Xi distributions", 900, 300); 
  c2->Divide(3,1);
  c2->cd(1);
  y7TeVFit->GetYaxis()->SetTitle("Candidates");
  y7TeVFit->GetXaxis()->SetTitle("#Xi^{-} rapidity");
  y7TeVFit->SetMinimum(0);
  y7TeVFit->Draw();
  cmsPre->Draw("same");
  c2->cd(2);
  effY7TeV->GetYaxis()->SetTitle("Efficiency");
  effY7TeV->GetXaxis()->SetTitle("#Xi^{-} rapidity");  
  effY7TeV->Draw();
  c2->Update();
  c2->cd(3);
  y7TeV->GetYaxis()->SetTitle("Candidates/unit Y");
  y7TeV->GetXaxis()->SetTitle("#Xi^{-} rapidity");  
  y7TeV->SetMinimum(0);
  y7TeV->Draw();
  c2->SaveAs(outputDir+"XiEffY_data_7.eps");

  c2->cd(1);
  y900Fit->GetYaxis()->SetTitle("Candidates");
  y900Fit->GetXaxis()->SetTitle("#Xi^{-} rapidity");
  y900Fit->SetMinimum(0);
  y900Fit->Draw();
  cmsPre->Draw("same");
  c2->cd(2);
  effY900->GetYaxis()->SetTitle("Efficiency");
  effY900->GetXaxis()->SetTitle("#Xi^{-} rapidity");
  effY900->Draw();
  c2->Update();
  c2->cd(3);
  y900->GetYaxis()->SetTitle("Candidates/unit Y");
  y900->GetXaxis()->SetTitle("#Xi^{-} rapidity");  
  y900->SetMinimum(0);
  y900->Draw();
  c2->SaveAs(outputDir+"XiEffY_data_900.eps");

  c2->cd(1);
  pT7TeVFit->GetYaxis()->SetTitle("Candidates");
  pT7TeVFit->GetXaxis()->SetTitle("#Xi^{-} p_{T} (GeV)");
  pT7TeVFit->SetMinimum(0);
  pT7TeVFit->Draw();
  cmsPre->Draw("same");
  c2->cd(2);
  effpT7TeV->GetYaxis()->SetTitle("Efficiency");
  effpT7TeV->GetXaxis()->SetTitle("#Xi^{-} p_{T} (GeV)");  
  effpT7TeV->Draw();
  c2->Update();
  c2->cd(3);
  pT7TeV->GetYaxis()->SetTitle("Candidates/GeV");
  pT7TeV->GetXaxis()->SetTitle("#Xi^{-} p_{T} (GeV)");  
  pT7TeV->SetMinimum(0);
  pT7TeV->Draw();
  c2->SaveAs(outputDir+"XiEffpT_data_7.eps");

  c2->cd(1);
  pT900Fit->GetYaxis()->SetTitle("Candidates");
  pT900Fit->GetXaxis()->SetTitle("#Xi^{-} p_{T} (GeV)");
  pT900Fit->SetMinimum(0);
  pT900Fit->Draw();
  cmsPre->Draw("same");
  c2->cd(2);
  effpT900->GetYaxis()->SetTitle("Efficiency");
  effpT900->GetXaxis()->SetTitle("#Xi^{-} p_{T} (GeV)");
  effpT900->Draw();
  c2->Update();
  c2->cd(3);
  pT900->GetYaxis()->SetTitle("Candidates/GeV");
  pT900->GetXaxis()->SetTitle("#Xi^{-} p_{T} (GeV)");  
  pT900->SetMinimum(0);
  pT900->Draw();
  c2->SaveAs(outputDir+"XiEffpT_data_900.eps");
  
  c2->cd(1);
  ct7TeVFit->GetYaxis()->SetTitle("Candidates");
  ct7TeVFit->GetXaxis()->SetTitle("#Xi^{-} ct (cm)");
  ct7TeVFit->SetMinimum(0);
  ct7TeVFit->Draw();
  cmsPre->Draw("same");
  c2->cd(2);
  effCtau7TeV->GetYaxis()->SetTitle("Efficiency");
  effCtau7TeV->GetXaxis()->SetTitle("#Xi^{-} ct (cm)");  
  effCtau7TeV->SetMaximum( effCtau7TeV->GetMaximum()*1.2 );
  effCtau7TeV->Draw();
  c2->cd(3);
  ct7TeV->GetYaxis()->SetTitle("Candidates/cm");
  ct7TeV->GetXaxis()->SetTitle("#Xi^{-} ct (cm)");  
  ct7TeV->SetMinimum(0);
  ct7TeV->Draw();
  c2->Update();
  TPaveStats *st1 = (TPaveStats*)ct7TeV->FindObject("stats");
  st1->SetX1NDC(0.5);
  st1->SetTextSize(0.035);
  st1->Draw();
  c2->SaveAs(outputDir+"XiEffCtau_data_7.eps");

  c2->cd(1);
  ct900Fit->GetYaxis()->SetTitle("Candidates");
  ct900Fit->GetXaxis()->SetTitle("#Xi^{-} ct (cm)");
  ct900Fit->SetMinimum(0);
  ct900Fit->Draw();
  cmsPre->Draw("same");
  c2->cd(2);
  effCtau900->GetYaxis()->SetTitle("Efficiency");
  effCtau900->GetXaxis()->SetTitle("#Xi^{-} ct (cm)");
  effCtau900->SetMaximum( effCtau900->GetMaximum()*1.2 );
  effCtau900->Draw();
  c2->cd(3);
  ct900->GetYaxis()->SetTitle("Candidates/cm");
  ct900->GetXaxis()->SetTitle("#Xi^{-} ct (cm)");  
  ct900->SetMinimum(0);
  ct900->Draw();
  c2->Update();
  TPaveStats *st1 = (TPaveStats*)ct900->FindObject("stats");
  st1->SetX1NDC(0.5);
  st1->SetTextSize(0.035);
  st1->Draw();
  c2->SaveAs(outputDir+"XiEffCtau_data_900.eps");
  
//////////////////////////////////////////////
// now do some cross section plots with syst errors
  TCanvas* c1 = new TCanvas("c1", "V0 pT distribution");
  
  TLegend* sigmaPtLegend = new TLegend(0.65, 0.7, 0.9, 0.8);
  sigmaPtLegend->SetTextSize(0.035);
  sigmaPtLegend->SetFillColor(0);
  TH1F* pT7TeV_sigma = pT7TeV->Clone("pT7TeV_sigma");
  TH1F* pT7TeV_sigma2 = pT7TeV->Clone("pT7TeV_sigma2");
  TH1F* pT7TeV_sigma3 = pT7TeV->Clone("pT7TeV_sigma3");
  TH1F* pT7TeV_sigma4 = pT7TeV->Clone("pT7TeV_sigma4");

  for (int i = 1; i<=pT7TeV_sigma->GetNbinsX(); ++i) {
    float stat = pT7TeV->GetBinError(i);
    float systUncorOnly = pT7TeV->GetBinContent(i)*systUncorpt;
    float systUncor = sqrt( stat*stat + systUncorOnly*systUncorOnly );
    float systCor = pT7TeV->GetBinContent(i)*systCorpt;
    float systRatOnly = pT7TeV->GetBinContent(i)*systRatpt;
    float systRat = sqrt( systRatOnly*systRatOnly + stat*stat );
    pT7TeV_sigma->SetBinContent(i,pT7TeV->GetBinContent(i)/lumi7);  
    pT7TeV_sigma->SetBinError(i,stat/lumi7);  
    pT7TeV_sigma2->SetBinContent(i,pT7TeV->GetBinContent(i)/lumi7);  
    pT7TeV_sigma2->SetBinError(i,systUncor/lumi7);  
    pT7TeV_sigma3->SetBinContent(i,pT7TeV->GetBinContent(i)/lumi7);  
    pT7TeV_sigma3->SetBinError(i,systCor/lumi7);  
    pT7TeV_sigma3->SetBinContent(i,pT7TeV->GetBinContent(i)/lumi7);  
    pT7TeV_sigma3->SetBinError(i,systRat/lumi7);  
  }
  pT7TeV_sigma3->SetMaximum(1.2*pT7TeV->GetMaximum()/lumi7);
  pT7TeV_sigma3->SetMinimum(0);
  pT7TeV_sigma3->GetYaxis()->SetTitle("#sigma(pp#rightarrow #Xi^{-}X) (mb/GeV)");
  pT7TeV_sigma3->SetStats(0);
  pT7TeV_sigma3->SetLineColor(4);
  pT7TeV_sigma3->SetFillColor(38);
  pT7TeV_sigma3->Draw("e3");
  pT7TeV_sigma->SetMarkerColor(4);
  pT7TeV_sigma->SetLineColor(4);
  pT7TeV_sigma->SetMarkerStyle(20);
  sigmaPtLegend->AddEntry(pT7TeV_sigma, "#sqrt{s} = 7 TeV", "pl");
  pT7TeV_sigma->Draw("e1,same");
  pT7TeV_sigma2->SetMarkerColor(4);
  pT7TeV_sigma2->SetLineColor(4);
  pT7TeV_sigma2->SetMarkerStyle(20);
  pT7TeV_sigma2->Draw("same");
  
  TH1F* pT900_sigma = pT900->Clone("pT900_sigma");
  TH1F* pT900_sigma2 = pT900->Clone("pT900_sigma2");
  TH1F* pT900_sigma3 = pT900->Clone("pT900_sigma3");
  TH1F* pT900_sigma4 = pT900->Clone("pT900_sigma4");  
  
  for (int i = 1; i<=pT900_sigma->GetNbinsX(); ++i) {
    float stat = pT900->GetBinError(i);
    float systUncorOnly = pT900->GetBinContent(i)*systUncorpt;
    float systUncor = sqrt( stat*stat + systUncorOnly*systUncorOnly );
    float systCor = pT900->GetBinContent(i)*systCorpt;
    float systRatOnly = pT900->GetBinContent(i)*systRatpt;
    float systRat = sqrt( systRatOnly*systRatOnly + stat*stat );
    pT900_sigma->SetBinContent(i,pT900->GetBinContent(i)/lumi900);  
    pT900_sigma->SetBinError(i,stat/lumi900);  
    pT900_sigma2->SetBinContent(i,pT900->GetBinContent(i)/lumi900);  
    pT900_sigma2->SetBinError(i,systUncor/lumi900);  
    pT900_sigma3->SetBinContent(i,pT900->GetBinContent(i)/lumi900);  
    pT900_sigma3->SetBinError(i,systCor/lumi900);  
    pT900_sigma4->SetBinContent(i,pT900->GetBinContent(i)/lumi900);  
    pT900_sigma4->SetBinError(i,systRat/lumi900); 
  }
  pT900_sigma3->SetStats(0);
  pT900_sigma3->SetLineColor(2);
  pT900_sigma3->SetFillColor(kRed-9);
  pT900_sigma3->Draw("same,e3");
  pT900_sigma->SetMarkerColor(2);
  pT900_sigma->SetLineColor(2);
  pT900_sigma->SetMarkerStyle(21);
  sigmaPtLegend->AddEntry(pT900_sigma, "#sqrt{s} = 900 GeV", "pl");
  pT900_sigma->Draw("e1,same");
  pT900_sigma2->SetMarkerColor(2);
  pT900_sigma2->SetLineColor(2);
  pT900_sigma2->SetMarkerStyle(21);
  pT900_sigma2->Draw("same");
  
  sigmaPtLegend->Draw("f");
  cmsPre->Draw();
  c1->SaveAs(outputDir+"pt_sigma.eps");
  
  TLegend* sigmaYLegend = new TLegend(0.65, 0.38, 0.9, 0.48);
  sigmaYLegend->SetTextSize(0.035);
  sigmaYLegend->SetFillColor(0);
  TH1F* y7TeV_sigma = y7TeV->Clone("y7TeV_sigma");
  TH1F* y7TeV_sigma2 = y7TeV->Clone("y7TeV_sigma2");
  TH1F* y7TeV_sigma3 = y7TeV->Clone("y7TeV_sigma3");
  TH1F* y7TeV_sigma4 = y7TeV->Clone("y7TeV_sigma4");
  
  for (int i = 1; i<=y7TeV_sigma->GetNbinsX(); ++i) {   
    float stat = y7TeV->GetBinError(i);
    float systUncorOnly = y7TeV->GetBinContent(i)*systUncory;
    float systUncor = sqrt( stat*stat + systUncorOnly*systUncorOnly );
    float systCor = y7TeV->GetBinContent(i)*systCory7;
    float systRatOnly = y7TeV->GetBinContent(i)*systRaty7;
    float systRat = sqrt( systRatOnly*systRatOnly + stat*stat );
    y7TeV_sigma->SetBinContent(i,y7TeV->GetBinContent(i)/lumi7);  
    y7TeV_sigma->SetBinError(i,stat/lumi7);  
    y7TeV_sigma2->SetBinContent(i,y7TeV->GetBinContent(i)/lumi7);  
    y7TeV_sigma2->SetBinError(i,systUncor/lumi7);  
    y7TeV_sigma3->SetBinContent(i,y7TeV->GetBinContent(i)/lumi7);  
    y7TeV_sigma3->SetBinError(i,systCor/lumi7);  
    y7TeV_sigma4->SetBinContent(i,y7TeV->GetBinContent(i)/lumi7);  
    y7TeV_sigma4->SetBinError(i,systRat/lumi7);  
  }
  y7TeV_sigma3->SetMaximum(1.5*y7TeV->GetMaximum()/lumi7);
  y7TeV_sigma3->SetMinimum(0);
  y7TeV_sigma3->GetYaxis()->SetTitle("#sigma(pp#rightarrow #Xi^{-}X) (mb)");
  y7TeV_sigma3->GetXaxis()->SetTitle("#Xi^{-} rapidity");
  y7TeV_sigma3->SetStats(0);
  y7TeV_sigma3->SetLineColor(4);
  y7TeV_sigma3->SetFillColor(38);
  y7TeV_sigma3->Draw("e3");
  y7TeV_sigma->SetMarkerColor(4);
  y7TeV_sigma->SetLineColor(4);
  y7TeV_sigma->SetMarkerStyle(20);
  sigmaYLegend->AddEntry(y7TeV_sigma, "#sqrt{s} = 7 TeV", "pl");
  y7TeV_sigma->Draw("e1,same");
  y7TeV_sigma2->SetMarkerColor(4);
  y7TeV_sigma2->SetLineColor(4);
  y7TeV_sigma2->SetMarkerStyle(20);
  y7TeV_sigma2->Draw("same");

  TH1F* y900_sigma = y900->Clone("y900_sigma");
  TH1F* y900_sigma2 = y900->Clone("y900_sigma2");
  TH1F* y900_sigma3 = y900->Clone("y900_sigma3");
  TH1F* y900_sigma4 = y900->Clone("y900_sigma4");
  
  for (int i = 1; i<=y900_sigma->GetNbinsX(); ++i) {
    float stat = y900->GetBinError(i);
    float systUncorOnly = y900->GetBinContent(i)*systUncory;
    float systUncor = sqrt( stat*stat + systUncorOnly*systUncorOnly );
    float systCor = y900->GetBinContent(i)*systCory9;
    float systRatOnly = y900->GetBinContent(i)*systRaty9;
    float systRat = sqrt( systRatOnly*systRatOnly + stat*stat );
    y900_sigma->SetBinContent(i,y900->GetBinContent(i)/lumi900);  
    y900_sigma->SetBinError(i,stat/lumi900);  
    y900_sigma2->SetBinContent(i,y900->GetBinContent(i)/lumi900);  
    y900_sigma2->SetBinError(i,systUncor/lumi900);  
    y900_sigma3->SetBinContent(i,y900->GetBinContent(i)/lumi900);  
    y900_sigma3->SetBinError(i,systCor/lumi900);  
    y900_sigma4->SetBinContent(i,y900->GetBinContent(i)/lumi900);  
    y900_sigma4->SetBinError(i,systRat/lumi900);  
  }
  y900_sigma3->SetStats(0);
  y900_sigma3->SetLineColor(4);
  y900_sigma3->SetFillColor(kRed-9);
  y900_sigma3->Draw("same,e3");
  y900_sigma->SetMarkerColor(2);
  y900_sigma->SetLineColor(2);
  y900_sigma->SetMarkerStyle(21);
  sigmaYLegend->AddEntry(y900_sigma, "#sqrt{s} = 900 GeV", "pl");
  y900_sigma->Draw("e1,same");
  y900_sigma2->SetMarkerColor(2);
  y900_sigma2->SetLineColor(2);
  y900_sigma2->SetMarkerStyle(21);
  y900_sigma2->Draw("same");
  sigmaYLegend->Draw("f");
  //cmsPre->SetX(1.2);
  //cmsPre->SetY(1.1*y7TeV->GetMaximum()/lumi7);
  cmsPre->Draw();
  c1->SaveAs(outputDir+"y_sigma.eps");


//////////////////////////////////////////////////////////////////////
// next do some dN/dX plots, first pT

//first get 7TeV MC truth plots for plotting range reasons
  TH1F* pT7TeVMCGen_NSD = pT7TeVMCGen->Clone("pT7TeVMCGen_NSD");
  TH1F* pT7TeVMCGenp0_NSD = pT7TeVMCGenp0->Clone("pT7TeVMCGenp0_NSD");
  TH1F* pT7TeVMCGenp8_NSD = pT7TeVMCGenp8->Clone("pT7TeVMCGenp8_NSD");
  for (int i = 1; i<=pT7TeVMCGen_NSD->GetNbinsX(); ++i) {
    pT7TeVMCGen_NSD->SetBinContent(i,pT7TeVMCGen->GetBinContent(i)/nEvt7MC);  
    pT7TeVMCGenp0_NSD->SetBinContent(i,pT7TeVMCGenp0->GetBinContent(i)/nEvt7MCp0);  
    pT7TeVMCGenp8_NSD->SetBinContent(i,pT7TeVMCGenp8->GetBinContent(i)/nEvt7MCp8);  
  }
  pT7TeVMCGen_NSD->SetMaximum(1.2*pT7TeVMCGen_NSD->GetMaximum());
  pT7TeVMCGen_NSD->GetYaxis()->SetTitle("(1/N) dN / dp_{T}  (GeV/c)^{-1}");
  pT7TeVMCGen_NSD->GetXaxis()->SetTitle("#Xi^{-} p_{T} (GeV)");
  pT7TeVMCGen_NSD->SetStats(0);
  pT7TeVMCGen_NSD->SetLineColor(4);
  pT7TeVMCGenp0_NSD->SetLineColor(4);
  pT7TeVMCGenp0_NSD->SetLineStyle(3);
  pT7TeVMCGenp8_NSD->SetLineColor(4);
  pT7TeVMCGenp8_NSD->SetLineStyle(2);
  
  TLegend* nsdLegend = new TLegend(0.35, 0.60, 0.90, 0.85);
  nsdLegend->SetFillColor(0);
  nsdLegend->SetTextSize(0.032);
  TH1F* pT7TeV_NSD = pT7TeV->Clone("pT7TeV_NSD");
  TH1F* pT7TeV_NSD2 = pT7TeV->Clone("pT7TeV_NSD2");
  TH1F* pT7TeV_NSD3 = pT7TeV->Clone("pT7TeV_NSD3");
  TH1F* pT7TeV_NSD4 = pT7TeV->Clone("pT7TeV_NSD4");
  TH1F* pT7TeV_NSD5 = pT7TeV->Clone("pT7TeV_NSD5");
  for (int i = 1; i<=pT7TeV_NSD->GetNbinsX(); ++i) {
    float stat = pT7TeV->GetBinError(i);
    float systUncorOnly = pT7TeV->GetBinContent(i)*systUncorpt;
    float systUncor = sqrt( stat*stat + systUncorOnly*systUncorOnly );
    float systCor = pT7TeV->GetBinContent(i)*systCorpt;
    float systRatOnly = pT7TeV->GetBinContent(i)*systRatpt;
    float systRat = sqrt( systRatOnly*systRatOnly + stat*stat );
    pT7TeV_NSD->SetBinContent(i,pT7TeV->GetBinContent(i)/nEvt7);  
    pT7TeV_NSD->SetBinError(i,stat/nEvt7);  
    pT7TeV_NSD2->SetBinContent(i,pT7TeV->GetBinContent(i)/nEvt7);  
    pT7TeV_NSD2->SetBinError(i,systUncor/nEvt7);  
    pT7TeV_NSD3->SetBinContent(i,pT7TeV->GetBinContent(i)/nEvt7);  
    pT7TeV_NSD3->SetBinError(i,systCor/nEvt7);  
    pT7TeV_NSD4->SetBinContent(i,pT7TeV->GetBinContent(i)/nEvt7);  
    pT7TeV_NSD4->SetBinError(i,systRat/nEvt7);  
    pT7TeV_NSD5->SetBinContent(i,pT7TeV->GetBinContent(i)/nEvt7);  
    pT7TeV_NSD5->SetBinError(i,sqrt(systUncor*systUncor+systCor*systCor)/nEvt7);  
  }
  pT7TeV_NSD3->SetMaximum(1.2*pT7TeV_NSD3->GetMaximum());
  pT7TeV_NSD3->SetMinimum(0);
  pT7TeV_NSD3->GetYaxis()->SetTitle("(1/N) dN / dp_{T}  (GeV/c)^{-1}");
  pT7TeV_NSD3->GetXaxis()->SetTitle("#Xi^{-} p_{T} (GeV)");
  pT7TeV_NSD3->SetStats(0);
  pT7TeV_NSD3->SetLineColor(4);
  pT7TeV_NSD3->SetFillColor(38);
  pT7TeV_NSD3->Draw("e3");
  pT7TeV_NSD->SetMarkerColor(4);
  pT7TeV_NSD->SetLineColor(4);
  pT7TeV_NSD->SetMarkerStyle(20);
  nsdLegend->AddEntry(pT7TeV_NSD, "#sqrt{s} = 7 TeV", "pl");
  pT7TeV_NSD->Draw("e1,same");
  pT7TeV_NSD2->SetMarkerColor(4);
  pT7TeV_NSD2->SetLineColor(4);
  pT7TeV_NSD2->SetMarkerStyle(20);
  pT7TeV_NSD2->Draw("same");

  TH1F* pT900_NSD = pT900->Clone("pT900_NSD");
  TH1F* pT900_NSD2 = pT900->Clone("pT900_NSD2");
  TH1F* pT900_NSD3 = pT900->Clone("pT900_NSD3");
  TH1F* pT900_NSD4 = pT900->Clone("pT900_NSD4");  
  TH1F* pT900_NSD5 = pT900->Clone("pT900_NSD5");  
  for (int i = 1; i<=pT900_NSD->GetNbinsX(); ++i) {
    float stat = pT900->GetBinError(i);
    float systUncorOnly = pT900->GetBinContent(i)*systUncory;
    float systUncor = sqrt( stat*stat + systUncorOnly*systUncorOnly );
    float systCor = pT900->GetBinContent(i)*systCory9;
    float systRatOnly = pT900->GetBinContent(i)*systRaty9;
    float systRat = sqrt( systRatOnly*systRatOnly + stat*stat );
    pT900_NSD->SetBinContent(i,pT900->GetBinContent(i)/nEvt900);  
    pT900_NSD->SetBinError(i,stat/nEvt900);  
    pT900_NSD2->SetBinContent(i,pT900->GetBinContent(i)/nEvt900);  
    pT900_NSD2->SetBinError(i,systUncor/nEvt900); 
    pT900_NSD3->SetBinContent(i,pT900->GetBinContent(i)/nEvt900);  
    pT900_NSD3->SetBinError(i,systCor/nEvt900); 
    pT900_NSD4->SetBinContent(i,pT900->GetBinContent(i)/nEvt900);  
    pT900_NSD4->SetBinError(i,systRat/nEvt900); 
    pT900_NSD5->SetBinContent(i,pT900->GetBinContent(i)/nEvt900);  
    pT900_NSD5->SetBinError(i,sqrt(systUncor*systUncor+systCor*systCor)/nEvt900); 
  }
  pT900_NSD3->SetMarkerColor(2);
  pT900_NSD3->SetLineColor(2);
  pT900_NSD3->SetMarkerStyle(21);
  pT900_NSD3->SetStats(0);
  pT900_NSD3->SetFillColor(kRed-9);
  pT900_NSD3->Draw("same,e3");
  pT900_NSD->SetMarkerColor(2);
  pT900_NSD->SetLineColor(2);
  pT900_NSD->SetMarkerStyle(21);
  pT900_NSD->Draw("e1,same");
  pT900_NSD2->SetMarkerColor(2);
  pT900_NSD2->SetLineColor(2);
  pT900_NSD2->SetMarkerStyle(21);
  pT900_NSD2->Draw("same");
  nsdLegend->Draw("f");
  pT900_NSD->Draw("same");
  cmsPre->Draw();
 

  nsdLegend->AddEntry(pT7TeVMCGen_NSD, "     Pythia6 D6T", "l");
  pT7TeVMCGen_NSD->Draw("same,hist");
//  nsdLegend->AddEntry(pT7TeVMCGenp0_NSD, "     TeV P0", "l");
//  pT7TeVMCGenp0_NSD->Draw("same,hist");
//  nsdLegend->AddEntry(pT7TeVMCGenp8_NSD, "     Pythia8", "l");
//  pT7TeVMCGenp8_NSD->Draw("same,hist");
 
  TH1F* pT900MCGen_NSD = pT900MCGen->Clone("pT900MCGen_NSD");
  TH1F* pT900MCGenNoEn_NSD = pT900MCGenNoEn->Clone("pT900MCGenNoEn_NSD");
  TH1F* pT900MCGenD6T_NSD = pT900MCGenD6T->Clone("pT900MCGenD6T_NSD");
  TH1F* pT900MCGenp0_NSD = pT900MCGenp0->Clone("pT900MCGenp0_NSD");
  TH1F* pT900MCGenp8_NSD = pT900MCGenp8->Clone("pT900MCGenp8_NSD");
  for (int i = 1; i<=pT900MCGen_NSD->GetNbinsX(); ++i) {
    pT900MCGen_NSD->SetBinContent(i,pT900MCGen->GetBinContent(i)/nEvt900MC);  
    pT900MCGenNoEn_NSD->SetBinContent(i,pT900MCGenNoEn->GetBinContent(i)/nEvt900MCnoEn);  
    pT900MCGenD6T_NSD->SetBinContent(i,pT900MCGenD6T->GetBinContent(i)/nEvt900MCD6T);  
    pT900MCGenp0_NSD->SetBinContent(i,pT900MCGenp0->GetBinContent(i)/nEvt900MCp0);  
    pT900MCGenp8_NSD->SetBinContent(i,pT900MCGenp8->GetBinContent(i)/nEvt900MCp8);  
  }
  nsdLegend->AddEntry(pT900_NSD, "  #sqrt{s} = 900 GeV", "pl");
  pT900MCGenNoEn_NSD->SetMaximum(1.2*pT900MCGenNoEn_NSD->GetMaximum());
  pT900MCGenNoEn_NSD->GetYaxis()->SetTitle("(1/N) dN / dp_{T}  (GeV/c)^{-1}");
  pT900MCGenNoEn_NSD->GetXaxis()->SetTitle("#Xi^{-} p_{T} (GeV)");
  pT900MCGenNoEn_NSD->SetStats(0);
  pT900MCGenNoEn_NSD->SetLineColor(2);
  nsdLegend->AddEntry(pT900MCGenNoEn_NSD, "     Pythia6 D6T", "l");
  pT900MCGenNoEn_NSD->Draw("same,hist");
  
//  pT900MCGenD6T_NSD->SetLineColor(2);
//  pT900MCGenD6T_NSD->SetLineStyle(1);
//  nsdLegend->AddEntry(pT900MCGenD6T_NSD, "     Pythia6 D6T", "l");
//  pT900MCGenD6T_NSD->Draw("same,hist");
//  pT900MCGenp0_NSD->SetLineColor(2);
//  pT900MCGenp0_NSD->SetLineStyle(3);
//  nsdLegend->AddEntry(pT900MCGenp0_NSD, "     Pythia6 P0", "l");
//  pT900MCGenp0_NSD->Draw("same,hist");
//  pT900MCGenp8_NSD->SetLineColor(2);
//  pT900MCGenp8_NSD->SetLineStyle(2);
//  nsdLegend->AddEntry(pT900MCGenp8_NSD, "     Pythia8", "l");
//  pT900MCGenp8_NSD->Draw("same,hist");

  c1->SaveAs(outputDir+"pt_NSD_dataMC.eps");

  TF1 *tsfunc7D = new TF1("tsfunc7D",tsallis,0.0,6.0,3);
  tsfunc7D->SetParNames("Amplitude","T","n");
  tsfunc7D->SetParameters(0.12,0.31,8.4); 

  pT7TeV_NSD->Fit("tsfunc7D","IEM0");
  tsfunc7D->SetNpx(50000);
  tsfunc7D->SetLineColor(4);
  tsfunc7D->GetXaxis()->SetTitle("#Xi^{-} p_{T} (GeV)");
  tsfunc7D->GetYaxis()->SetTitle("(1/N) dN / dp_{T}  (GeV/c)^{-1}");
  tsfunc7D->SetMaximum(pT7TeV_NSD3->GetMaximum());
  tsfunc7D->Draw();

  outTextFile << "Tsallis paramter values for 7 TeV data" << endl;
  outTextFile << "==================================================" << endl;
  outTextFile << "$\\Xi^-$ & " << tsfunc7D->GetParameter(1) << " \\pm " << tsfunc7D->GetParError(1) << 
  " & " << tsfunc7D->GetParameter(2) << " \\pm " << tsfunc7D->GetParError(2) << " & " << endl;

  TF1 *tsmean7D = new TF1("tsmean7D",tsallismean,0.0,6.0,3);
  tsmean7D->SetParNames("Amplitude","A","n");

  tsmean7D->SetParameters(0.12.,0.96,8.4); 
//  pT7TeV_NSD->Fit("tsmean7D","IEM0");
   pT7TeV_NSD->Fit("tsmean7D","IE0");
  
  outTextFile << "Tsallis mean pT (mass=0) = " << tsmean7D->GetParameter(1) << "\pm" << 
  tsmean7D->GetParError(1)<< " rel error: " << tsmean7D->GetParError(1)/tsmean7D->GetParameter(1) << endl;

  pT7TeV_NSD3->Draw("e3,same");
  pT7TeV_NSD->Draw("e1,same");
  pT7TeV_NSD2->Draw("same");
  tsfunc7D->Draw("same");

  tsfunc7D->SetParameters(0.136,0.30,8.0); 
  pT7TeV_NSD5->Fit("tsfunc7D","IEM0");
  outTextFile << "Tsallis parameter values with systematics: " << endl;
  outTextFile << "$\Xi^-$ & " << tsfunc7D->GetParameter(1) << " \\pm " << tsfunc7D->GetParError(1) << 
  " & " << tsfunc7D->GetParameter(2) << " \\pm " << tsfunc7D->GetParError(2) << " & " << endl;

  tsmean7D->SetParameters(0.136.,0.95,8.0); 
//  pT7TeV_NSD5->Fit("tsmean7D","IEM0");
    pT7TeV_NSD5->Fit("tsmean7D","IE0");
  
  outTextFile << "Tsallis mean pT (mass=0) = " << tsmean7D->GetParameter(1) << "\pm" << 
  tsmean7D->GetParError(1)<< " rel error: " << tsmean7D->GetParError(1)/tsmean7D->GetParameter(1) << endl;
  outTextFile << endl;


  
  TF1 *tsfunc = new TF1("tsfunc",tsallis,0.0,6.0,3);
  tsfunc->SetParNames("Amplitude","T","n");
  tsfunc->SetParameters(1000000.0,0.2,7.2);  

  pT7TeV->Fit("tsfunc","IEM0");
  float tsallisIntFit = tsallisInt(mXi, tsfunc->GetParameter(2), tsfunc->GetParameter(1), tsfunc->GetParameter(0) );
  outTextFile << "Integral over pT from formula = " << tsallisIntFit << ", sigma = " << tsallisIntFit/lumi7 << ", dXi/dNSD = " <<
  tsallisIntFit/nEvt7 << endl;

  float tsallisIntCur = tsfunc->Integral(0.0,10000.0);
  outTextFile << "Integral from pT curve function = " << tsallisIntCur << ", sigma = " << tsallisIntCur/lumi7 << ", dXi/dNSD = " <<
  tsallisIntCur/nEvt7 << endl;

  tsallisIntCur = tsfunc->Integral(0.0,6.0);
  outTextFile << "Integral from pT from 0.0 to 6.0 = " << tsallisIntCur << ", sigma = " << tsallisIntCur/lumi7 << ", dXi/dNSD = " <<
  tsallisIntCur/nEvt7 << endl;
  //these values are scaled by bin width, but it's constant bin width so the ratios are fine.
  float npT7TeVlt6 = pT7TeVMCGen->Integral(1,60);
  float npT7TeVgt6 = pT7TeVMCGen->Integral(61,101);
  float fracpT7TeV = npT7TeVlt6/(npT7TeVlt6+npT7TeVgt6);
  outTextFile << "From MC truth fraction from 0.0 to 6.0 = " << fracpT7TeV << "  two values are " << npT7TeVlt6 << 
  ", " << npT7TeVgt6 << endl;

  float totalpT7TeV = 0.0;
  float totalpT7TeVEsq = 0.0;
  for (int i=1; i<=pT7TeVFit->GetNbinsX(); i++) {
    float scaled = (pT7TeVFit->GetBinContent(i))/(effpT7TeV->GetBinContent(i));
    totalpT7TeV += scaled;
    totalpT7TeVEsq += ((pT7TeVFit->GetBinError(i))/(effpT7TeV->GetBinContent(i))*(pT7TeVFit->GetBinError(i))/(effpT7TeV->GetBinContent(i)));
  }
  outTextFile << "Sum of corrected yields in pT (0.0 to 6.0) = " << totalpT7TeV << "  scaled to full range = " << totalpT7TeV/fracpT7TeV  <<
  "  sigma = " <<   totalpT7TeV/(lumi7*fracpT7TeV) << "   dXi/dNSD = " << totalpT7TeV/(nEvt7*fracpT7TeV) << endl;
  outTextFile << "  integral 0-6.0 = " << tsfunc->Integral(0.0,6.0) << 
                 "  integral 6.0-10000 = " << tsfunc->Integral(6.0,10000.0) << endl;
  float totalpT7TeV_fitExtrap = totalpT7TeV + tsfunc->Integral(6.0,10000.0);
  outTextFile << "Sum of corrected yields in pT (0.0 to 6.0) + fit integral 6+ = " << totalpT7TeV_fitExtrap << 
  "  sigma = " <<   totalpT7TeV_fitExtrap/lumi7 << "   dXi/dNSD = " << totalpT7TeV_fitExtrap/nEvt7 << endl;
  outTextFile << "fractional error on sum of yields from bins = " << sqrt(totalpT7TeVEsq)/totalpT7TeV << endl;
  outTextFile << "  VALUE TO USE: " << totalpT7TeV_fitExtrap << " +- " << (sqrt(totalpT7TeVEsq)/totalpT7TeV)*totalpT7TeV_fitExtrap << " +- " << 
                 sqrt(systUncory*systUncory + systCory7*systCory7)*totalpT7TeV_fitExtrap << 
		 "   Xi/event =  " << totalpT7TeV_fitExtrap/nEvt7 << " +- " << (sqrt(totalpT7TeVEsq)/totalpT7TeV)*totalpT7TeV_fitExtrap/nEvt7 << " +- " << 
                 sqrt(systUncory*systUncory + systCory7*systCory7)*totalpT7TeV_fitExtrap/nEvt7 << endl;


      
  float totalct7TeV = 0.0;
  for (int i=1; i<=ct7TeVFit->GetNbinsX(); i++) {
    float scaled = (ct7TeVFit->GetBinContent(i))/(effCtau7TeV->GetBinContent(i));
    totalct7TeV += scaled;
  }
  outTextFile << "Sum of corrected yields in ct (0.0 to 10.0) = " << totalct7TeV << 
  "  sigma = " <<   totalct7TeV/lumi7 << "   dXi/dNSD = " << totalct7TeV/nEvt7 << endl;

  float totalY7TeV = 0.0;
  for (int i=1; i<=y7TeVFit->GetNbinsX(); i++) {
    float scaled = (y7TeVFit->GetBinContent(i))/(effY7TeV->GetBinContent(i));
    totalY7TeV += scaled;
  }
  outTextFile << "Sum of corrected yields in Y = " << totalY7TeV << 
  "  sigma = " <<   totalY7TeV/lumi7 << "   dXi/dNSD = " << totalY7TeV/nEvt7 << endl;


   



  TF1 *tsfunc9D = new TF1("tsfunc9D",tsallis,0.0,6.0,3);
  tsfunc9D->SetParNames("Amplitude","T","n");
  tsfunc9D->SetParameters(0.096,0.23,8.6); 

  pT900_NSD->Fit("tsfunc9D","IEM0");
  tsfunc9D->SetNpx(50000);
  tsfunc9D->SetLineColor(2);
  tsfunc9D->Draw("same");
  pT900_NSD3->Draw("e3,same");
  pT900_NSD->Draw("e1,same");
  pT900_NSD2->Draw("same");
  tsfunc9D->Draw("same");
  
  cmsPre->Draw();
  TLegend* nsdDFLegend = new TLegend(0.65, 0.65, 0.9, 0.8);
  nsdDFLegend->SetTextSize(0.035);
  nsdDFLegend->SetFillColor(0);
  nsdDFLegend->AddEntry(pT7TeV_NSD, "#sqrt{s} = 7 TeV", "lp");
  nsdDFLegend->AddEntry(pT900_NSD, "#sqrt{s} = 900 GeV", "lp");
  nsdDFLegend->Draw("same");
  
  c1->SaveAs(outputDir+"pt_NSD_dataFits.eps");
  outTextFile << endl; outTextFile << endl; outTextFile << endl; 
  outTextFile << "Tsallis paramter values for 900 GeV data" << endl;
  outTextFile << "====================================================" << endl;
  outTextFile << "$\Xi^-$ & " << tsfunc9D->GetParameter(1) << " \\pm " << tsfunc9D->GetParError(1) << 
  " & " << tsfunc9D->GetParameter(2) << " \\pm " << tsfunc9D->GetParError(2) << " & " << endl;


  TF1 *tsmean9D = new TF1("tsmean9D",tsallismean,0.0,6.0,3);
  tsmean9D->SetParNames("Amplitude","A","n");
  tsmean9D->SetParameters(0.096,0.72,8.6); 
  pT900_NSD->Fit("tsmean9D","IEM0");

  outTextFile << "Tsallis mean pT (mass=0) = " << tsmean9D->GetParameter(1) << "\pm" << 
  tsmean9D->GetParError(1)<< " rel error: " << tsmean9D->GetParError(1)/tsmean9D->GetParameter(1) << endl;

  outTextFile << "Next do Tsallis fit with systematics:" << endl;
  pT900_NSD5->Fit("tsfunc9D","IEM0");
  outTextFile << "$\Xi^-$ & " << tsfunc9D->GetParameter(1) << " \\pm " << tsfunc9D->GetParError(1) << 
  " & " << tsfunc9D->GetParameter(2) << " \\pm " << tsfunc9D->GetParError(2) << " & " << endl;



  pT900_NSD5->Fit("tsmean9D","IEM0");
  outTextFile << "Tsallis mean pT (mass=0) = " << tsmean9D->GetParameter(1) << "\pm" << 
  tsmean9D->GetParError(1)<< " rel error: " << tsmean9D->GetParError(1)/tsmean9D->GetParameter(1) << endl;
  outTextFile << endl;



  pT900->Fit("tsfunc","IEM0");
  tsallisIntFit = tsallisInt(mXi, tsfunc->GetParameter(2), tsfunc->GetParameter(1), tsfunc->GetParameter(0) );
  outTextFile << "Integral over pT from formula = " << tsallisIntFit << ", sigma = " << tsallisIntFit/lumi900 << ", dXi/dNSD = " <<
  tsallisIntFit/nEvt900 << endl;
  tsallisIntCur = tsfunc->Integral(0.0,10000.0);
  outTextFile << "Integral from pT curve function = " << tsallisIntCur << ", sigma = " << tsallisIntCur/lumi900 << ", dXi/dNSD = " <<
  tsallisIntCur/nEvt900 << endl;
  tsallisIntCur = tsfunc->Integral(0.0,6.0);
  outTextFile << "Integral from pT from 0.0 to 6.0 = " << tsallisIntCur << ", sigma = " << tsallisIntCur/lumi900 << ", dXi/dNSD = " <<
  tsallisIntCur/nEvt900 << endl;
  //these values are scaled by bin width, but it's constant bin width so the ratios are fine.
  float npT900lt6 = pT900MCGenNoEn->Integral(1,60);
  float npT900gt6 = pT900MCGenNoEn->Integral(61,101);
  float fracpT900 = npT900lt6/(npT900lt6+npT900gt6);
  outTextFile << "From MC truth fraction from 0.0 to 6.0 = " << fracpT900 << "  two values are " << npT900lt6 << 
  ", " << npT900gt6 << endl;

  float totalpT900 = 0.0;
  float totalpT900Esq = 0.0;
  for (int i=1; i<=pT900Fit->GetNbinsX(); i++) {
    float scaled = (pT900Fit->GetBinContent(i))/(effpT900->GetBinContent(i));
    totalpT900 += scaled;
    totalpT900Esq += ((pT900Fit->GetBinError(i))/(effpT900->GetBinContent(i))*(pT900Fit->GetBinError(i))/(effpT900->GetBinContent(i)));
  }
  outTextFile << "Sum of corrected yields in pT (0.0 to 6.0) = " << totalpT900 << " scaled to full range = " << totalpT900/fracpT900 <<
  "  sigma = " <<   totalpT900/(lumi900*fracpT900) << "   dXi/dNSD = " << totalpT900/(nEvt900*fracpT900) << endl;
  outTextFile << "  integral 0-6.0 = " << tsfunc->Integral(0.0,6.0) <<
                 "  integral 6.0-10000 = " << tsfunc->Integral(6.0,10000.0) << endl;
  float totalpT900_fitExtrap = totalpT900 + tsfunc->Integral(6.0,10000.0);
  outTextFile << "Sum of corrected yields in pT (0.0 to 6.0) + fit integral 6+ = " << totalpT900_fitExtrap << 
  "  sigma = " <<   totalpT900_fitExtrap/lumi900 << "   dXi/dNSD = " << totalpT900_fitExtrap/nEvt900 << endl; 
  outTextFile << "fractional error on fit sum = " << sqrt(totalpT900Esq)/totalpT900 << endl;

  outTextFile << "  VALUE TO USE: " << totalpT900_fitExtrap << " +- " << (sqrt(totalpT900Esq)/totalpT900)*totalpT900_fitExtrap << " +- " << 
                 sqrt(systUncory*systUncory + systCory9*systCory9)*totalpT900_fitExtrap <<  
		 "   Xi/event = " << totalpT900_fitExtrap/nEvt900 << " +- " << (sqrt(totalpT900Esq)/totalpT900)*totalpT900_fitExtrap/nEvt900 << " +- " << 
                 sqrt(systUncory*systUncory + systCory9*systCory9)*totalpT900_fitExtrap/nEvt900 <<endl;
  outTextFile << "  FOR TEX TABLE:  " <<  "   $\Xi^-$ & " << totalpT900_fitExtrap/nEvt900 << " & $\\pm$ & " << 
                 (sqrt(totalpT900Esq)/totalpT900)*totalpT900_fitExtrap/nEvt900 << " & $\\pm$ & " << 
                 sqrt(systUncory*systUncory + systCory9*systCory9)*totalpT900_fitExtrap/nEvt900 <<
		 "   &  " << totalpT7TeV_fitExtrap/nEvt7 << " & $\\pm$ " << (sqrt(totalpT7TeVEsq)/totalpT7TeV)*totalpT7TeV_fitExtrap/nEvt7 << " & $\\pm$ " << 
                 sqrt(systUncory*systUncory + systCory7*systCory7)*totalpT7TeV_fitExtrap/nEvt7 << endl;

  float totalct900 = 0.0;
  for (int i=1; i<=ct900Fit->GetNbinsX(); i++) {
    float scaled = (ct900Fit->GetBinContent(i))/(effCtau900->GetBinContent(i));
    totalct900 += scaled;
  }
  outTextFile << "Sum of corrected yields in ct (0.0 to 10.0) = " << totalct900 << 
  "  sigma = " <<   totalct900/lumi900 << "   dXi/dNSD = " << totalct900/nEvt900 << endl;
  
  float totalY900 = 0.0;
  for (int i=1; i<=y900Fit->GetNbinsX(); i++) {
    float scaled = (y900Fit->GetBinContent(i))/(effY900->GetBinContent(i));
    totalY900 += scaled;
  }
  outTextFile << "Sum of corrected yields in Y = " << totalY900 << 
  "  sigma = " <<   totalY900/lumi900 << "   dXi/dNSD = " << totalY900/nEvt900 << endl;


  
  TF1 *tsfunc7M = new TF1("tsfunc7M",tsallis,0.0,6.0,3);
  tsfunc7M->SetParNames("Amplitude","T","n");
  tsfunc7M->SetParameters(0.065,0.21,4.7); 
    
//  pT7TeVMCGen_NSD->Fit("tsfunc7M","IEM0");
  pT7TeVMCGen_NSD->Fit("tsfunc7M","I0");
  tsfunc7M->SetNpx(50000);
  tsfunc7M->SetLineColor(4);
  tsfunc7M->SetLineStyle(2);
  tsfunc7M->GetXaxis()->SetTitle("#Xi^{-} p_{T} (GeV)");
  tsfunc7M->GetYaxis()->SetTitle("(1/N) dN / dp_{T}  (GeV/c)^{-1}");
  tsfunc7M->Draw();
  pT7TeVMCGen_NSD->Draw("same,hist");  

  outTextFile << endl; outTextFile << endl; outTextFile << endl; 
  outTextFile << "Tsallis paramter values for 7 TeV D6T MC" << endl;
  outTextFile << "================================================" << endl;
  outTextFile << "$\Xi^-$ & " << tsfunc7M->GetParameter(1) << " \\pm " << tsfunc7M->GetParError(1) << 
  " & " << tsfunc7M->GetParameter(2) << " \\pm " << tsfunc7M->GetParError(2) << " & " << endl;


  TF1 *tsmean7M = new TF1("tsmean7M",tsallismean,0.0,6.0,3);
  tsmean7M->SetParNames("Amplitude","A","n");
  tsmean7M->SetParameters(0.065.,1.18,4.7); 
  //pT7TeVMCGen_NSD->Fit("tsmean7M","IEM0");
  pT7TeVMCGen_NSD->Fit("tsmean7M","I0");
  
  outTextFile << "Tsallis mean pT (mass=0) = " << tsmean7M->GetParameter(1) << "\pm" << 
  tsmean7M->GetParError(1)<< " rel error: " << tsmean7M->GetParError(1)/tsmean7M->GetParameter(1) << endl;
  outTextFile << "mean pT directly from generated events = " << pT7TeVMCGen->GetMean() << "\\pm" << pT7TeVMCGen->GetMeanError()<< endl;
  outTextFile << "total entries in MC = " << y7TeVMCGen->GetEntries() <<
  //" or integral = " << y7TeVMCGen->Integral() << 
  //"  sigma = " << y7TeVMCGen->GetEntries()/lumi7 << 
  "  dXi/dNSD = " << y7TeVMCGen->GetEntries()/nEvt7MC << endl;
  
  TH1F* pT7TeVMC_NSD = pT7TeVMC->Clone("pT7TeVMC_NSD");
  for (int i = 1; i<=pT7TeVMC_NSD->GetNbinsX(); ++i) {
    pT7TeVMC_NSD->SetBinContent(i,pT7TeVMC->GetBinContent(i)/nEvt7MC);  
    pT7TeVMC_NSD->SetBinError(i,pT7TeVMC->GetBinError(i)/nEvt7MC);  
  }
  pT7TeVMC_NSD->SetMaximum(1.2*pT7TeV->GetMaximum()/nEvt7);
  pT7TeVMC_NSD->GetYaxis()->SetTitle("(1/N) dN / dp_{T}  (GeV/c)^{-1}");
  pT7TeVMC_NSD->GetXaxis()->SetTitle("#Xi^{-} p_{T} (GeV)");
  pT7TeVMC_NSD->SetStats(0);
  pT7TeVMC_NSD->SetMarkerColor(4);
  pT7TeVMC_NSD->SetLineColor(4);
  pT7TeVMC_NSD->SetMarkerStyle(24);
  nsdLegend->AddEntry(pT7TeVMC_NSD, "#sqrt{s} = 7 TeV D6T MC Reco", "pl");
  pT7TeVMC_NSD->Draw("same");  


  
  TF1 *tsfunc9M = new TF1("tsfunc9M",tsallis,0.0,6.0,3);
  tsfunc9M->SetParNames("Amplitude","T","n");
  tsfunc9M->SetParameters(0.051,0.19,6.9);  

//  pT900MCGen_NSD->Fit("tsfunc9M","IEM0");
  pT900MCGen_NSD->Fit("tsfunc9M","I0");
  tsfunc9M->SetNpx(50000);
  tsfunc9M->SetLineColor(2);
  tsfunc9M->SetLineStyle(2);
  tsfunc9M->Draw("same");  
  pT900MCGen_NSD->SetLineColor(2);
  pT900MCGen_NSD->Draw("same,hist");

  
  outTextFile << endl; outTextFile << endl; outTextFile << endl; 
  outTextFile << "Tsallis paramter values for 900 GeV MC" << endl;
  outTextFile << "=================================================" << endl;
  outTextFile << "$\Xi^-$ & " << tsfunc9M->GetParameter(1) << " \\pm " << tsfunc9M->GetParError(1) << 
  " & " << tsfunc9M->GetParameter(2) << " \\pm " << tsfunc9M->GetParError(2) << " & " << endl;


  TF1 *tsmean9M = new TF1("tsmean9M",tsallismean,0.0,6.0,3);
  tsmean9M->SetParNames("Amplitude","A","n");
  tsmean9M->SetParameters(0.051,0.66,6.9); 

//  pT900MCGen_NSD->Fit("tsmean9M","IEM0");
  pT900MCGen_NSD->Fit("tsmean9M","I0");
  outTextFile << "Tsallis mean pT (mass=0) = " << tsmean9M->GetParameter(1) << "\pm" << 
  tsmean9M->GetParError(1)<< " rel error: " << tsmean9M->GetParError(1)/tsmean9M->GetParameter(1) << endl;
  outTextFile << "mean pT directly from generated events = " << pT900MCGenNoEn->GetMean() << "\\pm" << pT900MCGenNoEn->GetMeanError()<< endl;
  outTextFile << "total entries in (no enhanced) MC = " << y900MCGenNoEn->GetEntries() <<
  "  dXi/dNSD = " << y900MCGenNoEn->GetEntries()/nEvt900MCnoEn << endl;
  outTextFile << endl; outTextFile << endl; outTextFile << endl; 

  
  TH1F* pT900MC_NSD = pT900MC->Clone("pT900MC_NSD");
  for (int i = 1; i<=pT900MC_NSD->GetNbinsX(); ++i) {
    pT900MC_NSD->SetBinContent(i,pT900MC->GetBinContent(i)/nEvt900MC);  
    pT900MC_NSD->SetBinError(i,pT900MC->GetBinError(i)/nEvt900MC);  
  }
  pT900MC_NSD->SetMaximum(1.2*pT900MC->GetMaximum()/nEvt900MC);
  pT900MC_NSD->GetYaxis()->SetTitle("(1/N) dN / dp_{T}  (GeV/c)^{-1}");
  pT900MC_NSD->GetXaxis()->SetTitle("#Xi^{-} p_{T} (GeV)");
  pT900MC_NSD->SetStats(0);
  pT900MC_NSD->SetMarkerColor(2);
  pT900MC_NSD->SetLineColor(2);
  pT900MC_NSD->SetMarkerStyle(25);
  nsdLegend->AddEntry(pT900MC_NSD, "#sqrt{s} = 900 GeV MC Reco", "pl");
  pT900MC_NSD->Draw("same");
  
  cmsPre->Draw();
  TLegend* nsdMFLegend = new TLegend(0.45, 0.6, 0.9, 0.8);
  nsdMFLegend->SetTextSize(0.035);
  nsdMFLegend->SetFillColor(0);
  nsdMFLegend->AddEntry(pT7TeVMCGen_NSD, "#sqrt{s} = 7 TeV D6T MC Gen", "l");
  nsdMFLegend->AddEntry(pT7TeVMC_NSD, "#sqrt{s} = 7 TeV D6T MC Reco", "lp");
  nsdMFLegend->AddEntry(pT900MCGen_NSD, "#sqrt{s} = 900 GeV MC Gen", "l");
  nsdMFLegend->AddEntry(pT900MC_NSD, "#sqrt{s} = 900 GeV MC Reco", "lp");
  nsdMFLegend->Draw("same");
  
  c1->SaveAs(outputDir+"pt_NSD_MCFits.eps");

  for (int i = 1; i<=pt7TeVMCforRatio->GetNbinsX(); ++i){
    pt7TeVMCforRatio->SetBinContent(i,pt7TeVMCforRatio->GetBinContent(i)/(pt7TeVMCforRatio->GetBinWidth(i)*nEvt7MC));
    pt7TeVMCforRatio->SetBinError(i,pt7TeVMCforRatio->GetBinError(i)/(pt7TeVMCforRatio->GetBinWidth(i)*nEvt7MC));
  }
  
  TH1F* pt7TeVMCforRatio2 = pt7TeVMCforRatio->Clone("pt7TeVMCforRatio2");
  pt7TeVMCforRatio->Divide(pt7TeVMCforRatio, pT7TeV_NSD);
  pt7TeVMCforRatio->GetXaxis()->SetTitle("#Xi^{-} p_{T} (GeV)");
  pt7TeVMCforRatio->GetYaxis()->SetTitle("Ratio MC/data");
  pt7TeVMCforRatio->SetLineColor(4);
  pt7TeVMCforRatio->SetMarkerColor(4);
  pt7TeVMCforRatio->SetMinimum(0);
  pt7TeVMCforRatio->SetMaximum(1.5);
  pt7TeVMCforRatio->Draw("e1");
  pt7TeVMCforRatio2->Divide(pt7TeVMCforRatio2, pT7TeV_NSD4);
  pt7TeVMCforRatio2->SetLineColor(4);
  pt7TeVMCforRatio2->SetMarkerColor(4);
  pt7TeVMCforRatio2->Draw("same");

  for (int i = 1; i<=pt900MCforRatio->GetNbinsX(); ++i){
    pt900MCforRatio->SetBinContent(i,pt900MCforRatio->GetBinContent(i)/(pt900MCforRatio->GetBinWidth(i)*nEvt900MC));
    pt900MCforRatio->SetBinError(i,pt900MCforRatio->GetBinError(i)/(pt900MCforRatio->GetBinWidth(i)*nEvt900MC));
  }

  TH1F* pt900MCforRatio2 = pt900MCforRatio->Clone("pt900MCforRatio2");  
  pt900MCforRatio->Divide(pt900MCforRatio, pT900_NSD);
  pt900MCforRatio->SetLineColor(2);
  pt900MCforRatio->SetMarkerColor(2);
  pt900MCforRatio->Draw("e1,same");
  pt900MCforRatio2->Divide(pt900MCforRatio2, pT900_NSD4);
  pt900MCforRatio2->SetLineColor(2);
  pt900MCforRatio2->SetMarkerColor(2);
  pt900MCforRatio2->Draw("same");
  
  cmsPreUL->Draw();
  TLegend* nsdRatLegend = new TLegend(0.35, 0.75, 0.6, 0.85);
  nsdRatLegend->SetTextSize(0.035);
  nsdRatLegend->SetFillColor(0);
  nsdRatLegend->AddEntry(pt7TeVMCforRatio, "#sqrt{s} = 7 TeV", "l");
  nsdRatLegend->AddEntry(pt900MCforRatio, "#sqrt{s} = 900 GeV", "l");
  nsdRatLegend->Draw("same");
  
  c1->SaveAs(outputDir+"pt_NSD_ratio.eps");
  
  for (int i = 1; i<=pt900MCforFitVal->GetNbinsX(); ++i){
    pt900MCforFitVal->SetBinContent(i,pt900MCforFitVal->GetBinContent(i)/(pt900MCforFitVal->GetBinWidth(i)));
    pt900MCforFitVal->SetBinError(i,pt900MCforFitVal->GetBinError(i)/(pt900MCforFitVal->GetBinWidth(i)));
  }
  
  pt900MCforFitVal->Divide(pT900MC, pt900MCforFitVal);
  pt900MCforFitVal->GetXaxis()->SetTitle("#Xi^{-} p_{T} (GeV)");
  pt900MCforFitVal->GetYaxis()->SetTitle("Ratio MC fit/MC truth");
  pt900MCforFitVal->SetLineColor(2);
  pt900MCforFitVal->SetMarkerColor(2);
  pt900MCforFitVal->Draw();

  for (int i = 1; i<=pt7TeVMCforFitVal->GetNbinsX(); ++i){
    pt7TeVMCforFitVal->SetBinContent(i,pt7TeVMCforFitVal->GetBinContent(i)/pt7TeVMCforFitVal->GetBinWidth(i));
    pt7TeVMCforFitVal->SetBinError(i,pt7TeVMCforFitVal->GetBinError(i)/(pt7TeVMCforFitVal->GetBinWidth(i)));
  }
  
  pt7TeVMCforFitVal->Divide(pT7TeVMC, pt7TeVMCforFitVal);
  pt7TeVMCforFitVal->SetLineColor(4);
  pt7TeVMCforFitVal->SetMarkerColor(4);
  pt7TeVMCforFitVal->Draw("same");
  
  cmsPreUL->Draw();
  TLegend* nsdRatLegend = new TLegend(0.35, 0.75, 0.6, 0.85);
  nsdRatLegend->SetTextSize(0.035);
  nsdRatLegend->SetFillColor(0);
  nsdRatLegend->AddEntry(pt7TeVMCforFitVal, "#sqrt{s} = 7 TeV", "l");
  nsdRatLegend->AddEntry(pt900MCforFitVal, "#sqrt{s} = 900 GeV", "l");
  nsdRatLegend->Draw("same");
  
  c1->SaveAs(outputDir+"pt_fitval_ratio.eps");
  
  /////////////////////////////////////////////////////////////////////////
  //  Moving on to rapidity
  
  
  TLegend* nsdLegend = new TLegend(0.62, 0.67, 0.95, 0.92);
  nsdLegend->SetFillColor(0);
  nsdLegend->SetTextSize(0.032);
  TH1F* y7TeV_NSD = y7TeV->Clone("y7TeV_NSD");
  TH1F* y7TeV_NSD2 = y7TeV->Clone("y7TeV_NSD2");
  TH1F* y7TeV_NSD3 = y7TeV->Clone("y7TeV_NSD3");
  TH1F* y7TeV_NSD4 = y7TeV->Clone("y7TeV_NSD4");
  
  for (int i = 1; i<=y7TeV_NSD->GetNbinsX(); ++i) {
    float stat = y7TeV->GetBinError(i);
    float systUncorOnly = y7TeV->GetBinContent(i)*systUncory;
    float systUncor = sqrt( stat*stat + systUncorOnly*systUncorOnly );
    float systCor = y7TeV->GetBinContent(i)*systCory7;
    float systRatOnly = y7TeV->GetBinContent(i)*systRaty7;
    float systRat = sqrt( systRatOnly*systRatOnly + stat*stat );
    y7TeV_NSD->SetBinContent(i,y7TeV->GetBinContent(i)/nEvt7);  
    y7TeV_NSD->SetBinError(i,stat/nEvt7);  
    y7TeV_NSD2->SetBinContent(i,y7TeV->GetBinContent(i)/nEvt7);  
    y7TeV_NSD2->SetBinError(i,systUncor/nEvt7);  
    y7TeV_NSD3->SetBinContent(i,y7TeV->GetBinContent(i)/nEvt7);  
    y7TeV_NSD3->SetBinError(i,systCor/nEvt7);  
    y7TeV_NSD4->SetBinContent(i,y7TeV->GetBinContent(i)/nEvt7);  
    y7TeV_NSD4->SetBinError(i,systRat/nEvt7);  
  }

  y7TeV_NSD3->SetMaximum(1.7*y7TeV->GetMaximum()/nEvt7);
  y7TeV_NSD3->SetMinimum(0);
  y7TeV_NSD3->GetYaxis()->SetTitle("(1/N) dN / dy ");
  y7TeV_NSD3->GetXaxis()->SetTitle("#Xi^{-} rapidity");
  y7TeV_NSD3->SetStats(0);
  y7TeV_NSD3->SetLineColor(4);
  y7TeV_NSD3->SetFillColor(38);
  y7TeV_NSD3->Draw("e3");
  y7TeV_NSD->SetMarkerColor(4);
  y7TeV_NSD->SetLineColor(4);
  y7TeV_NSD->SetMarkerStyle(20);
  nsdLegend->AddEntry(y7TeV_NSD, " #sqrt{s} = 7 TeV", "pl");
  y7TeV_NSD->Draw("e1,same");
  y7TeV_NSD2->SetMarkerColor(4);
  y7TeV_NSD2->SetLineColor(4);
  y7TeV_NSD2->SetMarkerStyle(20);
  y7TeV_NSD2->Draw("same");
  
  outTextFile << endl;
  outTextFile << "Rapidity result in first bin 7 TeV data -0.2 - 0.2: dXi/dY = " << y7TeV_NSD->GetBinContent(1) << 
  "\\pm" << y7TeV_NSD->GetBinError(1) << "  d(sigma)/dY = " << y7TeV_sigma->GetBinContent(1) << endl;

  TH1F* y900_NSD = y900->Clone("y900_NSD");
  TH1F* y900_NSD2 = y900->Clone("y900_NSD2");
  TH1F* y900_NSD3 = y900->Clone("y900_NSD3");
  TH1F* y900_NSD4 = y900->Clone("y900_NSD4");

  for (int i = 1; i<=y900_NSD->GetNbinsX(); ++i) {
    float stat = y900->GetBinError(i);
    float systUncorOnly = y900->GetBinContent(i)*systUncory;
    float systUncor = sqrt( stat*stat + systUncorOnly*systUncorOnly );
    float systCor = y900->GetBinContent(i)*systCory9;
    float systRatOnly = y900->GetBinContent(i)*systRaty9;
    float systRat = sqrt( systRatOnly*systRatOnly + stat*stat );
    y900_NSD->SetBinContent(i,y900->GetBinContent(i)/nEvt900);  
    y900_NSD->SetBinError(i,stat/nEvt900);  
    y900_NSD2->SetBinContent(i,y900->GetBinContent(i)/nEvt900);  
    y900_NSD2->SetBinError(i,systUncor/nEvt900);  
    y900_NSD3->SetBinContent(i,y900->GetBinContent(i)/nEvt900);  
    y900_NSD3->SetBinError(i,systCor/nEvt900);  
    y900_NSD4->SetBinContent(i,y900->GetBinContent(i)/nEvt900);  
    y900_NSD4->SetBinError(i,systRat/nEvt900);  
  }
  y900_NSD3->SetMarkerColor(2);
  y900_NSD3->SetLineColor(2);
  y900_NSD3->SetMarkerStyle(21);
  y900_NSD3->SetStats(0);
  y900_NSD3->SetFillColor(kRed-9);
  y900_NSD3->Draw("same,e3");
  y900_NSD->SetMarkerColor(2);
  y900_NSD->SetLineColor(2);
  y900_NSD->SetMarkerStyle(21);
  y900_NSD->Draw("e1,same");
  y900_NSD2->SetMarkerColor(2);
  y900_NSD2->SetLineColor(2);
  y900_NSD2->SetMarkerStyle(21);
  y900_NSD2->Draw("same");
  nsdLegend->Draw("f");
  y900_NSD->Draw("same");

  TH1F* y7TeVMCGen_NSD = y7TeVMCGen->Clone("y7TeVMCGen_NSD");
  TH1F* y7TeVMCGenp0_NSD = y7TeVMCGenp0->Clone("y7TeVMCGenp0_NSD");
  TH1F* y7TeVMCGenp8_NSD = y7TeVMCGenp8->Clone("y7TeVMCGenp8_NSD");
  for (int i = 1; i<=y7TeVMCGen_NSD->GetNbinsX(); ++i) {
    y7TeVMCGen_NSD->SetBinContent(i,y7TeVMCGen->GetBinContent(i)/(nEvt7MC));  
    y7TeVMCGen_NSD->SetBinError(i,y7TeVMCGen->GetBinError(i)/(nEvt7MC));  
    y7TeVMCGenp0_NSD->SetBinContent(i,y7TeVMCGenp0->GetBinContent(i)/(nEvt7MCp0));   
    y7TeVMCGenp0_NSD->SetBinError(i,y7TeVMCGenp0->GetBinError(i)/(nEvt7MCp0));  
    y7TeVMCGenp8_NSD->SetBinContent(i,y7TeVMCGenp8->GetBinContent(i)/(nEvt7MCp8));   
    y7TeVMCGenp8_NSD->SetBinError(i,y7TeVMCGenp8->GetBinError(i)/(nEvt7MCp8));  
  }
  y7TeVMCGen_NSD->SetMaximum(1.5*y7TeVMCGen->GetMaximum()/nEvt7MC);
  y7TeVMCGen_NSD->GetYaxis()->SetTitle("(1/N) dN / dy ");
  y7TeVMCGen_NSD->GetXaxis()->SetTitle("#Xi^{-} rapidity");
  y7TeVMCGen_NSD->SetStats(0);
  y7TeVMCGen_NSD->SetLineColor(4);
  nsdLegend->AddEntry(y7TeVMCGen_NSD, "     Pythia6 D6T", "l");
  y7TeVMCGen_NSD->Draw("same,hist");
  y7TeVMCGenp0_NSD->SetLineColor(4);
  y7TeVMCGenp0_NSD->SetLineStyle(3);
//  y7TeVMCGenp0_NSD->Draw("same,hist");
//  nsdLegend->AddEntry(y7TeVMCGenp0_NSD, "     Pythia6 P0", "l");
  y7TeVMCGenp8_NSD->SetLineColor(4);
  y7TeVMCGenp8_NSD->SetLineStyle(2);
//  y7TeVMCGenp8_NSD->Draw("same,hist");
//  nsdLegend->AddEntry(y7TeVMCGenp8_NSD, "     Pythia8", "l");
 
  TH1F* y900MCGen_NSD = y900MCGen->Clone("y900MCGen_NSD");
  TH1F* y900MCGenNoEn_NSD = y900MCGenNoEn->Clone("y900MCGenNoEn_NSD");
  TH1F* y900MCGenD6T_NSD = y900MCGenD6T->Clone("y900MCGenD6T_NSD");
  TH1F* y900MCGenp0_NSD = y900MCGenp0->Clone("y900MCGenp0_NSD");
  TH1F* y900MCGenp8_NSD = y900MCGenp8->Clone("y900MCGenp8_NSD");
  for (int i = 1; i<=y900MCGen_NSD->GetNbinsX(); ++i) {
    y900MCGen_NSD->SetBinContent(i,y900MCGen->GetBinContent(i)/(nEvt900MC));  
    y900MCGenNoEn_NSD->SetBinContent(i,y900MCGenNoEn->GetBinContent(i)/(nEvt900MCnoEn));
    y900MCGenD6T_NSD->SetBinContent(i,y900MCGenD6T->GetBinContent(i)/(nEvt900MCD6T));   
    y900MCGenp0_NSD->SetBinContent(i,y900MCGenp0->GetBinContent(i)/(nEvt900MCp0));   
    y900MCGenp8_NSD->SetBinContent(i,y900MCGenp8->GetBinContent(i)/(nEvt900MCp8));   
  }
  y900MCGenNoEn_NSD->SetMaximum(1.2*y900MCGenNoEn_NSD->GetMaximum());
  y900MCGenNoEn_NSD->GetYaxis()->SetTitle("(1/N) dN / dy");  
  y900MCGenNoEn_NSD->GetXaxis()->SetTitle("#Xi^{-} rapidity");
  y900MCGenNoEn_NSD->SetStats(0);
  y900MCGenNoEn_NSD->SetLineColor(2);
  y900MCGenNoEn_NSD->SetLineStyle(1);
  y900MCGenNoEn_NSD->Draw("same,hist");
  nsdLegend->AddEntry(y900_NSD, "  #sqrt{s} = 900 GeV", "pl");
  nsdLegend->AddEntry(y900MCGenNoEn_NSD, "     Pythia6 D6T", "l");

  y900MCGen_NSD->SetLineColor(2);
  //nsdLegend->AddEntry(y900MCGen_NSD, "#sqrt{s} = 900 GeV MC", "l");
  //y900MCGen_NSD->Draw("same,hist");

  y900MCGenD6T_NSD->SetLineColor(2);
  y900MCGenD6T_NSD->SetLineStyle(1);
//  y900MCGenD6T_NSD->Draw("same,hist");
//  nsdLegend->AddEntry(y900MCGenD6T_NSD, "     Pythia6 D6T", "l");
  y900MCGenp0_NSD->SetLineColor(2);
  y900MCGenp0_NSD->SetLineStyle(3);
//  y900MCGenp0_NSD->Draw("same,hist");
//  nsdLegend->AddEntry(y900MCGenp0_NSD, "     Pythia6 P0", "l");
  y900MCGenp8_NSD->SetLineColor(2);
  y900MCGenp8_NSD->SetLineStyle(2);
//  y900MCGenp8_NSD->Draw("same,hist");
//  nsdLegend->AddEntry(y900MCGenp8_NSD, "     Pythia8", "l");

  outTextFile << "Rapidity result in first bin 900 GeV data -0.2 - 0.2: dXi/dY = " << y900_NSD->GetBinContent(1) << 
  "\\pm" << y900_NSD->GetBinError(1) << "  d(sigma)/dY = " << y900_sigma->GetBinContent(1) << endl;

  outTextFile << "Rapidity result in first bin 7 TeV D6T MC -0.1 - 0.1: dXi/dY = " << y7TeVMCGen_NSD->GetBinContent(1) << 
  "\\pm" << y7TeVMCGen_NSD->GetBinError(1) << 
  endl;
  
  outTextFile << "Rapidity result in first bin 900 GeV MC -0.1 - 0.1: dXi/dY = " << y900MCGenNoEn_NSD->GetBinContent(1) << 
  "\\pm" << y900MCGenNoEn_NSD->GetBinError(1) << 
  endl;
  
  cmsPreUL2->Draw();

  c1->SaveAs(outputDir+"y_NSD_dataMC.eps");

  y7TeV_NSD3->SetMaximum(y7TeV_NSD3->GetMaximum()*1.2);
  y7TeV_NSD3->Draw("e3");
  y7TeV_NSD->Draw("e1,same");
  y7TeV_NSD2->Draw("same");

  y900_NSD3->Draw("e3,same");
  y900_NSD->Draw("e1,same");
  y900_NSD2->Draw("same");
  cmsPre->Draw();
  TLegend* nsd2Legend = new TLegend(0.70, 0.15, 0.95, 0.25);
  nsd2Legend->SetTextSize(0.035);
  nsd2Legend->SetFillColor(0);
  nsd2Legend->AddEntry(y7TeV_NSD, "#sqrt{s} = 7 TeV", "pl");
  nsd2Legend->AddEntry(y900_NSD, "#sqrt{s} = 900 GeV", "pl");
  nsd2Legend->Draw("f");

  c1->SaveAs(outputDir+"y_NSD_data.eps");
  
  TH1F* y7TeVMC_NSD = y7TeVMC->Clone("y7TeVMC_NSD");
  for (int i = 1; i<=y7TeVMC_NSD->GetNbinsX(); ++i) {
    y7TeVMC_NSD->SetBinContent(i,y7TeVMC->GetBinContent(i)/(nEvt7MC));  
    y7TeVMC_NSD->SetBinError(i,y7TeVMC->GetBinError(i)/(nEvt7MC));  
  }
  y7TeVMC_NSD->SetMaximum(2.2*y7TeVMC->GetMaximum()/nEvt7MC);
  y7TeVMC_NSD->GetYaxis()->SetTitle("(1/N) dN / dy ");
  y7TeVMC_NSD->GetXaxis()->SetTitle("#Xi^{-} rapidity");
  y7TeVMC_NSD->SetStats(0);
  y7TeVMC_NSD->SetMarkerColor(4);
  y7TeVMC_NSD->SetLineColor(4);
  y7TeVMC_NSD->SetMarkerStyle(24);
  nsdLegend->AddEntry(y7TeVMC_NSD, "#sqrt{s} = 7 TeV D6T MC Reco", "pl");
  y7TeVMC_NSD->Draw();  

  y7TeVMCGen_NSD->Draw("same, hist");  
  
  TH1F* y900MC_NSD = y900MC->Clone("y900MC_NSD");
  for (int i = 1; i<=y900MC_NSD->GetNbinsX(); ++i) {
    y900MC_NSD->SetBinContent(i,y900MC->GetBinContent(i)/(nEvt900MC));  
    y900MC_NSD->SetBinError(i,y900MC->GetBinError(i)/(nEvt900MC));  
  }
  y900MC_NSD->GetYaxis()->SetTitle("(1/N) dN / dy ");
  y900MC_NSD->GetXaxis()->SetTitle("#Xi^{-} rapidity");
  y900MC_NSD->SetStats(0);
  y900MC_NSD->SetMarkerColor(2);
  y900MC_NSD->SetLineColor(2);
  y900MC_NSD->SetMarkerStyle(25);
  nsdLegend->AddEntry(y900MC_NSD, "#sqrt{s} = 900 GeV MC Reco", "pl");
  y900MC_NSD->Draw("same");
  
  y900MCGen_NSD->Draw("same,hist");
  
  cmsPre->Draw();
  TLegend* nsdMFLegend = new TLegend(0.50, 0.65, 0.75, 0.85);
  nsdMFLegend->SetTextSize(0.035);
  nsdMFLegend->SetFillColor(0);
  nsdMFLegend->AddEntry(y7TeVMCGen_NSD, "#sqrt{s} = 7 TeV D6T MC Gen", "l");
  nsdMFLegend->AddEntry(y7TeVMC_NSD, "#sqrt{s} = 7 TeV D6T MC Reco", "lp");
  nsdMFLegend->AddEntry(y900MCGen_NSD, "#sqrt{s} = 900 GeV MC Gen", "l");
  nsdMFLegend->AddEntry(y900MC_NSD, "#sqrt{s} = 900 GeV MC Reco", "lp");
  nsdMFLegend->Draw("same");
  
  c1->SaveAs(outputDir+"y_NSD_MCFits.eps");
  
  for (int i = 1; i<=y7TeVMCforRatio->GetNbinsX(); ++i){
    y7TeVMCforRatio->SetBinContent(i,y7TeVMCforRatio->GetBinContent(i)/(y7TeVMCforRatio->GetBinWidth(i)*2*nEvt7MC)); // factor of two
    y7TeVMCforRatio->SetBinError(i,y7TeVMCforRatio->GetBinError(i)/(y7TeVMCforRatio->GetBinWidth(i)*2*nEvt7MC));     // from uncor yield
  }
  
  TH1F* y7TeVMCforRatio2 = y7TeVMCforRatio->Clone("y7TeVMCforRatio2");
  y7TeVMCforRatio->Divide(y7TeVMCforRatio, y7TeV_NSD);
  y7TeVMCforRatio->GetXaxis()->SetTitle("#Xi^{-} rapidity");
  y7TeVMCforRatio->GetYaxis()->SetTitle("Ratio MC/data");
  y7TeVMCforRatio->SetMarkerColor(4);
  y7TeVMCforRatio->SetLineColor(4);
  y7TeVMCforRatio->SetMinimum(0);
  y7TeVMCforRatio->SetMaximum(1.5);
  y7TeVMCforRatio->Draw("e1");
  y7TeVMCforRatio2->Divide(y7TeVMCforRatio2, y7TeV_NSD4);
  y7TeVMCforRatio2->SetMarkerColor(4);
  y7TeVMCforRatio2->SetLineColor(4);
  y7TeVMCforRatio2->Draw("same");
  
  for (int i = 1; i<=y900MCforRatio->GetNbinsX(); ++i){
    y900MCforRatio->SetBinContent(i,y900MCforRatio->GetBinContent(i)/(y900MCforRatio->GetBinWidth(i)*2*nEvt900MC)); // factor of two
    y900MCforRatio->SetBinError(i,y900MCforRatio->GetBinError(i)/(y900MCforRatio->GetBinWidth(i)*2*nEvt900MC));     // from uncor yield
  }
  
  TH1F* y900MCforRatio2 = y900MCforRatio->Clone("y900MCforRatio2");
  y900MCforRatio->Divide(y900MCforRatio, y900_NSD);
  y900MCforRatio->SetMarkerColor(2);
  y900MCforRatio->SetLineColor(2);
  y900MCforRatio->Draw("e1,same");
  y900MCforRatio2->Divide(y900MCforRatio2, y900_NSD4);
  y900MCforRatio2->SetMarkerColor(2);
  y900MCforRatio2->SetLineColor(2);
  y900MCforRatio2->Draw("same");
  
  cmsPreUL->Draw();
  nsdRatLegend->Draw("same");
  
  c1->SaveAs(outputDir+"y_NSD_ratio.eps");
  
  for (int i = 1; i<=y900MCforFitVal->GetNbinsX(); ++i){
    y900MCforFitVal->SetBinContent(i,y900MCforFitVal->GetBinContent(i)/(y900MCforFitVal->GetBinWidth(i)*2));  // factor of two
    y900MCforFitVal->SetBinError(i,y900MCforFitVal->GetBinError(i)/(y900MCforFitVal->GetBinWidth(i)*2));      // from uncor yield
  }
  
  y900MCforFitVal->Divide(y900MC, y900MCforFitVal);
  y900MCforFitVal->GetXaxis()->SetTitle("#Xi^{-} rapidity");
  y900MCforFitVal->GetYaxis()->SetTitle("Ratio MC fit/MC truth");
  y900MCforFitVal->SetLineColor(2);
  y900MCforFitVal->SetMarkerColor(2);
  //y900MCforFitVal->SetMinimum(0);
  //y900MCforFitVal->SetMaximum(1.5);
  y900MCforFitVal->Draw();

  for (int i = 1; i<=y7TeVMCforFitVal->GetNbinsX(); ++i){
    y7TeVMCforFitVal->SetBinContent(i,y7TeVMCforFitVal->GetBinContent(i)/(y7TeVMCforFitVal->GetBinWidth(i)*2));  // factor of two
    y7TeVMCforFitVal->SetBinError(i,y7TeVMCforFitVal->GetBinError(i)/(y7TeVMCforFitVal->GetBinWidth(i)*2));      // from uncor yield
  }
  
  y7TeVMCforFitVal->Divide(y7TeVMC, y7TeVMCforFitVal);
  y7TeVMCforFitVal->SetLineColor(4);
  y7TeVMCforFitVal->SetMarkerColor(4);
  y7TeVMCforFitVal->Draw("same");
  
  cmsPreUL->Draw();
  TLegend* nsdRatLegend = new TLegend(0.35, 0.75, 0.6, 0.85);
  nsdRatLegend->SetTextSize(0.035);
  nsdRatLegend->SetFillColor(0);
  nsdRatLegend->AddEntry(y7TeVMCforFitVal, "#sqrt{s} = 7 TeV", "l");
  nsdRatLegend->AddEntry(y900MCforFitVal, "#sqrt{s} = 900 GeV", "l");
  nsdRatLegend->Draw("same");
  
  c1->SaveAs(outputDir+"y_fitval_ratio.eps");
  
  outTextFile << endl;

  exp1->SetParLimits(1,-0.3,0.0);
  //exp1->SetParameters(65000,-0.16);
  //exp1->SetParLimits(0,60000,70000);


  exp1->SetParameters(ct7TeV->GetMaximum(), -0.16);
  ct7TeV->Fit("exp1","IEM");
  ct7TeV->Draw();

  outTextFile << "integral ctau from 0.0 to 10.0 7TeV = " << exp1->Integral(0.0,10.0) << "  from 0.0 to 10000.0 = " << exp1->Integral(0.0,10000.0) <<
  " value calculated from fit results = " << -exp1->GetParameter(0)/exp1->GetParameter(1) << endl;
  
  exp1->SetParameters(ct900->GetMaximum(), -0.16);
  ct900->Fit("exp1","IEM");
  ct900->Draw();

  outTextFile << "integral ctau from 0.0 to 10.0 900 GeV = " << exp1->Integral(0.0,10.0) << "  from 0.0 to 10000.0 = " << exp1->Integral(0.0,10000.0) <<
  " value calculated from fit results = " << -exp1->GetParameter(0)/exp1->GetParameter(1) << endl;

  for (int i = 1; i<=ct900MCforFitVal->GetNbinsX(); ++i){
    ct900MCforFitVal->SetBinContent(i,ct900MCforFitVal->GetBinContent(i)/(ct900MCforFitVal->GetBinWidth(i)));
    ct900MCforFitVal->SetBinError(i,ct900MCforFitVal->GetBinError(i)/(ct900MCforFitVal->GetBinWidth(i)));
  }
  
  ct900MCforFitVal->Divide(ct900MC, ct900MCforFitVal);
  ct900MCforFitVal->GetXaxis()->SetTitle("#Xi^{-} ct (cm)");
  ct900MCforFitVal->GetYaxis()->SetTitle("Ratio MC fit/MC truth");
  ct900MCforFitVal->SetLineColor(2);
  ct900MCforFitVal->SetMarkerColor(2);
  //ct900MCforFitVal->SetMinimum(0);
  //ct900MCforFitVal->SetMaximum(1.5);
  ct900MCforFitVal->Draw();

  for (int i = 1; i<=ct7TeVMCforFitVal->GetNbinsX(); ++i){
    ct7TeVMCforFitVal->SetBinContent(i,ct7TeVMCforFitVal->GetBinContent(i)/ct7TeVMCforFitVal->GetBinWidth(i));
    ct7TeVMCforFitVal->SetBinError(i,ct7TeVMCforFitVal->GetBinError(i)/(ct7TeVMCforFitVal->GetBinWidth(i)));
  }
  
  ct7TeVMCforFitVal->Divide(ct7TeVMC, ct7TeVMCforFitVal);
  ct7TeVMCforFitVal->SetLineColor(4);
  ct7TeVMCforFitVal->SetMarkerColor(4);
  ct7TeVMCforFitVal->Draw("same");
  
  cmsPreUL->Draw();
  TLegend* nsdRatLegend = new TLegend(0.35, 0.75, 0.6, 0.85);
  nsdRatLegend->SetTextSize(0.035);
  nsdRatLegend->SetFillColor(0);
  nsdRatLegend->AddEntry(ct7TeVMCforFitVal, "#sqrt{s} = 7 TeV", "l");
  nsdRatLegend->AddEntry(ct900MCforFitVal, "#sqrt{s} = 900 GeV", "l");
  nsdRatLegend->Draw("same");
  
  c1->SaveAs(outputDir+"ct_fitval_ratio.eps");










  // Get Asymmetry plot here
  TCanvas* theAsymPtCanvas = new TCanvas("theAsymPtCanvas","theAsymPtCanvas",600,600);
  TLegend* asymPtLegend = new TLegend(0.2, 0.2, 0.85, 0.4);
  asymPtLegend->SetFillColor(0);
  asymPtLegend->SetTextSize(0.042);
  asymPtLegend->SetBorderSize(0);
  asymPtLegend->SetFillStyle(0);
  theAsymPtCanvas->cd();

  TGraphAsymmErrors* theAsymmPt7TeVGraph  = new TGraphAsymmErrors();
  TGraphAsymmErrors* theCorrAsymmPt7TeVGraph  = new TGraphAsymmErrors();
  TGraphAsymmErrors* theCorrAsymmPt7TeVGraph2  = new TGraphAsymmErrors();
  TGraphAsymmErrors* theCorrAsymmPt7TeVGraph3  = new TGraphAsymmErrors();
  TGraphAsymmErrors* theAsymmPt900GeVGraph  = new TGraphAsymmErrors();
  TGraphAsymmErrors* theCorrAsymmPt900GeVGraph  = new TGraphAsymmErrors();
  TGraphAsymmErrors* theCorrAsymmPt900GeVGraph2  = new TGraphAsymmErrors();
  TGraphAsymmErrors* theCorrAsymmPt900GeVGraph3  = new TGraphAsymmErrors();
  TGraphAsymmErrors* theAsymmY7TeVGraph  = new TGraphAsymmErrors();
  TGraphAsymmErrors* theAsymmY900GeVGraph  = new TGraphAsymmErrors();
  TGraphAsymmErrors* theCorrAsymmY7TeVGraph  = new TGraphAsymmErrors();
  TGraphAsymmErrors* theCorrAsymmY7TeVGraph2  = new TGraphAsymmErrors();
  TGraphAsymmErrors* theCorrAsymmY7TeVGraph3  = new TGraphAsymmErrors();
  TGraphAsymmErrors* theCorrAsymmY900GeVGraph  = new TGraphAsymmErrors();
  TGraphAsymmErrors* theCorrAsymmY900GeVGraph2  = new TGraphAsymmErrors();
  TGraphAsymmErrors* theCorrAsymmY900GeVGraph3  = new TGraphAsymmErrors();

  
  TAxis* theXAxis = thePosYieldsPt7TeVHisto->GetXaxis();
  unsigned int theBins = theXAxis->GetNbins();
  float lowLimit = 0.;
  float hiLimit = 0.;
  
  for(unsigned int vIt = 0; vIt < theBins; ++vIt){
    // first get values for normal binning
    if(vIt == 0) lowLimit = theXAxis->GetBinLowEdge(vIt + 1);
    double asymmetryMinus = theNegCorrYieldsPt7TeVHisto->GetBinContent(vIt + 1) - thePosCorrYieldsPt7TeVHisto->GetBinContent(vIt + 1);							  
    double asymmetryPlus = thePosCorrYieldsPt7TeVHisto->GetBinContent(vIt + 1) + theNegCorrYieldsPt7TeVHisto->GetBinContent(vIt + 1);							  
    double asymmetry = asymmetryMinus / asymmetryPlus;													  
    theCorrAsymmPt7TeVGraph->SetPoint(vIt, theXAxis->GetBinCenter(vIt + 1), asymmetry);							  
    theCorrAsymmPt7TeVGraph->SetPointEXlow(vIt,theXAxis->GetBinCenter(vIt + 1)-theXAxis->GetBinLowEdge(vIt + 1));					  
    theCorrAsymmPt7TeVGraph->SetPointEXhigh(vIt,theXAxis->GetBinLowEdge(vIt + 1)-theXAxis->GetBinCenter(vIt + 1)+theXAxis->GetBinWidth(vIt + 1)); 
    double plusminusErrorSqrd = thePosCorrYieldsPt7TeVHisto->GetBinError(vIt + 1) * thePosCorrYieldsPt7TeVHisto->GetBinError(vIt + 1) + 						  
    			       theNegCorrYieldsPt7TeVHisto->GetBinError(vIt + 1) * theNegCorrYieldsPt7TeVHisto->GetBinError(vIt + 1);							  
    double asymmError = (sqrt(plusminusErrorSqrd)/asymmetryPlus)*sqrt(1. + asymmetry*asymmetry);							  
    theCorrAsymmPt7TeVGraph->SetPointEYhigh(vIt,asymmError);													  
    theCorrAsymmPt7TeVGraph->SetPointEYlow(vIt,asymmError);
    hiLimit = theXAxis->GetBinUpEdge(vIt + 1);											  
    
    double plusminusErrorSystSqrd = thePosCorrYieldsPt7TeVHisto->GetBinContent(vIt + 1)*systRatpt*thePosCorrYieldsPt7TeVHisto->GetBinContent(vIt + 1)*systRatpt +
                                theNegCorrYieldsPt7TeVHisto->GetBinContent(vIt + 1)*systRatpt*theNegCorrYieldsPt7TeVHisto->GetBinContent(vIt + 1)*systRatpt;
    double asymmErrorSyst = (sqrt(plusminusErrorSystSqrd)/asymmetryPlus)*sqrt(1. + asymmetry*asymmetry);
    theCorrAsymmPt7TeVGraph2->SetPoint(vIt, theXAxis->GetBinCenter(vIt + 1), asymmetry);
    theCorrAsymmPt7TeVGraph2->SetPointEYhigh(vIt,asymmErrorSyst);
    theCorrAsymmPt7TeVGraph2->SetPointEYlow(vIt,asymmErrorSyst); 
    
    double plusminusErrorSystUncorSqrd = thePosCorrYieldsPt7TeVHisto->GetBinContent(vIt + 1)*systUncorpt*thePosCorrYieldsPt7TeVHisto->GetBinContent(vIt + 1)*systUncorpt +
                                    theNegCorrYieldsPt7TeVHisto->GetBinContent(vIt + 1)*systUncorpt*theNegCorrYieldsPt7TeVHisto->GetBinContent(vIt + 1)*systUncorpt;
    double asymmErrorSystUncor = (sqrt(plusminusErrorSystUncorSqrd)/asymmetryPlus)*sqrt(1. + asymmetry*asymmetry);
    theCorrAsymmPt7TeVGraph3->SetPoint(vIt, theXAxis->GetBinCenter(vIt + 1), asymmetry);
    theCorrAsymmPt7TeVGraph3->SetPointEYhigh(vIt,asymmErrorSystUncor);
    theCorrAsymmPt7TeVGraph3->SetPointEYlow(vIt,asymmErrorSystUncor); 
  }
  theCorrAsymmPt7TeVGraph->GetXaxis()->SetLimits(lowLimit, hiLimit); 
  theCorrAsymmPt7TeVGraph->SetMaximum(1); 					     
  theCorrAsymmPt7TeVGraph->SetMinimum(-1);					     
  theCorrAsymmPt7TeVGraph->SetMarkerColor(4);
  theCorrAsymmPt7TeVGraph->SetMarkerStyle(20);
  theCorrAsymmPt7TeVGraph->SetLineColor(4);			     
  theCorrAsymmPt7TeVGraph->SetMarkerSize(0.8);					     
  theCorrAsymmPt7TeVGraph->GetYaxis()->SetTitle("(#Xi^{-} - #Xi^{+})/(#Xi^{-} + #Xi^{+})");
  theCorrAsymmPt7TeVGraph->GetXaxis()->SetTitle("#Xi p_{T} (GeV/c)");
  theCorrAsymmPt7TeVGraph->Draw("ap");	
  asymPtLegend->AddEntry(theCorrAsymmPt7TeVGraph, "#sqrt{s} = 7 TeV Corrected", "pl");
  
  for(unsigned int vIt = 0; vIt < theBins; ++vIt){
    if(vIt == 0) lowLimit = theXAxis->GetBinLowEdge(vIt + 1);
    double asymmetryMinus = theNegYieldsPt7TeVHisto->GetBinContent(vIt + 1) - thePosYieldsPt7TeVHisto->GetBinContent(vIt + 1);							  
    double asymmetryPlus = thePosYieldsPt7TeVHisto->GetBinContent(vIt + 1) + theNegYieldsPt7TeVHisto->GetBinContent(vIt + 1);							  
    double asymmetry = asymmetryMinus / asymmetryPlus;													  
    theAsymmPt7TeVGraph->SetPoint(vIt, theXAxis->GetBinCenter(vIt + 1), asymmetry);							  
    theAsymmPt7TeVGraph->SetPointEXlow(vIt,theXAxis->GetBinCenter(vIt + 1)-theXAxis->GetBinLowEdge(vIt + 1));					  
    theAsymmPt7TeVGraph->SetPointEXhigh(vIt,theXAxis->GetBinLowEdge(vIt + 1)-theXAxis->GetBinCenter(vIt + 1)+theXAxis->GetBinWidth(vIt + 1)); 
    double plusminusErrorSqrd = thePosYieldsPt7TeVHisto->GetBinError(vIt + 1) * thePosYieldsPt7TeVHisto->GetBinError(vIt + 1) + 						  
    			       theNegYieldsPt7TeVHisto->GetBinError(vIt + 1) * theNegYieldsPt7TeVHisto->GetBinError(vIt + 1);							  
    double asymmError = (sqrt(plusminusErrorSqrd)/asymmetryPlus)*sqrt(1. + asymmetry*asymmetry);							  
    theAsymmPt7TeVGraph->SetPointEYhigh(vIt,asymmError);													  
    theAsymmPt7TeVGraph->SetPointEYlow(vIt,asymmError);
    hiLimit = theXAxis->GetBinUpEdge(vIt + 1);											  
  }
  theAsymmPt7TeVGraph->GetXaxis()->SetLimits(lowLimit, hiLimit); 
  theAsymmPt7TeVGraph->SetMaximum(1); 					     
  theAsymmPt7TeVGraph->SetMinimum(-1);	
  theAsymmPt7TeVGraph->SetMarkerStyle(23);			     
  theAsymmPt7TeVGraph->SetMarkerColor(7);
  theAsymmPt7TeVGraph->SetLineColor(7);			     
  theAsymmPt7TeVGraph->SetMarkerSize(0.8);					     
  theAsymmPt7TeVGraph->GetYaxis()->SetTitle("(#Xi^{-} - #Xi^{+})/(#Xi^{-} + #Xi^{+})");
  theAsymmPt7TeVGraph->GetXaxis()->SetTitle("#Xi^{-} p_{T} (GeV/c)");
  theAsymmPt7TeVGraph->Draw("p same");	
  asymPtLegend->AddEntry(theAsymmPt7TeVGraph, "#sqrt{s} = 7 TeV Uncorrected", "pl");
  
  TAxis* theXAxis = thePosYieldsPt900GeVHisto->GetXaxis();
  unsigned int theBins = theXAxis->GetNbins();
  float lowLimit = 0.;
  float hiLimit = 0.;
  for(unsigned int vIt = 0; vIt < theBins; ++vIt){
    if(vIt == 0) lowLimit = theXAxis->GetBinLowEdge(vIt + 1);
    double asymmetryMinus = theNegCorrYieldsPt900GeVHisto->GetBinContent(vIt + 1) - thePosCorrYieldsPt900GeVHisto->GetBinContent(vIt + 1);							  
    double asymmetryPlus = thePosCorrYieldsPt900GeVHisto->GetBinContent(vIt + 1) + theNegCorrYieldsPt900GeVHisto->GetBinContent(vIt + 1);							  
    double asymmetry = asymmetryMinus / asymmetryPlus;													  
    theCorrAsymmPt900GeVGraph->SetPoint(vIt, theXAxis->GetBinCenter(vIt + 1), asymmetry);							  
    theCorrAsymmPt900GeVGraph->SetPointEXlow(vIt,theXAxis->GetBinCenter(vIt + 1)-theXAxis->GetBinLowEdge(vIt + 1));					  
    theCorrAsymmPt900GeVGraph->SetPointEXhigh(vIt,theXAxis->GetBinLowEdge(vIt + 1)-theXAxis->GetBinCenter(vIt + 1)+theXAxis->GetBinWidth(vIt + 1)); 
    double plusminusErrorSqrd = thePosCorrYieldsPt900GeVHisto->GetBinError(vIt + 1) * thePosCorrYieldsPt900GeVHisto->GetBinError(vIt + 1) + 						  
    			       theNegCorrYieldsPt900GeVHisto->GetBinError(vIt + 1) * theNegCorrYieldsPt900GeVHisto->GetBinError(vIt + 1);							  
    double asymmError = (sqrt(plusminusErrorSqrd)/asymmetryPlus)*sqrt(1. + asymmetry*asymmetry);							  
    theCorrAsymmPt900GeVGraph->SetPointEYhigh(vIt,asymmError);													  
    theCorrAsymmPt900GeVGraph->SetPointEYlow(vIt,asymmError);
    hiLimit = theXAxis->GetBinUpEdge(vIt + 1);	
    
    double plusminusErrorSystSqrd = thePosCorrYieldsPt900GeVHisto->GetBinContent(vIt + 1)*systRatpt*thePosCorrYieldsPt900GeVHisto->GetBinContent(vIt + 1)*systRatpt +
                                theNegCorrYieldsPt900GeVHisto->GetBinContent(vIt + 1)*systRatpt*theNegCorrYieldsPt900GeVHisto->GetBinContent(vIt + 1)*systRatpt;
    double asymmErrorSyst = (sqrt(plusminusErrorSystSqrd)/asymmetryPlus)*sqrt(1. + asymmetry*asymmetry);
    theCorrAsymmPt900GeVGraph2->SetPoint(vIt, theXAxis->GetBinCenter(vIt + 1), asymmetry);
    theCorrAsymmPt900GeVGraph2->SetPointEYhigh(vIt,asymmErrorSyst);
    theCorrAsymmPt900GeVGraph2->SetPointEYlow(vIt,asymmErrorSyst); 
    
    double plusminusErrorSystUncorSqrd = thePosCorrYieldsPt900GeVHisto->GetBinContent(vIt + 1)*systUncorpt*thePosCorrYieldsPt900GeVHisto->GetBinContent(vIt + 1)*systUncorpt +
                                    theNegCorrYieldsPt900GeVHisto->GetBinContent(vIt + 1)*systUncorpt*theNegCorrYieldsPt900GeVHisto->GetBinContent(vIt + 1)*systUncorpt;
    double asymmErrorSystUncor = (sqrt(plusminusErrorSystUncorSqrd)/asymmetryPlus)*sqrt(1. + asymmetry*asymmetry);
    theCorrAsymmPt900GeVGraph3->SetPoint(vIt, theXAxis->GetBinCenter(vIt + 1), asymmetry);
    theCorrAsymmPt900GeVGraph3->SetPointEYhigh(vIt,asymmErrorSystUncor);
    theCorrAsymmPt900GeVGraph3->SetPointEYlow(vIt,asymmErrorSystUncor); 
    										  
  }
  theCorrAsymmPt900GeVGraph->GetXaxis()->SetLimits(lowLimit, hiLimit); 
  theCorrAsymmPt900GeVGraph->SetMaximum(1); 					     
  theCorrAsymmPt900GeVGraph->SetMinimum(-1);					     
  theCorrAsymmPt900GeVGraph->SetMarkerStyle(21);
  theCorrAsymmPt900GeVGraph->SetMarkerColor(2);
  theCorrAsymmPt900GeVGraph->SetLineColor(2);			     
  theCorrAsymmPt900GeVGraph->SetMarkerSize(0.8);					     
  theCorrAsymmPt900GeVGraph->GetYaxis()->SetTitle("(#Xi^{-} - #Xi^{+})/(#Xi^{-} + #Xi^{+})");
  theCorrAsymmPt900GeVGraph->GetXaxis()->SetTitle("#Xi^{-} p_{T} (GeV/c)");
  theCorrAsymmPt900GeVGraph->Draw("p same");	
  asymPtLegend->AddEntry(theCorrAsymmPt900GeVGraph, "#sqrt{s} = 900GeV Corrected", "pl");
  
  for(unsigned int vIt = 0; vIt < theBins; ++vIt){
    if(vIt == 0) lowLimit = theXAxis->GetBinLowEdge(vIt + 1);
    double asymmetryMinus = theNegYieldsPt900GeVHisto->GetBinContent(vIt + 1) - thePosYieldsPt900GeVHisto->GetBinContent(vIt + 1);							  
    double asymmetryPlus = thePosYieldsPt900GeVHisto->GetBinContent(vIt + 1) + theNegYieldsPt900GeVHisto->GetBinContent(vIt + 1);							  
    double asymmetry = asymmetryMinus / asymmetryPlus;													  
    theAsymmPt900GeVGraph->SetPoint(vIt, theXAxis->GetBinCenter(vIt + 1), asymmetry);							  
    theAsymmPt900GeVGraph->SetPointEXlow(vIt,theXAxis->GetBinCenter(vIt + 1)-theXAxis->GetBinLowEdge(vIt + 1));					  
    theAsymmPt900GeVGraph->SetPointEXhigh(vIt,theXAxis->GetBinLowEdge(vIt + 1)-theXAxis->GetBinCenter(vIt + 1)+theXAxis->GetBinWidth(vIt + 1)); 
    double plusminusErrorSqrd = thePosYieldsPt900GeVHisto->GetBinError(vIt + 1) * thePosYieldsPt900GeVHisto->GetBinError(vIt + 1) + 						  
    			       theNegYieldsPt900GeVHisto->GetBinError(vIt + 1) * theNegYieldsPt900GeVHisto->GetBinError(vIt + 1);							  
    double asymmError = (sqrt(plusminusErrorSqrd)/asymmetryPlus)*sqrt(1. + asymmetry*asymmetry);							  
    theAsymmPt900GeVGraph->SetPointEYhigh(vIt,asymmError);													  
    theAsymmPt900GeVGraph->SetPointEYlow(vIt,asymmError);
    hiLimit = theXAxis->GetBinUpEdge(vIt + 1);											  
  }
  theAsymmPt900GeVGraph->SetMaximum(1); 					     
  theAsymmPt900GeVGraph->SetMinimum(-1);					     
  theAsymmPt900GeVGraph->SetMarkerStyle(22);					     
  theAsymmPt900GeVGraph->SetMarkerColor(6);					     
  theAsymmPt900GeVGraph->SetLineColor(6);					     
  theAsymmPt900GeVGraph->SetMarkerSize(0.8);					     
  theAsymmPt900GeVGraph->GetYaxis()->SetTitle("(#Xi^{-} - #Xi^{+})/(#Xi^{-} + #Xi^{+})");
  theAsymmPt900GeVGraph->GetXaxis()->SetTitle("#Xi^{-} p_{T} (GeV/c)");
  asymPtLegend->AddEntry(theAsymmPt900GeVGraph, "#sqrt{s} = 900 GeV Uncorrected", "pl");
  asymPtLegend->Draw("F");
  theAsymmPt900GeVGraph->Draw("p same");	
  
  
  TLine* zeroLine = new TLine(1,0,6,0);
  zeroLine->Draw("same");
  
  cmsPre->Draw();				     
  theAsymPtCanvas->Print(outputDir+"xiAsymPtUncor.eps");
  
  TLegend* asymPtLegend2 = new TLegend(0.2, 0.2, 0.85, 0.4);  
  asymPtLegend2->SetFillColor(0);
  asymPtLegend2->SetTextSize(0.042);
  asymPtLegend2->SetBorderSize(0);
  asymPtLegend2->SetFillStyle(0);
  asymPtLegend2->AddEntry(theCorrAsymmPt900GeVGraph, "#sqrt{s} = 900 GeV Corrected", "pl");
  asymPtLegend2->AddEntry(theCorrAsymmPt7TeVGraph, "#sqrt{s} = 7 TeV Corrected", "pl");
  theCorrAsymmPt7TeVGraph2->SetLineColor(4);
  theCorrAsymmPt7TeVGraph2->SetFillColor(38);
  theCorrAsymmPt7TeVGraph2->SetMaximum(1); 					     
  theCorrAsymmPt7TeVGraph2->SetMinimum(-1);					     
  theCorrAsymmPt7TeVGraph2->SetMarkerColor(4);
  theCorrAsymmPt7TeVGraph2->SetMarkerStyle(20);
  theCorrAsymmPt7TeVGraph2->SetLineColor(4);			     
  theCorrAsymmPt7TeVGraph2->SetMarkerSize(0.8);					     
  theCorrAsymmPt7TeVGraph2->GetYaxis()->SetTitle("(#Xi^{-} - #Xi^{+})/(#Xi^{-} + #Xi^{+})");
  theCorrAsymmPt7TeVGraph2->GetXaxis()->SetTitle("#Xi^{-} p_{T} (GeV/c)");
  theCorrAsymmPt7TeVGraph2->Draw("e3,ap");
  theCorrAsymmPt7TeVGraph->Draw("e1,p");
  theCorrAsymmPt7TeVGraph3->SetLineColor(4);
  theCorrAsymmPt7TeVGraph3->SetMarkerColor(4);
  theCorrAsymmPt7TeVGraph3->Draw("p");
  
  theCorrAsymmPt900GeVGraph2->SetLineColor(2);
  theCorrAsymmPt900GeVGraph2->SetFillColor(kRed-9);				     
  theCorrAsymmPt900GeVGraph2->SetMarkerColor(2);
  theCorrAsymmPt900GeVGraph2->SetMarkerStyle(21);
  theCorrAsymmPt900GeVGraph2->SetLineColor(2);			     
  theCorrAsymmPt900GeVGraph2->SetMarkerSize(0.8);			
  theCorrAsymmPt900GeVGraph2->Draw("e3,p");
  theCorrAsymmPt900GeVGraph->Draw("e1,p");
  theCorrAsymmPt900GeVGraph3->SetLineColor(4);
  theCorrAsymmPt900GeVGraph3->SetMarkerColor(4);
  theCorrAsymmPt900GeVGraph3->Draw("p");
  	
  theCorrAsymmPt900GeVGraph->Draw("p same");	
  asymPtLegend2->Draw("F");
  cmsPre->Draw();	
  theAsymPtCanvas->Print(outputDir+"xiAsymPt.eps");

  TCanvas* theAsymYCanvas = new TCanvas("theAsymYCanvas","theAsymYCanvas",600,600);  
  theAsymYCanvas->cd();
  TLegend* asymYLegend = new TLegend(0.2, 0.2, 0.85, 0.4);  
  asymYLegend->SetFillColor(0);
  asymYLegend->SetTextSize(0.042);
  asymYLegend->SetBorderSize(0);
  asymYLegend->SetFillStyle(0);
  TAxis* theXAxis = thePosYieldsY7TeVHisto->GetXaxis();
  unsigned int theBins = theXAxis->GetNbins();
  float lowLimit = 0.;
  float hiLimit = 0.;
  
  for(unsigned int vIt = 0; vIt < theBins; ++vIt){
    if(vIt == 0) lowLimit = theXAxis->GetBinLowEdge(vIt + 1);
    double asymmetryMinus = theNegCorrYieldsY7TeVHisto->GetBinContent(vIt + 1) - thePosCorrYieldsY7TeVHisto->GetBinContent(vIt + 1);							  
    double asymmetryPlus = thePosCorrYieldsY7TeVHisto->GetBinContent(vIt + 1) + theNegCorrYieldsY7TeVHisto->GetBinContent(vIt + 1);							  
    double asymmetry = asymmetryMinus / asymmetryPlus;													  
    theCorrAsymmY7TeVGraph->SetPoint(vIt, theXAxis->GetBinCenter(vIt + 1), asymmetry);							  
    theCorrAsymmY7TeVGraph->SetPointEXlow(vIt,theXAxis->GetBinCenter(vIt + 1)-theXAxis->GetBinLowEdge(vIt + 1));					  
    theCorrAsymmY7TeVGraph->SetPointEXhigh(vIt,theXAxis->GetBinLowEdge(vIt + 1)-theXAxis->GetBinCenter(vIt + 1)+theXAxis->GetBinWidth(vIt + 1)); 
    double plusminusErrorSqrd = thePosCorrYieldsY7TeVHisto->GetBinError(vIt + 1) * thePosCorrYieldsY7TeVHisto->GetBinError(vIt + 1) + 						  
    			       theNegCorrYieldsY7TeVHisto->GetBinError(vIt + 1) * theNegCorrYieldsY7TeVHisto->GetBinError(vIt + 1);							  
    double asymmError = (sqrt(plusminusErrorSqrd)/asymmetryPlus)*sqrt(1. + asymmetry*asymmetry);							  
    theCorrAsymmY7TeVGraph->SetPointEYhigh(vIt,asymmError);													  
    theCorrAsymmY7TeVGraph->SetPointEYlow(vIt,asymmError);
    hiLimit = theXAxis->GetBinUpEdge(vIt + 1);											  

    double plusminusErrorSystSqrd = thePosCorrYieldsY7TeVHisto->GetBinContent(vIt + 1)*systRaty7*thePosCorrYieldsY7TeVHisto->GetBinContent(vIt + 1)*systRaty7 +
                                theNegCorrYieldsY7TeVHisto->GetBinContent(vIt +	1)*systRaty7*theNegCorrYieldsY7TeVHisto->GetBinContent(vIt + 1)*systRaty7;
    double asymmErrorSyst = (sqrt(plusminusErrorSystSqrd)/asymmetryPlus)*sqrt(1. + asymmetry*asymmetry);
    theCorrAsymmY7TeVGraph2->SetPoint(vIt, theXAxis->GetBinCenter(vIt + 1), asymmetry);
    theCorrAsymmY7TeVGraph2->SetPointEYhigh(vIt,asymmErrorSyst);
    theCorrAsymmY7TeVGraph2->SetPointEYlow(vIt,asymmErrorSyst); 
    
    double plusminusErrorSystUncorSqrd = thePosCorrYieldsY7TeVHisto->GetBinContent(vIt + 1)*systUncory*thePosCorrYieldsY7TeVHisto->GetBinContent(vIt + 1)*systUncory +
                                    theNegCorrYieldsY7TeVHisto->GetBinContent(vIt + 1)*systUncory*theNegCorrYieldsY7TeVHisto->GetBinContent(vIt + 1)*systUncory;
    double asymmErrorSystUncor = (sqrt(plusminusErrorSystUncorSqrd)/asymmetryPlus)*sqrt(1. + asymmetry*asymmetry);
    theCorrAsymmY7TeVGraph3->SetPoint(vIt, theXAxis->GetBinCenter(vIt + 1), asymmetry);
    theCorrAsymmY7TeVGraph3->SetPointEYhigh(vIt,asymmErrorSystUncor);
    theCorrAsymmY7TeVGraph3->SetPointEYlow(vIt,asymmErrorSystUncor); 

  }
  theCorrAsymmY7TeVGraph->GetXaxis()->SetLimits(lowLimit, hiLimit); 
  theCorrAsymmY7TeVGraph->SetMaximum(1); 					     
  theCorrAsymmY7TeVGraph->SetMinimum(-1);					     
  theCorrAsymmY7TeVGraph->SetMarkerColor(4);
  theCorrAsymmY7TeVGraph->SetMarkerStyle(20);
  theCorrAsymmY7TeVGraph->SetLineColor(4);			     
  theCorrAsymmY7TeVGraph->SetMarkerSize(0.8);					     
  theCorrAsymmY7TeVGraph->GetYaxis()->SetTitle("(#Xi^{-} - #Xi^{+})/(#Xi^{-} + #Xi^{+})");
  theCorrAsymmY7TeVGraph->GetXaxis()->SetTitle("#Xi^{-} rapidity");
  theCorrAsymmY7TeVGraph->Draw("ap");	
  asymYLegend->AddEntry(theCorrAsymmY7TeVGraph, "#sqrt{s} = 7 TeV Corrected", "pl");

  for(unsigned int vIt = 0; vIt < theBins; ++vIt){
    if(vIt == 0) lowLimit = theXAxis->GetBinLowEdge(vIt + 1);
    double asymmetryMinus = theNegYieldsY7TeVHisto->GetBinContent(vIt + 1) - thePosYieldsY7TeVHisto->GetBinContent(vIt + 1);							  
    double asymmetryPlus = thePosYieldsY7TeVHisto->GetBinContent(vIt + 1) + theNegYieldsY7TeVHisto->GetBinContent(vIt + 1);							  
    double asymmetry = asymmetryMinus / asymmetryPlus;													  
    theAsymmY7TeVGraph->SetPoint(vIt, theXAxis->GetBinCenter(vIt + 1), asymmetry);							  
    theAsymmY7TeVGraph->SetPointEXlow(vIt,theXAxis->GetBinCenter(vIt + 1)-theXAxis->GetBinLowEdge(vIt + 1));					  
    theAsymmY7TeVGraph->SetPointEXhigh(vIt,theXAxis->GetBinLowEdge(vIt + 1)-theXAxis->GetBinCenter(vIt + 1)+theXAxis->GetBinWidth(vIt + 1)); 
    double plusminusErrorSqrd = thePosYieldsY7TeVHisto->GetBinError(vIt + 1) * thePosYieldsY7TeVHisto->GetBinError(vIt + 1) + 						  
    			       theNegYieldsY7TeVHisto->GetBinError(vIt + 1) * theNegYieldsY7TeVHisto->GetBinError(vIt + 1);							  
    double asymmError = (sqrt(plusminusErrorSqrd)/asymmetryPlus)*sqrt(1. + asymmetry*asymmetry);							  
    theAsymmY7TeVGraph->SetPointEYhigh(vIt,asymmError);													  
    theAsymmY7TeVGraph->SetPointEYlow(vIt,asymmError);
    hiLimit = theXAxis->GetBinUpEdge(vIt + 1);											  
  }
  theAsymmY7TeVGraph->GetXaxis()->SetLimits(lowLimit, hiLimit); 
  theAsymmY7TeVGraph->SetMaximum(1); 					     
  theAsymmY7TeVGraph->SetMinimum(-1);					     
  theAsymmY7TeVGraph->SetMarkerColor(7);
  theAsymmY7TeVGraph->SetLineColor(7);			     
  theAsymmY7TeVGraph->SetMarkerSize(0.8);					     
  theAsymmY7TeVGraph->SetMarkerStyle(23);
  theAsymmY7TeVGraph->GetYaxis()->SetTitle("(#Xi^{-} - #Xi^{+})/(#Xi^{-} + #Xi^{+})");
  theAsymmY7TeVGraph->GetXaxis()->SetTitle("#Xi^{-} p_{T} (GeV/c)");
  theAsymmY7TeVGraph->Draw("p same");	
  asymYLegend->AddEntry(theAsymmY7TeVGraph, "#sqrt{s} = 7 TeV Uncorrected", "pl");
  
  TAxis* theXAxis = thePosYieldsY900GeVHisto->GetXaxis();
  unsigned int theBins = theXAxis->GetNbins();
  float lowLimit = 0.;
  float hiLimit = 0.;
  
  for(unsigned int vIt = 0; vIt < theBins; ++vIt){
    if(vIt == 0) lowLimit = theXAxis->GetBinLowEdge(vIt + 1);
    double asymmetryMinus = theNegCorrYieldsY900GeVHisto->GetBinContent(vIt + 1) - thePosCorrYieldsY900GeVHisto->GetBinContent(vIt + 1);							  
    double asymmetryPlus = thePosCorrYieldsY900GeVHisto->GetBinContent(vIt + 1) + theNegCorrYieldsY900GeVHisto->GetBinContent(vIt + 1);							  
    double asymmetry = asymmetryMinus / asymmetryPlus;													  
    theCorrAsymmY900GeVGraph->SetPoint(vIt, theXAxis->GetBinCenter(vIt + 1), asymmetry);							  
    theCorrAsymmY900GeVGraph->SetPointEXlow(vIt,theXAxis->GetBinCenter(vIt + 1)-theXAxis->GetBinLowEdge(vIt + 1));					  
    theCorrAsymmY900GeVGraph->SetPointEXhigh(vIt,theXAxis->GetBinLowEdge(vIt + 1)-theXAxis->GetBinCenter(vIt + 1)+theXAxis->GetBinWidth(vIt + 1)); 
    double plusminusErrorSqrd = thePosCorrYieldsY900GeVHisto->GetBinError(vIt + 1) * thePosCorrYieldsY900GeVHisto->GetBinError(vIt + 1) + 						  
    			       theNegCorrYieldsY900GeVHisto->GetBinError(vIt + 1) * theNegCorrYieldsY900GeVHisto->GetBinError(vIt + 1);							  
    double asymmError = (sqrt(plusminusErrorSqrd)/asymmetryPlus)*sqrt(1. + asymmetry*asymmetry);							  
    theCorrAsymmY900GeVGraph->SetPointEYhigh(vIt,asymmError);													  
    theCorrAsymmY900GeVGraph->SetPointEYlow(vIt,asymmError);
    hiLimit = theXAxis->GetBinUpEdge(vIt + 1);											  

    double plusminusErrorSystSqrd = thePosCorrYieldsY900GeVHisto->GetBinContent(vIt + 1)*systRaty9*thePosCorrYieldsY900GeVHisto->GetBinContent(vIt + 1)*systRaty9 +
                                theNegCorrYieldsY900GeVHisto->GetBinContent(vIt +	1)*systRaty9*theNegCorrYieldsY900GeVHisto->GetBinContent(vIt + 1)*systRaty9;
    double asymmErrorSyst = (sqrt(plusminusErrorSystSqrd)/asymmetryPlus)*sqrt(1. + asymmetry*asymmetry);
    theCorrAsymmY900GeVGraph2->SetPoint(vIt, theXAxis->GetBinCenter(vIt + 1), asymmetry);
    theCorrAsymmY900GeVGraph2->SetPointEYhigh(vIt,asymmErrorSyst);
    theCorrAsymmY900GeVGraph2->SetPointEYlow(vIt,asymmErrorSyst); 
    
    double plusminusErrorSystUncorSqrd = thePosCorrYieldsY900GeVHisto->GetBinContent(vIt + 1)*systUncory*thePosCorrYieldsY900GeVHisto->GetBinContent(vIt + 1)*systUncory +
                                    theNegCorrYieldsY900GeVHisto->GetBinContent(vIt + 1)*systUncory*theNegCorrYieldsY900GeVHisto->GetBinContent(vIt + 1)*systUncory;
    double asymmErrorSystUncor = (sqrt(plusminusErrorSystUncorSqrd)/asymmetryPlus)*sqrt(1. + asymmetry*asymmetry);
    theCorrAsymmY900GeVGraph3->SetPoint(vIt, theXAxis->GetBinCenter(vIt + 1), asymmetry);
    theCorrAsymmY900GeVGraph3->SetPointEYhigh(vIt,asymmErrorSystUncor);
    theCorrAsymmY900GeVGraph3->SetPointEYlow(vIt,asymmErrorSystUncor); 
  }
  theCorrAsymmY900GeVGraph->GetXaxis()->SetLimits(lowLimit, hiLimit); 
  theCorrAsymmY900GeVGraph->SetMaximum(1); 					     
  theCorrAsymmY900GeVGraph->SetMinimum(-1);					     
  theCorrAsymmY900GeVGraph->SetMarkerColor(2);
  theCorrAsymmY900GeVGraph->SetMarkerStyle(21);
  theCorrAsymmY900GeVGraph->SetLineColor(2);			     
  theCorrAsymmY900GeVGraph->SetMarkerSize(0.8);					     
  theCorrAsymmY900GeVGraph->GetYaxis()->SetTitle("(#Xi^{-} - #Xi^{+})/(#Xi^{-} + #Xi^{+})");
  theCorrAsymmY900GeVGraph->GetXaxis()->SetTitle("#Xi^{-} rapidity");
  theCorrAsymmY900GeVGraph->Draw("p same");	
  asymYLegend->AddEntry(theCorrAsymmY900GeVGraph, "#sqrt{s} = 900 GeV Corrected", "pl");

  for(unsigned int vIt = 0; vIt < theBins; ++vIt){
    if(vIt == 0) lowLimit = theXAxis->GetBinLowEdge(vIt + 1);
    double asymmetryMinus = theNegYieldsY900GeVHisto->GetBinContent(vIt + 1) - thePosYieldsY900GeVHisto->GetBinContent(vIt + 1);							  
    double asymmetryPlus = thePosYieldsY900GeVHisto->GetBinContent(vIt + 1) + theNegYieldsY900GeVHisto->GetBinContent(vIt + 1);							  
    double asymmetry = asymmetryMinus / asymmetryPlus;													  
    theAsymmY900GeVGraph->SetPoint(vIt, theXAxis->GetBinCenter(vIt + 1), asymmetry);							  
    theAsymmY900GeVGraph->SetPointEXlow(vIt,theXAxis->GetBinCenter(vIt + 1)-theXAxis->GetBinLowEdge(vIt + 1));					  
    theAsymmY900GeVGraph->SetPointEXhigh(vIt,theXAxis->GetBinLowEdge(vIt + 1)-theXAxis->GetBinCenter(vIt + 1)+theXAxis->GetBinWidth(vIt + 1)); 
    double plusminusErrorSqrd = thePosYieldsY900GeVHisto->GetBinError(vIt + 1) * thePosYieldsY900GeVHisto->GetBinError(vIt + 1) + 						  
    			       theNegYieldsY900GeVHisto->GetBinError(vIt + 1) * theNegYieldsY900GeVHisto->GetBinError(vIt + 1);							  
    double asymmError = (sqrt(plusminusErrorSqrd)/asymmetryPlus)*sqrt(1. + asymmetry*asymmetry);							  
    theAsymmY900GeVGraph->SetPointEYhigh(vIt,asymmError);													  
    theAsymmY900GeVGraph->SetPointEYlow(vIt,asymmError);
    hiLimit = theXAxis->GetBinUpEdge(vIt + 1);											  
  }
  theAsymmY900GeVGraph->SetMaximum(1); 					     
  theAsymmY900GeVGraph->SetMinimum(-1);					     
  theAsymmY900GeVGraph->SetMarkerColor(6);					     
  theAsymmY900GeVGraph->SetLineColor(6);					     
  theAsymmY900GeVGraph->SetMarkerSize(0.8);					     
  theAsymmY900GeVGraph->SetMarkerStyle(22);
  theAsymmY900GeVGraph->GetYaxis()->SetTitle("(#Xi^{-} - #Xi^{+})/(#Xi^{-} + #Xi^{+})");
  theAsymmY900GeVGraph->GetXaxis()->SetTitle("#Xi^{-} rapidity");
  asymYLegend->AddEntry(theAsymmY900GeVGraph, "#sqrt{s} = 900 GeV Uncorrected", "pl");
  asymYLegend->Draw("F");
  theAsymmY900GeVGraph->Draw("p same");	
  TLine* zeroLine = new TLine(0,0,2.2,0);
  zeroLine->Draw("same");
  
  cmsPre->Draw();				     
  theAsymYCanvas->Print(outputDir+"xiAsymYUncor.eps");



  TLegend* asymYLegend2 = new TLegend(0.2, 0.2, 0.85, 0.4);  
  asymYLegend2->SetFillColor(0);
  asymYLegend2->SetTextSize(0.042);
  asymYLegend2->SetBorderSize(0);
  asymYLegend2->SetFillStyle(0);
  asymYLegend2->AddEntry(theCorrAsymmY900GeVGraph, "#sqrt{s} = 900 GeV Corrected", "pl");
  asymYLegend2->AddEntry(theCorrAsymmY7TeVGraph, "#sqrt{s} = 7 TeV Corrected", "pl");

  theCorrAsymmY900GeVGraph2->SetLineColor(2);
  theCorrAsymmY900GeVGraph2->SetFillColor(kRed-9);				     
  theCorrAsymmY900GeVGraph2->SetMarkerColor(2);
  theCorrAsymmY900GeVGraph2->SetMaximum(1); 					     
  theCorrAsymmY900GeVGraph2->SetMinimum(-1);	
  theCorrAsymmY900GeVGraph2->SetMarkerStyle(21);
  theCorrAsymmY900GeVGraph2->SetLineColor(2);			     
  theCorrAsymmY900GeVGraph2->SetMarkerSize(0.8);			
  theCorrAsymmY900GeVGraph2->GetYaxis()->SetTitle("(#Xi^{-} - #Xi^{+})/(#Xi^{-} + #Xi^{+})");
  theCorrAsymmY900GeVGraph2->GetXaxis()->SetTitle("#Xi^{-} rapidity");
  theCorrAsymmY900GeVGraph2->Draw("e3,ap");
  theCorrAsymmY900GeVGraph->SetMarkerColor(2);
  theCorrAsymmY900GeVGraph->Draw("e1,p");
  theCorrAsymmY900GeVGraph3->SetLineColor(2);
  theCorrAsymmY900GeVGraph3->SetMarkerColor(2);
  theCorrAsymmY900GeVGraph3->Draw("p");

  theCorrAsymmY7TeVGraph2->SetLineColor(4);
  theCorrAsymmY7TeVGraph2->SetFillColor(38);				     
  theCorrAsymmY7TeVGraph2->SetMarkerColor(4);
  theCorrAsymmY7TeVGraph2->SetMarkerStyle(20);
  theCorrAsymmY7TeVGraph2->SetLineColor(4);			     
  theCorrAsymmY7TeVGraph2->SetMarkerSize(0.8);			
  theCorrAsymmY7TeVGraph2->Draw("e3,p");
  theCorrAsymmY7TeVGraph->SetMarkerColor(4);
  theCorrAsymmY7TeVGraph->Draw("e1,p");
  theCorrAsymmY7TeVGraph3->SetLineColor(4);
  theCorrAsymmY7TeVGraph3->SetMarkerColor(4);
  theCorrAsymmY7TeVGraph3->Draw("p");

  asymYLegend2->Draw("F");
  cmsPre->Draw();	
  theAsymYCanvas->Print(outputDir+"xiAsymY.eps");

  outTextFile << " Should really calculate data/MC and 900/7 results ratios here for inclusion in the note" << endl;

  outTextFile.close();



/*
Plots needed for AN:


convert XiEffCtau_data_7.eps XiEffCtau_data_7.pdf	    
convert systVar_900GeVData_looseDelR_nominal/cTauBinsFit.eps XiMassCtauBins_data_900.pdf   
convert systVar_900GeVData_looseDelR_nominal/rapidityBinsFit.eps XiMassYBins_data_900.pdf
convert systVar_900GeVData_looseDelR_nominal/ptBinsFit.eps XiMasspTBins_data_900.pdf
convert pt_NSD_dataMC.eps pt_NSD_dataMC.pdf    
convert y_NSD_dataMC.eps y_NSD_dataMC.pdf
convert XiEffCtau_data_900.eps XiEffCtau_data_900.pdf     
convert systVar_7TeVMC_looseDelR_nominal/XiMass_selected_MC_7TeV.eps XiMass_selected_MC_7.pdf      
convert pt_NSD_ratio.eps pt_NSD_ratio.pdf     
convert y_NSD_ratio.eps y_NSD_ratio.pdf
convert XiEffY_data_7.eps XiEffY_data_7.pdf	    
convert systVar_900GeVMC_looseDelR_nominal/XiMass_selected_MC_900GeV.eps XiMass_selected_MC_900.pdf   
convert pt_fitval_ratio.eps pt_fitval_ratio.pdf  
convert y_fitval_ratio.eps y_fitval_ratio.pdf
convert XiEffY_data_900.eps XiEffY_data_900.pdf	    
convert systVar_7TeVData_looseDelR_nominal/XiMassSelected.eps XiMass_selected_data_7.pdf    
convert pt_sigma.eps pt_sigma.pdf	      
convert y_sigma.eps y_sigma.pdf
convert XiEffpT_data_7.eps XiEffpT_data_7.pdf	    
convert systVar_900GeVData_looseDelR_nominal/XiMassSelected.eps XiMass_selected_data_900.pdf  
convert xiAsymPtUncor.eps xiAsymPtUncor.pdf
convert xiAsymPt.eps xiAsymPt.pdf
convert xiAsymY.eps xiAsymY.pdf
convert xiAsymYUncor.eps xiAsymYUncor.pdf
convert XiEffpT_data_900.eps XiEffpT_data_900.pdf	    
convert ct_fitval_ratio.eps ct_fitval_ratio.pdf	       
convert systVar_7TeVData_looseDelR_nominal/cTauBinsFit.eps XiMassCtauBins_data_7.pdf  
convert systVar_7TeVData_looseDelR_nominal/rapidityBinsFit.eps XiMassYBins_data_7.pdf
convert systVar_7TeVData_looseDelR_nominal/ptBinsFit.eps XiMasspTBins_data_7.pdf
convert pt_NSD_dataFits.eps pt_NSD_dataFits.pdf	       
convert y_NSD_data.eps y_NSD_data.pdf
convert systVar_7TeVMC_looseDelR_nominal/XiTruthComp_MC_7TeV.eps XiTruthComp_MC_7.pdf
convert systVar_900GeVMC_looseDelR_nominal/XiTruthComp_MC_900GeV.eps XiTruthComp_MC_900.pdf

convert XiCuts_pt.eps XiCuts_pt.pdf
convert XiCuts_y.eps XiCuts_y.pdf
convert xiAllCandsComp.eps xiAllCandsComp.pdf

*/


}
