#include <iostream>
#include "run16dAupc3matching.h"

void fitwitherr2() {

  fetchpc3dphidz();

  float pt[50];
  float errx[50];

  for (int ipt=0; ipt<50; ipt++) {
    pt[ipt] = ipt*0.1 + 0.05;
    errx[ipt] = 0.0;
  }

  gStyle->SetOptStat(kFALSE);
  //gStyle->SetOptFit(1101);

  ofstream fout("run16dAupc3matchingsmooth.h");
  ofstream fout2("run16dAupc3matchingsmooth2.h");


  for (int iarm=0; iarm<2; iarm++) {
    for (int ich=0; ich<2; ich++) {
    for(int ivz=0; ivz<5; ivz++){
      string arm;
      string ch;
      if (iarm == 0) arm = "east";
      else if (iarm == 1) arm = "west";
      else arm = "err";

      if (ich == 0) ch = "pos";
      else if (ich == 1) ch = "neg";
      else ch = "err";
        string bbcz;
	if (ivz==0) bbcz = "-10 to -6";
	else if (ivz==1) bbcz = "-6 to -2";
	else if (ivz==2) bbcz = "-2 to 2";
	else if (ivz==3) bbcz = "2 to 6";
	else if (ivz==4) bbcz = "6 to 10";
	else bbcz = "err";

      TCanvas *c1 = new TCanvas("c1","c1",500,500);
      c1->SetGridx();
      TGraphErrors* dphisigma = new TGraphErrors(50,pt,pc3dphisigma[iarm][ich][ivz],errx,pc3dphisigmaerr[iarm][ich][ivz]);
      dphisigma->SetTitle(Form("dphisigma_%s_%s_%s",arm,ch,bbcz));
      dphisigma->SetLineColor(1);
      dphisigma->SetMarkerStyle(20);
      dphisigma->SetMarkerSize(0.8);
      dphisigma->GetXaxis()->SetRangeUser(0.2,5.0);
      dphisigma->GetHistogram()->SetMaximum(0.01);
      dphisigma->GetHistogram()->SetMinimum(-0.01);
      dphisigma->GetXaxis()->SetTitle("p_{T}");
      dphisigma->GetYaxis()->SetTitle("dphi sigma");
      dphisigma->Draw("AP");

      TPaveText *p1 = new TPaveText(0.2,0.8,0.9,0.9,"NDC");
      p1->AddText("Fit function: a+bx+cx^{2}+dx^{3}+ex^{4}+fx^{5}+#frac{g}{#sqrt{x}}+#frac{h}{x^{2}}");
      p1->Draw("same");

      TF1 *fdphisigma = new TF1("fdphisigma","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x+[6]/TMath::Sqrt(x)+[7]/x/x",0.3,5.0);
      //TF1 *fdphisigma2 = new TF1("fdphisigma2","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x+[6]/TMath::Sqrt(x)+[7]/x/x",0.2,5.0);
      dphisigma->Fit("fdphisigma","Q0","",0.4,4.0);
      //double sigma_par[8];
      //fdphisigma->GetParameters(sigma_par);
      //fdphisigma->SetParameters(sigma_par);
      //dphisigma->Fit("fdphisigma2","Q0","",0.3,5.0);
      fdphisigma->Draw("same");
      fout << "fdphisigma->SetParameters(" << fdphisigma->GetParameter(0) << "," << fdphisigma->GetParameter(1) << "," << fdphisigma->GetParameter(2) << "," << fdphisigma->GetParameter(3) << "," << fdphisigma->GetParameter(4) << "," << fdphisigma->GetParameter(5) << "," << fdphisigma->GetParameter(6) << "," << fdphisigma->GetParameter(7) << ");" << endl;
      for(int ipar=0;ipar<8;ipar++){
      fout2 << "PC3_dphisigma[" << iarm <<"][" << ich << "][" << ivz <<"][" << ipar << "] = " << fdphisigma->GetParameter(ipar) << ";" << endl;
      }
      c1->Print(Form("smooth/dphisigma_%d_%d_%d.png",iarm,ich,ivz));
      delete c1;


      TCanvas *c3 = new TCanvas("c3","c3",500,500);
      c3->SetGridx();
      TGraphErrors* dphimean = new TGraphErrors(50,pt,pc3dphimean[iarm][ich][ivz],errx,pc3dphimeanerr[iarm][ich][ivz]);
      dphimean->SetTitle(Form("dphimean_%s_%s_%s",arm,ch,bbcz));
      dphimean->SetLineColor(1);
      dphimean->SetMarkerStyle(20);
      dphimean->SetMarkerSize(0.8);
      dphimean->GetXaxis()->SetRangeUser(0.2,5.0);
      dphimean->GetHistogram()->SetMaximum(0.01);
      dphimean->GetHistogram()->SetMinimum(-0.01);
      dphimean->GetXaxis()->SetTitle("p_{T}");
      dphimean->GetYaxis()->SetTitle("dphi mean");
      dphimean->Draw("AP");

      TPaveText *p3 = new TPaveText(0.2,0.8,0.9,0.9,"NDC");
      p3->AddText("Fit function: a+bx+#frac{c}{x}+#frac{d}{#sqrt{x}}+#frac{e}{x^{2}}+#frac{f}{x^{3}}+#frac{g}{x^{4}}");
      p3->Draw("same");
      TF1 *fdphimean = new TF1("fdphimean","[0]+[1]*x+[2]/x+[3]/TMath::Sqrt(x)+[4]/x/x+[5]/x/x/x+[6]/x/x/x/x",0.3,5.0);
      //TF1 *fdphimean2 = new TF1("fdphimean2","[0]+[1]*x+[2]/x+[3]/TMath::Sqrt(x)+[4]/x/x+[5]/x/x/x+[6]/x/x/x/x",0.2,5.0);
      dphimean->Fit("fdphimean","Q0","",0.4,4.0);
      double mean_par[7];
      //fdphimean->GetParameters(mean_par);
      //fdphimean2->SetParameters(mean_par);
      //dphimean->Fit("fdphimean2","Q0","",0.3,5.0);
      fdphimean->Draw("same");
      fout << "fdphimean->SetParameters(" << fdphimean->GetParameter(0) << "," << fdphimean->GetParameter(1) << "," << fdphimean->GetParameter(2) << "," << fdphimean->GetParameter(3) << "," << fdphimean->GetParameter(4) << "," << fdphimean->GetParameter(5) << "," << fdphimean->GetParameter(6) << ");" << endl;
      for(int ipar=0;ipar<7;ipar++){
      fout2 << "PC3_dphimean[" << iarm <<"][" << ich << "][" << ivz <<"][" << ipar << "] = " << fdphimean->GetParameter(ipar) << ";" << endl;
      }
      c3->Print(Form("smooth/dphimean_%d_%d_%d.png",iarm,ich,ivz));
      delete c3;

      TCanvas *c2 = new TCanvas("c2","c2",500,500);
      c2->SetGridx();
      TGraphErrors* dzsigma = new TGraphErrors(50,pt,pc3dzsigma[iarm][ich][ivz],errx,pc3dzsigmaerr[iarm][ich][ivz]);
      dzsigma->SetTitle(Form("dzsigma_%s_%s_%s",arm,ch,bbcz));
      dzsigma->SetLineColor(1);
      dzsigma->SetMarkerStyle(20);
      dzsigma->SetMarkerSize(0.8);
      dzsigma->GetXaxis()->SetRangeUser(0.2,5.0);
      dzsigma->GetHistogram()->SetMaximum(5);
      dzsigma->GetHistogram()->SetMinimum(0);
      dzsigma->GetXaxis()->SetTitle("p_{T}");
      dzsigma->GetYaxis()->SetTitle("dz sigma");
      dzsigma->Draw("AP");

      TPaveText *p2 = new TPaveText(0.2,0.8,0.9,0.9,"NDC");
      p2->AddText("Fit function: a+bx+cx^{2}+dx^{3}+ex^{4}+fx^{5}+#frac{g}{#sqrt{x}}+#frac{h}{x^{2}}");
      p2->Draw("same");

      //TF1 *fdzsigma = new TF1("fdzsigma","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]/TMath::Sqrt(x)+[6]/x/x",0.3,5.0);
      TF1 *fdzsigma = new TF1("fdzsigma","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x+[6]/TMath::Sqrt(x)+[7]/x/x",0.3,5.0);
      //TF1 *fdzsigma2 = new TF1("fdzsigma2","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]/TMath::Sqrt(x)+[6]/x/x",0.2,5.0);
      dzsigma->Fit("fdzsigma","Q0","",0.4,4.0);
      //fdzsigma->GetParameters(sigma_par);
      //fdzsigma2->SetParameters(sigma_par);
      //dzsigma->Fit("fdzsigma2","Q0","",0.3,5.0);
      fdzsigma->Draw("same");
      fout << "fdzsigma->SetParameters(" << fdzsigma->GetParameter(0) << "," << fdzsigma->GetParameter(1) << "," << fdzsigma->GetParameter(2) << "," << fdzsigma->GetParameter(3) << "," << fdzsigma->GetParameter(4) << "," << fdzsigma->GetParameter(5) << "," << fdzsigma->GetParameter(6) << ");" << endl;
      for(int ipar=0;ipar<8;ipar++){
      fout2 << "PC3_dzsigma[" << iarm <<"][" << ich << "][" << ivz <<"][" << ipar << "] = " << fdzsigma->GetParameter(ipar) <<  ";" << endl;
      }
      c2->Print(Form("smooth/dzsigma_%d_%d_%d.png",iarm,ich,ivz));
      delete c2;





      TCanvas *c4 = new TCanvas("c4","c4",500,500);
      c4->SetGridx();
      TGraphErrors* dzmean = new TGraphErrors(50,pt,pc3dzmean[iarm][ich][ivz],errx,pc3dzmeanerr[iarm][ich][ivz]);
      dzmean->SetTitle(Form("dzmean_%s_%s_%s",arm,ch,bbcz));
      dzmean->SetLineColor(1);
      dzmean->SetMarkerStyle(20);
      dzmean->SetMarkerSize(0.8);
      dzmean->GetXaxis()->SetRangeUser(0.2,5.0);
      dzmean->GetHistogram()->SetMaximum(3);
      dzmean->GetHistogram()->SetMinimum(-2);
      dzmean->GetXaxis()->SetTitle("p_{T}");
      dzmean->GetYaxis()->SetTitle("dz mean");
      dzmean->Draw("AP");

      TPaveText *p4 = new TPaveText(0.2,0.8,0.9,0.9,"NDC");
      p4->AddText("Fit function: a+bx+#frac{c}{x}+#frac{d}{#sqrt{x}}+#frac{e}{x^{2}}+#frac{f}{x^{3}}+#frac{g}{x^{4}}");
      p4->Draw("same");

      TF1 *fdzmean = new TF1("fdzmean","[0]+[1]*x+[2]/x+[3]/TMath::Sqrt(x)+[4]/x/x+[5]/x/x/x+[6]/x/x/x/x",0.3,5.0);
      //TF1 *fdzmean2 = new TF1("fdzmean2","[0]+[1]*x+[2]/x+[3]/TMath::Sqrt(x)+[4]/x/x+[5]/x/x/x+[6]/x/x/x/x",0.2,5.0);
      dzmean->Fit("fdzmean","Q0","",0.4,4.0);
      //fdzmean->GetParameters(mean_par);
      //fdzmean2->SetParameters(mean_par);
      //dzmean->Fit("fdzmean2","Q0","",0.3,5.0);
      fdzmean->Draw("same");
      fout << "fdzmean->SetParameters(" << fdzmean->GetParameter(0) << "," << fdzmean->GetParameter(1) << "," << fdzmean->GetParameter(2) << "," << fdzmean->GetParameter(3) << "," << fdzmean->GetParameter(4) << "," << fdzmean->GetParameter(5) << "," << fdzmean->GetParameter(6) << ");" << endl;
      for(int ipar=0;ipar<7;ipar++){
      fout2 << "PC3_dzmean[" << iarm <<"][" << ich << "][" << ivz <<"][" << ipar << "] = " << fdzmean->GetParameter(ipar) <<  ";" << endl;
      }
      c4->Print(Form("smooth/dzmean_%d_%d_%d.png",iarm,ich,ivz));
      delete c4;
      }
    }
  }
  

}








