#include <iostream>
#include "run16dAupc3matching.h"

void fitwitherr() {

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


  for (int iarm=0; iarm<2; iarm++) {
    for (int ich=0; ich<2; ich++) {



      TCanvas *c1 = new TCanvas("c1","c1",500,500);
      c1->SetGridx();
      TGraphErrors* dphisigma[10];
      TLegend *leg1 = new TLegend(0.15,0.1,0.4,0.35);
      leg1->SetBorderSize(0);
      leg1->SetTextSize(0.04);
      for (int ivz=0; ivz<10; ivz++) {
      
        string arm;
        string ch;
	string bbcz;
        if (iarm == 0) arm = "east";
        else if (iarm == 1) arm = "west";
        else arm = "err";

        if (ich == 0) ch = "pos";
        else if (ich == 1) ch = "neg";
        else ch = "err";

	if (ivz==0) bbcz = "-10 to -8";
        else if (ivz==1) bbcz = "-8 to -6";
	else if (ivz==2) bbcz = "-6 to -4";
	else if (ivz==3) bbcz = "-4 to -2";
	else if (ivz==4) bbcz = "-2 to 0";
	else if (ivz==5) bbcz = "0 to 2";
	else if (ivz==6) bbcz = "2 to 4";
	else if (ivz==7) bbcz = "4 to 6";
	else if (ivz==8) bbcz = "6 to 8";
	else if (ivz==9) bbcz = "8 to 10";
	else bbcz = "err";

        dphisigma[ivz] = new TGraphErrors(50,pt,pc3dphisigma[iarm][ich][ivz],errx,pc3dphisigmaerr[iarm][ich][ivz]);
	dphisigma[ivz]->SetName(Form("dphisigma_%s_%s_%d",arm,ch,ivz));
        dphisigma[ivz]->SetTitle(Form("dphisigma_%s_%s_%d",arm,ch,ivz));
        dphisigma[ivz]->SetLineColor(0);
        dphisigma[ivz]->SetMarkerStyle(1);
	dphisigma[ivz]->SetMarkerColor(0);
        dphisigma[ivz]->SetMarkerSize(0.8);
        dphisigma[ivz]->GetXaxis()->SetRangeUser(0.3,5.0);
        dphisigma[ivz]->GetHistogram()->SetMaximum(0.01);
        dphisigma[ivz]->GetHistogram()->SetMinimum(-0.01);
        dphisigma[ivz]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        dphisigma[ivz]->GetYaxis()->SetTitle("dphi sigma");

        if (ivz==0) {
	  dphisigma[ivz]->Draw("AP");
	}
	else if {
	  dphisigma[ivz]->Draw("CP");
	}

        TPaveText *p1 = new TPaveText(0.2,0.8,0.9,0.9,"NDC");
        p1->AddText("Fit function: a+bx+cx^{2}+dx^{3}+ex^{4}+fx^{5}+#frac{g}{#sqrt{x}}+#frac{h}{x^{2}}");
        p1->Draw("same");

        TF1 *fdphisigma = new TF1("fdphisigma","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x+[6]/TMath::Sqrt(x)+[7]/x/x",0.3,5.0);
        //TF1 *fdphisigma2 = new TF1("fdphisigma2","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x+[6]/TMath::Sqrt(x)+[7]/x/x",0.2,5.0);
        dphisigma[ivz]->Fit("fdphisigma","Q0","",0.4,.0);
        double sigma_par[10];
        //fdphisigma->GetParameters(sigma_par);
        //fdphisigma2->SetParameters(sigma_par);
        //dphisigma->Fit("fdphisigma2","Q0","",0.3,5.0);
	fdphisigma->SetLineColor(ivz+1);
	leg1->AddEntry(fdphisigma,Form("vtex zed bin %s",bbcz),"l");
        fdphisigma->Draw("same");
        fout << "fdphisigma->SetParameters(" << fdphisigma->GetParameter(0) << "," << fdphisigma->GetParameter(1) << "," << fdphisigma->GetParameter(2) << "," << fdphisigma->GetParameter(3) << "," << fdphisigma->GetParameter(4) << "," << fdphisigma->GetParameter(5) << "," << fdphisigma->GetParameter(6) << "," << fdphisigma->GetParameter(7) << ");" << endl;
      }
      leg1->Draw("same");
      c1->Print(Form("smooth/dphisigma_%d_%d.png",iarm,ich));
      delete c1;




      TCanvas *c3 = new TCanvas("c3","c3",500,500);
      c3->SetGridx();
      TGraphErrors* dphimean[10];
      TLegend *leg3 = new TLegend(0.15,0.1,0.4,0.35);
      leg3->SetBorderSize(0);
      leg3->SetTextSize(0.04);
      for (int ivz=0; ivz<10; ivz++) {
      
        string arm;
        string ch;
	string bbcz;
        if (iarm == 0) arm = "east";
        else if (iarm == 1) arm = "west";
        else arm = "err";

        if (ich == 0) ch = "pos";
        else if (ich == 1) ch = "neg";
        else ch = "err";

	if (ivz==0) bbcz = "-10 to -8";
        else if (ivz==1) bbcz = "-8 to -6";
	else if (ivz==2) bbcz = "-6 to -4";
	else if (ivz==3) bbcz = "-4 to -2";
	else if (ivz==4) bbcz = "-2 to 0";
	else if (ivz==5) bbcz = "0 to 2";
	else if (ivz==6) bbcz = "2 to 4";
	else if (ivz==7) bbcz = "4 to 6";
	else if (ivz==8) bbcz = "6 to 8";
	else if (ivz==9) bbcz = "8 to 10";
	else bbcz = "err";
        dphimean[ivz] = new TGraphErrors(50,pt,pc3dphimean[iarm][ich][ivz],errx,pc3dphimeanerr[iarm][ich][ivz]);
        dphimean[ivz]->SetTitle(Form("dphimean_%s_%s",arm,ch));
        dphimean[ivz]->SetLineColor(0);
        dphimean[ivz]->SetMarkerStyle(1);
        dphimean[ivz]->SetMarkerSize(0.8);
        dphimean[ivz]->GetXaxis()->SetRangeUser(0.3,5.0);
        dphimean[ivz]->GetHistogram()->SetMaximum(0.01);
        dphimean[ivz]->GetHistogram()->SetMinimum(-0.01);
        dphimean[ivz]->GetXaxis()->SetTitle("p_{T}/GeV/c");
        dphimean[ivz]->GetYaxis()->SetTitle("dphi mean");

	if (ivz==0) {
	  dphimean[ivz]->Draw("AP");
	}
	else {
	  dphimean[ivz]->Draw("CP");
	}

        TPaveText *p3 = new TPaveText(0.2,0.8,0.9,0.9,"NDC");
        p3->AddText("Fit function: a+bx+#frac{c}{x}+#frac{d}{#sqrt{x}}+#frac{e}{x^{2}}+#frac{f}{x^{3}}+#frac{g}{x^{4}}");
        p3->Draw("same");
        TF1 *fdphimean = new TF1("fdphimean","[0]+[1]*x+[2]/x+[3]/TMath::Sqrt(x)+[4]/x/x+[5]/x/x/x+[6]/x/x/x/x",0.3,5.0);
        //TF1 *fdphimean2 = new TF1("fdphimean2","[0]+[1]*x+[2]/x+[3]/TMath::Sqrt(x)+[4]/x/x+[5]/x/x/x+[6]/x/x/x/x",0.2,5.0);
        dphimean[ivz]->Fit("fdphimean","Q0","",0.4,.0);
        //double mean_par[7];
        //fdphimean->GetParameters(mean_par);
        //fdphimean2->SetParameters(mean_par);
        //dphimean->Fit("fdphimean2","Q0","",0.3,5.0);
	fdphimean->SetLineColor(ivz+1);
	leg3->AddEntry(fdphimean,Form("dz bin %s",bbcz),"l");
        fdphimean->Draw("same");
        fout << "fdphimean->SetParameters(" << fdphimean->GetParameter(0) << "," << fdphimean->GetParameter(1) << "," << fdphimean->GetParameter(2) << "," << fdphimean->GetParameter(3) << "," << fdphimean->GetParameter(4) << "," << fdphimean->GetParameter(5) << "," << fdphimean->GetParameter(6) << ");" << endl;
      }
      leg3->Draw("same");
      c3->Print(Form("smooth/dphimean_%d_%d.png",iarm,ich));
      delete c3;

      
      TCanvas *c2 = new TCanvas("c2","c2",500,500);
      c2->SetGridx();
      TGraphErrors* dzsigma[10];
      TLegend *leg2 = new TLegend(0.15,0.1,0.4,0.35);
      leg2->SetBorderSize(0);
      leg2->SetTextSize(0.04);
      for (int ivz=0; ivz<10; ivz++) {
      
        string arm;
        string ch;
	string bbcz;
        if (iarm == 0) arm = "east";
        else if (iarm == 1) arm = "west";
        else arm = "err";

        if (ich == 0) ch = "pos";
        else if (ich == 1) ch = "neg";
        else ch = "err";

	if (ivz==0) bbcz = "-10 to -8";
        else if (ivz==1) bbcz = "-8 to -6";
	else if (ivz==2) bbcz = "-6 to -4";
	else if (ivz==3) bbcz = "-4 to -2";
	else if (ivz==4) bbcz = "-2 to 0";
	else if (ivz==5) bbcz = "0 to 2";
	else if (ivz==6) bbcz = "2 to 4";
	else if (ivz==7) bbcz = "4 to 6";
	else if (ivz==8) bbcz = "6 to 8";
	else if (ivz==9) bbcz = "8 to 10";
	else bbcz = "err";
        dzsigma[ivz] = new TGraphErrors(50,pt,pc3dzsigma[iarm][ich][ivz],errx,pc3dzsigmaerr[iarm][ich][ivz]);
        dzsigma[ivz]->SetTitle(Form("dzsigma_%s_%s",arm,ch));
        dzsigma[ivz]->SetLineColor(0);
        dzsigma[ivz]->SetMarkerStyle(1);
        dzsigma[ivz]->SetMarkerSize(0.8);
        dzsigma[ivz]->GetXaxis()->SetRangeUser(0.3,5.0);
        dzsigma[ivz]->GetHistogram()->SetMaximum(5);
        dzsigma[ivz]->GetHistogram()->SetMinimum(0);
        dzsigma[ivz]->GetXaxis()->SetTitle("p_{T}");
        dzsigma[ivz]->GetYaxis()->SetTitle("dz sigma");
        if (ivz==0) {
	  dzsigma[ivz]->Draw("AP");
	}
	else {
	  dzsigma[ivz]->Draw("CP");
	}

        TPaveText *p2 = new TPaveText(0.2,0.8,0.9,0.9,"NDC");
        p2->AddText("Fit function: a+bx+cx^{2}+dx^{3}+ex^{4}+fx^{5}+#frac{g}{#sqrt{x}}+#frac{h}{x^{2}}");
        p2->Draw("same");

        TF1 *fdzsigma = new TF1("fdzsigma","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]/TMath::Sqrt(x)+[6]/x/x",0.3,5.0);
        //TF1 *fdzsigma2 = new TF1("fdzsigma2","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]/TMath::Sqrt(x)+[6]/x/x",0.2,5.0);
        dzsigma[ivz]->Fit("fdzsigma","Q0","",0.4,.0);
	fdzsigma->SetLineColor(ivz+1);
	leg2->AddEntry(fdzsigma,Form("vtx zed bin %s",bbcz),"l");
        //fdzsigma->GetParameters(sigma_par);
        //fdzsigma2->SetParameters(sigma_par);
        //dzsigma->Fit("fdzsigma2","Q0","",0.3,5.0);
        fdzsigma->Draw("same");
        fout << "fdzsigma->SetParameters(" << fdzsigma->GetParameter(0) << "," << fdzsigma->GetParameter(1) << "," << fdzsigma->GetParameter(2) << "," << fdzsigma->GetParameter(3) << "," << fdzsigma->GetParameter(4) << "," << fdzsigma->GetParameter(5) << "," << fdzsigma->GetParameter(6) << ");" << endl;
      }
      leg2->Draw("same");
      c2->Print(Form("smooth/dzsigma_%d_%d.png",iarm,ich));
      delete c2;





      TCanvas *c4 = new TCanvas("c4","c4",500,500);
      c4->SetGridx();
      TGraphErrors* dzmean[10];
      TLegend *leg4 = new TLegend(0.15,0.1,0.4,0.35);
      leg4->SetBorderSize(0);
      leg4->SetTextSize(0.04);
      for (int ivz=0; ivz<10; ivz++) {
      
        string arm;
        string ch;
	string bbcz;
        if (iarm == 0) arm = "east";
        else if (iarm == 1) arm = "west";
        else arm = "err";

        if (ich == 0) ch = "pos";
        else if (ich == 1) ch = "neg";
        else ch = "err";

	if (ivz==0) bbcz = "-10 to -8";
        else if (ivz==1) bbcz = "-8 to -6";
	else if (ivz==2) bbcz = "-6 to -4";
	else if (ivz==3) bbcz = "-4 to -2";
	else if (ivz==4) bbcz = "-2 to 0";
	else if (ivz==5) bbcz = "0 to 2";
	else if (ivz==6) bbcz = "2 to 4";
	else if (ivz==7) bbcz = "4 to 6";
	else if (ivz==8) bbcz = "6 to 8";
	else if (ivz==9) bbcz = "8 to 10";
	else bbcz = "err";

        dzmean[ivz] = new TGraphErrors(50,pt,pc3dzmean[iarm][ich][ivz],errx,pc3dzmeanerr[iarm][ich][ivz]);
        dzmean[ivz]->SetTitle(Form("dzmean_%s_%s",arm,ch));
        dzmean[ivz]->SetLineColor(0);
        dzmean[ivz]->SetMarkerStyle(1);
        dzmean[ivz]->SetMarkerSize(0.8);
        dzmean[ivz]->GetXaxis()->SetRangeUser(0.3,5.0);
        dzmean[ivz]->GetHistogram()->SetMaximum(3);
        dzmean[ivz]->GetHistogram()->SetMinimum(-2);
        dzmean[ivz]->GetXaxis()->SetTitle("p_{T}/GeV/c");
        dzmean[ivz]->GetYaxis()->SetTitle("dz mean");
        if(ivz==0) dzmean[ivz]->Draw("AP");
	else dzmean[ivz]->Draw("CP");

        TPaveText *p4 = new TPaveText(0.2,0.8,0.9,0.9,"NDC");
        p4->AddText("Fit function: a+bx+#frac{c}{x}+#frac{d}{#sqrt{x}}+#frac{e}{x^{2}}+#frac{f}{x^{3}}+#frac{g}{x^{4}}");
        p4->Draw("same");

        TF1 *fdzmean = new TF1("fdzmean","[0]+[1]*x+[2]/x+[3]/TMath::Sqrt(x)+[4]/x/x+[5]/x/x/x+[6]/x/x/x/x",0.3,5.0);
        //TF1 *fdzmean2 = new TF1("fdzmean2","[0]+[1]*x+[2]/x+[3]/TMath::Sqrt(x)+[4]/x/x+[5]/x/x/x+[6]/x/x/x/x",0.2,5.0);
        dzmean[ivz]->Fit("fdzmean","Q0","",0.4,.0);
	fdzmean->SetLineColor(ivz+1);
	leg4->AddEntry(fdzmean,Form("vtx zed bin %s",bbcz),"l");
        //fdzmean->GetParameters(mean_par);
        //fdzmean2->SetParameters(mean_par);
        //dzmean->Fit("fdzmean2","Q0","",0.3,5.0);
        fdzmean->Draw("same");
        fout << "fdzmean->SetParameters(" << fdzmean->GetParameter(0) << "," << fdzmean->GetParameter(1) << "," << fdzmean->GetParameter(2) << "," << fdzmean->GetParameter(3) << "," << fdzmean->GetParameter(4) << "," << fdzmean->GetParameter(5) << "," << fdzmean->GetParameter(6) << ");" << endl;
      }
      leg4->Draw("same");
      c4->Print(Form("smooth/dzmean_%d_%d.png",iarm,ich));
      delete c4;
    }
  }
}








