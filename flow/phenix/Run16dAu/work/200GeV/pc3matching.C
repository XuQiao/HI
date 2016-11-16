#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "TAttMarker.h"
#include "TFile.h"

using namespace std;

void pc3matching() {

  TFile *f = TFile::Open("output_perform.root");

  ofstream fout("run16dAupc3matching.h");
  ofstream fout2("run16dAupc3matchingfirst.h");
  fout << "float pc3dphimean[2][2][10][50];" << endl;
  fout << "float pc3dphisigma[2][2][10][50];" << endl;
  fout << "float pc3dzmean[2][2][10][50];" << endl;
  fout << "float pc3dzsigma[2][2][10][50];" << endl;
  fout << "float pc3dphimeanerr[2][2][10][50];" << endl;
  fout << "float pc3dphisigmaerr[2][2][10][50];" << endl;
  fout << "float pc3dzmeanerr[2][2][10][50];" << endl;
  fout << "float pc3dzsigmaerr[2][2][10][50];" << endl;
  fout << " " << endl;
  fout << "void fetchpc3dphidz();" << endl;
  fout << " " << endl;
  fout << "void fetchpc3dphidz() {" << endl;


  float max = 0.0;
  float sigma = 0.0;
  float mean = 0.0;

  for (Int_t iarm = 0; iarm < 2; iarm++) {
    for (Int_t ich = 0; ich < 2; ich++) {
      for (Int_t ipt = 0; ipt < 50; ipt++) {
        for(Int_t ivz = 0; ivz < 10; ivz++) {
           // if(ipt!=3 || iarm!=0 || ich!=0 || ivz!=0)continue;
            cout<<iarm<<" "<<ich<<" "<<ipt<<" "<<ivz<<endl;
          double sigmaerr=0.0;
	  double meanerr=0.0;
          TString ch = "";
            if(ich == 0) ch = "pos";
            if(ich == 1) ch = "neg";
	  TString histname = Form("pc3dphidz_arm%d_%s_z%d_%d",iarm,ch.Data(),ivz,ipt);
          TH2D *hist = (TH2D*) f->Get(histname);

          TH1D *dphi = (TH1D*) hist->ProjectionX(Form("pc3dphi_%d_%d_%d_%d",iarm,ich,ipt,ivz));
          TH1D *dz = (TH1D*) hist->ProjectionY(Form("pc3dz_%d_%d_%d_%d",iarm,ich,ipt,ivz));

	  dphi->GetXaxis()->SetRangeUser(-0.1,0.1);

          gStyle->SetOptFit(1101);

          TF1 *fphi1 = new TF1("fphi1","gaus",-0.1,0.1);
          TF1 *fz1 = new TF1("fz1","gaus",-10,10);
          TF1 *fphi2 = new TF1("fphi2","gaus(0)+gaus(3)",-0.1,0.1);
          TF1 *fz2 = new TF1("fz2","gaus(0)+gaus(3)",-10,10);
          TF1 *phi_gaus1 = new TF1("phi_gaus1","gaus",-0.1,0.1);
          TF1 *phi_gaus2 = new TF1("phi_gaus2","gaus",-0.1,0.1);
          TF1 *z_gaus1 = new TF1("z_gaus1","gaus",-10,10);
          TF1 *z_gaus2 = new TF1("z_gaus2","gaus",-10,10);


          Float_t Xbins = dphi->GetNbinsX();
          Float_t Xmin = dphi->GetXaxis()->GetXmin();
          Float_t Xmax = dphi->GetXaxis()->GetXmax();
          mean = (dphi->GetMaximumBin() * (Xmax-Xmin))/Xbins + Xmin;
          max = dphi->GetMaximum();
          fphi1->SetRange(mean-0.01,mean+0.01);
	  fphi1->SetParameters(max,mean);

          TCanvas *c = new TCanvas("c","c",500,500);
	  //dphi->Scale(1./dphi->Integral());
	  dphi->SetTitle("dphi matching");
	  dphi->GetXaxis()->SetTitle("dphi");
	  dphi->GetYaxis()->SetTitle("# tracks");
	  if (iarm == 1 && ich ==0) dphi->Rebin(2);
	//  dphi->Rebin(2);
	  dphi->SetMarkerSize(1);
	  dphi->SetMarkerStyle(kFullCircle);
          dphi->Draw("P");
          dphi->Fit("fphi1","RQ0");
          double dphi_par[6];
          fphi1->GetParameters(dphi_par);
          dphi_par[3] = 0.1*dphi_par[0];
          dphi_par[4] = dphi_par[1];
          dphi_par[5] = 8*dphi_par[2];
          fphi2->SetParameters(dphi_par);
	  fphi2->SetParLimits(3,0,5*dphi_par[3]);
	  fphi2->SetParLimits(4,-1,1);
	  fphi2->SetParLimits(5,0,100*dphi_par[2]);
          dphi->Fit("fphi2","RQ0");
	  fphi2->Draw("same");
          fphi2->GetParameters(dphi_par);
	  meanerr = fphi2->GetParError(1);
	  sigmaerr = fphi2->GetParError(2);

          fout << "pc3dphimean[" << iarm << "][" << ich << "][" << ivz << "][" << ipt << "]=" << dphi_par[1] << ";" << endl;
          fout << "pc3dphisigma[" << iarm << "][" << ich << "][" << ivz << "][" << ipt << "]=" <<  dphi_par[2] << ";" << endl;
          fout << "pc3dphimeanerr[" << iarm << "][" << ich << "][" << ivz << "][" << ipt << "]=" << meanerr << ";" << endl;
          fout << "pc3dphisigmaerr[" << iarm << "][" << ich << "][" << ivz << "][" << ipt << "]=" <<  sigmaerr << ";" << endl;
          if(ipt == 2){
            fout2 << "PC3_dphifirst_mean[" << iarm << "][" << ich << "][" << ivz << "][" << ipt << "] = " << dphi_par[1] << ";" << endl;
            fout2 << "PC3_dphifirst_sigma[" << iarm << "][" << ich << "][" << ivz << "][" << ipt << "] = " <<  dphi_par[2] << ";" << endl;
          }

          phi_gaus1->SetParameters(dphi_par[0],dphi_par[1],dphi_par[2]);
          phi_gaus1->SetLineColor(1);
          phi_gaus1->Draw("SAME");
          phi_gaus2->SetParameters(dphi_par[3],dphi_par[4],dphi_par[5]);
          phi_gaus2->SetLineColor(6);
          phi_gaus2->Draw("SAME");

          //c->Print(Form("matching/pc3dphi_%d_%d_%d_%d.png",iarm,ich,ipt,ivz));
          delete c;



          Float_t Xbins = dz->GetNbinsX();
          Float_t Xmin = dz->GetXaxis()->GetXmin();
          Float_t Xmax = dz->GetXaxis()->GetXmax();
          mean = (dz->GetMaximumBin() * (Xmax-Xmin))/Xbins + Xmin;
          max = dz->GetMaximum();
          fz1->SetRange(mean-4,mean+4);
	  fz1->SetParameters(max,mean);

          TCanvas *c = new TCanvas("c","c",500,500);
	  //dz->Scale(1./dz->Integral());
	  dz->SetTitle("dz matching");
	  dz->GetXaxis()->SetTitle("dz");
	  dz->GetYaxis()->SetTitle("# tracks");
	  dz->Rebin(5);
          dz->SetMarkerSize(1);
	  dz->SetMarkerStyle(kFullCircle);
	  dz->Draw("P");
          dz->Fit("fz1","RQ0");
          double dz_par[6];
          fz1->GetParameters(dz_par);
          dz_par[3] = 0.1*dz_par[0];
          dz_par[4] = dz_par[1];
          dz_par[5] = 7*dz_par[2];
          fz2->SetParameters(dz_par);
	  fz2->SetParLimits(3,0,5*dz_par[3]);
	  fz2->SetParLimits(4,-1,1);
	  fz2->SetParLimits(5,2*dz_par[2],8*dz_par[2]);
          dz->Fit("fz2","RQ0");
	  fz2->Draw("same");
          fz2->GetParameters(dz_par);
	  meanerr = fz2->GetParError(1);
	  sigmaerr = fz2->GetParError(2);


          fout << "pc3dzmean[" << iarm << "][" << ich << "][" << ivz << "][" << ipt << "]=" << dz_par[1] << ";" << endl;
          fout << "pc3dzsigma[" << iarm << "][" << ich << "][" << ivz << "][" << ipt << "]=" <<  dz_par[2] << ";" << endl;
          fout << "pc3dzmeanerr[" << iarm << "][" << ich << "][" << ivz << "][" << ipt << "]=" << meanerr << ";" << endl;
          fout << "pc3dzsigmaerr[" << iarm << "][" << ich << "][" << ivz << "][" << ipt << "]=" <<  sigmaerr << ";" << endl;
          if(ipt == 2){
            fout2 << "PC3_dzfirst_mean[" << iarm << "][" << ich << "][" << ivz << "][" << ipt << "] = " << dz_par[1] << ";" << endl;
            fout2 << "PC3_dzfirst_sigma[" << iarm << "][" << ich << "][" << ivz << "][" << ipt << "] = " <<  dz_par[2] << ";" << endl;
          }

          z_gaus1->SetParameters(dz_par[0],dz_par[1],dz_par[2]);
          z_gaus1->SetLineColor(1);
          z_gaus1->Draw("SAME");
          z_gaus2->SetParameters(dz_par[3],dz_par[4],dz_par[5]);
          z_gaus2->SetLineColor(6);
          z_gaus2->Draw("SAME");

          //c->Print(Form("matching/pc3dz_%d_%d_%d_%d.png",iarm,ich,ipt,ivz));
          delete c;
	}
      }
    }
  }
  fout << "}" << endl;

}


