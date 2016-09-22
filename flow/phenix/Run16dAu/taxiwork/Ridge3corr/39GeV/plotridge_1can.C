void plotridge_1can(){
//TString dire = "north";
TString dire = "south";
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetStripDecimals(0);
  gStyle->SetErrorX(0);
  float const PI = acos(-1.0);
const int ncent = 8;
const int npt = 10;
  TFile *f=TFile::Open("merged_AnapAumbcentral.root");
TH2F* kforebbcw[ncent][npt];
TH2F* kbackbbcw[ncent][npt];
TH2F* kbackbbcw2[ncent][npt];
TH2F* kforebbcw_ptIn[ncent];
TH2F* kbackbbcw_ptIn[ncent];
TH2F* kbackbbcw2_ptIn[ncent];
double ptbin[npt+1] = {0.2,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0};
double centbin[ncent] = {224,170,146,120,89,65,0};
/*double ptmin = 1.0; double ptmax = 3.0;//has to be boundary
for(int ipt=0; ipt<npt; ipt++){
if(ptmin >= ptbin[ipt] && ptmin < ptbin[ipt+1]){int xptmin = ipt; continue;}
if(ptmax >= ptbin[ipt] && ptmax < ptbin[ipt+1]){int xptmax = ipt; continue;}
}*/
for(int ipt_a=0;ipt_a<1;ipt_a++){
int xptmin = ipt_a*4+2;
int xptmax = (ipt_a+1)*4+2;
double ptmin = ptbin[xptmin];
double ptmax = ptbin[xptmax];

for(int icent=0; icent<ncent; icent++){
for(int ipt=0; ipt<npt; ipt++){
kforebbcw[icent][ipt] = (TH2F*)f->Get(Form("kfore%setabbc_%d_%d",dire.Data(),icent,ipt));
kbackbbcw[icent][ipt] = (TH2F*)f->Get(Form("kback%setabbc_%d_%d",dire.Data(),icent,ipt));
kbackbbcw2[icent][ipt] = (TH2F*)f->Get(Form("kback%setabbc2_%d_%d",dire.Data(),icent,ipt));
}
kforebbcw_ptIn[icent] = (TH2F*)kforebbcw[icent][xptmin]->Clone();
kbackbbcw_ptIn[icent] = (TH2F*)kbackbbcw[icent][xptmin]->Clone();
kbackbbcw2_ptIn[icent] = (TH2F*)kbackbbcw2[icent][xptmin]->Clone();
for(int ipt=xptmin+1; ipt<xptmax; ipt++){
kforebbcw_ptIn[icent]->Add(kforebbcw[icent][ipt]);
kbackbbcw_ptIn[icent]->Add(kbackbbcw[icent][ipt]);
kbackbbcw2_ptIn[icent]->Add(kbackbbcw2[icent][ipt]);
}
//kforebbcw_ptIn[icent]->Rebin(2);
//kbackbbcw_ptIn[icent]->Rebin(2);
//kbackbbcw2_ptIn[icent]->Rebin(2);
}
ofstream fout("out.dat");

TH2F* hpp[ncent];
TH2F* hbackpp[ncent];
  c1=new TCanvas("c1","c1", 600,600);
  c1->SetFillColor(10);
  c1->SetFillColor(10);
  c1->SetFillColor(10);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
for(int icent=0; icent<1; icent++){
/*	if(icent==5){
for(int icent_add=6;icent_add<ncent;icent_add++){
kforebbcw_ptIn[icent]->Add(kforebbcw_ptIn[icent_add]);
kbackbbcw_ptIn[icent]->Add(kbackbbcw_ptIn[icent_add]);
kbackbbcw2_ptIn[icent]->Add(kbackbbcw2_ptIn[icent_add]);
}
centbin[icent+1] = 1.;
}*/
  hpp[icent] = (TH2F*)kforebbcw_ptIn[icent]->Clone();
  hbackpp[icent] = (TH2F*)kbackbbcw2_ptIn[icent]->Clone();
  float nbackpp = hbackpp[icent]->Integral()/2.0/PI;
  float nforepp = hpp[icent]->Integral()/2.0/PI;
  //float ntrig0 = ptforedis_0->Integral(11,30);
   for(int i=0; i<100*40; i++){
     float pp_cont = 1.0*hpp[icent]->GetBinContent(i+1);
     //float pp0_err = 1.0*hpp[icent]->GetBinError(i+1);
     float weight2 = 0;//sqrt(1.0*hbackbbcw_ptIn[icent]->GetBinContent(i+1));

     float backpp_cont = 1.0*hbackpp[icent]->GetBinContent(i+1);
    if(backpp_cont==0){
	 float con = 0;
	float err = 0;
    }
    else{
     float con = pp_cont/backpp_cont*nbackpp/nforepp;
     float err = weight2/backpp_cont*nbackpp/nforepp;
	}
     hpp[icent]->SetBinContent(i+1, con);
     hpp[icent]->SetBinError(i+1, err);
   }

  c1->cd();
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetFrameFillColor(0);
  c1->SetFrameBorderMode(0);
  c1->SetFrameBorderMode(0);
  gPad->SetLeftMargin(0.130);
  gPad->SetRightMargin(0.02);
  gPad->SetTopMargin(0.01);
  gPad->SetTicks(1,1);
  
  float ymax = 1.147;//hpp[icent]->GetMaximum()*1.1;
  float ymin = 0.868;//hpp[icent]->GetMinimum()*0.9;

  hpp[icent]->SetMinimum(ymin);
  hpp[icent]->SetMaximum(ymax);
  hpp[icent]->GetXaxis()->SetRangeUser(2.8,3.6);
  hpp[icent]->GetYaxis()->SetRangeUser(-1,4);

  hpp[icent]->SetMarkerStyle(20);
  hpp[icent]->SetMarkerSize(1.3);
  
  hpp[icent]->SetMarkerColor(4);
  hpp[icent]->GetYaxis()->SetLabelSize(0.08);
  hpp[icent]->GetXaxis()->SetLabelSize(0.08);
  //hpp[icent]->GetYaxis()->SetTitleSize(0.6);
  //hpp[icent]->GetXaxis()->SetTitleSize(0.0);
  hpp[icent]->Draw("surf1");
  
TLatex *t=new TLatex(-1.0,0.88*(ymax-ymin)+ymin, Form("%d)",icent+1));
  t->SetTextSize(0.08);
  t->Draw();

  TMarker *m0 = new TMarker(-0.3, 0.91*(ymax-ymin)+ymin, 20);
  m0->SetMarkerSize(1.3);
  m0->SetMarkerColor(4);
  m0->Draw();

  TLatex *t=new TLatex(-0.0,0.88*(ymax-ymin)+ymin, Form("p+p %.f - %.f (vtx z)",centbin[icent],centbin[icent+1]));
  t->SetTextSize(0.08);
  t->Draw();

  TLatex *t=new TLatex(-1.2,0.78*(ymax-ymin)+ymin, Form("%.1f<p_{T,trig}<%.1f GeV/c",ptmin,ptmax));
  t->SetTextSize(0.08);
  if(icent==0)t->Draw();

  TLatex *t=new TLatex(-1.0,0.66*(ymax-ymin)+ymin, "|#eta_{trig}|<0.35");
  t->SetTextSize(0.08);
  if(icent==0)t->Draw();

if(dire=="north")
  TLatex *t=new TLatex(-1.2,0.14*(ymax-ymin)+ymin, "1.5<|#eta_{asso}|<3.0");
else
  TLatex *t=new TLatex(-1.2,0.14*(ymax-ymin)+ymin, "1.0<|#eta_{asso}|<3.0");
  t->SetTextSize(0.08);
  t->SetTextColor(2);
  if(icent==0)t->Draw();

}
  c1->cd();

  c1->Print(Form("fig/ridge2D_1can_%s_%.1f_%.1f.png",dire.Data(),ptmin,ptmax));
}

/*  
  cout<<"*************** v2 ***********"<<endl;
  float v2_0 = funvn0->GetParameter(1)/(zym_pp0 + funvn0->GetParameter(0));
  float v2_1 = funvn1->GetParameter(1)/(zym_pp1 + funvn1->GetParameter(0));
  float v2_2 = funvn2->GetParameter(1)/(zym_pp2 + funvn2->GetParameter(0));
  float v2_3 = funvn3->GetParameter(1)/(zym_pp3 + funvn3->GetParameter(0));
  

  cout<<v2_0<<" "<<v2_1<<" "<<v2_2<<" "<<v2_3<<endl;
 
  
  cout<<funvn0->GetParameter(0)*funvn0->GetParameter(1)<<" "
      <<funvn1->GetParameter(0)*funvn1->GetParameter(1)<<" "
      <<funvn2->GetParameter(0)*funvn2->GetParameter(1)<<" "
      <<funvn3->GetParameter(0)*funvn3->GetParameter(1)<<endl;
  

  cout<<"*************** v2 ***********"<<endl;
  cout<<funvn0->GetParameter(2)<<" "<<funvn1->GetParameter(2)<<" "<<funvn2->GetParameter(2)<<" "<<funvn3->GetParameter(2)<<endl;

  cout<<"*************** v3 ***********"<<endl;
  cout<<funvn0->GetParameter(3)<<" "<<funvn1->GetParameter(3)<<" "<<funvn2->GetParameter(3)<<" "<<funvn3->GetParameter(3)<<endl;
*/
}
