#include "/afs/cern.ch/user/q/qixu/CMSSW_6_2_5/src/Centrality/NBD_Glauber_fit/NBD/parameter.h"
#include "par.h"
void Drawsys(){ 
 TCanvas *c1 = new TCanvas("c1","c1",1,1,550,460);
  c1->SetFillColor(10);
  c1->SetFrameFillColor(0);
  c1->SetFrameBorderSize(0);
  c1->SetFrameBorderMode(0);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.15);
  c1->SetTopMargin(0.02);
  c1->SetRightMargin(0.02);
  gStyle->SetOptStat(0);
  c1->SetTicks(-1);
TString Var[nGlau+2]={"standard","Depth04914","Depth06006","dmin00","dmin08","Radius675","Radius649","sNN65mb","sNN75mb",Form("bin-%.f",binshift),Form("bin+%.f",binshift)};

 N=N-1;
 TString str="Ncoll";
 TString dirname="std";
 TH1D* hist = new TH1D("","",N,0,N);
 hist->GetXaxis()->SetNdivisions(502);
if(method==0)
 hist->SetXTitle("Centrality");
else 
 hist->SetXTitle("Multiplicity");
 hist->SetYTitle(Form("#frac{Different Glauber Parameter <%s>}{Default <%s>}",str.Data(),str.Data()));
 hist->SetMinimum(0.80);
 hist->SetMaximum(1.15);
 hist->GetXaxis()->CenterTitle(0);
 hist->GetYaxis()->CenterTitle(1);
 hist->GetYaxis()->SetTitleOffset(1.5);
 hist->GetXaxis()->SetTitleOffset(1.5);
 hist->GetXaxis()->SetTitleSize(0.04);
 hist->GetYaxis()->SetTitleSize(0.04);
 hist->GetXaxis()->SetLabelSize(0.05);
 hist->GetYaxis()->SetLabelSize(0.05);
 hist->GetXaxis()->SetLabelOffset(99);
 hist->Draw();

	double xbin[100];
	TFile *f = TFile::Open("G.root");
	TGraph *graph[nGlau+2];
        TVectorD* NcollAver0 = (TVectorD*)f->Get(Form("%s/G0/NcollAver",dirname.Data()));
	TVectorD *centbin = (TVectorD*)f->Get(Form("%s/G0/centbin",dirname.Data()));
        TVectorD *kpoint = (TVectorD*)f->Get(Form("%s/G0/kpoint",dirname.Data()));

	std::vector<double> ratio;
	TString name;
	for(int iGlau=0;iGlau<nGlau+2;iGlau++){
   		if(iGlau==0)
                	name = "G0";
                else if(iGlau<nGlau && iGlau>0)
                	name = Form("Glau_%d",iGlau);
               	else
                	name = Form("bin_%d",iGlau-nGlau+1);

               	TVectorD* NcollAver = (TVectorD*)f->Get(Form("%s/%s/NcollAver",dirname.Data(),name.Data()));
		ratio.clear();
		for(int i=0;i<N;i++){
        	xbin[i]=0.29+i;
		ratio.push_back((*NcollAver)[i]/(*NcollAver0)[i]);
		}
        	double* aratio  = &ratio[0];
		graph[iGlau] = new TGraph(N,xbin,aratio);
		graph[iGlau]->SetTitle("");
		graph[iGlau]->SetMarkerStyle(marker[iGlau]);
		graph[iGlau]->SetMarkerColor(color[iGlau]);
		graph[iGlau]->SetLineColor(1);
		graph[iGlau]->SetLineWidth(2);
		graph[iGlau]->SetMarkerSize(1.2);
		graph[iGlau]->Draw("Psameez");
	}

//char* label[9]={"60-90%","40-60%","30-40%","20-30%","10-20%","5-10%","1-5%","0-1%","0-90%"};
//char* label[N]={"59.7-100%", "27.4-59.7%","16.6-27.4%","9.41-16.6%", "2.46-9.41%", "0-2.46%"};
TString label[N];
for(int i=0;i<N;i++)
        if(method==0)label[i] = Form("%.2f-%.2f%%",(*centbin)[i]*100,(*centbin)[i+1]*100);
        else label[i] = Form("%.2f-%.2f",(*kpoint)[i],(*kpoint)[i+1]);

    TLatex *tex1= new TLatex(0.6,1.12,"CMS Preliminary pPb \ #sqrt{s_{NN}} = 5.02 TeV");
    tex1->SetTextColor(1);
    tex1->SetTextFont(42);
    tex1->SetTextSize(0.05);
    tex1->Draw();

TLine *l = new TLine(hist->GetXaxis()->GetXmin(),1.,hist->GetXaxis()->GetXmax(),1.);
l->SetLineStyle(2);
l->Draw();
double y = hist->GetMinimum();
   TText t;
   t.SetTextAngle(45);
   t.SetTextSize(0.03);
   t.SetTextAlign(33);
   for (int i=0;i<N;i++) {
      double x = hist->GetXaxis()->GetBinCenter(i+1);
      t.DrawText(x,y,label[i]);
   }
TLegend *leg0 = new TLegend(0.42,0.15,0.67,0.40);
TLegend *leg1 = new TLegend(0.75,0.15,0.83,0.40);
TLegend *leg2 = new TLegend(0.20,0.15,0.35,0.40);
    leg0->SetFillColor(10);
    leg1->SetFillColor(10);
    leg2->SetFillColor(10);
    leg0->SetBorderSize(0.035);
    leg1->SetBorderSize(0.035);
    leg2->SetBorderSize(0.035);
    leg0->SetTextFont(42);
    leg1->SetTextFont(42);
    leg2->SetTextFont(42);
    leg0->SetTextSize(0.035);
    leg1->SetTextSize(0.035);
    leg2->SetTextSize(0.035);
  for(int iGlau=0; iGlau<nGlau+2; iGlau++){
	if(iGlau==0 || iGlau>=nGlau)
  		leg2->AddEntry(graph[iGlau],Var[iGlau].Data(),"p");
 	else{
		if(iGlau%2==1)
  			leg0->AddEntry(graph[iGlau],Var[iGlau].Data(),"p");
		else
  			leg1->AddEntry(graph[iGlau],Var[iGlau].Data(),"p");
	}
   }

   leg0->Draw();	
   leg1->Draw();	
   leg2->Draw();	
c1->SaveAs(Form("%ssys.png",dirname.Data()));

}
