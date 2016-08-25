#include "/afs/cern.ch/user/q/qixu/CMSSW_6_2_5/src/Centrality/NBD_Glauber_fit/NBD/parameter.h"
#include <vector>
#include "par.h"

std::vector<double> calcsys(int sth, int Gth){
	if (!TClass::GetDict("NBD")) {
                gROOT->ProcessLine(".L /afs/cern.ch/user/q/qixu/CMSSW_6_2_5/src/Centrality/NBD_Glauber_fit/NBD/NBDclass.C+");
        }
	TFile *f0 = TFile::Open(Form("G%d%d.root",sth,Gth));
                TList *tl = f0->GetListOfKeys();
                key = (TKey*) tl->At(0);
                NBD* temp = (NBD*)key->ReadObj();
		return temp->NcollAver;
}

void Drawsys(){ 
/*	   if (!TClass::GetDict("NBD")) {
                gROOT->ProcessLine(".L /afs/cern.ch/user/q/qixu/CMSSW_6_2_5/src/Centrality/NBD_Glauber_fit/NBD/NBDclass.C++");
        }*/
/* TCanvas *c1 = new TCanvas("c1","c1",1,1,550,460);
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
*/
TString GlauVar[nGlau]={"standard","Depth04914","Depth06006","dmin00","dmin08","Radius675","Radius649","sNN65mb","sNN75mb"};

 N=N-1;
 TString str="Ncoll";
 TString dirname="std";
 int sth=0;
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
// hist->Draw();

	double xbin[100];
	ofstream fout("sys.txt");
	TGraph *graph[nGlau];
/*	TFile *f = TFile::Open("G.root");
	dir = f->GetDirectory(dirname);
        def = new TList();
        tl = dir->GetListOfKeys();
        for(int iGlau=0; iGlau<nGlau; iGlau++){
                key = ((TKey*) tl->At(iGlau));
                NBD* temp = (NBD*)key->ReadObj();
                def->Add(temp);
        }
*/
        for(int i=0;i<N;i++)
                xbin[i]=0.29+i;
		vector<double> NcollAver0 = calcsys(sth,0);
        for(int iGlau=0; iGlau<nGlau; iGlau++){
		/*fout<<dirname<<"\t"<<iGlau<<"\t"<<n->Glaubername<<endl;
		fout<<"mu="<<n->mubest<<"\t"<<"k="<<n->kbest<<endl;
		fout<<"chis="<<n->chis<"\t"<<"Ndf="<<n->Ndf<<endl;
		if(iGlau==0) fout<<"kpoint"<<"\t"<<"centbin"<<"\t"<<"Ncoll"<<"\t"<<"Npart"<<"\t"<<"B"<<endl;*/
//		for(int j=0;j<N;j++)
//		fout<<n->kpoint[j]<<"-"<<n->kpoint[j+1]<<"\t"<<n->centbin[j]<<"-"<<centbin[j+1]<<"\t"<<n->NcollAver[j]<<"\t"<<n->NpartAver[j]<<"\t"<<n->BAver[j]<<endl;
		vector<double> NcollAver= calcsys(sth,iGlau);
		vector<double> ratio;
		for(int j=0;j<N;j++){
		fout<<NcollAver[j]<<"\t";
		ratio.push_back(NcollAver[j]/NcollAver0[j]);
		}
		fout<<endl;
        	double* aratio  = &ratio[0];
		graph[iGlau] = new TGraph(N,xbin,aratio);
		graph[iGlau]->SetTitle("");
		graph[iGlau]->SetMarkerStyle(marker[iGlau]);
		graph[iGlau]->SetMarkerColor(color[iGlau]);
		graph[iGlau]->SetLineColor(1);
		graph[iGlau]->SetLineWidth(2);
		graph[iGlau]->SetMarkerSize(1.2);
//		graph[iGlau]->Draw("Psameez");
	}

//char* label[9]={"60-90%","40-60%","30-40%","20-30%","10-20%","5-10%","1-5%","0-1%","0-90%"};
//char* label[N]={"59.7-100%", "27.4-59.7%","16.6-27.4%","9.41-16.6%", "2.46-9.41%", "0-2.46%"};
TString label[N];
for(int i=0;i<N;i++){
	if(method==0)label[i] = Form("%.2f-%.2f%%",centbin[N-i-1]*100,centbin[N-i]*100);
	else label[i] = Form("%.2f-%.2f",centbin[N-i-1],centbin[N-i]);
}
    TLatex *tex1= new TLatex(0.6,1.12,"CMS Preliminary pPb \ #sqrt{s_{NN}} = 5.02 TeV");
    tex1->SetTextColor(1);
    tex1->SetTextFont(42);
    tex1->SetTextSize(0.05);
 //   tex1->Draw();

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
  //    t.DrawText(x,y,label[i]);
   }
TLegend *leg0 = new TLegend(0.32,0.15,0.50,0.40);
TLegend *leg1 = new TLegend(0.60,0.15,0.78,0.40);
    leg0->SetFillColor(10);
    leg1->SetFillColor(10);
    leg0->SetBorderSize(0.035);
    leg1->SetBorderSize(0.035);
    leg0->SetTextFont(42);
    leg1->SetTextFont(42);
    leg0->SetTextSize(0.035);
    leg1->SetTextSize(0.035);
  for(int iGlau=0; iGlau<5; iGlau++)
  leg0->AddEntry(graph[iGlau],GlauVar[iGlau].Data(),"p");
  for(int iGlau=5; iGlau<nGlau; iGlau++)
  leg1->AddEntry(graph[iGlau],GlauVar[iGlau].Data(),"p");
//   leg0->Draw();	
//   leg1->Draw();	
//c1->SaveAs(Form("%ssys.png",dirname.Data()));

}
