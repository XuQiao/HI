void plotspflowvsNtrk(){
gStyle->SetOptFit(kFALSE);
gStyle->SetOptStat(kFALSE);
TCanvas *c1 = new TCanvas;
TCanvas *c2 = new TCanvas;
TLegend *leg = new TLegend(0.2,0.1,0.4,0.4);
leg->SetTextSize(0.05);
leg->SetBorderSize(0);
leg->SetFillColor(0);
const int ntotbin=12; 
//const int ntotbin=6; 
const int multbin[ntotbin]={30,40,50,65,80,100,120,150,185,220,260,300};
//const int multbin[ntotbin]={50,65,100,120,150,185};
int xbin=0;
TString dir;
double V2_Prod[ntotbin],V2err_Prod[ntotbin],Mult[ntotbin],V2sp[ntotbin],V2spError[ntotbin];

TF1 *PtDistr  = new TF1("PtDistr","0.03*(exp (-(x/0.594540))+0.00499506*exp (-(x/1.89391)))", 0.3,6.0);       //Real Data
TF1 *V2vsPt = new TF1("V2vsPt","((x/3.31699)^2.35142/(1+(x/3.49188)^3.54429))*(.00005+(1/x)^1.50600)",0.3,6.0);
TF1 *V2 = new TF1("V2","0.03*(exp (-(x/0.594540))+0.00499506*exp (-(x/1.89391)))*((x/3.31699)^2.35142/(1+(x/3.49188)^3.54429))*(.00005+(1/x)^1.50600)",0.3,6.0);
double InV2 = V2->Integral(0.3,6.0)/PtDistr->Integral(0.3,6.0);  

for(int i=0;i<ntotbin;i++){
	dir=Form("pPbDataV205m%d",multbin[i]);
	TFile *mergedV_Prod = TFile::Open(Form("%s/mergedV_Prod.root",dir.Data()));
  //      if(multbin[i]==185)
	//TFile *mergedV_Prod = TFile::Open(Form("%s/mergedV_Sum.root",dir.Data()));
	TVectorD *vecMult = (TVectorD*)mergedV_Prod->Get("totmultall");
	TVectorD *vecNevent = (TVectorD*)mergedV_Prod->Get("Nevent");
	TVectorD *vecV2_Prod=(TVectorD*)mergedV_Prod->Get(Form("Vmean",xbin));
	TVectorD *vecV2err_Prod=(TVectorD*)mergedV_Prod->Get(Form("deltaVmean",xbin));
	Mult[i]=(*vecMult)[xbin]/(*vecNevent)[xbin];
	V2_Prod[i]=(*vecV2_Prod)[xbin];
	V2err_Prod[i]=(*vecV2err_Prod)[xbin];
	V2sp[i]=InV2/V2_Prod[i];
	V2spError[i]=InV2/V2_Prod[i]/V2_Prod[i]*(*vecV2err_Prod)[xbin];
	//cout<<V2sp[i]<<"pm"<<V2spError[i]<<"\t";
	cout<<V2sp[i]<<",";
	}
//    V2_Prod[1]*= 1.25;        
//    V2_Prod[2] = 0.0411;
cout<<endl<<InV2<<endl;
c1->cd();
TH1D* hFrame = new TH1D("","",400,0,400);
hFrame->GetXaxis()->SetTitle("Ntrkoffline");
hFrame->GetYaxis()->SetTitle("Real V2/Calc V2");
hFrame->GetXaxis()->SetRangeUser(25,320);
hFrame->GetYaxis()->SetRangeUser(0.5,1.1);
hFrame->Draw();
TGraphErrors *grProd=new TGraphErrors(ntotbin,Mult,V2_Prod,0,V2err_Prod);
TGraphErrors *grsp = new TGraphErrors(ntotbin,Mult,V2sp,0,V2spError);
TLine *l = new TLine(30,1,320,1);
l->SetLineStyle(2);
grsp->SetMarkerSize(1.6);
grsp->SetMarkerStyle(20);
grProd->SetMarkerColor(1);
grProd->SetMarkerSize(1.6);
grProd->SetLineColor(1);
grProd->SetMarkerStyle(20);
grsp->Draw("Psame");
l->Draw("same");
c1->SaveAs("spflow.gif");
c2->cd();
hFrame->GetYaxis()->SetTitle("Calc V2");
hFrame->GetXaxis()->SetRangeUser(25,320);
hFrame->GetYaxis()->SetRangeUser(0.042,0.1);
hFrame->Draw();
grProd->Draw("Psame");
TLine *l = new TLine(25,InV2,320,InV2);
l->SetLineStyle(2);
l->Draw("same");
c2->SaveAs("V205_Ntrk.gif");
}

