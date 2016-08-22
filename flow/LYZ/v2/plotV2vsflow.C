void plotV2vsflow(){
gStyle->SetOptFit(kFALSE);
gStyle->SetOptStat(kFALSE);
TCanvas *c1 = new TCanvas("c1","c1",600,600);
TLegend *leg = new TLegend(0.2,0.1,0.4,0.4);
leg->SetTextSize(0.05);
leg->SetBorderSize(0);
leg->SetFillColor(0);
const int ninflow = 6;
int xbin=0;
TString dir;
double V2_Prod[ninflow],V2err_Prod[ninflow];
const int mult = 150;

TF1 *PtDistr  = new TF1("PtDistr","0.03*(exp (-(x/0.594540))+0.00499506*exp (-(x/1.89391)))", 0.3,6.0);       //Real Data
TF1 *V2vsPt = new TF1("V2vsPt","((x/3.31699)^2.35142/(1+(x/3.49188)^3.54429))*(.00005+(1/x)^1.50600)",0.3,6.0);
TF1 *V2 = new TF1("V2","0.03*(exp (-(x/0.594540))+0.00499506*exp (-(x/1.89391)))*((x/3.31699)^2.35142/(1+(x/3.49188)^3.54429))*(.00005+(1/x)^1.50600)",0.3,6.0);
double InV2 = V2->Integral(0.3,6.0)/PtDistr->Integral(0.3,6.0);
double inV2[ninflow]={};
for(int i=0;i<ninflow;i++){
        inV2[i] = InV2 * i / 5 ; 
	dir=Form("pPbDataV20%dm%d",i,mult);
	TFile *mergedV_Prod = TFile::Open(Form("%s/mergedV_Prod.root",dir.Data()));
        if(i==)
	TFile *mergedV_Prod = TFile::Open(Form("%s/mergedV_Sum.root",dir.Data()));
	TVectorD *vecMult = (TVectorD*)mergedV_Prod->Get("totmultall");
	TVectorD *vecNevent = (TVectorD*)mergedV_Prod->Get("Nevent");
	TVectorD *vecV2_Prod=(TVectorD*)mergedV_Prod->Get(Form("Vmean",xbin));
	TVectorD *vecV2err_Prod=(TVectorD*)mergedV_Prod->Get(Form("deltaVmean",xbin));
	V2_Prod[i]=(*vecV2_Prod)[xbin];
	V2err_Prod[i]=(*vecV2err_Prod)[xbin];
        cout<<inV2[i]<<endl;
//	V2sp[i]=InV2/V2_Prod[i];
//	V2spError[i]=InV2/V2_Prod[i]/V2_Prod[i]*(*vecV2err_Prod)[xbin];
//	cout<<V2sp[i]<<"pm"<<V2spError[i]<<"\t";
	}
        V2_Prod[1]=0;
cout<<endl<<InV2<<endl;
TH1D* hFrame = new TH1D("","",2000,-1,1);
hFrame->GetXaxis()->SetRangeUser(0., 0.08);
hFrame->GetYaxis()->SetRangeUser(0., 0.08);
hFrame->Draw();
TGraphErrors *grProd=new TGraphErrors(ninflow,inV2,V2_Prod,0,V2err_Prod);
grProd->SetMarkerColor(1);
grProd->SetMarkerSize(1.6);
grProd->SetLineColor(1);
grProd->SetMarkerStyle(20);
hFrame->GetXaxis()->SetTitle("Input V2");
hFrame->GetYaxis()->SetTitle("Calc V2");
hFrame->Draw();
grProd->Draw("Psame");
TLine *l = new TLine(0,0,0.06,0.06);
l->SetLineStyle(2);
l->Draw("same");
c1->SaveAs("V2vsflow.gif");
}

