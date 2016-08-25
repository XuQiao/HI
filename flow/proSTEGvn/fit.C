void fit(){
	TCanvas *c1 = new TCanvas;
	c1->SetLogy();
	TFile *f = TFile::Open("/home/xuq7/HI/flow/LYZ/v2/tracknormcpt03to6/M150120/histomerged.root");
	TH1D* hpt = (TH1D*)f->Get("D_0/hpt");
	hpt =(TH1D*)hpt->DrawNormalized();
	hpt->SetMarkerSize(0.9);
	hpt->SetMarkerStyle(20);
        hpt->GetXaxis()->SetTitle("p_{T}");
        hpt->GetYaxis()->SetTitle("Probability");
	hpt->Draw("P");
   	//TF1 *PtDistr  = new TF1("PtDistr","[3]*(exp (-(x/[0]))+[1]*exp (-(x/[2])))", 0.3,6);
   	 TF1 *PtDistr  = new TF1("PtDistr","0.03*(exp (-(x/0.594540))+0.00499506*exp (-(x/1.89391)))", 0.3,6.0);
	//TF1 *PtDistr = new TF1("PtDistr","[0]*x*TMath::Power(1+(sqrt(x*x+[1]*[1])-[1]/[2]),-[3])",0.2,10);
	//TF1 *PtDistr = new TF1("PtDistr","[0]*TMath::Power(x,[1])",0.2,10);
	//PtDistr->SetParameters(1.7,0.45,10000,6.5,);
	//PtDistr->SetParameters(0.594,0.005,1.9,1.0);
	//hpt->Fit(PtDistr,"","same",0.3,6);
	//cout<<PtDistr->GetParameter(0)<<","<<PtDistr->GetParameter(1)<<","<<PtDistr->GetParameter(2)<<","<<PtDistr->GetParameter(3)<<endl;
	//cout<<PtDistr->GetChisquare()<<"\t"<<PtDistr->GetNDF()<<endl;
	PtDistr->Draw("same");
	c1->Print("ptfit.pdf");
//	TF1 *V2vsPt = new TF1("V2vsPt","((x/4.3)^1.8/(1+(x/3)^1.8))*(.00005+(1/x)^0.8)",0.2,10);      //Real data V~0.06
	TCanvas *c2 = new TCanvas();
	//TF1 *V2vsPt = new TF1("V2vsPt","((x/[0])^[1]/(1+(x/[2])^[3]))*(.00005+(1/x)^[4])",0.3,6);      //Real data V~0.06
	TF1 *V2vsPt = new TF1("V2vsPt","((x/3.31699)^2.35142/(1+(x/3.49188)^3.54429))*(.00005+(1/x)^1.50600)",0.3,6.0);
	const double pt[9]={0.4,0.7,1.2,1.7,2.2,2.7,3.4,4.4,5.4};
	//const double v2[9]={0.02582,0.04324,0.06644,0.08396,0.09725,0.1032,0.09603,0.07668,0.0605};
	const double v2[9]={0.02682,0.04447,0.06771,0.08586,0.09932,0.09782,0.08714,0.06139,0.05068};
	const int N=9;
	TGraph *gr = new TGraph(N,pt,v2);
        gr->SetTitle("hv2");
        gr->GetXaxis()->SetTitle("p_{T}");
        gr->GetYaxis()->SetTitle("v_{2}");
	gr->SetMarkerSize(1.3);
	gr->SetMarkerStyle(20);
	gr->SetMarkerColor(1);
	gr->Draw("AP");
	//V2vsPt->SetParameters(4.3,1.8,3,1.8,0.8);
	//gr->Fit(V2vsPt,"","same",0.3,6);
	//cout<<V2vsPt->GetParameter(0)<<","<<V2vsPt->GetParameter(1)<<","<<V2vsPt->GetParameter(2)<<","<<V2vsPt->GetParameter(3)<<","<<V2vsPt->GetParameter(4)<<endl;
	V2vsPt->Draw("same");
	c2->Print("v2fit.pdf");

}
