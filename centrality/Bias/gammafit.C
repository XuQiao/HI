double function(double *v, double *par){
	TF1 *gamma = new TF1("gamma","TMath::GammaDist(x,[0],0,[1])",0,100);
//	TF1 *gaus = new TF1("gaus","[0]*TMath::Gaus(x,[1]-[2],[3])",0,100);
	gamma->SetParameters(par[0],par[1]);
//	gaus->SetParameters(par[2],v[0],par[3],par[4]);
//	TF1 *conv = new TF1("conv","gamma*gaus",0,100);
//	return conv->Integral(0,v[0]);
	return gamma->Eval(v[0]);
}

void gammafit(){
	gStyle->SetOptFit(kTRUE);
	double xmin=5;
	double xmax=60;
	ofstream fout("fitresultMC.dat");
	TFile *f = TFile::Open("PythiaMBhistos.root");
	TH1D* hHFp4 = (TH1D*)f->Get("hHFplus4");
	TH1D* hHFp = (TH1D*)f->Get("hHFplus");
	TH1D* hHF = (TH1D*)f->Get("hHF");
	TH1D* hHF4 = (TH1D*)f->Get("hHF4");
	hHFp4->Scale(1.0/hHFp4->Integral("width"));
	hHFp->Scale(1.0/hHFp->Integral("width"));
	hHF->Scale(1.0/hHF->Integral("width"));
	hHF4->Scale(1.0/hHF4->Integral("width"));
	hHFp4->GetXaxis()->SetRangeUser(0,50);
	hHFp->GetXaxis()->SetRangeUser(0,50);
	hHF->GetXaxis()->SetRangeUser(0,100);
	TCanvas *c1= new TCanvas;
	c1->SetLogy();
	TF1 *g = new TF1("g",function,0,100,2);
	TF1 *NBD_fun = new TF1("NBD_fun","[0]*TMath::Gamma(x+[1])/(TMath::Gamma(x+1)*TMath::Gamma([1]))*TMath::Power([2]/[1],x)/TMath::Power([2]/[1]+1,x+[1])",0,100);
	NBD_fun->FixParameter(0,1);
	NBD_fun->SetParameter(1,0.5);
	NBD_fun->SetParameter(2,1.0);
	g->SetParameters(0.8,2);
	hHF4->Fit("g","","",xmin,xmax);
//	hHFp->Fit("g","","",xmin,xmax);
//	hHFp4->Fit("g","","",xmin,xmax);
//	g->Draw("same");
	TLatex *tl = new TLatex(0.2,0.6,Form("fit range %.1f to %.1f",xmin,xmax));
	tl->Draw("same");
	c1->Print("hHFp3FitGamma.png");
	fout<<g->GetParameter(0)<<"\t"<<g->GetParError(0)<<endl;
	fout<<g->GetParameter(1)<<"\t"<<g->GetParError(1)<<endl;
//	cout<<g->GetParameter(2)<<"\t"<<g->GetParError(2)<<endl;
//	cout<<g->GetParameter(3)<<"\t"<<g->GetParError(3)<<endl;
//	cout<<g->GetParameter(4)<<"\t"<<g->GetParError(4)<<endl;
	fout<<g->GetNDF()<<"\t"<<g->GetChisquare()<<endl;
}
