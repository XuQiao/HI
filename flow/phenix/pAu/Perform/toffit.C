#include "/phenix/u/xuq/util/SimplifyLife.C"
void toffit(){
	gStyle->SetErrorX(0);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	TFile *f = TFile::Open("mergedFull.root");
	TString charge;
	ofstream fout("Run15tofdphidzcalib.dat");
	ofstream fouttmp("Run15tofdphidzcalibtmp.dat");
	for(int iarm=0; iarm<2; iarm++){
	for(int ich=0; ich<2; ich++){
	if(iarm==0)
	TH2F *tofdphidz = (TH2F*)f->Get(Form("tofdphidz_%d",ich));
	else
	TH2F *tofdphidz = (TH2F*)f->Get(Form("tofwdphidz_%d",ich));
	if(tofdphidz->Integral()==0) continue;
	TH1F* tofdphi = (TH1F*)tofdphidz->ProjectionX(Form("tofdphi_%d_%d",iarm,ich),0,-1);
	TH1F* tofdz = (TH1F*)tofdphidz->ProjectionY(Form("tofdz_%d_%d",iarm,ich),0,-1);
	//tofdphi->Rebin(4);
	tofdz->Rebin(5);
	tofdphi->Scale(1./tofdphi->Integral());
	tofdz->Scale(1./tofdz->Integral());
	TF1 *fphi = new TF1("fphi","gaus(0)+gaus(3)",-0.1,0.1);
	fphi->SetNpx(10000);
	//TF1 *fphi = new TF1("fphi","[0]*([1]*exp(-(x-[2])**2/[3]/[3])+(1-[1])*exp(-(x-[4])**2/[5]/[5]))",-0.1,0.1);
	fphi->SetParameters(1e-1,0,1e-3,1e-3,0,1e-2);
	//fphi->SetParameters(4.79253e-02,9.25529e-01,-8.56488e-05,-7.46701e-03,-5.37828e-04,5.99178e-0);
	SetStyle(*tofdphi,1.2,1,20,0,0);
	tofdphi->Fit("fphi","RQ");
	double par[6];
	fphi->GetParameters(par);
	TF1 *fphi1 = new TF1("fphi1","gaus",-0.1,0.1);
	fphi1->SetNpx(10000);
	fphi1->SetParameters(par);
	TF1 *fphi2 = new TF1("fphi2","gaus",-0.1,0.1);
	fphi2->SetNpx(10000);
	fphi2->SetParameters(&par[3]);
	fphi1->SetLineColor(4);
	//fphi1->Draw("same");
	fphi2->SetLineColor(6);
	//fphi2->Draw("same");
	fout<<par[1]<<",";
	if(par[1]>1) cout<<"dphi "<<iarm<<"\t"<<ich<<"\t"<<ipt<<endl;

	TF1 *fz = new TF1("fz","gaus(0)+gaus(3)",-10,10);
	//fz->SetParameters(3.02349e+07,8.32462e-01,2.39187e+00,2.17631e+07,6.35460e-01,8.09821e+00);
	fz->SetParameters(1e-1,0,2,1e-3,0,8);
	SetStyle(*tofdz,1.2,1,20,0,0);
	tofdz->Fit("fz","RQ");
	fz->GetParameters(par);
	TF1 *fz1 = new TF1("fz1","gaus",-10,10);
	fz1->SetNpx(10000);
	fz1->SetParameters(par);
	TF1 *fz2 = new TF1("fz2","gaus",-10,10);
	fz2->SetNpx(10000);
	fz2->SetParameters(&par[3]);
	fz1->SetLineColor(4);
	//fz1->Draw("same");
	fz2->SetLineColor(6);
	//fz2->Draw("same");
	if(par[1]>10) cout<<"dz"<<iarm<<"\t"<<ich<<"\t"<<ipt<<endl;
	fouttmp<<par[1]<<",";
	}
	fout<<endl;
	fouttmp<<endl;
}

}
