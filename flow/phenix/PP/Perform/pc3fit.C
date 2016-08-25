//#include "/phenix/u/xuq/util/SimplifyLife.C"
void pc3fit(){
	gStyle->SetErrorX(0);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	TFile *f = TFile::Open("merged_AnapAlmb.root");
	TString charge;
	ofstream fout("Run15pc3dphidzcalib.dat");
	fout<<"double dphishift[2][2][50];"<<endl;
	fout<<"double dzshift[2][2][50];"<<endl;
	fout<<"double dphisigma[2][2][50];"<<endl;
	fout<<"double dzsigma[2][2][50];"<<endl;
	TCanvas *c1 = new TCanvas("c1","c1",500,500);
	TCanvas *c2 = new TCanvas("c2","c2",500,500);
	for(int ipt=0; ipt<50; ipt++){
	for(int iarm=0; iarm<2; iarm++){
	for(int ich=0; ich<2; ich++){
	if(ich==0) charge = "pos";
	else if(ich==1) charge = "neg";
	//if(iarm!=0 || ich!=1 || ipt !=1) continue;
	TH2F *pc3dphidz = (TH2F*)f->Get(Form("pc3dphidz_arm%d_%s_%d",iarm,charge.Data(),ipt));
	if(pc3dphidz->Integral()==0){
		fout<<"dphishift["<<iarm<<"]["<<ich<<"]["<<ipt<<"] = "<<-0.0<<";  ";
		fout<<"dzshift["<<iarm<<"]["<<ich<<"]["<<ipt<<"] = "<<-0.0<<";  ";
		fout<<"dphisigma["<<iarm<<"]["<<ich<<"]["<<ipt<<"] = "<<-0.0<<";  ";
		fout<<"dzsigma["<<iarm<<"]["<<ich<<"]["<<ipt<<"] = "<<-0.0<<";  ";
		continue;
	}
	TH1F* pc3dphi = (TH1F*)pc3dphidz->ProjectionX(Form("pc3dphi_%d_%d_%d",iarm,ich,ipt),0,-1);
	TH1F* pc3dz = (TH1F*)pc3dphidz->ProjectionY(Form("pc3dz_%d_%d_%d",iarm,ich,ipt),0,-1);
	//pc3dphi->Rebin(4);
	pc3dz->Rebin(5);
	pc3dphi->Scale(1./pc3dphi->Integral());
	pc3dz->Scale(1./pc3dz->Integral());
	TF1 *fphi = new TF1("fphi","gaus(0)+gaus(3)",-0.1,0.1);
	fphi->SetNpx(10000);
	//TF1 *fphi = new TF1("fphi","[0]*([1]*exp(-(x-[2])**2/[3]/[3])+(1-[1])*exp(-(x-[4])**2/[5]/[5]))",-0.1,0.1);
	fphi->SetParameters(1e-1,0,1e-3,1e-3,0,1e-2);
	//fphi->SetParameters(4.79253e-02,9.25529e-01,-8.56488e-05,-7.46701e-03,-5.37828e-04,5.99178e-0);
	//SetStyle(*pc3dphi,1.2,1,20,0,0);
	c1->cd();
	pc3dphi->Fit("fphi","RQ");
	double par[6];
	fphi->GetParameters(par);
	TF1 *fphi1 = new TF1("fphi1","gaus",-0.1,0.1);
	fphi1->SetNpx(10000);
	fphi1->SetParameters(par);
	TF1 *fphi2 = new TF1("fphi2","gaus",-0.1,0.1);
	fphi2->SetNpx(10000);
	fphi2->SetParameters(&par[3]);
	fphi1->SetLineColor(4);
	fphi1->Draw("same");
	fphi2->SetLineColor(6);
	fphi2->Draw("same");
	fout<<"dphishift["<<iarm<<"]["<<ich<<"]["<<ipt<<"] = "<<par[1]<<";  ";
	fout<<"dphisigma["<<iarm<<"]["<<ich<<"]["<<ipt<<"] = "<<par[2]<<";  ";
	if(par[1]>1) cout<<"Problem!: dphi: "<<"iarm = "<<iarm<<" ich = "<<ich<<" ipt = "<<ipt<<endl;
	c1->Print(Form("fig/pc3matching/dphi_%d_%d_%d.png",iarm,ich,ipt));
		
	TF1 *fz = new TF1("fz","gaus(0)+gaus(3)",-10,10);
	//fz->SetParameters(3.02349e+07,8.32462e-01,2.39187e+00,2.17631e+07,6.35460e-01,8.09821e+00);
	fz->SetParameters(1e-1,0,2,1e-3,0,8);
	//SetStyle(*pc3dz,1.2,1,20,0,0);
	c2->cd();
	pc3dz->Fit("fz","RQ");
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
	fout<<"dzshift["<<iarm<<"]["<<ich<<"]["<<ipt<<"] = "<<par[1]<<";  ";
	fout<<"dzsigma["<<iarm<<"]["<<ich<<"]["<<ipt<<"] = "<<par[2]<<";  ";
	if(par[1]>10) cout<<"Problem! dz: "<<"iarm = "<<iarm<<" ich = "<<ich<<" ipt = "<<ipt<<endl;
	c2->Print(Form("fig/pc3matching/dz_%d_%d_%d.png",iarm,ich,ipt));

	}
	}
	fout<<endl;
}
}
