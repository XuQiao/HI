void GetRandomTest(){
	double k0=1.39;
	double k1 = 0.425;
	double theta0 = 3.41;
	double theta1 = 1.30;
	int iNpart=2;
	double k_=k0+k1*(iNpart-2);
        double theta_=theta0+theta1*TMath::Log(iNpart-1);
	TF1 *f = new TF1("f","TMath::GammaDist(x,[0],0,[1])",0,200);
	f->SetParameters(k1,theta_);
	cout<<"Value at 0 = "<<f->Eval(0)<<endl;
	cout<<"Value at 1e-11 = "<<f->Eval(1e-11)<<endl;
	cout<<"Integral 1= "<<f->Integral(f->GetXmin(),f->GetXmax())<<endl;
	f->SetRange(1e-12,200);
	cout<<"Integral 2= "<<f->Integral(f->GetXmin(),f->GetXmax())<<endl;
	cout<<"fXmin = "<<f->GetXmin()<<"\tfXmax = "<<f->GetXmax()<<"\tfNpx = "<<f->GetNpx()<<endl;
	f->SetNpx(1e5);
	TCanvas *c1 = new TCanvas();
	c1->SetLogy();
	cout<<"f mean = "<<f->Mean(0,200)<<endl;
	cout<<"math mean = "<<f->GetParameter(0)*f->GetParameter(1)<<endl;
	TH1D* h = new TH1D("h","h",1000,0,200);
	for(int i=0;i<1e6;i++){
		double para = f->GetRandom();
		h->Fill(para);
	}
	h->Scale(1.0/h->Integral()*1000/200);
	h->GetYaxis()->SetRangeUser(1e-10,1);
	h->SetMarkerStyle(24);
	h->SetMarkerColor(4);
	h->SetMarkerSize(1.1);
	TLegend *leg = new TLegend(0.6,0.7,0.8,0.9);
        leg->SetFillColor(0);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetTextFont(42);
        leg->SetTextSize(0.03);
        leg->AddEntry(f,"function","lp");
        leg->AddEntry(h,"filled histogram","lp");
	h->Draw("P");
	f->Draw("same");
        leg->Draw("same");
	cout<<"h mean = "<<h->GetMean(1)<<endl;
	c1->Print("TestGetRandom.png");
	}
	
