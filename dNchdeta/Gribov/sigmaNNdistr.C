void sigmaNNdistr(){
	double omega[2] = {0.55,1.01};
	double sigNN = 70;
	double sigtot[2] = {86,94.8};
	double sig0[2] = {78.6,72.5};
	TF1 *Phsig[2];
	TCanvas *c1 = new TCanvas();
	for(int i=0;i<2;i++){
		double lambda = sigNN/sigtot[i];
		cout<<"lambda="<<lambda<<endl;
		Phsig[i] = new TF1("Phsig","[0]*x/[1]/(x/[1]+[2])*exp(-((x/[1]/[2]-1)/[3])^2)",0,300);
		Phsig[i]->SetParameters(1,lambda,sig0[i],omega[i]);
		Phsig[i]->SetTitle("");
		Phsig[i]->SetLineColor(i+1);
		double norm = Phsig[i]->Integral(0,300);
		Phsig[i]->SetParameter(0,1.0/norm);
	}
		Phsig[1]->GetYaxis()->SetRangeUser(0,0.02);
		Phsig[1]->GetYaxis()->SetTitle("P_{H}(#sigma_{NN})");
		Phsig[1]->GetXaxis()->SetTitle("#sigma_{NN} (mb)");
		TLegend *tl = new TLegend(0.6,0.7,0.8,0.8);
		tl->AddEntry(Phsig[0],Form("Gribov #Omega = %.2f",omega[0]),"l");
		tl->AddEntry(Phsig[1],Form("Gribov #Omega = %.2f",omega[1]),"l");
		tl->SetBorderSize(0);
		tl->SetFillColor(0);
		tl->SetTextSize(0.05);
		Phsig[1]->GetYaxis()->SetTitleOffset(1.2);
		Phsig[1]->Draw();
		Phsig[0]->Draw("same");
		tl->Draw("same");
		c1->SaveAs("Gribov_sigmaDis.pdf");
		c1->SaveAs("Gribov_sigmaDis.png");
}
