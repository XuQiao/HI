void plotnsvsV2sub(){
	TCanvas *c1 = new TCanvas();
	int isSum=1;
	if(isSum==0)	TFile *f = TFile::Open("V2_nsvssub_Prod.root");
	else TFile *f = TFile::Open("V2_nsvssub_Sum.root");
	int xbin=0;
	TVectorD* vecnsamples = (TVectorD*)f->Get(Form("D_%d/nsamples",xbin));
	TVectorD* vecVmeanmean = (TVectorD*)f->Get(Form("D_%d/Vmeanmean",xbin));
	TVectorD* vecdeltaVmeanmean = (TVectorD*)f->Get(Form("D_%d/deltaVmeanmean",xbin));
	int N = vecnsamples->GetNrows();
	double *nsamples = vecnsamples->GetMatrixArray();
	double *Vmeanmean = vecVmeanmean->GetMatrixArray();
	double *deltaVmeanmean = vecdeltaVmeanmean->GetMatrixArray();
	TGraphErrors* gr = new TGraphErrors(N,nsamples,Vmeanmean,0,deltaVmeanmean);
	gr->SetTitle("");
	gr->GetXaxis()->SetRangeUser(0,50);
	gr->GetXaxis()->SetTitle("Nsamples Div");
	gr->GetYaxis()->SetTitle("V2 mean");
	gr->SetMarkerSize(1.2);
	gr->SetMarkerStyle(20);
	gr->SetMarkerColor(1);
	gr->SetLineColor(1);
	gr->Draw("AP");
	if(isSum==0)c1->SaveAs("grnsvsV2_Prod.gif");
	else c1->SaveAs("grnsvsV2_Sum.gif");
	}
	
