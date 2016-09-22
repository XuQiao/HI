void SetStyle(TH1& h, double size, int color, int style, int fillstyle=0, int linestyle=1){
	h.SetMarkerSize(size);
	h.SetMarkerColor(color);
	h.SetLineColor(color);
	h.SetMarkerStyle(style);
	h.SetFillStyle(fillstyle);
	h.SetLineStyle(linestyle);
	h.GetXaxis()->SetTitleFont(42);
	h.GetYaxis()->SetTitleFont(42);
	h.GetXaxis()->SetTitleSize(0.048);
	h.GetYaxis()->SetTitleSize(0.048);
	h.GetXaxis()->CenterTitle();
	h.GetYaxis()->CenterTitle();
}

void SetStyle(TGraph& g, double size, int color, int style, int fillstyle=0, int linestyle=1){
	g.SetMarkerSize(size);
	g.SetMarkerColor(color);
	g.SetLineColor(color);
	g.SetMarkerStyle(style);
	g.SetFillStyle(fillstyle);
	g.SetLineStyle(linestyle);
	g.GetXaxis()->SetTitleFont(42);
	g.GetYaxis()->SetTitleFont(42);
	g.GetXaxis()->SetTitleSize(0.048);
	g.GetYaxis()->SetTitleSize(0.048);
	g.GetXaxis()->CenterTitle();
	g.GetYaxis()->CenterTitle();
}

void SetTitle(TH1& h, TString Xtitle, TString Ytitle, TString title){
	h.GetXaxis()->SetTitle(Xtitle);
	h.GetYaxis()->SetTitle(Ytitle);
	h.SetTitle(title);
}

void SetTitle(TGraph& g, TString Xtitle, TString Ytitle, TString title){
	g.GetXaxis()->SetTitle(Xtitle);
	g.GetYaxis()->SetTitle(Ytitle);
	g.SetTitle(title);
}

void SetTitle(TH1& h, TString Xtitle, TString Ytitle, TString title){
	h.GetXaxis()->SetTitle(Xtitle);
	h.GetYaxis()->SetTitle(Ytitle);
	h.SetTitle(title);
}

void SetTitle(TGraph& g, TString Xtitle, TString Ytitle, TString title){
	g.GetXaxis()->SetTitle(Xtitle);
	g.GetYaxis()->SetTitle(Ytitle);
	g.SetTitle(title);
}

void SetRange(TGraph &g, double xmin, double ymin, double xmax, double ymax){
	g.GetXaxis()->SetRangeUser(xmin,xmax);
	g.GetYaxis()->SetRangeUser(ymin,ymax);
}

void SetRange(TH1 &h, double xmin, double ymin, double xmax, double ymax){
	h.GetXaxis()->SetRangeUser(xmin,xmax);
	h.GetYaxis()->SetRangeUser(ymin,ymax);
}
void SetXRange(TGraph &g, double xmin, double xmax){
	g.GetXaxis()->SetRangeUser(xmin,xmax);
}

void SetXRange(TH1 &h, double xmin, double xmax){
	h.GetXaxis()->SetRangeUser(xmin,xmax);
}

void SetYRange(TGraph &g, double ymin, double ymax){
	g.GetYaxis()->SetRangeUser(ymin,ymax);
}

void SetYRange(TH1 &h, double ymin, double ymax){
	h.GetYaxis()->SetRangeUser(ymin,ymax);
}

TCanvas* CompareTwoHisto(TH1 *h1, TH1 *h2){
	TCanvas *c1 = new TCanvas();
	c1->cd();
	h1->Draw("P");
	h2->Draw("Psame");
	return c1;
}

TGraph* DivideTwoGraphs(TGraph *g1, TGraph *g2){
        int N = g1->GetN();
        double ratio[100], err[100], x[100];
	for(int i=1;i<N;i++){
            if(g2->GetY()[i] * g1->GetY()[i]==0) continue;
            ratio[i] = g1->GetY()[i]/g2->GetY()[i];
            x[i] = g1->GetX()[i];
            err[i] = ratio[i]*sqrt(TMath::Power((g1->GetEY()[i]/g1->GetY()[i]),2)+TMath::Power((g2->GetEY()[i]/g2->GetY()[i]),2));
        }
        TGraphErrors *gr = new TGraphErrors(N,x,ratio,0,err);
        SetStyle(*gr,g1->GetMarkerSize(),g1->GetMarkerColor(),g1->GetMarkerStyle());
	return gr;
}
