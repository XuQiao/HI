void plotV2vsNtrk(){

const int nbin=5;
const int nbin24 = 12;
const double trkbin[nbin]={125.5,167.5,202.5,240,280};
const double avgtrkbin[nbin24]={44.36,54.33,68.63,88.39,108.2,131.3,162.1,196.6,227.5,247.2,269.2,301.2};
const double V24[nbin24]={0.02965,0.03913,0.04832,0.04941,0.04822,0.04955,0.049,0.04805,0.04709,0.04665,0.04772,0.04797};
const double V24err[nbin24]={0.005967,0.004262,0.001496,0.001351,0.001599,0.0003922,0.0003065,0.0002939,0.0004568,0.0008684,0.001471,0.004329};
const double V2[nbin]={0.0465,0.0498,0.0430,0.0447,0.0462};
const double V2Prod[nbin]={0.0536,0.0514,0.0491,0.0482,0.0462};
const double V2Proderr[nbin]={0.0013,0.0015,0.0015,0.0070,0.0027};
TCanvas *c1 = new TCanvas;
TGraphErrors *gr=new TGraphErrors(nbin,trkbin,V2,0,0);
TGraphErrors *grProd=new TGraphErrors(nbin,trkbin,V2Prod,0,V2Proderr);
TGraphErrors *grv24=new TGraphErrors(nbin24,avgtrkbin,V24,0,V24err);
TLatex *tl = new TLatex(0.2,0.2,Form("%.1f < p_{T} < %.1f (GeV/c)",0.3,6.0));
tl->SetNDC();
grv24->SetTitle("V_{2} vs Ntrkoffline");
grv24->GetXaxis()->SetTitle("Ntrkoffline");
grv24->GetYaxis()->SetTitle("V_{2}");
grv24->SetMaximum(0.08);
grv24->SetMinimum(0.0);
gr->SetMarkerSize(1);
grv24->SetMarkerSize(1);
gr->SetMarkerColor(1);
grProd->SetMarkerColor(4);
grv24->SetMarkerColor(2);
gr->SetMarkerStyle(20);
grProd->SetMarkerStyle(29);
grv24->SetMarkerStyle(20);
TLegend *leg = new TLegend(0.2,0.7,0.4,0.8);
leg->SetTextSize(0.05);
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->AddEntry(gr,"LYZ Sum method","P");
leg->AddEntry(grProd,"LYZ Prod method","P");
leg->AddEntry(grv24,"4-particle cumulant","P");
grv24->Draw("AP");
gr->Draw("Psame");
grProd->Draw("Psame");
tl->Draw("same");
leg->Draw("same");
c1->Print("V2vsNtrk.png");
tl->Draw();

}

