void plot(){
int xmin=110;
int xmax=750;
TFile *f_atlas_fulljet=new TFile("atlas_fulljet.root");
TFile *f_cms=new TFile("cms.root");
TGraphAsymmErrors *Ratio_atlas_fulljet=(TGraphAsymmErrors*)f_atlas_fulljet->Get("Ratio_0_03");
TF1 *fConstant_1 = Ratio_atlas_fulljet->GetFunction("pol0");
TGraphAsymmErrors *Ratio_cms=(TGraphAsymmErrors*)f_cms->Get("Ratio_0_05");
TF1 *fConstant_2 = Ratio_cms->GetFunction("pol0");
TCanvas *c1=new TCanvas("c1");
Ratio_atlas_fulljet->SetMarkerStyle(22);
Ratio_atlas_fulljet->SetMarkerSize(1);
Ratio_atlas_fulljet->SetMarkerColor(kRed);
Ratio_atlas_fulljet->GetXaxis()->SetLimits(xmin,xmax);
Ratio_atlas_fulljet->SetMinimum(0.6);
Ratio_atlas_fulljet->SetMaximum(1.0);
fConstant_1->SetLineColor(kRed);
fConstant_1->Draw();
Ratio_atlas_fulljet->Draw("AP");
Ratio_cms->SetMarkerStyle(20);
Ratio_cms->SetMarkerSize(1);
Ratio_cms->SetMarkerColor(kBlack);
fConstant_2->SetLineColor(kBlack);
fConstant_2->Draw("same");
Ratio_cms->Draw("P");
TLegend *leg = new TLegend(0.8,0.8,0.95,0.95);
leg->AddEntry(Ratio_atlas_fulljet,"Ratio Plot for Atlas fulljet","lp");
leg->AddEntry(Ratio_cms,"Ratio Plot for CMS","lp");
leg->Draw("same");
c1->Print("plots/combined/ratio_atlasfull_vs_cms.png");
}

