F_center(){
ofstream foutput;
foutput.open("F_cent_fit.txt");
const int N=16;
char* crystal[N]={"LiF","LiCl","LiBr","LiI","NaF","NaCl","NaBr","NaI","KF","KCl","KBr","KI","RbF","RbCl","RbBr","RbI"};
double a[N]={0.201,0.257,0.275,0.300,0.231,0.281,0.298,0.323,0.267,0.314,0.329,0.353,0.282,0.327,0.343,0.366};	//lattice constant: nm
double F_band[N]={257,385,452,531,335,462,540,594,455,562,627,695,474,616,707,766};	//maximum F_band wavelength: nm
gStyle->SetOptFit(0);
TCanvas *c1=new TCanvas("c1");
c1->SetGridx();
c1->SetGridy();
TGraph *gr=new TGraph(N,a,F_band);
gr->SetTitle("F_band simple fit of #lambda_{max} vs a");
gr->GetXaxis()->SetTitle("lattice constant a [nm]");
gr->GetYaxis()->SetTitle("F-band Maximum #lambda [nm]");
gr->SetMarkerSize(0.8);
gr->SetMarkerStyle(20);
TF1 *fitfun=new TF1("fitfun","[0]*TMath::Power(x,[1])");
fitfun->SetParameter(1,2);
gr->Fit("fitfun","Q");
gr->Draw("AP");

TPaveText *text=new TPaveText(0.25,0.7,0.30,0.85,"brNDC");
text->AddText(Form("fit function:#lambda=%.3f*a^{%.3f}",fitfun->GetParameter(0),fitfun->GetParameter(1)));
text->SetTextSize(0.035);
text->Draw("same");

c1->SaveAs("./plots/F_band.png");

foutput<<"Constant C ="<<fitfun->GetParameter(0)<<"#pm"<<fitfun->GetParError(0)<<endl;
foutput<<"Power n ="<<fitfun->GetParameter(1)<<"#pm"<<fitfun->GetParError(1)<<endl;

foutput.close();
}

