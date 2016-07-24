#include <iomanip>      // std::setprecision
/*
Centrality:
0:MB  1:0-5%  2:5-10%  3:10-15%  4:15-20%  5:20-30%  6:30-40%  7:40-50%  8:50-60%  9:60-70%  10:70-80%  11:80-92%  12:60-92%

Particle:
0:pi+  1:pi-  2:K+  3:K-  4:proton  5:anti-pronton
*/

void convert_cent(int cent,char *buf){
if(cent==0){sprintf(buf,"Min. Bias");}
else if(cent==1){sprintf(buf,"0-5%%");}
else if(cent==2){sprintf(buf,"5-10%%");}
else if(cent==3){sprintf(buf,"10-15%%");}
else if(cent==4){sprintf(buf,"15-20%%");}
else if(cent==5){sprintf(buf,"20-30%%");}
else if(cent==6){sprintf(buf,"30-40%%");}
else if(cent==7){sprintf(buf,"40-50%%");}
else if(cent==8){sprintf(buf,"50-60%%");}
else if(cent==9){sprintf(buf,"60-70%%");}
else if(cent==10){sprintf(buf,"70-80%%");}
else if(cent==11){sprintf(buf,"80-92%%");}
else if(cent==12){sprintf(buf,"60-92%%");}
else;
}
void convert_part(int part,char *buf){
if(part==0){sprintf(buf,"#pi +");}
else if(part==1){sprintf(buf,"#pi -");}
else if(part==2){sprintf(buf,"K +");}
else if(part==3){sprintf(buf,"K -");}
else if(part==4){sprintf(buf,"p");}
else if(part==5){sprintf(buf,"#bar{p}");}
else;
}

TGraphErrors* MakeGraph(double ptmin, double m0,const int N,double particle_invy[],double particle_invy_err[],int MarkerStyle,int cent, int part){

double pt[N],mt_m0[N];
double x_err[N]={0};
char buf[1024],buf1[1024],buf2[1024];
convert_cent(cent,buf2);
if(part%2==0) sprintf(buf1,"Positive"); else sprintf(buf1,"Negative");

if(cent==1)
sprintf(buf,"%s Charged Particles %s Invariable Yield",buf1,buf2);
else
sprintf(buf,"%s Charged Particles %s Centrality Invariable Yield",buf1,buf2);

if(N!=22)
for(int i=0;i<N;i++){
pt[i]=ptmin+i*0.10;
mt_m0[i]=TMath::Sqrt(TMath::Power(pt[i],2)+TMath::Power(m0,2))-m0;
}
else {
for(int i=0;i<19;i++){
if(i<14)
pt[i]=ptmin+i*0.10;
pt[14]=2.10;
if(i>14&&i<19)
pt[i]=pt[14]+(i-14)*0.20;
}
pt[19]=3.25;pt[20]=3.75;pt[21]=4.25;

for(int i=0;i<N;i++)
mt_m0[i]=TMath::Sqrt(TMath::Power(pt[i],2)+TMath::Power(m0,2))-m0;
}

TGraphErrors *mass_vs_invy= new TGraphErrors(N,mt_m0,particle_invy,x_err,particle_invy_err);
mass_vs_invy->SetMarkerStyle(MarkerStyle);
mass_vs_invy->SetMarkerSize(0.8);
mass_vs_invy->GetXaxis()->SetLimits(0,4);
mass_vs_invy->SetMaximum(1e3);
if(cent<9) mass_vs_invy->SetMinimum(1e-5);
else mass_vs_invy->SetMinimum(1e-6);
mass_vs_invy->SetTitle(buf);
mass_vs_invy->GetYaxis()->SetTitle("invariant yield (1/2#pi m_{T}) d^{2}N/dm_{T}dy [GeV^{-2}]");
mass_vs_invy->GetXaxis()->SetTitle("m_{T}-m_{0} [GeV]");
return mass_vs_invy;

}

const int cent_all=13;

void read(char *filename,const int N,TArrayD* &array, TArrayD* &array_err) {

int cent=4;
ifstream fs(filename);
double tmp,tmp_err;

for(int j=0;j<cent_all;j++){
array[j]=new TArrayD(N);
array_err[j]=new TArrayD(N);
}

for(int row=0;row<cent_all/cent+1;row++){
for(int k=0;k<N;k++){
for(int j=cent*row;j<(cent*(row+1)>cent_all?cent_all:cent*(row+1));j++){
fs>>tmp;
fs>>tmp_err;
array[j]->AddAt(tmp,k);
array_err[j]->AddAt(tmp_err,k);
}
}
}

fs.close();

}

double fitf(double *x,double *par){
double arg=0;
if(par[0]!=0) arg=-x[0]/par[0];
double fitval=1/(2*TMath::Pi()*par[0]*(par[0]+par[2]))*par[1]*TMath::Exp(arg);
return fitval;
}

double linf(double *x,double *par){
double fitval=x[0]*par[0]+par[1];
return fitval;
}

void Iden_charged(){
const int N_pi=28,N_K=16,N_pro=22,part_all=6,k_max=3;
const double pimass=139.57018/1000, Kmass=493.677/1000, promass=938.272/1000;
double x_err[part_all]={0};
int Np[part_all]={N_pi,N_K,N_pro};
int markerstyle[part_all]={20,26,30};
char *input[part_all]={"piplus.txt","piminus.txt","Kplus.txt","Kminus.txt","proplus.txt","prominus.txt"};
double mass[part_all]={pimass,Kmass,promass};
double ptmin[part_all]={0.25,0.45,0.65};
double fitmin[part_all]={0.2,0.1,0.1};
double fitfixPar0[part_all]={0.2,0.3,0.4};
TCanvas *c[cent_all],*c_com[k_max];
TGraphErrors *gr_err[part_all][cent_all],*gr_inv_slope_plus[cent_all],*gr_inv_slope_minus[cent_all];
TF1 *fitfun[part_all][cent_all];
char fitname[1024][part_all][cent_all];
double T[part_all][cent_all],T_err[part_all][cent_all],T_cent[part_all],T_cent_err[part_all];
ofstream foutput;
foutput.open("results.txt");
char buf1[1024][part_all],buf2[1024][cent_all];

//------------------------------Miscellaneous---------------------------------------

foutput<<fixed<<setprecision(5);
foutput<<"Inverse Slope T"<<endl;
for(int i=0;i<part_all;i++)
convert_part(i,buf1[i]);
for(int j=0;j<cent_all;j++)
convert_cent(j,buf2[j]);

//---------------------------------read from files------------------------------------

TArrayD *array[part_all][cent_all],*array_err[part_all][cent_all];
double *arr[part_all][cent_all],*arr_err[part_all][cent_all];

for(int i=0;i<part_all;i++)
read(input[i],Np[i/2],array+cent_all*i,array_err+cent_all*i);

//-------------------------------fit data and get inverse slope---------------

for(int i=0;i<part_all;i++){		//particle loop
foutput<<buf1[i]<<":"<<endl;
for(int j=0;j<cent_all;j++){		//centrality bin loop
sprintf(fitname[i][j],"fit_%d_cent_%d",i,j);
fitfun[i][j]=new TF1(fitname[i][j],fitf,fitmin[i/2],1.0,3);
fitfun[i][j]->SetParameter(0,fitfixPar0[i/2]);
fitfun[i][j]->FixParameter(2,mass[i/2]);
fitfun[i][j]->SetParNames("T","A","m0");

arr[i][j]=array[i][j]->GetArray();
arr_err[i][j]=array_err[i][j]->GetArray();

gr_err[i][j]=MakeGraph(ptmin[i/2],mass[i/2],Np[i/2],arr[i][j],arr_err[i][j],markerstyle[i/2],j,i);
gr_err[i][j]->Fit(fitname[i][j],"QRM","");
T[i][j]=gr_err[i][j]->GetFunction(fitname[i][j])->GetParameter(0);
T_err[i][j]=gr_err[i][j]->GetFunction(fitname[i][j])->GetParError(0);
foutput<<buf2[j]<<":";
foutput<<T[i][j]<<"pm";
foutput<<T_err[i][j]<<'\t';
}
foutput<<endl;
}

//--------------------------------get flow velocity and T0---------------------------

foutput<<endl<<"flow velocity and T0:"<<endl;
TF1 *lin=new TF1("linfun",linf,0,1,2);

for(j=0;j<cent_all;j++){

TLegend *leg1 = new TLegend(0.7,0.7,0.85,0.85);
TLegend *leg2 = new TLegend(0.7,0.7,0.85,0.85);

for(i=0;i<part_all;i++){
if(i<3){
T_cent[i]=T[i*2][j];
T_cent_err[i]=T_err[i*2][j];
}
else{
T_cent[i]=T[i*2-part_all+1][j];
T_cent_err[i]=T_err[i*2-part_all+1][j];
}
}

gr_inv_slope_plus[j]=new TGraphErrors(part_all/2,mass,T_cent,x_err,T_cent_err);
gr_inv_slope_plus[j]->Fit("linfun","Q");
gr_inv_slope_plus[j]->SetTitle(Form("Positive Charged Particles in %s",buf2[j]));
gr_inv_slope_plus[j]->GetXaxis()->SetTitle("Mass [GeV]");
gr_inv_slope_plus[j]->GetYaxis()->SetTitle("Inverse Slope T [GeV]");
gr_inv_slope_plus[j]->SetMinimum(0.1);
if(j<8)
gr_inv_slope_plus[j]->SetMaximum(0.5);
else
gr_inv_slope_plus[j]->SetMaximum(0.4);
gr_inv_slope_plus[j]->SetMarkerStyle(j%8+20);
gr_inv_slope_plus[j]->SetMarkerSize(0.8);
gr_inv_slope_plus[j]->SetLineStyle(2);
gr_inv_slope_plus[j]->SetLineColor(j%8+1);
gr_inv_slope_plus[j]->GetFunction("linfun")->SetLineStyle(2);
gr_inv_slope_plus[j]->GetFunction("linfun")->SetLineColor(j%8+1);
gr_inv_slope_plus[j]->GetFunction("linfun")->SetLineWidth(0.8);

foutput<<buf2[j]<<": positive flow:"<<TMath::Sqrt(gr_inv_slope_plus[j]->GetFunction("linfun")->GetParameter(0))<<'\t';
foutput<<buf2[j]<<": positive T0:"<<gr_inv_slope_plus[j]->GetFunction("linfun")->GetParameter(1)<<endl;

gr_inv_slope_minus[j]=new TGraphErrors(part_all/2,mass,T_cent+part_all/2,x_err,T_cent_err+part_all/2);
gr_inv_slope_minus[j]->Fit("linfun","Q");
gr_inv_slope_minus[j]->SetTitle(Form("Negative Charged Particles in %s",buf2[j]));
gr_inv_slope_minus[j]->GetXaxis()->SetTitle("Mass [GeV]");
gr_inv_slope_minus[j]->GetYaxis()->SetTitle("Inverse Slope T");
gr_inv_slope_minus[j]->SetMinimum(0.1);
if(j<8)
gr_inv_slope_minus[j]->SetMaximum(0.5);
else
gr_inv_slope_minus[j]->SetMaximum(0.4);
gr_inv_slope_minus[j]->SetMarkerStyle(j%8+20);
gr_inv_slope_minus[j]->SetMarkerSize(0.8);
gr_inv_slope_minus[j]->SetLineStyle(2);
gr_inv_slope_minus[j]->SetLineColor(j%8+1);
gr_inv_slope_minus[j]->GetFunction("linfun")->SetLineStyle(2);
gr_inv_slope_minus[j]->GetFunction("linfun")->SetLineColor(j%8+1);
gr_inv_slope_minus[j]->GetFunction("linfun")->SetLineWidth(0.8);
foutput<<buf2[j]<<": negative flow:"<<TMath::Sqrt(gr_inv_slope_minus[j]->GetFunction("linfun")->GetParameter(0))<<'\t';
foutput<<buf2[j]<<": negative T0:"<<gr_inv_slope_minus[j]->GetFunction("linfun")->GetParameter(1)<<endl;
}

//------------------------------------plots----------------------------------------------------------------------------

for(int j=0;j<cent_all;j++){           //centrality bin loop

c[j]=new TCanvas();
c[j]->Divide(2);

TLegend *leg1 = new TLegend(0.7,0.7,0.85,0.85);
TLegend *leg2 = new TLegend(0.7,0.7,0.85,0.85);

for(int i=0;i<part_all;i++){           //particle loop

if(i%2==0){
c[j]->cd(1)->SetLogy();

if(i==0){
gr_err[i][j]->Draw("AP");

}
else
gr_err[i][j]->Draw("P");

leg1->AddEntry(gr_err[i][j],buf1[i],"p");
leg1->Draw("same");
}

else{
c[j]->cd(2)->SetLogy();

if(i==1){
gr_err[i][j]->Draw("AP");
}
else
gr_err[i][j]->Draw("P");

leg2->AddEntry(gr_err[i][j],buf1[i],"p");
leg2->Draw("same");
}


}

c[j]->Print(Form("./plots/inv_yield_cent_%s.png",buf2[j]));
delete c[j];

}

for(j=0;j<cent_all;j++){
c[j]=new TCanvas();
c[j]->Divide(2);
c[j]->cd(1);
gr_inv_slope_plus[j]->Draw("AP");
c[j]->cd(2);
gr_inv_slope_minus[j]->Draw("AP");
c[j]->Print(Form("./plots/inv_slope_vs_mass_cent_%s.png",buf2[j]));
delete c[j];
}

for(int k=0;k<k_max;k++){
c_com[k]=new TCanvas();
TLegend *leg3=new TLegend(0.15,0.65,0.35,0.85);
TLegend *leg4=new TLegend(0.15,0.65,0.35,0.85);
leg3->SetTextSize(0.04);
leg4->SetTextSize(0.04);
c_com[k]->Divide(2);
for(j=0;j<cent_all;j++){
if(j==k+1||j==k+4||j==k+7||j==k+10){
c_com[k]->cd(1);
gr_inv_slope_plus[j]->SetTitle("Postive Charged Particles");
if(j==k+1) gr_inv_slope_plus[j]->Draw("AP");
else  gr_inv_slope_plus[j]->Draw("P");
leg3->AddEntry(gr_inv_slope_plus[j],buf2[j],"lp");
leg3->Draw("same");
c_com[k]->cd(2);
gr_inv_slope_minus[j]->SetTitle("Negative Charged Particles");
if(j==k+1) gr_inv_slope_minus[j]->Draw("AP");
else gr_inv_slope_minus[j]->Draw("P");
leg4->AddEntry(gr_inv_slope_minus[j],buf2[j],"lp");
leg4->Draw("same");
}

}
c_com[k]->Print(Form("./plots/inv_slope_vs_mass_combined_%d.png",k));
}

}
