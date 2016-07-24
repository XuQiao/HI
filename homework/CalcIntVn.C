TF1 *fun_fit_v2;
TF1 *fun_fit_v3;
TF1 *fun_fit_v4;
TF1 *N_pt;
TF1 *fun_pt;
TF1 *fun_fit;

void convert(int cent,char *buf){
sprintf(buf,"%d-%d%%",cent*10,(cent+1)*10);
}

void convert_ex(int cent,char *buf){
if(cent==0)     {sprintf(buf,"0-92%%");}
else if(cent==1)        {sprintf(buf,"0-5%%");}
else if(cent==2)        {sprintf(buf,"0-10%%");}
else if(cent==3)        {sprintf(buf,"10-15%%");}
else if(cent==4)        {sprintf(buf,"15-20%%");}
else if(cent==5)        {sprintf(buf,"10-20%%");}
else if(cent==6)        {sprintf(buf,"20-30%%");}
else if(cent==7)        {sprintf(buf,"30-40%%");}
else if(cent==8)        {sprintf(buf,"40-50%%");}
else if(cent==9)        {sprintf(buf,"50-60%%");}
else if(cent==10)       {sprintf(buf,"60-70%%");}
else if(cent==11)       {sprintf(buf,"70-80%%");}
else if(cent==12)       {sprintf(buf,"80-92%%");}
else if(cent==13)       {sprintf(buf,"60-92%%");}
}

void read_ex(double **pt,double **v,double **v_err,char *input,const int cent,const int* pt_,const int col,const int ptbin_max){
ifstream f;
double tmp;
f.open(input);
for(int k=0;k<cent;k++){
for(int i=0;i<pt_[k];i++){
f>>(*((double*)pt+k*ptbin_max+i));
f>>(*((double*)v+k*ptbin_max+i));
f>>(*((double*)v_err+k*ptbin_max+i));
(*((double*)v_err+k*ptbin_max+i))=(*((double*)v+k*ptbin_max+i))*(*((double*)v_err+k*ptbin_max+i))/100;
for(int j=0;j<col-3;j++)
f>>tmp;
}
}

f.close();
}

double fitf_ex(double *x,double *par){
double arg=0;
if(par[0]!=0) arg=1+x[0]/par[1];
double fitval=par[0]*TMath::Power(arg,-par[2]);
//*if(par[0]!=0) arg=1+x[0]/par[0];
//double fitval=1/(2*TMath::Pi()*TMath::Power(par[0],2))*(par[2]-1)*(par[2]-2)*par[1]*TMath::Power(arg,-par[2]);
return fitval;
}

double multiply_1(double *x,double *par){
double mult_result;
double fXmin,fXmax,gXmin,gXmax;
fun_pt->GetRange(fXmin,fXmax);
fun_fit->GetRange(gXmin,gXmax);
if(fXmin!=gXmin || fXmax!=gXmax){
cout<<"not same function range!"<<endl;
return 0;
}
else{
mult_result=fun_pt->Eval(x[0]) * fun_fit->Eval(x[0]);
return mult_result;
}

}

double multiply_2(double *x,double *par){
double mult_result;
double fXmin,fXmax,gXmin,gXmax;
N_pt->GetRange(fXmin,fXmax);
fun_fit_v2->GetRange(gXmin,gXmax);
if(fXmin!=gXmin || fXmax!=gXmax){
cout<<"not same function range!"<<endl;
return 0;
}
else{
mult_result=N_pt->Eval(x[0]) * fun_fit_v2->Eval(x[0]);
return mult_result;
}
}

double multiply_3(double *x,double *par){
double mult_result;
double fXmin,fXmax,gXmin,gXmax;
N_pt->GetRange(fXmin,fXmax);
fun_fit_v3->GetRange(gXmin,gXmax);
if(fXmin!=gXmin || fXmax!=gXmax){
cout<<"not same function range!"<<endl;
return 0;
}
else{
mult_result=N_pt->Eval(x[0]) * fun_fit_v3->Eval(x[0]);
return mult_result;
}
}

double multiply_4(double *x,double *par){
double mult_result;
double fXmin,fXmax,gXmin,gXmax;
N_pt->GetRange(fXmin,fXmax);
fun_fit_v4->GetRange(gXmin,gXmax);
if(fXmin!=gXmin || fXmax!=gXmax){
cout<<"not same function range!"<<endl;
return 0;
}
else{
mult_result=N_pt->Eval(x[0]) * fun_fit_v4->Eval(x[0]);
return mult_result;
}
}

double pol3_(double *x,double *par){
double y;
//y=par[0]*TMath::Power(x[0],3)+par[1]*TMath::Power(x[0],2)+par[2]*x[0]+par[3];
y=par[0]*TMath::Power(x[0],4)+par[1]*TMath::Power(x[0],3)+par[2]*TMath::Power(x[0],2)+par[3]*x[0]+par[4];
//y=TMath::Power((x[0]/par[0]),par[1])/TMath::Power((1+x[0]/par[2]),par[3])*(.00005+TMath::Power(1/x[0],par[4]));
return y;
}

void read(double **pt,double **v,double **v_err,char *input,const int cent,const int pt_,const int col_){
ifstream f;
double tmp;
f.open(input);
f>>tmp;
for(int k=0;k<cent;k++){
for(int i=0;i<pt_;i++){
f>>(*((double*)pt+k*pt_+i));
f>>(*((double*)v+k*pt_+i));
f>>(*((double*)v_err+k*pt_+i));
for(int j=0;j<col_-3;j++)
f>>tmp;
}
}

f.close();
}

void CalcIntVn(){
ofstream foutput ;
foutput.open("IntVn.txt");
const int col=6,ptbinm=10,cent_all=6;
const int cent=14,col_ex=5;
double pt_v2[cent_all][ptbinm],pt_v3[cent_all][ptbinm],pt_v4[cent_all][ptbinm];
double v2[cent_all][ptbinm],v2_err[cent_all][ptbinm];
double v3[cent_all][ptbinm],v3_err[cent_all][ptbinm];
double v4[cent_all][ptbinm],v4_err[cent_all][ptbinm];
TCanvas *c1[cent_all],*c2[cent_all];
TF1 *Int_V2,*Int_V3,*Int_V4;
double fitmin=0.3,fitmax=3;
fun_fit_v2=new TF1("fun_fit_v2",pol3_,fitmin,fitmax,5);
fun_fit_v3=new TF1("fun_fit_v3",pol3_,fitmin,fitmax,5);
fun_fit_v4=new TF1("fun_fit_v4",pol3_,fitmin,fitmax,5);
/*fun_fit_v2->SetParameters(3.4,1.8,3,1.8,0.8);
fun_fit_v3->SetParameters(3.4,1.8,3,1.8,0.8);
fun_fit_v4->SetParameters(3.4,1.8,3,1.8,0.8);*/
fun_fit_v2->SetParameters(0.0008,0.0466,-0.0058,0.0015,-0.00067);
fun_fit_v3->SetParameters(0.0008,0.0466,-0.0058,0.0015,-0.00067);
fun_fit_v4->SetParameters(0.0008,0.0466,-0.0058,0.0015,-0.00067);
TGraphErrors *gr_v2[cent_all],*gr_v3[cent_all],*gr_v4[cent_all];
TGraphErrors *gr_err[cent];

const int ptbin[cent]={34,34,34,33,33,34,34,32,32,32,32,31,30,33};
const int ptbin_max=34;
double pt[cent][ptbin_max],y[cent][ptbin_max],y_err[cent][ptbin_max];

double pt_err[ptbin_max]={0};
char buf[cent_all][1024],buf_ex[cent][1024];
double v2Int[cent_all],v3Int[cent_all],v4Int[cent_all];

read_ex(pt,y,y_err,"charged.txt",cent,ptbin,col_ex,ptbin_max);

fun_fit=new TF1("fun_fit",fitf_ex,fitmin,fitmax,3);
fun_pt=new TF1("fun_pt","x",fitmin,fitmax);
fun_fit->SetParameters(80,3,14);
fun_fit->SetLineColor(2);
fun_fit->SetLineWidth(1);
fun_fit->SetLineStyle(1);
read(pt_v2,v2,v2_err,"v2.txt",cent_all,ptbinm,col);
read(pt_v3,v3,v3_err,"v3.txt",cent_all,ptbinm,col);
read(pt_v4,v4,v4_err,"v4.txt",cent_all-1,ptbinm,col);

for(int i=0;i<cent_all;i++)
convert_ex(i,buf_ex[i]);

for(i=0;i<cent_all;i++){

if(i==0)
int j=2;
else 
j=i+4;

convert(i,buf[i]);
c2[i]=new TCanvas();
c2[i]->cd();
c2[i]->SetLogy();
gr_err[i]=new TGraphErrors(ptbin[j],pt[j],y[j],pt_err,y_err[j]);
gr_err[i]->Fit("fun_fit","QR");
gr_err[i]->SetTitle(Form("Invariant Yield vs p_{T} in Centrality %s",buf[i]));
gr_err[i]->GetXaxis()->SetTitle("p_{T} [GeV]");
gr_err[i]->GetYaxis()->SetTitle("invariant yield (1/2#pi p_{T}) d^{2}N/dp_{T}dy [GeV^{-2}]");
gr_err[i]->GetXaxis()->SetLimits(fitmin,fitmax);
gr_err[i]->SetMaximum(3e2);
gr_err[i]->SetMinimum(1e-3);
gr_err[i]->SetMarkerStyle(20);
gr_err[i]->Draw("AP");
N_pt=new TF1("N_pt",multiply_1,fitmin,fitmax,0);
foutput<<"N_integral="<<N_pt->Integral(0.3,3)<<endl;
TLegend *leg2 = new TLegend(0.53,0.78,0.88,0.85);
leg2->AddEntry(fun_fit,Form("fit:y(p_{T})=%.3f*(1+#frac{p_{T}}{%.3f})^{-%.3f}",fun_fit->GetParameter(0),fun_fit->GetParameter(1),fun_fit->GetParameter(2)),"lp");
leg2->SetTextSize(0.03);
leg2->SetBorderSize(0);
leg2->Draw("same");
c2[i]->Print(Form("./plots/Inclusive_charged_yield_cent_%s.png",buf[i]));

//----------------------------------------------------------------------------------------------

c1[i]=new TCanvas();
c1[i]->cd();
fun_fit_v2->SetLineStyle(2);
fun_fit_v2->SetLineWidth(0.5);
fun_fit_v2->SetLineColor(1);
gr_v2[i]=new TGraphErrors(ptbinm,pt_v2[i],v2[i],pt_err,v2_err[i]);
gr_v2[i]->Fit("fun_fit_v2","QR","");
gr_v2[i]->SetTitle(Form("v_{n} vs p_{T} in Centrality %s",buf[i]));
gr_v2[i]->GetXaxis()->SetTitle("p_{T} [GeV]");
gr_v2[i]->GetYaxis()->SetTitle("v_{n}");
gr_v2[i]->SetMinimum(0);
gr_v2[i]->GetXaxis()->SetLimits(fitmin,fitmax);
gr_v2[i]->SetMarkerStyle(20);
gr_v2[i]->SetLineColor(1);
gr_v2[i]->Draw("AP");
Int_V2=new TF1("Int_V2",multiply_2,fitmin,fitmax,0);
v2Int[i]=Int_V2->Integral(0.3,3)/N_pt->Integral(0.3,3);
foutput<<"Integrated V2 in Centrality "<<buf[i]<<" is: "<<v2Int[i]<<endl;

fun_fit_v3->SetLineStyle(2);
fun_fit_v3->SetLineWidth(0.5);
fun_fit_v3->SetLineColor(2);
gr_v3[i]=new TGraphErrors(ptbinm,pt_v3[i],v3[i],pt_err,v3_err[i]);
gr_v3[i]->Fit("fun_fit_v3","QR","");
gr_v3[i]->SetMarkerStyle(22);
gr_v3[i]->SetLineColor(2);
gr_v3[i]->Draw("P");
Int_V3=new TF1("Int_V3",multiply_3,fitmin,fitmax,0);
v3Int[i]=Int_V3->Integral(0.3,3)/N_pt->Integral(0.3,3);
foutput<<"Integrated V3 in Centrality "<<buf[i]<<" is: "<<v3Int[i]<<endl;

TLegend *leg=new TLegend(0.12,0.68,0.27,0.88);
leg->AddEntry(gr_v2[i],"v_{2}{#Phi_{2}}","lp");
leg->AddEntry(gr_v3[i],"v_{3}{#Phi_{3}}","lp");
if(i!=cent_all-1){
fun_fit_v4->SetLineStyle(2);
fun_fit_v4->SetLineWidth(0.5);
fun_fit_v4->SetLineColor(3);
gr_v4[i]=new TGraphErrors(ptbinm,pt_v4[i],v4[i],pt_err,v4_err[i]);
gr_v4[i]->Fit("fun_fit_v4","QR","");
gr_v4[i]->SetMarkerStyle(24);
gr_v4[i]->SetLineColor(3);
gr_v4[i]->Draw("P");
Int_V4=new TF1("Int_V4",multiply_4,fitmin,fitmax,0);
v4Int[i]=Int_V4->Integral(0.3,3)/N_pt->Integral(0.3,3);
foutput<<"Integrated V4 in Centrality "<<buf[i]<<" is: "<<v4Int[i]<<endl;
leg->AddEntry(gr_v4[i],"v_{4}{#Phi_{4}}","lp");

}
TPaveText *text=new TPaveText(0.40,0.80,0.60,0.88,"brNDC");
text->AddText("Use Fourth Polynomial to fit v_{n}(p_{T})");
text->SetTextSize(0.03);
text->Draw("same");

leg->Draw("same");
c1[i]->Update();
c1[i]->Print(Form("./plots/vn_pt_cent_%s.png",buf[i]));

}
TCanvas *c3=new TCanvas("c3");
TCanvas *c4=new TCanvas("c4");
TCanvas *c5=new TCanvas("c5");

double cent_bin[cent_all];
for(i=0;i<cent_all;i++){
cent_bin[i]=i*10+5;
}

c3->cd();
TGraph *IntV2_cent=new TGraph(cent_all,cent_bin,v2Int);
IntV2_cent->SetTitle("Integrate v_{2} vs Centrality");
IntV2_cent->GetXaxis()->SetTitle("Centrality %");
IntV2_cent->GetYaxis()->SetTitle("Int. v_{2}");
IntV2_cent->SetMarkerStyle(20);
IntV2_cent->SetMarkerSize(0.8);
IntV2_cent->Draw("AP");
c3->SetGridy();
c3->Print("./plots/v2_vs_cent.png");

c4->cd();
TGraph *IntV3_cent=new TGraph(cent_all,cent_bin,v3Int);
IntV3_cent->SetTitle("Integrate v_{3} vs Centrality");
IntV3_cent->GetXaxis()->SetTitle("Centrality %");
IntV3_cent->GetYaxis()->SetTitle("Int. v_{3}");
IntV3_cent->SetMarkerStyle(22);
IntV3_cent->SetMarkerSize(0.8);
IntV3_cent->Draw("AP");
c4->SetGridy();
c4->Print("./plots/v3_vs_cent.png");

c5->cd();
TGraph *IntV4_cent=new TGraph(cent_all-1,cent_bin,v4Int);
IntV4_cent->SetTitle("Integrate v_{4} vs Centrality");
IntV4_cent->GetXaxis()->SetTitle("Centrality %");
IntV4_cent->GetYaxis()->SetTitle("Int. v_{4}");
IntV4_cent->SetMarkerStyle(24);
IntV4_cent->SetMarkerSize(0.8);
IntV4_cent->Draw("AP");
c5->SetGridy();
c5->Print("./plots/v4_vs_cent.png");

}
