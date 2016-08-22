#include <TGraph.h>
#include "fstream"
#include "iostream"
void plotgvsr(){

ifstream fstr("G2.txt");
vector <double> gvec;
vector <double> rvec;
double n;

for(int i=0;i<500;i++){
//while(!fstr.eof()){
	fstr>>n;
    rvec.push_back(n);
	fstr>>n;
    gvec.push_back(n);
}
if(gvec.size()==rvec.size()){
const int N=gvec.size()-1;
double g[N],r[N];
//cout<<N<<endl;
for(int i=0;i<N;i++){
r[i]=rvec[i];cout<<r[i]<<"\t";
g[i]=gvec[i];cout<<g[i]<<endl;
}
TGraph *gr=new TGraph(N,r,g);
gr->GetXaxis()->SetTitle("r");
gr->GetYaxis()->SetTitle("|G^{#theta}(ir)|");
gr->SetTitle("Title");
gr->SetMarkerSize(1);
gr->SetMarkerColor(1);
gr->SetMarkerStyle(20);
gr->Draw("AP");
//c1->Print("gvsr.png");
}
fstr.close();

}

