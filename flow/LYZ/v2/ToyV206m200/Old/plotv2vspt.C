void plotv2vspt(){
const double ptbin[]={0.2,0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0,10.0};
//const double ptbin[]={0.2,0.3,0.5,0.8,1.0,1.2,1.5,1.8,2.0,2.2,2.5,2.8,3.0,3.3,3.8,4.0,4.5,5.0,5.4,6.0,7.0,8.0,9.0,10.0};
const int npt=10;

ifstream fstr("v.txt",ios::in);
vector <double> v2vec;
vector <double> v2errvec;
vector <double> ptvec;
double n;double nerr;

while(!fstr.eof()){
	fstr>>n;
	fstr>>nerr;
    v2vec.push_back(n);
    v2errvec.push_back(nerr);
}
const int N=v2vec.size()-1;
double v2[N],v2err[N],pt[N];
//cout<<N<<endl;
for(int i=0;i<npt;i++){
v2[i]=v2vec[i];//cout<<v2[i]<<"\t";
v2err[i]=v2errvec[i];//cout<<v2[i]<<"\t";
pt[i]=(ptbin[i]+ptbin[i+1])/2;//cout<<pt[i]<<endl;
}
TGraphErrors *gr=new TGraphErrors(N,pt,v2,0,v2err);
gr->SetTitle("v_{2} vs momentum");
gr->GetXaxis()->SetTitle("p_{T} (GeV/c)");
gr->GetYaxis()->SetTitle("v_{2}");
gr->SetMarkerSize(1);
gr->SetMarkerColor(1);
gr->SetMarkerStyle(20);
gr->Draw("AP");

TF1 *V2vsPt = new TF1("V2vsPt","((x/3.4)^1.8/(1+(x/3)^1.8))*(.00005+(1/x)^0.8)",0.2,10);
V2vsPt->SetLineColor(2);
V2vsPt->Draw("same");
c1->Print("v2vspt.png");

fstr.close();

}

