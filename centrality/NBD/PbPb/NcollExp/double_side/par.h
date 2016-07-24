//const double bin[]={0,12.61,18.99,23.23,29.50,150};	//G1.root
double bin[]={1.0,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0};	//G0.root
//const double bin[]={1.0,0};	//G2.root
int N=(int)(sizeof(bin)/sizeof(double));
int method=0; //method=0: use bin; method=1: use kpoint

TString datafile="/scratch/xuq7/Centrality/PbPbHistPlusEta4.root";

TString histoname="hHFEnergy";

double binshift = 20;

struct para1 var1={40,4000,	0,10,1,    0,5,0.5, 0.5,1.5,0.1};  //standard Glauber
struct para2 bestlist1[nGlau+2]={
    {3.1,2.42,0.62},{4,1.02,0.8},{3.08,1.58,0.62},{3.08,1.8,0.64},{3.02,1.98,0.6},{3.02,0.58,0.68},{3.9,1.54,0.78},{3.1,2.44,0.62},{3.02,1.98,0.6},{3.1,2.42,0.62},{3.1,2.42,0.62}
};

	
struct para1 var2= {15,70,	1.00,2.01,0.50,	2.01,5.01,0.2, 0.9,1.1,0.05};   //Gribov Omega = 0.55
struct para2 bestlist2[nGlau+2]={
    {1.5,3.52},{1.5,3.916},{1.503,3.008},{1.5,3.52},{1.5,3.52},{1.402,2.604},{1.5,2.42},{1.601,3.117},{1.4,3.913},{1.5,3.12},{1.5,3.52}
};


struct para1 var3={15,70,	1.00,2.01,0.50,	2.01,5.01,0.2, 0.9,1.1,0.05};   //Gribov Omega = 1.01
struct para2 bestlist3[nGlau+2]={
    {1.3,2.419},{1.3,2.918},{1.4,3.919},{1.3,2.419},{1.3,2.419},{1.3,3.42},{1.4,3.418},{1.5,3.92},{1.2,2.22},{1.4,3.82},{1.3,2.419}
};


TString outG = Form("G%d.root",method);
