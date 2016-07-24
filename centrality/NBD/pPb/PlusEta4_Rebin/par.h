const double bin[]={0,6.63,11.27,14.62,19.82,100};	//G0.root
//const double bin[]={1.0,0.5,0.3,0.2,0.1,0};	//G1.root
//const double bin[]={1.0,0};	//G2.root
int N=(int)(sizeof(bin)/sizeof(double));
int method=1;	//method=0 using percent method=1 using boundary

const TString datafile="/scratch/xuq7/Centrality/pPbHistPlusEta4.root";

const TString histoname="hHFEnergyPlus4_Rebin";

//struct para1 var={8,55, 0.20,0.25,0.01,	0.05,0.08,0.01};
struct para1 var1={12,70, 0.20,0.40,0.01,      	0.01,0.20,0.01};

struct para1 var2={12,70, 0.25,0.45,0.01,       0.10,0.25,0.01};

struct para1 var3={12,70, 0.35,0.55,0.01,       0.15,0.35,0.01};

double binshift = 3;

struct para2 bestlist1[nGlau+2]={{0.3,0.085},{0.31,0.09},{0.32,0.09},{0.31,0.09},{0.30,0.08},{0.30,0.09},{0.35,0.1},{0.36,0.1},{0.3,0.09},{0.31,0.09},{0.34,0.1}};

struct para2 bestlist2[nGlau+2]={{0.375,0.146},{0.384,0.156},{0.412,0.165},{0.375,0.146},{0.375,0.146},{0.331,0.125},{0.415,0.163},{0.435,0.166},{0.301,0.114},{0.373,0.145},{0.401,0.163}};

struct para2 bestlist3[nGlau+2]={{0.485,0.286},{0.499,0.319},{0.488,0.28},{0.485,0.286},{0.485,0.286},{0.472,0.292},{0.521,0.31},{0.537,0.299},{0.465,0.312},{0.489,0.29},{0.488,0.29}};

const TString outG =  "G0.root";
