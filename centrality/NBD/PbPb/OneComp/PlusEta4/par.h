//const double bin[]={0,6.63,11.27,14.62,19.82,100};	//G0.root
const double bin[]={1.0,0.5,0.3,0.2,0.1,0};	//G1.root
//const double bin[]={1.0,0};	//G2.root
int N=(int)(sizeof(bin)/sizeof(double));
int method=1;	//method=0;

const TString datafile="/scratch/xuq7/Centrality/Correction/pPbHistPlusEta4.root";

const TString histoname="hHFEnergyPlus4";

//struct para1 var={8,55, 0.20,0.25,0.01,	0.05,0.08,0.01};
struct para1 var1={8,55, 0.15,0.40,0.01,	0.01,0.20,0.01};

struct para1 var2={8,55, 0.20,0.40,0.01,	0.01,0.20,0.01};

struct para1 var3={8,55, 0.20,0.40,0.01,	0.01,0.20,0.01};

double binshift = 2;

struct para2 bestlist1[nGlau+2]={{0.23,0.06},{0.23,0.06},{0.29,0.08},{0.23,0.06},{0.24,0.06},{0.25,0.07},{0.29,0.08},{0.33,0.09},{0.22,0.06}
			      	,{0.20,0.05},{0.31,0.09}};

struct para2 bestlist2[nGlau+2]={{0.35,0.13},{0.34,0.13},{0.37,0.14},{0.35,0.13},{0.35,0.13},{0.35,0.14},{0.34,0.12},{0.39,0.14},{0.39,0.18}
				,{0.35,0.13},{0.35,0.13}};

struct para2 bestlist3[nGlau+2]={{0.38,0.18},{0.39,0.19},{0.36,0.16},{0.38,0.18},{0.38,0.18},{0.38,0.19},{0.38,0.17},{0.39,0.17},{0.35,0.18}
				,{0.39,0.18},{0.38,0.18}};

const TString outG =  "G1.root";	//"G0.root";
