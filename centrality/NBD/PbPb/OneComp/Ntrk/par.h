const double bin[]={1,35,54,67,86,310};	//G0.root
//const double bin[]={1.0,0.5,0.3,0.2,0.1,0};	//G1.root
//const double bin[]={1.0,0};	//G2.root
int N=(int)(sizeof(bin)/sizeof(double));
int method=1;

const TString datafile="/scratch/xuq7/Centrality/pPbHist.root";

const TString histoname="hNtrack";

//struct para1 var={8,55, 0.20,0.25,0.01,	0.05,0.08,0.01};
struct para1 var1={20,200, 3.80,4.20,0.01,	0.30,0.60,0.01};

struct para1 var2={20,200, 4.9,5.3,0.01,	2.80,3.40,0.01};

struct para1 var3={20,200, 4.5,5.0,0.01,	1.90,2.30,0.01};

double binshift = 5;

struct para2 bestlist1[nGlau+2]=
//{{3.95,0.42},{3.87,0.41},{4.20,0.47},{4.09,0.46},{3.93,0.38},{3.85,0.42},{4.19,0.45},{4.15,0.40},{3.81,0.44},{4.01,0.43},{3.89,0.40}};
{{3.942,0.417},{3.861,0.409},{4.205,0.469},{4.084,0.457},{3.939,0.376},{3.842,0.414},{4.182,0.441},{4.158,0.399},{3.803,0.434},{4.017,0.428},{3.882,0.397}};

struct para2 bestlist2[nGlau+2]=
//{{5.19,3.10},{5.06,3.00},{5.27,3.26},{5.19,3.10},{5.19,3.10},{5.04,3.29},{5.29,3.01},{5.30,2.80},{4.90,3.03},{5.23,3.31},{5.13,2.95}};
{{5.18,3.11},{5.061,2.99},{5.271,3.251},{5.18,3.11},{5.18,3.11},{5.043,3.293},{5.299,3.009},{5.309,2.791},{4.89,3.04},{5.236,3.315},{5.136,2.955}};

struct para2 bestlist3[nGlau+2]=
//{{4.76,2.26},{4.63,2.14},{4.85,2.28},{4.76,2.26},{4.76,2.26},{4.62,2.30},{4.91,2.25},{5.00,1.94},{4.51,2.28},{4.85,2.22},{4.66,2.25}};
{{4.755,2.268},{4.635,2.145},{4.846,2.286},{4.755,2.268},{4.755,2.268},{4.61,2.304},{4.918,2.257},{5.009,1.933},{4.502,2.289},{4.847,2.22},{4.667,2.242}};

TString outG = "G0.root";	//G0.root
