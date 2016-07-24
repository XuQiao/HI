const double bin[]={0.9,0.6,0.4,0.3,0.2,0.1,0};
int N=(int)(sizeof(bin)/sizeof(double));
int method=0;

const TString datafile="/afs/cern.ch/user/q/qixu/CMSSW_6_2_5/src/Centrality/Correction/MBCen_weight_RecoHI.root";

const TString histoname="pfjet_HFPlus4_QiaoHFPlus4Wt";

struct para1 var1={1,50, 0.15,0.40,0.01,	0.05,0.05,0.01};
//struct para1 var1={1,50, 0.15,0.35,0.01,	0.01,0.20,0.01};

struct para1 var2={1,50, 1.20,1.20,0.03,	1.20,1.20,0.05};

struct para1 var3={1,50, 1.20,1.20,0.03,	1.20,1.20,0.05};

struct para2 bestlist1[nGlau]={{0.20,0.05},{0.22,0.06},{0.23,0.06},{0.22,0.06},{0.20,0.05},{0.24,0.07},{0.26,0.07},{0.21,0.05},{0.21,0.06}};

struct para2 bestlist2[nGlau]={{0.35,0.13},{0.34,0.13},{0.37,0.14},{0.35,0.13},{0.35,0.13},{0.35,0.14},{0.34,0.12},{0.39,0.14},{0.39,0.18}};

struct para2 bestlist3[nGlau]={{0.38,0.18},{0.39,0.19},{0.36,0.16},{0.38,0.18},{0.38,0.18},{0.38,0.19},{0.38,0.17},{0.39,0.17},{0.35,0.18}};
