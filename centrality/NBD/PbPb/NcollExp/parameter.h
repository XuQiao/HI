const int nGlau=9;
//double centbin[]={0,0.01,0.05,0.10,0.20,0.30,0.40,0.60,0.90,0,0.90};
//double centbin[]={0,0.0246,0.0941, 0.166, 0.274,0.597,1.0};
//int method=0;	//0 for centrality percent, 1 for boundary
//int N=(int)(sizeof(centbin)/sizeof(double));
TString cenvarlist[3]={"ds","ss","Ntrk"};
TString cenvardir[3]={"double_side","PlusEta4_Rebin","Ntrk"};
TString cenvariable[3]={"HF Energy |#eta|>4","HF Energy #eta>4","N_{trk}^{offline}"};
TString Grilist[3]={"Standard Glauber","Glauber Gribov #Omega=0.55","Glauber Gribov #Omega=1.01"};
struct para1{
	double xmin;	double xmax;
	double mumin;	double mumax;	double mustep;
        double kmin; 	double kmax;	double kstep;
        double fmin;    double fmax;    double fstep;
};
struct para2{
	double mubest; double kbest; double fbest;
};

const int color[nGlau+2]={2,1,1,2,2,4,4,6,6,3,3};
const int marker[nGlau+2]={29,20,24,21,25,27,33,28,34,20,24};
const int color_[10]={2,3,4,5,6,7,8,9,2,3};

TString stdGlaulist[nGlau]={
/* "/lio/lfs/cms/store/user/qixu/Glauber/PbPb/glauberPbPb_Default_1M.root",
 "/lio/lfs/cms/store/user/qixu/Glauber/PbPb/glauberPbPb_SkinDepth536fm_1M.root",
 "/lio/lfs/cms/store/user/qixu/Glauber/PbPb/glauberPbPb_SkinDepth556fm_1M.root",
 "/lio/lfs/cms/store/user/qixu/Glauber/PbPb/glauberPbPb_dMin00fm_1M.root",
 "/lio/lfs/cms/store/user/qixu/Glauber/PbPb/glauberPbPb_dMin08fm_1M.root",
 "/lio/lfs/cms/store/user/qixu/Glauber/PbPb/glauberPbPb_NuclearRadius656fm_1M.root",
 "/lio/lfs/cms/store/user/qixu/Glauber/PbPb/glauberPbPb_NuclearRadius668fm_1M.root",
 "/lio/lfs/cms/store/user/qixu/Glauber/PbPb/glauberPbPb_SigmaInelNN65mb_1M.root",
 "/lio/lfs/cms/store/user/qixu/Glauber/PbPb/glauberPbPb_SigmaInelNN75mb_1M.root"*/
"/store/user/tuos/GlauberMC/PbPb2760GeV/glauberPbPb_Default_1M.root",
"/store/user/tuos/GlauberMC/PbPb2760GeV/glauberPbPb_R649_1M.root",
"/store/user/tuos/GlauberMC/PbPb2760GeV/glauberPbPb_R675_1M.root",
"/store/user/tuos/GlauberMC/PbPb2760GeV/glauberPbPb_SD4914_1M.root",
"/store/user/tuos/GlauberMC/PbPb2760GeV/glauberPbPb_SD6006_1M.root",
"/store/user/tuos/GlauberMC/PbPb2760GeV/glauberPbPb_SigmaNN59_1M.root",
"/store/user/tuos/GlauberMC/PbPb2760GeV/glauberPbPb_SigmaNN69_1M.root",
"/store/user/tuos/GlauberMC/PbPb2760GeV/glauberPbPb_dmin00_1M.root",
"/store/user/tuos/GlauberMC/PbPb2760GeV/glauberPbPb_dmin08_1M.root"
};

TString Gri055Glaulist[nGlau]={
 "/lio/lfs/cms/store/user/qixu/Glauber/glaubergribovpPbv1omega055_default.root",
 "/lio/lfs/cms/store/user/qixu/Glauber/glaubergribovpPbv1omega055_D04914.root",
 "/lio/lfs/cms/store/user/qixu/Glauber/glaubergribovpPbv1omega055_D06006.root",
 "/lio/lfs/cms/store/user/qixu/Glauber/glaubergribovpPbv1omega055_dmin00.root",
 "/lio/lfs/cms/store/user/qixu/Glauber/glaubergribovpPbv1omega055_dmin08.root",
 "/lio/lfs/cms/store/user/qixu/Glauber/glaubergribovpPbv1omega055_R649.root",
 "/lio/lfs/cms/store/user/qixu/Glauber/glaubergribovpPbv1omega055_R675.root",
 "/lio/lfs/cms/store/user/qixu/Glauber/glaubergribovpPbv1omega055_sigmaNN65.root",
 "/lio/lfs/cms/store/user/qixu/Glauber/glaubergribovpPbv1omega055_sigma75.root"
};


TString Gri101Glaulist[nGlau]={
 "/lio/lfs/cms/store/user/qixu/Glauber/glaubergribovpPbv1omega101_default.root",
 "/lio/lfs/cms/store/user/qixu/Glauber/glaubergribovpPbv1omega101_D04914.root",
 "/lio/lfs/cms/store/user/qixu/Glauber/glaubergribovpPbv1omega101_D06006.root",
 "/lio/lfs/cms/store/user/qixu/Glauber/glaubergribovpPbv1omega101_dmin00.root",
 "/lio/lfs/cms/store/user/qixu/Glauber/glaubergribovpPbv1omega101_dmin08.root",
 "/lio/lfs/cms/store/user/qixu/Glauber/glaubergribovpPbv1omega101_R649.root",
 "/lio/lfs/cms/store/user/qixu/Glauber/glaubergribovpPbv1omega101_R675.root",
 "/lio/lfs/cms/store/user/qixu/Glauber/glaubergribovpPbv1omega101_sigmaNN65.root",
 "/lio/lfs/cms/store/user/qixu/Glauber/glaubergribovpPbv1omega101_sigma75.root"
};
