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
};
struct para2{
	double mubest; double kbest;
};

const int color[nGlau+2]={2,1,1,2,2,4,4,6,6,3,3};
const int marker[nGlau+2]={29,20,24,21,25,27,33,28,34,20,24};
const int color_[10]={2,3,4,5,6,7,8,9,10};

TString stdGlaulist[nGlau]={
 "/lio/lfs/cms/store/user/qixu/centrality/Glauber/pPb/systematics/Phob_Glau_pPb_sNN70mb_v15_1M_dmin04.root",
 "/lio/lfs/cms/store/user/qixu/centrality/Glauber/pPb/systematics/Phob_Glau_pPb_sNN70mb_v15_1M_D04914.root",
 "/lio/lfs/cms/store/user/qixu/centrality/Glauber/pPb/systematics/Phob_Glau_pPb_sNN70mb_v15_1M_D06006.root",
 "/lio/lfs/cms/store/user/qixu/centrality/Glauber/pPb/systematics/Phob_Glau_pPb_sNN70mb_v15_1M_dmin00.root",
 "/lio/lfs/cms/store/user/qixu/centrality/Glauber/pPb/systematics/Phob_Glau_pPb_sNN70mb_v15_1M_dmin08.root",
 "/lio/lfs/cms/store/user/qixu/centrality/Glauber/pPb/systematics/Phob_Glau_pPb_sNN70mb_v15_1M_R649.root",
 "/lio/lfs/cms/store/user/qixu/centrality/Glauber/pPb/systematics/Phob_Glau_pPb_sNN70mb_v15_1M_R675.root",
 "/lio/lfs/cms/store/user/qixu/centrality/Glauber/pPb/systematics/Phob_Glau_pPb_sNN65mb_v15_1M_dmin04.root",
 "/lio/lfs/cms/store/user/qixu/centrality/Glauber/pPb/systematics/Phob_Glau_pPb_sNN75mb_v15_1M_dmin04.root"
};

TString Gri055Glaulist[nGlau]={
 "/lio/lfs/cms/store/user/qixu/centrality/Glauber/pPb/systematics/glaubergribovpPbv1omega055_default.root",
 "/lio/lfs/cms/store/user/qixu/centrality/Glauber/pPb/systematics/glaubergribovpPbv1omega055_D04914.root",
 "/lio/lfs/cms/store/user/qixu/centrality/Glauber/pPb/systematics/glaubergribovpPbv1omega055_D06006.root",
 "/lio/lfs/cms/store/user/qixu/centrality/Glauber/pPb/systematics/glaubergribovpPbv1omega055_dmin00.root",
 "/lio/lfs/cms/store/user/qixu/centrality/Glauber/pPb/systematics/glaubergribovpPbv1omega055_dmin08.root",
 "/lio/lfs/cms/store/user/qixu/centrality/Glauber/pPb/systematics/glaubergribovpPbv1omega055_R649.root",
 "/lio/lfs/cms/store/user/qixu/centrality/Glauber/pPb/systematics/glaubergribovpPbv1omega055_R675.root",
 "/lio/lfs/cms/store/user/qixu/centrality/Glauber/pPb/systematics/glaubergribovpPbv1omega055_sigmaNN65.root",
 "/lio/lfs/cms/store/user/qixu/centrality/Glauber/pPb/systematics/glaubergribovpPbv1omega055_sigma75.root"
};


TString Gri101Glaulist[nGlau]={
 "/lio/lfs/cms/store/user/qixu/centrality/Glauber/pPb/systematics/glaubergribovpPbv1omega101_default.root",
 "/lio/lfs/cms/store/user/qixu/centrality/Glauber/pPb/systematics/glaubergribovpPbv1omega101_D04914.root",
 "/lio/lfs/cms/store/user/qixu/centrality/Glauber/pPb/systematics/glaubergribovpPbv1omega101_D06006.root",
 "/lio/lfs/cms/store/user/qixu/centrality/Glauber/pPb/systematics/glaubergribovpPbv1omega101_dmin00.root",
 "/lio/lfs/cms/store/user/qixu/centrality/Glauber/pPb/systematics/glaubergribovpPbv1omega101_dmin08.root",
 "/lio/lfs/cms/store/user/qixu/centrality/Glauber/pPb/systematics/glaubergribovpPbv1omega101_R649.root",
 "/lio/lfs/cms/store/user/qixu/centrality/Glauber/pPb/systematics/glaubergribovpPbv1omega101_R675.root",
 "/lio/lfs/cms/store/user/qixu/centrality/Glauber/pPb/systematics/glaubergribovpPbv1omega101_sigmaNN65.root",
 "/lio/lfs/cms/store/user/qixu/centrality/Glauber/pPb/systematics/glaubergribovpPbv1omega101_sigma75.root"
};
