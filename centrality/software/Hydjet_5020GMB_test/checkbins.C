#include "iostream"
#include "iomanip"
void checkbins(){
    TString fs = "/cms/store/user/qixu/centrality/software/Hydjet_5020GMB_test/eventtree.root";
    const double HFbins[]={0, 0.59639, 1.29968, 1.97784, 2.75471, 3.45017, 4.25126, 5.11264, 5.90352, 6.68281, 7.50329, 8.34573, 9.28051, 10.1751, 11.0425, 11.9863, 13.0763, 14.1645, 15.2643, 16.4405, 17.476, 18.5394, 19.7623, 21.1482, 22.5326, 23.9879, 25.5085, 26.8893, 28.399, 29.9031, 31.6868, 33.3283, 35.1382, 36.8985, 38.9402, 40.86, 43.0312, 45.2068, 47.444, 49.9495, 52.3375, 54.9166, 57.4942, 60.2328, 63.2836, 66.2615, 69.2931, 72.5929, 76.0682, 79.7105, 83.4692, 87.247, 91.0883, 95.4537, 99.3284, 103.87, 108.218, 112.768, 117.628, 122.395, 127.284, 132.54, 138.527, 144.089, 150.588, 156.714, 162.969, 169.273, 176.176, 183.024, 190.035, 197.182, 204.326, 211.773, 219.763, 227.411, 235.715, 244.734, 253.971, 263.454, 272.179, 281.496, 290.787, 300.438, 310.689, 320.978, 330.659, 341.554, 352.199, 363.815, 375.502, 387.705, 399.861, 412.892, 426.11, 438.371, 450.723, 464.662, 477.878, 491.297, 505.658, 519.494, 535.056, 550.44, 564.311, 578.085, 593.454, 610.702, 628.058, 646.714, 663.501, 681.378, 698.18, 715.534, 733.211, 750.997, 769.845, 790.143, 810.741, 831.265, 851.627, 871.813, 892.845, 916.292, 939.527, 960.58, 983.788, 1009.2, 1032.86, 1058.49, 1082.36, 1107.23, 1131.6, 1158.65, 1185.81, 1210.9, 1238.42, 1265.31, 1293.52, 1320.97, 1349.39, 1376.55, 1407.17, 1435.37, 1465.2, 1496.71, 1525.76, 1557.54, 1591.16, 1623.29, 1657.58, 1693, 1733.09, 1767.15, 1804.91, 1843.26, 1879.61, 1915.43, 1954.76, 1994.7, 2037.45, 2079.48, 2123.4, 2166.54, 2209.46, 2255.51, 2300.83, 2342.17, 2389.45, 2434.96, 2483.89, 2532.47, 2582.19, 2633.2, 2684.11, 2737.95, 2794.76, 2847.79, 2903.33, 2958.13, 3013, 3072.54, 3134.47, 3197.52, 3259.3, 3322.45, 3388.86, 3454.35, 3523.45, 3596.01, 3671.27, 3746.02, 3816.44, 3896.17, 3975.58, 4061.46, 4152.57, 4240.07, 4329.56, 4443.85, 4999.26};
    int nbin=sizeof(HFbins)/sizeof(double)-1;
    ofstream fout("eventsinbins.txt");
    fout<<"nbins = "<<nbin<<endl;
    int count[200]={};
    Float_t hiHFplus;
    Float_t hiHFminus;
    Float_t hiHFminusEta4;
    Float_t hiHFplusEta4;
    Float_t hiHF;
    TFile *f = new TFile(fs,"ReadOnly");
    TTree *tree = (TTree*)f->Get("hiEvtAnalyzer/HiTree");
    tree->SetBranchAddress("hiHFplusEta4",&hiHFplusEta4);
    tree->SetBranchAddress("hiHFplus",&hiHFplus);
    tree->SetBranchAddress("hiHFminusEta4",&hiHFminusEta4);
    tree->SetBranchAddress("hiHFminus",&hiHFminus);
    tree->SetBranchAddress("hiHF",&hiHF);
    for(int i=0;i<tree->GetEntries();i++){
        tree->GetEntry(i);
 //       Float_t par = hiHFplus+hiHFminus;
        Float_t par = hiHF;
        for(int ibin=0;ibin<nbin;ibin++){
            if(par>=HFbins[ibin] && par < HFbins[ibin+1])
                break;
        }
        count[ibin]++;
    }
    int totalevent=0;
    fout<<"event falls into each bin"<<endl;
    fout<<"ibin"<<"\t"<<"nevent"<<endl;
    for(int ibin=0;ibin<nbin;ibin++){
        fout<<ibin<<"\t"<<count[ibin]<<endl;
        totalevent += count[ibin];
    }
    fout<<"totalevent = "<<totalevent<<endl;
}

