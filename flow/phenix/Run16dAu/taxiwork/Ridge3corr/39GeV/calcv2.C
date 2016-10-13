void calcv2(){
    const int ncent = 6;
    const int npt = 8;
    int pporper = 1; // 0 is pp, 1 is per
    int centbin[ncent+1] = {0,5,10,20,40,60,100};
    gStyle->SetOptFit(kFALSE);
    gStyle->SetOptStat(kFALSE);

    ifstream fcntbbc("c1_c2_central_ptfiner_south.dat");
    ifstream fcntbbcIn("c1_c2_central_ptIn_south.dat");
    ifstream fcntfvtx("c1_c2_central_ptfiner_north.dat");
    ifstream fcntfvtxIn("c1_c2_central_ptIn_north.dat");
    ifstream fbbcfvtx("c1_c2_central_ptIn_sn.dat");

    float c1cntbbc[ncent][npt], c2cntbbc[ncent][npt], c3cntbbc[ncent][npt];
    float c1cntbbcIn[ncent][npt], c2cntbbcIn[ncent][npt], c3cntbbcIn[ncent][npt];
    float c1cntfvtx[ncent][npt], c2cntfvtx[ncent][npt], c3cntfvtx[ncent][npt];
    float c1cntfvtxIn[ncent][npt], c2cntfvtxIn[ncent][npt], c3cntfvtxIn[ncent][npt];
    float c1bbcfvtx[ncent][npt], c2bbcfvtx[ncent][npt], c3bbcfvtx[ncent][npt];
    float c1errcntbbc[ncent][npt], c2errcntbbc[ncent][npt], c3errcntbbc[ncent][npt];
    float c1errcntbbcIn[ncent][npt], c2errcntbbcIn[ncent][npt], c3errcntbbcIn[ncent][npt];
    float c1errcntfvtx[ncent][npt], c2errcntfvtx[ncent][npt], c3errcntfvtx[ncent][npt];
    float c1errcntfvtxIn[ncent][npt], c2errcntfvtxIn[ncent][npt], c3errcntfvtxIn[ncent][npt];
    float c1errbbcfvtx[ncent][npt], c2errbbcfvtx[ncent][npt], c3errbbcfvtx[ncent][npt];
    if(pporper ==0 ){
    ifstream PPfcntbbc("c1_c2_PP_centIn_south.dat");
    ifstream PPfcntbbcIn("c1_c2_PP_ptIn_south.dat");
    ifstream PPfcntfvtx("c1_c2_PP_centIn_north.dat");
    ifstream PPfcntfvtxIn("c1_c2_PP_ptIn_north.dat");
    ifstream PPfbbcfvtx("c1_c2_bbcfvtx_PP_centIn.dat");
    }
    else{
    ifstream PPfcntbbc("c1_c2_central_per_south.dat");
    ifstream PPfcntbbcIn("c1_c2_central_per_ptIn_south.dat");
    ifstream PPfcntfvtx("c1_c2_central_per_north.dat");
    ifstream PPfcntfvtxIn("c1_c2_central_per_ptIn_north.dat");
    ifstream PPfbbcfvtx("c1_c2_central_per_ptIn_sn.dat");
    }
    float c1PPcntbbc[ncent][npt], c2PPcntbbc[ncent][npt], c3PPcntbbc[ncent][npt];
    float c1PPcntbbcIn[ncent][npt], c2PPcntbbcIn[ncent][npt], c3PPcntbbcIn[ncent][npt];
    float c1PPcntfvtx[ncent][npt], c2PPcntfvtx[ncent][npt], c3PPcntfvtx[ncent][npt];
    float c1PPcntfvtxIn[ncent][npt], c2PPcntfvtxIn[ncent][npt], c3PPcntfvtxIn[ncent][npt];
    float c1PPbbcfvtx[ncent][npt], c2PPbbcfvtx[ncent][npt], c3PPbbcfvtx[ncent][npt];
    float c1PPerrcntbbc[ncent][npt], c2PPerrcntbbc[ncent][npt], c3PPerrcntbbc[ncent][npt];
    float c1PPerrcntbbcIn[ncent][npt], c2PPerrcntbbcIn[ncent][npt], c3PPerrcntbbcIn[ncent][npt];
    float c1PPerrcntfvtx[ncent][npt], c2PPerrcntfvtx[ncent][npt], c3PPerrcntfvtx[ncent][npt];
    float c1PPerrcntfvtxIn[ncent][npt], c2PPerrcntfvtxIn[ncent][npt], c3PPerrcntfvtxIn[ncent][npt];
    float c1PPerrbbcfvtx[ncent][npt], c2PPerrbbcfvtx[ncent][npt], c3PPerrbbcfvtx[ncent][npt];

    float tmp;
    int scale = 1; //0: use c1 as scale factor 1: use multiplicity as scale factor

//--------------Multiplicity-----------------------
//float MPPbbc[ncent] = {7.43,7.43,7.43,7.43,7.43,7.43};
    if(pporper ==0 ){
float MPPbbc[ncent] = {3.8,3.8,3.8,3.8,3.8,3.8};
float MPPfvtx[ncent] = {0.69,0.69,0.69,0.69,0.69,0.69};
    }
    else{
float MPPbbc[ncent] = {4.21,4.21,4.21,4.21,4.21,4.21};
float MPPfvtx[ncent] = {2.34,2.34,2.34,2.34,2.34,2.34};
    }
float Mbbc[ncent] = {58.21,44.28,35.55,24.85,8.87,4.16};
float Mfvtx[ncent] = {8.58,7.51,6.71,5.58,3.52,2.80};
//MPPNpart[0] = 2;
//MNpart[0] = 11;
    if(pporper ==0 ){
float MPPNpart[ncent] = {3.8,3.8,3.8,3.8,3.8,3.8};
float MNpart[ncent] = {63.93,47.41,38.03,25.73,13.63,4.21};
    }
    else{
float MPPNpart[ncent] = {4.21,4.21,4.21,4.21,4.21,4.21};
float MNpart[ncent] = {63.93,47.41,38.03,25.73,13.63,4.21};
    }

//--------------read in parameters--------------------------------------------
for(int icent=0;icent<ncent;icent++){
    //output txt
    ofstream of1(Form("v2cntbbcs_cent%d.dat",icent));
    ofstream of2(Form("v2cntfvtxs_cent%d.dat",icent));
    ofstream of(Form("v2_cent%d.dat",icent));
    ofstream ofsub1(Form("v2cntbbcs_cent%d_scale%d.dat",icent,scale));
    ofstream ofsub2(Form("v2cntfvtxs_cent%d_scale%d.dat",icent,scale));
    ofstream ofsub(Form("v2_cent%d_scale%d.dat",icent,scale));

      fcntbbcIn>>c1cntbbcIn[icent][0]>>c1errcntbbcIn[icent][0]>>c2cntbbcIn[icent][0]>>c2errcntbbcIn[icent][0]>>c3cntbbcIn[icent][0]>>c3errcntbbcIn[icent][0];
      fcntfvtxIn>>c1cntfvtxIn[icent][0]>>c1errcntfvtxIn[icent][0]>>c2cntfvtxIn[icent][0]>>c2errcntfvtxIn[icent][0]>>c3cntfvtxIn[icent][0]>>c3errcntfvtxIn[icent][0];
      fbbcfvtx>>c1bbcfvtx[icent][0]>>c1errbbcfvtx[icent][0]>>c2bbcfvtx[icent][0]>>c2errbbcfvtx[icent][0]>>c3bbcfvtx[icent][0]>>c3errbbcfvtx[icent][0];

    PPfcntbbcIn>>c1PPcntbbcIn[icent][0]>>c1PPerrcntbbcIn[icent][0]>>c2PPcntbbcIn[icent][0]>>c2PPerrcntbbcIn[icent][0]>>c3PPcntbbcIn[icent][0]>>c3PPerrcntbbcIn[icent][0];
    PPfcntfvtxIn>>c1PPcntfvtxIn[icent][0]>>c1PPerrcntfvtxIn[icent][0]>>c2PPcntfvtxIn[icent][0]>>c2PPerrcntfvtxIn[icent][0]>>c3PPcntfvtxIn[icent][0]>>c3PPerrcntfvtxIn[icent][0];
    PPfbbcfvtx>>c1PPbbcfvtx[icent][0]>>c1PPerrbbcfvtx[icent][0]>>c2PPbbcfvtx[icent][0]>>c2PPerrbbcfvtx[icent][0]>>c3PPbbcfvtx[icent][0]>>c3PPerrbbcfvtx[icent][0];
    //    fcntbbcIn>>tmp>>tmp>>c1cntbbcIn[icent][0]>>c1errcntbbcIn[icent][0]>>c2cntbbcIn[icent][0]>>c2errcntbbcIn[icent][0];
    //    fcntfvtxIn>>tmp>>tmp>>c1cntfvtxIn[icent][0]>>c1errcntfvtxIn[icent][0]>>c2cntfvtxIn[icent][0]>>c2errcntfvtxIn[icent][0];
    //    fbbcfvtx>>tmp>>tmp>>c1bbcfvtx[icent][0]>>c1errbbcfvtx[icent][0]>>c2bbcfvtx[icent][0]>>c2errbbcfvtx[icent][0];
    
for(int ipt=0;ipt<npt;ipt++){
        fcntbbc>>c1cntbbc[icent][ipt]>>c1errcntbbc[icent][ipt]>>c2cntbbc[icent][ipt]>>c2errcntbbc[icent][ipt]>>c3cntbbc[icent][ipt]>>c3errcntbbc[icent][ipt];
        fcntfvtx>>c1cntfvtx[icent][ipt]>>c1errcntfvtx[icent][ipt]>>c2cntfvtx[icent][ipt]>>c2errcntfvtx[icent][ipt]>>c3cntfvtx[icent][ipt]>>c3errcntfvtx[icent][ipt];
        PPfcntbbc>>c1PPcntbbc[icent][ipt]>>c1PPerrcntbbc[icent][ipt]>>c2PPcntbbc[icent][ipt]>>c2PPerrcntbbc[icent][ipt]>>c3PPcntbbc[icent][ipt]>>c3PPerrcntbbc[icent][ipt];
        PPfcntfvtx>>c1PPcntfvtx[icent][ipt]>>c1PPerrcntfvtx[icent][ipt]>>c2PPcntfvtx[icent][ipt]>>c2PPerrcntfvtx[icent][ipt]>>c3PPcntfvtx[icent][ipt]>>c3PPerrcntfvtx[icent][ipt];
    }
    if(pporper ==0 ){
//Only PP c2 are fitted!
    const int nptpp = 10;
    float ptppmean[nptpp] = {0.25,0.75,1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75};
    TGraphErrors *grPPcntbbc = new TGraphErrors(nptpp,ptppmean,c2PPcntbbc[icent],0,c2PPerrcntbbc[icent]);
    TGraphErrors *grPPcntfvtx = new TGraphErrors(nptpp,ptppmean,c2PPcntfvtx[icent],0,c2PPerrcntfvtx[icent]);
    TF1 *f1 = new TF1("f1","pol3",0,3.5);
    TF1 *f2 = new TF1("f2","pol3",0,3.5);
    TFitResultPtr r1 = grPPcntbbc->Fit("f1", "S");
    TFitResultPtr r2 = grPPcntfvtx->Fit("f2", "S");
    double x[npt]={0.3,0.7,1.2,1.6,2.2,2.8,3.5,4.5};
    double ci1[npt],ci2[npt];
    double cl = 0.683;  // for 1 sigma error
    r1->GetConfidenceIntervals(npt,1,1,x,ci1,cl);
    r2->GetConfidenceIntervals(npt,1,1,x,ci2,cl);
    for(int ipt=0;ipt<npt;ipt++){
        c2PPcntbbc[icent][ipt] = f1->Eval(x[ipt]);
        c2PPerrcntbbc[icent][ipt] = ci1[ipt];
        c2PPcntfvtx[icent][ipt] = f2->Eval(x[ipt]);
        c2PPerrcntfvtx[icent][ipt] = ci2[ipt];
        }
    }

        cout<<"icent: "<<centbin[icent] << "\% to " << centbin[icent+1]<< "\%"<<endl;
        cout<<"cnt - bbc "<<c2cntbbcIn[icent][0]<<" "<<c2errcntbbcIn[icent][0]<<endl;
        cout<<"cnt - fvtx "<<c2cntfvtxIn[icent][0]<<" "<<c2errcntfvtxIn[icent][0]<<endl;
        cout<<"bbc - fvtx "<<c2bbcfvtx[icent][0]<<" "<<c2errbbcfvtx[icent][0]<<endl;

        cout<<"cnt - bbc PP"<<c2PPcntbbcIn[icent][0]<<" "<<c2PPerrcntbbcIn[icent][0]<<endl;
        cout<<"cnt - fvtx PP"<<c2PPcntfvtxIn[icent][0]<<" "<<c2PPerrcntfvtxIn[icent][0]<<endl;
        cout<<"bbc - fvtx PP"<<c2PPbbcfvtx[icent][0]<<" "<<c2PPerrbbcfvtx[icent][0]<<endl;

//---------------calculating cnt/bbc/fvtx inclusive v2------------------------------

        if(c2cntbbcIn[icent][0]*c2cntfvtxIn[icent][0]/c2bbcfvtx[icent][0]<=0) {cout<<"negative?"<<endl; continue;}
        float v2cnt = sqrt(c2cntbbcIn[icent][0]*c2cntfvtxIn[icent][0]/c2bbcfvtx[icent][0]);
        //float v2errcnt = 1/2.*v2cnt*sqrt(TMath::Power(c2errcntbbcIn[icent][0]/c2cntbbcIn[icent][0],2)+TMath::Power(c2errcntfvtxIn[icent][0]/c2cntfvtxIn[icent][0],2)+TMath::Power(c2errbbcfvtx[icent][0]/c2bbcfvtx[icent][0],2));
        float v2errcnt = get3sqerr(c2cntbbcIn[icent][0],c2errcntbbcIn[icent][0],c2cntfvtxIn[icent][0],c2errcntfvtxIn[icent][0],c2bbcfvtx[icent][0],c2errbbcfvtx[icent][0]);

        float v2bbc = c2cntbbcIn[icent][0]/v2cnt;
        //float v2errbbc = v2bbc * sqrt(TMath::Power(c2errcntbbcIn[icent][0]/c2cntbbcIn[icent][0],2)+TMath::Power(v2errcnt/v2cnt,2));
        float v2errbbc = get2derr(c2cntbbcIn[icent][0],c2errcntbbcIn[icent][0],v2cnt,v2errcnt);
        
        float v2fvtx = c2cntfvtxIn[icent][0]/v2cnt;
        //float v2errfvtx = v2fvtx * sqrt(TMath::Power(c2errcntfvtxIn[icent][0]/c2cntfvtxIn[icent][0],2)+TMath::Power(v2errcnt/v2cnt,2));
        float v2errfvtx = get2derr(c2cntfvtxIn[icent][0],c2errcntfvtxIn[icent][0],v2cnt,v2errcnt);

        cout<<"Integral v2 cnt = "<<v2cnt<<" "<<v2errcnt<<endl;
        cout<<"Integral v2 bbc = "<<v2bbc<<" "<<v2errbbc<<endl;
        cout<<"Integral v2 fvtx = "<<v2fvtx<<" "<<v2errfvtx<<endl;

//------------calculating cnt pt dependence v2--------------------method1----------------
        
        float v2cntpt1[ncent][npt], v2errcntpt1[ncent][npt];
        float v2cntpt2[ncent][npt], v2errcntpt2[ncent][npt];
        float ptmean[npt]={0.3,0.7,1.2,1.6,2.2,2.8,3.5,4.5};
    for(int ipt=0;ipt<npt;ipt++){
        v2cntpt1[icent][ipt] = c2cntbbc[icent][ipt]/v2bbc;
        //v2errcntpt1[icent][ipt] = v2cntpt1[icent][ipt]*sqrt(TMath::Power(c2errcntbbc[icent][ipt]/c2cntbbc[icent][ipt],2)+TMath::Power(v2errbbc/v2bbc,2));
        v2errcntpt1[icent][ipt] = get2derr(c2cntbbc[icent][ipt],c2errcntbbc[icent][ipt],v2bbc,v2errbbc);

        v2cntpt2[icent][ipt] = c2cntfvtx[icent][ipt]/v2fvtx;
        //v2errcntpt2[icent][ipt] = v2cntpt2[icent][ipt]*sqrt(TMath::Power(c2errcntfvtx[icent][ipt]/c2cntfvtx[icent][ipt],2)+TMath::Power(v2errfvtx/v2fvtx,2));
        v2errcntpt2[icent][ipt] = get2derr(c2cntfvtx[icent][ipt],c2errcntfvtx[icent][ipt],v2bbc,v2errbbc);
    }

//------------calculating cnt pt dependence v2--------------------method2----------------
        
        float v2cntpt[ncent][npt];
        float v2errcntpt[ncent][npt];
    for(int ipt=0;ipt<npt;ipt++){
        v2cntpt[icent][ipt] = sqrt(c2cntbbc[icent][ipt]*c2cntfvtx[icent][ipt]/c2bbcfvtx[icent][0]);
        //v2errcntpt[icent][ipt] = v2cntpt[icent][ipt]*sqrt(TMath::Power(c2errcntbbc[icent][ipt]/c2cntbbc[icent][ipt],2)+TMath::Power(c2errcntfvtx[icent][ipt]/c2cntfvtx[icent][ipt],2)+TMath::Power(c2errbbcfvtx[icent][0]/c2bbcfvtx[icent][0],2));
        v2errcntpt[icent][ipt] = get3sqerr(c2cntbbc[icent][ipt],c2errcntbbc[icent][ipt],c2cntfvtx[icent][ipt],c2errcntfvtx[icent][ipt],c2bbcfvtx[icent][0],c2errbbcfvtx[icent][0]);

        of1<<ptmean[ipt]<<" "<<v2cntpt1[icent][ipt]<<" "<<v2errcntpt1[icent][ipt]<<endl;
        of2<<ptmean[ipt]<<" "<<v2cntpt2[icent][ipt]<<" "<<v2errcntpt2[icent][ipt]<<endl;
        of<<ptmean[ipt]<<" "<<v2cntpt[icent][ipt]<<" "<<v2errcntpt[icent][ipt]<<endl;
    }
    of1.close();
    of2.close();
    of.close();
        
//----subtract---------calculating cnt/bbc/fvtx inclusive v2------------------------------

    if(scale==0){
        c2cntbbcIn[icent][0] = c2cntbbcIn[icent][0] - c2PPcntbbcIn[icent][0]*c1cntbbcIn[icent][0]/c1PPcntbbcIn[icent][0];
        c2cntfvtxIn[icent][0] = c2cntfvtxIn[icent][0] - c2PPcntfvtxIn[icent][0]*c1cntfvtxIn[icent][0]/c1PPcntfvtxIn[icent][0];
        c2bbcfvtx[icent][0] = c2bbcfvtx[icent][0] - c2PPbbcfvtx[icent][0]*c1bbcfvtx[icent][0]/c1PPbbcfvtx[icent][0];
        c2errcntbbcIn[icent][0] = get2mierr(c2cntbbcIn[icent][0],c2errcntbbcIn[icent][0],c2PPcntbbcIn[icent][0]*c1cntbbcIn[icent][0]/c1PPcntbbcIn[icent][0], get3err(c2PPcntbbcIn[icent][0],c2PPerrcntbbcIn[icent][0],c1cntbbcIn[icent][0],c1errcntbbcIn[icent][0],c1PPcntbbcIn[icent][0],c1PPerrcntbbcIn[icent][0]));
        c2errcntfvtxIn[icent][0] = get2mierr(c2cntfvtxIn[icent][0],c2errcntfvtxIn[icent][0],c2PPcntfvtxIn[icent][0]*c1cntfvtxIn[icent][0]/c1PPcntfvtxIn[icent][0], get3err(c2PPcntfvtxIn[icent][0],c2PPerrcntfvtxIn[icent][0],c1cntfvtxIn[icent][0],c1errcntfvtxIn[icent][0],c1PPcntfvtxIn[icent][0],c1PPerrcntfvtxIn[icent][0]));
        c2errbbcfvtx[icent][0] = get2mierr(c2bbcfvtx[icent][0],c2errbbcfvtx[icent][0],c2PPbbcfvtx[icent][0]*c1bbcfvtx[icent][0]/c1PPbbcfvtx[icent][0], get3err(c2PPbbcfvtx[icent][0],c2PPerrbbcfvtx[icent][0],c1bbcfvtx[icent][0],c1errbbcfvtx[icent][0],c1PPbbcfvtx[icent][0],c1PPerrbbcfvtx[icent][0]));
    }
    else{
        c2cntbbcIn[icent][0] = c2cntbbcIn[icent][0] - c2PPcntbbcIn[icent][0]*MPPbbc[icent]/Mbbc[icent];
        c2cntfvtxIn[icent][0] = c2cntfvtxIn[icent][0] - c2PPcntfvtxIn[icent][0]*MPPfvtx[icent]/Mfvtx[icent];
        c2bbcfvtx[icent][0] = c2bbcfvtx[icent][0] - c2PPbbcfvtx[icent][0]*(MPPbbc[icent]*MPPfvtx[icent]/MPPNpart[icent])/(Mbbc[icent]*Mfvtx[icent]/MNpart[icent]);
        c2errcntbbcIn[icent][0] = get2mierr(c2cntbbcIn[icent][0],c2errcntbbcIn[icent][0],c2PPcntbbcIn[icent][0]*MPPbbc[icent]/Mbbc[icent],c2PPerrcntbbcIn[icent][0]*MPPbbc[icent]/Mbbc[icent]);
        c2errcntfvtxIn[icent][0] = get2mierr(c2cntfvtxIn[icent][0],c2errcntfvtxIn[icent][0],c2PPcntfvtxIn[icent][0]*MPPfvtx[icent]/Mfvtx[icent],c2PPerrcntfvtxIn[icent][0]*MPPfvtx[icent]/Mfvtx[icent]);
        c2errbbcfvtx[icent][0] = get2mierr(c2bbcfvtx[icent][0],c2errbbcfvtx[icent][0],c2PPbbcfvtx[icent][0]*(MPPbbc[icent]*MPPfvtx[icent]/MPPNpart[icent])/(Mbbc[icent]*Mfvtx[icent]/MNpart[icent]),c2PPerrbbcfvtx[icent][0]*(MPPbbc[icent]*MPPfvtx[icent]/MPPNpart[icent])/(Mbbc[icent]*Mfvtx[icent]/MNpart[icent]));
    }
        cout<<"cnt - bbc subtract ="<<c2cntbbcIn[icent][0]<<" "<<c2errcntbbcIn[icent][0]<<endl;
        cout<<"cnt - fvtx subtract ="<<c2cntfvtxIn[icent][0]<<" "<<c2errcntfvtxIn[icent][0]<<endl;
        cout<<"bbc - fvtx subtract ="<<c2bbcfvtx[icent][0]<<" "<<c2errbbcfvtx[icent][0]<<endl;

        float v2cnt = sqrt(c2cntbbcIn[icent][0]*c2cntfvtxIn[icent][0]/c2bbcfvtx[icent][0]);
      //  float v2errcnt = 1/2.*v2cnt*sqrt(TMath::Power(c2errcntbbcIn[icent][0]/c2cntbbcIn[icent][0],2)+TMath::Power(c2errcntfvtxIn[icent][0]/c2cntfvtxIn[icent][0],2)+TMath::Power(c2errbbcfvtx[icent][0]/c2bbcfvtx[icent][0],2));
        float v2errcnt = get3sqerr(c2cntbbcIn[icent][0],c2errcntbbcIn[icent][0],c2cntfvtxIn[icent][0],c2errcntfvtxIn[icent][0],c2bbcfvtx[icent][0],c2errbbcfvtx[icent][0]);

        float v2bbc = c2cntbbcIn[icent][0]/v2cnt;
        //float v2errbbc = v2bbc * sqrt(TMath::Power(c2errcntbbcIn[icent][0]/c2cntbbcIn[icent][0],2)+TMath::Power(v2errcnt/v2cnt,2));
        float v2errbbc = get2derr(c2cntbbcIn[icent][0],c2errcntbbcIn[icent][0],v2cnt,v2errcnt);
        
        float v2fvtx = c2cntfvtxIn[icent][0]/v2cnt;
        //float v2errfvtx = v2fvtx * sqrt(TMath::Power(c2errcntfvtxIn[icent][0]/c2cntfvtxIn[icent][0],2)+TMath::Power(v2errcnt/v2cnt,2));
        float v2errfvtx = get2derr(c2cntfvtxIn[icent][0],c2errcntfvtxIn[icent][0],v2cnt,v2errcnt);

        cout<<"Integral v2 subtract cnt = "<<v2cnt<<" "<<v2errcnt<<endl;
        cout<<"Integral v2 subtract bbc = "<<v2bbc<<" "<<v2errbbc<<endl;
        cout<<"Integral v2 subtract fvtx = "<<v2fvtx<<" "<<v2errfvtx<<endl;

//----subtract--------calculating cnt pt dependence v2--------------------method1----------------
        
        float v2cntpt1[ncent][npt], v2errcntpt1[ncent][npt];
        float v2cntpt2[ncent][npt], v2errcntpt2[ncent][npt];
        float ptmean[npt]={0.3,0.7,1.2,1.6,2.2,2.8,3.5,4.5};
    for(int ipt=0;ipt<npt;ipt++){
    if(scale==0){
        c2cntbbc[icent][ipt] = c2cntbbc[icent][ipt] - c2PPcntbbc[icent][ipt]*c1cntbbc[icent][ipt]/c1PPcntbbc[icent][ipt];
        c2cntfvtx[icent][ipt] = c2cntfvtx[icent][ipt] - c2PPcntfvtx[icent][ipt]*c1cntfvtx[icent][ipt]/c1PPcntfvtx[icent][ipt];
        c2bbcfvtx[icent][ipt] = c2bbcfvtx[icent][ipt] - c2PPbbcfvtx[icent][ipt]*c1bbcfvtx[icent][ipt]/c1PPbbcfvtx[icent][ipt];
        c2errcntbbc[icent][ipt] = get2mierr(c2cntbbc[icent][ipt],c2errcntbbc[icent][ipt],c2PPcntbbc[icent][ipt]*c1cntbbc[icent][ipt]/c1PPcntbbc[icent][ipt], get3err(c2PPcntbbc[icent][ipt],c2PPerrcntbbc[icent][ipt],c1cntbbc[icent][ipt],c1errcntbbc[icent][ipt],c1PPcntbbc[icent][ipt],c1PPerrcntbbc[icent][ipt]));
        c2errcntfvtx[icent][ipt] = get2mierr(c2cntfvtx[icent][ipt],c2errcntfvtx[icent][ipt],c2PPcntfvtx[icent][ipt]*c1cntfvtx[icent][ipt]/c1PPcntfvtx[icent][ipt], get3err(c2PPcntfvtx[icent][ipt],c2PPerrcntfvtx[icent][ipt],c1cntfvtx[icent][ipt],c1errcntfvtx[icent][ipt],c1PPcntfvtx[icent][ipt],c1PPerrcntfvtx[icent][ipt]));
        c2errbbcfvtx[icent][ipt] = get2mierr(c2bbcfvtx[icent][ipt],c2errbbcfvtx[icent][ipt],c2PPbbcfvtx[icent][ipt]*c1bbcfvtx[icent][ipt]/c1PPbbcfvtx[icent][ipt], get3err(c2PPbbcfvtx[icent][ipt],c2PPerrbbcfvtx[icent][ipt],c1bbcfvtx[icent][ipt],c1errbbcfvtx[icent][ipt],c1PPbbcfvtx[icent][ipt],c1PPerrbbcfvtx[icent][ipt]));
    }
    else{
        c2cntbbc[icent][ipt] = c2cntbbc[icent][ipt] - c2PPcntbbc[icent][ipt]*MPPbbc[icent]/Mbbc[icent];
        c2cntfvtx[icent][ipt] = c2cntfvtx[icent][ipt] - c2PPcntfvtx[icent][ipt]*MPPfvtx[icent]/Mfvtx[icent];
      //  c2bbcfvtx[icent][ipt] = c2bbcfvtx[icent][ipt] - c2PPbbcfvtx[icent][ipt]*(MPPbbc[icent]*MPPfvtx[icent]/MPPNpart[icent])/(Mbbc[icent]*Mfvtx[icent]/MNpart[icent]);
if(icent==0&&ipt>=2 &&ipt<=3){
        cout<<c2cntbbc[icent][ipt]<<" "<<c2errcntbbc[icent][ipt]<<" "<<c2PPcntbbc[icent][ipt]*MPPbbc[icent]/Mbbc[icent]<<" "<<c2PPerrcntbbc[icent][ipt]*MPPbbc[icent]/Mbbc[icent]<<endl;
        cout<<c2cntfvtx[icent][ipt]<<" "<<c2errcntfvtx[icent][ipt]<<" "<<c2PPcntfvtx[icent][ipt]*MPPfvtx[icent]/Mfvtx[icent]<<endl;
        cout<<c2bbcfvtx[icent][0]<<" "<<c2errbbcfvtx[icent][0]<<" "<<c2PPbbcfvtx[icent][0]*MPPfvtx[icent]/Mfvtx[icent]<<endl;
        cout<<sqrt(c2cntbbc[icent][ipt]*c2cntfvtx[icent][ipt]/c2bbcfvtx[icent][0])<<endl;
        cout<<get3sqerr(c2cntbbc[icent][ipt],c2errcntbbc[icent][ipt],c2cntfvtx[icent][ipt],c2errcntfvtx[icent][ipt],c2bbcfvtx[icent][0],c2errbbcfvtx[icent][0])<<endl;;
}
        c2errcntbbc[icent][ipt] = get2mierr(c2cntbbc[icent][ipt],c2errcntbbc[icent][ipt],c2PPcntbbc[icent][ipt]*MPPbbc[icent]/Mbbc[icent],c2PPerrcntbbc[icent][ipt]*MPPbbc[icent]/Mbbc[icent]);
        c2errcntfvtx[icent][ipt] = get2mierr(c2cntfvtx[icent][ipt],c2errcntfvtx[icent][ipt],c2PPcntfvtx[icent][ipt]*MPPfvtx[icent]/Mfvtx[icent],c2PPerrcntfvtx[icent][ipt]*MPPfvtx[icent]/Mfvtx[icent]);
        //c2errbbcfvtx[icent][ipt] = get2mierr(c2bbcfvtx[icent][ipt],c2errbbcfvtx[icent][ipt],c2PPbbcfvtx[icent][ipt]*(MPPbbc[icent]*MPPfvtx[icent]/MPPNpart[icent])/(Mbbc[icent]*Mfvtx[icent]/MNpart[icent]),c2PPerrbbcfvtx[icent][ipt]*(MPPbbc[icent]*MPPfvtx[icent]/MPPNpart[icent])/(Mbbc[icent]*Mfvtx[icent]/MNpart[icent]));
if(icent==0&&ipt>=2&&ipt<=3){
        cout<<c2cntbbc[icent][ipt]<<" "<<c2errcntbbc[icent][ipt]<<" "<<c2PPcntbbc[icent][ipt]*MPPbbc[icent]/Mbbc[icent]<<" "<<c2PPerrcntbbc[icent][ipt]*MPPbbc[icent]/Mbbc[icent]<<endl;
        cout<<c2cntfvtx[icent][ipt]<<" "<<c2errcntfvtx[icent][ipt]<<" "<<c2PPcntfvtx[icent][ipt]*MPPfvtx[icent]/Mfvtx[icent]<<endl;
        cout<<c2bbcfvtx[icent][0]<<" "<<c2errbbcfvtx[icent][0]<<" "<<c2PPbbcfvtx[icent][0]*MPPfvtx[icent]/Mfvtx[icent]<<endl;
        cout<<sqrt(c2cntbbc[icent][ipt]*c2cntfvtx[icent][ipt]/c2bbcfvtx[icent][0])<<endl;
        cout<<get3sqerr(c2cntbbc[icent][ipt],c2errcntbbc[icent][ipt],c2cntfvtx[icent][ipt],c2errcntfvtx[icent][ipt],c2bbcfvtx[icent][0],c2errbbcfvtx[icent][0])<<endl;;
}
    }
    }
    for(int ipt=0;ipt<npt;ipt++){
        v2cntpt1[icent][ipt] = c2cntbbc[icent][ipt]/v2bbc;
        //v2errcntpt1[icent][ipt] = v2cntpt1[icent][ipt]*sqrt(TMath::Power(c2errcntbbc[icent][ipt]/c2cntbbc[icent][ipt],2)+TMath::Power(v2errbbc/v2bbc,2));
        v2errcntpt1[icent][ipt] = get2derr(c2cntbbc[icent][ipt],c2errcntbbc[icent][ipt],v2bbc,v2errbbc);

        v2cntpt2[icent][ipt] = c2cntfvtx[icent][ipt]/v2fvtx;
        //v2errcntpt2[icent][ipt] = v2cntpt2[icent][ipt]*sqrt(TMath::Power(c2errcntfvtx[icent][ipt]/c2cntfvtx[icent][ipt],2)+TMath::Power(v2errfvtx/v2fvtx,2));
        v2errcntpt2[icent][ipt] = get2derr(c2cntfvtx[icent][ipt],c2errcntfvtx[icent][ipt],v2fvtx,v2errfvtx);
    }

//----subtract--------calculating cnt pt dependence v2--------------------method2----------------
        
        float v2cntpt[ncent][npt];
        float v2errcntpt[ncent][npt];
    for(int ipt=0;ipt<npt;ipt++){
        v2cntpt[icent][ipt] = sqrt(c2cntbbc[icent][ipt]*c2cntfvtx[icent][ipt]/c2bbcfvtx[icent][0]);
        //v2errcntpt[icent][ipt] = v2cntpt[icent][ipt]*sqrt(TMath::Power(c2errcntbbc[icent][ipt]/c2cntbbc[icent][ipt],2)+TMath::Power(c2errcntfvtx[icent][ipt]/c2cntfvtx[icent][ipt],2)+TMath::Power(c2errbbcfvtx[icent][0]/c2bbcfvtx[icent][0],2));
        v2errcntpt[icent][ipt] = get3sqerr(c2cntbbc[icent][ipt],c2errcntbbc[icent][ipt],c2cntfvtx[icent][ipt],c2errcntfvtx[icent][ipt],c2bbcfvtx[icent][0],c2errbbcfvtx[icent][0]);
        
        ofsub1<<ptmean[ipt]<<" "<<v2cntpt1[icent][ipt]<<" "<<v2errcntpt1[icent][ipt]<<endl;
        ofsub2<<ptmean[ipt]<<" "<<v2cntpt2[icent][ipt]<<" "<<v2errcntpt2[icent][ipt]<<endl;
        ofsub<<ptmean[ipt]<<" "<<v2cntpt[icent][ipt]<<" "<<v2errcntpt[icent][ipt]<<endl;
    }
    ofsub1.close();
    ofsub2.close();
    ofsub.close();
    
if(pporper == 0 ){
    ifstream PPfcntbbc("c1_c2_PP_centIn_south.dat");
    ifstream PPfcntbbcIn("c1_c2_PP_ptIn_south.dat");
    ifstream PPfcntfvtx("c1_c2_PP_centIn_north.dat");
    ifstream PPfcntfvtxIn("c1_c2_PP_ptIn_north.dat");
    ifstream PPfbbcfvtx("c1_c2_bbcfvtx_PP_centIn.dat");
}
else{
    ifstream PPfcntbbc("c1_c2_central_per_south.dat");
    ifstream PPfcntbbcIn("c1_c2_central_per_ptIn_south.dat");
    ifstream PPfcntfvtx("c1_c2_central_per_north.dat");
    ifstream PPfcntfvtxIn("c1_c2_central_per_ptIn_north.dat");
    ifstream PPfbbcfvtx("c1_c2_central_per_ptIn_sn.dat");
    }
}
}

float get3sqerr(float a, float ea, float b, float eb, float c, float ec){
    float d = sqrt(a*b/c);
    if(a==0 || b==0 || c==0) return 0;
    float ed = 0.5 * d * sqrt(ea/a*ea/a+eb/b*eb/b+ec/c*ec/c);
    return ed;
}

float get3err(float a, float ea, float b, float eb, float c, float ec){
    float d = a*b/c;
    if(a==0 || b==0 || c==0) return 0;
    float ed = d * sqrt(ea/a*ea/a+eb/b*eb/b+ec/c*ec/c);
    return ed;
}

float get2merr(float a, float ea, float b, float eb){
    float c = a*b;
    if(a==0 || b==0) return 0;
    float ec = c * sqrt(ea/a*ea/a+eb/b*eb/b);
    return ec;
}

float get2derr(float a, float ea, float b, float eb){
    float c = a/b;
    if(a==0 || b==0) return 0;
    float ec = c * sqrt(ea/a*ea/a+eb/b*eb/b);
    return ec;
}

float get2aerr(float a, float ea, float b, float eb){
    float c = a+b;
    if(a==0 || b==0) return 0;
    float ec = c * sqrt(ea/a*ea/a+eb/b*eb/b);
    return ec;
}

float get2mierr(float a, float ea, float b, float eb){
    float c = a-b;
    if(a==0 || b==0) return 0;
    float ec = c * sqrt(ea/a*ea/a+eb/b*eb/b);
    return ec;
}
