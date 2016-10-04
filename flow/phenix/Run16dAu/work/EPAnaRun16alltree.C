#include "EPAnaRun16alltree.h"

#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h" 
#include "TString.h"
#include "TVectorD.h"
#include <fstream>
#include <string>
#include <TROOT.h>
#include "TSystem.h"
#include <iostream>

#include "Run16dAupc3dphidzcalibsmoothpass1.h"
#include "func.h"

using namespace std;
int get_fvtx_layer(float);
void initialize_pmt_position();
float d_pmt_z = -1443.5; // same for all tubes

//_____________________________________________________________________________________________________________________________
EPAnaRun16alltree::EPAnaRun16alltree(std::vector<TString> input, std::vector<TString> input1, const char* output) :
  OutputFileName(output), InputFileName(input), InputFileName1(input1), ievent(0), jevent(0), calFlag(0), RunNumber(0)
{
  d_outfile=NULL;

  for(int icent=0;icent<ncent;icent++){
    for(int ibbcz=0;ibbcz<nbbcz;ibbcz++){
      for(int ihar=0;ihar<nhar;ihar++){
        for(int isub=0;isub<nsub;isub++){
            psi[icent][ibbcz][ihar][isub]=NULL;
            psiFla[icent][ibbcz][ihar][isub]=NULL;
            phiweight[icent][ibbcz][ihar][isub]=NULL;
            phiweightbbc[icent][ibbcz][ihar][isub]=NULL;
          for(int ixy=0;ixy<nxy;ixy++){
            q[icent][ibbcz][ihar][isub][ixy]=NULL;
            qRec[icent][ibbcz][ihar][isub][ixy]=NULL;
          }
          for(int iord=0;iord<nord;iord++){
            sinflt[icent][ibbcz][ihar][isub][iord]=NULL;
            cosflt[icent][ibbcz][ihar][isub][iord]=NULL;
          }
        }
      }
    }
  }
    
  for(int icent=0;icent<ncent;icent++){
      for(int ihar=0;ihar<nhar;ihar++){
        /*
          EPRCNTBBCS[icent][ihar]=NULL;
          EPRCNTFVTX1S[icent][ihar]=NULL;
          EPRCNTFVTX1S[icent][ihar]=NULL;
          EPRCNTFVTX2S[icent][ihar]=NULL;
          EPRCNTFVTX1LS[icent][ihar]=NULL;
          EPRCNTFVTX2LS[icent][ihar]=NULL;
          EPRCNTFVTX3LS[icent][ihar]=NULL;
          EPRCNTFVTX4LS[icent][ihar]=NULL;
          EPRBBCSFVTX1S[icent][ihar]=NULL;
          EPRBBCSFVTX2S[icent][ihar]=NULL;
          */
          EPRBBCSFVTX1S[icent][ihar]=NULL;
          EPRBBCSFVTX1LS[icent][ihar]=NULL;
          EPRBBCSFVTX2LS[icent][ihar]=NULL;
          EPRBBCSFVTX3LS[icent][ihar]=NULL;
          EPRBBCSFVTX4LS[icent][ihar]=NULL;
          EPRBBCSFVTX1p2p3LS[icent][ihar]=NULL;
          EPRBBCSFVTX1p2p4LS[icent][ihar]=NULL;
          EPRFVTX1NFVTX1S[icent][ihar]=NULL;
          EPRBBCSFVTXtrkS[icent][ihar]=NULL;
          /*
          EPRCNTFVTX1trkS[icent][ihar]=NULL;
          EPRCNTFVTX2trkS[icent][ihar]=NULL;

          EPRCNTcBBCSc[icent][ihar]=NULL;
          EPRCNTcFVTX1S[icent][ihar]=NULL;
          EPRBBCScFVTX1S[icent][ihar]=NULL;

          EPRBBCSFVTX1trkS[icent][ihar]=NULL;
          EPRBBCSFVTX2trkS[icent][ihar]=NULL;
          EPRCNTcBBCSc[icent][ihar]=NULL;
          EPRCNTcFVTX1Sc[icent][ihar]=NULL;
          EPRBBCScFVTX1Sc[icent][ihar]=NULL;
          */
          
        for(int iphi=0;iphi<nphi;iphi++){
          /*
          vobsBBCS[icent][ihar][iphi]=NULL;
          vobsFVTX1S[icent][ihar][iphi]=NULL;
          vobsFVTX2S[icent][ihar][iphi]=NULL;
          */
          /*
          vobsBBCSsq[icent][ihar][iphi]=NULL;
          vobsFVTX1Ssq[icent][ihar][iphi]=NULL;
          vobsFVTX2Ssq[icent][ihar][iphi]=NULL;
          */
          vobsFVTX1LS[icent][ihar][iphi]=NULL;
          vobsFVTX2LS[icent][ihar][iphi]=NULL;
          vobsFVTX3LS[icent][ihar][iphi]=NULL;
          vobsFVTX4LS[icent][ihar][iphi]=NULL;
          vobsFVTX1p2p3LS[icent][ihar][iphi]=NULL;
          vobsFVTX1p2p4LS[icent][ihar][iphi]=NULL;
          vobsBBCS[icent][ihar][iphi]=NULL;
          vobsFVTXtrkS[icent][ihar][iphi]=NULL;
          vobsFVTX1N[icent][ihar][iphi]=NULL;
          vobsFVTX1S[icent][ihar][iphi]=NULL;

          vFVTX1LS[icent][ihar][iphi]=NULL;
          vFVTX2LS[icent][ihar][iphi]=NULL;
          vFVTX3LS[icent][ihar][iphi]=NULL;
          vFVTX4LS[icent][ihar][iphi]=NULL;
          vFVTX1p2p3LS[icent][ihar][iphi]=NULL;
          vFVTX1p2p4LS[icent][ihar][iphi]=NULL;
          vBBCS[icent][ihar][iphi]=NULL;
          vFVTXtrkS[icent][ihar][iphi]=NULL;
          vFVTX1N[icent][ihar][iphi]=NULL;
          vFVTX1S[icent][ihar][iphi]=NULL;

          vnFVTX1LS[icent][ihar][iphi]=NULL;
          vnFVTX2LS[icent][ihar][iphi]=NULL;
          vnFVTX3LS[icent][ihar][iphi]=NULL;
          vnFVTX4LS[icent][ihar][iphi]=NULL;
          vnFVTX1p2p3LS[icent][ihar][iphi]=NULL;
          vnFVTX1p2p4LS[icent][ihar][iphi]=NULL;
          vnBBCS[icent][ihar][iphi]=NULL;
          vnFVTXtrkS[icent][ihar][iphi]=NULL;
          vnFVTX1N[icent][ihar][iphi]=NULL;
          vnFVTX1S[icent][ihar][iphi]=NULL;
          /*
          vobsFVTX2S[icent][ihar][iphi]=NULL;
          vobsFVTX2Ssq[icent][ihar][iphi]=NULL;
          vobscBBCSc[icent][ihar][iphi]=NULL;
          vobscBBCScsq[icent][ihar][iphi]=NULL;
          vobscFVTX1Sc[icent][ihar][iphi]=NULL;
          vobscFVTX1Scsq[icent][ihar][iphi]=NULL;
          vobsFVTX1trkS[icent][ihar][iphi]=NULL;
          vobsFVTX2trkS[icent][ihar][iphi]=NULL;
          vobsFVTX1trkSsq[icent][ihar][iphi]=NULL;
          vobsFVTX2trkSsq[icent][ihar][iphi]=NULL;
          */
        }
        }
  }


for(int istation=0;istation<4;istation++){
    hfvtx[istation]=NULL;
}
    hbbc=NULL;
    hnfvtxtrk=NULL;
    hbbcsnfvtxtrk=NULL;
    hfvtxtrk=NULL;
}

//_____________________________________________________________________________________________________________________________
EPAnaRun16alltree::~EPAnaRun16alltree()
{
//  delete bbccalib;
//  delete bbcgeo;
  cout << " EPAnaRun16alltree::~EPAnaRun16alltree " << endl;
}

//_____________________________________________________________________________________________________________________________

int EPAnaRun16alltree::Init()
{
  cout << " EPAnaRun16alltree::Init " << endl;
  /* initialize random seed: */

  initialize_pmt_position();
  ievent = 0;
  jevent = 0;

  cout<<"finish of initialize"<<endl;

  cout << " EPAnaRun16alltree::Init Histos " << endl;
  char name[200];
  float pi = acos(-1.0);
  d_outfile = new TFile(OutputFileName.c_str(),"recreate");

  for(int icent=0;icent<ncent;icent++){
    for(int ibbcz=0;ibbcz<nbbcz;ibbcz++){
      for(int ihar=0;ihar<nhar;ihar++){
        for(int isub=0;isub<nsub;isub++){
            sprintf(name,"psi_%d_%d_%d_%d",icent,ibbcz,ihar,isub);
            psi[icent][ibbcz][ihar][isub]=new TH1F(name,name,100,-pi,pi);
            sprintf(name,"psiFla_%d_%d_%d_%d",icent,ibbcz,ihar,isub);
            psiFla[icent][ibbcz][ihar][isub]=new TH1F(name,name,100,-pi,pi);
            sprintf(name,"phiweight_%d_%d_%d_%d",icent,ibbcz,ihar,isub);
            phiweight[icent][ibbcz][ihar][isub]=new TH1F(name,name,50,-pi,pi);
            sprintf(name,"phiweightbbc_%d_%d_%d_%d",icent,ibbcz,ihar,isub);
            phiweightbbc[icent][ibbcz][ihar][isub]=new TProfile(name,name,200,0,200,0,100);
          for(int ixy=0;ixy<nxy;ixy++){
            sprintf(name,"q_%d_%d_%d_%d_%d",icent,ibbcz,ihar,isub,ixy);
            q[icent][ibbcz][ihar][isub][ixy]=new TH1F(name,name,820,-4.1,4.1);
            sprintf(name,"qRec_%d_%d_%d_%d_%d",icent,ibbcz,ihar,isub,ixy);
            qRec[icent][ibbcz][ihar][isub][ixy]=new TH1F(name,name,820,-4.1,4.1);
          }   
          for(int iord=0;iord<nord;iord++){
            sprintf(name,"fltsin_%d_%d_%d_%d_%d",icent,ibbcz,ihar,isub,iord);
            sinflt[icent][ibbcz][ihar][isub][iord]=new TH1F(name,name,220,-1.1,1.1);
            sprintf(name,"fltcos_%d_%d_%d_%d_%d",icent,ibbcz,ihar,isub,iord);
            cosflt[icent][ibbcz][ihar][isub][iord]=new TH1F(name,name,220,-1.1,1.1);
            }
      }
    }
  }
  }

  for(int icent=0;icent<ncent;icent++){
      for(int ihar=0;ihar<nhar;ihar++){
          /*
        sprintf(name,"EPRCNTBBCS_%d_%d",icent,ihar);
        EPRCNTBBCS[icent][ihar] = new TH1F(name,name,220,-1.1,1.1);
        sprintf(name,"EPRCNTFVTX1S_%d_%d",icent,ihar);
        EPRCNTFVTX1S[icent][ihar] = new TH1F(name,name,220,-1.1,1.1);
        sprintf(name,"EPRBBCSFVTX1S_%d_%d",icent,ihar);
        EPRBBCSFVTX1S[icent][ihar] = new TH1F(name,name,220,-1.1,1.1);
        sprintf(name,"EPRCNTFVTX1S_%d_%d",icent,ihar);
        EPRCNTFVTX2S[icent][ihar] = new TH1F(name,name,220,-1.1,1.1);
        sprintf(name,"EPRCNTFVTX1LS_%d_%d",icent,ihar);
        EPRCNTFVTX1S[icent][ihar] = new TH1F(name,name,220,-1.1,1.1);
        sprintf(name,"EPRCNTFVTX2S_%d_%d",icent,ihar);
        EPRCNTFVTX1LS[icent][ihar] = new TH1F(name,name,220,-1.1,1.1);
        sprintf(name,"EPRCNTFVTX2LS_%d_%d",icent,ihar);
        EPRCNTFVTX2LS[icent][ihar] = new TH1F(name,name,220,-1.1,1.1);
        sprintf(name,"EPRCNTFVTX3LS_%d_%d",icent,ihar);
        EPRCNTFVTX3LS[icent][ihar] = new TH1F(name,name,220,-1.1,1.1);
        sprintf(name,"EPRCNTFVTX4LS_%d_%d",icent,ihar);
        EPRCNTFVTX4LS[icent][ihar] = new TH1F(name,name,220,-1.1,1.1);
        sprintf(name,"EPRBBCSFVTX2S_%d_%d",icent,ihar);
        EPRBBCSFVTX2S[icent][ihar] = new TH1F(name,name,220,-1.1,1.1);
        */
        sprintf(name,"EPRBBCSFVTX1S_%d_%d",icent,ihar);
        EPRBBCSFVTX1S[icent][ihar] = new TH1F(name,name,220,-1.1,1.1);
        sprintf(name,"EPRFVTX1NFVTX1S_%d_%d",icent,ihar);
        EPRFVTX1NFVTX1S[icent][ihar] = new TH1F(name,name,220,-1.1,1.1);
        sprintf(name,"EPRBBCSFVTX1LS_%d_%d",icent,ihar);
        EPRBBCSFVTX1LS[icent][ihar] = new TH1F(name,name,220,-1.1,1.1);
        sprintf(name,"EPRBBCSFVTX2LS_%d_%d",icent,ihar);
        EPRBBCSFVTX2LS[icent][ihar] = new TH1F(name,name,220,-1.1,1.1);
        sprintf(name,"EPRBBCSFVTX3LS_%d_%d",icent,ihar);
        EPRBBCSFVTX3LS[icent][ihar] = new TH1F(name,name,220,-1.1,1.1);
        sprintf(name,"EPRBBCSFVTX4LS_%d_%d",icent,ihar);
        EPRBBCSFVTX4LS[icent][ihar] = new TH1F(name,name,220,-1.1,1.1);
        sprintf(name,"EPRBBCSFVTX1p2p3LS_%d_%d",icent,ihar);
        EPRBBCSFVTX1p2p3LS[icent][ihar] = new TH1F(name,name,220,-1.1,1.1);
        sprintf(name,"EPRBBCSFVTX1p2p4LS_%d_%d",icent,ihar);
        EPRBBCSFVTX1p2p4LS[icent][ihar] = new TH1F(name,name,220,-1.1,1.1);
        sprintf(name,"EPRBBCSFVTXtrkS_%d_%d",icent,ihar);
        EPRBBCSFVTXtrkS[icent][ihar] = new TH1F(name,name,220,-1.1,1.1);
        /*
        sprintf(name,"EPRCNTcBBCSc_%d_%d",icent,ihar);
        EPRCNTcBBCSc[icent][ihar] = new TH1F(name,name,220,-1.1,1.1);
        sprintf(name,"EPRBBCScFVTX1S_%d_%d",icent,ihar);
        EPRBBCScFVTX1S[icent][ihar] = new TH1F(name,name,220,-1.1,1.1);
        sprintf(name,"EPRCNTcFVTX1S_%d_%d",icent,ihar);
        EPRCNTcFVTX1S[icent][ihar] = new TH1F(name,name,220,-1.1,1.1);

        sprintf(name,"EPRCNTFVTX1trkS_%d_%d",icent,ihar);
        EPRCNTFVTX1trkS[icent][ihar] = new TH1F(name,name,220,-1.1,1.1);
        sprintf(name,"EPRCNTFVTX2trkS_%d_%d",icent,ihar);
        EPRCNTFVTX2trkS[icent][ihar] = new TH1F(name,name,220,-1.1,1.1);

        sprintf(name,"EPRBBCSFVTX1trkS_%d_%d",icent,ihar);
        EPRBBCSFVTX1trkS[icent][ihar] = new TH1F(name,name,220,-1.1,1.1);
        sprintf(name,"EPRBBCSFVTX2trkS_%d_%d",icent,ihar);
        EPRBBCSFVTX2trkS[icent][ihar] = new TH1F(name,name,220,-1.1,1.1);
        sprintf(name,"EPRCNTcBBCSc_%d_%d",icent,ihar);
        EPRCNTcBBCSc[icent][ihar] = new TH1F(name,name,220,-1.1,1.1);
        sprintf(name,"EPRBBCScFVTX1Sc_%d_%d",icent,ihar);
        EPRBBCScFVTX1Sc[icent][ihar] = new TH1F(name,name,220,-1.1,1.1);
        sprintf(name,"EPRCNTcFVTX1Sc_%d_%d",icent,ihar);
        EPRCNTcFVTX1Sc[icent][ihar] = new TH1F(name,name,220,-1.1,1.1);
        */
        for(int iphi=0;iphi<nphi;iphi++){
        /*
        sprintf(name,"vobsBBCS_%d_%d_%d",icent,ihar,iphi);
        vobsBBCS[icent][ihar][iphi] = new TProfile(name,name,60,0,6,-1.1,1.1);
        sprintf(name,"vobsFVTX1S_%d_%d_%d",icent,ihar,iphi);
        vobsFVTX1S[icent][ihar][iphi] = new TProfile(name,name,60,0,6,-1.1,1.1);
        sprintf(name,"vobsFVTX2S_%d_%d_%d",icent,ihar,iphi);
        vobsFVTX2S[icent][ihar][iphi] = new TProfile(name,name,60,0,6,-1.1,1.1);
        */

        /*
        sprintf(name,"vobsBBCSsq_%d_%d_%d",icent,ihar,iphi);
        vobsBBCSsq[icent][ihar][iphi] = new TProfile(name,name,60,0,6,-1.1,1.1);
        sprintf(name,"vobsFVTX1Ssq_%d_%d_%d",icent,ihar,iphi);
        vobsFVTX1Ssq[icent][ihar][iphi] = new TProfile(name,name,60,0,6,-1.1,1.1);
        sprintf(name,"vobsFVTX2Ssq_%d_%d_%d",icent,ihar,iphi);
        vobsFVTX2Ssq[icent][ihar][iphi] = new TProfile(name,name,60,0,6,-1.1,1.1);
        */
        sprintf(name,"vobsFVTX1LS_%d_%d_%d",icent,ihar,iphi);
        vobsFVTX1LS[icent][ihar][iphi] = new TH2F(name,name,60,0,6,220,-1.1,1.1);
        sprintf(name,"vobsFVTX2LS_%d_%d_%d",icent,ihar,iphi);
        vobsFVTX2LS[icent][ihar][iphi] = new TH2F(name,name,60,0,6,220,-1.1,1.1);
        sprintf(name,"vobsFVTX3LS_%d_%d_%d",icent,ihar,iphi);
        vobsFVTX3LS[icent][ihar][iphi] = new TH2F(name,name,60,0,6,220,-1.1,1.1);
        sprintf(name,"vobsFVTX4LS_%d_%d_%d",icent,ihar,iphi);
        vobsFVTX4LS[icent][ihar][iphi] = new TH2F(name,name,60,0,6,220,-1.1,1.1);
        sprintf(name,"vobsFVTX1p2p3LS_%d_%d_%d",icent,ihar,iphi);
        vobsFVTX1p2p3LS[icent][ihar][iphi] = new TH2F(name,name,60,0,6,220,-1.1,1.1);
        sprintf(name,"vobsFVTX1p2p4LS_%d_%d_%d",icent,ihar,iphi);
        vobsFVTX1p2p4LS[icent][ihar][iphi] = new TH2F(name,name,60,0,6,220,-1.1,1.1);
        sprintf(name,"vobsFVTX1S_%d_%d_%d",icent,ihar,iphi);
        vobsFVTX1S[icent][ihar][iphi] = new TH2F(name,name,60,0,6,220,-1.1,1.1);
        sprintf(name,"vobsFVTX1N_%d_%d_%d",icent,ihar,iphi);
        vobsFVTX1N[icent][ihar][iphi] = new TH2F(name,name,60,0,6,220,-1.1,1.1);
        sprintf(name,"vobsBBCS_%d_%d_%d",icent,ihar,iphi);
        vobsBBCS[icent][ihar][iphi] = new TH2F(name,name,60,0,6,220,-1.1,1.1);
        sprintf(name,"vobsFVTXtrkS_%d_%d_%d",icent,ihar,iphi);
        vobsFVTXtrkS[icent][ihar][iphi] = new TH2F(name,name,60,0,6,220,-1.1,1.1);

        sprintf(name,"vFVTX1LS_%d_%d_%d",icent,ihar,iphi);
        vFVTX1LS[icent][ihar][iphi] = new TH2F(name,name,60,0,6,200,-4,4);
        sprintf(name,"vFVTX2LS_%d_%d_%d",icent,ihar,iphi);
        vFVTX2LS[icent][ihar][iphi] = new TH2F(name,name,60,0,6,200,-4,4);
        sprintf(name,"vFVTX3LS_%d_%d_%d",icent,ihar,iphi);
        vFVTX3LS[icent][ihar][iphi] = new TH2F(name,name,60,0,6,200,-4,4);
        sprintf(name,"vFVTX4LS_%d_%d_%d",icent,ihar,iphi);
        vFVTX4LS[icent][ihar][iphi] = new TH2F(name,name,60,0,6,200,-4,4);
        sprintf(name,"vFVTX1p2p3LS_%d_%d_%d",icent,ihar,iphi);
        vFVTX1p2p3LS[icent][ihar][iphi] = new TH2F(name,name,60,0,6,200,-4,4);
        sprintf(name,"vFVTX1p2p4LS_%d_%d_%d",icent,ihar,iphi);
        vFVTX1p2p4LS[icent][ihar][iphi] = new TH2F(name,name,60,0,6,200,-4,4);
        sprintf(name,"vFVTX1N_%d_%d_%d",icent,ihar,iphi);
        vFVTX1N[icent][ihar][iphi] = new TH2F(name,name,60,0,6,200,-4,4);
        sprintf(name,"vFVTX1S_%d_%d_%d",icent,ihar,iphi);
        vFVTX1S[icent][ihar][iphi] = new TH2F(name,name,60,0,6,200,-4,4);
        sprintf(name,"vBBCS_%d_%d_%d",icent,ihar,iphi);
        vBBCS[icent][ihar][iphi] = new TH2F(name,name,60,0,6,200,-4,4);
        sprintf(name,"vFVTXtrkS_%d_%d_%d",icent,ihar,iphi);
        vFVTXtrkS[icent][ihar][iphi] = new TH2F(name,name,60,0,6,200,-4,4);

        sprintf(name,"vnFVTX1LS_%d_%d_%d",icent,ihar,iphi);
        vnFVTX1LS[icent][ihar][iphi] = new TH2F(name,name,60,0,6,100,-2,2);
        sprintf(name,"vnFVTX2LS_%d_%d_%d",icent,ihar,iphi);
        vnFVTX2LS[icent][ihar][iphi] = new TH2F(name,name,60,0,6,100,-2,2);
        sprintf(name,"vnFVTX3LS_%d_%d_%d",icent,ihar,iphi);
        vnFVTX3LS[icent][ihar][iphi] = new TH2F(name,name,60,0,6,100,-2,2);
        sprintf(name,"vnFVTX4LS_%d_%d_%d",icent,ihar,iphi);
        vnFVTX4LS[icent][ihar][iphi] = new TH2F(name,name,60,0,6,100,-2,2);
        sprintf(name,"vnFVTX1p2p3LS_%d_%d_%d",icent,ihar,iphi);
        vnFVTX1p2p3LS[icent][ihar][iphi] = new TH2F(name,name,60,0,6,100,-2,2);
        sprintf(name,"vnFVTX1p2p4LS_%d_%d_%d",icent,ihar,iphi);
        vnFVTX1p2p4LS[icent][ihar][iphi] = new TH2F(name,name,60,0,6,100,-2,2);
        sprintf(name,"vnFVTX1N_%d_%d_%d",icent,ihar,iphi);
        vnFVTX1N[icent][ihar][iphi] = new TH2F(name,name,60,0,6,100,-2,2);
        sprintf(name,"vnFVTX1S_%d_%d_%d",icent,ihar,iphi);
        vnFVTX1S[icent][ihar][iphi] = new TH2F(name,name,60,0,6,100,-2,2);
        sprintf(name,"vnBBCS_%d_%d_%d",icent,ihar,iphi);
        vnBBCS[icent][ihar][iphi] = new TH2F(name,name,60,0,6,100,-2,2);
        sprintf(name,"vnFVTXtrkS_%d_%d_%d",icent,ihar,iphi);
        vnFVTXtrkS[icent][ihar][iphi] = new TH2F(name,name,60,0,6,100,-2,2);
        /*
        sprintf(name,"vobsFVTX2S_%d_%d_%d",icent,ihar,iphi);
        vobsFVTX2S[icent][ihar][iphi] = new TProfile(name,name,60,0,6,-1.1,1.1);
        sprintf(name,"vobsFVTX2Ssq_%d_%d_%d",icent,ihar,iphi);
        vobsFVTX2Ssq[icent][ihar][iphi] = new TProfile(name,name,60,0,6,-1.1,1.1);
        sprintf(name,"vobscBBCSc_%d_%d_%d",icent,ihar,iphi);
        vobscBBCSc[icent][ihar][iphi] = new TProfile(name,name,60,0,6,-1.1,1.1);
        sprintf(name,"vobscBBCScsq_%d_%d_%d",icent,ihar,iphi);
        vobscBBCScsq[icent][ihar][iphi] = new TProfile(name,name,60,0,6,-1.1,1.1);
        sprintf(name,"vobscFVTX1Sc_%d_%d_%d",icent,ihar,iphi);
        vobscFVTX1Sc[icent][ihar][iphi] = new TProfile(name,name,60,0,6,-1.1,1.1);
        sprintf(name,"vobscFVTX1Scsq_%d_%d_%d",icent,ihar,iphi);
        vobscFVTX1Scsq[icent][ihar][iphi] = new TProfile(name,name,60,0,6,-1.1,1.1);
        sprintf(name,"vobsFVTX1trkS_%d_%d_%d",icent,ihar,iphi);
        vobsFVTX1trkS[icent][ihar][iphi] = new TProfile(name,name,60,0,6,-1.1,1.1);
        sprintf(name,"vobsFVTX1trkSsq_%d_%d_%d",icent,ihar,iphi);
        vobsFVTX1trkSsq[icent][ihar][iphi] = new TProfile(name,name,60,0,6,-1.1,1.1);
        sprintf(name,"vobsFVTX2trkS_%d_%d_%d",icent,ihar,iphi);
        vobsFVTX2trkS[icent][ihar][iphi] = new TProfile(name,name,60,0,6,-1.1,1.1);
        sprintf(name,"vobsFVTX2trkSsq_%d_%d_%d",icent,ihar,iphi);
        vobsFVTX2trkSsq[icent][ihar][iphi] = new TProfile(name,name,60,0,6,-1.1,1.1);
        */
        }
  }
  }

for(int istation=0;istation<4;istation++){
    sprintf(name,"hfvtx_%d",istation);
    hfvtx[istation] = new TH2F(name,name,200,-20,20,200,-20,20);
}
  sprintf(name,"hbbc");
  hbbc = new TH2F(name,name,200,-200,200,200,-200,200);
  sprintf(name,"hnfvtxtrk");
  hnfvtxtrk = new TH1F(name,name,100,0,100);
  sprintf(name,"hbbcsnfvtxtrk");
  hbbcsnfvtxtrk = new TH2F(name,name,200,0,200,100,0,100);
  sprintf(name,"hfvtxtrk");
  hfvtxtrk = new TH2F(name,name,100,-4,4,50,-pi,pi);
  return 0;
}

int EPAnaRun16alltree::Inittree(){
  //------------------------------------------------------------//
  //               Initializing Tree Variables                  //
  //------------------------------------------------------------//
  tree = NULL;
  tree1 = NULL;
  tree = new TChain("tree");
  tree1 = new TChain("tree");
  for(unsigned int itree=0;itree<InputFileName.size();itree++){
     cout<<InputFileName[itree]<<endl;
     try{
         tree->Add(InputFileName[itree]);}
     catch(int e){ cout << "An exception occurred. Exception Nr. " << e << '\n';}
  }
  for(unsigned int itree=0;itree<InputFileName1.size();itree++){
     cout<<InputFileName1[itree]<<endl;
     try{
         tree1->Add(InputFileName1[itree]);}
     catch(int e){ cout << "An exception occurred. Exception Nr. " << e << '\n';}
  }
  cout << "Now getting ready to read in the tree branch addresses and stuff...." << endl;

  // List of branches
  TBranch* b_run;
//  TBranch* b_event;   //!
  TBranch* b_bbc_z;   //!
  TBranch* b_centrality;   //!
  TBranch* b_bbc_qn;   //!
  TBranch* b_bbc_qs;   //!
  TBranch* b_npc1;   //!
  TBranch* b_trigger_scaled;   //!
//  TBranch* b_trigger_live;   //!
//  TBranch* b_d_Qx;   //!
//  TBranch* b_d_Qy;   //!
//  TBranch* b_d_Qw;   //!
//  TBranch* b_bc_x;   //!
//  TBranch* b_bc_y;   //!
//  TBranch* b_vtx_z;   //!
  TBranch* b_fvtx_x;   //!
  TBranch* b_fvtx_y;   //!
  TBranch* b_fvtx_z;   //!
  TBranch* b_nbbc;   //!
  TBranch* b_d_BBC_charge;   //!
  TBranch* b_d_BBC_time0;   //!
  TBranch* b_d_BBC_valid;   //!
  TBranch* b_d_nFVTX_clus;   //!
//  TBranch* b_d_nFVTXN_clus;   //!
//  TBranch* b_d_nFVTXS_clus;   //!
  TBranch* b_d_FVTX_x;   //!
  TBranch* b_d_FVTX_y;   //!
  TBranch* b_d_FVTX_z;   //!
  TBranch* b_ntrk;   //!
  TBranch* b_mom;   //!
  TBranch* b_phi0;   //!
  TBranch* b_the0;   //!
  TBranch* b_charge;   //!
  TBranch* b_pc3dphi;   //!
  TBranch* b_pc3dz;   //!
  // TBranch* b_d_ntrk;   //!
  // TBranch* b_d_cntpx;   //!
  // TBranch* b_d_cntpy;   //!
  // TBranch* b_d_cntpz;   //!
  TBranch* b_d_nfvtxtrk;   //!
  TBranch* b_d_nfvtxtrk1;  //!
  TBranch* b_fvtxchi;
  TBranch* b_farm;
  TBranch* b_fnhits;
  TBranch* b_feta;
  TBranch* b_fphi;
  TBranch* b_fvtxX;
  TBranch* b_fvtxY;
  TBranch* b_fvtxZ;


  tree->SetBranchAddress("run",&RunNumber,&b_run);
  tree->SetBranchAddress("bbcv",&d_bbcz,&b_bbc_z);
  tree->SetBranchAddress("cent",&centrality,&b_centrality);
  tree->SetBranchAddress("bbc_n",&bbc_qn,&b_bbc_qn);
  tree->SetBranchAddress("bbc_s",&bbc_qs,&b_bbc_qs);
  tree->SetBranchAddress("npc1hits",&npc1,&b_npc1);
  tree->SetBranchAddress("trig",&trigger_scaled,&b_trigger_scaled);

  tree->SetBranchAddress("fvtxx",&eventfvtx_x,&b_fvtx_x);
  tree->SetBranchAddress("fvtxy",&eventfvtx_y,&b_fvtx_y);
  tree->SetBranchAddress("fvtxz",&eventfvtx_z,&b_fvtx_z);

  tree->SetBranchAddress("nbbc", &d_nbbc,&b_nbbc);
  tree->SetBranchAddress("bbccharge",d_BBC_charge,&b_d_BBC_charge);
  tree->SetBranchAddress("bbct0",d_BBC_time0,&b_d_BBC_time0);
  tree1->SetBranchAddress("bbcvalid",d_BBC_valid,&b_d_BBC_valid);

  tree->SetBranchAddress("nclus",&d_nFVTX_clus,&b_d_nFVTX_clus);
  tree->SetBranchAddress("fclusX",d_FVTX_x,&b_d_FVTX_x);
  tree->SetBranchAddress("fclusY",d_FVTX_y,&b_d_FVTX_y);
  tree->SetBranchAddress("fclusZ",d_FVTX_z,&b_d_FVTX_z);

  tree->SetBranchAddress("ntrack",&d_ntrk,&b_ntrk);
  tree->SetBranchAddress("mom",d_mom,&b_mom);
  tree->SetBranchAddress("phi0",d_phi0,&b_phi0);
  tree->SetBranchAddress("the0",d_the0,&b_the0);
  tree->SetBranchAddress("charge",d_charge,&b_charge);
  tree->SetBranchAddress("pc3dphi",d_pc3dphi,&b_pc3dphi);
  tree->SetBranchAddress("pc3dz",d_pc3dz,&b_pc3dz);
  
  tree->SetBranchAddress("nfvtxtrack",&d_nfvtxtrk,&b_d_nfvtxtrk);
  tree1->SetBranchAddress("nfvtxtrack",&d_nfvtxtrk1,&b_d_nfvtxtrk1);
  tree1->SetBranchAddress("fvtxchi2",d_fvtxchi,&b_fvtxchi);
  tree->SetBranchAddress("farm",d_farm,&b_farm);
  tree->SetBranchAddress("fnhits",d_fnhits,&b_fnhits);
  tree->SetBranchAddress("feta",d_feta,&b_feta);
  tree->SetBranchAddress("fphi",d_fphi,&b_fphi);
  tree->SetBranchAddress("fvtxX",d_fvtxX,&b_fvtxX);
  tree->SetBranchAddress("fvtxY",d_fvtxY,&b_fvtxY);
  tree->SetBranchAddress("fvtxZ",d_fvtxZ,&b_fvtxZ);


  return 0;
}

//_____________________________________________________________________________________________________________________________
int EPAnaRun16alltree::process_event()
{
  int nEvent = tree->GetEntries();
  cout<<nEvent<<endl;
  for(ievent=0;ievent < nEvent; ievent++){
      tree->GetEntry(ievent);
      tree1->GetEntry(ievent);
  if(ievent%10000==0) {
    cout<<"EPall calFlag = "<< calFlag << "************* ievent= "<<ievent<<"    *************"<<endl;
  }
//  float pi = acos(-1.0);
  
//event check
   if(d_nfvtxtrk != d_nfvtxtrk1) {cout << d_nfvtxtrk << " EVENT NOT SYNC ! " << d_nfvtxtrk1 << endl; continue;}

//global
    unsigned int trigger_FVTXNSBBCScentral = 0x00100000;
    unsigned int trigger_FVTXNSBBCS        = 0x00400000;
    unsigned int trigger_BBCLL1narrowcent  = 0x00000008;
    unsigned int trigger_BBCLL1narrow      = 0x00000010;
    unsigned int accepted_triggers = 0;
      // --- Run16dAu200                                                                                                    
   if ( RunNumber >= 454774 && RunNumber <= 455639 ) accepted_triggers = trigger_BBCLL1narrowcent | trigger_BBCLL1narrow;
      // --- Run16dAu62                                                                                                                     
   if ( RunNumber >= 455792 && RunNumber <= 456283 ) accepted_triggers = trigger_BBCLL1narrowcent | trigger_BBCLL1narrow;
      // --- Run16dAu20                                                                   
   if ( RunNumber >= 456652 && RunNumber <= 457298 ) accepted_triggers = trigger_FVTXNSBBCScentral | trigger_FVTXNSBBCS;
      // --- Run16dAu39                                                                        
   if ( RunNumber >= 457634 && RunNumber <= 458167 ) accepted_triggers = trigger_FVTXNSBBCScentral | trigger_FVTXNSBBCS;

  unsigned int passes_trigger = trigger_scaled & accepted_triggers;
  if ( passes_trigger == 0 )      continue;
  float bbcv = d_bbcz;
  int cent = centrality;
//  float bbc_s = bbc_qs;
  float vvertex = bbcv;
  
//  float fvtxx = eventfvtx_x;
//  float fvtxy = eventfvtx_y;
  float fvtxz = eventfvtx_z;
  if(fvtxz != fvtxz){
      fvtxz = -9999;
  }
  
  if(RunNumber>=455792 && RunNumber<=458167){ // low energies use fvtx as vertex if there is one
  if(fvtxz != -9999)
    vvertex = fvtxz;
  }

 int ibbcz = -1;
 ibbcz = nbbcz*(vvertex+10)/20;
 if(ibbcz<0||ibbcz>=nbbcz) continue;

 int icent = -9999;
 if(cent<0) continue;
 if ( RunNumber >= 454774 && RunNumber <= 455639 && cent <= 5) icent = 0;
 if ( RunNumber >= 455792 && RunNumber <= 456283 && cent <= 10) icent = 0;
 if ( RunNumber >= 456652 && RunNumber <= 457298 && cent <= 20) icent = 0;
 if ( RunNumber >= 457634 && RunNumber <= 458167 && cent <= 20) icent = 0;
 if(icent<0) continue;

 // --- all numbers from Darren 2016-06-23
      const float x_off = 0.3;
      const float beam_angle = 0.001;
      float vtx_z = vvertex;
      float vtx_x = x_off + atan(beam_angle)*vtx_z;
      float vtx_y = 0.02;
    
 if ( RunNumber >= 456652 && RunNumber <= 457298 && d_nFVTX_clus > 300) continue;
 if ( RunNumber >= 457634 && RunNumber <= 458167 && d_nFVTX_clus > 500) continue;

  jevent++;

 float Qx[nsub][nhar];
 float Qy[nsub][nhar];
 float Qw[nsub][nhar];
 float RPPhi[nsub][nhar];
 bool GoodPsi[nsub][nhar];

 int iharE=0;
 if(nhar == 1 || nhar == 2) iharE=1;

 for(int ihar=0;ihar<nhar;ihar++){

    for(int isub=0;isub<nsub;isub++){
      Qx[isub][ihar] = 0;
      Qy[isub][ihar] = 0;
      Qw[isub][ihar] = 0;
      RPPhi[isub][ihar] = -9999;
      GoodPsi[isub][ihar] = false;
    }

  int n = ihar+1.0+iharE;

//FVTX tracks
    for(int iftrk = 0; iftrk < d_nfvtxtrk; iftrk++){
      int fvtx_arm = d_farm[iftrk];
      int nhits = d_fnhits[iftrk];
      float fvtx_eta = d_feta[iftrk];
      float fvtx_the = 2.*atan(exp(-fvtx_eta));
      float fvtx_phi = d_fphi[iftrk];

      float     fvtx_x = d_fvtxX[iftrk];
      float     fvtx_y = d_fvtxY[iftrk];
      float     fvtx_z = d_fvtxZ[iftrk];
      float     chi2   = d_fvtxchi[iftrk];
      if(fvtx_arm) continue;
      if(fvtx_the==0) continue;
      if(chi2>4) continue;
      float weight = 1.;
 
      float DCA_x = fvtx_x + tan(fvtx_the)*cos(fvtx_phi)*(fvtxz - fvtx_z);
      float DCA_y = fvtx_y + tan(fvtx_the)*sin(fvtx_phi)*(fvtxz - fvtx_z);

     bool dcacut = fabs(DCA_x)<2.0 && fabs(DCA_y)<2.0; //fabs(sigma_dcax)<2.0 && fabs(sigma_dcay)<2.0;
      if(fvtx_phi<10 && fvtx_phi > -10 && fabs(fvtx_eta)<3.5 && dcacut && nhits>=3){
            Qx[9][ihar] += weight * cos(n*fvtx_phi);
            Qy[9][ihar] += weight * sin(n*fvtx_phi);
            Qw[9][ihar] += weight;
        }
}

    for(int iclus = 0; iclus < d_nFVTX_clus; iclus++){
        float fvtx_x      = d_FVTX_x[iclus];
        float fvtx_y      = d_FVTX_y[iclus];
        float fvtx_z      = d_FVTX_z[iclus];
    int istation = get_fvtx_layer(fvtx_z);
    if(istation<0) continue;
    fvtx_x = fvtx_x-vtx_x;
    fvtx_y = fvtx_y-vtx_y;
    fvtx_z = fvtx_z-vtx_z;
    fvtx_x = fvtx_z*sin(-beam_angle) + fvtx_x*cos(-beam_angle);
    double weight = 1.0;
//    int iarm = 0;
//    if(fvtx_z>0) iarm = 1; 
    float fvtx_r = sqrt(pow(fvtx_x,2.0)+pow(fvtx_y,2.0));
    if ( RunNumber >= 456652 && RunNumber <= 458167 && fvtx_r < 5.2 ) continue; //radius cut
    if( (fabs(fvtx_x)>999) ||(fabs(fvtx_y)>999) || (fabs(fvtx_z)>999)) continue;
    float fvtx_the = atan2(fvtx_r,fvtx_z);
    float fvtx_phi = atan2(fvtx_y,fvtx_x);
    float fvtx_eta = -log(tan(0.5*fvtx_the));
    if(!(fabs(fvtx_eta)<3.5)) continue;
    if(calFlag == 0){
        phiweight[icent][ibbcz][ihar][istation]->Fill(fvtx_phi);
      //  phiweight[icent][ibbcz][ihar][5]->Fill(fvtx_phi);
    }
    else{
        int ibin = phiweight[icent][ibbcz][ihar][istation]->FindBin(fvtx_phi);
        float binc = phiweight[icent][ibbcz][ihar][istation]->GetBinContent(ibin);
        if(binc!=0){
        if(ibin > 0 && ibin<=50)
         weight = phiweight[icent][ibbcz][ihar][istation]->Integral()/phiweight[icent][ibbcz][ihar][istation]->GetNbinsX()/binc;
        //if(fabs(weight-1.)>0.2) weight = 0.;
        else weight = 0.;
        }
        //else weight = 0.;
    if(istation < 4){
      Qx[istation][ihar] += weight * cos(n*fvtx_phi);
      Qy[istation][ihar] += weight * sin(n*fvtx_phi);
      Qw[istation][ihar] += weight;
    }
      /*
        ibin = phiweight[icent][ibbcz][ihar][5]->FindBin(fvtx_phi);
        binc = phiweight[icent][ibbcz][ihar][5]->GetBinContent(ibin);
        if(binc!=0){
        if(ibin > 0 && ibin<=50)
        weight = phiweight[icent][ibbcz][ihar][5]->Integral()/phiweight[icent][ibbcz][ihar][5]->GetNbinsX()/binc;
      //if(fabs(weight-1.)>0.2) weight = 0.;
        else weight = 0.;
        }
    else weight = 0.;
    */
    if(istation == 0 || istation == 1 || istation == 2){
     Qx[6][ihar] += weight * cos(n*fvtx_phi);
     Qy[6][ihar] += weight * sin(n*fvtx_phi);
     Qw[6][ihar] += weight;
    }
    if(istation == 0 || istation == 1 || istation == 3){
     Qx[7][ihar] += weight * cos(n*fvtx_phi);
     Qy[7][ihar] += weight * sin(n*fvtx_phi);
     Qw[7][ihar] += weight;
    }
    if(istation < 4){
      Qx[5][ihar] += weight * cos(n*fvtx_phi);
      Qy[5][ihar] += weight * sin(n*fvtx_phi);
      Qw[5][ihar] += weight;
    }
    if(istation >= 4){
      Qx[8][ihar] += weight * cos(n*fvtx_phi);
      Qy[8][ihar] += weight * sin(n*fvtx_phi);
      Qw[8][ihar] += weight;
    }

    }//calflag > 0 
    }//fvtx cluster loop

/*
// Tracks 
    for(int itrk=0; itrk< d_ntrk; itrk++){
      float mom    = d_mom[itrk];
      float phi0    = d_phi0[itrk];
      float the0    = d_the0[itrk];
      float pt        = mom*sin(the0);
      float pz        = mom*cos(the0);
      float px = pt * cos(phi0);
      float py = pt * sin(phi0);
      px = pz*sin(-beam_angle) + px*cos(-beam_angle);
      float phi = atan2(py,px);
      float weight = 1.0;
      phi = atan2(sin(phi),cos(phi));
     // float mass = 0.1396;//assume charged pion mass
      //float eta       = -log(tan(0.5*d_scnt->get_the0()));
      weight = pt;
      if(pt<0.0 || pt>2.0) weight = 0.0;

        Qx[6][ihar] += weight * cos(n*phi);
        Qy[6][ihar] += weight * sin(n*phi);
        Qw[6][ihar] += weight;
    }//track loop
*/

  // beam beam counter (bbc) r.p. // -----------------------------------
    int ibbc = 0;
    for(int ipmt = 0; ipmt < 128; ipmt++){
        float bbcx, bbcy, bbcz;
        if(!d_BBC_valid[ipmt]) continue;
    if(ipmt < 64){
    bbcx      = d_pmt_x[ipmt];
    bbcy      = d_pmt_y[ipmt];
    bbcz      = d_pmt_z;
    }
    else{
    bbcx      = d_pmt_x[ipmt-64];
    bbcy      = d_pmt_y[ipmt-64];
    bbcz      = -d_pmt_z;
    }
    float charge = d_BBC_charge[ibbc];
    float time0 = d_BBC_time0[ibbc];
    ibbc++;
    bbcx = bbcx - vtx_x*10.0;
    bbcy = bbcy - vtx_y*10.0;
    bbcz = bbcz - vtx_z*10.0;
    // --- rotation
    bbcx = bbcz*sin(-beam_angle) + bbcx*cos(-beam_angle);

    if (time0>0 && charge>0 && bbcz<0) {
      float phi=atan2(bbcy,bbcx);
      float weight = charge;
      if(calFlag == 0){
          phiweightbbc[icent][ibbcz][ihar][4]->Fill(ipmt,charge);
      }
      else{
           //        phiweight[icent][ibbcz][ihar][4] = GetPPbbcpmt(ibbcz);
           //        if(phiweight[icent][ibbcz][ihar][4]->GetBinContent(phiweight[icent][ibbcz][ihar][4]->FindBin(phi))!=0)
           //        weight = charge * phiweight[icent][ibbcz][ihar][4]->Integral()/phiweight[icent][ibbcz][ihar][4]->GetNbinsX()/phiweight[icent][ibbcz][ihar][4]->GetBinContent(phiweight[icent][ibbcz][ihar][4]->FindBin(phi));
          //        else weight = charge * 1.;
        int ibin = phiweightbbc[icent][ibbcz][ihar][4]->FindBin(ipmt);
        float binc = phiweightbbc[icent][ibbcz][ihar][4]->GetBinContent(ibin);
        if(binc!=0)
        weight = charge * PPbbcpmtweight[ibbcz][ipmt]/binc;
      //  else weight = 0;
      Qx[4][ihar] += weight * cos(n*phi);
      Qy[4][ihar] += weight * sin(n*phi);
      Qw[4][ihar] += weight;
    }
    }
    }//bbc loop

//fvtx tracks


if(calFlag>0){
    for(int isub=0;isub<nsub;isub++){

        if(Qw[isub][ihar]>0){
          float qx = Qx[isub][ihar]/Qw[isub][ihar];
          float qy = Qy[isub][ihar]/Qw[isub][ihar];
          if(calFlag==1){
            q[icent][ibbcz][ihar][isub][0]->Fill(qx);
            q[icent][ibbcz][ihar][isub][1]->Fill(qy);
          }
          else{
            qx=(qx-meanx[icent][ibbcz][ihar][isub])/rmsx[icent][ibbcz][ihar][isub];
            qy=(qy-meany[icent][ibbcz][ihar][isub])/rmsy[icent][ibbcz][ihar][isub];
            float Psi = atan2(qy,qx)/n;
            if(calFlag==2){
              qRec[icent][ibbcz][ihar][isub][0]->Fill(qx);
              qRec[icent][ibbcz][ihar][isub][1]->Fill(qy);
              psi[icent][ibbcz][ihar][isub]->Fill(Psi*n);  //make sure the Psi angle range (0,2*pi)
              for (int iord=0; iord<nord; iord++) {
                float cosPsi=cos((iord+1.0)*n*Psi);
                float sinPsi=sin((iord+1.0)*n*Psi);
                sinflt[icent][ibbcz][ihar][isub][iord]->Fill(sinPsi);
                cosflt[icent][ibbcz][ihar][isub][iord]->Fill(cosPsi);
              }
            }
            else{
              float dPsi = 0;
              for (int iord=0; iord<nord; iord++) {
                float cosPsi=cos((iord+1.0)*n*Psi);
                float sinPsi=sin((iord+1.0)*n*Psi);
                //dPsi+=(fltcos[icent][ibbcz][ihar][isub][iord]*sinPsi-fltsin[icent][ibbcz][ihar][isub][iord]*cosPsi)*2.0/(iord+1)/n;
                dPsi+=(cosfltarr[icent][ibbcz][ihar][isub][iord]*sinPsi-sinfltarr[icent][ibbcz][ihar][isub][iord]*cosPsi)*2.0/(iord+1)/n;
              }
            Psi+=dPsi;
            Psi=atan2(sin(n*Psi),cos(n*Psi))/n;
            RPPhi[isub][ihar] = Psi;
            GoodPsi[isub][ihar] = (fabs(Psi)<4);
            psiFla[icent][ibbcz][ihar][isub]->Fill(Psi*n);  //make sure the Psi angle range (0,2*pi)
            }
          }
        }
    }
    
    if(calFlag==3){
        /*
        if(GoodPsi[1][ihar] && GoodPsi[5][ihar])
          EPRCNTBBCS[icent][ihar]->Fill(cos(n*(RPPhi[1][ihar]-RPPhi[5][ihar])));
        if(GoodPsi[5][ihar] && GoodPsi[2][ihar])
          EPRCNTFVTX1S[icent][ihar]->Fill(cos(n*(RPPhi[5][ihar]-RPPhi[2][ihar])));
        if(GoodPsi[1][ihar] && GoodPsi[2][ihar])
          EPRBBCSFVTX1S[icent][ihar]->Fill(cos(n*(RPPhi[1][ihar]-RPPhi[2][ihar])));
        if(GoodPsi[2][ihar] && GoodPsi[5][ihar])
          EPRCNTFVTX1trkS[icent][ihar]->Fill(cos(n*(RPPhi[2][ihar]-RPPhi[5][ihar])));
        if(GoodPsi[3][ihar] && GoodPsi[5][ihar])
          EPRCNTFVTX2trkS[icent][ihar]->Fill(cos(n*(RPPhi[3][ihar]-RPPhi[5][ihar])));
        if(GoodPsi[1][ihar] && GoodPsi[2][ihar])
          EPRBBCSFVTX1trkS[icent][ihar]->Fill(cos(n*(RPPhi[1][ihar]-RPPhi[2][ihar])));
        if(GoodPsi[1][ihar] && GoodPsi[3][ihar])
          EPRBBCSFVTX2trkS[icent][ihar]->Fill(cos(n*(RPPhi[1][ihar]-RPPhi[3][ihar])));

        if(GoodPsi[0][ihar] && GoodPsi[4][ihar])
          EPRCNTcBBCSc[icent][ihar]->Fill(cos(n*(RPPhi[0][ihar]-RPPhi[4][ihar])));
        if(GoodPsi[4][ihar] && GoodPsi[3][ihar])
          EPRCNTcFVTX1Sc[icent][ihar]->Fill(cos(n*(RPPhi[4][ihar]-RPPhi[3][ihar])));
        if(GoodPsi[0][ihar] && GoodPsi[3][ihar])
          EPRBBCScFVTX1Sc[icent][ihar]->Fill(cos(n*(RPPhi[0][ihar]-RPPhi[3][ihar])));

        if(GoodPsi[2][ihar] && GoodPsi[5][ihar])
            EPRCNTBBCS[icent][ihar]->Fill(cos(n*(RPPhi[5][ihar]-RPPhi[2][ihar])));//cnt-bbcs
        if(GoodPsi[6][ihar] && GoodPsi[5][ihar])
            EPRCNTFVTX1S[icent][ihar]->Fill(cos(n*(RPPhi[5][ihar]-RPPhi[6][ihar])));//cnt-fvtx1s
        if(GoodPsi[7][ihar] && GoodPsi[5][ihar])
            EPRCNTFVTX2S[icent][ihar]->Fill(cos(n*(RPPhi[5][ihar]-RPPhi[7][ihar])));//cnt-fvtx2s
        if(GoodPsi[2][ihar] && GoodPsi[6][ihar])
            EPRBBCSFVTX1S[icent][ihar]->Fill(cos(n*(RPPhi[2][ihar]-RPPhi[6][ihar])));//bbcs-fvtx1s
        if(GoodPsi[2][ihar] && GoodPsi[7][ihar])
            EPRBBCSFVTX2S[icent][ihar]->Fill(cos(n*(RPPhi[2][ihar]-RPPhi[7][ihar])));//bbcs-fvtx2s
        */
        if(GoodPsi[4][ihar] && GoodPsi[0][ihar])
            EPRBBCSFVTX1LS[icent][ihar]->Fill(cos(n*(RPPhi[4][ihar]-RPPhi[0][ihar])));//bbcs-fvtx1ls
        if(GoodPsi[4][ihar] && GoodPsi[1][ihar])
            EPRBBCSFVTX2LS[icent][ihar]->Fill(cos(n*(RPPhi[4][ihar]-RPPhi[1][ihar])));//bbcs-fvtx2ls
        if(GoodPsi[4][ihar] && GoodPsi[2][ihar])
            EPRBBCSFVTX3LS[icent][ihar]->Fill(cos(n*(RPPhi[4][ihar]-RPPhi[2][ihar])));//bbcs-fvtx1ls
        if(GoodPsi[4][ihar] && GoodPsi[3][ihar])
            EPRBBCSFVTX4LS[icent][ihar]->Fill(cos(n*(RPPhi[4][ihar]-RPPhi[3][ihar])));//bbcs-fvtx1ls
        if(GoodPsi[4][ihar] && GoodPsi[5][ihar])
            EPRBBCSFVTX1S[icent][ihar]->Fill(cos(n*(RPPhi[4][ihar]-RPPhi[5][ihar])));//bbcs-fvtx1ls
        if(GoodPsi[4][ihar] && GoodPsi[6][ihar])
            EPRBBCSFVTX1p2p3LS[icent][ihar]->Fill(cos(n*(RPPhi[4][ihar]-RPPhi[6][ihar])));//bbcs-fvtx1ls
        if(GoodPsi[4][ihar] && GoodPsi[7][ihar])
            EPRBBCSFVTX1p2p4LS[icent][ihar]->Fill(cos(n*(RPPhi[4][ihar]-RPPhi[7][ihar])));//bbcs-fvtx1ls
        if(GoodPsi[5][ihar] && GoodPsi[8][ihar])
            EPRFVTX1NFVTX1S[icent][ihar]->Fill(cos(n*(RPPhi[5][ihar]-RPPhi[8][ihar])));//bbcs-fvtx1ls
        if(GoodPsi[4][ihar] && GoodPsi[9][ihar])
            EPRBBCSFVTXtrkS[icent][ihar]->Fill(cos(n*(RPPhi[4][ihar]-RPPhi[9][ihar])));//bbcs-fvtx1ls
        /*
        if(GoodPsi[5][ihar] && GoodPsi[8][ihar])
            EPRCNTFVTX1LS[icent][ihar]->Fill(cos(n*(RPPhi[5][ihar]-RPPhi[8][ihar])));//cnt-fvtx1ls
        if(GoodPsi[5][ihar] && GoodPsi[9][ihar])
            EPRCNTFVTX2LS[icent][ihar]->Fill(cos(n*(RPPhi[5][ihar]-RPPhi[9][ihar])));//cnt-fvtx2ls
        if(GoodPsi[5][ihar] && GoodPsi[10][ihar])
            EPRCNTFVTX3LS[icent][ihar]->Fill(cos(n*(RPPhi[5][ihar]-RPPhi[10][ihar])));//cnt-fvtx3ls
        if(GoodPsi[5][ihar] && GoodPsi[11][ihar])
            EPRCNTFVTX4LS[icent][ihar]->Fill(cos(n*(RPPhi[5][ihar]-RPPhi[11][ihar])));//cnt-fvtx4ls
*/
// Tracks 
    for(int itrk=0; itrk< d_ntrk; itrk++){
      float mom    = d_mom[itrk];
      float phi0    = d_phi0[itrk];
      float the0    = d_the0[itrk];
      int charge = d_charge[itrk];
      float pc3dphi   = d_pc3dphi[itrk];
      float pc3dz   = d_pc3dz[itrk];
      float pt        = mom*sin(the0);
      float pz        = mom*cos(the0);
      float px = pt * cos(phi0);
      float py = pt * sin(phi0);
      px = pz*sin(-beam_angle) + px*cos(-beam_angle);
      pt = sqrt(px*px+py*py);
      float phi = atan2(py,px);
      int dcarm=0;
      if(px>0) dcarm=1;

      double sdphi = calcsdphi(pc3dphi,dcarm,charge,mom,RunNumber);
      double sdz =  calcsdz(pc3dz,dcarm,charge,mom,RunNumber);
      int iphi = 0;
      if(fabs(sdphi)<2.0 && fabs(sdz)<2.0){
          if(dcarm==0) iphi = 0;
          else iphi = 1;
        if(GoodPsi[4][ihar]){
          vobsBBCS[icent][ihar][iphi]->Fill(pt, cos(n*(phi-RPPhi[4][ihar])));
          float dphi = phi-RPPhi[4][ihar];
          vBBCS[icent][ihar][iphi]->Fill(pt, atan2(sin(dphi),cos(dphi)));
          vnBBCS[icent][ihar][iphi]->Fill(pt, atan2(sin(n*dphi),cos(n*dphi))/n);
        }
        if(GoodPsi[0][ihar]){
          vobsFVTX1LS[icent][ihar][iphi]->Fill(pt, cos(n*(phi-RPPhi[0][ihar])));
          float dphi = phi-RPPhi[0][ihar];
          vFVTX1LS[icent][ihar][iphi]->Fill(pt, atan2(sin(dphi),cos(dphi)));
          vnFVTX1LS[icent][ihar][iphi]->Fill(pt, atan2(sin(n*dphi),cos(n*dphi))/n);
        }
        if(GoodPsi[1][ihar]){
          vobsFVTX2LS[icent][ihar][iphi]->Fill(pt, cos(n*(phi-RPPhi[1][ihar])));
          float dphi = phi-RPPhi[1][ihar];
          vFVTX2LS[icent][ihar][iphi]->Fill(pt, atan2(sin(dphi),cos(dphi)));
          vnFVTX2LS[icent][ihar][iphi]->Fill(pt, atan2(sin(n*dphi),cos(n*dphi))/n);
        }
        if(GoodPsi[2][ihar]){
          vobsFVTX3LS[icent][ihar][iphi]->Fill(pt, cos(n*(phi-RPPhi[2][ihar])));
          float dphi = phi-RPPhi[2][ihar];
          vFVTX3LS[icent][ihar][iphi]->Fill(pt, atan2(sin(dphi),cos(dphi)));
          vnFVTX3LS[icent][ihar][iphi]->Fill(pt, atan2(sin(n*dphi),cos(n*dphi))/n);
        }
        if(GoodPsi[3][ihar]){
          vobsFVTX4LS[icent][ihar][iphi]->Fill(pt, cos(n*(phi-RPPhi[3][ihar])));
          float dphi = phi-RPPhi[3][ihar];
          vFVTX4LS[icent][ihar][iphi]->Fill(pt, atan2(sin(dphi),cos(dphi)));
          vnFVTX4LS[icent][ihar][iphi]->Fill(pt, atan2(sin(n*dphi),cos(n*dphi))/n);
        }
        if(GoodPsi[5][ihar]){
          vobsFVTX1S[icent][ihar][iphi]->Fill(pt, cos(n*(phi-RPPhi[5][ihar])));
          float dphi = phi-RPPhi[5][ihar];
          vFVTX1S[icent][ihar][iphi]->Fill(pt, atan2(sin(dphi),cos(dphi)));
          vnFVTX1S[icent][ihar][iphi]->Fill(pt, atan2(sin(n*dphi),cos(n*dphi))/n);
        }
        if(GoodPsi[6][ihar]){
          vobsFVTX1p2p3LS[icent][ihar][iphi]->Fill(pt, cos(n*(phi-RPPhi[6][ihar])));
          float dphi = phi-RPPhi[6][ihar];
          vFVTX1p2p3LS[icent][ihar][iphi]->Fill(pt, atan2(sin(dphi),cos(dphi)));
          vnFVTX1p2p3LS[icent][ihar][iphi]->Fill(pt, atan2(sin(n*dphi),cos(n*dphi))/n);
        }
        if(GoodPsi[7][ihar]){
          vobsFVTX1p2p4LS[icent][ihar][iphi]->Fill(pt, cos(n*(phi-RPPhi[7][ihar])));
          float dphi = phi-RPPhi[7][ihar];
          vFVTX1p2p4LS[icent][ihar][iphi]->Fill(pt, atan2(sin(dphi),cos(dphi)));
          vnFVTX1p2p4LS[icent][ihar][iphi]->Fill(pt, atan2(sin(n*dphi),cos(n*dphi))/n);
        }
        if(GoodPsi[8][ihar]){
          vobsFVTX1N[icent][ihar][iphi]->Fill(pt, cos(n*(phi-RPPhi[8][ihar])));
          float dphi = phi-RPPhi[8][ihar];
          vFVTX1N[icent][ihar][iphi]->Fill(pt, atan2(sin(dphi),cos(dphi)));
          vnFVTX1N[icent][ihar][iphi]->Fill(pt, atan2(sin(n*dphi),cos(n*dphi))/n);
        }
        if(GoodPsi[9][ihar]){
          vobsFVTXtrkS[icent][ihar][iphi]->Fill(pt, cos(n*(phi-RPPhi[9][ihar])));
          float dphi = phi-RPPhi[9][ihar];
          vFVTXtrkS[icent][ihar][iphi]->Fill(pt, atan2(sin(dphi),cos(dphi)));
          vnFVTXtrkS[icent][ihar][iphi]->Fill(pt, atan2(sin(n*dphi),cos(n*dphi))/n);
        }
        /*
        if(GoodPsi[2][ihar]){
          vobsFVTX1trkS[icent][ihar][iphi]->Fill(pt, cos(n*(phi-RPPhi[2][ihar])));
          vobsFVTX1trkSsq[icent][ihar][iphi]->Fill(pt, cos(n*(phi-RPPhi[2][ihar]))*cos(n*(phi-RPPhi[2][ihar])));
        }
        if(GoodPsi[3][ihar]){
          vobsFVTX2trkS[icent][ihar][iphi]->Fill(pt, cos(n*(phi-RPPhi[3][ihar])));
          vobsFVTX2trkSsq[icent][ihar][iphi]->Fill(pt, cos(n*(phi-RPPhi[3][ihar]))*cos(n*(phi-RPPhi[3][ihar])));
        }
        */
/*
          float angleblue = 0.0016+0.0008-iangle1*0.0002;
          float angleyellow = pi+0.0036+0.0016-iangle2*0.0004;
          TLorentzVector particle_vec = RotateAndBoost(angleblue,angleyellow,px,py,pz,mass);
          float pt = TMath::Sqrt(particle_vec.Py()*particle_vec.Py()+particle_vec.Px()*particle_vec.Px());                   
          float phi = TMath::ATan2(particle_vec.Py(),particle_vec.Px());
        if(GoodPsi[0][ihar]){
          vobscBBCSc[icent][ihar][iphi]->Fill(pt, cos(n*(phi-RPPhi[0][ihar])));
          vobscBBCScsq[icent][ihar][iphi]->Fill(pt, cos(n*(phi-RPPhi[0][ihar]))*cos(n*(phi-RPPhi[0][ihar])));
        }
        if(GoodPsi[3][ihar]){
          vobscFVTX1Sc[icent][ihar][iphi]->Fill(pt, cos(n*(phi-RPPhi[3][ihar])));
          vobscFVTX1Scsq[icent][ihar][iphi]->Fill(pt, cos(n*(phi-RPPhi[3][ihar]))*cos(n*(phi-RPPhi[3][ihar])));
        }
        */
        /*
        if(GoodPsi[6][ihar]){
          vobsFVTX1S[icent][ihar][iphi]->Fill(pt, cos(n*(phi-RPPhi[6][ihar])));
          vobsFVTX1Ssq[icent][ihar][iphi]->Fill(pt, cos(n*(phi-RPPhi[6][ihar]))*cos(n*(phi-RPPhi[6][ihar])));
        }
        if(GoodPsi[7][ihar]){
          vobsFVTX2S[icent][ihar][iphi]->Fill(pt, cos(n*(phi-RPPhi[7][ihar])));
          vobsFVTX2Ssq[icent][ihar][iphi]->Fill(pt, cos(n*(phi-RPPhi[7][ihar]))*cos(n*(phi-RPPhi[7][ihar])));
        }
        
        if(GoodPsi[8][ihar]){
          vobsFVTX1LS[icent][ihar][iphi]->Fill(pt, cos(n*(phi-RPPhi[8][ihar])));
          vobsFVTX1LSsq[icent][ihar][iphi]->Fill(pt, cos(n*(phi-RPPhi[8][ihar]))*cos(n*(phi-RPPhi[8][ihar])));
        }
        if(GoodPsi[9][ihar]){
          vobsFVTX2LS[icent][ihar][iphi]->Fill(pt, cos(n*(phi-RPPhi[9][ihar])));
          vobsFVTX2LSsq[icent][ihar][iphi]->Fill(pt, cos(n*(phi-RPPhi[9][ihar]))*cos(n*(phi-RPPhi[9][ihar])));
        }
        if(GoodPsi[10][ihar]){
          vobsFVTX3LS[icent][ihar][iphi]->Fill(pt, cos(n*(phi-RPPhi[10][ihar])));
          vobsFVTX3LSsq[icent][ihar][iphi]->Fill(pt, cos(n*(phi-RPPhi[10][ihar]))*cos(n*(phi-RPPhi[10][ihar])));
        }
        if(GoodPsi[11][ihar]){
          vobsFVTX4LS[icent][ihar][iphi]->Fill(pt, cos(n*(phi-RPPhi[11][ihar])));
          vobsFVTX4LSsq[icent][ihar][iphi]->Fill(pt, cos(n*(phi-RPPhi[11][ihar]))*cos(n*(phi-RPPhi[11][ihar])));
        }*/
      }
    } //track loop
    } //calFlag == 3
    } //calFlag > 0
 }//ihar

} //event loop

  return 0;
}   

//_____________________________________________________________________________________________________________________________
int EPAnaRun16alltree::End()
{
  cout << "End of EPAnaRun16alltree for Run " << RunNumber << endl;
  cout << "Total # of events = " << ievent << ", Total # of passed events = "<< jevent << endl;
  cout << "OutputFileName = " << OutputFileName << endl;

 if(calFlag == 3){
  if(d_outfile && d_outfile->IsOpen()) {
    d_outfile->cd();
    /*
    for(int icent=0;icent<ncent;icent++){
    for(int ibbcz=0;ibbcz<nbbcz;ibbcz++){
      for(int ihar=0;ihar<nhar;ihar++){
        for(int isub=0;isub<nsub;isub++){
        psi[icent][ibbcz][ihar][isub]->Write();
        psiFla[icent][ibbcz][ihar][isub]->Write();
          for(int ixy=0;ixy<nxy;ixy++){
            q[icent][ibbcz][ihar][isub][ixy]->Write();
            qRec[icent][ibbcz][ihar][isub][ixy]->Write();
          }
          for(int iord=0;iord<nord;iord++){
            sinflt[icent][ibbcz][ihar][isub][iord]->Write();
            cosflt[icent][ibbcz][ihar][isub][iord]->Write();
          }
        }
      }
    }
  }
  */
  
    for(int icent=0;icent<ncent;icent++){
    for(int ibbcz=0;ibbcz<nbbcz;ibbcz++){
      for(int ihar=0;ihar<nhar;ihar++){
        for(int isub=0;isub<nsub;isub++){
          for(int ixy=0;ixy<nxy;ixy++){
            q[icent][ibbcz][ihar][isub][ixy]->Write();
            qRec[icent][ibbcz][ihar][isub][ixy]->Write();
            }
        psi[icent][ibbcz][ihar][isub]->Write();
        psiFla[icent][ibbcz][ihar][isub]->Write();
        phiweight[icent][ibbcz][ihar][isub]->Write();
        phiweightbbc[icent][ibbcz][ihar][isub]->Write();

        }
      }
    }
    }

  for(int icent=0;icent<ncent;icent++){
      for(int ihar=0;ihar<nhar;ihar++){
          /*
          EPRCNTBBCS[icent][ihar]->Write();
          EPRCNTFVTX1S[icent][ihar]->Write();
          EPRCNTFVTX2S[icent][ihar]->Write();
          EPRBBCSFVTX2S[icent][ihar]->Write();
          */
          EPRBBCSFVTX1S[icent][ihar]->Write();
          EPRFVTX1NFVTX1S[icent][ihar]->Write();
          EPRBBCSFVTX1LS[icent][ihar]->Write();
          EPRBBCSFVTX2LS[icent][ihar]->Write();
          EPRBBCSFVTX3LS[icent][ihar]->Write();
          EPRBBCSFVTX4LS[icent][ihar]->Write();
          EPRBBCSFVTX1p2p3LS[icent][ihar]->Write();
          EPRBBCSFVTX1p2p4LS[icent][ihar]->Write();
          EPRBBCSFVTXtrkS[icent][ihar]->Write();
        /*
          EPRCNTBBCS[icent][ihar]->Write();
          EPRCNTFVTX1S[icent][ihar]->Write();
          EPRBBCSFVTX1S[icent][ihar]->Write();
          EPRCNTFVTX1trkS[icent][ihar]->Write();
          EPRCNTFVTX2trkS[icent][ihar]->Write();
          EPRBBCSFVTX1trkS[icent][ihar]->Write();
          EPRBBCSFVTX2trkS[icent][ihar]->Write();

          EPRCNTcBBCSc[icent][ihar]->Write();
          EPRCNTcFVTX1Sc[icent][ihar]->Write();
          EPRBBCScFVTX1Sc[icent][ihar]->Write();
          */
        for(int iphi=0;iphi<nphi;iphi++){
          /*
          vobsFVTX2S[icent][ihar][iphi]->Write();
          vobsFVTX2Ssq[icent][ihar][iphi]->Write();
          */
          vobsFVTX1LS[icent][ihar][iphi]->Write();
          vobsFVTX2LS[icent][ihar][iphi]->Write();
          vobsFVTX3LS[icent][ihar][iphi]->Write();
          vobsFVTX4LS[icent][ihar][iphi]->Write();
          vobsFVTX1S[icent][ihar][iphi]->Write();
          vobsFVTX1N[icent][ihar][iphi]->Write();
          vobsFVTX1p2p3LS[icent][ihar][iphi]->Write();
          vobsFVTX1p2p4LS[icent][ihar][iphi]->Write();
          vobsBBCS[icent][ihar][iphi]->Write();
          vobsFVTXtrkS[icent][ihar][iphi]->Write();

          vFVTX1LS[icent][ihar][iphi]->Write();
          vFVTX2LS[icent][ihar][iphi]->Write();
          vFVTX3LS[icent][ihar][iphi]->Write();
          vFVTX4LS[icent][ihar][iphi]->Write();
          vFVTX1S[icent][ihar][iphi]->Write();
          vFVTX1N[icent][ihar][iphi]->Write();
          vFVTX1p2p3LS[icent][ihar][iphi]->Write();
          vFVTX1p2p4LS[icent][ihar][iphi]->Write();
          vBBCS[icent][ihar][iphi]->Write();
          vFVTXtrkS[icent][ihar][iphi]->Write();

          vnFVTX1LS[icent][ihar][iphi]->Write();
          vnFVTX2LS[icent][ihar][iphi]->Write();
          vnFVTX3LS[icent][ihar][iphi]->Write();
          vnFVTX4LS[icent][ihar][iphi]->Write();
          vnFVTX1N[icent][ihar][iphi]->Write();
          vnFVTX1S[icent][ihar][iphi]->Write();
          vnFVTX1p2p3LS[icent][ihar][iphi]->Write();
          vnFVTX1p2p4LS[icent][ihar][iphi]->Write();
          vnBBCS[icent][ihar][iphi]->Write();
          vnFVTXtrkS[icent][ihar][iphi]->Write();
         /* 
          vobscBBCSc[icent][ihar][iphi]->Write();
          vobscBBCScsq[icent][ihar][iphi]->Write();
          vobscFVTX1Sc[icent][ihar][iphi]->Write();
          vobscFVTX1Scsq[icent][ihar][iphi]->Write();
          */
         /*
          vobsBBCS[icent][ihar][iphi]->Write();
          vobsBBCSsq[icent][ihar][iphi]->Write();
          vobsFVTX1S[icent][ihar][iphi]->Write();
          vobsFVTX1Ssq[icent][ihar][iphi]->Write();
          vobsFVTX1trkS[icent][ihar][iphi]->Write();
          vobsFVTX1trkSsq[icent][ihar][iphi]->Write();
          vobsFVTX2trkS[icent][ihar][iphi]->Write();
          vobsFVTX2trkSsq[icent][ihar][iphi]->Write();
          */
          
          }
        }
        }

/*
for(int istation=0;istation<4;istation++){
    hfvtx[istation]->Write();
}
    hbbc->Write();
    hfvtxtrk->Write();
    hbbcsnfvtxtrk->Write();
    */
//    hnfvtxtrk->Write();
    d_outfile->Close();
  }else{

    cout<<"ERROR: No output file set!"<<endl;

  }
 }
  return 0;
}

int EPAnaRun16alltree::SetcalFlag(int _flag){
    calFlag = _flag;
    ievent = 0;
    jevent = 0;
    return 0;
}

void EPAnaRun16alltree::Getrec(){
    for(int icent=0;icent<ncent;icent++){
    for(int ibbcz=0;ibbcz<nbbcz;ibbcz++){
      for(int ihar=0;ihar<nhar;ihar++){
        for(int isub=0;isub<nsub;isub++){
    meanx[icent][ibbcz][ihar][isub] = q[icent][ibbcz][ihar][isub][0]->GetMean();
    rmsx[icent][ibbcz][ihar][isub] = q[icent][ibbcz][ihar][isub][0]->GetRMS();
    meany[icent][ibbcz][ihar][isub] = q[icent][ibbcz][ihar][isub][1]->GetMean();
    rmsy[icent][ibbcz][ihar][isub] = q[icent][ibbcz][ihar][isub][1]->GetRMS();
                }
        }
        }
    }
}

void EPAnaRun16alltree::Getflt(){
    for(int icent=0;icent<ncent;icent++){
    for(int ibbcz=0;ibbcz<nbbcz;ibbcz++){
      for(int ihar=0;ihar<nhar;ihar++){
        for(int isub=0;isub<nsub;isub++){
            for(int iord=0; iord<nord; iord++) {
                cosfltarr[icent][ibbcz][ihar][isub][iord] = cosflt[icent][ibbcz][ihar][isub][iord]->GetMean();
                sinfltarr[icent][ibbcz][ihar][isub][iord] = sinflt[icent][ibbcz][ihar][isub][iord]->GetMean();
        }
      }
    }
    }
    }
}
