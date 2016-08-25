#include <iomanip>
void EvSel(){
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);
ofstream fstr("EventSelPbP.txt");
const int Ncolor=7; const int N=5; const int NSel=9;
static const int Color[Ncolor]={
   1, 2, 4, 46,6, 7,8
};
static const int Style[Ncolor] = {
    20, 34, 33, 25,27, 28,24
};
//double vzbin[]={-30,}
//TString TrigName[N]={"Jet20","Jet40","Jet60","Jet80","Jet100","MB"};
TString TrigName[N]={"Jet20","Jet40","Jet60","Jet80","Jet100"};
//TString histoName[NSel]={"Ev","EvRun","EvRunHLT","EvRunHLTEvt","EvRunHLTEvtppr","EvRunHLTEvtpprNoi","EvRunHLTEvtpprNoivz","EvRunHLTEvtpprNoivzGps"};
TString histoName[NSel]={"Ev","EvRun","EvRunHLT","EvRunHLTBeam","EvRunHLTBeamppr","EvRunHLTBeampprHFp","EvRunHLTBeampprHFpHFn","EvRunHLTBeampprHFpHFnNoi","EvRunHLTBeampprHFpHFnNoivz"};
TString CutName[NSel]={"0. No cut","1. Run >211256\t","2. HLT_PAJet*_NoJetID_v1","3. pBeamScrapingFilter\t","4. pprimaryvertexFilter\t","5. phfPosFilter1\t","6. phfNegFilter1\t","7. pHBHENoiseFilter\t","8. |vertex.z| < 15\t"};
TFile *f[N];	TH1F* hEvSel[N][NSel];

const int Nbase=0;
fstr<<"\tCutName\t\t"<<"Number of Events"<<"\t"<<"fraction to the total\t"<<"fraction to the previous"<<endl;
for(int i=0;i<N;i++){
f[i]=TFile::Open(Form("/scratch/xuq7/RpA/GlobalEvent/EventSel/Kurtfile/PbP%sEvSel.root",TrigName[i].Data()));
fstr<<TrigName[i]<<":"<<endl;
for(int j=0;j<NSel;j++){
hEvSel[i][j]=(TH1F*)f[i]->Get(Form("h%s",histoName[j].Data()));
if(j>=Nbase){
if(j==Nbase){
if(i==N-1)
fstr<<CutName[j]<<"\t"<<setprecision(0)<<fixed<<hEvSel[i][j]->GetEntries()<<"\t\t"<<setprecision(4)<<fixed<<(double)hEvSel[i][j]->GetEntries()/hEvSel[i][Nbase]->GetEntries()*100<<"%"<<"\t\t\t"<<endl;
else
fstr<<CutName[j]<<"\t"<<setprecision(0)<<fixed<<hEvSel[i][j]->GetEntries()<<"\t\t\t"<<setprecision(4)<<fixed<<(double)hEvSel[i][j]->GetEntries()/hEvSel[i][Nbase]->GetEntries()*100<<"%"<<"\t\t\t"<<endl;
}
else
fstr<<CutName[j]<<"\t"<<setprecision(0)<<fixed<<hEvSel[i][j]->GetEntries()<<"\t\t\t"<<setprecision(4)<<fixed<<(double)hEvSel[i][j]->GetEntries()/hEvSel[i][Nbase]->GetEntries()*100<<"%"<<"\t\t\t"<<(double)hEvSel[i][j]->GetEntries()/hEvSel[i][j-1]->GetEntries()*100<<"%"<<endl;
}
else
fstr<<CutName[j]<<"\t"<<setprecision(0)<<fixed<<hEvSel[i][j]->GetEntries()<<endl;
}
fstr<<setprecision(4)<<fixed<<(double)hEvSel[i][NSel-1]->GetEntries()/hEvSel[i][2]->GetEntries()*100<<"%"<<endl;
fstr<<endl;
}


}
