#include "../../RpA/Quality/root_setting.h"
void ComparePaperppcross(){
    gStyle->SetOptStat(kFALSE);
    //gStyle->SetOptStat(100);
    gStyle->SetErrorX(0);
    gStyle->SetLabelFont(70);
    const double deta[]={0., 0.5, 1.0, 1.5, 2.0, 2.5} ;
    int nptpoint[]={33, 30, 27, 24, 19};
    int netabin = sizeof(deta)/sizeof(Double_t)-1 ;
    const Double_t jetPtBin[]={56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790, 846};
    const int nJetPtBin = sizeof(jetPtBin)/sizeof(Double_t)-1 ;
    const int N  = netabin ;
    double RptCen[N][100], R7ptCen[N][100], xerrminus[N][100], xerrplus[N][100], Ryval[N][100],Rystatup[N][100], Rystatdown[N][100], Rysysup[N][100], Rysysdown[N][100],R7yval[N][100],R7ystatup[N][100], R7ystatdown[N][100], R7ysysup[N][100], R7ysysdown[N][100];
    ifstream infile;
    
    std::string inname ;
    std::string inname2 ;
    
    bool Save=1;
    int conesize=7;
if(conesize==5) TString ref1 = "QCD-10-011";//    "PRL107(2011)132001"
if(conesize==7) TString ref1 = "QCD-11-004";//    "PRD87 (2012)112002"
TString ref2 =  "SMP-13-002"; //"PRD90 (2014)072006 

    TH1F *Rjet[N];
    TH1F *Rjet1[N];
    TH1F *Rratio[N];
    for(int ieta=0; ieta<N; ieta++)
    {
        Rjet[ieta] = new TH1F(Form("RjetcsAbsEta%.f",deta[ieta+1]*10),Form("RjetcsAbsEta%.f",deta[ieta+1]*10), nJetPtBin,jetPtBin);
        Rjet[ieta]->Sumw2();
        Rjet1[ieta] = new TH1F(Form("RjetcsAbsEta%.f_1",deta[ieta+1]*10),Form("RjetcsAbsEta%.f_1",deta[ieta+1]*10), nJetPtBin,jetPtBin);
        Rjet1[ieta]->Sumw2();
        Rratio[ieta] = new TH1F(Form("RjetcsratioAbsEta%.f_1",deta[ieta+1]*10),Form("RjetcsratioAbsEta%.f_1",deta[ieta+1]*10), nJetPtBin,jetPtBin);
        Rratio[ieta]->Sumw2();
    }

    inname = Form("InclusiveJets_Table_AK%d_new.txt",conesize);
    //cout <<"input = " <<inname.c_str() <<endl ;
    //  FILE* infile = fopen(inname.c_str(), "r");
    infile.open(inname.c_str());
    std::string line;
    int np =0;
    //  infile.open(inname.c_str());
    // string line;
    Int_t row = 0;
    double v[1000][100];
    // cout << "going to reading files !!!" << endl ;
    //   cout << "infile = " << infile <<endl ;
    if(infile.fail()) {
        cerr << "unable to open file for reading" << endl;
        exit(1);
    }
    
    while(std::getline(infile, line) && row < 1000 )
    {
        int column = 0;
        //     cout << "inside loop now 2222222 !!!" << endl ;
        
        istringstream lineStream(line);
        while( lineStream >> v[row][column] && column < 100 ) column++;
        row++;
        //      cout << "row =" << row  << "column = " << column <<endl ;
    }
    infile.close();
    Int_t ieta=0;
    for(int irow = 0 ; irow < row ; irow++){
        if(v[irow][0]>1.5 || v[irow][3]> jetPtBin[nJetPtBin]) continue ;
        for(int ibin = 0 ; ibin < netabin ; ibin++){
            if(v[irow][0]==deta[ibin] && v[irow][1]==deta[ibin+1]) ieta=ibin ;
        }
        if(v[irow][2]==56) np=0 ;
        //cout << " ieta = " << ieta << " np = " << np <<endl ;
        RptCen[ieta][np]=v[irow][2]+(v[irow][3]-v[irow][2])/2.;
        Ryval[ieta][np]=v[irow][4];
        Rystatup[ieta][np]=(v[irow][5])*v[irow][4];
        Rystatdown[ieta][np]=(v[irow][6])*v[irow][4];
        Rysysup[ieta][np]=TMath::Sqrt(v[irow][7]*v[irow][7]+v[irow][10]*v[irow][10]+v[irow][12]*v[irow][12]+v[irow][14]*v[irow][14]+v[irow][16]*v[irow][16]+v[irow][18]*v[irow][18]+v[irow][20]*v[irow][20]+v[irow][22]*v[irow][22]+v[irow][24]*v[irow][24]+v[irow][26]*v[irow][26]+v[irow][28]*v[irow][28]+v[irow][30]*v[irow][30]+v[irow][32]*v[irow][32]+v[irow][34]*v[irow][34]+v[irow][36]*v[irow][36]+v[irow][38]*v[irow][38]+v[irow][40]*v[irow][40]+v[irow][42]*v[irow][42]+v[irow][44]*v[irow][44]);
        
        int point = -1 ;
        for(int ibin = 0 ; ibin < nJetPtBin ; ibin++){
            if(v[irow][2]==jetPtBin[ibin] && v[irow][3]==jetPtBin[ibin+1]) point=ibin ;
        }
        if(point <0) continue ;
        int hbin = Rjet[ieta]->FindBin(RptCen[ieta][np]);
//        if(sys=="SysUp")
//            Rjet[ieta]->SetBinContent(hbin, Ryval[ieta][point]+Rysysup[ieta][point]*Ryval[ieta][point]);
//        else if(sys=="SysDown")
//            Rjet[ieta]->SetBinContent(hbin, Ryval[ieta][point]-Rysysdown[ieta][point]*Ryval[ieta][point]);
//        else
            Rjet[ieta]->SetBinContent(hbin, Ryval[ieta][point]);
        Rjet[ieta]->SetBinError(hbin, Rystatup[ieta][point]);
        
        //      sys[np]=v2*TMath::Sqrt(2);
        np++;
        
    }

//---------------------------------------------------------------------------------------------------------
    
    for(int ieta=0; ieta<4; ieta++){
    inname = Form("R%dY%.fp%.f_%.fp%.f.txt",conesize,TMath::Floor(deta[ieta]),(deta[ieta]-TMath::Floor(deta[ieta]))*10,TMath::Floor(deta[ieta+1]),(deta[ieta+1]-TMath::Floor(deta[ieta+1]))*10);
    //cout <<"input = " <<inname.c_str() <<endl ;
    //  FILE* infile = fopen(inname.c_str(), "r");
    infile.open(inname.c_str());
    std::string line;
    int np =0;
    //  infile.open(inname.c_str());
    // string line;
    Int_t row = 0;
    double v[1000][100];
    double RptCen[N][100], R7ptCen[N][100], xerrminus[N][100], xerrplus[N][100], Ryval[N][100],Rystatup[N][100], Rystatdown[N][100], Rysysup[N][100], Rysysdown[N][100],R7yval[N][100],R7ystatup[N][100], R7ystatdown[N][100], R7ysysup[N][100], R7ysysdown[N][100];
    // cout << "going to reading files !!!" << endl ;
    //   cout << "infile = " << infile <<endl ;
    if(infile.fail()) {
        cerr << "unable to open file for reading" << endl;
        exit(1);
    }
    
    while(std::getline(infile, line) && row < 1000 )
    {
        int column = 0;
        //     cout << "inside loop now 2222222 !!!" << endl ;
        
        istringstream lineStream(line);
        while( lineStream >> v[row][column] && column < 100 ) column++;
        row++;
        //      cout << "row =" << row  << "column = " << column <<endl ;
    }
    infile.close();
    for(int irow = 0 ; irow < row ; irow++){
        if(v[irow][2] > jetPtBin[nJetPtBin]) continue ;
        if(conesize == 5 && v[irow][1]==56) np=0 ;
        if(conesize == 7 && v[irow][1]==114) np=5 ;
        //cout << " ieta = " << ieta << " np = " << np <<endl ;
        RptCen[ieta][np]=v[irow][1]+(v[irow][2]-v[irow][1])/2.;
        Ryval[ieta][np]=v[irow][3];
        Rystatup[ieta][np]=(v[irow][4]);
        Rystatdown[ieta][np]=(v[irow][5]);
        Rysysup[ieta][np]=TMath::Sqrt(v[irow][7]*v[irow][7]+v[irow][10]*v[irow][10]+v[irow][12]*v[irow][12]+v[irow][14]*v[irow][14]+v[irow][16]*v[irow][16]+v[irow][18]*v[irow][18]+v[irow][20]*v[irow][20]+v[irow][22]*v[irow][22]+v[irow][24]*v[irow][24]+v[irow][26]*v[irow][26]+v[irow][28]*v[irow][28]+v[irow][30]*v[irow][30]+v[irow][32]*v[irow][32]+v[irow][34]*v[irow][34]+v[irow][36]*v[irow][36]+v[irow][38]*v[irow][38]+v[irow][40]*v[irow][40]+v[irow][42]*v[irow][42]+v[irow][44]*v[irow][44]);
        int point = -1 ;
        for(int ibin = 0 ; ibin < nJetPtBin ; ibin++){
            if(v[irow][1]==jetPtBin[ibin] && v[irow][2]==jetPtBin[ibin+1]) point=ibin ;
        }
        if(point <0) continue ;
        int hbin = Rjet1[ieta]->FindBin(RptCen[ieta][np]);
//        if(sys=="SysUp")
//            Rjet[ieta]->SetBinContent(hbin, Ryval[ieta][point]+Rysysup[ieta][point]*Ryval[ieta][point]);
//        else if(sys=="SysDown")
//            Rjet[ieta]->SetBinContent(hbin, Ryval[ieta][point]-Rysysdown[ieta][point]*Ryval[ieta][point]);
//        else
            Rjet1[ieta]->SetBinContent(hbin, Ryval[ieta][point]);
            Rjet1[ieta]->SetBinError(hbin, Rystatup[ieta][point]);       
        //      sys[np]=v2*TMath::Sqrt(2);
        np++;
    }   
         //   Rratio[ieta]=(TH1F*)Rjet[ieta]->Clone(Form("R_%d",ieta));
         //   Rratio[ieta]->Divide(Rjet1[ieta]);
        for(int ibin = Rjet[ieta]->FindBin(56); ibin < Rjet[ieta]->FindBin(569); ibin ++){
            double bincenter = Rjet[ieta]->GetBinCenter(ibin);
            Rratio[ieta]->SetBinContent(ibin, Rjet1[ieta]->GetBinContent(Rjet1[ieta]->FindBin(bincenter))/Rjet[ieta]->GetBinContent(ibin));
            Rratio[ieta]->SetBinError(ibin,0);
        }
}
 
c1 = new TCanvas("c1"," ",1200,400);
c2 = new TCanvas("c2"," ",1200,400);
makeMultiPanelCanvas(c1,4,1,0,0,0.22,0.18,0.03);
makeMultiPanelCanvas(c2,4,1,0,0,0.22,0.18,0.03);
TH1F* hFrame1=new TH1F("","",1000,0,1000);
fixedFontHist(hFrame1,1.4,1.6);
hFrame1->SetTitle("");
TLegend *leg1=new TLegend(0.26,0.26,0.60,0.45);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetLineWidth(0);
leg1->SetTextSize(0.045);
TLatex T1;
T1.SetNDC();
T1.SetTextSize(0.065);
T1.SetTextFont(42);
TLine *l = new TLine(hFrame1->GetXaxis()->GetXmin(),1,hFrame1->GetXaxis()->GetXmax(),1);
l->SetLineStyle(2);

for(int ieta=0;ieta<4;ieta++){
	double ybase=0.15;	double xbase=0.28;
    if(ieta==0){
        hFrame1->GetYaxis()->SetTitle("d^{2}#sigma/dp_{T}dy");
        hFrame1->SetLabelSize(18,"Y");
    }
    else{
        hFrame1->SetLabelSize(0,"Y");
        hFrame1->GetYaxis()->SetTitle("");
    }
        hFrame1->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
    c1->cd(ieta+1)->SetLogy();	
   hFrame1->GetXaxis()->SetNdivisions(510);
   hFrame1->GetXaxis()->SetRangeUser(32,597);
   hFrame1->SetMinimum(1e-4);
   hFrame1->SetMaximum(1e9);
   hFrame1->DrawCopy();
        Rjet[ieta]->SetMarkerSize(1.5);
        Rjet[ieta]->SetMarkerStyle(20);
        Rjet[ieta]->SetMarkerColor(1);
        Rjet[ieta]->SetLineColor(1);
        Rjet1[ieta]->SetMarkerSize(1.5);
        Rjet1[ieta]->SetMarkerStyle(24);
        Rjet1[ieta]->SetMarkerColor(2);
        Rjet1[ieta]->SetLineColor(2);

        Rjet[ieta]->Draw("Psame");
        Rjet1[ieta]->Draw("Psame");
        if(ieta==0)
        T1.DrawLatex(0.3,0.85,Form("|y| < %.1f",deta[ieta+1]));
        else
        T1.DrawLatex(0.3,0.85,Form("%.1f < |y| < %.1f",deta[ieta],deta[ieta+1]));
        if(ieta==0){
        leg1->AddEntry(Rjet[ieta],Form("%s",ref1.Data()),"lp");
        leg1->AddEntry(Rjet1[ieta],Form("%s",ref2.Data()),"lp");
        leg1->Draw("same");
        T1.DrawLatex(0.7,0.25,Form("R = 0.%d",conesize));
        }
    c2->cd(ieta+1);
    hFrame1->GetYaxis()->SetTitle("pp 7TeV measurement");
    if(conesize==5){
   hFrame1->SetMinimum(0.80);
   hFrame1->SetMaximum(1.20);
    }
    if(conesize==7){
   hFrame1->SetMinimum(0.95);
   hFrame1->SetMaximum(1.05);
    }
   hFrame1->DrawCopy();
   Rratio[ieta]->SetMarkerSize(1.5);
   Rratio[ieta]->SetMarkerStyle(20);
   Rratio[ieta]->SetMarkerColor(1);
   Rratio[ieta]->SetLineColor(1);
   //getRidYError(Rratio[ieta]);
   Rratio[ieta]->Draw("Psame");
   TLegend *leg= new TLegend(*leg1);
   leg->Clear();
        leg->AddEntry(Rratio[ieta],Form("#frac{%s}{%s}",ref1.Data(),ref2.Data()),"lp");
        if(ieta==0)
        T1.DrawLatex(0.3,0.85,Form("|y| < %.1f",deta[ieta+1]));
        else
        T1.DrawLatex(0.3,0.85,Form("%.1f < |y| < %.1f",deta[ieta],deta[ieta+1]));
        if(ieta==0){
        leg->Draw("same");
        T1.DrawLatex(0.7,0.25,Form("R = 0.%d",conesize));
        }
   l->Draw("same");
    }

if(Save){
c1->Print(Form("ppR%dcs.png",conesize));
c1->Print(Form("ppR%dcs.pdf",conesize));
c2->Print(Form("ppR%dcsratio.png",conesize));
c2->Print(Form("ppR%dcsratio.pdf",conesize));
}
}
