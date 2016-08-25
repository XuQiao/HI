#include <iostream>
#include "Run15pc3dphidzcalib.dat"
void Drawpc3calibhisto(){
    double pt[50];
    for(int ipt = 0; ipt<50; ipt++){
        pt[ipt] = ipt * 0.1 + 0.05;
    }
    TGraph* pc3dphipt_arm0_pos = new TGraph(50,pt,dphishift[0][0]);
    TGraph* pc3dphipt_arm0_neg = new TGraph(50,pt,dphishift[0][1]);
    TGraph* pc3dphipt_arm1_pos = new TGraph(50,pt,dphishift[1][0]);
    TGraph* pc3dphipt_arm1_neg = new TGraph(50,pt,dphishift[1][1]);
    pc3dphipt_arm0_pos -> SetMarkerStyle(20);
    pc3dphipt_arm0_pos -> SetMarkerSize(1.2);
    pc3dphipt_arm0_pos -> Draw("AP");
    
    TGraph* pc3dzpt_arm0_pos = new TGraph(50,pt,dzshift[0][0]);
    TGraph* pc3dzpt_arm0_neg = new TGraph(50,pt,dzshift[0][1]);
    TGraph* pc3dzpt_arm1_pos = new TGraph(50,pt,dzshift[1][0]);
    TGraph* pc3dzpt_arm1_neg = new TGraph(50,pt,dzshift[1][1]);
    pc3dzpt_arm0_pos -> SetMarkerStyle(20);
    pc3dzpt_arm0_pos -> SetMarkerSize(1.2);
    pc3dzpt_arm0_pos -> Draw("AP");
    
    TGraph* pc3dphispt_arm0_pos = new TGraph(50,pt,dphisigma[0][0]);
    TGraph* pc3dphispt_arm0_neg = new TGraph(50,pt,dphisigma[0][1]);
    TGraph* pc3dphispt_arm1_pos = new TGraph(50,pt,dphisigma[1][0]);
    TGraph* pc3dphispt_arm1_neg = new TGraph(50,pt,dphisigma[1][1]);
    pc3dphispt_arm0_pos -> SetMarkerStyle(20);
    pc3dphispt_arm0_pos -> SetMarkerSize(1.2);
    pc3dphispt_arm0_pos -> Draw("AP");
    
    TGraph* pc3dzspt_arm0_pos = new TGraph(50,pt,dzsigma[0][0]);
    TGraph* pc3dzspt_arm0_neg = new TGraph(50,pt,dzsigma[0][1]);
    TGraph* pc3dzspt_arm1_pos = new TGraph(50,pt,dzsigma[1][0]);
    TGraph* pc3dzspt_arm1_neg = new TGraph(50,pt,dzsigma[1][1]);
    pc3dzspt_arm0_pos -> SetMarkerStyle(20);
    pc3dzspt_arm0_pos -> SetMarkerSize(1.2);
    pc3dzspt_arm0_pos -> Draw("AP");
    
}
