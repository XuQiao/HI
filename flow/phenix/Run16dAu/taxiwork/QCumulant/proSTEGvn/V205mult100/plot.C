void plot(){
  TFile *fp = TFile::Open("outP.root");

  int xcorr = 0;
  int iscorr = 0;


  TH1D* h = new TH1D("","",100,0,10);
  int color[12] = {1,2,4,7,8,5,1,2,4,7,8,5};
  int style[12] = {20,21,24,25,26,27,29,30,31,32,33,34};
  TCanvas *c1 = new TCanvas();
  c1->cd();
  int xcent = 0;
  int xhar = 1;
  int n = xhar + 1 + iharE;
  SetTitle(*h,"n-particle correlation", Form("<d%d>",n),"");
  SetXRange(*h,0,8);
  h->DrawCopy();
  TLegend *leg = new TLegend(0.5,0.7,0.7,0.9);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  for(int isub=0;isub<nsub;isub++){
    SetStyle(*grwd[xcent][isub][xhar],1.2,color[isub],style[isub]);
    grwd[xcent][isub][xhar]->Draw("Psame");
    leg->AddEntry(grwd[xcent][isub][xhar],Form("<d%d> cent:%d-%d, %s",n,0,5,choosesub(isub).Data()),"P");
  }
  leg->Draw("same");
  c1->Print(Form("d%d_cent%d.png",n,xcent));

  TCanvas *c2 = new TCanvas();
  c2->cd();
  int xcent = 0;
  int xhar = 1;
  int n = xhar + 1 + iharE;
  SetTitle(*h,"n-particle correlation", Form("<c%d>",n),"");
  SetXRange(*h,0,8);
  h->DrawCopy();
  TLegend *leg = new TLegend(0.5,0.7,0.7,0.9);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  for(int isub=0;isub<nsub;isub++){
    SetStyle(*grwc[xcent][isub][xhar],1.2,color[isub],style[isub]);
    grwc[xcent][isub][xhar]->Draw("Psame");
    leg->AddEntry(grwd[xcent][isub][xhar],Form("<c%d> cent:%d-%d, %s",n,0,5,choosesub(isub).Data()),"P");
  }
  leg->Draw("same");
  c2->Print(Form("c%d_cent%d.png",n,xcent));

  TCanvas *c3 = new TCanvas();
  c3->cd();
  int xcent = 0;
  int xhar = 1;
  int n = xhar + 1 + iharE;
  SetTitle(*h,"n-particle correlation", Form("<v%d>",n),"");
  SetXRange(*h,0,8);
  h->DrawCopy();
  TLegend *leg = new TLegend(0.5,0.7,0.7,0.9);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  for(int isub=0;isub<nsub;isub++){
    SetStyle(*grwv[xcent][isub][xhar],1.2,color[isub],style[isub]);
    grwv[xcent][isub][xhar]->Draw("Psame");
    leg->AddEntry(grwv[xcent][isub][xhar],Form("<v%d> cent:%d-%d, %s",n,0,5,choosesub(isub).Data()),"P");
  }
  leg->Draw("same");
  c3->Print(Form("v%d_cent%d.png",n,xcent));

  TCanvas *c4 = new TCanvas();
  c4->cd();
  int xcent = 0;
  int xhar = 1;
  int xcorr = 0;
  int n = xhar + 1 + iharE;
  TH1D* h = new TH1D("","",50,0,5);
  SetTitle(*h,"p_{T} (GeV/c)", Form("<d%d>",n),"");
  SetXRange(*h,0,5);
  h->DrawCopy();
  TLegend *leg = new TLegend(0.5,0.7,0.7,0.9);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  TLatex *tl = new TLatex();
  tl->SetNDC();
  tl->DrawLatex(0.4,0.7,Form("%d-particle correlation differential",2*(xcorr+1)));
  for(int isub=0;isub<nsub;isub++){
    SetStyle(*grwdpr[xcent][isub][xhar][xcorr],1.2,color[isub],style[isub]);
    grwdpr[xcent][isub][xhar][xcorr]->Draw("Psame");
    leg->AddEntry(grwdpr[xcent][isub][xhar][xcorr],Form("<d%d> cent:%d-%d, %s",n,0,5,choosesub(isub).Data()),"P");
  }
  leg->Draw("same");
  c4->Print(Form("dpr%d_cent%d.png",n,xcent));
  
  TCanvas *c5 = new TCanvas();
  c5->cd();
  int xcent = 0;
  int xhar = 1;
  int xcorr = 0;
  int n = xhar + 1 + iharE;
  SetTitle(*h,"p_{T} (GeV/c)", Form("<c%d>",n),"");
  SetXRange(*h,0,5);
  h->DrawCopy();
  TLegend *leg = new TLegend(0.5,0.7,0.7,0.9);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  TLatex *tl = new TLatex();
  tl->SetNDC();
  tl->DrawLatex(0.4,0.7,Form("%d-particle correlation differential",2*(xcorr+1)));
  for(int isub=0;isub<nsub;isub++){
    SetStyle(*grwcpr[xcent][isub][xhar][xcorr],1.2,color[isub],style[isub]);
    grwcpr[xcent][isub][xhar][xcorr]->Draw("Psame");
    leg->AddEntry(grwcpr[xcent][isub][xhar][xcorr],Form("<c%d> cent:%d-%d, %s",n,0,5,choosesub(isub).Data()),"P");
  }
  leg->Draw("same");
  c5->Print(Form("cpr%d_cent%d.png",n,xcent));
  
  TCanvas *c6 = new TCanvas();
  c6->cd();
  int xcent = 0;
  int xhar = 1;
  int xcorr = 0;
  int n = xhar + 1 + iharE;
  SetTitle(*h,"p_{T} (GeV/c)", Form("<v%d>",n),"");
  SetXRange(*h,0,5);
  h->DrawCopy();
  TLegend *leg = new TLegend(0.1,0.6,0.45,0.9);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetFillColor(0);
  TLatex *tl = new TLatex();
  tl->SetTextSize(0.04);
  tl->SetNDC();
  tl->DrawLatex(0.5,0.75,Form("cent: %d-%d",0,5));
  tl->DrawLatex(0.5,0.8,Form("%d-particle correlation differential",2*(xcorr+1)));
  for(int isub=0;isub<nsub;isub++){
    SetStyle(*grwvpr[xcent][isub][xhar][xcorr],1.2,color[isub],style[isub]);
    grwvpr[xcent][isub][xhar][xcorr]->Draw("Psame");
    leg->AddEntry(grwvpr[xcent][isub][xhar][xcorr],Form("<v%d>, %s",n,0,5,choosesub(isub).Data()),"P");
  }
  leg->Draw("same");
  c6->Print(Form("vpr%d_cent%d.png",n,xcent));
}
