#define lam_efficiency_cxx
#include "lam_efficiency.h"
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
#include "constants.h"
#include "ntrk_wgt.C"
#include "lampt_wgt.C"
#include "lamy_wgt.C"

void lam_efficiency::Loop(int weight_choice,
			  TH1F* veematch,
			  TH1F* lam_dNdy_gen,
			  TH1F* lam_dNdy_rec,
			  TH1F* lam_dNdpT_gen,
			  TH1F* lam_dNdpT_rec,
			  TH1F* lam_dNdlife_gen,
			  TH1F* lam_dNdlife_rec,
			  TH1F* lamReg_dNdy_gen,
			  TH1F* lamReg_dNdy_rec,
			  TH1F* lamReg_dNdpT_gen,
			  TH1F* lamReg_dNdpT_rec,
			  TH1F* lamReg_dNdlife_gen,
			  TH1F* lamReg_dNdlife_rec,
			  TH1F* lamBar_dNdy_gen,
			  TH1F* lamBar_dNdy_rec,
			  TH1F* lamBar_dNdpT_gen,
			  TH1F* lamBar_dNdpT_rec,
			  TH1F* lamBar_dNdlife_gen,
			  TH1F* lamBar_dNdlife_rec)
{
//   In a ROOT session, you can do:
//      Root > .L lam_efficiency.C
//      Root > lam_efficiency t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

  const Double_t lamMassConst = 1.115683;
  const Double_t lamMassConst2 = lamMassConst * lamMassConst;
  int trk_weight_choice, kin_weight_choice, part_weight_choice, anti_weight_choice;
  if ( weight_choice >= 100 ) {
    trk_weight_choice = weight_choice - 100;
    kin_weight_choice = trk_weight_choice;
    part_weight_choice = 1000 + kin_weight_choice;
    anti_weight_choice = 2000 + kin_weight_choice;
  } else {
    trk_weight_choice = weight_choice;
    kin_weight_choice = 0;
    part_weight_choice = 0;
    anti_weight_choice = 0;
  }

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   // HARDCODE for testing
   //   Long64_t nentries = 500000;
   double p, pt, y, life, mass, p2, pz, energy, wgt, 
     y_wgt, yReg_wgt, yBar_wgt, pt_wgt, ptReg_wgt, ptBar_wgt;
   int ncand;
   bool lambda;
   //   cout << "Check incoming values: " << lamEtaXmin << ", " << lamEtaBinWidth << ", " << lamptNbins << endl;
   cout << "nentries = " << nentries << endl;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      ncand = (*lamv0MatchStatus).size();
      veematch->Fill(ncand);
      wgt = ntrk_wgt(nLooseTracks, trk_weight_choice);

      for ( int i = 0; i<ncand; i++ ) {
	if ( (*lamgenMass)[i] > 1.115682 ) continue;
	  //eta = fabs((*lamgenEta)[i]);
	  //if (eta > 2.5) continue;
	  // Calculate rapidity
	p2 = ((*lamgenp)[i]) * ((*lamgenp)[i]);
	pt = (*lamgenpT)[i];
	pz = sqrt( p2 - pt*pt );
	energy = sqrt( p2 + lamMassConst2 );
	y = 0.5*log((energy + pz) / (energy - pz));
	if ( y > 2.3 ) continue;
	lambda = (*lamgenPDG)[i] > 0;
	life = (*lamgenCtau)[i];
	y_wgt = lampt_wgt( pt, kin_weight_choice );
	yReg_wgt = lampt_wgt( pt, part_weight_choice );
	yBar_wgt = lampt_wgt( pt, anti_weight_choice );
	pt_wgt = lamy_wgt( y, kin_weight_choice );
	ptReg_wgt = lamy_wgt( y, part_weight_choice );
	ptBar_wgt = lamy_wgt( y, anti_weight_choice );
	if (processType != 92 && processType != 93 && processType != 103 && processType != 104 ) {
	  lam_dNdy_gen->Fill(y,wgt*y_wgt);
	  lam_dNdpT_gen->Fill(pt,wgt*pt_wgt);
	  lam_dNdlife_gen->Fill(life,wgt*y_wgt);
	  if ( lambda ) {
	    lamReg_dNdy_gen->Fill(y,wgt*yReg_wgt);
	    lamReg_dNdpT_gen->Fill(pt,wgt*ptReg_wgt);
	    lamReg_dNdlife_gen->Fill(life,wgt*yReg_wgt);
	  } else {
	    lamBar_dNdy_gen->Fill(y,wgt*yBar_wgt);
	    lamBar_dNdpT_gen->Fill(pt,wgt*ptBar_wgt);
	    lamBar_dNdlife_gen->Fill(life,wgt*yBar_wgt);
	  }
	}

	if ( (*lamv0MatchStatus)[i] != 1 ) continue;
	if ( !lamtrigTech34 || !lamtrigScraping 
	     || !lamtrigHF || !lamtrigTech36_39 || !trigHLTminBias
	     || priVtxIsFake ) continue;
	if ( fabs((*lam3dIpWrtPrimarySig)[i]) >= 3.0 ) continue;
	if ( fabs((*lamrecoOtherCandMass)[i]-0.497614)<0.020 ) continue;
	if ( (*lamrecoVtxSig)[i] <= 10.0 ) continue;
	if ( (*lamrecoNegDauNormChi2)[i] < 0.0 || (*lamrecoPosDauNormChi2)[i] < 0.0 ) continue;
	//if ( fabs((*lamrecoEta)[i]) > 2.5 ) continue;

	//eta = fabs((*lamrecoEta)[i]);
	p2 = ((*lamrecop)[i]) * ((*lamrecop)[i]);
	pt = (*lamrecopT)[i];
	pz = (*lamrecopZ)[i];
	energy = sqrt( p2 + lamMassConst2 );
	y = fabs(0.5*log((energy + pz) / (energy - pz)));
	if( y > 2.3 ) continue;
	life = (*lamrecoCtau)[i];
	lam_dNdy_rec->Fill(y,wgt*y_wgt);
	lam_dNdpT_rec->Fill(pt,wgt*pt_wgt);
	lam_dNdlife_rec->Fill(life,wgt*y_wgt);
	if ( lambda ) {
	  lamReg_dNdy_rec->Fill(y,wgt*y_wgt);
	  lamReg_dNdpT_rec->Fill(pt,wgt*pt_wgt);
	  lamReg_dNdlife_rec->Fill(life,wgt*y_wgt);
	} else {
	  lamBar_dNdy_rec->Fill(y,wgt*y_wgt);
	  lamBar_dNdpT_rec->Fill(pt,wgt*pt_wgt);
	  lamBar_dNdlife_rec->Fill(life,wgt*y_wgt);
	}
      }
   }
}
