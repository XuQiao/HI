#define lam_efficiency_cxx
#include "lam_efficiency.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void lam_efficiency::Loop(int weight_choice)
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
  int trk_weight_choice, kin_weight_choice;
  if ( weight_choice >= 100 ) {
    trk_weight_choice = weight_choice - 100;
    kin_weight_choice = trk_weight_choice;
  } else {
    trk_weight_choice = weight_choice;
    kin_weight_choice = 0;
  }

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   // HARDCODE for testing
   //Long64_t nentries = 500000
   double p, pt, y, life, mass, p2, pz, energy;
   //   cout << "Check incoming values: " << lamEtaXmin << ", " << lamEtaBinWidth << ", " << lamptNbins << endl;
   cout << "nentries = " << nentries << endl;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      int ncand = (*lamv0MatchStatus).size();
      veematch->Fill(ncand);
      double wgt = ntrk_wgt(nLooseTracks, trk_weight_choice);

      for ( int i = 0; i<ncand; i++ ) {
	  //eta = fabs((*lamgenEta)[i]);
	  //if (eta > 2.5) continue;
	  // Calculate rapidity
	p2 = ((*lamgenp)[i]) * ((*lamgenp)[i]);
	pt = (*lamgenpT)[i];
	pz = sqrt( p2 - pt*pt );
	energy = sqrt( p2 + lamMassConst2 );
	y = 0.5*log((energy + pz) / (energy - pz));
	if ( y > 2.3 ) continue;
	life = (*lamgenCtau)[i];
	double y_wgt = lampt_wgt( pt, kin_weight_choice );
	double pt_wgt = lamy_wgt( y, kin_weight_choice );
	if (processType != 92 && processType != 93 && (*lamgenMass)[i] < 1.115682 ) {
	  lam_dNdy_gen->Fill(y,wgt*y_wgt);
	  lam_dNdpT_gen->Fill(pt,wgt*pt_wgt);
	  lam_dNdlife_gen->Fill(life,wgt*y_wgt);
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

      }
   }
}
