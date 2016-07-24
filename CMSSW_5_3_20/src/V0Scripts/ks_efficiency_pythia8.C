#define ks_efficiency_cxx
#include "ks_efficiency.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void ks_efficiency::Loop(int weight_choice)
{
//   In a ROOT session, you can do:
//      Root > .L ks_efficiency.C
//      Root > ks_efficiency t
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

  const Double_t ksMassConst = 0.497614;
  const Double_t ksMassConst2 = ksMassConst * ksMassConst;
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
   // HARDCODE to test
   //Long64_t nentries = 500000;
   double p, pt, y, life, mass, p2, pz, energy;
   //   cout << "Check incoming values: " << ksEtaXmin << ", " << ksEtaBinWidth << ", " << ksptNbins << endl;
   cout << "nentries = " << nentries << endl;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      int ncand = (*ksv0MatchStatus).size();
      veematch->Fill(ncand);
      double wgt = ntrk_wgt(nLooseTracks, trk_weight_choice);

      for ( int i = 0; i<ncand; i++ ) {
	if ( (*ksgenMass)[i] < 0.49763 ) continue;
	// Switching to rapidity
	//eta = fabs((*ksgenEta)[i]);
	//if (eta > 2.5) continue;
	p2 = ((*ksgenp)[i]) * ((*ksgenp)[i]);
	pt = (*ksgenpT)[i];
	pz = sqrt( p2 - pt*pt );
	energy = sqrt( p2 + ksMassConst2 );
	y = 0.5*log((energy + pz) / (energy - pz));
	if ( y > 2.3 ) continue;
	life = (*ksgenCtau)[i];
	double y_wgt = kspt_wgt( pt, kin_weight_choice );
	double pt_wgt = ksy_wgt( y, kin_weight_choice );
	if (processType != 103 && processType != 104 ) {
	  ks_dNdy_gen->Fill(y,wgt*y_wgt);
	  ks_dNdpT_gen->Fill(pt,wgt*pt_wgt);
	  ks_dNdlife_gen->Fill(life,wgt*y_wgt);
	}

	if ( (*ksv0MatchStatus)[i] != 1 ) continue;
	if ( !kstrigTech34 || !kstrigScraping 
	     || !kstrigHF || !kstrigTech36_39 || !trigHLTminBias
	     || priVtxIsFake) continue;
	if ( fabs((*ks3dIpWrtPrimarySig)[i]) >= 3.0 ) continue;
	if ( fabs((*ksrecoOtherCandMass)[i]-1.115683)<0.0075 ) continue;
	if ( (*ksrecoVtxSig)[i] <= 10.0 ) continue;
	if ( (*ksrecoNegDauNormChi2)[i] < 0.0 || (*ksrecoPosDauNormChi2)[i] < 0.0 ) continue;
	//if ( fabs((*ksrecoEta)[i]) > 2.5 ) continue;

	p2 = ((*ksrecop)[i]) * ((*ksrecop)[i]);
	pt = (*ksrecopT)[i];
	pz = (*ksrecopZ)[i];
	energy = sqrt( p2 + ksMassConst2 );
	y = fabs(0.5*log((energy + pz) / (energy - pz)));
	if( y > 2.3 ) continue;

	//eta = fabs((*ksrecoEta)[i]);
	life = (*ksrecoCtau)[i];
	ks_dNdy_rec->Fill(y,wgt*y_wgt);
	ks_dNdpT_rec->Fill(pt,wgt*pt_wgt);
	ks_dNdlife_rec->Fill(life,wgt*y_wgt);

      }
   }
}
