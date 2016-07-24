#define kshortrec_masspeak_cxx
#include "kshortrec_masspeak.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void kshortrec_masspeak::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L kshortrec.C
//      Root > kshortrec t
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
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   int ietabin, iptbin, ipbin, ilifebin;
   float p, pt, eta, lifetime, mass;
   //   cout << "Check incoming values: " << ksEtaXmin << ", " << ksEtaBinWidth << ", " << ksptNbins << endl;
   cout << "nentires = " << nentries << endl;


   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      //      cout << "nb, nbytes = " << nb << ", " << nbytes << endl;
      //      if ( nV0s > 0 ) cout << "number of v0s = " << nV0s << endl;
      nV0sH->Fill(nV0s);
      //      if ( priVtxChi2 <= 0.0 ) continue;
      int ncand = (*v0CandMass).size();
      numvee->Fill(ncand);
      if ( !trigTech34 || !trigHLTminBias || !trigScraping || !trigHF || !trigTech36_39 ) continue;
      ntrk->Fill(nLooseTracks);
      for ( int i = 0; i<ncand; i++ ) {
	mass = (*v0CandMass)[i];
	vee_nocuts->Fill( mass );
	if ( fabs((*v3dIpWrtPrimarySig)[i]) >= 3.0 ) continue;
	if ( kshort && fabs((*v0OtherCandMass)[i]-1.115683)<0.0075 ) continue;
	if ( !kshort && fabs((*v0OtherCandMass)[i]-0.497614)<0.020 ) continue;
	eta = fabs((*v0Eta)[i]);
	if ( eta >= 2.5 ) continue;
	vee_cut1->Fill( mass );
	if ( kshort ) {
	  if ( (*v0VtxSig)[i] > 10.0 ) vee_cut2->Fill(mass);
	  if ( (*thePosDau3DIpSig)[i] > 2.0 || (*theNegDau3DIpSig)[i] > 2.0 ) vee_cut3->Fill(mass);
	  if ( (*thePosDau3DIpSig)[i] > 1.0 && (*theNegDau3DIpSig)[i] > 1.0 ) vee_cut4->Fill(mass);
	  if ( (*posDauNhits)[i] > 3.5 && (*negDauNhits)[i] > 3.5 && (*v0VtxSig)[i] > 10.0) vee_cut5->Fill(mass);
	} else {
	  if ( (*v0VtxSig)[i] > 10.0 ) vee_cut2->Fill(mass);
	  if ( (*v0VtxSig)[i] > 10.0 && ( (*thePosDau3DIpSig)[i] > 2.0 || (*theNegDau3DIpSig)[i] > 2.0 ) ) vee_cut3->Fill(mass);
	  if ( (*v0VtxSig)[i] > 10.0 && (*thePosDau3DIpSig)[i] > 1.0 && (*theNegDau3DIpSig)[i] > 1.0 ) vee_cut4->Fill(mass);
	  if ( ( (*posDauNhits)[i] + (*negDauNhits)[i] ) > 10 && (*v0VtxSig)[i] > 10.0) vee_cut5->Fill(mass);
	}
      }
   }
}
