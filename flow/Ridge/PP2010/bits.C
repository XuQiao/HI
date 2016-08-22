    void bits() {
       TFile *file = TFile::Open("/store/user/qixu/flow/NewSTEG/pPbDataV205m100/vndata_50k_0.root");
       TTree *T = (TTree*)file->Get("tree");
       Long64_t nentries = T->GetEntries();
       printf("tree has %d entries\n",nentries);
       TBits *bits = new TBits(nentries); //see http://root.cern.ch/root/html534/TBits.html
       //now assume you want to read N=1000 entries out of nentries
       Int_t N = 100;
       Int_t i=0;
       TRandom3 rand(0);
       while(i<N) {
          Int_t j = rand.Uniform(0,nentries);
          if (bits->TestBitNumber(j)) continue;  //we have already seen this entry
          bits->SetBitNumber(j);
          i++;
       }
       //now loop on the N selected entries
       Int_t fbit = 0;
       for (i=0;i<N;i++) {
          j = bits->FirstSetBit(fbit);
          //read this Tree entry
          T->GetEntry(j);
          //process entry
          printf("reading entry %d\n",j);
          //.....
          fbit = j+1;
       }
    }
     


