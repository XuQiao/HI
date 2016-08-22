void getevents(){
	int totevent=0;
	for(int i=0;i<221;i++){
		TFile *f = TFile::Open(Form("/lio/lfs/cms/store/user/qixu/flow/PACorrHM/skim/PFcandi/tracktest/pbs/M300260/FlowHMpPbReReco/pPbskimAnaTree_%d.root",i));
		TTree *tr = (TTree*)f->Get("demo/PFTree");
		int nevent=tr->GetEntries();
		totevent+=nevent;
		cout<<i<<" file: "<<nevent<<endl;
	}
	cout<<totevent<<endl;
}
