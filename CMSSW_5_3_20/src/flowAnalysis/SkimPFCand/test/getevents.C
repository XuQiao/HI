void getevents(){
	int file;	TString name;
	for(int dataset=0;dataset<3;dataset++){
	int totevent=0;	
	if(dataset==0){	file = 221; name = "pPbReReco";}
	else if(dataset==1){ file = 556; name = "pPb";}
	else{	file = 703; name = "Pbp";}

	for(int i=0;i<file;i++){
		TFile *f = TFile::Open(Form("/lio/lfs/cms/store/user/qixu/flow/PACorrHM/skim/PFcandi/pbs/M300260/FlowHM%s/pPbskimAnaTree_%d.root",name.Data(),i));
		TTree *tr = (TTree*)f->Get("demo/PFTree");
		int nevent=tr->GetEntries();
		totevent+=nevent;
		f->Close();
		delete f;
	}

	TFile *fskim = TFile::Open(Form("/cms/store/user/qixu/flow/PACorrHM/skim/tracknormcpt03to6/multiM300260/FlowLYZHM%s/mergedTree.root",name.Data()));
	TTree *tree = (TTree*)fskim->Get("demo/TrackTree");
	int neventfromtrack = tree->GetEntries();
	fskim->Close();		

	int totevent1=0;
	for(int i=0;i<getnumber(Form("/home/xuq7/HI/flow/LYZ/v2/PFcandpt01to10/M300260/filelistcrab%s.dat",name.Data()));i++){
		string crabfilename = readline(Form("/home/xuq7/HI/flow/LYZ/v2/PFcandpt01to10/M300260/filelistcrab%s.dat",name.Data()),i);
		TFile *fcrab = TFile::Open(crabfilename.c_str());
		TTree *tr1 = (TTree*)fcrab->Get("demo/PFTree");
		int nevent1 = tr1->GetEntries();
		totevent1+=nevent1;
		fcrab->Close();
		delete fcrab;
	}

	cout<<"Pbs result for "<<name<<" = "<<totevent<<";Crab result = "<<totevent1<<"; track result = "<<neventfromtrack<<endl;
	}

	TFile *ftot = TFile::Open("/home/xuq7/HI/flow/LYZ/v2/tracknormcpt03to6/M300260/mergedV_Sum.root");
	TVectorD* Nevent = (TVectorD*)ftot->Get("Nevent");
	cout<<"totevent = "<<(*Nevent)[0]<<endl;
	ftot->Close();
}

int getnumber(char* name){
        std::ifstream backstory (name);
        std::string line;
        int nlines=0;
        if (backstory.is_open())
        {
                while (std::getline(backstory, line))
                        nlines++;
                backstory.close();
                return nlines;
        }
        else
        {
                std::cout << "Unable to open file" << std::endl << std::endl;
                return 0;
        }

}

std::string readline(char* name, int iline){
        std::ifstream backstory(name);
        std::string line;
        if (backstory.is_open())
                if(backstory.good()){
                        for(int i = 0; i < iline+1; ++i)
                           getline(backstory, line);
                }
        return line;
}

