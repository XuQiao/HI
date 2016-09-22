void split(char *inFileList= "tree.lst"){
  //split file to different run #
  char cc[100];
  ofstream lstfile[10000];
  int tmprun[10000];
  
  for(int irun=0;irun<10000;irun++){
    tmprun[irun]=0;
  }

  //get the list of run
  ofstream runlist("Run16dAu20GeV.lst");

  ifstream filelist(inFileList);
  char cntfile[1000];

  int nrun=0;
  int ifile=0;
  while (filelist.getline(cntfile, 1000)) {    
    int run=0;
    ifile++;
    //if(ifile>10) continue;

    cout<<cntfile<<endl;

    int index=0;
    //get run number
    for(int i=0; i<200; i++){
      if(cntfile[i]=='r'&&cntfile[i+1]=='o'&&cntfile[i+2]=='o'&&cntfile[i+3]=='t'){
	cout<<"i= "<<i<<endl;
	index=53;
      }
    }

    //int index=104-12;

    for(int j=index;j<index+6;j++){
      run=run*10+(int)cntfile[j]-48;
    }
    //}
    //}
    cout<<run<<endl;

    bool runflag=true;
    for(int irun=0;irun<nrun;irun++){
      if(run==tmprun[irun]) {
	runflag=false;
	lstfile[irun]<<cntfile<<endl;
      }
    }
    if(runflag) {
      //split
      tmprun[nrun]=run;
      sprintf(cc,"%s%s%d%s","run-by-run/run","_",run,".lst");
      cout<<cc<<endl;
      lstfile[nrun].open(cc);
      lstfile[nrun]<<cntfile<<endl;
      //store the run #
      runlist<<run<<endl;
      nrun++;
    }

  }//end of while
  runlist.close();
}
