{
TString cmsswbase = getenv("CMSSW_BASE");
    if(cmsswbase.Length() > 0){
        cout << "Loading FW Lite Setup" << endl;
        gSystem->Load("libFWCoreFWLite.so");
        AutoLibraryLoader::enable();
        gSystem->Load("libDataFormatsFWLite.so");
        gSystem->Load("libDataFormatsPatCandidates.so");
    }

}
