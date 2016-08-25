const double ptmin = 0.2;
const double ptmax = 10.0; 
const double ptbin[]={0.3,0.5,0.7,1.0,1.5,2.0,3.0};
const int npt= 6;
const int trkbin[]={185,180,175,170,165,160,155,150};
const int nbin = 7;
const int nn=2;
const int mm=1;
const double j01=2.404826;
const int ntheta = 5;
const int nstepr=200;
   TF1 *V2vsPt = new TF1("V2vsPt","((x/[0])^[1]/(1+(x/[2])^[3]))*(.00005+(1/x)^[4])",0.2,10);	
	
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
