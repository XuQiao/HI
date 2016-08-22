const double ptmin = 0.3;
const double ptmax = 6.0; 
const double ptbin[] = {0.3,0.4,0.7,1.2,1.7,2.2,2.7,3.4,4.4,5.4,6.0};
const int npt= 10;
const int trkbin[]={260,220};
const int nbin = 1;
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
