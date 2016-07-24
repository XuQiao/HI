double partition(double array[],int low, int high){

double pivot=array[low];
while(low<high){
while(low<high&&array[high]>=pivot)
high--;
array[low]=array[high];
while(low<high&&array[low]<=pivot)
low++;
array[high]=array[low];
}
array[low]=pivot;
return low;
}

double findkth(double array[],int low,int high, int k){

double kth=0;
if(low == high) kth = array[low];
else{
 double pivot = array[low];
 int splitpoint = partition(array,low,high);
if(k-1 == splitpoint)   kth = pivot;
else if(k-1<splitpoint) kth = findkth(array,low,splitpoint-1,k);
else if(k-1>splitpoint) kth = findkth(array,splitpoint+1,high,k);
}


return kth;
}

double findpoint(TH1* h, double cent){
	int ibin;
	int binNX = -1;
	int binNY = -1;
	double goal=cent*h->Integral(0,binNX);
	//double goal=cent*h->Integral(0,binNX,0,binNY);
	//int xbinmin=1;//for some histograms there are bins in the negative multiplicity area due to binning issue
	for(ibin=1;ibin<h->GetNbinsX();ibin++){
		//if(h->GetBinLowEdge(ibin+1)<0) {xbinmin++; continue;}
		if(h->Integral(ibin,binNX)<=goal){
		//if(h->Integral(ibin,binNX,ibin,binNY)<=goal){
			//cout<<h->GetXaxis()->GetBinLowEdge(ibin)<<endl;
		break;
		}
	}
	//if(ibin==h->GetNbinsX()) return -1;
	if(fabs(h->Integral(ibin,binNX)-goal)>fabs(h->Integral(ibin-1,binNX)-goal))
	//if(fabs(h->Integral(ibin,binNX,ibin,binNY)-goal)>fabs(h->Integral(ibin-1,binNX,ibin-1,binNY)-goal))
		return h->GetXaxis()->GetBinLowEdge(ibin-1);
	else
		return h->GetXaxis()->GetBinLowEdge(ibin);
}
