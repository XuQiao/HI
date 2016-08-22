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
	for(ibin=1;ibin<h->GetNbinsX();ibin++){
		double goal=cent*h->Integral();
		if(h->Integral(ibin,h->GetNbinsX())<=goal){
			//cout<<h->GetXaxis()->GetBinLowEdge(ibin)<<endl;
		break;
		}
	}
	return h->GetXaxis()->GetBinLowEdge(ibin);
}
