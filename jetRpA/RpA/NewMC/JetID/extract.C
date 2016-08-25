TVectorD* extract(double x){
        int N=6;
	x=x-0.5;
	if(x<0) exit();
	else if(x==0) int n=N;
	else int n  = (int)TMath::Log2(x)+1;
	n=(n>N?n:N);
	TVectorD* x2 = new TVectorD(); TVectorD* x2_ = new TVectorD();
	x2->ResizeTo(n);x2_->ResizeTo(n);
	x2->Zero();
	(*x2)[n-1] = (int)x%2;
	(*x2_)[n-1] = (int)x/2;
	for(int i=1; i<n; i++){
	(*x2)[n-i-1] = (int)(*x2_)[n-i]%2;
	(*x2_)[n-i-1] = (int)(*x2_)[n-i]/2;
	}
	return x2;
}
