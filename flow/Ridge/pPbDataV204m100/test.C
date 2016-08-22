#include "vector"
void test(){
vector<vector<double> >pVect;
for(int i=0;i<10;i++){
     //   pVector.clear();
vector<double> pVector;

    for(int j=0;j<5;j++){
        pVector.push_back(i*10+j);
    }
    cout<<"original container size = "<<pVector.size()<<endl;
    pVect.push_back(pVector);
}
cout<<"large container size = "<<pVect.size()<<endl;
for(int i=0;i<10;i++){
    vector<double> pVector_t = pVect[i];
    cout<<"small container size = "<<pVector_t.size()<<endl;
    for(int j=0;j<5;j++){
        cout<<j<<"\t"<<pVector_t[j]<<"\t";
    }
    cout<<endl;
}
}
