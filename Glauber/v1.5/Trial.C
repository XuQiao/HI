Trial(){
class Object{
private:
	int x;
public:
	Object (int a1=0, int b1=1, int c1=2);
	virtual ~Object(){};
	void assign(int A){
		x=A;
		}
	int getvalue(){return x;}
};
Object::Object(int a1, int b1, int c1) : x(a1);	{}

Object obj;
//obj.assign();
cout<<"x="<<obj.getvalue()<<endl;
//<<"a="<<obj.a<<endl<<"b="<<obj.b<<endl<<"c="<<obj.c<<endl;
}