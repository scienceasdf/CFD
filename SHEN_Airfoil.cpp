#include<cmath>#include<fstream>#include<iostream>#include<cstdlib>#include<iomanip>#include<string>#include<sstream>/*Chase Method, solving tridiagnol linear equations. */double* root(double* a, double* b, double* c, double* d, int n){    if(n==1){        double* x=new double[1];        x[0]=d[0]/b[0];        return x;    }    else{        double* alpha=new double[n-1];        double* beta;        beta=new double[n];        beta[0]=b[0];        for(int i=0;i<=n-2;++i){            alpha[i]=a[i]/beta[i];            beta[i+1]=b[i+1]-alpha[i]*(c[i]);        }        double* x=new double[n],*y=new double[n];        y[0]=d[0];        for(int i=1;i<=n-1;++i) y[i]=d[i]-alpha[i-1]*(y[i-1]);        x[n-1]=beta[n-1];        for(int i=n-1;i>=0;--i) x[i]=(y[i]-c[i]*x[i+1])/beta[i];        delete [] y;        delete [] alpha;        delete [] beta;        return x;    }}inline double aux1(double x){    return (.29690*pow(x,.5)-.35160*x*x+.284330*x*x*x-.10363*x*x*x*x);}inline double daux1(double x){    return (.14845/pow(x,.5)-.7032*x+.85299*x*x-.41452*x*x*x);}class airfoil{public:    airfoil(double var1, double var2, double var3);    double uPos(double x);    double lPos(double x);private:    double thick(double x);    double camber(double x);    double tk,tf,tb;    double p,t,q;};airfoil::airfoil(double var1, double var2, double var3){    p=var1;    t=var2;    q=var3;    std::cout<<p<<"\t"<<t<<"\t"<<q<<"\n";    double p1[]={5.0*t*q,1.0};    double p2[]={-5.0*t,-1.0,2.0*(p-q)/p/p};    double p3[]={1.0,(2.0*p*q-q*q)/p/p};    //double p4[]={-5.0*t*aux1(1.0),5.0*t*aux1(q),5.0*t*(daux1(q)-aux1(1.0))};    double p4[]={5.0*t*aux1(1.0),-5.0*t*aux1(q),-5.0*t*(daux1(q)-aux1(1.0))};    std::cout<<p2[0]<<"\t"<<p3[0]<<"\n";    std::cout<<p1[0]<<"\t"<<p2[1]<<"\t"<<p3[1]<<"\n";    std::cout<<.0<<"\t"<<p1[1]<<"\t"<<p2[2]<<"\n";    std::cout<<p4[0]<<"\t"<<p4[1]<<"\t"<<p4[2]<<"\n";    double* x;    x=root(p1,p2,p3,p4,3);    tk=x[0];    tb=x[1];    tf=x[2];    std::cout<<tk<<"\n"<<tb<<"\n"<<tf<<"\n\na";    x[0]=(5.0*t*aux1(1.0)+5.0*t*tk-tb);    std::cout<<x[0]<<"a\n";    delete [] x;}double airfoil::thick(double x){    return (5.0*t*(aux1(x)-.126*x));}double airfoil::camber(double x){    double y;    if(x<q){        y=tf/p/p*(2.0*p*x-x*x);    }    else{        y=-5.0*t*(aux1(x)+tk*x)+tb;    }    return y;}double airfoil::uPos(double x){    return (camber(x)+thick(x));}double airfoil::lPos(double x){    return (camber(x)-thick(x));}int main(){    double t,p,q;    //...    /*std::cin>>p;    std::cin>>t;    std::cin>>q;*/    int IP,IT,IQ;    p=20.0; t=8.0; q=25.0;    IP=p;    IT=t;    IQ=q;    std::string name("SHEN");    std::stringstream ss;    ss<<name;    ss<<"-";    ss<<IP;    ss<<"-"<<IT<<"-"<<IQ;    ss>>name;    std::cout<<name<<"\n";    std::cout<<aux1(1)<<"\n";    p*=.1; t*=.01; q*=.01;    airfoil foil1(p,t,q);    std::ofstream fout("E://data//"+name+".dat");    if(!fout.is_open()){         std::cout << "Error opening file";         exit (1);     }    fout<<name<<"\n";    fout<<std::setprecision(6);    double x=1.0, y=.0;    for(int i=0;i<100;++i){        x=exp(-i/10.0);        y=foil1.uPos(x);        fout<<x<<"\t"<<y<<"\n";    }    x=.0;    y=.0;    fout<<x<<"\t"<<y<<"\n";    for(int i=99;i>-1;--i){        x=exp(-i/10.0);        y=foil1.lPos(x);        fout<<x<<"\t"<<y<<"\n";    }    fout.close();    return 0;}
