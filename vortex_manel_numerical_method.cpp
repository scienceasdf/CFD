#include <iostream>
#include<vector>
#include<cmath>
#include<fstream>
#include<string>
#include<E:\eigen\Eigen/Dense>
#include<E:\eigen\Eigen/LU>


int main()
{

    const double tPi=6.28319;

    /*std::vector<double> px;
    std::vector<double> py;
    px.reserve(100);
    py.reserve(100);*/


    double vInf=10.0;
    double alpha=.0*tPi/360.0;
    double airDen=1.225;

    double px[100],py[100],ctrlpx[100],ctrlpy[100],nx[100],ny[100],s[100];
    double beta[100];



    //At present here I use array instead of std::vector
    Eigen::MatrixXd A;
    Eigen::VectorXd gamma;
    Eigen::VectorXd v;

    double cj;
    int n=0;    //n means the number of the points, as points are from zero to (n-1)
    std::string foilName;
    std::ifstream in("e://SD7037.dat");
    if(!in.is_open()){
        std::cout << "Error opening file";
        exit (1);
    }
    in>>foilName;
    std::cout<<foilName<<"\n";
    while (!in.eof() )
    {
        in>>px[n]>>py[n];
        std::cout<<px[n]<<"\t"<<py[n]<<"\n";
        ++n;
    }
    --n;
    std::cout<<px[n-1]<<"\t"<<py[n-1]<<"\n";

    //std::vector<double> px,py,ctrlpx,ctrlpy,nx,ny,beta,gamma;

    A.conservativeResize(n,n);
    gamma.conservativeResize(n,1);
    v.conservativeResize(n,1);

    double dx=.0,dy=.0;
    for(int i=0;i<n-1;++i){
        ctrlpx[i]=(px[i]+px[i+1])/2.0;
        ctrlpy[i]=(py[i]+py[i+1])/2.0;
        nx[i]=py[i]-py[i+1];
        ny[i]=px[i+1]-px[i];        //nx, ny represents the x and y components of the normal vector
        beta[i]=std::atan2(nx[i],.0-ny[i])+alpha;
        v[i]=vInf*beta[i];
        s[i]=pow((py[i]-py[i+1])*(py[i]-py[i+1])+(px[i+1]-px[i])*(px[i+1]-px[i]),.5);

    }
    ctrlpx[n-1]=(px[0]+px[n-1])/2.0;
    ctrlpy[n-1]=(py[0]+py[n-1])/2.0;
    s[n-1]=pow((py[0]-py[n-2])*(py[0]-py[n-2])+(px[0]-px[n-2])*(px[0]-px[n-2]),.5);

    for(int i=0;i<n;++i){
        for(int j=0;j<n;++j){
            dx=ctrlpx[i]-ctrlpx[j];
            dy=ctrlpy[i]-ctrlpy[j];
            cj=1.0/(1.0+dy*dy/dx/dx);
            cj=cj-(dy/dx/dx*nx[i]+1.0/dx*ny[i]);
            A(i,j)=(i==j)?0:cj*s[j]/tPi;            //Maybe here A(i,j) use the reload of the operator ()
        }
    }
    for(int i=1;i<n;++i){
        A(10,i)=.0;
    }
    A(10,0)=1.0;    A(10,n-1)=1.0;

    gamma=A.inverse()*v;

    double L=.0;
    for(int i=0;i<n;++i){
        L+=gamma[i]*s[i];
        std::cout<<"\n"<<gamma[i]<<"\t"<<s[i];
    }
    L*=airDen*vInf;
    std::cout<<"\n\n"<<L;
    return 0;
}

