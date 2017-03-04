#include<iostream>
#include<cmath>
#include<algorithm>
#include<cstdio>

const double dx=.1;
const double gamma=1.4;

int main()
{

    double pro_pt1[31];
    double pV_pt1[31];
    double pT_pt1[31];
    double pro_pt2[31];
    double pV_pt2[31];
    double pT_pt2[31];
    double pro_pt[31];
    double pV_pt[31];
    double pT_pt[31];
    double ro_[31];
    double V_[31];
    double T_[31];
    double ro[31];
    double V[31];
    double T[31];
    double step[31];
    double A[31];
    double dt;

    for(int i=0;i<31;++i){
        A[i]=1.0+2.2*(.1*i-1.5)*(.1*i-1.5);
        ro[i]=1.0-.03146*i;
        T[i]=1.0-.02314*i;
        V[i]=(.1+.109*i)*pow(T[i],.5);
    }

    /*for(int i=0;i<31;++i){
        std::cout<<A[i]<<"\t"<<ro[i]<<"\t"<<V[i]<<"\t"<<T[i]<<"\n";
    }*/
    // above blocks are checked

    for(int counter=0;counter<=1400;++counter){

        for(int i=0;i<31;++i){
            step[i]=.5*dx/(V[i]+pow(T[i],.5));
        }
        dt=*std::min_element(step,step+31);     // checked

        for(int i=0;i<30;++i){
            pro_pt1[i]=-V[i]*(ro[i+1]-ro[i])/dx-ro[i]*(V[i+1]-V[i])/dx-ro[i]*V[i]*(log(A[i+1])-log(A[i]))/dx;
            pV_pt1[i]=-V[i]*(V[i+1]-V[i])/dx-1.0/gamma*((T[i+1]-T[i])/dx+T[i]*(ro[i+1]-ro[i])/ro[i]/dx);
            pT_pt1[i]=-V[i]*(T[i+1]-T[i])/dx-(gamma-1.0)*T[i]*((V[i+1]-V[i])/dx+V[i]*(log(A[i+1])-log(A[i]))/dx);

            ro_[i]=ro[i]+pro_pt1[i]*dt;
            V_[i]=V[i]+pV_pt1[i]*dt;
            T_[i]=T[i]+pT_pt1[i]*dt;
        }

        //std::cout<<pro_pt1[15]<<"\t"<<pV_pt1[15]<<"\t"<<pT_pt1[15]<<"\t"<<ro_[15]<<"\t"<<V_[15]<<"\t"<<T_[15]<<"\n";
        // above blocks are checked

        for(int i=1;i<31;++i){
            pro_pt2[i]=-V_[i]*(ro_[i]-ro_[i-1])/dx-ro_[i]*(V_[i]-V_[i-1])/dx-ro_[i]*V_[i]*(log(A[i])-log(A[i-1]))/dx;
            pV_pt2[i]=-V_[i]*(V_[i]-V_[i-1])/dx-1.0/gamma*((T_[i]-T_[i-1])/dx+T_[i]*(ro_[i]-ro_[i-1])/ro_[i]/dx);
            pT_pt2[i]=-V_[i]*(T_[i]-T_[i-1])/dx-(gamma-1.0)*T_[i]*((V_[i]-V_[i-1])/dx+V_[i]*(log(A[i])-log(A[i-1]))/dx);

            pro_pt[i]=.5*(pro_pt1[i]+pro_pt2[i]);
            pV_pt[i]=.5*(pV_pt1[i]+pV_pt2[i]);
            pT_pt[i]=.5*(pT_pt1[i]+pT_pt2[i]);

            ro[i]=ro[i]+pro_pt[i]*dt;
            V[i]=V[i]+pV_pt[i]*dt;
            T[i]=T[i]+pT_pt[i]*dt;
        }
        //std::cout<<pro_pt2[15]<<"\t"<<pV_pt2[15]<<"\t"<<pT_pt2[15]<<"\n";//<<"\t"<<ro_[15]<<"\t"<<V_[15]<<"\t"<<T_[15]<<"\n";
        //std::cout<<pro_pt[15]<<"\t"<<pV_pt[15]<<"\t"<<pT_pt[15]<<"\t"<<ro[15]<<"\t"<<V[15]<<"\t"<<T[15]<<"\n";

        V[0]=2.0*V[1]-V[2];

        V[30]=2.0*V[29]-V[28];
        ro[30]=2.0*ro[29]-ro[28];
        T[30]=2.0*T[29]-T[28];
    }

    for(int i=0;i<31;++i){
        std::cout<<A[i]<<"\t"<<ro[i]<<"\t"<<V[i]<<"\t"<<T[i]<<"\n";
    }

    return 0;
}
