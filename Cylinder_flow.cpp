/// \author SHEN Weihong
/// \date 2017/10/11
/// These codes are for solving laplace PDE,
/// here is for solving the flow function of a flow over a cylinder


#include <iostream>
#include <cmath>
#include <fstream>
#include <random>
#include <vector>

const int r_slices = 350;
const int theta_slices = 100;
const double r_ub = 350.0;
const double dtheta = 2.0 * 3.141592653589793 / theta_slices;
const double dr = (r_ub - 1.0) / r_slices;

int main()
{
    //double psi[r_slices+1][theta_slices+1];   //psi = psi(r,theta)
    std::vector<std::vector<double>> psi(r_slices+1);
    for(int i=0; i<r_slices+1; ++i){
        psi[i].resize(theta_slices+1);
    }
    double r = .0;

    //double theta = .0;

    std::random_device rd;
    std::default_random_engine gen(rd());
    std::uniform_real_distribution<> dis(-30.0,30.0);

    for(int i=0; i<r_slices+1; ++i){
        for(int j=0; j<theta_slices+1; ++j){
            psi[i][j] = dis(gen);
        }
    }

    /*for(int i=0; i<101; ++i){
        psi[0][i] = .0;
        psi[i][0] = .0;
        psi[i][100] = .0;
        r = i * dr + 1.0;
        psi[100][i] = 500.0 * sin (static_cast<double>(i) * dtheta);
    }*/

    for(int i=0; i<r_slices+1; ++i){
        psi[i][0]=.0;
        psi[i][theta_slices]=.0;
    }

    for(int i=0; i<theta_slices+1; ++i){
        psi[0][i] = .0;
        psi[r_slices][i] = 5.0 * r_ub * sin (static_cast<double>(i) * dtheta);
    }

    double dmax = 5.0;

    while(dmax>1e-2){
        dmax = .0;
        for(int i=1; i<r_slices; ++i){
            for(int j=1; j<theta_slices; ++j){
                double storage = psi[i][j];
                r = i * dr + 1.0;
                //theta = i * dtheta;
                psi[i][j] = (psi[i+1][j] + psi[i-1][j]) / dr / dr +
                        (psi[i][j+1] + psi[i][j-1]) / dtheta / dtheta / r / r
                        + (psi
                           [i+1][j]-psi[i-1][j]) / 2.0 / r / dr;
                psi[i][j] /= 2.0 / dr / dr + 2.0 / dtheta / dtheta / r / r;
                psi[i][j] = storage + 1.4 * (psi[i][j] - storage);
                double d = fabs(psi[i][j] - storage);
                dmax = (d>dmax)? d : dmax;
            }
        }

        //std::cout<<dmax<<"\n";
    }

    int i1 = 15, i2 = 64;
    std::cout << psi[i1][i2]<<"\n";
    double theta = i2 * dtheta;
    r = i1 * dr + 1.0;
    std::cout << 5.0 * r * sin(theta) * (1.0 - 1.0 / r / r);

    return 0;

}
