#include "Lattice.hpp"

void LatticeBoltzmann::Advection(void){
    int ixnew, iynew, iznew;
    for(int ix=0; ix<Lx;ix++){
        for(int iy=0;iy<Ly;iy++){
            for(int iz=0;iz<Lz;iz++){
                for(int r=0;r<2;r++){
                    f_0[nf0(ix,iy,iz,r)] = f_0_new[nf0(ix,iy,iz,r)];
                }
                for(int p=0;p<3;p++){
                    for(int i=0; i<4;i++){
                        ixnew=(ix+V[0][p][i]+Lx)%Lx; iynew=(iy+V[1][p][i]+Ly)%Ly; iznew=(iz+V[2][p][i]+Lz)%Lz;
                        for(int j=0;j<2;j++){
                            for(int r=0;r<2;r++){
                                f[nf(ixnew,iynew,iznew,p,i,j,r)] = f_new[nf(ix,iy,iz,p,i,j,r)];
                            }
                        }
                    }
                }
            }
        }
    }
}

void LatticeBoltzmann::Collision(void){
    vector3D Emed0, B0, Jmed0;
    double FEQ, FEQ0, rhoc, eps, mu;
    for(int ix=0; ix<Lx ; ix++){
        for(int iy=0; iy<Ly ; iy++){
            for(int iz=0; iz<Lz ; iz++){
                Emed0 = LatticeBoltzmann::E_med(ix,iy,iz,false);
                B0 = LatticeBoltzmann::B(ix,iy,iz,false);
                Jmed0 = LatticeBoltzmann::J_med(ix,iy,iz,false);
                rhoc = LatticeBoltzmann::rho_c(ix,iy,iz,false);
                eps = epsr(ix,iy,iz); mu = mur(ix,iy,iz);

                FEQ0 = LatticeBoltzmann::feq0(rhoc);
                for(int r=0;r<2;r++){
                    f_0_new[nf0(ix,iy,iz,r)] = f_0[nf0(ix,iy,iz,r)] - (1.0/tau)*(f_0[nf0(ix,iy,iz,r)] - FEQ0);
                }
                for(int p=0;p<3;p++){
                    for(int i=0; i<4;i++){
                        for(int j=0;j<2;j++){
                            for(int r=0;r<2;r++){
                                FEQ = LatticeBoltzmann::feq(Jmed0,Emed0,B0,eps,mu,p,i,j,r);
                                f_new[nf(ix,iy,iz,p,i,j,r)] = f[nf(ix,iy,iz,p,i,j,r)] - (1.0/tau)*(f[nf(ix,iy,iz,p,i,j,r)]-FEQ);
                            }
                        }
                    }
                }
            }
        }
    }
}

void LatticeBoltzmann::Start(void)
{
    double B0norm = E0/c_light;
    vector3D Emed0, B0, Jmed0;
    Jmed0.load(0,0,0);
    for (int ix = 0; ix < Lx; ix++) //for each cell
        for (int iy = 0; iy < Ly; iy++)
            for (int iz = 0; iz < Lz; iz++){
                double pulse = exp(-alpha*(iz-z0)*(iz-z0));
                double eps = epsr(ix,iy,iz), mu = mur(ix,iy,iz);
                Emed0.load(E0*pulse,0,0); B0.load(0,B0norm*pulse,0);
                for (int r = 0; r < 2; r++){
                    f_0[nf0(ix,iy,iz,r)] = feq0(0);
                }
                for (int p = 0; p < 3; p++){ //in each plane
                    for (int i = 0; i < 4; i++){// at each direcction
                        for (int j = 0; j < 2; j++){
                            for (int r = 0; r < 2; r++){
                                f[nf(ix,iy,iz,p,i,j,r)] = feq(Jmed0,Emed0,B0,eps,mu,p,i,j,r);
                            }
                        }
                    }
                }
        }
}

void LatticeBoltzmann::Print(void)
{
    int ix,iy,iz;
    for (ix = 0; ix < Lx; ix++)
    {
        for (iy = 0; iy < Ly; iy++){

        for (iz=0;iz<Lz;iz++){
            vector3D E0 = E_med(ix,iy,iz,false);
            vector3D B0 = B(ix,iy,iz,false);

            std::cout<<iz<<" "<<E0.x()<<" "<<B0.y()<<" "<<U(ix,iy,iz,false)<<" "<<epsr(ix,iy,iz)<<std::endl;

            }}
        std::cout<<std::endl;
    }
}

void LatticeBoltzmann::Peaks(double &Et, double &Er)
{
    Et = 0; Er = 0;
    for (int iz = 0; iz < Lz; iz++){
        double Ex = E_med(0,0,iz,false).x();
        if (iz > Lz/2 && Ex > Et) Et = Ex;
        if (iz < Lz/2 && Ex < Er) Er = Ex;
    }
}
