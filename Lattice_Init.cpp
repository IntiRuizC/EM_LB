#include "Lattice.hpp"
LatticeBoltzmann::LatticeBoltzmann(void){
    v0.load(0,0,0);

    V[0][0][0]=1;  V[0][1][0]=1;  V[0][2][0]=0;
    V[1][0][0]=1;  V[1][1][0]=0;  V[1][2][0]=1;
    V[2][0][0]=0;  V[2][1][0]=1;  V[2][2][0]=1;

    V[0][0][1]=-1; V[0][1][1]=-1; V[0][2][1]=0;
    V[1][0][1]=1;  V[1][1][1]=0;  V[1][2][1]=-1;
    V[2][0][1]=0;  V[2][1][1]=1;  V[2][2][1]=1;

    V[0][0][2]=-1; V[0][1][2]=-1; V[0][2][2]=0;
    V[1][0][2]=-1; V[1][1][2]=0;  V[1][2][2]=-1;
    V[2][0][2]=0;  V[2][1][2]=-1; V[2][2][2]=-1;

    V[0][0][3]=1;  V[0][1][3]=1;  V[0][2][3]=0;
    V[1][0][3]=-1; V[1][1][3]=0;  V[1][2][3]=1;
    V[2][0][3]=0;  V[2][1][3]=-1; V[2][2][3]=-1;

    for(int p=0;p<3;p++){
        for(int i=0;i<4;i++){
            v[p][i].load(V[0][p][i],V[1][p][i],V[2][p][i]);
        }
    }

    e0.load(0,0,0);
    for(int p=0;p<3;p++){
        for(int i=0; i<4;i++){
            e[p][i][0] = 0.5*v[p][(i+1)%4];
            e[p][i][1] = 0.5*v[p][(i+3)%4];
        }
    }

    b0.load(0,0,0);
    for(int p=0;p<3;p++){
        for(int i=0;i<4;i++){
            for(int j=0;j<2;j++){
	            b[p][i][j]=(v[p][i]^e[p][i][j]);
            }
        }
    }
}

double LatticeBoltzmann::epsr(int ix, int iy, int iz){
    return 0.5*(eps1+eps2) + 0.5*(eps2-eps1)*tanh(iz-Lz/2.0);
}

double LatticeBoltzmann::mur(int ix, int iy, int iz){
    return 1.0;
}

double LatticeBoltzmann::sigma(int ix, int iy, int iz){
    return 0.0;
}

double LatticeBoltzmann::rho_c(int ix, int iy, int iz, bool Usenew){
    double sum = 0.0;
    if (Usenew){
        sum += f_0_new[nf0(ix,iy,iz,0)];
        for (int p=0; p<3;p++){
            for (int i=0; i<4;i++){
                for (int j=0; j<2;j++){
                    sum += f_new[nf(ix,iy,iz,p,i,j,0)];
                }
            }
        }
    }
    else{
        sum += f_0[nf0(ix,iy,iz,0)];
        for (int p=0; p<3;p++){
            for (int i=0; i<4;i++){
                for (int j=0; j<2;j++){
                    sum += f[nf(ix,iy,iz,p,i,j,0)];
                }
            }
        }
    }
    return sum;
}

vector3D LatticeBoltzmann::D(int ix, int iy, int iz, bool Usenew){
    vector3D sum;
    sum.load(0,0,0);
    if (Usenew){
        for(int p=0;p<3;p++){
            for(int i=0; i<4;i++){
                for(int j=0;j<2;j++){
                    sum+=f_new[nf(ix,iy,iz,p,i,j,0)]*e[p][i][j];
                }
            }
        }
    }
    else{
        for(int p=0;p<3;p++){
            for(int i=0; i<4;i++){
                for(int j=0;j<2;j++){
                    sum+=f[nf(ix,iy,iz,p,i,j,0)]*e[p][i][j];
                }
            }
        }
    }
    return sum;
}

vector3D LatticeBoltzmann::E(int ix, int iy, int iz, bool Usenew){
    return LatticeBoltzmann::D(ix,iy,iz,Usenew)/epsr(ix,iy,iz);
}

vector3D LatticeBoltzmann::B(int ix, int iy, int iz, bool Usenew){
    vector3D sum;
    sum.load(0,0,0);
    if (Usenew){
        for(int p=0;p<3;p++){
            for(int i=0; i<4;i++){
                for(int j=0;j<2;j++){
                    sum+=f_new[nf(ix,iy,iz,p,i,j,1)]*b[p][i][j];
                }
            }
        }
    }
    else{
        for(int p=0;p<3;p++){
            for(int i=0; i<4;i++){
                for(int j=0;j<2;j++){
                    sum+=f[nf(ix,iy,iz,p,i,j,1)]*b[p][i][j];
                }
            }
        }
    }
    return sum;
}

vector3D LatticeBoltzmann::J_med(int ix, int iy, int iz, bool use_new){
    double sig = sigma(ix,iy,iz), eps = epsr(ix,iy,iz);
    return (sig/(1+0.25*mu0*sig/eps))*LatticeBoltzmann::E(ix,iy,iz,use_new);
}

vector3D LatticeBoltzmann::E_med(int ix, int iy, int iz, bool use_new){
    vector3D E0, J0;
    E0 = LatticeBoltzmann::E(ix,iy,iz,use_new);
    J0 = LatticeBoltzmann::J_med(ix,iy,iz,use_new);
    return E0 - (0.25*mu0/epsr(ix,iy,iz))*J0;
}

double LatticeBoltzmann::U(int ix, int iy, int iz, bool use_new){
    vector3D E0 = LatticeBoltzmann::E_med(ix,iy,iz,use_new), B0 = LatticeBoltzmann::B(ix,iy,iz,use_new);
    return 0.5*(epsr(ix,iy,iz)*E0.norm2() + B0.norm2()/mur(ix,iy,iz));
}

double LatticeBoltzmann::feq(vector3D &Jmed, vector3D &Emed, vector3D &B, double eps, double mu, int p, int i, int j, int r){
    if (r==0)
        return 0.0625*v[p][i]*Jmed + 0.25*eps*Emed*e[p][i][j] + 0.125/mu*B*b[p][i][j];
    return 0.0625*v[p][i]*Jmed + 0.25*Emed*e[p][i][j] + 0.125*B*b[p][i][j];
}

double LatticeBoltzmann::feq0(double rhoc){
    return rhoc;
}
