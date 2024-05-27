#include "Lattice.hpp"

void LatticeBoltzmann::ExchangeGhosts(void){
    // Post-collision boundary planes go to the neighbours' ghost planes
    int plane = Lx*Ly*3*4*2*2;
    MPI_Sendrecv(&f_new[nf(0,0,Nz,0,0,0,0)], plane, MPI_DOUBLE, right, 0,
                 &f_new[nf(0,0,0,0,0,0,0)], plane, MPI_DOUBLE, left, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&f_new[nf(0,0,1,0,0,0,0)], plane, MPI_DOUBLE, left, 1,
                 &f_new[nf(0,0,Nz+1,0,0,0,0)], plane, MPI_DOUBLE, right, 1,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

void LatticeBoltzmann::Advection(void){
    // Pull scheme: every owned cell takes its populations from the upstream
    // neighbour, which may sit in a ghost plane
    ExchangeGhosts();
    int ixold, iyold, izold;
    for(int ix=0; ix<Lx;ix++){
        for(int iy=0;iy<Ly;iy++){
            for(int iz=1;iz<=Nz;iz++){
                for(int r=0;r<2;r++){
                    f_0[nf0(ix,iy,iz,r)] = f_0_new[nf0(ix,iy,iz,r)];
                }
                for(int p=0;p<3;p++){
                    for(int i=0; i<4;i++){
                        ixold=(ix-V[0][p][i]+Lx)%Lx; iyold=(iy-V[1][p][i]+Ly)%Ly; izold=iz-V[2][p][i];
                        for(int j=0;j<2;j++){
                            for(int r=0;r<2;r++){
                                f[nf(ix,iy,iz,p,i,j,r)] = f_new[nf(ixold,iyold,izold,p,i,j,r)];
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
            for(int iz=1; iz<=Nz ; iz++){
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
            for (int iz = 1; iz <= Nz; iz++){
                double pulse = exp(-alpha*(zg(iz)-z0)*(zg(iz)-z0));
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
    // Rank 0 gathers Ex, By, U, eps of every slab and writes them in z order
    std::vector<int> counts(size), displs(size);
    for (int k = 0; k < size; k++){
        counts[k] = 4*(Lz/size + (k < Lz % size ? 1 : 0));
        displs[k] = k == 0 ? 0 : displs[k-1] + counts[k-1];
    }
    std::vector<double> local(4*Nz), global(rank == 0 ? 4*Lz : 0);
    for (int ix = 0; ix < Lx; ix++)
    {
        for (int iy = 0; iy < Ly; iy++){

        for (int iz = 1; iz <= Nz; iz++){
            vector3D E0 = E_med(ix,iy,iz,false);
            vector3D B0 = B(ix,iy,iz,false);
            local[4*(iz-1)] = E0.x(); local[4*(iz-1)+1] = B0.y();
            local[4*(iz-1)+2] = U(ix,iy,iz,false); local[4*(iz-1)+3] = epsr(ix,iy,iz);
            }
        MPI_Gatherv(local.data(), 4*Nz, MPI_DOUBLE, global.data(), counts.data(), displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (rank == 0)
            for (int iz = 0; iz < Lz; iz++)
                std::cout<<iz<<" "<<global[4*iz]<<" "<<global[4*iz+1]<<" "<<global[4*iz+2]<<" "<<global[4*iz+3]<<std::endl;
        }
        if (rank == 0) std::cout<<std::endl;
    }
}

void LatticeBoltzmann::Peaks(double &Et, double &Er)
{
    double Et_loc = 0, Er_loc = 0;
    for (int iz = 1; iz <= Nz; iz++){
        double Ex = E_med(0,0,iz,false).x();
        if (zg(iz) > Lz/2 && Ex > Et_loc) Et_loc = Ex;
        if (zg(iz) < Lz/2 && Ex < Er_loc) Er_loc = Ex;
    }
    MPI_Allreduce(&Et_loc, &Et, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&Er_loc, &Er, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
}
