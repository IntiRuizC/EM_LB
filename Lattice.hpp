#include "constants.hpp"
#include "vector.hpp"
#include <vector>
#include <mpi.h>

#ifndef LATTICE_H
#define LATTICE_H

// The lattice is split in slabs along z, one per MPI rank. Each rank stores its
// Nz owned planes plus one ghost plane on each side (local iz = 0 and Nz+1);
// owned planes are iz = 1..Nz, global z = z_start + iz - 1.
class LatticeBoltzmann{
    private:
        int V[3][3][4]; /*V[xyz][p][i]*/  vector3D v[3][4], v0; //v[p][i]
        vector3D e[3][4][2]; vector3D e0;//e[p][i][j]
        vector3D b[3][4][2]; vector3D b0;//b[p][i][j]
        int rank, size, left, right; //MPI neighbours along z (periodic)
        int Nz, z_start; //owned planes and global z of the first one
        std::vector<double> f, f_new; // f[z][y][x][p][i][j][r], z slowest so a plane is contiguous
        std::vector<double> f_0, f_0_new; //f[z][y][x][r]

    public:
        LatticeBoltzmann(void);
        int nf(int ix, int iy, int iz, int p, int i, int j, int r) {
        return r + 2 * j + 2 * 2 * i + 2 * 2 * 4 * p + 2 * 2 * 4 * 3 * ix + 2 * 2 * 4 * 3 * Lx * iy + 2 * 2 * 4 * 3 * Lx * Ly * iz;}
        int nf0(int ix, int iy, int iz, int r)
        {return r + 2 * ix + 2 * Lx * iy + 2 * Lx * Ly * iz;};
        int zg(int iz) {return z_start + iz - 1;}
        //Media
        double epsr(int ix, int iy, int iz);
        double mur(int ix, int iy, int iz);
        double sigma(int ix, int iy, int iz);
        //Fields and magnituds
        double rho_c(int ix, int iy, int iz, bool use_new);
        vector3D D(int ix, int iy, int iz, bool use_new);
        vector3D E(int ix, int iy, int iz, bool use_new);
        vector3D B(int ix, int iy, int iz, bool use_new);
        vector3D J_med(int ix, int iy, int iz, bool use_new);
        vector3D E_med(int ix, int iy, int iz, bool use_new);
        double U(int ix, int iy, int iz, bool use_new);
        double feq(vector3D &Jmed, vector3D &Emed, vector3D &B, double eps, double mu, int p, int i, int j, int r);
        double feq0(double rhoc);
        //Evolution dynamics
        void Start(void);
        void ExchangeGhosts(void);
        void Advection(void);
        void Collision(void);

        //Extract data
        int Rank(void) {return rank;}
        void Print(void);
        void Peaks(double &Et, double &Er);

};

#endif //LATTICE_H
