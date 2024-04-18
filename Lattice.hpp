#include "constants.hpp"
#include "vector.hpp"

#ifndef LATTICE_H
#define LATTICE_H

class LatticeBoltzmann{
    private:
        int V[3][3][4]; /*V[xyz][p][i]*/  vector3D v[3][4], v0; //v[p][i]
        vector3D e[3][4][2]; vector3D e0;//e[p][i][j]
        vector3D b[3][4][2]; vector3D b0;//b[p][i][j]
        double f[Lx*Ly*Lz*3*4*2*2]; // f[x][y][z][p][i][j][r]
        double f_new[Lx*Ly*Lz*3*4*2*2];
        double f_0[Lx*Ly*Lz*2], f_0_new[Lx*Ly*Lz*2]; //f[x][y][z][r]

    public:
        LatticeBoltzmann(void);
        int nf(int ix, int iy, int iz, int p, int i, int j, int r) {
        return r + 2 * j + 2 * 2 * i + 2 * 2 * 4 * p + 2 * 2 * 4 * 3 * iz + 2 * 2 * 4 * 3 * Lz * iy + 2 * 2 * 4 * 3 * Lz * Ly * ix;}
        int nf0(int ix, int iy, int iz, int r)
        {return r + 2 * iz + 2 * Lz * iy + 2 * Lz * Ly * ix;};
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
        void Advection(void);
        void Collision(void);

        //Extract data
        void Print(void);
        void Peaks(double &Et, double &Er);

};

#endif //LATTICE_H
