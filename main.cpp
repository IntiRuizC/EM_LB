#include "Lattice.hpp"


int main(void){
    LatticeBoltzmann pulse;

    int t;

    pulse.Start();

    for ( t = 0; t < tmax; t++)
    {
        pulse.Collision();

        pulse.Advection();
    }

    pulse.Print();

    double Et, Er;
    pulse.Peaks(Et, Er);
    double n = sqrt(eps2/eps1);
    std::cerr<<"transmitted: "<<Et/E0<<" theory: "<<2/(n+1)<<std::endl;
    std::cerr<<"reflected:   "<<-Er/E0<<" theory: "<<(n-1)/(n+1)<<std::endl;
    return 0;
}
