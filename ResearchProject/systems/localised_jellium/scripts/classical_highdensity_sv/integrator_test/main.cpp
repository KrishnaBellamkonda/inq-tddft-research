// Phase 2 Test B: velocity-Verlet integrator (inqkit::dynamics::Projectile) vs an
// independent Newton-ODE. Pure host (no INQ): fire the projectile head-on at a fixed
// repulsive Gaussian source (analytic two-Gaussian force), integrate, dump trajectory.
#include <inqkit/dynamics/projectile.hpp>
#include <cmath>
#include <fstream>
using inqkit::dynamics::Projectile;
using inqkit::detail::Vec3;
int main(){
  const double sp=0.35355, ss=0.5, s12=std::sqrt(sp*sp+ss*ss), a=1.0/(std::sqrt(2.0)*s12);
  const double m=1.0, dt=0.01, Q=1.0;
  auto Uof=[&](double d){ return Q*std::erf(a*d)/d; };
  auto Fmag=[&](double d){ return Q*( std::erf(a*d)/(d*d) - (2*a/std::sqrt(M_PI))*std::exp(-a*a*d*d)/d ); };
  double z0=8.0, v0=-0.8;                      // launched toward the source at origin
  Projectile proj(m, -1.0, Vec3{0,0,z0}, Vec3{0,0,v0});
  std::ofstream csv("results/traj.csv"); csv.precision(12);
  csv<<"step,t,z,v,KE,U,E\n";
  const int N=3000;
  for(int n=0;n<=N;n++){
    double zn=proj.R().z;
    double Fz=(zn>=0? +1.0:-1.0)*Fmag(std::fabs(zn));   // repulsive: pushes away from origin
    proj.advance(Vec3{0,0,Fz}, dt);
    double vn=proj.V().z, KE=proj.ke(), U=Uof(std::fabs(zn));
    csv<<n<<","<<n*dt<<","<<zn<<","<<vn<<","<<KE<<","<<U<<","<<(KE+U)<<"\n";
  }
  return 0;
}
