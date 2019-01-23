#include "EM.hxx"
#include <cmath>
#include <iostream>
//----------------------------------
using namespace std;
//----------------------------------
class particle{
public:
	particle(){x=0; y=0; z=0;  // location at t = 0
		         px=0; py=0; pz=0;  // momentum at t = -dt
						 gam=0; }; // gamma at t= -dt};
	// Print particle position, momentum u and gamma-factor
	void print_info(const double t) const;
    void boris_push(double Ex,double Ey, double Ez,double Bx, double By,double Bz,const double dt);
	
  double get_x() const{return x;}
private:
	double x,y,z;
	double px,py,pz;
    double gam;
    void calcgam();
	// Calculate gamma-factor from current momentum
    double bla(){return sqrt(1+(px*px+py*py+pz*pz));};
};
//----------------------------------
//----------------------------------
int main(){
  const double FWHM = 10;
  const double a0 = 1.44;
  const double T0 = 200;
  const double Tend = 800;
  const double dt = 0.001;
  const int N = int(Tend/dt + 0.5);

  double Ex,Ey,Ez, Bx, By, Bz;
  double t=0;

  EM em(a0, FWHM,T0);
  particle p;

  for(int i=0; i<N; i++){
	  em.fields(Ex,Ey,Ez,Bx,By,Bz, t + dt*0.5, p.get_x());
	  p.boris_push(Ex, Ey, 0, Bx, 0, Bz, dt);
      //p.gam();
	  if (i%5 == 0 ) p.print_info(t);
	  t += dt;

  }


  return 0;
}
//----------------------------------
 void particle::boris_push(double Ex,double Ey,double Ez,double Bx,double By,double Bz,const double dt){
     double px1 = px + dt/2.0 * Ex;
     double py1 = py + dt /2.0 * Ey;
     double pz1 = pz + dt /2.0 * Ez;
     
     calcgam();
     
     double tx = Bx* dt*0.5/gam;
     double ty = By*dt*0.5/gam;
     double tz = Bz*dt*0.5 / gam;
     double psx = py1*tz -pz1 * ty;
     double psy = pz1*tx -px1*tz;
     double psz = px1 * ty - py1*tx;
     
     psx += px1;
     psy += py1;
     psz += pz1;
     
     double n = 1+ tx*tx + ty*ty +tz*tz;
     double sx = 2*tx /n;
     double sy = 2*ty /n;
     double sz = 2*tz/n;
     
     double ppx = psy *sz -psz*sy;
     double ppy = psz * sx - psx *sz;
     double ppz = psx * sy - psy *sx;
     
     ppx += px1;
     ppy += py1;
     ppz += pz1;
     
     px = ppx + Ex * 0.5*dt;
     py = ppy + Ey* 0.5 *dt;
     pz = ppz + Ez * 0.5 *dt;
     
    calcgam();
    x += px /gam *dt;
    y += py /gam *dt;
    z += pz / gam *dt;
     
}
//----------------------------------
void particle::calcgam(){
    gam = sqrt(1.0 + px*px + py*py + pz * pz);
}
//--------------------------------------
void particle::print_info(const double t) const
{
	cout << t << "\t" << x << "\t" << y << "\t" << z
			 << "\t" << px << "\t" << py << "\t" << pz << "\t" << gam << endl;
}
