#ifndef KINEMATICS_H
#define KINEMATICS_H
#include <TMath.h>
#include <iostream>
#include <fstream>
#include <string>
#include <TVector3.h>
const double u_MeV = 931.49410372; //u in MeV

class Kinematics {
public:
	/*
		A(d,p)B is used as template, with A being beam, and p being ejectile.
		Always, make m3 the thing you detect, and m1 the beam
	*/
    double m_A, m_d, m_p, m_B;
    double E_beam;

/**
 * @{ \name List of funny impossible default values
 */
/**
 * \brief Default values used for all the physics values
 */
    double Q0=-9999, Qx=-9999;
    double P4=-9999, E4=-9999, T4=-9999; //heavy recoil momentum, totE, KE
    double P3=-9999, E3=-9999, T3=-9999; //light recoil
    double ET=-9999;
    double gamma4=-9999, beta4=-9999, theta4=-9999; //theta=heavy-recoil lab angle
    double brho=-9999;
/**
 * @}
 */
    Kinematics(double m1, double m2, double m3, double m4, double ebeam) {
	/*
		A(d,p)B is used as template, with A being beam, and p being ejectile.
		Always, make m3 the thing you detect, m2 the target, and m1 the beam
		
		ebeam is in MeV/u, all others are in amu
	*/
        m_A = m1;
        m_d = m2;
        m_p = m3, m_B = m4, E_beam = ebeam*m_A;
        Q0 = (- m_B - m_p + m_d + m_A)*u_MeV;
    }
    Kinematics() {}

    void setValues(double m1, double m2, double m3, double m4, double ebeam) {
	/*
		Can be used to 'live update' say the beam energy in the case of active target detectors.
	*/

        m_A = m1;
        m_d = m2;
        m_p = m3, m_B = m4, E_beam = ebeam*m_A;
        Q0 = (- m_B - m_p + m_d + m_A)*u_MeV;
        //std::cout << "Q0 MeV: " << Q0 << std::endl;
    }

    void setEBeam(double Ebeam) {E_beam = Ebeam;}

    double getBeta4(double t3, double angle3);
    double getTheta4(double t3, double angle3);
    double getBrho(double t3, double angle3, double charge_state);
    double getExc(double t3, double angle3); //t3 is proton energy detected in ORRUBA, angle3 is proton angle in degrees

    double getBeta4_fromvec(double t3, const TVector3 &pos, const TVector3 &origin) {
        TVector3 local = pos - origin; //position w.r.t origin
        float angle = local.Theta()*180./M_PI;
        return getBeta4(t3, angle);
    }

    double getExc_fromvec(double t3, const TVector3 &pos, const TVector3 &origin) {
        TVector3 local = pos - origin; //position w.r.t origin
        float angle = local.Theta()*180./M_PI;
        return getExc(t3, angle);
    }

    void setValuesFromFile(const std::string& filename) {
    	(void) filename;
        /*std::ifstream in;
        in.open(filename);
        if(!in) {
            std::cerr<< "File not open at " << filename << std::endl;
            return;
        }
        for(std::string line; std::getline(in, line); ) {
            if(line.size()!=0 && line[0]=='#')
                ; //don't do anything with '#' lines
            else {
                std::stringstream ss(line);
                ss>>m_A>>m_d>>m_p>>m_B>>E_beam;
            }
        }
        in.close();*/

    }
};

//double Kinematics::getQval(double m1, double m2, double m3, double t1, double t3, double angle3)
double Kinematics::getExc(double t3, double angle3)
/*
  \brief Follows convention in Marion, 2013: (1 - beam, 2- target, 3-ejectile, 4-recoil)

  m1 is beam, (typically heavy nucleus)
  m2 is 'd', (light target)
  m3 is 'p', (light ejectile mass)
  t1 is beam kinetic energy,
  
  All calculations are done here, and other wrapper functions written make derived quantities from stuff calculated here.
  
  \param t3 lab-kinetic-energy in MeV
  \param angle3 lab-angle in deg of detected proton (if d,p) or other charged particle
  \return Excitation energy of the heavy-recoil nucleus in MeV
*/
{
    double m1 = m_A;
    double m2 = m_d;
    double m3 = m_p;
    double m4 = m_B;
    double t1 = E_beam; 

    m1 *= u_MeV;
    m2 *= u_MeV;
    m3 *= u_MeV;
    m4 *= u_MeV;

    double e1 = m1 + t1;
    double e3 = m3 + t3;
    ET = t1 + m1 + m2;
    double p1 = TMath::Sqrt(t1*t1 + 2*m1*t1);
    double p3 = TMath::Sqrt(t3*t3 + 2*m3*t3);
    double cosTheta = TMath::Cos(angle3*TMath::Pi()/180.);
//    return m1+m2-m3-TMath::Sqrt(m1*m1 + m2*m2 + m3*m3 + 2.*m2*e1 - 2.*e3*(e1+m2)+ 2.*p1*p3*cosTheta);
    double Q =  m1+m2-m3-TMath::Sqrt(m1*m1 + m2*m2 + m3*m3 + 2.*m2*e1 - 2.*e3*(e1+m2)+ 2.*p1*p3*cosTheta);

    Qx = Q;
    //Recoil properties just in case it's useful
    T4 = ET - e3 - (m1+m2-m3-Q);
    P4 = TMath::Sqrt(T4*T4 + 2*m4*T4);

    //this angle will not be affected by eloss
    theta4 = (180./M_PI)*TMath::ASin((p3/P4)*TMath::Sin(angle3*M_PI/180.));

    //recalculate everything other than angle with lowered Kinetic energy
    P4 = TMath::Sqrt(T4*T4 + 2*m4*T4);
    gamma4 = T4/m4+1.;
    beta4 = TMath::Sqrt(1. - 1./(gamma4*gamma4));
    theta4 = (180./M_PI)*TMath::ASin((p3/P4)*TMath::Sin(angle3*M_PI/180.));

    return Q0 - Q;//Q0  = Q + Exc

}

double Kinematics::getBeta4(double t3, double angle3) 
/*
  \brief Follows convention in Marion, 2013: (1 - beam, 2- target, 3-ejectile, 4-recoil)

  m1 is beam, (typically heavy nucleus)
  m2 is 'd', (light target)
  m3 is 'p', (light ejectile mass)
  t1 is beam kinetic energy,
  \param t3 lab-kinetic-energy in MeV
  \param angle3 lab-angle in deg of detected proton (if d,p) or other charged particle
  \return doppler-shift beta (=v/c) of the heavy-recoil nucleus, calls getExc() to fill the value
*/
{
    getExc(t3, angle3);
    return beta4;
}

double Kinematics::getBrho(double t3, double angle3, double charge_state) {
/*
  \brief Follows convention in Marion, 2013: (1 - beam, 2- target, 3-ejectile, 4-recoil)

  m1 is beam, (typically heavy nucleus)
  m2 is 'd', (light target)
  m3 is 'p', (light ejectile mass)
  t1 is beam kinetic energy,
  \param t3 lab-kinetic-energy in MeV
  \param angle3 lab-angle in deg of detected proton (if d,p) or other charged particle
  \param charge_state charge state of the intended nucleus, in units of elementary charge (= +1 for H+, +2 for He2+ etc).
  \return b-rho value, generated from P4*3.3359e-3/charge_state where P4 is 4momentum of heavy-recoil calculated from orruba kinematics
*/
	getExc(t3,angle3);
	return P4*3.3359e-3/charge_state;
}

double Kinematics::getTheta4(double t3, double angle3) {
/*
  \brief Follows convention in Marion, 2013: (1 - beam, 2- target, 3-ejectile, 4-recoil)

  m1 is beam, (typically heavy nucleus)
  m2 is 'd', (light target)
  m3 is 'p', (light ejectile mass)
  t1 is beam kinetic energy,
  \param t3 lab-kinetic-energy in MeV
  \param angle3 lab-angle in deg of detected proton (if d,p) or other charged particle
  \param charge_state charge state of the intended nucleus, in units of elementary charge (= +1 for H+, +2 for He2+ etc).
  \return lab theta value of heavy recoil in degrees
*/
	getExc(t3,angle3);
	return theta4;
}


#endif
