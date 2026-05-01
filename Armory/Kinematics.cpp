#include "Kinematics.h"
//#include "elastcaller.h"
//double Kinematics::getQval(double m1, double m2, double m3, double t1, double t3, double angle3)
double Kinematics::getExc(double t3, double angle3)
/*
  \brief Follows convention in Marion, 2013: (1 - beam, 2- target, 3-ejectile, 4-recoil)

  m1 is beam, (typically heavy nucleus)
  m2 is 'd', (light target)
  m3 is 'p', (light ejectile mass)
  t1 is beam kinetic energy,
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
    //t1 = slowitdown("75Ga",t1,"1(12C)2(2H)",1.2);

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
    T4 = ET - e3 - (m1+m2-m3-Q);
    //T4 = slowitdown("75Ga",T4,"1(12C)2(2H)",1.2);
    P4 = TMath::Sqrt(T4*T4 + 2*m4*T4);

    //this angle will not be affected by eloss
    theta4 = (180./M_PI)*TMath::ASin((p3/P4)*TMath::Sin(angle3*M_PI/180.));

//    T4-=16.5; //eloss in about 1.4 mg CD2
    T4-=15.5; //eloss in 1.35 mg CD2 //TODO: actually degrade recoil/ejectiles in target
//    T4-=14.1; //eloss in about 1.2 mg CD2
//    T4-=31.7; //eloss in 2.7 mg CD2

    //recalculate everything other than angle with lowered Kinetic energy
    P4 = TMath::Sqrt(T4*T4 + 2*m4*T4);
    gamma4 = T4/m4+1.;
    beta4 = TMath::Sqrt(1. - 1./(gamma4*gamma4));
    //beta4 = TMath::Sqrt((P4*P4)/(P4*P4 + m4*m4));
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
