#include "/home/sud/Desktop/Software2/propagator/elastcaller.h"
void make_eloss_table_protons() {
    double einput = 20.0, estepnow;
    double target_thickness_unit = 4e-2; //mg/cm2. 
    double density = 0.0711;//mg/cm3
    long i=0;
    while(einput > 0.001) {
        std::cout << "After " << i << " steps, 1H is at " << einput << " MeV after penetrating " << i*target_thickness_unit << " mg/cm2 " << i*target_thickness_unit/density << " cm of HeCO2" << std::endl;
        estepnow = slowmedown("1H",einput,"3(12C)6(16O)97(4He)",target_thickness_unit);

        einput = estepnow;
        i+=1;
    }
}

void make_eloss_table() {
    double einput = 20.0, estepnow;
    double target_thickness_unit = 1e-3; //mg/cm2. 
    double density = 0.0711;//mg/cm3
    long i=0;
    while(einput > 0.001) {
        std::cout << "After " << i << " steps, 4He is at " << einput << " MeV after penetrating " << i*target_thickness_unit << " mg/cm2 " << i*target_thickness_unit/density << " cm of HeCO2" << std::endl;
        estepnow = slowmedown("4He",einput,"3(12C)6(16O)97(4He)",target_thickness_unit);


        einput = estepnow;
        i+=1;
    }
}
