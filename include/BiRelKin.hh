////////////////////////////////////
////        BiRelKin        ////
////////////////////////////////////
//
//  A calculator for Binary Relativistic Kinematics
//  Kevin C.W. Li, kcwli@sun.ac.za
//  10/02/15

//  Based upon the Derivations for two-body kinematics (both relativistic and non-relativistic), C. Wheldon

#include "iostream"

////////////////////////
////    Constants   ////
////////////////////////

////  Speed of Light
double c2 = 931.494;     // MeV/u, c^2
double c4 = c2*c2;  // (MeV/u)^2, c^4

void BiRelKin(double *m, double *T, double *E, double *p, double ThetaSCAT_ejectile,  double &ThetaSCAT_recoil, double Ex) {
    
    ////    Q-value calculation
    double Q = (m[0] + m[1])*c2 - (m[2] + m[3])*c2; // MeV
    
    ////    Conversion of ejectile scattering angle from degrees to radians
    double theta_lab_ejectile = ThetaSCAT_ejectile*0.017453292; // radians
    //double theta_lab_recoil = ThetaSCAT_recoil*0.017453292; // radians
    double theta_lab_recoil = 0.0; // radians
    
    ////    Initial Total Energy Calculation
    E[0] = T[0] + (m[0]*c2);
    E[1] = T[1] + (m[1]*c2);
    double Etot = E[0] + E[1] + Q - Ex;
    
    ////    Initial Momentum Calculation
    p[0] = sqrt((E[0]*E[0]) - (m[0]*m[0]*c4))*(1.0/pow(c2, 0.5));
    p[1] = sqrt((E[1]*E[1]) - (m[1]*m[1]*c4))*(1.0/pow(c2, 0.5));
    
    ////////////////////////////////////////////////////////////////////////////////
    ////    E[2] is now solved via the solution of a quadratic equation
    double a = 4*p[0]*p[0]*c2*cos(theta_lab_ejectile)*cos(theta_lab_ejectile) - 4*Etot*Etot;
    double b = (4*Etot*Etot*Etot) - (4*p[0]*p[0]*c2*Etot) + (4*m[2]*m[2]*c4*Etot) - (4*m[3]*m[3]*c4*Etot);
    double c = (2*p[0]*p[0]*c2*Etot*Etot) - (2*m[2]*m[2]*c4*Etot*Etot) + (2*m[2]*m[2]*c4*p[0]*p[0]*c2) + (2*m[3]*m[3]*c4*Etot*Etot) - (2*m[3]*m[3]*c4*p[0]*p[0]*c2) + (2*m[3]*m[3]*c4*m[2]*m[2]*c4) - (Etot*Etot*Etot*Etot) - (p[0]*p[0]*p[0]*p[0]*c4) - (m[2]*m[2]*m[2]*m[2]*c4*c4) - (m[3]*m[3]*m[3]*m[3]*c4*c4) - (4*m[2]*m[2]*c4*p[0]*p[0]*c2*cos(theta_lab_ejectile)*cos(theta_lab_ejectile));
    
    ////    Total Energy, Kinetic Energy and Momentum of Ejectile
    if(ThetaSCAT_ejectile<90.0)
    {
        E[2] = (- b - sqrt((b*b) - (4*a*c)))/(2*a);
    }
    else
    {
        E[2] = (- b + sqrt((b*b) - (4*a*c)))/(2*a);
    }
    
    T[2] = E[2] - m[2]*c2;
    p[2] = sqrt((E[2]*E[2]) - (m[2]*m[2]*c4))*(1.0/pow(c2, 0.5));
    
    ////    Total Energy, Kinetic Energy and Momentum of Recoil
    E[3] = Etot - E[2];
    T[3] = E[3] - m[3]*c2;
    p[3] = sqrt((E[3]*E[3]) - (m[3]*m[3]*c4))*(1.0/pow(c2, 0.5));
    
    ////    Angle between beam axis and Recoil velocity
    theta_lab_recoil = asin((p[2]/p[3])*sin(theta_lab_ejectile)); // radians
    ThetaSCAT_recoil = theta_lab_recoil/0.0174532925; // deg
    
    ////    For backward-angled calculations
    if(p[2]>p[0] && ThetaSCAT_ejectile<90.0)
    {
        ThetaSCAT_recoil = 180 - ThetaSCAT_recoil;
    }
    
    bool printResults = false;
    
    if(printResults)
    {
        
        std::cout << "///////////////////////////////" << std::endl;
        std::cout << "////   Particle Massses    ////" << std::endl;
        
        std::cout << "m[0]:    " << m[0] << " u"  << std::endl;
        std::cout << "m[1]:    " << m[1] << " u"  << std::endl;
        std::cout << "m[2]:    " << m[2] << " u"  << std::endl;
        std::cout << "m[3]:    " << m[3] << " u"  << std::endl;
        
        std::cout << "/////////////////////////" << std::endl;
        std::cout << "////   Projectile    ////" << std::endl;
        
        std::cout << "Total Energy, E[0]:    " << E[0] << " MeV"  << std::endl;
        std::cout << "Kinetic Energy, T[0]:    " << T[0] << " MeV"  << std::endl;
        std::cout << "Momentum, p[0]:    " << p[0] << " MeV/c"  << std::endl;
        
        
        std::cout << "///////////////////////" << std::endl;
        std::cout << "////   Ejectile    ////" << std::endl;
        
        std::cout << "Total Energy, E[2]:    " << E[2] << " MeV"  << std::endl;
        std::cout << "Kinetic Energy, T[2]:    " << T[2] << " MeV"  << std::endl;
        std::cout << "Momentum, p[2]:    " << p[2] << " MeV/c"  << std::endl;
        
        
        std::cout << "///////////////////////" << std::endl;
        std::cout << "////   Recoil    ////" << std::endl;
        std::cout << "Total Energy, E[3]:    " << E[3] << " MeV"  << std::endl;
        std::cout << "Kinetic Energy, T[3]:    " << T[3] << " MeV"  << std::endl;
        std::cout << "Momentum, p[3]:    " << p[3] << " MeV/c"  << std::endl;
        
        
        std::cout << "   " << std::endl;
        std::cout << "ThetaSCAT_ejectile (angle of ejectile w.r.t. beam axis):    " << ThetaSCAT_ejectile << " deg"  << std::endl;
        std::cout << "ThetaSCAT_recoil (angle of recoil w.r.t. beam axis):    " << ThetaSCAT_recoil << " deg"  << std::endl;
        
        
        //std::cout << "///////////////////////" << std::endl;
        //std::cout << "////   Ejectile    ////" << std::endl;
        /*
         std::cout << "Total Energy, E[2]:    " << E[2] << " MeV"  << std::endl;
         std::cout << "Momentum, p[2]:    " << p[2] << " MeV/c"  << std::endl;
         std::cout << "Excitation Energy, Ex:    " << Ex << " MeV"  << std::endl;
         std::cout << "Kinetic Energy, T[2]:    " << T[2] << " MeV"  << std::endl;
         */
        
        //std::cout << "Kinetic Energy, T[3]:    " << T[3] << " MeV"  << std::endl;
        //std::cout << "Q:    " << Q << " MeV"  << std::endl;
        
        //std::cout << "Momentum, E[2]:    " << E[2] << " MeV"  << std::endl;
        //std::cout << "Momentum, E[3]:    " << E[3] << " MeV"  << std::endl;
        //std::cout << "Momentum, Etot:    " << Etot << " MeV"  << std::endl;
        
        //Etot = E[0] + E[1] + Q - Ex;
        //std::cout << "Calculated Ex:    " << Ex << " MeV"  << std::endl;
    }
}

double CalculateEjectileEnergy(double mass_projectile_amu, double mass_target_amu, double mass_ejectile_amu, double mass_recoil_amu, double beamEnergy_MeV, double recoilExcitationEnergy_MeV, double ejectilePolarAngle_labFrame_deg)
{
    double m[4];
    double T[4];
    double E[4];
    double p[4];
    
    ////    Excitation of the Recoil Nucleus
    double Ex; // MeV
    
    ////    theta_ejectile is the angle of the Ejectile with respect to the beam axis
    double theta_ejectile; // deg
    
    ////    theta_recoil is the angle of the Recoil with respect to the beam axis
    double theta_recoil; // deg
    
    for(int i=0; i<4; i++)
    {
        m[i] = 0.0;
        T[i] = 0.0;
        E[i] = 0.0;
        p[i] = 0.0;
    }
    
    ////    These are the input variables
    m[0] = mass_projectile_amu;
    m[1] = mass_target_amu;
    m[2] = mass_ejectile_amu;
    m[3] = mass_recoil_amu;
    
    ////    Projectile/Beam Energy (Lab Frame)
    T[0] = beamEnergy_MeV;
    
    ////    Target Energy (Lab Frame)
    T[1] = 0.0;
    
    theta_ejectile = ejectilePolarAngle_labFrame_deg; // deg
    
    Ex = recoilExcitationEnergy_MeV; // MeV
    
    BiRelKin(m, T, E, p, theta_ejectile, theta_recoil, Ex);
    
    return T[2];
}

double CalculateRecoilEnergy(double mass_projectile_amu, double mass_target_amu, double mass_ejectile_amu, double mass_recoil_amu, double beamEnergy_MeV, double recoilExcitationEnergy_MeV, double ejectilePolarAngle_labFrame_deg)
{
    double m[4];
    double T[4];
    double E[4];
    double p[4];
    
    ////    Excitation of the Recoil Nucleus
    double Ex; // MeV
    
    ////    theta_ejectile is the angle of the Ejectile with respect to the beam axis
    double theta_ejectile; // deg
    
    ////    theta_recoil is the angle of the Recoil with respect to the beam axis
    double theta_recoil; // deg
    
    for(int i=0; i<4; i++)
    {
        m[i] = 0.0;
        T[i] = 0.0;
        E[i] = 0.0;
        p[i] = 0.0;
    }
    
    ////    These are the input variables
    m[0] = mass_projectile_amu;
    m[1] = mass_target_amu;
    m[2] = mass_ejectile_amu;
    m[3] = mass_recoil_amu;
    
    ////    Projectile/Beam Energy (Lab Frame)
    T[0] = beamEnergy_MeV;
    
    ////    Target Energy (Lab Frame)
    T[1] = 0.0;
    
    theta_ejectile = ejectilePolarAngle_labFrame_deg; // deg
    
    Ex = recoilExcitationEnergy_MeV; // MeV
    
    BiRelKin(m, T, E, p, theta_ejectile, theta_recoil, Ex);
    
    return T[3];
}

void CalcEx(double *m, double *T, double *E, double *p, double Xpos, double theta, double &Ex, std::vector<double> momentumCalPars) {
    
    ////    Q-value calculation
    double Q = (m[0] + m[1])*c2 - (m[2] + m[3])*c2; // MeV
    
    ////    Conversion of ejectile scattering angle from degrees to radians
    double theta_lab = theta*0.017453292; // radians
    
    ////    Initial Total Energy Calculation
    E[0] = T[0] + (m[0]*c2);
    E[1] = T[1] + (m[1]*c2);
    
    ////    Initial Momentum Calculation
    p[0] = sqrt((E[0]*E[0]) - (m[0]*m[0]*c4))*(1.0/pow(c2, 0.5));
    p[1] = sqrt((E[1]*E[1]) - (m[1]*m[1]*c4))*(1.0/pow(c2, 0.5));
    
    p[2] = 0.0;
    for(int i=0; i<(int) momentumCalPars.size(); i++)
    {
        p[2] += momentumCalPars[i]*(1.0/pow(c2, 0.5))*pow(Xpos, i);
    }
    E[2] = sqrt(p[2]*p[2]*c2 + (m[2]*m[2]*c4));
    T[2] = E[2] - m[2]*c2;
    
    p[3] = sqrt((p[0]*p[0]) + (p[2]*p[2]) - 2*p[0]*p[2]*cos(theta_lab));
    E[3] = sqrt(p[3]*p[3]*c2 + (m[3]*m[3]*c4));
    T[3] = E[3] - m[3]*c2;
    
    Ex = E[0] + E[1] + Q - E[2] - E[3];
}




