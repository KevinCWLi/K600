////////////////////////////////////
////        BiRelKin        ////
////////////////////////////////////
//
//  A calculator for Binary Relativistic Kinematics
//  Kevin C.W. Li, kcwli@sun.ac.za
//  10/02/15

//  Based upon the Derivations for two-body kinematics (both relativistic and non-relativistic), C. Wheldon

#include <cmath>
#include <iomanip>
#include <sstream>
#include "iostream"

////////////////////////
////    Constants   ////
////////////////////////

////  Speed of Light
long double c2 = 931.49410242;     // MeV/u, c^2
long double c4 = c2*c2;  // (MeV/u)^2, c^4

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
    //double Etot = E[0] + E[1] + Q - Ex;
    double Etot = E[0] + E[1] - Ex;
    
    ////    Initial Momentum Calculation
    p[0] = std::sqrt((E[0]*E[0]) - (m[0]*m[0]*c4))*(1.0/std::pow(c2, 0.5));
    p[1] = std::sqrt((E[1]*E[1]) - (m[1]*m[1]*c4))*(1.0/std::pow(c2, 0.5));
    
    ////////////////////////////////////////////////////////////////////////////////
    ////    E[2] is now solved via the solution of a quadratic equation
    double a = 4*p[0]*p[0]*c2*std::cos(theta_lab_ejectile)*std::cos(theta_lab_ejectile) - 4*Etot*Etot;
    double b = (4*Etot*Etot*Etot) - (4*p[0]*p[0]*c2*Etot) + (4*m[2]*m[2]*c4*Etot) - (4*m[3]*m[3]*c4*Etot);
    double c = (2*p[0]*p[0]*c2*Etot*Etot) - (2*m[2]*m[2]*c4*Etot*Etot) + (2*m[2]*m[2]*c4*p[0]*p[0]*c2) + (2*m[3]*m[3]*c4*Etot*Etot) - (2*m[3]*m[3]*c4*p[0]*p[0]*c2) + (2*m[3]*m[3]*c4*m[2]*m[2]*c4) - (Etot*Etot*Etot*Etot) - (p[0]*p[0]*p[0]*p[0]*c4) - (m[2]*m[2]*m[2]*m[2]*c4*c4) - (m[3]*m[3]*m[3]*m[3]*c4*c4) - (4*m[2]*m[2]*c4*p[0]*p[0]*c2*std::cos(theta_lab_ejectile)*std::cos(theta_lab_ejectile));
    
    ////    Total Energy, Kinetic Energy and Momentum of Ejectile
    if(ThetaSCAT_ejectile<90.0)
    {
        E[2] = (- b - std::sqrt((b*b) - (4*a*c)))/(2*a);
    }
    else
    {
        E[2] = (- b + std::sqrt((b*b) - (4*a*c)))/(2*a);
        //
        //        if(std::abs(ThetaSCAT_ejectile-90)<0.2)
        //        {
        //            std::cout << std::setprecision(20) << "(b*b) - (4*a*c): " << (b*b) - (4*a*c) << std::endl;
        //            std::cout << std::setprecision(20) << "(b*b): " << (b*b) << std::endl;
        //            std::cout << std::setprecision(20) << "(4*a*c): " << (4*a*c) << std::endl;
        //            std::cout << "E[2]: " << E[2] << std::endl;
        //        }
    }
    
    if((b*b) - (4*a*c)<0.0)
    {
        E[2] = (- b + std::sqrt(0))/(2*a);
    }
    
    
    //    std::cout << "E[2] = (- b - std::sqrt((b*b) - (4*a*c)))/(2*a): " << (- b - std::sqrt((b*b) - (4*a*c)))/(2*a) << std::endl;
    //    std::cout << "E[2] = (- b + std::sqrt((b*b) - (4*a*c)))/(2*a): " << (- b - std::sqrt((b*b) + (4*a*c)))/(2*a) << std::endl;
    
    //    std::cout << "Q: " << Q << std::endl;
    //    std::cout << "E[2]: " << E[2] << std::endl;
    //    E[2] = 0.0;
    //    std::cout << "E[2]: " << E[2] << std::endl;
    
    T[2] = E[2] - m[2]*c2;
    p[2] = std::sqrt((E[2]*E[2]) - (m[2]*m[2]*c4))*(1.0/std::pow(c2, 0.5));
    
    //    if(std::abs(ThetaSCAT_ejectile-90)<0.2)
    //    {
    //        std::cout << "a: " << a << std::endl;
    //        std::cout << "b: " << b << std::endl;
    //        std::cout << "c: " << c << std::endl;
    //
    //        std::cout << "(b*b) - (4*a*c): " << (b*b) - (4*a*c) << std::endl;
    //        std::cout << "ThetaSCAT_ejectile: " << ThetaSCAT_ejectile << std::endl;
    //        std::cout << "E[2] = (- b - std::sqrt((b*b) - (4*a*c)))/(2*a): " << (- b - std::sqrt((b*b) - (4*a*c)))/(2*a) << std::endl;
    //        std::cout << "E[2] = (- b + std::sqrt((b*b) - (4*a*c)))/(2*a): " << (- b + std::sqrt((b*b) - (4*a*c)))/(2*a) << std::endl;
    //        std::cout << "T[2]: " << T[2] << std::endl;
    //
    //        std::cout << std::endl;
    //    }
    
    ////    Total Energy, Kinetic Energy and Momentum of Recoil
    E[3] = Etot - E[2];
    T[3] = E[3] - m[3]*c2;
    p[3] = std::sqrt((E[3]*E[3]) - (m[3]*m[3]*c4))*(1.0/std::pow(c2, 0.5));
    
    ////    Angle between beam axis and Recoil velocity
    theta_lab_recoil = std::asin((p[2]/p[3])*std::sin(theta_lab_ejectile)); // radians
    ThetaSCAT_recoil = theta_lab_recoil/0.0174532925; // deg
    
    //    std::cout << "p[2]: " << p[2] << std::endl;
    //    std::cout << "p[3]: " << p[3] << std::endl;
    
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

void BiRelKin(long double *m, long double *T, long double *E, long double *p, long double ThetaSCAT_ejectile,  long double &ThetaSCAT_recoil, long double Ex) {
    
    ////    Q-value calculation
    long double Q = (m[0] + m[1])*c2 - (m[2] + m[3])*c2; // MeV
    
    ////    Conversion of ejectile scattering angle from degrees to radians
    long double theta_lab_ejectile = ThetaSCAT_ejectile*0.017453292; // radians
    //long double theta_lab_recoil = ThetaSCAT_recoil*0.017453292; // radians
    long double theta_lab_recoil = 0.0; // radians
    
    ////    Initial Total Energy Calculation
    E[0] = T[0] + (m[0]*c2);
    E[1] = T[1] + (m[1]*c2);
    //long double Etot = E[0] + E[1] + Q - Ex;
    long double Etot = E[0] + E[1] - Ex;
    
    ////    Initial Momentum Calculation
    p[0] = std::sqrt((E[0]*E[0]) - (m[0]*m[0]*c4))*(1.0/std::pow(c2, 0.5));
    p[1] = std::sqrt((E[1]*E[1]) - (m[1]*m[1]*c4))*(1.0/std::pow(c2, 0.5));
    
    ////////////////////////////////////////////////////////////////////////////////
    ////    E[2] is now solved via the solution of a quadratic equation
    long double a = 4*p[0]*p[0]*c2*std::cos(theta_lab_ejectile)*std::cos(theta_lab_ejectile) - 4*Etot*Etot;
    long double b = (4*Etot*Etot*Etot) - (4*p[0]*p[0]*c2*Etot) + (4*m[2]*m[2]*c4*Etot) - (4*m[3]*m[3]*c4*Etot);
    long double c = (2*p[0]*p[0]*c2*Etot*Etot) - (2*m[2]*m[2]*c4*Etot*Etot) + (2*m[2]*m[2]*c4*p[0]*p[0]*c2) + (2*m[3]*m[3]*c4*Etot*Etot) - (2*m[3]*m[3]*c4*p[0]*p[0]*c2) + (2*m[3]*m[3]*c4*m[2]*m[2]*c4) - (Etot*Etot*Etot*Etot) - (p[0]*p[0]*p[0]*p[0]*c4) - (m[2]*m[2]*m[2]*m[2]*c4*c4) - (m[3]*m[3]*m[3]*m[3]*c4*c4) - (4*m[2]*m[2]*c4*p[0]*p[0]*c2*std::cos(theta_lab_ejectile)*std::cos(theta_lab_ejectile));
    
    ////    Total Energy, Kinetic Energy and Momentum of Ejectile
    if(ThetaSCAT_ejectile<90.0)
    {
        E[2] = (- b - std::sqrt((b*b) - (4*a*c)))/(2*a);
    }
    else
    {
        E[2] = (- b + std::sqrt((b*b) - (4*a*c)))/(2*a);
        //
        //        if(std::abs(ThetaSCAT_ejectile-90)<0.2)
        //        {
        //            std::cout << std::setprecision(20) << "(b*b) - (4*a*c): " << (b*b) - (4*a*c) << std::endl;
        //            std::cout << std::setprecision(20) << "(b*b): " << (b*b) << std::endl;
        //            std::cout << std::setprecision(20) << "(4*a*c): " << (4*a*c) << std::endl;
        //            std::cout << "E[2]: " << E[2] << std::endl;
        //        }
    }
    
    if((b*b) - (4*a*c)<0.0)
    {
        E[2] = (- b + std::sqrt(0))/(2*a);
    }
    
    
    //    std::cout << "E[2] = (- b - std::sqrt((b*b) - (4*a*c)))/(2*a): " << (- b - std::sqrt((b*b) - (4*a*c)))/(2*a) << std::endl;
    //    std::cout << "E[2] = (- b + std::sqrt((b*b) - (4*a*c)))/(2*a): " << (- b - std::sqrt((b*b) + (4*a*c)))/(2*a) << std::endl;
    
    //    std::cout << "Q: " << Q << std::endl;
    //    std::cout << "E[2]: " << E[2] << std::endl;
    //    E[2] = 0.0;
    //    std::cout << "E[2]: " << E[2] << std::endl;
    
    T[2] = E[2] - m[2]*c2;
    p[2] = std::sqrt((E[2]*E[2]) - (m[2]*m[2]*c4))*(1.0/std::pow(c2, 0.5));
    
    //    if(std::abs(ThetaSCAT_ejectile-90)<0.2)
    //    {
    //        std::cout << "a: " << a << std::endl;
    //        std::cout << "b: " << b << std::endl;
    //        std::cout << "c: " << c << std::endl;
    //
    //        std::cout << "(b*b) - (4*a*c): " << (b*b) - (4*a*c) << std::endl;
    //        std::cout << "ThetaSCAT_ejectile: " << ThetaSCAT_ejectile << std::endl;
    //        std::cout << "E[2] = (- b - std::sqrt((b*b) - (4*a*c)))/(2*a): " << (- b - std::sqrt((b*b) - (4*a*c)))/(2*a) << std::endl;
    //        std::cout << "E[2] = (- b + std::sqrt((b*b) - (4*a*c)))/(2*a): " << (- b + std::sqrt((b*b) - (4*a*c)))/(2*a) << std::endl;
    //        std::cout << "T[2]: " << T[2] << std::endl;
    //
    //        std::cout << std::endl;
    //    }
    
    ////    Total Energy, Kinetic Energy and Momentum of Recoil
    E[3] = Etot - E[2];
    T[3] = E[3] - m[3]*c2;
    p[3] = std::sqrt((E[3]*E[3]) - (m[3]*m[3]*c4))*(1.0/std::pow(c2, 0.5));
    
    ////    Angle between beam axis and Recoil velocity
    theta_lab_recoil = std::asin((p[2]/p[3])*std::sin(theta_lab_ejectile)); // radians
    ThetaSCAT_recoil = theta_lab_recoil/0.0174532925; // deg
    
    //    std::cout << "p[2]: " << p[2] << std::endl;
    //    std::cout << "p[3]: " << p[3] << std::endl;
    
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

void DetermineCOMDynamics(double *Masses, double TBeam, double Ex)
{
    //Function to calculate the recoil TLorentzVector in the lab frame
    
    double RecoilMass = Masses[3] + Ex; //Convert the recoil mass to the right value including the excitation-energy dependence
    
    double s = std::pow(Masses[0],2.) + std::pow(Masses[1],2.) + 2*Masses[1]*(TBeam + Masses[0]);
    
    std::cout << "s: " << s << std::endl;
    
    double TBeam_inverseKinematics = (s - std::pow(Masses[0],2.) - std::pow(Masses[1],2.))/(2.0*Masses[0]) - Masses[1];
    
    std::cout << "TBeam_inverseKinematics: " << TBeam_inverseKinematics << std::endl;
    
    //TLorentzVector f4MomentumCentreOfMass(0,0,0,std::sqrt(s));
    
    
    /*
     double ECM0 = (s + std::pow(Masses[0],2.) - std::pow(Masses[1],2.))/(2*std::sqrt(s));
     double ECM1 = (s + std::pow(Masses[1],2.) - std::pow(Masses[0],2.))/(2*std::sqrt(s));
     double ECM2 = (s + std::pow(Masses[2],2.) - std::pow(RecoilMass,2.))/(2*std::sqrt(s));
     double ECM3 = (s + std::pow(RecoilMass,2.) - std::pow(Masses[2],2.))/(2*std::sqrt(s));
     
     double PCM0 = std::sqrt(std::pow(ECM0,2.) - std::pow(Masses[0],2.));
     double PCM1 = std::sqrt(std::pow(ECM1,2.) - std::pow(Masses[1],2.));
     double PCM2 = std::sqrt(std::pow(ECM2,2.) - std::pow(Masses[2],2.));
     double PCM3 = std::sqrt(std::pow(ECM3,2.) - std::pow(RecoilMass,2.));
     
     TVector3 Lab3Momentum0(0,0,std::sqrt(std::pow(TBeam,2.) + 2 * TBeam * Masses[0]));
     TVector3 Lab3Momentum1(0,0,0);
     
     TLorentzVector Lab4Momentum0(Lab3Momentum0,Masses[0] + TBeam);
     TLorentzVector Lab4Momentum1(Lab3Momentum1,Masses[1]);
     
     double BetaCM = (Lab4Momentum0+Lab4Momentum1).Beta();
     //std::cout << "BetaCM = " << BetaCM << std::endl;
     
     TLorentzVector CoM4Momentum0 = Lab4Momentum0;
     CoM4Momentum0.Boost(TVector3(0,0,-BetaCM));
     
     TLorentzVector CoM4Momentum1 = Lab4Momentum1;
     CoM4Momentum1.Boost(TVector3(0,0,-BetaCM));
     
     TLorentzVector CoM4Momentum2 = TLorentzVector(PCM2 * std::sin(ThetaAlphaCM * M_PI/180.) * std::cos(PhiAlphaCM * M_PI/180.),
     PCM2 * std::sin(ThetaAlphaCM * M_PI/180.) * std::sin(PhiAlphaCM * M_PI/180.),
     PCM2 * std::cos(ThetaAlphaCM * M_PI/180.),
     ECM2);
     TLorentzVector CoM4Momentum3 = f4MomentumCentreOfMass - CoM4Momentum2;
     
     TLorentzVector Lab4Momentum2 = CoM4Momentum2;
     Lab4Momentum2.Boost(0,0,BetaCM);
     TLorentzVector Lab4Momentum3 = CoM4Momentum3;
     Lab4Momentum3.Boost(0,0,BetaCM);
     
     //if(VerboseFlag)Lab4Momentum2.Print();
     //if(VerboseFlag)Lab4Momentum3.Print();
     
     if((Lab4Momentum0+Lab4Momentum1-Lab4Momentum2-Lab4Momentum3).Mag()>1.e-6)
     {
     std::cout << "Problem with two-body kinematics" << std::endl;
     (Lab4Momentum0+Lab4Momentum1-Lab4Momentum2-Lab4Momentum3).Print();
     }
     */
    
    //    result.push_back(Lab4Momentum2);
    //    result.push_back(Lab4Momentum3);
    //
    //    return result;
}

//================================================================================

//double SolidAngleScalingFactor_COMtoLAB(double TBeam, double Ex, double polarAngle_COM)
//{
//    //------------------------------------------------
//    double Masses[4];
//
//    for(int i=0; i<4; i++)
//    {
//        Masses[i] = m[i]*931.49410242;
//    }
//
//    //------------------------------------------------
//    double RecoilMass = Masses[3] + Ex; //Convert the recoil mass to the right value including the excitation-energy dependence
//
//    double s = std::pow(Masses[0],2.) + std::pow(Masses[1],2.) + 2 * Masses[1] * (TBeam + Masses[0]);
//
//    TLorentzVector f4MomentumCentreOfMass(0,0,0,std::sqrt(s));
//
//    double ECM0 = (s + std::pow(Masses[0],2.) - std::pow(Masses[1],2.))/(2*std::sqrt(s));
//    double ECM1 = (s + std::pow(Masses[1],2.) - std::pow(Masses[0],2.))/(2*std::sqrt(s));
//    double ECM2 = (s + std::pow(Masses[2],2.) - std::pow(RecoilMass,2.))/(2*std::sqrt(s));
//    double ECM3 = (s + std::pow(RecoilMass,2.) - std::pow(Masses[2],2.))/(2*std::sqrt(s));
//
//    double PCM0 = std::sqrt(std::pow(ECM0,2.) - std::pow(Masses[0],2.));
//    double PCM1 = std::sqrt(std::pow(ECM1,2.) - std::pow(Masses[1],2.));
//    double PCM2 = std::sqrt(std::pow(ECM2,2.) - std::pow(Masses[2],2.));
//    double PCM3 = std::sqrt(std::pow(ECM3,2.) - std::pow(RecoilMass,2.));
//
//    TVector3 Lab3Momentum0(0,0,std::sqrt(std::pow(TBeam,2.) + 2 * TBeam * Masses[0]));
//    TVector3 Lab3Momentum1(0,0,0);
//
//    TLorentzVector Lab4Momentum0(Lab3Momentum0,Masses[0] + TBeam);
//    TLorentzVector Lab4Momentum1(Lab3Momentum1,Masses[1]);
//
//    double BetaCM = (Lab4Momentum0+Lab4Momentum1).Beta();
//
//    double gamma = 1.0/std::sqrt(1-std::pow(BetaCM,2.0));
//
//    //------------------------------------------------
//    double solidAngleScalingFactor_COMtoLAB = pow(1 + pow(gamma, 2.0) + (2.0*gamma*cos(polarAngle_COM*0.0174533)), 3.0/2.0);
//    solidAngleScalingFactor_COMtoLAB *= (1.0/(1 + (gamma*cos(polarAngle_COM*0.0174533))));
//
//    //std::cout << "solidAngleScalingFactor_COMtoLAB:    " << solidAngleScalingFactor_COMtoLAB << std::endl;
//
//    return solidAngleScalingFactor_COMtoLAB;
//}

void BiRelKin_PrintResults(long double *m, long double *T, long double *E, long double *p, long double ThetaSCAT_ejectile,  long double &ThetaSCAT_recoil, long double Ex) {
    
    ////    Q-value calculation
    long double Q = (m[0] + m[1])*c2 - (m[2] + m[3])*c2; // MeV
    
    ////    Conversion of ejectile scattering angle from degrees to radians
    long double theta_lab_ejectile = ThetaSCAT_ejectile*0.017453292; // radians
    //long double theta_lab_recoil = ThetaSCAT_recoil*0.017453292; // radians
    long double theta_lab_recoil = 0.0; // radians
    
    ////    Initial Total Energy Calculation
    E[0] = T[0] + (m[0]*c2);
    E[1] = T[1] + (m[1]*c2);
    long double Etot = E[0] + E[1] + Q - Ex;
    
    ////    Initial Momentum Calculation
    p[0] = std::sqrt((E[0]*E[0]) - (m[0]*m[0]*c4))*(1.0/std::pow(c2, 0.5));
    p[1] = std::sqrt((E[1]*E[1]) - (m[1]*m[1]*c4))*(1.0/std::pow(c2, 0.5));
    
    ////////////////////////////////////////////////////////////////////////////////
    ////    E[2] is now solved via the solution of a quadratic equation
    long double a = 4*p[0]*p[0]*c2*std::cos(theta_lab_ejectile)*std::cos(theta_lab_ejectile) - 4*Etot*Etot;
    long double b = (4*Etot*Etot*Etot) - (4*p[0]*p[0]*c2*Etot) + (4*m[2]*m[2]*c4*Etot) - (4*m[3]*m[3]*c4*Etot);
    long double c = (2*p[0]*p[0]*c2*Etot*Etot) - (2*m[2]*m[2]*c4*Etot*Etot) + (2*m[2]*m[2]*c4*p[0]*p[0]*c2) + (2*m[3]*m[3]*c4*Etot*Etot) - (2*m[3]*m[3]*c4*p[0]*p[0]*c2) + (2*m[3]*m[3]*c4*m[2]*m[2]*c4) - (Etot*Etot*Etot*Etot) - (p[0]*p[0]*p[0]*p[0]*c4) - (m[2]*m[2]*m[2]*m[2]*c4*c4) - (m[3]*m[3]*m[3]*m[3]*c4*c4) - (4*m[2]*m[2]*c4*p[0]*p[0]*c2*std::cos(theta_lab_ejectile)*std::cos(theta_lab_ejectile));
    
    ////    Total Energy, Kinetic Energy and Momentum of Ejectile
    if(ThetaSCAT_ejectile<90.0)
    {
        E[2] = (- b - std::sqrt((b*b) - (4*a*c)))/(2*a);
    }
    else
    {
        E[2] = (- b + std::sqrt((b*b) - (4*a*c)))/(2*a);
    }
    
    T[2] = E[2] - m[2]*c2;
    p[2] = std::sqrt((E[2]*E[2]) - (m[2]*m[2]*c4))*(1.0/std::pow(c2, 0.5));
    
    ////    Total Energy, Kinetic Energy and Momentum of Recoil
    E[3] = Etot - E[2];
    T[3] = E[3] - m[3]*c2;
    p[3] = std::sqrt((E[3]*E[3]) - (m[3]*m[3]*c4))*(1.0/std::pow(c2, 0.5));
    
    ////    Angle between beam axis and Recoil velocity
    theta_lab_recoil = std::asin((p[2]/p[3])*std::sin(theta_lab_ejectile)); // radians
    ThetaSCAT_recoil = theta_lab_recoil/0.0174532925; // deg
    
    ////    For backward-angled calculations
    if(p[2]>p[0] && ThetaSCAT_ejectile<90.0)
    {
        ThetaSCAT_recoil = 180 - ThetaSCAT_recoil;
    }
    
    bool printResults = true;
    
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

void BiRelKin_PrintResults(long double mass_projectile_amu, long double mass_target_amu, long double mass_ejectile_amu, long double mass_recoil_amu, long double beamEnergy_MeV, long double recoilExcitationEnergy_MeV, long double ejectilePolarAngle_labFrame_deg) {
    
    long double m[4];
    long double T[4];
    long double E[4];
    long double p[4];
    
    ////    Excitation of the Recoil Nucleus
    long double Ex; // MeV
    
    ////    theta_ejectile is the angle of the Ejectile with respect to the beam axis
    long double theta_ejectile; // deg
    
    ////    theta_recoil is the angle of the Recoil with respect to the beam axis
    long double theta_recoil; // deg
    
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
    
    BiRelKin_PrintResults(m, T, E, p, theta_ejectile, theta_recoil, Ex);
}

long double CalculateEjectileEnergy(long double mass_projectile_amu, long double mass_target_amu, long double mass_ejectile_amu, long double mass_recoil_amu, long double beamEnergy_MeV, long double recoilExcitationEnergy_MeV, long double ejectilePolarAngle_labFrame_deg)
{
    long double m[4];
    long double T[4];
    long double E[4];
    long double p[4];
    
    ////    Excitation of the Recoil Nucleus
    long double Ex; // MeV
    
    ////    theta_ejectile is the angle of the Ejectile with respect to the beam axis
    long double theta_ejectile; // deg
    
    ////    theta_recoil is the angle of the Recoil with respect to the beam axis
    long double theta_recoil; // deg
    
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

long double CalculateRecoilEnergy(long double mass_projectile_amu, long double mass_target_amu, long double mass_ejectile_amu, long double mass_recoil_amu, long double beamEnergy_MeV, long double recoilExcitationEnergy_MeV, long double ejectilePolarAngle_labFrame_deg)
{
    long double m[4];
    long double T[4];
    long double E[4];
    long double p[4];
    
    ////    Excitation of the Recoil Nucleus
    long double Ex; // MeV
    
    ////    theta_ejectile is the angle of the Ejectile with respect to the beam axis
    long double theta_ejectile; // deg
    
    ////    theta_recoil is the angle of the Recoil with respect to the beam axis
    long double theta_recoil; // deg
    
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

long double CalculateRecoilPolarAngle(long double mass_projectile_amu, long double mass_target_amu, long double mass_ejectile_amu, long double mass_recoil_amu, long double beamEnergy_MeV, long double recoilExcitationEnergy_MeV, long double ejectilePolarAngle_labFrame_deg)
{
    long double m[4];
    long double T[4];
    long double E[4];
    long double p[4];
    
    ////    Excitation of the Recoil Nucleus
    long double Ex; // MeV
    
    ////    theta_ejectile is the angle of the Ejectile with respect to the beam axis
    long double theta_ejectile; // deg
    
    ////    theta_recoil is the angle of the Recoil with respect to the beam axis
    long double theta_recoil; // deg
    
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
    
    return theta_recoil;
}

void CalcEx(long double *m, long double *T, long double *E, long double *p, long double Xpos, long double theta, long double &Ex, std::vector<long double> momentumCalPars) {
    
    ////    Q-value calculation
    long double Q = (m[0] + m[1])*c2 - (m[2] + m[3])*c2; // MeV
    
    ////    Conversion of ejectile scattering angle from degrees to radians
    long double theta_lab = theta*0.017453292; // radians
    
    ////    Initial Total Energy Calculation
    E[0] = T[0] + (m[0]*c2);
    E[1] = T[1] + (m[1]*c2);
    
    ////    Initial Momentum Calculation
    p[0] = std::sqrt((E[0]*E[0]) - (m[0]*m[0]*c4))*(1.0/std::pow(c2, 0.5));
    p[1] = std::sqrt((E[1]*E[1]) - (m[1]*m[1]*c4))*(1.0/std::pow(c2, 0.5));
    
    p[2] = 0.0;
    for(int i=0; i<(int) momentumCalPars.size(); i++)
    {
        p[2] += momentumCalPars[i]*(1.0/std::pow(c2, 0.5))*std::pow(Xpos, i);
    }
    E[2] = std::sqrt(p[2]*p[2]*c2 + (m[2]*m[2]*c4));
    T[2] = E[2] - m[2]*c2;
    
    p[3] = std::sqrt((p[0]*p[0]) + (p[2]*p[2]) - 2*p[0]*p[2]*std::cos(theta_lab));
    E[3] = std::sqrt(p[3]*p[3]*c2 + (m[3]*m[3]*c4));
    T[3] = E[3] - m[3]*c2;
    
    Ex = E[0] + E[1] + Q - E[2] - E[3];
}




