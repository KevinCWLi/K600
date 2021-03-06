//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//      ----------------------------------------------------------------
//                      K600 Spectrometer (iThemba Labs)
//      ----------------------------------------------------------------
//
//      Github repository: https://www.github.com/KevinCWLi/K600
//
//      Main Author:    K.C.W. Li
//
//      email: likevincw@gmail.com
//

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "EventAction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include <TF1.h>
#include <TMath.h>

class G4ParticleGun;
class G4Event;
class EventAction;

/// The primary generator action class with particle gum.
///
/// It defines a single particle which hits the calorimeter
/// perpendicular to the input face. The type of the particle
/// can be changed via the G4 build-in commands of G4ParticleGun class
/// (see the macros provided with this example).

//static std::mutex mutex_particleN;  // protects particleN
static int particleN;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    PrimaryGeneratorAction(EventAction* eventAction);
    virtual ~PrimaryGeneratorAction();
    
    virtual void GeneratePrimaries(G4Event* event);
    
    // set methods
    void SetRandomFlag(G4bool value);
    
    G4double    *angDist_function;
    G4double    *angDist_theta;
    G4int       angDist_noOfPoints;
    
    
    void initialiseAngDist_interpolated(G4String name_angDist);
    
    G4double evaluateAngDist_interpolated(G4double chosenTheta);
    
    void CalculateBinaryRelativisticKinematics(long double theta_ejectile, long double phi_ejectile, long double &theta_recoil_LAB, long double &phi_recoil_LAB, long double &eX);

    void CalculateBinaryDecayKinematics(double eX, double eX_residual, double sepE, double initialScatteringReactionMomentum, double initialScatteringReactionEnergy, std::vector<double> decay_A, double theta_recoil_LAB, double phi_recoil_LAB, std::vector<double>& decay_theta_recoilCOM, std::vector<double>& decay_phi_recoilCOM, std::vector<double>& decay_theta_LAB, std::vector<double>& decay_phi_LAB, std::vector<double>& decay_T_LAB);
    void CalculateBinaryDecayKinematics2(long double eX, long double eX_residual, long double sepE, long double initialScatteringReactionMomentum, long double initialScatteringReactionEnergy, std::vector<long double> decay_A, long double theta_recoil_LAB, long double phi_recoil_LAB, std::vector<long double>& decay_theta_recoilCOM, std::vector<long double>& decay_phi_recoilCOM, std::vector<long double>& decay_theta_LAB, std::vector<long double>& decay_phi_LAB, std::vector<long double>& decay_T_LAB);

    std::vector<G4LorentzVector> CalculateBinaryDecayKinematics3(double eX, double eX_residual, double sepE, double initialScatteringReactionMomentum, double initialScatteringReactionEnergy, std::vector<double> decay_A, double theta_recoil_LAB, double phi_recoil_LAB, std::vector<double>& decay_theta_recoilCOM, std::vector<double>& decay_phi_recoilCOM, std::vector<double>& decay_theta_LAB, std::vector<double>& decay_phi_LAB, std::vector<double>& decay_T_LAB);

    std::vector<G4LorentzVector> CalculateReactionMomentumVectors(double *massEnergies, double TBeam, double ThetaAlphaCM, double PhiAlphaCM);
    
//    std::mutex angDist_mutex;

//    static G4ThreadLocal TF1* fAngCor_2plus_9870;
//    static G4ThreadLocal TF1* fAngDist_2plus_9870;

    TF1* fAngDist_1minus;

//    TF1* fAngCor_2plus_9870;
    TF1* fAngDist_2plus_9870;

//    TF1* fAngCor_2plus_16106;
//    TF1* fAngDist_2plus_16106;
//
    TF1* fAngDist_3minus;
    
    TF1* fAngDist_4plus;

//    TF1* fAngCor_4plus_13300;
//    TF1* fAngDist_4plus_13300;
//    
//    TF1* fAngCor_4plus_14079;
//    TF1* fAngDist_4plus_14079;
    
private:
    G4ParticleGun*  fParticleGun; // G4 particle gun
    EventAction*  fEventAction;
    
    long double m[4], T[4], E[4], p[4];
    long double eX;
    
    G4double    mx;
    G4double    my;
    G4double    mz;
    
    G4double    theta;
    G4double    phi;
    
    G4double    beamEnergy;
    G4double    energy;
    G4double    recoilExcitationEnergy;
    G4double    alphaSeperationEnergy;
    G4double    protonSeperationEnergy;
    G4double    daughterExcitationEnergy;
    G4int       alphaDecayMode;
    
    G4ThreeVector ejectileDirection;
    G4ThreeVector recoilDirection;
    
    //------------------------------------------------
    int nEnergies;
    int nParticlesPerEnergy;
    
    std::vector<G4double> nEventsPerKineticEnergy;
    std::vector<G4double> initialKineticEnergies;
    
    int GetParticleN();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline int PrimaryGeneratorAction::GetParticleN() {
    
//    std::lock_guard<std::mutex> lock(mutex_particleN);
    return particleN++;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double PrimaryGeneratorAction::evaluateAngDist_interpolated(G4double chosenTheta) {
    
    G4double result = 0.0;
    
    G4int index_L = (G4int) ((chosenTheta/180.0)*angDist_noOfPoints);
    G4int index_H = index_L + 1;
    
    G4double m, c;
    
    m = (angDist_function[index_H]-angDist_function[index_L])/(angDist_theta[index_H]-angDist_theta[index_L]);
    c = angDist_function[index_H] - (m*angDist_theta[index_H]);
    
    result = (m*chosenTheta) + c;
    
    ////    TEST
    //result = angDist_function[index_H];
    //result = sin(chosenTheta);
    //result = pow(result, 2.0);
    
    //result = angDist_function[index_H];
    
    
    /*
     G4cout << "" << G4endl;
     G4cout << "chosenTheta:    " << chosenTheta << G4endl;
     G4cout << "index_L:    " << index_L << G4endl;
     G4cout << "index_H:    " << index_H << G4endl;
     G4cout << "angDist_theta[index_L]:    " << angDist_theta[index_L] << G4endl;
     G4cout << "angDist_theta[index_H]:    " << angDist_theta[index_H] << G4endl;
     */
    
    //G4cout << "result:    " << result << G4endl;
    
    return result;
}

#endif


