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

#include "PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"

#include "G4IonTable.hh"

#include "BiRelKin.hh"

#include <TCanvas.h>

////    Distributions
#include "alpha0_0plus_0plus_12049_L0.h"
#include "alpha0_0plus_0plus_14032_L0.h"
#include "alpha0_0plus_0plus_15097_L0.h"
#include "alpha0_1minus_0plus_15097_L1.h"
#include "alpha0_2plus_0plus_11520_L2.h"
#include "alpha0_2plus_0plus_15097_L2.h"
#include "alpha0_3minus_0plus_15097_L3.h"
#include "alpha0_4plus_0plus_15097_L4.h"
#include "alpha0_5minus_0plus_15097_L5.h"
#include "alpha1_2plus_2plus_15097_L0.h"
#include "alpha1_3minus_2plus_15097_L1_L3_L5.h"
#include "alpha1_3minus_2plus_15097_L1.h"
#include "alpha1_3minus_2plus_15097_L3.h"
#include "alpha1_3minus_2plus_15097_L5.h"
#include "p0_2plus_12minus_14926_L1.h"
#include "p0_2plus_12minus_14926_L3.h"

#include <unistd.h>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//TF1* G4ThreadLocal fAngDist_2plus_9870;
//TF1* fAngDist_2plus_9870;

PrimaryGeneratorAction::PrimaryGeneratorAction(EventAction* eventAction): G4VUserPrimaryGeneratorAction(), fParticleGun(0), fEventAction(eventAction)
{
    
    ///////////////////////////////////////////////////////////////
    //          To generate radioactive decay - enabled particles
    ///////////////////////////////////////////////////////////////
    
    G4int n_particle = 1;
    fParticleGun  = new G4ParticleGun(n_particle);
    
    //G4ParticleDefinition* particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("He3");
    //G4ParticleDefinition* particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("neutron");
    //G4ParticleDefinition* particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("mu-");
//    G4ParticleDefinition* particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
    //G4ParticleDefinition* particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("e-");
//    G4ParticleDefinition* particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("alpha");
    //G4ParticleDefinition* particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("proton");
    G4ParticleDefinition* particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("geantino");
    
    fParticleGun->SetParticleDefinition(particleDefinition);
    
    //fParticleGun->SetParticleEnergy(0.*MeV);
    //fParticleGun->SetParticleEnergy(0.100*MeV);
    //fParticleGun->SetParticleEnergy(1.332*MeV);
    //fParticleGun->SetParticleEnergy(82.0*keV);
    
    //fParticleGun->SetParticleEnergy(200.*MeV);
    //fParticleGun->SetParticleEnergy(22.5*MeV);
    //fParticleGun->SetParticleEnergy(50.0*MeV);
    //fParticleGun->SetParticleEnergy(4*GeV);
    //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
    
    
    //fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.*m,1.90*m));
    //fParticleGun->SetParticlePosition(G4ThreeVector(0.,57.5*mm,1.6*m));
    //fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
    //fParticleGun->SetParticlePosition(G4ThreeVector(0.,57.5*mm,2.10*m));
    
    fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
    //fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,-1.*mm));
    //fParticleGun->SetParticlePosition(G4ThreeVector(2000.*mm,0.,3000.0*mm));
    
    //fParticleGun->SetParticlePosition(G4ThreeVector(-3.*m, 0., -3.8*m));
    
    
    /*
     ////////    4He, +1 charge
     G4int Z = 6, A = 12;
     //G4int Z = 2, A = 4;
     
     //G4double ionCharge   = 1.*eplus;
     G4double ionCharge   = 0.*eplus;
     G4double excitEnergy = 12.049*MeV;
     
     //G4ParticleDefinition* ion = G4ParticleTable::GetParticleTable()->GetIon(Z,A,excitEnergy);
     G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
     //ion->SetPDGLifeTime(1*ns);
     fParticleGun->SetParticleDefinition(ion);
     fParticleGun->SetParticleCharge(ionCharge);
     */
    
    /*
    ////////    4He, +1 charge
    G4int Z = 4, A = 8;
    G4double ionCharge   = 0.*eplus;
    G4double excitEnergy = 0.0*MeV;
    
    G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
    //ion->SetPDGLifeTime(1*ns);
    fParticleGun->SetParticleDefinition(ion);
    fParticleGun->SetParticleCharge(ionCharge);
    */
    
    
    ////========================================================////
    ////    Initialise a pre-calculated angular distribution
    G4String angDist_name = "a0_15097_16O";
    
    /*
    initialiseAngDist_interpolated(angDist_name);
    
    
    ////    Distributions - From AngCor
    initialiseAngCor_alpha0_0plus_0plus_12049_L0();
    initialiseAngCor_alpha0_0plus_0plus_14032_L0();
    initialiseAngCor_alpha0_0plus_0plus_15097_L0();
    initialiseAngCor_alpha0_1minus_0plus_15097_L1();
    initialiseAngCor_alpha0_2plus_0plus_11520_L2();
    initialiseAngCor_alpha0_2plus_0plus_15097_L2();
    initialiseAngCor_alpha0_3minus_0plus_15097_L3();
    initialiseAngCor_alpha0_4plus_0plus_15097_L4();
    initialiseAngCor_alpha0_5minus_0plus_15097_L5();
    initialiseAngCor_alpha1_2plus_2plus_15097_L0();
    initialiseAngCor_alpha1_3minus_2plus_15097_L1_L3_L5();
    initialiseAngCor_alpha1_3minus_2plus_15097_L1();
    initialiseAngCor_alpha1_3minus_2plus_15097_L3();
    initialiseAngCor_alpha1_3minus_2plus_15097_L5();
    initialiseAngCor_p0_2plus_12minus_14926_L1();
    initialiseAngCor_p0_2plus_12minus_14926_L3();
    */
    
    
    ////    Decay particle energy
    //G4double decayParticleEnergy = (11.520-7.16192)*(12.0/16.0);
    //fParticleGun->SetParticleEnergy(decayParticleEnergy*MeV);
    
//    std::cout << "CHECKPOINT, PrimaryGeneratorAction, 1" << std::endl;
    
    //----------------------------------------------------
    /*
    nEnergies += 1401;
    nParticlesPerEnergy = 100; // 20000000 particles to be simulated
    
    double energyMin = 6.0;
    double energyMax = 20.0;
    double energyDivision = ((energyMax-energyMin)/(nEnergies-1));
    
    for(int i=0; i<nEnergies; i++)
    {
        initialKineticEnergies.push_back(((i*energyDivision) + energyMin)*MeV);
    }
    */
    
    //----------------------------------------------------
    nEnergies = 40;
    nParticlesPerEnergy = 1.e+07;
    
    double energyMin = 500.0;
    double energyMax = 20000.0;
    double energyDivision = ((energyMax-energyMin)/(nEnergies-1));
    
    for(int i=0; i<nEnergies; i++)
    {
        initialKineticEnergies.push_back(((i*energyDivision) + energyMin)*keV);
    }

//    std::cout << "CHECKPOINT, PrimaryGeneratorAction, 2" << std::endl;

    double xMin = 0.0;
    double xMax = 180.0;

//    fAngDist_2plus_9870 = new TF1("fAngDist_2plus_9870", "[0]*TMath::Power(0.5*(3.0*TMath::Power(TMath::Cos(x*0.0174533), 2.0) - 1.0), 2.0)*TMath::Sin(x*0.0174533)", xMin, xMax);
//    fAngDist_2plus_9870->SetNpx(1000);
//    fAngDist_2plus_9870->SetParameter(0, 1.0);

    //------------------------------------------------
//    fAngCor_2plus_9870 = new TF1("fAngCor_2plus_9870",
//                                 [&](double*x, double *p)
//                                 {
//                                     return p[0]*TMath::Power(0.5*(3.0*TMath::Power(TMath::Cos(x[0]*0.0174533), 2.0) - 1.0), 2.0);
//                                 },
//                                 xMin, xMax, 1);
//    fAngCor_2plus_9870->SetNpx(1000);
//    fAngCor_2plus_9870->SetParameter(0, 1.0);

//    Double_t w = 600;
//    Double_t h = 600;
//    auto c1 = new TCanvas("c1", "c1", w, h);
//    
    //------------------------------------------------
    fAngDist_1minus = new TF1("fAngDist_1minus",
                                  [&](double*x, double *p)
                                  {
                                      //                                     return p[0]*fAngCor_2plus_9870->Eval(x[0])*TMath::Sin(x[0]*0.0174533);
                                      return p[0]*TMath::Power(TMath::Cos(x[0]*0.0174533), 2.0)*TMath::Sin(x[0]*0.0174533);
                                  },
                                  xMin, xMax, 1);
    fAngDist_1minus->SetNpx(1000);
    fAngDist_1minus->SetParameter(0, 1.0);

    G4cout << "fAngDist_1minus->GetMaximum(): " << fAngDist_1minus->GetMaximum() << G4endl;
    
//    fAngDist_1minus->Draw();
//    c1->Print("test.root");
    
//    fAngDist_1minus->SaveAs("fAngDist_1minus.root");
    
    //------------------------------------------------
    fAngDist_2plus_9870 = new TF1("fAngDist_2plus_9870",
                                 [&](double*x, double *p)
                                 {
//                                     return p[0]*fAngCor_2plus_9870->Eval(x[0])*TMath::Sin(x[0]*0.0174533);
                                     return p[0]*TMath::Power(0.5*(3.0*TMath::Power(TMath::Cos(x[0]*0.0174533), 2.0) - 1.0), 2.0)*TMath::Sin(x[0]*0.0174533);
                                 },
                                 xMin, xMax, 1);
    fAngDist_2plus_9870->SetNpx(1000);
    fAngDist_2plus_9870->SetParameter(0, 1.0);

    G4cout << "fAngDist_2plus_9870->GetMaximum(): " << fAngDist_2plus_9870->GetMaximum() << G4endl;

    //------------------------------------------------
    fAngDist_3minus = new TF1("fAngDist_3minus",
                              [&](double*x, double *p)
                              {
                                  //                                     return p[0]*fAngCor_2plus_9870->Eval(x[0])*TMath::Sin(x[0]*0.0174533);
                                  return p[0]*TMath::Power(0.5*(5.0*TMath::Power(TMath::Cos(x[0]*0.0174533), 2.0) - 3.0*TMath::Cos(x[0]*0.0174533)), 2.0)*TMath::Sin(x[0]*0.0174533);
                              },
                              xMin, xMax, 1);
    fAngDist_3minus->SetNpx(1000);
    fAngDist_3minus->SetParameter(0, 1.0);

    G4cout << "fAngDist_3minus->GetMaximum(): " << fAngDist_3minus->GetMaximum() << G4endl;

    //------------------------------------------------
    fAngDist_4plus = new TF1("fAngDist_4plus",
                              [&](double*x, double *p)
                              {
                                  //                                     return p[0]*fAngCor_2plus_9870->Eval(x[0])*TMath::Sin(x[0]*0.0174533);
                                  return p[0]*TMath::Power((1.0/8.0)*(35.0*TMath::Power(TMath::Cos(x[0]*0.0174533), 4.0) - 30.0*TMath::Power(TMath::Cos(x[0]*0.0174533), 2.0) + 3.0), 2.0)*TMath::Sin(x[0]*0.0174533);
                              },
                              xMin, xMax, 1);
    fAngDist_4plus->SetNpx(1000);
    fAngDist_4plus->SetParameter(0, 1.0);
    
    G4cout << "fAngDist_4plus->GetMaximum(): " << fAngDist_4plus->GetMaximum() << G4endl;

    
//    delete c1;
    
    //    fAngCor_2plus_9870->SaveAs("test1.root");
    //    fAngDist_2plus_9870->SaveAs("test2.root");

    //------------------------------------------------
//    fAngCor_2plus_16106 = new TF1("fAngCor_2plus_16106",
//                                  [&](double*x, double *p)
//                                  {
//                                      return p[0]*TMath::Power(0.5*(3.0*TMath::Power(TMath::Cos(x[0]*0.0174533), 2.0) - 1.0), 2.0);
//                                  },
//                                  xMin, xMax, 1);
//    fAngCor_2plus_16106->SetNpx(1000);
//    fAngCor_2plus_16106->SetParameter(0, 1.0);
//    
//    fAngDist_2plus_16106 = new TF1("fAngDist_2plus_16106",
//                                   [&](double*x, double *p)
//                                   {
//                                       return p[0]*fAngCor_2plus_16106->Eval(x[0])*TMath::Sin(x[0]*0.0174533);
//                                   },
//                                   xMin, xMax, 1);
//    fAngDist_2plus_16106->SetNpx(1000);
//    fAngDist_2plus_16106->SetParameter(0, 1.0);

    //------------------------------------------------
//    fAngCor_4plus_13300 = new TF1("fAngCor_4plus_13300",
//                                 [&](double*x, double *p)
//                                 {
//                                     return p[0]*TMath::Power((1.0/8.0)*(35.0*TMath::Power(TMath::Cos(x[0]*0.0174533), 4.0) - 30.0*TMath::Power(TMath::Cos(x[0]*0.0174533), 2.0) + 3.0), 2.0);
//                                 },
//                                 xMin, xMax, 1);
//    fAngCor_4plus_13300->SetNpx(1000);
//    fAngCor_4plus_13300->SetParameter(0, 1.0);
//    
//    fAngDist_4plus_13300 = new TF1("fAngDist_4plus_13300",
//                                  [&](double*x, double *p)
//                                  {
//                                      return p[0]*fAngCor_4plus_13300->Eval(x[0])*TMath::Sin(x[0]*0.0174533);
//                                  },
//                                  xMin, xMax, 1);
//    fAngDist_4plus_13300->SetNpx(1000);
//    fAngDist_4plus_13300->SetParameter(0, 1.0);
    
//    fAngCor_4plus_13300->SaveAs("test1.root");
//    fAngDist_4plus_13300->SaveAs("test2.root");

    //------------------------------------------------
//    fAngCor_4plus_14079 = new TF1("fAngCor_4plus_14079",
//                                  [&](double*x, double *p)
//                                  {
//                                      return p[0]*TMath::Power((1.0/8.0)*(35.0*TMath::Power(TMath::Cos(x[0]*0.0174533), 4.0) - 30.0*TMath::Power(TMath::Cos(x[0]*0.0174533), 2.0) + 3.0), 2.0);
//                                  },
//                                  xMin, xMax, 1);
//    fAngCor_4plus_14079->SetNpx(1000);
//    fAngCor_4plus_14079->SetParameter(0, 1.0);
//    
//    fAngDist_4plus_14079 = new TF1("fAngDist_4plus_14079",
//                                   [&](double*x, double *p)
//                                   {
//                                       return p[0]*fAngCor_4plus_14079->Eval(x[0])*TMath::Sin(x[0]*0.0174533);
//                                   },
//                                   xMin, xMax, 1);
//    fAngDist_4plus_14079->SetNpx(1000);
//    fAngDist_4plus_14079->SetParameter(0, 1.0);

//    fAngCor_4plus_14079->SaveAs("test1.root");
//    fAngDist_4plus_14079->SaveAs("test2.root");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
//    std::cout << "fAngCor_4plus_14079->Eval(10.0): " << fAngCor_4plus_14079->Eval(10.0) << std::endl;
    
    /*
    ////////    8B3, +1 charge
    G4int Z = 4, A = 8;
    G4double ionCharge   = 0.*eplus;
    G4double excitEnergy = 0.0*MeV;
    
    G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(Z,A,0);
    //ion->SetPDGLifeTime(1*ns);
    fParticleGun->SetParticleDefinition(ion);
    fParticleGun->SetParticleCharge(ionCharge);
    */
    
    ////////////////////////////////////////////////////////////////
    ////    Definition of states/resonances for 9Be(3He, t)9B
    /*
     G4int noOfStates_9Be = 5;
     G4double resonanceEnergy_9B[noOfStates_9Be]; // MeV
     G4double resonanceWidth_9B[noOfStates_9Be]; // MeV
     
     ////
     resonanceEnergy_9B[0] = 1.5;
     resonanceWidth_9B[0] = 1.2;
     
     ////    5/2-
     resonanceEnergy_9B[1] = 2.345;
     resonanceWidth_9B[1] = 0.081;
     
     ////    5/2+
     resonanceEnergy_9B[2] = 2.751;
     resonanceWidth_9B[2] = 0.614;
     
     ////    1/2-
     resonanceEnergy_9B[3] = 2.78;
     resonanceWidth_9B[3] = 3.13;
     
     ////
     resonanceEnergy_9B[4] = 4.8;
     resonanceWidth_9B[4] = 1.2;
     */
    
    
    ////////////////////////////////////////////////////////////////
    ////    Definition of states/resonances for 9Be(3He, t)9B
    /*
     G4int noOfStates_9Be = 5;
     G4double resonanceEnergy_9B[noOfStates_9Be]; // MeV
     G4double resonanceWidth_9B[noOfStates_9Be]; // MeV
     
     ////
     resonanceEnergy_9B[0] = 1.5;
     resonanceWidth_9B[0] = 1.2;
     
     ////    5/2-
     resonanceEnergy_9B[1] = 2.345;
     resonanceWidth_9B[1] = 0.081;
     
     ////    5/2+
     resonanceEnergy_9B[2] = 2.751;
     resonanceWidth_9B[2] = 0.614;
     
     ////    1/2-
     resonanceEnergy_9B[3] = 2.78;
     resonanceWidth_9B[3] = 3.13;
     
     ////
     resonanceEnergy_9B[4] = 4.8;
     resonanceWidth_9B[4] = 1.2;
     */
    
    
    
    ////    Different particles for Visualisation
    /*
     double test = G4UniformRand();
     
     G4ParticleDefinition* particleDefinitionNew;
     if(test<0.3) particleDefinitionNew = G4ParticleTable::GetParticleTable()->FindParticle("alpha");
     if(test>=0.3 && test<0.6) particleDefinitionNew = G4ParticleTable::GetParticleTable()->FindParticle("e-");
     if(test>=0.6 && test<=1.0) particleDefinitionNew = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
     fParticleGun->SetParticleDefinition(particleDefinitionNew);
     */
    
    ///////////////////////////////////////////////////////////////
    //          To generate radioactive decay - enabled particles
    ///////////////////////////////////////////////////////////////
    
    /*
     //G4int Z = 3, A = 11;
     //G4int Z = 6, A = 21;
     //G4int Z = 10, A = 18;
     //G4int Z = 11, A = 22;
     //G4int Z = 27, A = 60;
     G4int Z = 8, A = 16;
     //G4int Z = 63, A = 152;
     
     G4double ionCharge   = 0.*eplus;
     G4double excitEnergy = 12.049*MeV;
     
     //G4ParticleDefinition* ion = G4ParticleTable::GetParticleTable()->GetIon(Z,A,excitEnergy);
     G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
     //ion->SetPDGLifeTime(1*ns);
     ion->SetPDGLifeTime(0.5485101634774202e-10*ns);
     fParticleGun->SetParticleDefinition(ion);
     fParticleGun->SetParticleCharge(ionCharge);
     */
    
    
    ///////////////////////////////////////////////////////////
    //       Initial Position Distribution of the Particle
    ///////////////////////////////////////////////////////////
    /*
     G4double targetZoffset = -1000.0; // um
     G4double targetThickness = 2.42; // um
     G4double InitialPosition = G4RandFlat::shoot(-(targetThickness/2) + targetZoffset, (targetThickness/2) + targetZoffset); // um
     */
    
    /*
     G4double targetZoffset = 0.0; // um
     G4double targetThickness = 1.0; // um
     G4double InitialPosition = G4RandFlat::shoot(-(targetThickness/2) + targetZoffset, (targetThickness/2) + targetZoffset); // um
     fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,InitialPosition*um));
     */
    
    ///////////////////////////////////////////////////////////
    //       Initial Direction of the Particle
    ///////////////////////////////////////////////////////////
    /*
     G4double theta = 2*M_PI*G4UniformRand();
     
     mz = -1.0 + 2*G4UniformRand();
     
     ////    Backward Angles
     //mz = -1.0 + G4UniformRand();
     
     ////    Forward Angles
     //mz = 1.0 - G4UniformRand();
     
     ////    Very Forward Angles
     //mz = 0.001 - 0.001*G4UniformRand() + (1.0-0.001);
     
     G4double a = sqrt(1-(mz*mz));
     
     mx = a*cos(theta);
     my = a*sin(theta);
     
     fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
     */
    

    
    
    ////////////////////////////////////////////////////
    ////    ISOTROPIC - Inverse Transform Method
    
    /*
     G4double theta = acos(1 - (2.0*G4UniformRand()))/deg; // 0.0->180.0
     //G4double theta = acos(1 - (1.0*G4UniformRand()))/deg; // 0.0->90.0 deg
     //G4double theta = acos(1 - (1.0*G4UniformRand() + 1.0))/deg; // 90.0->180.0 deg
     
     ////    Collimator - 2.0 deg range
     //G4double theta = acos(1 - (6.09173503004822869e-04*G4UniformRand()))/deg; // deg
     
     G4double phi = 360.0*G4UniformRand(); // deg
     
     //------------------------------------------------
//     fEventAction->SetInitialTheta_COM(theta);
//     fEventAction->SetInitialPhi_COM(phi);
//     
//     fEventAction->SetInitialTheta_LAB(theta);
//     fEventAction->SetInitialPhi_LAB(phi);
    
     //------------------------------------------------
     mx = sin(theta*deg)*cos(phi*deg);
     my = sin(theta*deg)*sin(phi*deg);
     mz = cos(theta*deg);
     
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
    fParticleGun->SetParticleEnergy(1.0*MeV);

    fParticleGun->GeneratePrimaryVertex(anEvent);
    */
    
    
    /*
     void SetInitialTheta_LAB(G4double theta) {initialTheta_LAB = theta;};
     void SetInitialPhi_LAB(G4double phi) {initialPhi_LAB = phi;};
     void SetInitialTheta_COM(G4double theta) {initialTheta_COM = theta;};
     void SetInitialPhi_COM(G4double phi) {initialPhi_COM = phi;};
     */
    
    
//    std::vector<double> decay_theta_LAB;
//    
//    for(int i=0; i<5; i++)
//    {
//        G4double theta = acos(1 - (2.0*G4UniformRand()))/deg; // 0.0->180.0
//        //G4double theta = acos(1 - (1.0*G4UniformRand()))/deg; // 0.0->90.0 deg
//        //G4double theta = acos(1 - (1.0*G4UniformRand() + 1.0))/deg; // 90.0->180.0 deg
//     
//        decay_theta_LAB.push_back(theta);
//        
//        ////    Collimator - 2.0 deg range
//        //G4double theta = acos(1 - (6.09173503004822869e-04*G4UniformRand()))/deg; // deg
//        
//        G4double phi = 360.0*G4UniformRand(); // deg
//        
//        //------------------------------------------------
//        mx = sin(theta*deg)*cos(phi*deg);
//        my = sin(theta*deg)*sin(phi*deg);
//        mz = cos(theta*deg);
//        
//        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
//        fParticleGun->GeneratePrimaryVertex(anEvent);
//    }
//    
//    fEventAction->SetDecay_theta_LAB(decay_theta_LAB);

    
    ////////////////////////////////////////////////////
    ////    X17, gamma-gamma
    
//    G4double theta = acos(1 - (2.0*G4UniformRand()))/deg; // 0.0->180.0
//    G4double phi = 360.0*G4UniformRand(); // deg
//    
//    //------------------------------------------------
//    //     fEventAction->SetInitialTheta_COM(theta);
//    //     fEventAction->SetInitialPhi_COM(phi);
//    //
//    //     fEventAction->SetInitialTheta_LAB(theta);
//    //     fEventAction->SetInitialPhi_LAB(phi);
//    
//    //------------------------------------------------
//    mx = sin(theta*deg)*cos(phi*deg);
//    my = sin(theta*deg)*sin(phi*deg);
//    mz = cos(theta*deg);
//    
//    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
//    
//    //------------------------------------------------
//    std::vector<double> gammaRayEnergies;
//    double e1 = 18150.0*G4UniformRand();
//    double e2 = 18150.0 - e1;
//    
//    gammaRayEnergies.push_back(e1);
//    gammaRayEnergies.push_back(e2);
//    fEventAction->SetReaction_Ex(gammaRayEnergies);
//    
//    fParticleGun->SetParticleEnergy(e1*keV);
//    fParticleGun->GeneratePrimaryVertex(anEvent);
//
//    
//    //------------------------------------------------
//    theta = acos(1 - (2.0*G4UniformRand()))/deg; // 0.0->180.0
//    phi = 360.0*G4UniformRand(); // deg
//    
//    mx = sin(theta*deg)*cos(phi*deg);
//    my = sin(theta*deg)*sin(phi*deg);
//    mz = cos(theta*deg);
//    
//    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
//    
//    //------------------------------------------------
//    fParticleGun->SetParticleEnergy(e2*keV);
//    fParticleGun->GeneratePrimaryVertex(anEvent);

    
    
    ////////////////////////////////////////////////////
    ////    Efficiency
    /*
    static int eventN = 0;
    static G4double gammaEnergy = 0.0*keV;
    
    double theta = acos(1 - (2.0*G4UniformRand()))/deg; // 0.0->180.0
    double phi = 360.0*G4UniformRand(); // deg
    
    //------------------------------------------------
    //     fEventAction->SetInitialTheta_COM(theta);
    //     fEventAction->SetInitialPhi_COM(phi);
    //
    //     fEventAction->SetInitialTheta_LAB(theta);
    //     fEventAction->SetInitialPhi_LAB(phi);
    
    //------------------------------------------------
    mx = sin(theta*deg)*cos(phi*deg);
    my = sin(theta*deg)*sin(phi*deg);
    mz = cos(theta*deg);
    
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));

    //------------------------------------------------
//    if(eventN%(20000)==0)
//    {
//        gammaEnergy += 500.0*keV;
//    }

    G4double initialParticleKineticEnergy = 0.0*keV;
    
    int energyN = (int) GetParticleN()/nParticlesPerEnergy;
    
    if(energyN<((int) initialKineticEnergies.size()))
    {
        gammaEnergy = initialKineticEnergies[energyN];
    }
    else
    {
        gammaEnergy = 20000.0*keV;
    }

    //------------------------------------------------
    std::vector<double> gammaRayEnergies;
    gammaRayEnergies.push_back(gammaEnergy/keV);
    fEventAction->SetReaction_Ex(gammaRayEnergies);

    fParticleGun->SetParticleEnergy(gammaEnergy);
    fParticleGun->GeneratePrimaryVertex(anEvent);
    
    eventN++;
    */
    

    /*
    double targetThickness = 3.0*1.1197459; // um
    double reactionPosition = targetThickness*G4UniformRand() - targetThickness/2.0;
    fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,reactionPosition*um));

    //========================================================================================================================
    //========================================================================================================================
    //      X17
    
    int nReactionProducts = 4;
    std::vector<int>    reaction_Z{3, 1, 0, 4}; // e
    std::vector<long double> reaction_A{7.01600342665, 1.00782503223, 16.7/931.49410242, 8.005305102}; // amu
    auto reaction_P = std::vector<long double>(4, 0.0);
    std::vector<long double> reaction_T{7.13079, 0.0, 0.0, 0.0}; // MeV
    std::vector<long double> reaction_Ex{0.0, 0.0, 0.0, 0.0}; // MeV
    auto reaction_E = std::vector<long double>(4, 0.0);
    

    //------------------------------------------------------------------------------------
    //      Polar angles (LAB)
//    std::vector<double> reaction_theta_LAB{0.0, 0.0, acos(1 - (6.09173503004822869e-04*G4UniformRand()))/deg, 0.0}; // deg
    std::vector<long double> reaction_theta_LAB{0.0, 0.0, 180.0*G4UniformRand(), 0.0}; // deg
    auto reaction_theta_reactionCOM = std::vector<double>(4, 0.0);
    
    //------------------------------------------------------------------------------------
    //      Azimuthal angles (LAB)
    double reaction_ejectile_phi_LAB = 360.0*G4UniformRand(); // deg
    double reaction_recoil_phi_LAB = reaction_ejectile_phi_LAB - 180.0; // deg
    if(reaction_recoil_phi_LAB<0.0)
    {
        reaction_recoil_phi_LAB += 360.0;
    }
    
    std::vector<double> reaction_phi_LAB{0.0, 0.0, reaction_ejectile_phi_LAB, reaction_recoil_phi_LAB}; // deg
    auto reaction_phi_reactionCOM = std::vector<double>(4, 0.0);
    
    //------------------------------------------------------------------------------------
    //      Recoil Nucleux Excitation Energy
    //double recoilExcitationEnergy = 0.0; // MeV
    double recoilExcitationEnergy = 0.0; // MeV
    //double recoilExcitationEnergy = 7.36659; // MeV
    //double recoilExcitationEnergy = 7.36659 + 0.5; // MeV
    //double recoilExcitationEnergy = 12.0; // MeV
    //double recoilExcitationEnergy = 20.0*G4UniformRand(); // MeV
    
    double daughterNucleusFinalStateEnergy = 0.0; // MeV
    
    //------------------------------------------------------------------------------------
    //      Executing the Binary Relativistic Kinematics code
    BiRelKin(&reaction_A[0], &reaction_T[0], &reaction_E[0], &reaction_P[0], reaction_theta_LAB[2], reaction_theta_LAB[3], recoilExcitationEnergy+daughterNucleusFinalStateEnergy);

//    G4cout << "reaction_T[2]: " << reaction_T[2] << G4endl;
//    G4cout << "reaction_T[3]: " << reaction_T[3] << G4endl;
    
    reaction_Ex[3] = recoilExcitationEnergy;
    
    //------------------------------------------------------------------------------------
//    std::vector<double>(reaction_Ex.begin(), reaction_Ex.end());
    
    //      Filling the vectors in the EventAction object
    fEventAction->SetNReactionProducts(nReactionProducts);
//    fEventAction->SetReaction_Z(reaction_Z);
//    fEventAction->SetReaction_A(reaction_A);
//    fEventAction->SetReaction_P(reaction_P);
//    fEventAction->SetReaction_T(reaction_T);
//    fEventAction->SetReaction_Ex(reaction_Ex);
    fEventAction->SetReaction_Z(std::vector<int>(reaction_Z.begin(), reaction_Z.end()));
    fEventAction->SetReaction_A(std::vector<double>(reaction_A.begin(), reaction_A.end()));
    fEventAction->SetReaction_P(std::vector<double>(reaction_P.begin(), reaction_P.end()));
    fEventAction->SetReaction_T(std::vector<double>(reaction_T.begin(), reaction_T.end()));
    fEventAction->SetReaction_Ex(std::vector<double>(reaction_Ex.begin(), reaction_Ex.end()));
    
    fEventAction->SetReaction_theta_LAB(std::vector<double>(reaction_theta_LAB.begin(), reaction_theta_LAB.end()));
    fEventAction->SetReaction_phi_LAB(reaction_phi_LAB);
    fEventAction->SetReaction_theta_reactionCOM(reaction_theta_reactionCOM);
    fEventAction->SetReaction_phi_reactionCOM(reaction_phi_reactionCOM);
    
  
    //========================================================================================
    //      8Be Decay (0, 1)
    G4ParticleDefinition* particleDefinition_alpha = G4ParticleTable::GetParticleTable()->FindParticle("alpha");
    fParticleGun->SetParticleDefinition(particleDefinition_alpha);

    int nDecayProducts_8Be = 2;
    std::vector<int>    decay8Be_Z{2, 2}; // e
    std::vector<double> decay8Be_A{4.00260325413, 4.00260325413}; // amu
    auto decay8Be_P = std::vector<double>(2, 0.0);
    auto decay8Be_T = std::vector<double>(2, 0.0);
    
    auto decay8Be_theta_LAB = std::vector<double>(2, 0.0);
    auto decay8Be_phi_LAB = std::vector<double>(2, 0.0);
    auto decay8Be_theta_recoilCOM = std::vector<double>(2, 0.0);
    auto decay8Be_phi_recoilCOM = std::vector<double>(2, 0.0);
    auto decay8Be_theta_reactionCOM = std::vector<double>(2, 0.0);
    auto decay8Be_phi_reactionCOM = std::vector<double>(2, 0.0);

    //----------------------------------------
    //      Polar angles (recoilCOM)
    decay8Be_theta_recoilCOM[0] = acos(1 - (2.0*G4UniformRand()))/deg;
    decay8Be_theta_recoilCOM[1] = 180.0 - decay8Be_theta_recoilCOM[0];

    //----------------------------------------
    //      Azimuthal angles (LAB)
    decay8Be_phi_recoilCOM[0] = 360.0*G4UniformRand(); // deg
    decay8Be_phi_recoilCOM[1] = decay8Be_phi_recoilCOM[0] - 180.0; // deg
    
    if(decay8Be_phi_recoilCOM[1]<0.0)
    {
        decay8Be_phi_recoilCOM[1] += 360.0;
    }
    
    //----------------------------------------
    CalculateBinaryDecayKinematics(recoilExcitationEnergy, daughterNucleusFinalStateEnergy, -0.09184, reaction_P[3], reaction_E[3], decay8Be_A, reaction_theta_LAB[3], reaction_phi_LAB[3], decay8Be_theta_recoilCOM, decay8Be_phi_recoilCOM, decay8Be_theta_LAB, decay8Be_phi_LAB, decay8Be_T);

    double mx, my, mz;
    mx = sin(decay8Be_theta_LAB[0]*deg)*cos(decay8Be_phi_LAB[0]*deg);
    my = sin(decay8Be_theta_LAB[0]*deg)*sin(decay8Be_phi_LAB[0]*deg);
    mz = cos(decay8Be_theta_LAB[0]*deg);
    
    fEventAction->SetInitialParticleKineticEnergy(decay8Be_T[0]*MeV); // MeV
    fParticleGun->SetParticleEnergy(decay8Be_T[0]*MeV);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
//    G4cout << "decay8Be_T[0]: " << decay8Be_T[0] << G4endl;
//    G4cout << "decay8Be_theta_LAB[0]: " << decay8Be_theta_LAB[0] << G4endl;
    
    fParticleGun->GeneratePrimaryVertex(anEvent);
    
    //----------------------------------------
    mx = sin(decay8Be_theta_LAB[1]*deg)*cos(decay8Be_phi_LAB[1]*deg);
    my = sin(decay8Be_theta_LAB[1]*deg)*sin(decay8Be_phi_LAB[1]*deg);
    mz = cos(decay8Be_theta_LAB[1]*deg);
    
    fEventAction->SetInitialParticleKineticEnergy(decay8Be_T[1]*MeV); // MeV
    fParticleGun->SetParticleEnergy(decay8Be_T[1]*MeV);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
//    G4cout << "decay8Be_T[1]: " << decay8Be_T[1] << G4endl;
//    G4cout << "decay8Be_theta_LAB[1]: " << decay8Be_theta_LAB[1] << G4endl;
    
    fParticleGun->GeneratePrimaryVertex(anEvent);

    
    fEventAction->SetNDecayProducts(nDecayProducts_8Be);
    fEventAction->SetDecay_Z(decay8Be_Z);
    fEventAction->SetDecay_A(decay8Be_A);
    fEventAction->SetDecay_P(decay8Be_P);
    fEventAction->SetDecay_T(decay8Be_T);
    
    fEventAction->SetDecay_theta_LAB(decay8Be_theta_LAB);
    fEventAction->SetDecay_phi_LAB(decay8Be_phi_LAB);
    fEventAction->SetDecay_theta_reactionCOM(decay8Be_theta_reactionCOM);
    fEventAction->SetDecay_phi_reactionCOM(decay8Be_phi_reactionCOM);
    fEventAction->SetDecay_theta_recoilCOM(decay8Be_theta_recoilCOM);
    fEventAction->SetDecay_phi_recoilCOM(decay8Be_phi_recoilCOM);
    */
    
    
    
    
    
//    //========================================================================================
//    double m[4];
//    
//    m[0] = 7.01600342665;
//    m[1] = 1.00782503223;
//    m[2] = 8.005305102;
//    m[3] = 16.7/931.49410242;
//    
//    ////    Transforming the angles of the differential cross section from COM to LAB
//    double c2 = 931.49410242; // MeV/u, c^2
//    double Q = (m[0] + m[1])*c2 - (m[2] + m[3])*c2;
//    
//    //------------------------------------------------
//    double MassEnergy[4];
//    
//    for(int i=0; i<4; i++)
//    {
//        MassEnergy[i] = m[i]*931.49410242;
//    }
//    
//    double beamEnergy = 7.13079; // MeV
//    
//    auto reactionMomentumVectors = CalculateReactionMomentumVectors(MassEnergy, beamEnergy, 0.0, 2.0, 0.0);
//    double polarAngle_LAB = reactionMomentumVectors[0].theta()*180./M_PI;
//    double energy_LAB = reactionMomentumVectors[0].e();
//    double kineticEnergy_LAB = energy_LAB - MassEnergy[2];
//    
//    std::cout << "polarAngle_LAB: " << polarAngle_LAB << std::endl;
//    std::cout << "energy_LAB: " << energy_LAB << std::endl;
//    std::cout << "kineticEnergy_LAB: " << kineticEnergy_LAB << std::endl;

    
    
    //========================================================================================================================
    //      X17
    /*
    double c2 = 931.49410242; // MeV/u, c^2

//    double targetThickness = 3.0*1.1197459; // um
    double targetThickness = 50.0; // um
    double reactionPosition = targetThickness*G4UniformRand() - targetThickness/2.0;
    fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,reactionPosition*um));
    
    int nReactionProducts = 4;
    std::vector<int>    reaction_Z{3, 1, 0, 4}; // e
    std::vector<double> reaction_A{7.01600342665, 1.00782503223, 8.005305102, 16.7/931.49410242}; // amu
    std::vector<double> reaction_P{0.0, 0.0, 0.0, 0.0};
    std::vector<double> reaction_T{7.13079, 0.0, 0.0, 0.0}; // MeV
    std::vector<double> reaction_E{0.0, 0.0, 0.0, 0.0};
    std::vector<double> reaction_Ex{0.0, 0.0, 0.0, 0.0}; // MeV
    
    //------------------------------------------------------------------------------------
    //      Polar angles (LAB)
    double reaction_ejectile_theta_COM = acos(1 - (2.0*G4UniformRand()))/deg; // deg
    reaction_ejectile_theta_COM = 0.0; // deg
    
    std::vector<double> reaction_theta_reactionCOM{0.0, 0.0, reaction_ejectile_theta_COM, 0.0};
    std::vector<double> reaction_theta_LAB{0.0, 0.0, 0.0, 0.0}; // deg
    
    //------------------------------------------------------------------------------------
    //      Azimuthal angles (LAB)
    double reaction_ejectile_phi_COM = 360.0*G4UniformRand(); // deg
    double reaction_recoil_phi_COM = reaction_ejectile_phi_COM - 180.0; // deg
    if(reaction_recoil_phi_COM<0.0)
    {
        reaction_recoil_phi_COM += 360.0;
    }
    
    std::vector<double> reaction_phi_reactionCOM{0.0, 0.0, reaction_ejectile_phi_COM, reaction_recoil_phi_COM};
    std::vector<double> reaction_phi_LAB{0.0, 0.0, 0.0, 0.0}; // deg
    
    //------------------------------------------------------------------------------------
    //      Recoil Nucleux Excitation Energy
    //double recoilExcitationEnergy = 0.0; // MeV
    double recoilExcitationEnergy = 0.0; // MeV
    //double recoilExcitationEnergy = 7.36659; // MeV
    //double recoilExcitationEnergy = 7.36659 + 0.5; // MeV
    //double recoilExcitationEnergy = 12.0; // MeV
    //double recoilExcitationEnergy = 20.0*G4UniformRand(); // MeV
    
    double daughterNucleusFinalStateEnergy = 0.0; // MeV
    
    //========================================================================================
//    double m[4];
//    
//    m[0] = 7.01600342665;
//    m[1] = 1.00782503223;
//    m[2] = 8.005305102;
//    m[3] = 16.7/931.49410242;
    
    ////    Transforming the angles of the differential cross section from COM to LAB
//    double Q = (m[0] + m[1])*c2 - (m[2] + m[3])*c2;
    
    //------------------------------------------------
    double reactionMassEnergies[4];
    
    for(int i=0; i<4; i++)
    {
        reactionMassEnergies[i] = reaction_A[i]*931.49410242;
    }
    
    double beamEnergy = reaction_T[0]; // MeV
    
    auto reactionMomentumVectors = CalculateReactionMomentumVectors(reactionMassEnergies, beamEnergy, reaction_theta_reactionCOM[2], reaction_phi_reactionCOM[2]);
    
    for(int i=0; i<static_cast<int>(reactionMomentumVectors.size()); i++)
    {
        reaction_E[i] = reactionMomentumVectors[i].e();
        reaction_T[i] = reaction_E[i] - reactionMassEnergies[i];
        reaction_P[i] = std::sqrt(std::pow(reaction_E[i], 2.0) - std::pow(reactionMassEnergies[i], 2.0))*(1.0/std::pow(c2, 0.5));
//        std::cout << "reaction_P[i]: " << reaction_P[i] << std::endl;
//        std::cout << "reaction_E[i]: " << reaction_E[i] << std::endl;
//        std::cout << "reaction_T[i]: " << reaction_T[i] << std::endl;
        
//        double polarAngle_LAB = reactionMomentumVectors[2].theta()*180./M_PI;
//        std::vector<double> reaction_theta_LAB{0.0, 0.0, 0.0, 0.0}; // deg

        reaction_theta_LAB[i] = reactionMomentumVectors[i].theta()*180./M_PI;
    }
    
//    std::cout << std::endl;

    
//    double polarAngle_LAB = reactionMomentumVectors[2].theta()*180./M_PI;
//    double energy_LAB = reactionMomentumVectors[2].e();
//    double kineticEnergy_LAB = energy_LAB - MassEnergy[2];
    
//    std::cout << "polarAngle_LAB: " << polarAngle_LAB << std::endl;
//    std::cout << "energy_LAB: " << energy_LAB << std::endl;
//    std::cout << "kineticEnergy_LAB: " << kineticEnergy_LAB << std::endl;
    
    //------------------------------------------------------------------------------------
    //    std::vector<double>(reaction_Ex.begin(), reaction_Ex.end());
    
    //      Filling the vectors in the EventAction object
//    fEventAction->SetNReactionProducts(nReactionProducts);
//    
//    fEventAction->SetReaction_Z(reaction_Z);
//    fEventAction->SetReaction_A(reaction_A);
//    fEventAction->SetReaction_T(reaction_T);
//    
//    fEventAction->SetReaction_theta_LAB(reaction_theta_LAB);
//    fEventAction->SetReaction_phi_LAB(reaction_phi_LAB);
//    fEventAction->SetReaction_theta_reactionCOM(reaction_theta_reactionCOM);
//    fEventAction->SetReaction_phi_reactionCOM(reaction_phi_reactionCOM);

    //========================================================================================
    //      8Be Decay (0, 1)
    G4ParticleDefinition* particleDefinition_alpha = G4ParticleTable::GetParticleTable()->FindParticle("alpha");
    fParticleGun->SetParticleDefinition(particleDefinition_alpha);
    
    int nDecayProducts_8Be = 2;
    std::vector<int>    decay8Be_Z{2, 2}; // e
    std::vector<double> decay8Be_A{4.00260325413, 4.00260325413}; // amu
    auto decay8Be_P = std::vector<double>(2, 0.0);
    auto decay8Be_T = std::vector<double>(2, 0.0);
    
    auto decay8Be_theta_LAB = std::vector<double>(2, 0.0);
    auto decay8Be_phi_LAB = std::vector<double>(2, 0.0);
    auto decay8Be_theta_recoilCOM = std::vector<double>(2, 0.0);
    auto decay8Be_phi_recoilCOM = std::vector<double>(2, 0.0);
    auto decay8Be_theta_reactionCOM = std::vector<double>(2, 0.0);
    auto decay8Be_phi_reactionCOM = std::vector<double>(2, 0.0);
    
    //----------------------------------------
    //      Polar angles (recoilCOM)
    decay8Be_theta_recoilCOM[0] = acos(1 - (2.0*G4UniformRand()))/deg;
    decay8Be_theta_recoilCOM[1] = 180.0 - decay8Be_theta_recoilCOM[0];
    
    //----------------------------------------
    //      Azimuthal angles (LAB)
    decay8Be_phi_recoilCOM[0] = 360.0*G4UniformRand(); // deg
    decay8Be_phi_recoilCOM[1] = decay8Be_phi_recoilCOM[0] - 180.0; // deg
    
    if(decay8Be_phi_recoilCOM[1]<0.0)
    {
        decay8Be_phi_recoilCOM[1] += 360.0;
    }
    */
    
//    auto decayVectors_8Be_4He_4He = CalculateBinaryDecayMomentumVectors(decayMassEnergies, reactionMomentumVectors[2], 0.0, 0.0);

    
//    reaction_P[2] = 0.0;
//    reaction_P[3] = 0.0;
    
//    reaction_T[2] = 0.0;
//    reaction_theta_LAB[2] = 0.0;
//    std::cout << "reaction_T[2]: " << reaction_T[2] << std::endl;
    
//    //----------------------------------------
//    CalculateBinaryDecayKinematics(recoilExcitationEnergy, daughterNucleusFinalStateEnergy, -0.09184, reaction_P[2], reaction_E[2], decay8Be_A, reaction_theta_LAB[2], reaction_phi_LAB[2], decay8Be_theta_recoilCOM, decay8Be_phi_recoilCOM, decay8Be_theta_LAB, decay8Be_phi_LAB, decay8Be_T);
//    
//    double mx, my, mz;
//    mx = sin(decay8Be_theta_LAB[0]*deg)*cos(decay8Be_phi_LAB[0]*deg);
//    my = sin(decay8Be_theta_LAB[0]*deg)*sin(decay8Be_phi_LAB[0]*deg);
//    mz = cos(decay8Be_theta_LAB[0]*deg);
//    
//    fEventAction->SetInitialParticleKineticEnergy(decay8Be_T[0]*MeV); // MeV
//    fParticleGun->SetParticleEnergy(decay8Be_T[0]*MeV);
//    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
//    std::cout << "decay8Be_T[0]: " << decay8Be_T[0] << std::endl;
//    std::cout << "decay8Be_theta_LAB[0]: " << decay8Be_theta_LAB[0] << std::endl;
//    
//    fParticleGun->GeneratePrimaryVertex(anEvent);
//    
//    //----------------------------------------
//    mx = sin(decay8Be_theta_LAB[1]*deg)*cos(decay8Be_phi_LAB[1]*deg);
//    my = sin(decay8Be_theta_LAB[1]*deg)*sin(decay8Be_phi_LAB[1]*deg);
//    mz = cos(decay8Be_theta_LAB[1]*deg);
//    
//    fEventAction->SetInitialParticleKineticEnergy(decay8Be_T[1]*MeV); // MeV
//    fParticleGun->SetParticleEnergy(decay8Be_T[1]*MeV);
//    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
//    std::cout << "decay8Be_T[1]: " << decay8Be_T[1] << std::endl;
//    std::cout << "decay8Be_theta_LAB[1]: " << decay8Be_theta_LAB[1] << std::endl;
//    
//    fParticleGun->GeneratePrimaryVertex(anEvent);
//    
//    fEventAction->SetNDecayProducts(nDecayProducts_8Be);
//    fEventAction->SetDecay_Z(decay8Be_Z);
//    fEventAction->SetDecay_A(decay8Be_A);
//    fEventAction->SetDecay_P(decay8Be_P);
//    fEventAction->SetDecay_T(decay8Be_T);
//    
//    fEventAction->SetDecay_theta_LAB(decay8Be_theta_LAB);
//    fEventAction->SetDecay_phi_LAB(decay8Be_phi_LAB);
//    fEventAction->SetDecay_theta_reactionCOM(decay8Be_theta_reactionCOM);
//    fEventAction->SetDecay_phi_reactionCOM(decay8Be_phi_reactionCOM);
//    fEventAction->SetDecay_theta_recoilCOM(decay8Be_theta_recoilCOM);
//    fEventAction->SetDecay_phi_recoilCOM(decay8Be_phi_recoilCOM);
//    
//    //========================================================================================
//    //              ELASTIC
//    G4ParticleDefinition* Lithium7 = G4IonTable::GetIonTable()->GetIon(3, 7.01600342665, 0.0);
//    fParticleGun->SetParticleDefinition(Lithium7);
//
//    double m_7Li_1H_Elastic[4];
//    
//    m_7Li_1H_Elastic[0] = 7.01600342665;
//    m_7Li_1H_Elastic[1] = 1.00782503223;
//    m_7Li_1H_Elastic[2] = 7.01600342665;
//    m_7Li_1H_Elastic[3] = 1.00782503223;
//    
//    //------------------------------------------------
//    double MassEnergy_7Li_1H_Elastic[4];
//    
//    for(int i=0; i<4; i++)
//    {
//        MassEnergy_7Li_1H_Elastic[i] = m_7Li_1H_Elastic[i]*931.49410242;
//    }
//    
//    double reaction_theta_reactionCOM_7Li_1H_Elastic = acos(1 - (2.0*G4UniformRand()))/deg;
//    double reaction_phi_reactionCOM_7Li_1H_Elastic = 360.0*G4UniformRand();
//    
//    auto reactionMomentumVectors_7Li_1H_Elastic = CalculateReactionMomentumVectors(MassEnergy_7Li_1H_Elastic, beamEnergy, reaction_theta_reactionCOM_7Li_1H_Elastic, reaction_phi_reactionCOM_7Li_1H_Elastic);
//    
////    double reaction_E_7Li_1H_Elastic = reactionMomentumVectors_7Li_1H_Elastic[2].e();
////    double reaction_T_7Li_1H_Elastic = reaction_E_7Li_1H_Elastic - MassEnergy_7Li_1H_Elastic[2];
//    
//    //----------------------------------------
//    double polarAngle_LAB_7Li_1H_Elastic = reactionMomentumVectors_7Li_1H_Elastic[2].theta()*180./M_PI;
//    double azimuthalAngle_LAB_7Li_1H_Elastic = reactionMomentumVectors_7Li_1H_Elastic[2].phi()*180./M_PI;
//    double energy_LAB_7Li_1H_Elastic = reactionMomentumVectors_7Li_1H_Elastic[2].e();
//    double kineticEnergy_LAB_7Li_1H_Elastic = energy_LAB_7Li_1H_Elastic - MassEnergy_7Li_1H_Elastic[2];
//    
////    std::cout << "polarAngle_LAB_7Li_1H_Elastic: " << polarAngle_LAB_7Li_1H_Elastic << std::endl;
////    std::cout << "energy_LAB_7Li_1H_Elastic: " << energy_LAB_7Li_1H_Elastic << std::endl;
////    std::cout << "kineticEnergy_LAB_7Li_1H_Elastic: " << kineticEnergy_LAB_7Li_1H_Elastic << std::endl;
//
//    mx = sin(polarAngle_LAB_7Li_1H_Elastic*deg)*cos(azimuthalAngle_LAB_7Li_1H_Elastic*deg);
//    my = sin(polarAngle_LAB_7Li_1H_Elastic*deg)*sin(azimuthalAngle_LAB_7Li_1H_Elastic*deg);
//    mz = cos(polarAngle_LAB_7Li_1H_Elastic*deg);
//
//    //----------------------------------------
//    fEventAction->SetInitialParticleKineticEnergy(kineticEnergy_LAB_7Li_1H_Elastic*MeV); // MeV
//    fParticleGun->SetParticleEnergy(kineticEnergy_LAB_7Li_1H_Elastic*MeV);
//    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
//    
//    fParticleGun->GeneratePrimaryVertex(anEvent);
//
//    //----------------------------------------
//    nReactionProducts++;
//    reaction_Z.push_back(3);
//    reaction_A.push_back(7.01600342665);
//    reaction_T.push_back(kineticEnergy_LAB_7Li_1H_Elastic);
//    
//    reaction_theta_LAB.push_back(polarAngle_LAB_7Li_1H_Elastic);
//    reaction_phi_LAB.push_back(azimuthalAngle_LAB_7Li_1H_Elastic);
//    reaction_theta_reactionCOM.push_back(reaction_theta_reactionCOM_7Li_1H_Elastic);
//    reaction_phi_reactionCOM.push_back(reaction_phi_reactionCOM_7Li_1H_Elastic);
//    
//    //----------------------------------------
//    fEventAction->SetNReactionProducts(nReactionProducts);
//    
//    fEventAction->SetReaction_Z(reaction_Z);
//    fEventAction->SetReaction_A(reaction_A);
//    fEventAction->SetReaction_T(reaction_T);
//    
//    fEventAction->SetReaction_theta_LAB(reaction_theta_LAB);
//    fEventAction->SetReaction_phi_LAB(reaction_phi_LAB);
//    fEventAction->SetReaction_theta_reactionCOM(reaction_theta_reactionCOM);
//    fEventAction->SetReaction_phi_reactionCOM(reaction_phi_reactionCOM);

    
    
    //        reaction_T[i] = reaction_E[i] - MassEnergy[i];
    //        reaction_P[i] = std::sqrt(std::pow(reaction_E[i], 2.0) - std::pow(MassEnergy[i], 2.0))*(1.0/std::pow(c2, 0.5));

//    for(int i=0; i<static_cast<int>(reactionMomentumVectorsreactionMomentumVectors_7Li_1H_Elastic.size()); i++)
//    {
//        reaction_E[i] = reactionMomentumVectors[i].e();
//        reaction_T[i] = reaction_E[i] - MassEnergy[i];
//        reaction_P[i] = std::sqrt(std::pow(reaction_E[i], 2.0) - std::pow(MassEnergy[i], 2.0))*(1.0/std::pow(c2, 0.5));
//    }

//    //========================================================================================
//    //      X17 decay (2, 3)
//    G4ParticleDefinition* particleDefinition_positron = G4ParticleTable::GetParticleTable()->FindParticle("e+");
//    G4ParticleDefinition* particleDefinition_electron = G4ParticleTable::GetParticleTable()->FindParticle("e-");
//    
//    int nDecayProducts_X17 = 2;
//    std::vector<int>    decayX17_Z{+1, -1}; // e
////    std::vector<double> decayX17_A{0.000548756, 0.000548756}; // amu
//    std::vector<double> decayX17_A{0.100548756, 0.100548756}; // amu
////    std::vector<double> decayX17_A{4.00260325413, 4.00260325413}; // amu
//    auto decayX17_P = std::vector<double>(2, 0.0);
//    auto decayX17_T = std::vector<double>(2, 0.0);
//    
//    auto decayX17_theta_LAB = std::vector<double>(2, 0.0);
//    auto decayX17_phi_LAB = std::vector<double>(2, 0.0);
//    auto decayX17_theta_recoilCOM = std::vector<double>(2, 0.0);
//    auto decayX17_phi_recoilCOM = std::vector<double>(2, 0.0);
//    auto decayX17_theta_reactionCOM = std::vector<double>(2, 0.0);
//    auto decayX17_phi_reactionCOM = std::vector<double>(2, 0.0);
//    
//    //----------------------------------------
//    //      Polar angles (recoilCOM)
//    decayX17_theta_recoilCOM[0] = acos(1 - (2.0*G4UniformRand()))/deg;
//    decayX17_theta_recoilCOM[1] = 180.0 - decayX17_theta_recoilCOM[0];
//    
//    //----------------------------------------
//    //      Azimuthal angles (LAB)
//    decayX17_phi_recoilCOM[0] = 360.0*G4UniformRand(); // deg
//    decayX17_phi_recoilCOM[1] = decayX17_phi_recoilCOM[0] - 180.0; // deg
//    
//    if(decayX17_phi_recoilCOM[1]<0.0)
//    {
//        decayX17_phi_recoilCOM[1] += 360.0;
//    }
//    
//    //    std::cout << "reaction_P[0]: " << reaction_P[0] << std::endl;
//    //    std::cout << "reaction_P[1]: " << reaction_P[1] << std::endl;
//    //    std::cout << "reaction_P[2]: " << reaction_P[2] << std::endl;
//    //    std::cout << "reaction_P[3]: " << reaction_P[3] << std::endl;
//    //
//    //    std::cout << "reaction_T[0]: " << reaction_T[0] << std::endl;
//    //    std::cout << "reaction_T[1]: " << reaction_T[1] << std::endl;
//    //    std::cout << "reaction_T[2]: " << reaction_T[2] << std::endl;
//    //    std::cout << "reaction_T[3]: " << reaction_T[3] << std::endl;
//    
//    std::cout << std::endl;
//    std::cout << "reaction_T[3]: " << reaction_T[3] << std::endl;
//    std::cout << "reaction_P[3]: " << reaction_P[3] << std::endl;
//    std::cout << "reaction_theta_LAB[3]: " << reaction_theta_LAB[3] << std::endl;
//    std::cout << "decayX17_theta_recoilCOM[0]: " << decayX17_theta_recoilCOM[0] << std::endl;
//    std::cout << "decayX17_theta_recoilCOM[1]: " << decayX17_theta_recoilCOM[1] << std::endl;
//    std::cout << "decayX17_phi_recoilCOM[0]: " << decayX17_phi_recoilCOM[0] << std::endl;
//    std::cout << "decayX17_phi_recoilCOM[1]: " << decayX17_phi_recoilCOM[1] << std::endl;
//
//    reaction_P[3] = 0.0;
//    reaction_E[3] = 0.0;
//    //----------------------------------------
////    CalculateBinaryDecayKinematics(16.7, 0.0, (2.0*0.00054858*931.49410242), reaction_P[3], reaction_E[3], decayX17_A, reaction_theta_LAB[3], reaction_phi_LAB[3], decayX17_theta_recoilCOM, decayX17_phi_recoilCOM, decayX17_theta_LAB, decayX17_phi_LAB, decayX17_T);
//    CalculateBinaryDecayKinematics(16.7, 0.0, (0.0), reaction_P[3], reaction_E[3], decayX17_A, reaction_theta_LAB[3], reaction_phi_LAB[3], decayX17_theta_recoilCOM, decayX17_phi_recoilCOM, decayX17_theta_LAB, decayX17_phi_LAB, decayX17_T);
//    
//    mx = sin(decayX17_theta_LAB[0]*deg)*cos(decayX17_phi_LAB[0]*deg);
//    my = sin(decayX17_theta_LAB[0]*deg)*sin(decayX17_phi_LAB[0]*deg);
//    mz = cos(decayX17_theta_LAB[0]*deg);
//    
//    fEventAction->SetInitialParticleKineticEnergy(decayX17_T[0]*MeV); // MeV
//    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
//    std::cout << "decayX17_T[0]: " << decayX17_T[0] << std::endl;
//    std::cout << "decayX17_theta_LAB[0]: " << decayX17_theta_LAB[0] << std::endl;
//    
//    fParticleGun->SetParticleDefinition(particleDefinition_positron);
//    fParticleGun->GeneratePrimaryVertex(anEvent);
//    
//    //----------------------------------------
//    mx = sin(decayX17_theta_LAB[1]*deg)*cos(decayX17_phi_LAB[1]*deg);
//    my = sin(decayX17_theta_LAB[1]*deg)*sin(decayX17_phi_LAB[1]*deg);
//    mz = cos(decayX17_theta_LAB[1]*deg);
//    
//    fEventAction->SetInitialParticleKineticEnergy(decayX17_T[1]*MeV); // MeV
//    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
//    std::cout << "decayX17_T[1]: " << decayX17_T[1] << std::endl;
//    std::cout << "decayX17_theta_LAB[1]: " << decayX17_theta_LAB[1] << std::endl;
//    std::cout << std::endl;
//    
//    fParticleGun->SetParticleDefinition(particleDefinition_electron);
//    fParticleGun->GeneratePrimaryVertex(anEvent);

    
//    void PrimaryGeneratorAction::CalculateBinaryDecayKinematics(double eX, double eX_residual, double sepE, double initialScatteringReactionMomentum, double initialScatteringReactionEnergy, std::vector<double> decay_A, double theta_recoil_LAB, double phi_recoil_LAB, std::vector<double>& decay_theta_recoilCOM, std::vector<double>& decay_phi_recoilCOM, std::vector<double>& decay_theta_LAB, std::vector<double>& decay_phi_LAB, std::vector<double>& decay_T_LAB)

    
    /*
    //========================================================================================
    //      TEST with CalculateBinaryDecayKinematics2
    //      X17 decay (2, 3)
    G4ParticleDefinition* particleDefinition_positron = G4ParticleTable::GetParticleTable()->FindParticle("e+");
    G4ParticleDefinition* particleDefinition_electron = G4ParticleTable::GetParticleTable()->FindParticle("e-");
    
    int nDecayProducts_X17 = 2;
    std::vector<int>    decayX17_Z{+1, -1}; // e
        std::vector<long double> decayX17_A{0.000548756, 0.000548756}; // amu
//    std::vector<long double> decayX17_A{0.100548756, 0.100548756}; // amu
    //    std::vector<long double> decayX17_A{4.00260325413, 4.00260325413}; // amu
    auto decayX17_P = std::vector<long double>(2, 0.0);
    auto decayX17_T = std::vector<long double>(2, 0.0);
    
    auto decayX17_theta_LAB = std::vector<long double>(2, 0.0);
    auto decayX17_phi_LAB = std::vector<long double>(2, 0.0);
    auto decayX17_theta_recoilCOM = std::vector<long double>(2, 0.0);
    auto decayX17_phi_recoilCOM = std::vector<long double>(2, 0.0);
    auto decayX17_theta_reactionCOM = std::vector<long double>(2, 0.0);
    auto decayX17_phi_reactionCOM = std::vector<long double>(2, 0.0);
    
    //----------------------------------------
    //      Polar angles (recoilCOM)
    decayX17_theta_recoilCOM[0] = acos(1 - (2.0*G4UniformRand()))/deg;
    decayX17_theta_recoilCOM[1] = 180.0 - decayX17_theta_recoilCOM[0];
    
    //----------------------------------------
    //      Azimuthal angles (LAB)
    decayX17_phi_recoilCOM[0] = 360.0*G4UniformRand(); // deg
    decayX17_phi_recoilCOM[1] = decayX17_phi_recoilCOM[0] - 180.0; // deg
    
    if(decayX17_phi_recoilCOM[1]<0.0)
    {
        decayX17_phi_recoilCOM[1] += 360.0;
    }
    
    //    std::cout << "reaction_P[0]: " << reaction_P[0] << std::endl;
    //    std::cout << "reaction_P[1]: " << reaction_P[1] << std::endl;
    //    std::cout << "reaction_P[2]: " << reaction_P[2] << std::endl;
    //    std::cout << "reaction_P[3]: " << reaction_P[3] << std::endl;
    //
    //    std::cout << "reaction_T[0]: " << reaction_T[0] << std::endl;
    //    std::cout << "reaction_T[1]: " << reaction_T[1] << std::endl;
    //    std::cout << "reaction_T[2]: " << reaction_T[2] << std::endl;
    //    std::cout << "reaction_T[3]: " << reaction_T[3] << std::endl;
    
    std::cout << std::endl;
    std::cout << "reaction_T[3]: " << reaction_T[3] << std::endl;
    std::cout << "reaction_P[3]: " << reaction_P[3] << std::endl;
    std::cout << "reaction_theta_LAB[3]: " << reaction_theta_LAB[3] << std::endl;
    std::cout << "decayX17_theta_recoilCOM[0]: " << decayX17_theta_recoilCOM[0] << std::endl;
    std::cout << "decayX17_theta_recoilCOM[1]: " << decayX17_theta_recoilCOM[1] << std::endl;
    std::cout << "decayX17_phi_recoilCOM[0]: " << decayX17_phi_recoilCOM[0] << std::endl;
    std::cout << "decayX17_phi_recoilCOM[1]: " << decayX17_phi_recoilCOM[1] << std::endl;
    
    reaction_P[3] = 0.0;
    reaction_E[3] = 0.0;
    //----------------------------------------
    //    CalculateBinaryDecayKinematics(16.7, 0.0, (2.0*0.00054858*931.49410242), reaction_P[3], reaction_E[3], decayX17_A, reaction_theta_LAB[3], reaction_phi_LAB[3], decayX17_theta_recoilCOM, decayX17_phi_recoilCOM, decayX17_theta_LAB, decayX17_phi_LAB, decayX17_T);
    CalculateBinaryDecayKinematics2(16.7, 0.0, (2.0*0.00054858*931.49410242), (long double) reaction_P[3], (long double) reaction_E[3], decayX17_A, (long double) reaction_theta_LAB[3], (long double) reaction_phi_LAB[3], decayX17_theta_recoilCOM, decayX17_phi_recoilCOM, decayX17_theta_LAB, decayX17_phi_LAB, decayX17_T);
//    CalculateBinaryDecayKinematics3(16.7, 0.0, (2.0*0.00054858*931.49410242), (long double) reaction_P[3], (long double) reaction_E[3], decayX17_A, (long double) reaction_theta_LAB[3], (long double) reaction_phi_LAB[3], decayX17_theta_recoilCOM, decayX17_phi_recoilCOM, decayX17_theta_LAB, decayX17_phi_LAB, decayX17_T);
    
    mx = sin(decayX17_theta_LAB[0]*deg)*cos(decayX17_phi_LAB[0]*deg);
    my = sin(decayX17_theta_LAB[0]*deg)*sin(decayX17_phi_LAB[0]*deg);
    mz = cos(decayX17_theta_LAB[0]*deg);
    
    fEventAction->SetInitialParticleKineticEnergy(decayX17_T[0]*MeV); // MeV
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
    std::cout << "decayX17_T[0]: " << decayX17_T[0] << std::endl;
    std::cout << "decayX17_theta_LAB[0]: " << decayX17_theta_LAB[0] << std::endl;
    
    fParticleGun->SetParticleDefinition(particleDefinition_positron);
    fParticleGun->GeneratePrimaryVertex(anEvent);
    
    //----------------------------------------
    mx = sin(decayX17_theta_LAB[1]*deg)*cos(decayX17_phi_LAB[1]*deg);
    my = sin(decayX17_theta_LAB[1]*deg)*sin(decayX17_phi_LAB[1]*deg);
    mz = cos(decayX17_theta_LAB[1]*deg);
    
    fEventAction->SetInitialParticleKineticEnergy(decayX17_T[1]*MeV); // MeV
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
    std::cout << "decayX17_T[1]: " << decayX17_T[1] << std::endl;
    std::cout << "decayX17_theta_LAB[1]: " << decayX17_theta_LAB[1] << std::endl;
    std::cout << std::endl;
    
    fParticleGun->SetParticleDefinition(particleDefinition_electron);
    fParticleGun->GeneratePrimaryVertex(anEvent);
*/
    
/*
    //========================================================================================
    //      X17 decay (2, 3)
    G4ParticleDefinition* particleDefinition_positron = G4ParticleTable::GetParticleTable()->FindParticle("e+");
    G4ParticleDefinition* particleDefinition_electron = G4ParticleTable::GetParticleTable()->FindParticle("e-");

    int nDecayProducts_X17 = 2;
    std::vector<int>    decayX17_Z{+1, -1}; // e
    std::vector<double> decayX17_A{0.000548756, 0.000548756}; // amu
    auto decayX17_P = std::vector<double>(2, 0.0);
    auto decayX17_T = std::vector<double>(2, 0.0);
    
    auto decayX17_theta_LAB = std::vector<double>(2, 0.0);
    auto decayX17_phi_LAB = std::vector<double>(2, 0.0);
    auto decayX17_theta_recoilCOM = std::vector<double>(2, 0.0);
    auto decayX17_phi_recoilCOM = std::vector<double>(2, 0.0);
    auto decayX17_theta_reactionCOM = std::vector<double>(2, 0.0);
    auto decayX17_phi_reactionCOM = std::vector<double>(2, 0.0);
    
    //----------------------------------------
    //      Polar angles (recoilCOM)
    decayX17_theta_recoilCOM[0] = acos(1 - (2.0*G4UniformRand()))/deg;
    decayX17_theta_recoilCOM[1] = 180.0 - decayX17_theta_recoilCOM[0];
    
    //----------------------------------------
    //      Azimuthal angles (LAB)
    decayX17_phi_recoilCOM[0] = 360.0*G4UniformRand(); // deg
    decayX17_phi_recoilCOM[1] = decayX17_phi_recoilCOM[0] - 180.0; // deg
    
    if(decayX17_phi_recoilCOM[1]<0.0)
    {
        decayX17_phi_recoilCOM[1] += 360.0;
    }
    
//    G4cout << "reaction_P[0]: " << reaction_P[0] << G4endl;
//    G4cout << "reaction_P[1]: " << reaction_P[1] << G4endl;
//    G4cout << "reaction_P[2]: " << reaction_P[2] << G4endl;
//    G4cout << "reaction_P[3]: " << reaction_P[3] << G4endl;
//
//    G4cout << "reaction_T[0]: " << reaction_T[0] << G4endl;
//    G4cout << "reaction_T[1]: " << reaction_T[1] << G4endl;
//    G4cout << "reaction_T[2]: " << reaction_T[2] << G4endl;
//    G4cout << "reaction_T[3]: " << reaction_T[3] << G4endl;

    //----------------------------------------
    CalculateBinaryDecayKinematics(0.0, 0.0, -(17.0 - 2.0*0.000548756*931.49410242), reaction_P[2], reaction_E[2], decayX17_A, reaction_theta_LAB[2], reaction_phi_LAB[2], decayX17_theta_recoilCOM, decayX17_phi_recoilCOM, decayX17_theta_LAB, decayX17_phi_LAB, decayX17_T);
    
    mx = sin(decayX17_theta_LAB[0]*deg)*cos(decayX17_phi_LAB[0]*deg);
    my = sin(decayX17_theta_LAB[0]*deg)*sin(decayX17_phi_LAB[0]*deg);
    mz = cos(decayX17_theta_LAB[0]*deg);
    
    fEventAction->SetInitialParticleKineticEnergy(decayX17_T[0]*MeV); // MeV
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
//    G4cout << "decayX17_T[0]: " << decayX17_T[0] << G4endl;
//    G4cout << "decayX17_theta_LAB[0]: " << decayX17_theta_LAB[0] << G4endl;
    
    fParticleGun->SetParticleDefinition(particleDefinition_positron);
    fParticleGun->GeneratePrimaryVertex(anEvent);
    
    //----------------------------------------
    mx = sin(decayX17_theta_LAB[1]*deg)*cos(decayX17_phi_LAB[1]*deg);
    my = sin(decayX17_theta_LAB[1]*deg)*sin(decayX17_phi_LAB[1]*deg);
    mz = cos(decayX17_theta_LAB[1]*deg);
    
    fEventAction->SetInitialParticleKineticEnergy(decayX17_T[1]*MeV); // MeV
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
//    G4cout << "decayX17_T[1]: " << decayX17_T[1] << G4endl;
//    G4cout << "decayX17_theta_LAB[1]: " << decayX17_theta_LAB[1] << G4endl;
//    G4cout << G4endl;
    
    fParticleGun->SetParticleDefinition(particleDefinition_electron);
    fParticleGun->GeneratePrimaryVertex(anEvent);

    
    //------------------------------------------------------------------------------------
    int nDecayProducts = nDecayProducts_8Be + nDecayProducts_X17;
    std::vector<int>    decay_Z{decay8Be_Z[0], decay8Be_Z[1], decayX17_Z[0], decayX17_Z[1]}; // e
    std::vector<double> decay_A{decay8Be_A[0], decay8Be_A[1], decayX17_A[0], decayX17_A[1]}; // amu
    
    std::vector<double> decay_P{decay8Be_P[0], decay8Be_P[1], decayX17_P[0], decayX17_P[1]};
    std::vector<double> decay_T{decay8Be_T[0], decay8Be_T[1], decayX17_T[0], decayX17_T[1]};
    
    std::vector<double> decay_theta_LAB{decay8Be_theta_LAB[0], decay8Be_theta_LAB[1], decayX17_theta_LAB[0], decayX17_theta_LAB[1]};
    std::vector<double> decay_phi_LAB{decay8Be_phi_LAB[0], decay8Be_phi_LAB[1], decayX17_phi_LAB[0], decayX17_phi_LAB[1]};
    std::vector<double> decay_theta_recoilCOM{decay8Be_theta_recoilCOM[0], decay8Be_theta_recoilCOM[1], decayX17_theta_recoilCOM[0], decayX17_theta_recoilCOM[1]};
    std::vector<double> decay_phi_recoilCOM{decay8Be_phi_recoilCOM[0], decay8Be_phi_recoilCOM[1], decayX17_phi_recoilCOM[0], decayX17_phi_recoilCOM[1]};
    std::vector<double> decay_theta_reactionCOM{decay8Be_theta_reactionCOM[0], decay8Be_theta_reactionCOM[1], decayX17_theta_reactionCOM[0], decayX17_theta_reactionCOM[1]};
    std::vector<double> decay_phi_reactionCOM{decay8Be_phi_reactionCOM[0], decay8Be_phi_reactionCOM[1], decayX17_phi_reactionCOM[0], decayX17_phi_reactionCOM[1]};

    //----------------------------------------
//    fEventAction->SetNDecayProducts(nDecayProducts);
//    fEventAction->SetDecay_Z(decay_Z);
//    fEventAction->SetDecay_A(decay_A);
//    fEventAction->SetDecay_P(decay_P);
//    fEventAction->SetDecay_T(decay_T);
//
//    fEventAction->SetDecay_theta_LAB(decay_theta_LAB);
//    fEventAction->SetDecay_phi_LAB(decay_phi_LAB);
//    fEventAction->SetDecay_theta_reactionCOM(decay_theta_reactionCOM);
//    fEventAction->SetDecay_phi_reactionCOM(decay_phi_reactionCOM);
//    fEventAction->SetDecay_theta_recoilCOM(decay_theta_recoilCOM);
//    fEventAction->SetDecay_phi_recoilCOM(decay_phi_recoilCOM);
*/
    
    
    
    
//    //------------------------------------------------------------------------------------
//    int nDecayProducts_8Be = 2;
//    std::vector<int>    decay_Z_8Be{2, 2}; // e
//    std::vector<double> decay_A_8Be{4.00260325413, 4.00260325413}; // amu
//    auto decay_P_8Be = std::vector<double>(2, 0.0);
//    auto decay_T_8Be = std::vector<double>(2, 0.0);
//    
//    auto decay_theta_LAB_8Be = std::vector<double>(2, 0.0);
//    auto decay_phi_LAB_8Be = std::vector<double>(2, 0.0);
//    auto decay_theta_recoilCOM_8Be = std::vector<double>(2, 0.0);
//    auto decay_phi_recoilCOM_8Be = std::vector<double>(2, 0.0);
//    auto decay_theta_reactionCOM_8Be = std::vector<double>(2, 0.0);
//    auto decay_phi_reactionCOM_8Be = std::vector<double>(2, 0.0);
//
//    //------------------------------------------------------------------------------------
//    //      Polar angles (recoilCOM)
//    decay_theta_recoilCOM_8Be[0] = acos(1 - (2.0*G4UniformRand()))/deg;
//    decay_theta_recoilCOM_8Be[1] = 180.0 - decay_theta_recoilCOM_8Be[0];
//    
//    //std::cout << "decay_theta_recoilCOM_8Be[0]: " << decay_theta_recoilCOM_8Be[0] << std::endl;
//    //std::cout << "decay_theta_recoilCOM_8Be[1]: " << decay_theta_recoilCOM_8Be[1] << std::endl;
//    
//    //------------------------------------------------------------------------------------
//    //      Azimuthal angles (LAB)
//    decay_phi_recoilCOM_8Be[0] = 360.0*G4UniformRand(); // deg
//    decay_phi_recoilCOM_8Be[1] = decay_phi_recoilCOM_8Be[0] - 180.0; // deg
//    
//    if(decay_phi_recoilCOM_8Be[1]<0.0)
//    {
//        decay_phi_recoilCOM_8Be[1] += 360.0;
//    }
//    
//    //------------------------------------------------------------------------------------
//    //reaction_P = std::vector<double>(4, 0.0);
//
//    //G4cout << "decay_T[0]: " << decay_T[0] << G4endl;
//    //G4cout << "decay_T[1]: " << decay_T[1] << G4endl;
//    
//    double c2 = 931.49410242;     // MeV/u, c^2
//    double P_8Be = sqrt(2.0*decay_A[1]*decay_T[1]);
//    double E_8Be = sqrt(pow(P_8Be, 2.0)*c2 + pow(decay_A[1]*c2, 2.0));
//    
//    //double P_8Be = 0.0;
//    //double E_8Be = sqrt(pow(P_8Be*931.5, 2.0) + pow(decay_A[1]*931.5, 2.0));
//    
//    //std::cout << "P_8Be: " << P_8Be << std::endl;
//    //std::cout << "reaction_theta_LAB[3]: " << reaction_theta_LAB[3] << std::endl;
//    
//    //void PrimaryGeneratorAction::CalculateBinaryDecayKinematics(double eX, double eX_residual, double sepE, double initialScatteringReactionMomentum, double initialScatteringReactionEnergy, std::vector<double> decay_A, double theta_recoil_LAB, double phi_recoil_LAB, std::vector<double>& decay_theta_recoilCOM, std::vector<double>& decay_phi_recoilCOM, std::vector<double>& decay_theta_LAB, std::vector<double>& decay_phi_LAB, std::vector<double>& decay_T)
//    
//    CalculateBinaryDecayKinematics(0.0, 0.0, -0.09184, P_8Be, E_8Be, decay_A_8Be, reaction_theta_LAB[3], reaction_phi_LAB[3], decay_theta_recoilCOM_8Be, decay_phi_recoilCOM_8Be, decay_theta_LAB_8Be, decay_phi_LAB_8Be, decay_T_8Be);
//    
//    //decay_theta_LAB_8Be[0] = decay_theta_recoilCOM_8Be[0];
//    //decay_theta_LAB_8Be[1] = decay_theta_recoilCOM_8Be[1];
//    
//    //decay_phi_LAB_8Be[0] = decay_phi_recoilCOM_8Be[0];
//    //decay_phi_LAB_8Be[1] = decay_phi_recoilCOM_8Be[1];
//    
//    //------------------------------------------------------------------------------------
//    double mx, my, mz;
//    
//    mx = sin(decay_theta_LAB_8Be[0]*deg)*cos(decay_phi_LAB_8Be[0]*deg);
//    my = sin(decay_theta_LAB_8Be[0]*deg)*sin(decay_phi_LAB_8Be[0]*deg);
//    mz = cos(decay_theta_LAB_8Be[0]*deg);
//    
//    fEventAction->SetInitialParticleKineticEnergy(decay_T_8Be[0]*MeV); // MeV
//    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
//    //G4cout << "decay_T_8Be[0]: " << decay_T_8Be[0] << G4endl;
//    //G4cout << "decay_theta_LAB_8Be[0]: " << decay_theta_LAB_8Be[0] << G4endl;
//    
//    fParticleGun->GeneratePrimaryVertex(anEvent);
//    
//    //------------------------------------------------------------------------------------
//    mx = sin(decay_theta_LAB_8Be[1]*deg)*cos(decay_phi_LAB_8Be[1]*deg);
//    my = sin(decay_theta_LAB_8Be[1]*deg)*sin(decay_phi_LAB_8Be[1]*deg);
//    mz = cos(decay_theta_LAB_8Be[1]*deg);
//    
//    fEventAction->SetInitialParticleKineticEnergy(decay_T_8Be[1]*MeV); // MeV
//    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
//    //G4cout << "decay_T_8Be[1]: " << decay_T_8Be[1] << G4endl;
//    //G4cout << "decay_theta_LAB_8Be[1]: " << decay_theta_LAB_8Be[1] << G4endl;
//    
//    fParticleGun->GeneratePrimaryVertex(anEvent);

    
    
    
    
    /*
    //------------------------------------------------
    double reactionPosition = 1.1197459*G4UniformRand() - 1.1197459/2.0;
    fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,reactionPosition*um));

    for(G4int i=0; i<4; i++)
    {
        T[i] = 0.0;
        E[i] = 0.0;
        p[i] = 0.0;
    }
    
    //------------------------------------------------
    ////	X17
    m[0] = 7.01600342665;
    m[1] = 1.00782503223;
    m[2] = 16.7/931.49410242;
    m[3] = 8.005305102;
    
    ////    Projectile/Beam Kinetic Energy (Lab Frame)
    T[0] = 3.07003;
    
    ////    Target Kinetic Energy (Lab Frame)
    T[1] = 0.0;
    
    ////    Excitation Energy
    eX = 0.0;
    
    //------------------------------------------------
    ////    Executing BiRelKin
    theta_ejectile = G4UniformRand()*180.0;
    BiRelKin(m, T, E, p, theta_ejectile, theta_recoil_LAB, eX);
    
    //------------------------------------------------------------------------------------
    //      Polar angles (LAB)
    std::vector<double> reaction_theta_LAB{0.0, 0.0, theta_ejectile, 0.0}; // deg
    auto reaction_theta_reactionCOM = std::vector<double>(4, 0.0);
    
    //------------------------------------------------------------------------------------
    //      Azimuthal angles (LAB)
    double reaction_ejectile_phi_LAB = 360.0*G4UniformRand(); // deg
    double reaction_recoil_phi_LAB = reaction_ejectile_phi_LAB - 180.0; // deg
    if(reaction_recoil_phi_LAB<0.0)
    {
        reaction_recoil_phi_LAB += 360.0;
    }
    
    std::vector<double> reaction_phi_LAB{0.0, 0.0, reaction_ejectile_phi_LAB, reaction_recoil_phi_LAB}; // deg
    auto reaction_phi_reactionCOM = std::vector<double>(4, 0.0);

    //------------------------------------------------------------------------------------
    int nDecayProducts = 2;
    std::vector<int>    decay_Z{2, 4}; // e
    std::vector<double> decay_A{4.00260325413, 8.005305102}; // amu
    auto decay_P = std::vector<double>(2, 0.0);
    auto decay_T = std::vector<double>(2, 0.0);
    
    auto decay_theta_LAB = std::vector<double>(2, 0.0);
    auto decay_phi_LAB = std::vector<double>(2, 0.0);
    auto decay_theta_recoilCOM = std::vector<double>(2, 0.0);
    auto decay_phi_recoilCOM = std::vector<double>(2, 0.0);
    auto decay_theta_reactionCOM = std::vector<double>(2, 0.0);
    auto decay_phi_reactionCOM = std::vector<double>(2, 0.0);
    
    //------------------------------------------------------------------------------------
    //      Polar angles (recoilCOM)
    decay_theta_recoilCOM[0] = acos(1 - (2.0*G4UniformRand()))/deg;
    decay_theta_recoilCOM[1] = 180.0 - decay_theta_recoilCOM[0];
    
    //------------------------------------------------------------------------------------
    //      Azimuthal angles (LAB)
    decay_phi_recoilCOM[0] = 360.0*G4UniformRand(); // deg
    decay_phi_recoilCOM[1] = decay_phi_recoilCOM[0] - 180.0; // deg
    
    if(decay_phi_recoilCOM[1]<0.0)
    {
        decay_phi_recoilCOM[1] += 360.0;
    }
    
    //------------------------------------------------------------------------------------
    reaction_P = std::vector<double>(4, 0.0);
    CalculateBinaryDecayKinematics(0.0, 0.0, -0.09184, p[3], E[3], decay_A, reaction_theta_LAB[3], reaction_phi_LAB[3], decay_theta_recoilCOM, decay_phi_recoilCOM, decay_theta_LAB, decay_phi_LAB, decay_T);
    
//    CalculateBinaryDecayKinematics(recoilExcitationEnergy, daughterNucleusFinalStateEnergy, 7.36659, reaction_P[3], reaction_E[3], decay_A, reaction_theta_LAB[3], reaction_phi_LAB[3], decay_theta_recoilCOM, decay_phi_recoilCOM, decay_theta_LAB, decay_phi_LAB, decay_T);

//    void PrimaryGeneratorAction::CalculateBinaryDecayKinematics(double eX, double eX_residual, double sepE, double initialScatteringReactionMomentum, double initialScatteringReactionEnergy, std::vector<double> decay_A, double theta_recoil_LAB, double phi_recoil_LAB, std::vector<double>& decay_theta_recoilCOM, std::vector<double>& decay_phi_recoilCOM, std::vector<double>& decay_theta_LAB, std::vector<double>& decay_phi_LAB, std::vector<double>& decay_T_LAB)

    //------------------------------------------------------------------------------------
    fEventAction->SetNDecayProducts(nDecayProducts);
    fEventAction->SetDecay_Z(decay_Z);
    fEventAction->SetDecay_A(decay_A);
    fEventAction->SetDecay_P(decay_P);
    fEventAction->SetDecay_T(decay_T);

//    fEventAction->SetDecay_theta_reactionCOM(decay_theta_reactionCOM);
//    fEventAction->SetDecay_phi_reactionCOM(decay_phi_reactionCOM);

    fEventAction->SetDecay_theta_recoilCOM(decay_theta_recoilCOM);
    fEventAction->SetDecay_phi_recoilCOM(decay_phi_recoilCOM);

    fEventAction->SetDecay_theta_LAB(decay_theta_LAB);
    fEventAction->SetDecay_phi_LAB(decay_phi_LAB);
    
    std::cout << "decay_theta_recoilCOM[0]: " << decay_theta_recoilCOM[0] << std::endl;
    std::cout << "decay_theta_recoilCOM[1]: " << decay_theta_recoilCOM[1] << std::endl;
    
    std::cout << "decay_theta_LAB[0]: " << decay_theta_LAB[0] << std::endl;
    std::cout << "decay_theta_LAB[1]: " << decay_theta_LAB[1] << std::endl;
    
    //------------------------------------------------------------------------------------
    double mx, my, mz;
    
    mx = sin(decay_theta_LAB[0]*deg)*cos(decay_phi_LAB[0]*deg);
    my = sin(decay_theta_LAB[0]*deg)*sin(decay_phi_LAB[0]*deg);
    mz = cos(decay_theta_LAB[0]*deg);
    
    fEventAction->SetInitialParticleKineticEnergy(decay_T[0]); // MeV
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));

    fParticleGun->GeneratePrimaryVertex(anEvent);

    //------------------------------------------------------------------------------------
    mx = sin(decay_theta_LAB[1]*deg)*cos(decay_phi_LAB[1]*deg);
    my = sin(decay_theta_LAB[1]*deg)*sin(decay_phi_LAB[1]*deg);
    mz = cos(decay_theta_LAB[1]*deg);
    
    fEventAction->SetInitialParticleKineticEnergy(decay_T[1]); // MeV
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
    
    fParticleGun->GeneratePrimaryVertex(anEvent);
    */
    
    
    
    
    
    /*
    std::vector<double> decay_theta_LAB;
    
    G4double theta = acos(1 - (2.0*G4UniformRand()))/deg; // 0.0->180.0
    //G4double theta = acos(1 - (1.0*G4UniformRand()))/deg; // 0.0->90.0 deg
    //G4double theta = acos(1 - (1.0*G4UniformRand() + 1.0))/deg; // 90.0->180.0 deg
    
    decay_theta_LAB.push_back(theta);
    
    ////    Collimator - 2.0 deg range
    //G4double theta = acos(1 - (6.09173503004822869e-04*G4UniformRand()))/deg; // deg
    
    G4double phi = 360.0*G4UniformRand(); // deg
    
    //------------------------------------------------
    mx = sin(theta*deg)*cos(phi*deg);
    my = sin(theta*deg)*sin(phi*deg);
    mz = cos(theta*deg);
    
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
    
    fEventAction->SetDecay_theta_LAB(decay_theta_LAB);
    
    fParticleGun->GeneratePrimaryVertex(anEvent);

    */
    
    
    //========================================================================================================================
    //========================================================================================================================
    //      PR240: alpha0 decay, 0+ -> 0+, alpha decay from 12C
    
    int nReactionProducts = 4;
    std::vector<int>    reaction_Z{1, 6, 1, 6}; // e
    std::vector<double> reaction_A{1.00727647, 14.00324199, 3.01550071, 12}; // amu
    auto reaction_P = std::vector<double>(4, 0.0);
    std::vector<double> reaction_T{100.0, 0.0, 0.0, 0.0}; // MeV
    std::vector<double> reaction_Ex{0.0, 0.0, 0.0, 0.0}; // MeV
    auto reaction_E = std::vector<double>(4, 0.0);
    
    //------------------------------------------------------------------------------------
    //      Polar angles (LAB)
    std::vector<double> reaction_theta_LAB{0.0, 0.0, acos(1 - (6.09173503004822869e-04*G4UniformRand()))/deg, 0.0}; // deg
    //auto reaction_theta_reactionCOM = std::vector<double>(4, 0.0);
    
    //------------------------------------------------------------------------------------
    //      Azimuthal angles (LAB)
    double reaction_ejectile_phi_LAB = 360.0*G4UniformRand(); // deg
    double reaction_recoil_phi_LAB = reaction_ejectile_phi_LAB - 180.0; // deg
    if(reaction_recoil_phi_LAB<0.0)
    {
        reaction_recoil_phi_LAB += 360.0;
    }
    
    std::vector<double> reaction_phi_LAB{0.0, 0.0, reaction_ejectile_phi_LAB, reaction_recoil_phi_LAB}; // deg
    //auto reaction_phi_reactionCOM = std::vector<double>(4, 0.0);
    
    //------------------------------------------------------------------------------------
    //      Recoil Nucleux Excitation Energy
    //double recoilExcitationEnergy = 0.0; // MeV
    //double recoilExcitationEnergy = 7.654; // MeV
    //double recoilExcitationEnergy = 7.36659; // MeV
    //double recoilExcitationEnergy = 7.36659 + 0.5; // MeV
    //double recoilExcitationEnergy = 12.0; // MeV
    double recoilExcitationEnergy = 20.0*G4UniformRand(); // MeV
//    double recoilExcitationEnergy = 9.0; // MeV
//    double recoilExcitationEnergy = G4RandFlat::shoot(7.0, 20.0); // MeV

    
    //------------------------------------------------------------------------------------
    
    //G4double recoilExcitationEnergy = 0.0*MeV;
    
    //int energyN = (int) GetParticleN()/nParticlesPerEnergy;
    
    //if(energyN<((int) initialKineticEnergies.size()))
    //{
    //    recoilExcitationEnergy = initialKineticEnergies[energyN];
    //}
    //else
    //{
    //    recoilExcitationEnergy = 30.0*MeV;
    //}
    
    
    //------------------------------------------------------------------------------------
    //      Executing the Binary Relativistic Kinematics code
    BiRelKin(&reaction_A[0], &reaction_T[0], &reaction_E[0], &reaction_P[0], reaction_theta_LAB[2], reaction_theta_LAB[3], recoilExcitationEnergy);
    
    //G4cout << "reaction_T[3]: " << reaction_T[3] << G4endl;
    
    reaction_Ex[3] = recoilExcitationEnergy;
    
    //------------------------------------------------------------------------------------
    //      Filling the vectors in the EventAction object
    fEventAction->SetNReactionProducts(nReactionProducts);
    fEventAction->SetReaction_Z(reaction_Z);
    fEventAction->SetReaction_A(reaction_A);
    fEventAction->SetReaction_P(reaction_P);
    fEventAction->SetReaction_T(reaction_T);
    fEventAction->SetReaction_Ex(reaction_Ex);
    
    fEventAction->SetReaction_theta_LAB(reaction_theta_LAB);
    fEventAction->SetReaction_phi_LAB(reaction_phi_LAB);
    //fEventAction->SetReaction_theta_reactionCOM(reaction_theta_reactionCOM);
    //fEventAction->SetReaction_phi_reactionCOM(reaction_phi_reactionCOM);
    
    
    //------------------------------------------------------------------------------------
    //      Recoil Nucleus Decay
    
    int nDecayProducts = 2;
    std::vector<int>    decay_Z{2, 4}; // e
    std::vector<double> decay_A{4.00260325413, 8.005305102}; // amu
    auto decay_P = std::vector<double>(2, 0.0);
    auto decay_T = std::vector<double>(2, 0.0);
    
    auto decay_theta_LAB = std::vector<double>(2, 0.0);
    auto decay_phi_LAB = std::vector<double>(2, 0.0);
    auto decay_theta_recoilCOM = std::vector<double>(2, 0.0);
    auto decay_phi_recoilCOM = std::vector<double>(2, 0.0);
    auto decay_theta_reactionCOM = std::vector<double>(2, 0.0);
    auto decay_phi_reactionCOM = std::vector<double>(2, 0.0);
    
    //------------------------------------------------------------------------------------
    //      Polar angles (recoilCOM)

//    std::cout << "fAngDist_2plus_9870->Eval(10.0): " << fAngDist_2plus_9870->Eval(10.0) << std::endl;
//    std::cout << "fAngDist_2plus_9870->GetMaximum(): " << fAngDist_2plus_9870->GetMaximum() << std::endl;

//    fAngDist_2plus_9870->SaveAs("test.root");
//    fAngDist_4plus_14079->SaveAs("test.root");

    //------------
//    fAngDist_1minus
//    fAngDist_2plus_9870
//    fAngDist_3minus
//    fAngDist_4plus
    
//    fAngDist_1minus->GetMaximum(): 0.3849
//    fAngDist_2plus_9870->GetMaximum(): 0.25
//    fAngDist_3minus->GetMaximum(): 5.04066
//    fAngDist_4plus->GetMaximum(): 0.140625
    
    double xMax = 0.3849 + 0.001;   // fAngDist_1minus
//    double xMax = 0.25 + 0.001;     // fAngDist_2plus_9870
//    double xMax = 5.04066 + 0.001;  // fAngDist_3minus
//    double xMax = 0.140625 + 0.001; // fAngDist_4plus

//    double xMax = fAngDist_1minus->GetMaximum() + 0.001;

    double thetaTest = 0.0; // deg
    double functionTest = 0.0; // deg
    
    bool found = false;
    while(!found)
    {
        thetaTest = G4RandFlat::shoot(0.0, 180.0); // deg
        functionTest = G4RandFlat::shoot(xMax); // deg

        if(functionTest<=fAngDist_1minus->Eval(thetaTest))
        {
            found = true;
            decay_theta_recoilCOM[0] = thetaTest;
        }
    }
    
    
//    bool draw = true;
//    
//    if(draw)
//    {
//        std::lock_guard<std::mutex> guard(angDist_mutex);
//        decay_theta_recoilCOM[0] = fAngDist_2plus_9870->GetRandom(0.0, 180.0);
//    }

//    decay_theta_recoilCOM[0] = fAngDist_2plus_9870->GetRandom(0.0, 180.0);

    
    //------------
    /*
    double xMax = 0.14062500 + 0.001;
    double thetaTest = 0.0; // deg
    double functionTest = 0.0; // deg
    
    bool found = false;
    while(!found)
    {
        thetaTest = G4RandFlat::shoot(0.0, 180.0); // deg
        functionTest = G4RandFlat::shoot(xMax); // deg
        
        if(functionTest<=fAngDist_4plus_14079->Eval(thetaTest))
        {
            found = true;
            decay_theta_recoilCOM[0] = thetaTest;
        }
    }
    */

    
    
    //decay_theta_recoilCOM[0] = fAngDist_4plus_14079->GetRandom(0.0, 180.0);
    
    
    //------------
//    decay_theta_recoilCOM[0] = acos(1 - (2.0*G4UniformRand()))/deg;
//    decay_theta_recoilCOM[0] = fAngDist_2plus_9870->GetRandom();
//    decay_theta_recoilCOM[0] = fAngDist_4plus_13300->GetRandom(0.0, 180.0);
    decay_theta_recoilCOM[1] = 180.0 - decay_theta_recoilCOM[0];
    
    //------------------------------------------------------------------------------------
    //      Azimuthal angles (LAB)
    decay_phi_recoilCOM[0] = 360.0*G4UniformRand(); // deg
    decay_phi_recoilCOM[1] = decay_phi_recoilCOM[0] - 180.0; // deg
    
    if(decay_phi_recoilCOM[1]<0.0)
    {
        decay_phi_recoilCOM[1] += 360.0;
    }
    
    //------------------------------------------------------------------------------------
    //reaction_P = std::vector<double>(4, 0.0);
    CalculateBinaryDecayKinematics(recoilExcitationEnergy, 0.0, 7.36659, reaction_P[3], reaction_E[3], decay_A, reaction_theta_LAB[3], reaction_phi_LAB[3], decay_theta_recoilCOM, decay_phi_recoilCOM, decay_theta_LAB, decay_phi_LAB, decay_T);
    
    //------------------------------------------------------------------------------------
    double mx, my, mz;
    
    mx = sin(decay_theta_LAB[0]*deg)*cos(decay_phi_LAB[0]*deg);
    my = sin(decay_theta_LAB[0]*deg)*sin(decay_phi_LAB[0]*deg);
    mz = cos(decay_theta_LAB[0]*deg);
    
    fEventAction->SetInitialParticleKineticEnergy(decay_T[0]*MeV); // MeV
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
    fParticleGun->SetParticleEnergy(decay_T[0]*MeV);
    
    fEventAction->SetDecay_theta_LAB(decay_theta_LAB);
    fEventAction->SetDecay_phi_LAB(decay_phi_LAB);
    fEventAction->SetDecay_theta_reactionCOM(decay_theta_reactionCOM);
    fEventAction->SetDecay_phi_reactionCOM(decay_phi_reactionCOM);
    fEventAction->SetDecay_theta_recoilCOM(decay_theta_recoilCOM);
    fEventAction->SetDecay_phi_recoilCOM(decay_phi_recoilCOM);

    fParticleGun->GeneratePrimaryVertex(anEvent);
    
    
//    G4cout << "decay_T[0]: " << decay_T[0] << G4endl;
//    G4cout << "fParticleGun->GetParticleMomentum(): " << fParticleGun->GetParticleMomentum() << G4endl;
//    G4cout << "decay_theta_LAB[0]: " << decay_theta_LAB[0] << G4endl;

    //========================================================================================================================
    //========================================================================================================================
    //      8Be simulation SCATTERING REACTION
    /*
    int nReactionProducts = 4;
    std::vector<int>    reaction_Z{1, 6, 1, 6}; // e
    std::vector<double> reaction_A{1.00727647, 14.00324199, 3.01550071, 12}; // amu
    auto reaction_P = std::vector<double>(4, 0.0);
    std::vector<double> reaction_T{140.0, 0.0, 0.0, 0.0}; // MeV
    std::vector<double> reaction_Ex{0.0, 0.0, 0.0, 0.0}; // MeV
    auto reaction_E = std::vector<double>(4, 0.0);
    
    //------------------------------------------------------------------------------------
    //      Polar angles (LAB)
    std::vector<double> reaction_theta_LAB{0.0, 0.0, acos(1 - (6.09173503004822869e-04*G4UniformRand()))/deg, 0.0}; // deg
    auto reaction_theta_reactionCOM = std::vector<double>(4, 0.0);
    
    //------------------------------------------------------------------------------------
    //      Azimuthal angles (LAB)
    double reaction_ejectile_phi_LAB = 360.0*G4UniformRand(); // deg
    double reaction_recoil_phi_LAB = reaction_ejectile_phi_LAB - 180.0; // deg
    if(reaction_recoil_phi_LAB<0.0)
    {
        reaction_recoil_phi_LAB += 360.0;
    }
    
    std::vector<double> reaction_phi_LAB{0.0, 0.0, reaction_ejectile_phi_LAB, reaction_recoil_phi_LAB}; // deg
    auto reaction_phi_reactionCOM = std::vector<double>(4, 0.0);
    
    //------------------------------------------------------------------------------------
    //      Recoil Nucleux Excitation Energy
    //double recoilExcitationEnergy = 0.0; // MeV
    double recoilExcitationEnergy = 7.654; // MeV
    //double recoilExcitationEnergy = 7.36659; // MeV
    //double recoilExcitationEnergy = 7.36659 + 0.5; // MeV
    //double recoilExcitationEnergy = 12.0; // MeV
    //double recoilExcitationEnergy = 20.0*G4UniformRand(); // MeV

    double daughterNucleusFinalStateEnergy = 0.0; // MeV

    //------------------------------------------------------------------------------------
    //      Executing the Binary Relativistic Kinematics code
    BiRelKin(&reaction_A[0], &reaction_T[0], &reaction_E[0], &reaction_P[0], reaction_theta_LAB[2], reaction_theta_LAB[3], recoilExcitationEnergy+daughterNucleusFinalStateEnergy);
    
    //G4cout << "reaction_T[3]: " << reaction_T[3] << G4endl;

    reaction_Ex[3] = recoilExcitationEnergy;
    
    //------------------------------------------------------------------------------------
    //      Filling the vectors in the EventAction object
    fEventAction->SetNReactionProducts(nReactionProducts);
    fEventAction->SetReaction_Z(reaction_Z);
    fEventAction->SetReaction_A(reaction_A);
    fEventAction->SetReaction_P(reaction_P);
    fEventAction->SetReaction_T(reaction_T);
    fEventAction->SetReaction_Ex(reaction_Ex);
    
    fEventAction->SetReaction_theta_LAB(reaction_theta_LAB);
    fEventAction->SetReaction_phi_LAB(reaction_phi_LAB);
    fEventAction->SetReaction_theta_reactionCOM(reaction_theta_reactionCOM);
    fEventAction->SetReaction_phi_reactionCOM(reaction_phi_reactionCOM);
    
    
    //------------------------------------------------------------------------------------
    //      Recoil Nucleus Decay
    
    int nDecayProducts = 2;
    std::vector<int>    decay_Z{2, 4}; // e
    std::vector<double> decay_A{4.00260325413, 8.005305102}; // amu
    auto decay_P = std::vector<double>(2, 0.0);
    auto decay_T = std::vector<double>(2, 0.0);
    
    auto decay_theta_LAB = std::vector<double>(2, 0.0);
    auto decay_phi_LAB = std::vector<double>(2, 0.0);
    auto decay_theta_recoilCOM = std::vector<double>(2, 0.0);
    auto decay_phi_recoilCOM = std::vector<double>(2, 0.0);
    auto decay_theta_reactionCOM = std::vector<double>(2, 0.0);
    auto decay_phi_reactionCOM = std::vector<double>(2, 0.0);
    
    //------------------------------------------------------------------------------------
    //      Polar angles (recoilCOM)
    decay_theta_recoilCOM[0] = acos(1 - (2.0*G4UniformRand()))/deg;
    decay_theta_recoilCOM[1] = 180.0 - decay_theta_recoilCOM[0];
    
    //------------------------------------------------------------------------------------
    //      Azimuthal angles (LAB)
    decay_phi_recoilCOM[0] = 360.0*G4UniformRand(); // deg
    decay_phi_recoilCOM[1] = decay_phi_recoilCOM[0] - 180.0; // deg
    
    if(decay_phi_recoilCOM[1]<0.0)
    {
        decay_phi_recoilCOM[1] += 360.0;
    }
    
    //------------------------------------------------------------------------------------
    //reaction_P = std::vector<double>(4, 0.0);
    //std::cout << "reaction_P[3]: " << reaction_P[3] << std::endl;
    
    CalculateBinaryDecayKinematics(recoilExcitationEnergy, daughterNucleusFinalStateEnergy, 7.36659, reaction_P[3], reaction_E[3], decay_A, reaction_theta_LAB[3], reaction_phi_LAB[3], decay_theta_recoilCOM, decay_phi_recoilCOM, decay_theta_LAB, decay_phi_LAB, decay_T);

    //------------------------------------------------------------------------------------
    int nDecayProducts_8Be = 2;
    std::vector<int>    decay_Z_8Be{2, 2}; // e
    std::vector<double> decay_A_8Be{4.00260325413, 4.00260325413}; // amu
    auto decay_P_8Be = std::vector<double>(2, 0.0);
    auto decay_T_8Be = std::vector<double>(2, 0.0);
    
    auto decay_theta_LAB_8Be = std::vector<double>(2, 0.0);
    auto decay_phi_LAB_8Be = std::vector<double>(2, 0.0);
    auto decay_theta_recoilCOM_8Be = std::vector<double>(2, 0.0);
    auto decay_phi_recoilCOM_8Be = std::vector<double>(2, 0.0);
    auto decay_theta_reactionCOM_8Be = std::vector<double>(2, 0.0);
    auto decay_phi_reactionCOM_8Be = std::vector<double>(2, 0.0);

    //------------------------------------------------------------------------------------
    //      Polar angles (recoilCOM)
    decay_theta_recoilCOM_8Be[0] = acos(1 - (2.0*G4UniformRand()))/deg;
    decay_theta_recoilCOM_8Be[1] = 180.0 - decay_theta_recoilCOM_8Be[0];

    //std::cout << "decay_theta_recoilCOM_8Be[0]: " << decay_theta_recoilCOM_8Be[0] << std::endl;
    //std::cout << "decay_theta_recoilCOM_8Be[1]: " << decay_theta_recoilCOM_8Be[1] << std::endl;
    
    //------------------------------------------------------------------------------------
    //      Azimuthal angles (LAB)
    decay_phi_recoilCOM_8Be[0] = 360.0*G4UniformRand(); // deg
    decay_phi_recoilCOM_8Be[1] = decay_phi_recoilCOM_8Be[0] - 180.0; // deg
    
    if(decay_phi_recoilCOM_8Be[1]<0.0)
    {
        decay_phi_recoilCOM_8Be[1] += 360.0;
    }
    
    //------------------------------------------------------------------------------------
    //reaction_P = std::vector<double>(4, 0.0);
    
    //G4cout << "decay_T[0]: " << decay_T[0] << G4endl;
    //G4cout << "decay_T[1]: " << decay_T[1] << G4endl;
    
    double c2 = 931.49410242;     // MeV/u, c^2
    double P_8Be = sqrt(2.0*decay_A[1]*decay_T[1]);
    double E_8Be = sqrt(pow(P_8Be, 2.0)*c2 + pow(decay_A[1]*c2, 2.0));
    
    //double P_8Be = 0.0;
    //double E_8Be = sqrt(pow(P_8Be*931.5, 2.0) + pow(decay_A[1]*931.5, 2.0));

    //std::cout << "P_8Be: " << P_8Be << std::endl;
    //std::cout << "reaction_theta_LAB[3]: " << reaction_theta_LAB[3] << std::endl;
    
    //void PrimaryGeneratorAction::CalculateBinaryDecayKinematics(double eX, double eX_residual, double sepE, double initialScatteringReactionMomentum, double initialScatteringReactionEnergy, std::vector<double> decay_A, double theta_recoil_LAB, double phi_recoil_LAB, std::vector<double>& decay_theta_recoilCOM, std::vector<double>& decay_phi_recoilCOM, std::vector<double>& decay_theta_LAB, std::vector<double>& decay_phi_LAB, std::vector<double>& decay_T)
    
    CalculateBinaryDecayKinematics(0.0, 0.0, -0.09184, P_8Be, E_8Be, decay_A_8Be, reaction_theta_LAB[3], reaction_phi_LAB[3], decay_theta_recoilCOM_8Be, decay_phi_recoilCOM_8Be, decay_theta_LAB_8Be, decay_phi_LAB_8Be, decay_T_8Be);

    //decay_theta_LAB_8Be[0] = decay_theta_recoilCOM_8Be[0];
    //decay_theta_LAB_8Be[1] = decay_theta_recoilCOM_8Be[1];
    
    //decay_phi_LAB_8Be[0] = decay_phi_recoilCOM_8Be[0];
    //decay_phi_LAB_8Be[1] = decay_phi_recoilCOM_8Be[1];
    
    //------------------------------------------------------------------------------------
    double mx, my, mz;
    
    mx = sin(decay_theta_LAB_8Be[0]*deg)*cos(decay_phi_LAB_8Be[0]*deg);
    my = sin(decay_theta_LAB_8Be[0]*deg)*sin(decay_phi_LAB_8Be[0]*deg);
    mz = cos(decay_theta_LAB_8Be[0]*deg);
    
    fEventAction->SetInitialParticleKineticEnergy(decay_T_8Be[0]*MeV); // MeV
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
    //G4cout << "decay_T_8Be[0]: " << decay_T_8Be[0] << G4endl;
    //G4cout << "decay_theta_LAB_8Be[0]: " << decay_theta_LAB_8Be[0] << G4endl;
    
    fParticleGun->GeneratePrimaryVertex(anEvent);
    
    //------------------------------------------------------------------------------------
    mx = sin(decay_theta_LAB_8Be[1]*deg)*cos(decay_phi_LAB_8Be[1]*deg);
    my = sin(decay_theta_LAB_8Be[1]*deg)*sin(decay_phi_LAB_8Be[1]*deg);
    mz = cos(decay_theta_LAB_8Be[1]*deg);
    
    fEventAction->SetInitialParticleKineticEnergy(decay_T_8Be[1]*MeV); // MeV
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
    //G4cout << "decay_T_8Be[1]: " << decay_T_8Be[1] << G4endl;
    //G4cout << "decay_theta_LAB_8Be[1]: " << decay_theta_LAB_8Be[1] << G4endl;

    fParticleGun->GeneratePrimaryVertex(anEvent);
    */
    
    //------------------------------------------------------------------------------------
//    fEventAction->SetNDecayProducts(nDecayProducts);
//    fEventAction->SetDecay_Z(decay_Z);
//    fEventAction->SetDecay_A(decay_A);
//    fEventAction->SetDecay_P(decay_P);
//    fEventAction->SetDecay_T(decay_T);
//    
//    fEventAction->SetDecay_theta_LAB(decay_theta_LAB);
//    fEventAction->SetDecay_phi_LAB(decay_phi_LAB);
//    fEventAction->SetDecay_theta_reactionCOM(decay_theta_reactionCOM);
//    fEventAction->SetDecay_phi_reactionCOM(decay_phi_reactionCOM);
//    fEventAction->SetDecay_theta_recoilCOM(decay_theta_recoilCOM);
//    fEventAction->SetDecay_phi_recoilCOM(decay_phi_recoilCOM);
//    
//    double mx, my, mz;
//    
//    mx = sin(decay_theta_LAB[1]*deg)*cos(decay_phi_LAB[1]*deg);
//    my = sin(decay_theta_LAB[1]*deg)*sin(decay_phi_LAB[1]*deg);
//    mz = cos(decay_theta_LAB[1]*deg);
//    
//    fEventAction->SetInitialParticleKineticEnergy(decay_T[1]*MeV); // MeV
//    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));

    
    
    //========================================================================================================================
    //========================================================================================================================
    //      PHD DEFENCE VISUALISATION
    //      PR240: alpha0 decay, 0+ -> 0+, alpha decay from 12C
    /*
    static int eventN = 0;
    
    int nReactionProducts = 4;
    std::vector<int>    reaction_Z{1, 6, 1, 6}; // e
    std::vector<double> reaction_A{1.00727647, 14.00324199, 3.01550071, 12}; // amu
    auto reaction_P = std::vector<double>(4, 0.0);
    std::vector<double> reaction_T{140.0, 0.0, 0.0, 0.0}; // MeV
    std::vector<double> reaction_Ex{0.0, 0.0, 0.0, 0.0}; // MeV
    auto reaction_E = std::vector<double>(4, 0.0);
    
    //------------------------------------------------------------------------------------
    //      Polar angles (LAB)
    std::vector<double> reaction_theta_LAB{0.0, 0.0, acos(1 - (6.09173503004822869e-04*G4UniformRand()))/deg, 0.0}; // deg
    //auto reaction_theta_reactionCOM = std::vector<double>(4, 0.0);
    
    //------------------------------------------------------------------------------------
    //      Azimuthal angles (LAB)
    double reaction_ejectile_phi_LAB = 360.0*G4UniformRand(); // deg
    double reaction_recoil_phi_LAB = reaction_ejectile_phi_LAB - 180.0; // deg
    if(reaction_recoil_phi_LAB<0.0)
    {
        reaction_recoil_phi_LAB += 360.0;
    }
    
    std::vector<double> reaction_phi_LAB{0.0, 0.0, reaction_ejectile_phi_LAB, reaction_recoil_phi_LAB}; // deg
    //auto reaction_phi_reactionCOM = std::vector<double>(4, 0.0);
    
    //------------------------------------------------------------------------------------
    //      Recoil Nucleux Excitation Energy
    double recoilExcitationEnergy = 7.654; // MeV
    //double recoilExcitationEnergy = 9.0; // MeV
    
    if(eventN<2000)
    {
        recoilExcitationEnergy = 9.0;
    }
    else if(eventN<4000)
    {
        recoilExcitationEnergy = 7.654;
    }

    //------------------------------------------------------------------------------------
    
    //------------------------------------------------------------------------------------
    //      Executing the Binary Relativistic Kinematics code
    BiRelKin(&reaction_A[0], &reaction_T[0], &reaction_E[0], &reaction_P[0], reaction_theta_LAB[2], reaction_theta_LAB[3], recoilExcitationEnergy);
    
    //G4cout << "reaction_T[3]: " << reaction_T[3] << G4endl;
    
    reaction_Ex[3] = recoilExcitationEnergy;
    
    //------------------------------------------------------------------------------------
    //      Filling the vectors in the EventAction object
    fEventAction->SetNReactionProducts(nReactionProducts);
    fEventAction->SetReaction_Z(reaction_Z);
    fEventAction->SetReaction_A(reaction_A);
    fEventAction->SetReaction_P(reaction_P);
    fEventAction->SetReaction_T(reaction_T);
    fEventAction->SetReaction_Ex(reaction_Ex);
    
    fEventAction->SetReaction_theta_LAB(reaction_theta_LAB);
    fEventAction->SetReaction_phi_LAB(reaction_phi_LAB);
    //fEventAction->SetReaction_theta_reactionCOM(reaction_theta_reactionCOM);
    //fEventAction->SetReaction_phi_reactionCOM(reaction_phi_reactionCOM);
    
    
    //------------------------------------------------------------------------------------
    //      Recoil Nucleus Decay
    
    int nDecayProducts = 2;
    std::vector<int>    decay_Z{2, 4}; // e
    std::vector<double> decay_A{4.00260325413, 8.005305102}; // amu
    auto decay_P = std::vector<double>(2, 0.0);
    auto decay_T = std::vector<double>(2, 0.0);
    
    auto decay_theta_LAB = std::vector<double>(2, 0.0);
    auto decay_phi_LAB = std::vector<double>(2, 0.0);
    auto decay_theta_recoilCOM = std::vector<double>(2, 0.0);
    auto decay_phi_recoilCOM = std::vector<double>(2, 0.0);
    auto decay_theta_reactionCOM = std::vector<double>(2, 0.0);
    auto decay_phi_reactionCOM = std::vector<double>(2, 0.0);
    
    //------------------------------------------------------------------------------------
    //      Polar angles (recoilCOM)
    decay_theta_recoilCOM[0] = acos(1 - (2.0*G4UniformRand()))/deg;
    decay_theta_recoilCOM[1] = 180.0 - decay_theta_recoilCOM[0];
    
    //------------------------------------------------------------------------------------
    //      Azimuthal angles (LAB)
    decay_phi_recoilCOM[0] = 360.0*G4UniformRand(); // deg
    decay_phi_recoilCOM[1] = decay_phi_recoilCOM[0] - 180.0; // deg
    
    if(decay_phi_recoilCOM[1]<0.0)
    {
        decay_phi_recoilCOM[1] += 360.0;
    }
    
    //------------------------------------------------------------------------------------
    //reaction_P = std::vector<double>(4, 0.0);
    if((eventN>=0 && eventN<1000) || (eventN>=2000 && eventN<3000))
    {
        usleep(500);
        
        //------------------------------------------------------------------------------------
        double mx, my, mz;
        
        mx = sin(decay_theta_recoilCOM[0]*deg)*cos(decay_phi_recoilCOM[0]*deg);
        my = sin(decay_theta_recoilCOM[0]*deg)*sin(decay_phi_recoilCOM[0]*deg);
        mz = cos(decay_theta_recoilCOM[0]*deg);
        
        fEventAction->SetInitialParticleKineticEnergy((8.005305102/12.0)*(recoilExcitationEnergy-7.36659)*MeV); // MeV
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
        fParticleGun->SetParticleEnergy((8.005305102/12.0)*(recoilExcitationEnergy-7.36659)*MeV);
        
        fEventAction->SetDecay_theta_LAB(decay_theta_LAB);
        fEventAction->SetDecay_phi_LAB(decay_phi_LAB);
        fEventAction->SetDecay_theta_reactionCOM(decay_theta_reactionCOM);
        fEventAction->SetDecay_phi_reactionCOM(decay_phi_reactionCOM);
        fEventAction->SetDecay_theta_recoilCOM(decay_theta_recoilCOM);
        fEventAction->SetDecay_phi_recoilCOM(decay_phi_recoilCOM);
        
        fParticleGun->GeneratePrimaryVertex(anEvent);
    }
    else if((eventN>=1000 && eventN<2000) || (eventN>=3000 && eventN<4000))
    {
        usleep(500);

        CalculateBinaryDecayKinematics(recoilExcitationEnergy, 0.0, 7.36659, reaction_P[3], reaction_E[3], decay_A, reaction_theta_LAB[3], reaction_phi_LAB[3], decay_theta_recoilCOM, decay_phi_recoilCOM, decay_theta_LAB, decay_phi_LAB, decay_T);
        
        //------------------------------------------------------------------------------------
        double mx, my, mz;
        
        mx = sin(decay_theta_LAB[0]*deg)*cos(decay_phi_LAB[0]*deg);
        my = sin(decay_theta_LAB[0]*deg)*sin(decay_phi_LAB[0]*deg);
        mz = cos(decay_theta_LAB[0]*deg);
        
        fEventAction->SetInitialParticleKineticEnergy(decay_T[0]*MeV); // MeV
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
        fParticleGun->SetParticleEnergy(decay_T[0]*MeV);
        
        fEventAction->SetDecay_theta_LAB(decay_theta_LAB);
        fEventAction->SetDecay_phi_LAB(decay_phi_LAB);
        fEventAction->SetDecay_theta_reactionCOM(decay_theta_reactionCOM);
        fEventAction->SetDecay_phi_reactionCOM(decay_phi_reactionCOM);
        fEventAction->SetDecay_theta_recoilCOM(decay_theta_recoilCOM);
        fEventAction->SetDecay_phi_recoilCOM(decay_phi_recoilCOM);
        
        fParticleGun->GeneratePrimaryVertex(anEvent);
    }
    
    eventN++;
    */
    
    //========================================================================================
    
    
    //------------------------------------------------------------------------------------
    //      INITIAL SCATTERING REACTION
    /*
    int nReactionProducts = 4;
    std::vector<int>    reaction_Z{1, 6, 1, 6}; // e
    std::vector<double> reaction_A{1.00727647, 14.00324199, 3.01550071, 12}; // amu
    auto reaction_P = std::vector<double>(4, 0.0);
    std::vector<double> reaction_T{140.0, 0.0, 0.0, 0.0}; // MeV
    std::vector<double> reaction_Ex{0.0, 0.0, 0.0, 0.0}; // MeV
    auto reaction_E = std::vector<double>(4, 0.0);
    
    //------------------------------------------------------------------------------------
    //      Polar angles (LAB)
    std::vector<double> reaction_theta_LAB{0.0, 0.0, acos(1 - (6.09173503004822869e-04*G4UniformRand()))/deg, 0.0}; // deg
    auto reaction_theta_reactionCOM = std::vector<double>(4, 0.0);
    
    //------------------------------------------------------------------------------------
    //      Azimuthal angles (LAB)
    double reaction_ejectile_phi_LAB = 360.0*G4UniformRand(); // deg
    double reaction_recoil_phi_LAB = reaction_ejectile_phi_LAB - 180.0; // deg
    if(reaction_recoil_phi_LAB<0.0)
    {
        reaction_recoil_phi_LAB += 360.0;
    }
    
    std::vector<double> reaction_phi_LAB{0.0, 0.0, reaction_ejectile_phi_LAB, reaction_recoil_phi_LAB}; // deg
    auto reaction_phi_reactionCOM = std::vector<double>(4, 0.0);

    //------------------------------------------------------------------------------------
    //      Recoil Nucleux Excitation Energy
    double recoilExcitationEnergy = 7.654; // MeV
    
    //------------------------------------------------------------------------------------
    //      Executing the Binary Relativistic Kinematics code
    BiRelKin(&reaction_A[0], &reaction_T[0], &reaction_E[0], &reaction_P[0], reaction_theta_LAB[2], reaction_theta_LAB[3], recoilExcitationEnergy);
    
    reaction_Ex[3] = recoilExcitationEnergy;
    
    //------------------------------------------------------------------------------------
    //      Filling the vectors in the EventAction object
    fEventAction->SetNReactionProducts(nReactionProducts);
    fEventAction->SetReaction_Z(reaction_Z);
    fEventAction->SetReaction_A(reaction_A);
    fEventAction->SetReaction_P(reaction_P);
    fEventAction->SetReaction_T(reaction_T);
    fEventAction->SetReaction_Ex(reaction_Ex);
    
    fEventAction->SetReaction_theta_LAB(reaction_theta_LAB);
    fEventAction->SetReaction_phi_LAB(reaction_phi_LAB);
    fEventAction->SetReaction_theta_reactionCOM(reaction_theta_reactionCOM);
    fEventAction->SetReaction_phi_reactionCOM(reaction_phi_reactionCOM);
   
    
    //------------------------------------------------------------------------------------
    //      Recoil Nucleus Decay
    
    int nDecayProducts = 2;
    std::vector<int>    decay_Z{2, 4}; // e
    std::vector<double> decay_A{4.00260325413, 8.005305102}; // amu
    auto decay_P = std::vector<double>(2, 0.0);
    auto decay_T = std::vector<double>(2, 0.0);

    auto decay_theta_LAB = std::vector<double>(2, 0.0);
    auto decay_phi_LAB = std::vector<double>(2, 0.0);
    auto decay_theta_recoilCOM = std::vector<double>(2, 0.0);
    auto decay_phi_recoilCOM = std::vector<double>(2, 0.0);
    auto decay_theta_reactionCOM = std::vector<double>(2, 0.0);
    auto decay_phi_reactionCOM = std::vector<double>(2, 0.0);

    //------------------------------------------------------------------------------------
    //      Polar angles (recoilCOM)
    decay_theta_recoilCOM[0] = acos(1 - (2.0*G4UniformRand()))/deg;
    decay_theta_recoilCOM[1] = 180.0 - decay_theta_recoilCOM[0];
    
    //------------------------------------------------------------------------------------
    //      Azimuthal angles (LAB)
    decay_phi_recoilCOM[0] = 360.0*G4UniformRand(); // deg
    decay_phi_recoilCOM[1] = decay_phi_recoilCOM[0] - 180.0; // deg
    
    if(decay_phi_recoilCOM[1]<0.0)
    {
        decay_phi_recoilCOM[1] += 360.0;
    }
    
    //------------------------------------------------------------------------------------
    reaction_P = std::vector<double>(4, 0.0);
    CalculateBinaryDecayKinematics(recoilExcitationEnergy, 0.0, 7.36659, reaction_P[3], reaction_E[3], decay_A, reaction_theta_LAB[3], reaction_phi_LAB[3], decay_theta_recoilCOM, decay_phi_recoilCOM, decay_theta_LAB, decay_phi_LAB, decay_T);
    
    //------------------------------------------------------------------------------------
    fEventAction->SetNDecayProducts(nDecayProducts);
    fEventAction->SetDecay_Z(decay_Z);
    fEventAction->SetDecay_A(decay_A);
    fEventAction->SetDecay_P(decay_P);
    fEventAction->SetDecay_T(decay_T);
    
    fEventAction->SetDecay_theta_LAB(decay_theta_LAB);
    fEventAction->SetDecay_phi_LAB(decay_phi_LAB);
    fEventAction->SetDecay_theta_reactionCOM(decay_theta_reactionCOM);
    fEventAction->SetDecay_phi_reactionCOM(decay_phi_reactionCOM);
    fEventAction->SetDecay_theta_recoilCOM(decay_theta_recoilCOM);
    fEventAction->SetDecay_phi_recoilCOM(decay_phi_recoilCOM);
    
    double mx, my, mz;
    
    
    mx = sin(decay_theta_LAB[0]*deg)*cos(decay_phi_LAB[0]*deg);
    my = sin(decay_theta_LAB[0]*deg)*sin(decay_phi_LAB[0]*deg);
    mz = cos(decay_theta_LAB[0]*deg);
    
    fEventAction->SetInitialParticleKineticEnergy(decay_T[0]); // MeV
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
    */
    
    
    //====================================================================================
    //      SOLID ANGLE DETERMINATION (LAB FRAME)
    /*
    double mx, my, mz;

    auto decay_theta_LAB = std::vector<double>(1, 0.0);
    auto decay_phi_LAB = std::vector<double>(1, 0.0);

    //decay_theta_LAB[0] = acos(1 - (2.0*G4UniformRand()))/deg; // 0.0->180.0 deg
    decay_theta_LAB[0] = acos(1 - (1.0*G4UniformRand() + 1.0))/deg; // 90.0->180.0 deg
    decay_phi_LAB[0] = 360.0*G4UniformRand(); // deg

    fEventAction->SetDecay_theta_LAB(decay_theta_LAB);
    fEventAction->SetDecay_phi_LAB(decay_phi_LAB);
    
    mx = sin(decay_theta_LAB[0]*deg)*cos(decay_phi_LAB[0]*deg);
    my = sin(decay_theta_LAB[0]*deg)*sin(decay_phi_LAB[0]*deg);
    mz = cos(decay_theta_LAB[0]*deg);
    
    //fEventAction->SetInitialParticleKineticEnergy(1.0); // MeV
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
    
    fParticleGun->GeneratePrimaryVertex(anEvent);
    */
    
    /*
    G4double range;
    G4double mx, my, mz;
    G4bool acceptance = false;
    
    G4double angDist_theta = 0.0;
    G4double function = 0.0;
    
    while(!acceptance)
    {
        //angDist_theta = 180.0*G4UniformRand();
        angDist_theta = 90.0*G4UniformRand() + 90.0;
        range = 1.00*G4UniformRand();
        function = sin(angDist_theta*deg);
        
        if(range<function)
        {
            acceptance = true;
        }
    }
    
    phi = 360.*G4UniformRand();
    
    mx = sin(angDist_theta*deg)*cos(phi*deg);
    my = sin(angDist_theta*deg)*sin(phi*deg);
    mz = cos(angDist_theta*deg);
    
    auto decay_theta_LAB = std::vector<double>(1, 0.0);
    auto decay_phi_LAB = std::vector<double>(1, 0.0);
    
    decay_theta_LAB[0] = angDist_theta; // 90.0->180.0 deg
    decay_phi_LAB[0] = phi; // deg

    fEventAction->SetDecay_theta_LAB(decay_theta_LAB);
    fEventAction->SetDecay_phi_LAB(decay_phi_LAB);

    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
    */
    
    ////////////////////////////////////////////////////
    ////    NEW CAKE SIMULATIONS
    
    /*
    //------------------------------------------------
    G4double theta = acos(1 - (2.0*G4UniformRand()))/deg; // 0.0->180.0
    //G4double theta = acos(1 - (1.0*G4UniformRand()))/deg; // 0.0->90.0 deg
    //G4double theta = acos(1 - (1.0*G4UniformRand() + 1.0))/deg; // 90.0->180.0 deg
    ////    Collimator - 2.0 deg range
    //G4double theta = acos(1 - (6.09173503004822869e-04*G4UniformRand()))/deg; // deg
    
    G4double phi = 360.0*G4UniformRand(); // deg
    
    
    fEventAction->SetInitialTheta_COM(theta);
    fEventAction->SetInitialPhi_COM(phi);
    
    
    //------------------------------------------------
    double theta_ejectile, phi_ejectile, theta_recoil_LAB, phi_recoil_LAB, eX;
    ////    Uniformly distributed theta_ejectile
    //theta_ejectile = G4RandFlat::shoot(-2., 2.); // deg
    //theta_ejectile = 0.1; // deg
    theta_ejectile = acos(1 - (6.09173503004822869e-04*G4UniformRand()))/deg; // deg
    phi_ejectile = G4RandFlat::shoot(0.0, 360.0); // deg
    
    CalculateBinaryRelativisticKinematics(theta_ejectile, phi_ejectile, theta_recoil_LAB, phi_recoil_LAB, eX);
    
    //G4cout << "theta_recoil_LAB: " << theta_recoil_LAB << std::endl;
    fEventAction->SetRecoilExcitationEnergy(eX);
    fEventAction->SetInitialTheta_LAB(theta_recoil_LAB);
    fEventAction->SetInitialPhi_LAB(phi_recoil_LAB);
    
    //------------------------------------------------
    //void PrimaryGeneratorAction::CalculateBinaryDecayKinematics(double eX, double eX_residual, double sepE, double mass_decayParticle, double mass_residualParticle, double theta_recoil_LAB, double phi_recoil_LAB, double theta_decayParticle_LAB, double phi_decayParticle_LAB, double &theta_decayParticle_recoilCOM, double &phi_decayParticle_recoilCOM)
    
    double theta_decayParticle_recoilCOM, phi_decayParticle_recoilCOM, theta_decayParticle_LAB, phi_decayParticle_LAB;
    theta_decayParticle_recoilCOM = acos(1 - (2.0*G4UniformRand()))/deg;
    phi_decayParticle_recoilCOM = G4RandFlat::shoot(0.0, 360.0);
    
    double initialT_decayParticle_LAB = 0.0;
    
    CalculateBinaryDecayKinematics(eX, 0.0, 7.16192, 4.001506179127, 12.0, theta_recoil_LAB, phi_recoil_LAB, theta_decayParticle_recoilCOM, phi_decayParticle_recoilCOM, theta_decayParticle_LAB, phi_decayParticle_LAB, initialT_decayParticle_LAB);
    
    mx = sin(theta_decayParticle_LAB*deg)*cos(phi_decayParticle_LAB*deg);
    my = sin(theta_decayParticle_LAB*deg)*sin(phi_decayParticle_LAB*deg);
    mz = cos(theta_decayParticle_LAB*deg);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
    fParticleGun->SetParticleEnergy(initialT_decayParticle_LAB*MeV);
    
    fEventAction->SetInitialT_decayParticle_LAB(initialT_decayParticle_LAB);
    */
    
    //------------------------------------------------
    /*
     mx = sin(theta*deg)*cos(phi*deg);
     my = sin(theta*deg)*sin(phi*deg);
     mz = cos(theta*deg);
     fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
     */
    //------------------------------------------------
    
    
    
    
    
    ////////////////////////////////////////////////////////////
    ////    Distributions - Acceptance/Rejection Method
    
    /*
     G4double range;
     G4double mx, my, mz;
     G4bool acceptance = false;
     
     G4double angDist_theta = 0.0;
     G4double function = 0.0;
     
     while(!acceptance)
     {
     angDist_theta = 180.0*G4UniformRand();
     //angDist_theta = 90.0*G4UniformRand() + 90.0;
     range = 1.001*G4UniformRand();
     
     //  EvalAngCor_alpha0_2plus_0plus_11520_L2: 0.25
     //  EvalAngCor_alpha0_0plus_0plus_12049_L0: 0.95
     //  EvalAngCor_alpha0_0plus_0plus_15097_L0: 0.95
     //  EvalAngCor_alpha1_2plus_2plus_15097_L0: 0.92
     
     
     //function = sin(angDist_theta*deg);
     
     //function = EvalAngCor_alpha0_2plus_0plus_11520_L2(angDist_theta)*sin(angDist_theta*deg);
     //function = EvalAngCor_alpha0_0plus_0plus_12049_L0(angDist_theta)*sin(angDist_theta*deg);
     //function = EvalAngCor_alpha0_0plus_0plus_15097_L0(angDist_theta)*sin(angDist_theta*deg);
     //function = EvalAngCor_alpha1_2plus_2plus_15097_L0(angDist_theta)*sin(angDist_theta*deg);
     
     function = EvalAngCor_p0_2plus_12minus_14926_L1(angDist_theta)*sin(angDist_theta*deg);
     
     if(range<function)
     {
     acceptance = true;
     }
     
     //G4cout << "Testing" << G4endl;
     }
     
     phi = 360.*G4UniformRand();
     
     mx = sin(angDist_theta*deg)*cos(phi*deg);
     my = sin(angDist_theta*deg)*sin(phi*deg);
     mz = cos(angDist_theta*deg);
     
     fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
     */
    
    
    
    
    ////    Rangefinder
    /*
     G4double highestValue = 0.0;
     G4double angDist_theta = 0.0;
     G4double function;
     
     for(G4double i=0; i<180.0; i+=0.1)
     {
     angDist_theta = i;
     
     //function = EvalAngCor_alpha0_2plus_0plus_11520_L2(angDist_theta)*sin(angDist_theta*deg);
     //function = EvalAngCor_alpha0_0plus_0plus_12049_L0(angDist_theta)*sin(angDist_theta*deg);
     function = EvalAngCor_alpha0_0plus_0plus_15097_L0(angDist_theta)*sin(angDist_theta*deg);
     //function = EvalAngCor_alpha1_2plus_2plus_15097_L0(angDist_theta)*sin(angDist_theta*deg);
     
     
     if(function>highestValue) highestValue = function;
     
     //G4cout << "EvalAngCor_alpha0_2plus_0plus_11520_L2:    " << function << G4endl;
     }
     
     G4cout << "Here is the highestValue:    " << highestValue << G4endl;
     */
    
    
    
    ////////////////////////////////////
    ////    Alpha Decay
    ////    Associated Legendre Polynomial
    
    ////    Initialise a pre-calculated angular distribution
    /*
     G4String angDist_name = "a0_15097_16O";
     G4int noOfPoints;
     
     initialiseAngDist_interpolated(angDist_name, noOfPoints);
     G4cout << "noOfPoints:    " << noOfPoints << G4endl;
     */
    
    /*
     G4double range;
     G4double mx, my, mz;
     G4bool acceptance = false;
     
     G4double angDist_theta = 0.0;
     G4double function = 0.0;
     
     while(!acceptance)
     {
     
     //angDist_theta = 180.0*G4UniformRand();
     angDist_theta = 90.0*G4UniformRand() + 90.0;
     range = 1.95*G4UniformRand();
     
     //range = 180.*G4UniformRand();
     
     ////====================================////
     ////        NEW (Bohr's Theorem)
     
     ////    alpha0 decay: 2+ -> 0+
     //function += pow( (2.0*0.203684)*(1.0/4.0)*sqrt(15.0/(2.0*M_PI))*pow(sin(angDist_theta), 2.0), 2.0);
     //function += pow( (0.592631)*(1.0/4.0)*sqrt(5.0/M_PI)*(3*cos(angDist_theta)*cos(angDist_theta) - 1), 2.0);
     //function *= sin(angDist_theta);
     
     //function += pow( (2.0*0.204898)*(1.0/4.0)*sqrt(15.0/(2.0*M_PI))*pow(sin(angDist_theta), 2.0), 2.0);
     //function += pow( (0.590203)*(1.0/4.0)*sqrt(5.0/M_PI)*(3*cos(angDist_theta)*cos(angDist_theta) - 1), 2.0);
     //function *= sin(angDist_theta);
     
     
     ////    alpha0 decay: 2+ -> 0+
     //function = 0.0;
     //function += pow( (2.0*0.204898)*(1.0/4.0)*sqrt(15.0/(2.0*M_PI))*pow(sin(angDist_theta), 2.0), 2.0);
     //function += pow( (0.590203)*(1.0/4.0)*sqrt(5.0/M_PI)*(3*cos(angDist_theta)*cos(angDist_theta) - 1), 2.0);
     //function *= sin(angDist_theta);
     
     function = evaluateAngDist_interpolated(angDist_theta)*sin(angDist_theta*deg);
     //function = 0.01*angDist_theta;
     
     //if(angDist_theta>90.0) G4cout << "PROBLEM: angDist_theta>90.0" << G4endl;
     
     //function = sin(angDist_theta);
     //function = sin(angDist_theta) + 0.1*angDist_theta;
     
     //function = 0.2*angDist_theta;
     
     //function = angDist_theta;
     
     //G4cout << "function:    " << function << G4endl;
     
     
     if(range<function)
     {
     acceptance = true;
     }
     }
     
     phi = 2.0*M_PI*G4UniformRand();
     
     mx = sin(angDist_theta*deg)*cos(phi*deg);
     my = sin(angDist_theta*deg)*sin(phi*deg);
     mz = cos(angDist_theta*deg);
     
     fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
     */
    
    
    
    
    
    
    
    
    ////////////////////////////////////////////////////////////////
    ////    Range-finder for the acceptance-rejection method
    
    /*
     double highestValue = 0.0;
     
     for(G4double i=0; i<180.0; i+=0.1)
     {
     
     //angDist_theta = i/10000.0;
     //angDist_theta*=0.0174533;
     
     //function = 0.0;
     
     //function += pow( (2.0*0.203684)*(1.0/4.0)*sqrt(15.0/(2.0*M_PI))*pow(sin(angDist_theta), 2.0), 2.0);
     //function += pow( (0.592631)*(1.0/4.0)*sqrt(5.0/M_PI)*(3*cos(angDist_theta)*cos(angDist_theta) - 1), 2.0);
     //function *= sin(angDist_theta);
     
     
     angDist_theta = i;
     function = evaluateAngDist_interpolated(angDist_theta);
     
     //G4cout << "function, " << angDist_theta << ":    " << function << G4endl;
     
     
     if(function>highestValue) highestValue = function;
     }
     
     G4cout << "Here is the highestValue:    " << highestValue << G4endl;
     */
    
    
    
    ///////////////////////////////////////////////////
    //       Initial Energy Distribution of Particle
    ///////////////////////////////////////////////////
    
    //G4double InitialEnergy = G4RandGauss::shoot(((12.049-7.16192)*(0.75)), 0.012);
    //G4double InitialEnergy = ((12.049-7.16192)*(0.75)); // MeV
    
    //  Empirical
    //G4double InitialEnergy = G4RandGauss::shoot(2.96058e+00, ((1.5e-3)/2.3548));
    
    //fParticleGun->SetParticleEnergy(InitialEnergy*MeV);
    
    
    ////////////////////////
    ////    BiRelKin
    /*
     double m[4], T[4], E[4], p[4];
     double ThSCAT, ThSCAT_Recoil, Ex, Phi;
     
     for(G4int i=0; i<4; i++)
     {
     T[i] = 0.0;
     E[i] = 0.0;
     p[i] = 0.0;
     }
     
     ////	PR226: 16O(a, a')
     m[0] = 3.0160293191; // u
     m[1] = 9.0121822; // u
     m[2] = 3.0160492777; // u
     m[3] = 9.0133288; // u
     
     ////    Projectile Energy
     T[0] = 50.0; // MeV
     
     ////    Ejectile Energy
     T[1] = 0.0; // MeV
     
     
     ////    Uniformly distributed ThSCAT
     ThSCAT = G4RandFlat::shoot( -2., 2.); // deg
     */
    
    ////    Gaussian distributed ThSCAT
    //ThSCAT = 2.1;
    //while(pow(ThSCAT*ThSCAT, 0.5)>2.)
    //{
    //    ThSCAT = G4RandGauss::shoot(0., 1.0);
    //}
    //if(abs(ThSCAT)>2.) "while condition has been breached";
    
    
    ////    12.049 MeV 0+, 16O
    //Ex = G4RandGauss::shoot(12.049, (0.012/2.3548));
    
    ////    12.049 MeV 0+, 16O
    //Ex = G4RandGauss::shoot(15.097, (0.166/2.3548));
    
    
    
    ////////////////////////////////////////////
    ////    Selection of various states
    /*
     int decayPath = G4RandBit::shootBit();
     int decayProductNumber;
     int resonanceIndex;
     G4double Ex, decayProductEnergy;
     
     //G4ParticleDefinition* decayProduct_alpha = G4ParticleTable::GetParticleTable()->FindParticle("alpha");
     //G4ParticleDefinition* decayProduct_alpha = G4IonTable::GetIonTable()->GetIon(2,4,0.0);
     
     //G4ParticleDefinition* decayProduct_alpha = G4ParticleTable::GetParticleTable()->FindParticle("alpha");
     //G4ParticleDefinition* decayProduct_alpha = G4ParticleTable::GetParticleTable()->GetIon(2,4,0.0);
     
     G4ParticleDefinition* decayProduct;
     G4double usedEnergy;
     G4String decayModeName;
     
     if(decayPath==0)
     {
     ////////////////////////////
     ////    9B -> 4He + 5Li
     
     ////    Selecting the Decay Product of the Binary Decay
     decayProductNumber = G4RandBit::shootBit();
     
     if(decayProductNumber==0)
     {
     decayProduct = G4IonTable::GetIonTable()->GetIon(3,5,0.0);
     decayModeName = "9B->4He+5Li, particle=5Li";
     }
     if(decayProductNumber==1)
     {
     decayProduct = G4ParticleTable::GetParticleTable()->FindParticle("alpha");
     decayModeName = "9B->4He+5Li, particle=alpha";
     }
     
     
     ////    Override to use an alpha particle to approximate the energy loss for the 5Li particle
     //decayProduct = G4ParticleTable::GetParticleTable()->FindParticle("alpha");
     //decayProduct = G4IonTable::GetIonTable()->GetIon(3,5,0.0);
     fParticleGun->SetParticleDefinition(decayProduct);
     
     ////    Selecting the Resonance Energy
     resonanceIndex = (int) G4RandFlat::shoot( 0., noOfStates_9Be + 0.99999);
     
     
     Ex = CLHEP::RandBreitWigner::shoot(resonanceEnergy_9B[resonanceIndex], resonanceWidth_9B[resonanceIndex]);
     while(Ex>20.0)
     {
     Ex = CLHEP::RandBreitWigner::shoot(resonanceEnergy_9B[resonanceIndex], resonanceWidth_9B[resonanceIndex]);
     }
     
     
     //Ex = G4RandGauss::shoot(12.049, (0.012/2.3548));
     
     fEventAction->SetRecoilExcitationEnergy(Ex);
     
     ////    Setting the Decay Product Energy
     if(decayProductNumber==0) decayProductEnergy = (4.0/9.0)*Ex; // MeV
     if(decayProductNumber==1) decayProductEnergy = (5.0/9.0)*Ex; // MeV
     
     fParticleGun->SetParticleEnergy(decayProductEnergy*MeV);
     }
     
     
     if(decayPath==1)
     {
     ////////////////////////////
     ////    9B -> 8Be + p
     
     ////    Selecting the Decay Product of the Binary Decay
     decayProductNumber = G4RandBit::shootBit();
     
     if(decayProductNumber==0)
     {
     //decayProduct = G4IonTable::GetIonTable()->GetIon(1,1,0.0);
     decayProduct = G4ParticleTable::GetParticleTable()->FindParticle("proton");
     decayModeName = "9B->8Be+p, particle=proton";
     }
     if(decayProductNumber==1)
     {
     //decayProduct = G4IonTable::GetIonTable()->GetIon(4.0,8.0,0.0);
     decayProduct = G4ParticleTable::GetParticleTable()->FindParticle("alpha");
     decayModeName = "8Be->4He+4He, particle=proton";
     }
     
     
     ////    Override to use an alpha particle to approximate the energy loss for the 5Li particle
     //decayProduct = G4ParticleTable::GetParticleTable()->FindParticle("alpha");
     //decayProduct = G4IonTable::GetIonTable()->GetIon(3,5,0.0);
     fParticleGun->SetParticleDefinition(decayProduct);
     
     ////    Selecting the Resonance Energy
     resonanceIndex = (int) G4RandFlat::shoot( 0., noOfStates_9Be + 0.99999);
     
     
     Ex = CLHEP::RandBreitWigner::shoot(resonanceEnergy_9B[resonanceIndex], resonanceWidth_9B[resonanceIndex]);
     while(Ex>20.0)
     {
     Ex = CLHEP::RandBreitWigner::shoot(resonanceEnergy_9B[resonanceIndex], resonanceWidth_9B[resonanceIndex]);
     }
     
     fEventAction->SetRecoilExcitationEnergy(Ex);
     
     
     ////    Setting the Decay Product Energy
     if(decayProductNumber==0) decayProductEnergy = (8.0/9.0)*Ex; // MeV
     
     if(decayProductNumber==1)
     {
     decayProductEnergy = (1.0/9.0)*Ex; // MeV
     
     ////    Selecting a direction for the 8Be -> 4He + 4He decay
     mz = -1.0 + 2.0*G4UniformRand();
     a = sqrt(1-(mz*mz));
     mx = a*cos(theta);
     my = a*sin(theta);
     
     G4ThreeVector energy_COM = G4ThreeVector(mx, my, mz);
     energy_COM.setMag(0.120);
     
     ////    Selecting a direction for the kineitc energy of the 8Be decay product
     mz = -1.0 + 2.0*G4UniformRand();
     a = sqrt(1-(mz*mz));
     mx = a*cos(theta);
     my = a*sin(theta);
     
     G4ThreeVector energy_8Be = G4ThreeVector(mx, my, mz);
     energy_8Be.setMag(decayProductEnergy);
     
     ////    Calculating the possible energy of the
     G4ThreeVector energy_8Be_alpha = 0.5*energy_COM + energy_8Be;
     usedEnergy = energy_8Be_alpha.mag();
     
     fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
     
     }
     
     
     if(decayProductNumber==0) fParticleGun->SetParticleEnergy(decayProductEnergy*MeV);
     if(decayProductNumber==1) fParticleGun->SetParticleEnergy(usedEnergy*MeV);
     
     }
     
     
     
     fEventAction->SetDecayMode(decayModeName);
     
     fParticleGun->GeneratePrimaryVertex(anEvent);
     */
    
    
    ////    Inelastically scattered alpha
    /*
     decayProduct = G4ParticleTable::GetParticleTable()->FindParticle("alpha");
     //decayProduct = G4IonTable::GetIonTable()->GetIon(3,5,0.0);
     //decayProduct = G4IonTable::GetIonTable()->GetIon(2,4,0.0);
     fParticleGun->SetParticleDefinition(decayProduct);
     //fParticleGun->SetParticleEnergy(200.*MeV);
     */
    
    /*
     fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
     fParticleGun->GeneratePrimaryVertex(anEvent);
     */
    
    
    
    
    ////////////////////////
    ////    BiRelKin
    /*
     bool produceEjectile = true;
     bool produceRecoil = false;
     
     double m[4], T[4], E[4], p[4];
     double ThSCAT, ThSCAT_Recoil, Phi, Theta;
     //double Ex;
     
     for(G4int i=0; i<4; i++)
     {
     T[i] = 0.0;
     E[i] = 0.0;
     p[i] = 0.0;
     }
     
     ////	PR226: 16O(a, a')
     m[0] = 3.0160293191; // u
     m[1] = 9.0121822; // u
     m[2] = 3.0160492777; // u
     m[3] = 9.0133288; // u
     
     ////    Projectile Energy
     T[0] = 50.0; // MeV
     
     ////    Ejectile Energy
     T[1] = 0.0; // MeV
     
     
     ////    Uniformly distributed ThSCAT
     ThSCAT = G4RandFlat::shoot( 0., 2.); // deg
     
     
     ////    Gaussian distributed ThSCAT
     //ThSCAT = 2.1;
     //while(pow(ThSCAT*ThSCAT, 0.5)>2.)
     //{
     //    ThSCAT = G4RandGauss::shoot(0., 1.0);
     //}
     //if(abs(ThSCAT)>2.) "while condition has been breached";
     
     
     ////    12.049 MeV 0+, 16O
     //Ex = G4RandGauss::shoot(12.049, (0.012/2.3548));
     
     
     ////    Executing BiRelKin
     BiRelKin(m, T, E, p, ThSCAT, ThSCAT_Recoil, Ex);
     
     
     ////    Setting Recoil Energy
     if(produceEjectile)
     {
     fParticleGun->SetParticleEnergy(T[2]*MeV);
     Theta = ThSCAT;
     }
     
     if(produceRecoil)
     {
     fParticleGun->SetParticleEnergy(T[3]*MeV);
     Theta = ThSCAT_Recoil;
     }
     
     
     
     ////    Setting Recoil Momentum Direction
     Phi = G4RandFlat::shoot( 0., 360.);
     
     mx = sin(Theta*deg)*cos(Phi*deg);
     my = sin(Theta*deg)*sin(Phi*deg);
     mz = abs(cos(Theta*deg));
     fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
     fParticleGun->GeneratePrimaryVertex(anEvent);
     */
    
    
    
    
    //if(decayProductNumber==1) decayProduct = G4ParticleTable::GetParticleTable()->FindParticle("alpha");
    //fParticleGun->GeneratePrimaryVertex(anEvent);
    
    
    /*
     if(decayPath==0)
     {
     ////    Selecting the Decay Product of the Binary Decay
     decayProductNumber = (int) G4RandFlat::shoot( 0., 1 + 0.99999);
     
     ////    9B -> 4He +5Li
     if(decayProductNumber==0)
     {
     G4ParticleDefinition* decayProduct = G4ParticleTable::GetParticleTable()->FindParticle("alpha");
     fParticleGun->SetParticleDefinition(decayProduct);
     }
     
     if(decayProductNumber==1)
     {
     G4ParticleDefinition* decayProduct = G4IonTable::GetIonTable()->GetIon(3,5,0.0);
     fParticleGun->SetParticleDefinition(decayProduct);
     }
     
     ////    Selecting the Resonance Energy
     resonanceIndex = (int) G4RandFlat::shoot( 0., noOfStates_9Be + 0.99999);
     Ex = CLHEP::RandBreitWigner::shoot(resonanceEnergy_9B[resonanceIndex], resonanceWidth_9B[resonanceIndex]);
     fEventAction->SetRecoilExcitationEnergy(Ex);
     
     ////    Setting the Decay Product Energy
     if(decayProductNumber==0) decayProductEnergy = (5.0/4.0)*Ex; // MeV
     if(decayProductNumber==1) decayProductEnergy = (4.0/5.0)*Ex; // MeV
     
     fParticleGun->SetParticleEnergy(decayProductEnergy*MeV);
     }
     */
    
    
    
    
    
    
    
    /*
     G4double theta_cm, phi_cm;
     
     ////    Isotropic in Centre-Of-Mass Frame
     bool isotropy = true;
     
     
     theta_cm = G4RandFlat::shoot( 0., 360.); // deg
     //if()
     
     fParticleGun->SetParticleMomentumDirection(G4ThreeVector( mx, my, mz));
     */
    
    
    
    
    /*
     ////////    4He, +1 charge
     G4int Z = 2, A = 4;
     
     G4double ionCharge   = 1.*eplus;
     G4double excitEnergy = 0.*MeV;
     
     //G4ParticleDefinition* ion = G4ParticleTable::GetParticleTable()->GetIon(Z,A,excitEnergy);
     G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
     //ion->SetPDGLifeTime(1*ns);
     fParticleGun->SetParticleDefinition(ion);
     fParticleGun->SetParticleCharge(ionCharge);
     */
    
    
    ///////////////////////////////////////////////////
    //              TIME DISTRIBUTION                //
    ///////////////////////////////////////////////////
    
    /*
     initialclocktime = G4RandGauss::shoot( 0, 1.5);
     fParticleGun->SetParticleTime(initialclocktime*ns);
     */
    
    
    
    
    ///////////////////////////////////////////////////
    //       Initial Direction of Particle
    ///////////////////////////////////////////////////
    
    //      Flat Distribution
    //mx = G4RandFlat::shoot( -.15, .15);
    //my = G4RandFlat::shoot( -.15, .15);
    //mz = G4RandFlat::shoot( -1, 1);
    
    ////////////////////////////////
    ////    ISOTROPIC
    //      Exploiting the spherical symmetry of normal distributions
    //      Gaussian Distribution
    /*
     mx = G4RandGauss::shoot( 0, 1.);
     my = G4RandGauss::shoot( 0, 1.);
     mz = G4RandGauss::shoot( 0, 1.);
     
     if(mz<0) mz = -mz;
     //mz = abs(mz);
     //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,-1.));
     //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, 10.));
     fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, -mz));
     */
    
    ////    Alternative Method for ISOTROPY
    /*
     G4double theta = 2*M_PI*G4UniformRand();
     //G4double mz = -1.0 + 2*G4UniformRand();
     G4double mz = -1.0 + G4UniformRand();
     
     G4double a = sqrt(1-(mz*mz));
     
     mx = a*cos(theta);
     my = a*sin(theta);
     
     fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
     */
    
    
    /*
     theta = M_PI*G4UniformRand();
     phi = 2*M_PI*G4UniformRand();
     fParticleGun->SetParticleMomentumDirection(G4ThreeVector(sin(theta*rad)*cos(phi*rad), sin(theta*rad)*sin(phi*rad), cos(theta*rad)));
     */
    
    
    
    ///////////////////////////////////////////////////
    //                  16O(a,a')
    ///////////////////////////////////////////////////
    /*
     G4double test = G4UniformRand();
     
     //if(test<=0.5) beamEnergy = G4RandGauss::shoot( 200, 0.5); // MeV
     //if(test>0.5) beamEnergy = G4RandGauss::shoot( 180, 0.5); // MeV
     
     if(test>0.0 && test<=0.1) beamEnergy = G4RandGauss::shoot( 215, 0.5); // MeV
     if(test>0.1 && test<=0.2) beamEnergy = G4RandGauss::shoot( 210, 0.5); // MeV
     if(test>0.2 && test<=0.3) beamEnergy = G4RandGauss::shoot( 200, 0.5); // MeV
     if(test>0.3 && test<=0.4) beamEnergy = G4RandGauss::shoot( 190, 0.5); // MeV
     if(test>0.4 && test<=0.5) beamEnergy = G4RandGauss::shoot( 180, 0.5); // MeV
     if(test>0.5 && test<=0.6) beamEnergy = G4RandGauss::shoot( 185, 0.5); // MeV
     if(test>0.6 && test<=0.7) beamEnergy = G4RandGauss::shoot( 170, 0.5); // MeV
     if(test>0.7 && test<=0.8) beamEnergy = G4RandGauss::shoot( 160, 0.5); // MeV
     if(test>0.8 && test<=0.9) beamEnergy = G4RandGauss::shoot( 150, 0.5); // MeV
     if(test>0.9 && test<=1.0) beamEnergy = G4RandGauss::shoot( 140, 0.5); // MeV
     
     //beamEnergy = G4RandGauss::shoot( 200, 0.5); // MeV
     //beamEnergy = 200; // MeV
     
     recoilExcitationEnergy = 15.097;  // MeV
     alphaSeperationEnergy = 7.16192;    // MeV
     protonSeperationEnergy = 12.1274;  // MeV
     //  Elastically scattered outgoing alpha
     mx = G4RandGauss::shoot( 0, 0.00001);
     my = G4RandGauss::shoot( 0, .05);
     //mx = G4RandGauss::shoot( 0, 0.00001);
     //my = G4RandGauss::shoot( 0, 0.00001);
     //mx = 0;
     //my = 0;
     
     fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, 10.));
     energy = beamEnergy - recoilExcitationEnergy; // MeV
     fParticleGun->SetParticleEnergy(energy*MeV);
     fParticleGun->GeneratePrimaryVertex(anEvent);
     
     
     //////////////////////////////////////////////////////////////
     //                      16O -> 4He + 12C
     //////////////////////////////////////////////////////////////
     //  Alpha 0 isotropic decay:    15.097 MeV state of 16O -> to Ground state of 12C
     //  Alpha 1 anisotropic decay:  15.097 MeV state of 16O -> to 4.43891 MeV state of 12C
     //  Assuming recoil nucleus takes no kinetic energy from inelastic scatter
     
     G4double pAlphaDecayMode = G4UniformRand();
     
     if(pAlphaDecayMode<=0.5) alphaDecayMode = 0;
     if(pAlphaDecayMode>0.5) alphaDecayMode = 1;
     
     //alphaDecayMode = 1;
     
     ////    ALPHA 1
     if(alphaDecayMode == 0) daughterExcitationEnergy = 0; //    MeV
     
     ////    Alpha 0
     if(alphaDecayMode == 1) daughterExcitationEnergy = 4.43891; //    MeV
     
     mx = G4RandFlat::shoot( -1., 1.);
     my = G4RandFlat::shoot( -1., 1.);
     mz = G4RandFlat::shoot( -1., 1.);
     
     fParticleGun->SetParticleEnergy( (recoilExcitationEnergy - alphaSeperationEnergy - daughterExcitationEnergy)*(3/4)*MeV);
     
     fParticleGun->SetParticleMomentumDirection(ejectileDirection);
     
     fParticleGun->GeneratePrimaryVertex(anEvent);
     */
    
    
    ////    Recoil - 12C
    /*
     fParticleGun->SetParticleEnergy( (recoilExcitationEnergy - alphaSeperationEnergy - daughterExcitationEnergy)*(1/4)*MeV);
     
     recoilDirection = -ejectileDirection;
     fParticleGun->SetParticleMomentumDirection(recoilDirection);
     */
    
    /*
     Z = 6, A = 12;
     
     G4ParticleDefinition* recoil = G4IonTable::GetIonTable()->GetIon(Z, A, 4.43891*MeV);
     //G4ParticleDefinition* recoil = G4ParticleTable::GetParticleTable()->GetIon(Z, A, 0.*keV);
     fParticleGun->SetParticleCharge(0.*eplus);
     //recoil->SetPDGLifeTime(1*ns);
     //recoil->SetPDGStable(true);
     fParticleGun->SetParticleDefinition(recoil);
     
     fParticleGun->SetParticleEnergy( (recoilExcitationEnergy - alphaSeperationEnergy - daughterExcitationEnergy)*(1/4)*MeV);
     
     recoilDirection = -ejectileDirection;
     fParticleGun->SetParticleMomentumDirection(recoilDirection);
     */
    
    
    //fParticleGun->GeneratePrimaryVertex(anEvent);
    
    //G4double Lifetime = fParticleGun->GetParticleDefinition()->GetPDGLifeTime();
    //G4double Lifetime = recoil->GetPDGLifeTime();
    
    //G4cout << "Here is the alphaSeperationEnergy    "<< alphaSeperationEnergy << G4endl;
    
    
    
    ////     Gamma Decay - from Recoil Nucleus
    /*
     if(alphaDecayMode == 1)
     {
     G4ParticleDefinition* particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
     fParticleGun->SetParticleDefinition(particleDefinition);
     
     fParticleGun->SetParticleEnergy(daughterExcitationEnergy*MeV);
     
     mx = G4RandFlat::shoot( -1., 1.);
     my = G4RandFlat::shoot( -1., 1.);
     mz = G4RandFlat::shoot( -1., 1.);
     fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
     
     fParticleGun->GeneratePrimaryVertex(anEvent);
     }
     */
    
    
    
    
    
    
    
    ////////////////////////////////////////////////////////
    //          PARTICLE INFORMATION OUTPUT
    ////////////////////////////////////////////////////////
    /*
     G4double Lifetime = fParticleGun->GetParticleDefinition()->GetPDGLifeTime();
     G4double Spin = fParticleGun->GetParticleDefinition()->GetPDGSpin();
     G4double Isospin = fParticleGun->GetParticleDefinition()->GetPDGIsospin();
     G4double Parity = fParticleGun->GetParticleDefinition()->GetPDGiGParity();
     G4String ParticleName = fParticleGun->GetParticleDefinition()->GetParticleName();
     
     
     
     G4cout << "Here is the Lifetime    "<< Lifetime << G4endl;
     G4cout << "Here is the Spin    "<< Spin << G4endl;
     G4cout << "Here is the Isospin    "<< Isospin << G4endl;
     G4cout << "Here is the Parity    "<< Parity << G4endl;
     G4cout << "Here is the ParticleName    "<< ParticleName << G4endl;
     */
    
    
    //create vertex
    //
    //fParticleGun->GeneratePrimaryVertex(anEvent);
    
    
    
    /*
     // This function is called at the begining of event
     
     // In order to avoid dependence of PrimaryGeneratorAction
     // on DetectorConstruction class we get world volume
     // from G4LogicalVolumeStore
     //
     G4double worldZHalfLength = 0;
     G4LogicalVolume* worlLV
     = G4LogicalVolumeStore::GetInstance()->GetVolume("World");
     G4Box* worldBox = 0;
     if ( worlLV) worldBox = dynamic_cast< G4Box*>(worlLV->GetSolid());
     if ( worldBox ) {
     worldZHalfLength = worldBox->GetZHalfLength();
     }
     else  {
     G4ExceptionDescription msg;
     msg << "World volume of box not found." << G4endl;
     msg << "Perhaps you have changed geometry." << G4endl;
     msg << "The gun will be place in the center.";
     G4Exception("PrimaryGeneratorAction::GeneratePrimaries()",
     "MyCode0002", JustWarning, msg);
     }
     
     // Set gun position
     fParticleGun
     ->SetParticlePosition(G4ThreeVector(0., 0., -worldZHalfLength));
     
     fParticleGun->GeneratePrimaryVertex(anEvent);
     */
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::CalculateBinaryRelativisticKinematics(long double theta_ejectile, long double phi_ejectile, long double &theta_recoil_LAB, long double &phi_recoil_LAB, long double &eX)
{
    for(G4int i=0; i<4; i++)
    {
        T[i] = 0.0;
        E[i] = 0.0;
        p[i] = 0.0;
    }
    
//    ////	PR226: 16O(a, a')
//    m[0] = 1.00782503224; // u
//    m[1] = 14.003241989; // u
//    m[2] = 3.01604928199; // u
//    m[3] = 12.0; // u
//    
//    ////    Projectile Energy
//    T[0] = 100.0; // MeV
//    T[0] = G4RandGauss::shoot(100.0, (0.010/2.3548));
//
//    ////    Excitation Energy
//    eX = G4RandGauss::shoot(12.049, (0.012/2.3548));

    //------------------------------------------------
    ////	X17
    m[0] = 7.01600342665;
    m[1] = 1.00782503223;
    m[2] = 16.7/931.49410242;
    m[3] = 8.005305102;
    
    ////    Projectile/Beam Kinetic Energy (Lab Frame)
    T[0] = 3.07003;

    ////    Target Kinetic Energy (Lab Frame)
    T[1] = 0.0;
    
    ////    Excitation Energy
    eX = 0.0;

    //------------------------------------------------
    ////    Executing BiRelKin
    BiRelKin(m, T, E, p, theta_ejectile, theta_recoil_LAB, eX);
    
    ////    Setting Recoil Momentum Direction
    phi_recoil_LAB = phi_ejectile - 180.0; // deg
    
    if(phi_recoil_LAB<0.0)
    {
        phi_recoil_LAB += 360.0; // deg
    }
    
    /*
     mx = sin(theta*deg)*cos(phi*deg);
     my = sin(theta*deg)*sin(phi*deg);
     mz = abs(cos(theta*deg));
     */
    //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::CalculateBinaryDecayKinematics(double eX, double eX_residual, double sepE, double initialScatteringReactionMomentum, double initialScatteringReactionEnergy, std::vector<double> decay_A, double theta_recoil_LAB, double phi_recoil_LAB, std::vector<double>& decay_theta_recoilCOM, std::vector<double>& decay_phi_recoilCOM, std::vector<double>& decay_theta_LAB, std::vector<double>& decay_phi_LAB, std::vector<double>& decay_T_LAB)
{
    if(eX - eX_residual - sepE >= 0.0) // Energy conservation
    {
        //------------------------------------------------
        G4ThreeVector recoilCOM_P_v3 = initialScatteringReactionMomentum*G4ThreeVector(sin(theta_recoil_LAB*deg)*cos(phi_recoil_LAB*deg), sin(theta_recoil_LAB*deg)*sin(phi_recoil_LAB*deg), cos(theta_recoil_LAB*deg));
        G4LorentzVector recoilCOM_P_v4 = G4LorentzVector(recoilCOM_P_v3, initialScatteringReactionEnergy/pow(c2, 0.5));
        double betaRecoil = recoilCOM_P_v4.beta();
//        std::cout << "betaRecoil: " << betaRecoil << std::endl;
        //        betaRecoil = 0.0;
        G4ThreeVector boostVector_recoilCOMtoLAB = betaRecoil*(recoilCOM_P_v3.unit());
        
        //------------------------------------------------
        double sumOfIndividualMasses = decay_A[0] + decay_A[1]; // amu
        //decay_T[0] = (eX - sepE - eX_residual)*(decay_A[1]/sumOfIndividualMasses); // MeV
        //decay_T[1] = (eX - sepE - eX_residual)*(decay_A[0]/sumOfIndividualMasses); // MeV
        double decay_T_COM_0 = (eX - sepE - eX_residual)*(decay_A[1]/sumOfIndividualMasses); // MeV
        double decay_T_COM_1 = (eX - sepE - eX_residual)*(decay_A[0]/sumOfIndividualMasses); // MeV
        
//        std::cout << "decay_T_COM_0: " << decay_T_COM_0 << std::endl;
//        std::cout << "decay_T_COM_1: " << decay_T_COM_1 << std::endl;
        
        double decay_0_E = decay_T_COM_0 + (decay_A[0]*c2);
        double decay_0_P = (1.0/pow(c2, 0.5))*sqrt(pow(decay_0_E, 2.0) - pow(decay_A[0]*c2, 2.0));
        
        //        std::cout << "(decay_A[0]*c2): " << (decay_A[0]*c2) << std::endl;
        
        //        std::cout << "decay_0_E: " << decay_0_E << std::endl;
        //        std::cout << "decay_0_P: " << decay_0_P << std::endl;
        
        double decay_1_E = decay_T_COM_1 + (decay_A[1]*c2);
        double decay_1_P = (1.0/pow(c2, 0.5))*sqrt(pow(decay_1_E, 2.0) - pow(decay_A[1]*c2, 2.0));
        
        //        std::cout << "decay_1_E: " << decay_1_E << std::endl;
        //        std::cout << "decay_1_P: " << decay_1_P << std::endl;
        
        //------------------------------------------------
        G4ThreeVector decay_0_P_v3 = decay_0_P*G4ThreeVector(sin(decay_theta_recoilCOM[0]*deg)*cos(decay_phi_recoilCOM[0]*deg), sin(decay_theta_recoilCOM[0]*deg)*sin(decay_phi_recoilCOM[0]*deg), cos(decay_theta_recoilCOM[0]*deg));
        G4LorentzVector decay_0_P_v4 = G4LorentzVector(decay_0_P_v3, decay_0_E/pow(c2, 0.5));
        
        G4ThreeVector decay_1_P_v3 = decay_1_P*G4ThreeVector(sin(decay_theta_recoilCOM[1]*deg)*cos(decay_phi_recoilCOM[1]*deg), sin(decay_theta_recoilCOM[1]*deg)*sin(decay_phi_recoilCOM[1]*deg), cos(decay_theta_recoilCOM[1]*deg));
        G4LorentzVector decay_1_P_v4 = G4LorentzVector(decay_1_P_v3, decay_1_E/pow(c2, 0.5));
        
        //------------------------------------------------
        decay_0_P_v4.boost(boostVector_recoilCOMtoLAB);
        decay_theta_LAB[0] = decay_0_P_v4.vect().theta()/deg;
        
        decay_1_P_v4.boost(boostVector_recoilCOMtoLAB);
        decay_theta_LAB[1] = decay_1_P_v4.vect().theta()/deg;
        
        //------------------------------------------------
        decay_phi_LAB[0] = decay_0_P_v4.vect().phi()/deg;
        
        while(decay_phi_LAB[0]<0.0)
        {
            decay_phi_LAB[0] += 360.0;
        }
        
        decay_phi_LAB[1] = decay_phi_LAB[0] - 180.0;
        
        while(decay_phi_LAB[1]<0.0)
        {
            decay_phi_LAB[1] += 360.0;
        }
        
        //------------------------------------------------
        decay_T_LAB[0] = decay_0_P_v4.vect().mag2()/(2.0*decay_A[0]); // MeV
        decay_T_LAB[1] = decay_1_P_v4.vect().mag2()/(2.0*decay_A[1]); // MeV
        
//        std::cout << "decay_0_P_v4.vect().mag2(): " << decay_0_P_v4.vect().mag2() << std::endl;
//        std::cout << "2.0*decay_A[0]: " << 2.0*decay_A[0] << std::endl;
//        
//        std::cout << "decay_T_LAB[0]: " << decay_T_LAB[0] << std::endl;
//        std::cout << "decay_T_LAB[1]: " << decay_T_LAB[1] << std::endl;
        
        
        
        //        decay_T_LAB[0] = decay_0_P_v4.vect().mag2()/(2.0*decay_A[0]*c2); // MeV
        //        decay_T_LAB[1] = decay_1_P_v4.vect().mag2()/(2.0*decay_A[1]*c2); // MeV
        //        std::cout << "decay_T_LAB[0]: " << decay_T_LAB[0] << std::endl;
        //        std::cout << "decay_T_LAB[1]: " << decay_T_LAB[1] << std::endl;
        
        //        std::cout << std::endl;
    }
}

void PrimaryGeneratorAction::CalculateBinaryDecayKinematics2(long double eX, long double eX_residual, long double sepE, long double initialScatteringReactionMomentum, long double initialScatteringReactionEnergy, std::vector<long double> decay_A, long double theta_recoil_LAB, long double phi_recoil_LAB, std::vector<long double>& decay_theta_recoilCOM, std::vector<long double>& decay_phi_recoilCOM, std::vector<long double>& decay_theta_LAB, std::vector<long double>& decay_phi_LAB, std::vector<long double>& decay_T_LAB)
{
    if(eX - eX_residual - sepE >= 0.0) // Energy conservation
    {
        //------------------------------------------------
        G4ThreeVector recoilCOM_P_v3 = initialScatteringReactionMomentum*G4ThreeVector(sin(theta_recoil_LAB*deg)*cos(phi_recoil_LAB*deg), sin(theta_recoil_LAB*deg)*sin(phi_recoil_LAB*deg), cos(theta_recoil_LAB*deg));
        G4LorentzVector recoilCOM_P_v4 = G4LorentzVector(recoilCOM_P_v3, initialScatteringReactionEnergy/pow(c2, 0.5));
        long double betaRecoil = recoilCOM_P_v4.beta();
        std::cout << "betaRecoil: " << betaRecoil << std::endl;
        //        betaRecoil = 0.0;
        G4ThreeVector boostVector_recoilCOMtoLAB = betaRecoil*(recoilCOM_P_v3.unit());
        
        //------------------------------------------------
        long double sumOfIndividualMasses = decay_A[0] + decay_A[1]; // amu
        //decay_T[0] = (eX - sepE - eX_residual)*(decay_A[1]/sumOfIndividualMasses); // MeV
        //decay_T[1] = (eX - sepE - eX_residual)*(decay_A[0]/sumOfIndividualMasses); // MeV
        long double decay_T_COM_0 = (eX - sepE - eX_residual)*(decay_A[1]/sumOfIndividualMasses); // MeV
        long double decay_T_COM_1 = (eX - sepE - eX_residual)*(decay_A[0]/sumOfIndividualMasses); // MeV
        
        std::cout << "decay_T_COM_0: " << decay_T_COM_0 << std::endl;
        std::cout << "decay_T_COM_1: " << decay_T_COM_1 << std::endl;
        
        long double decay_0_E = decay_T_COM_0 + (decay_A[0]*c2);
        long double decay_0_P = (1.0/std::pow(c2, 0.5))*std::sqrt(std::pow(decay_0_E, 2.0) - std::pow(decay_A[0]*c2, 2.0));
        
        //        std::cout << "(decay_A[0]*c2): " << (decay_A[0]*c2) << std::endl;
        
        //        std::cout << "decay_0_E: " << decay_0_E << std::endl;
        //        std::cout << "decay_0_P: " << decay_0_P << std::endl;
        
        long double decay_1_E = decay_T_COM_1 + (decay_A[1]*c2);
        long double decay_1_P = (1.0/std::pow(c2, 0.5))*std::sqrt(std::pow(decay_1_E, 2.0) - std::pow(decay_A[1]*c2, 2.0));
        
        //        std::cout << "decay_1_E: " << decay_1_E << std::endl;
        //        std::cout << "decay_1_P: " << decay_1_P << std::endl;
        
        //------------------------------------------------
        G4ThreeVector decay_0_P_v3 = decay_0_P*G4ThreeVector(sin(decay_theta_recoilCOM[0]*deg)*cos(decay_phi_recoilCOM[0]*deg), sin(decay_theta_recoilCOM[0]*deg)*sin(decay_phi_recoilCOM[0]*deg), cos(decay_theta_recoilCOM[0]*deg));
        G4LorentzVector decay_0_P_v4 = G4LorentzVector(decay_0_P_v3, decay_0_E/std::pow(c2, 0.5));
        
        G4ThreeVector decay_1_P_v3 = decay_1_P*G4ThreeVector(sin(decay_theta_recoilCOM[1]*deg)*cos(decay_phi_recoilCOM[1]*deg), sin(decay_theta_recoilCOM[1]*deg)*sin(decay_phi_recoilCOM[1]*deg), cos(decay_theta_recoilCOM[1]*deg));
        G4LorentzVector decay_1_P_v4 = G4LorentzVector(decay_1_P_v3, decay_1_E/std::pow(c2, 0.5));
        
        //------------------------------------------------
        decay_0_P_v4.boost(boostVector_recoilCOMtoLAB);
        decay_theta_LAB[0] = decay_0_P_v4.vect().theta()/deg;
        
        decay_1_P_v4.boost(boostVector_recoilCOMtoLAB);
        decay_theta_LAB[1] = decay_1_P_v4.vect().theta()/deg;
        
        //------------------------------------------------
        decay_phi_LAB[0] = decay_0_P_v4.vect().phi()/deg;
        
        while(decay_phi_LAB[0]<0.0)
        {
            decay_phi_LAB[0] += 360.0;
        }
        
        decay_phi_LAB[1] = decay_phi_LAB[0] - 180.0;
        
        while(decay_phi_LAB[1]<0.0)
        {
            decay_phi_LAB[1] += 360.0;
        }
        
        //------------------------------------------------
        decay_T_LAB[0] = decay_0_P_v4.vect().mag2()/(2.0*decay_A[0]); // MeV
        decay_T_LAB[1] = decay_1_P_v4.vect().mag2()/(2.0*decay_A[1]); // MeV
        
        std::cout << "decay_0_P_v4.vect().mag2(): " << decay_0_P_v4.vect().mag2() << std::endl;
        std::cout << "2.0*decay_A[0]: " << 2.0*decay_A[0] << std::endl;
        
        std::cout << "decay_T_LAB[0]: " << decay_T_LAB[0] << std::endl;
        std::cout << "decay_T_LAB[1]: " << decay_T_LAB[1] << std::endl;
        
        
        
        //        decay_T_LAB[0] = decay_0_P_v4.vect().mag2()/(2.0*decay_A[0]*c2); // MeV
        //        decay_T_LAB[1] = decay_1_P_v4.vect().mag2()/(2.0*decay_A[1]*c2); // MeV
        //        std::cout << "decay_T_LAB[0]: " << decay_T_LAB[0] << std::endl;
        //        std::cout << "decay_T_LAB[1]: " << decay_T_LAB[1] << std::endl;
        
        //        std::cout << std::endl;
    }
}

std::vector<G4LorentzVector> PrimaryGeneratorAction::CalculateBinaryDecayKinematics3(double eX, double eX_residual, double sepE, double initialScatteringReactionMomentum, double initialScatteringReactionEnergy, std::vector<double> decay_A, double theta_recoil_LAB, double phi_recoil_LAB, std::vector<double>& decay_theta_recoilCOM, std::vector<double>& decay_phi_recoilCOM, std::vector<double>& decay_theta_LAB, std::vector<double>& decay_phi_LAB, std::vector<double>& decay_T_LAB)
{
    //Function to calculate the recoil G4LorentzVector in the lab frame
    std::vector<G4LorentzVector> result;
    
//    double RecoilMass = Masses[3] + recoilEx; //Convert the recoil mass to the right value including the excitation-energy dependence
//    
//    double s = std::pow(Masses[0],2.) + std::pow(Masses[1],2.) + 2 * Masses[1] * (TBeam + Masses[0]);
//    
//    G4LorentzVector f4MomentumCentreOfMass(0,0,0,std::sqrt(s));
//    
//    double ECM0 = (s + std::pow(Masses[0],2.) - std::pow(Masses[1],2.))/(2*std::sqrt(s));
//    double ECM1 = (s + std::pow(Masses[1],2.) - std::pow(Masses[0],2.))/(2*std::sqrt(s));
//    double ECM2 = (s + std::pow(Masses[2],2.) - std::pow(RecoilMass,2.))/(2*std::sqrt(s));
//    double ECM3 = (s + std::pow(RecoilMass,2.) - std::pow(Masses[2],2.))/(2*std::sqrt(s));
//    
//    //    std::cout << "ECM2: " << ECM2 << std::endl;
//    //    std::cout << "ECM3: " << ECM3 << std::endl;
//    
//    double PCM0 = std::sqrt(std::pow(ECM0,2.) - std::pow(Masses[0],2.));
//    double PCM1 = std::sqrt(std::pow(ECM1,2.) - std::pow(Masses[1],2.));
//    double PCM2 = std::sqrt(std::pow(ECM2,2.) - std::pow(Masses[2],2.));
//    double PCM3 = std::sqrt(std::pow(ECM3,2.) - std::pow(RecoilMass,2.));
//    
//    G4ThreeVector Lab3Momentum0(0,0,std::sqrt(std::pow(TBeam,2.) + 2 * TBeam * Masses[0]));
//    G4ThreeVector Lab3Momentum1(0,0,0);
//    
//    G4LorentzVector Lab4Momentum0(Lab3Momentum0,Masses[0] + TBeam);
//    G4LorentzVector Lab4Momentum1(Lab3Momentum1,Masses[1]);
//    
//    double BetaCM = (Lab4Momentum0+Lab4Momentum1).beta();
//    //std::cout << "BetaCM = " << BetaCM << std::endl;
//    
//    G4LorentzVector CoM4Momentum0 = Lab4Momentum0;
//    CoM4Momentum0.boost(G4ThreeVector(0,0,-BetaCM));
//    
//    G4LorentzVector CoM4Momentum1 = Lab4Momentum1;
//    CoM4Momentum1.boost(G4ThreeVector(0,0,-BetaCM));
//    
//    G4LorentzVector CoM4Momentum2 = G4LorentzVector(PCM2 * std::sin(ThetaAlphaCM * M_PI/180.) * std::cos(PhiAlphaCM * M_PI/180.),
//                                                    PCM2 * std::sin(ThetaAlphaCM * M_PI/180.) * std::sin(PhiAlphaCM * M_PI/180.),
//                                                    PCM2 * std::cos(ThetaAlphaCM * M_PI/180.),
//                                                    ECM2);
//    G4LorentzVector CoM4Momentum3 = f4MomentumCentreOfMass - CoM4Momentum2;
//    
//    G4LorentzVector Lab4Momentum2 = CoM4Momentum2;
//    Lab4Momentum2.boost(0,0,BetaCM);
//    G4LorentzVector Lab4Momentum3 = CoM4Momentum3;
//    Lab4Momentum3.boost(0,0,BetaCM);
//    
//    //if(VerboseFlag)Lab4Momentum2.Print();
//    //if(VerboseFlag)Lab4Momentum3.Print();
//    
//    if((Lab4Momentum0+Lab4Momentum1-Lab4Momentum2-Lab4Momentum3).mag()>1.e-6)
//    {
//        std::cout << "Problem with two-body kinematics" << std::endl;
//        //        (Lab4Momentum0+Lab4Momentum1-Lab4Momentum2-Lab4Momentum3).Print();
//    }
//    
//    result.push_back(Lab4Momentum0);
//    result.push_back(Lab4Momentum1);
//    result.push_back(Lab4Momentum2);
//    result.push_back(Lab4Momentum3);
    
    return result;
}

std::vector<G4LorentzVector> PrimaryGeneratorAction::CalculateReactionMomentumVectors(double *massEnergies, double TBeam, double ThetaAlphaCM, double PhiAlphaCM)
{
    //Function to calculate the recoil G4LorentzVector in the lab frame
    std::vector<G4LorentzVector> result;
    
//    double RecoilMass = Masses[3] + recoilEx; //Convert the recoil mass to the right value including the excitation-energy dependence
    
    double s = std::pow(massEnergies[0],2.) + std::pow(massEnergies[1],2.) + 2 * massEnergies[1] * (TBeam + massEnergies[0]);
    
    G4LorentzVector f4MomentumCentreOfMass(0,0,0,std::sqrt(s));
    
    double ECM0 = (s + std::pow(massEnergies[0],2.) - std::pow(massEnergies[1],2.))/(2*std::sqrt(s));
    double ECM1 = (s + std::pow(massEnergies[1],2.) - std::pow(massEnergies[0],2.))/(2*std::sqrt(s));
    double ECM2 = (s + std::pow(massEnergies[2],2.) - std::pow(massEnergies[3],2.))/(2*std::sqrt(s));
    double ECM3 = (s + std::pow(massEnergies[3],2.) - std::pow(massEnergies[2],2.))/(2*std::sqrt(s));
    
//    std::cout << "ECM2: " << ECM2 << std::endl;
//    std::cout << "ECM3: " << ECM3 << std::endl;
    
    double PCM0 = std::sqrt(std::pow(ECM0,2.) - std::pow(massEnergies[0],2.));
    double PCM1 = std::sqrt(std::pow(ECM1,2.) - std::pow(massEnergies[1],2.));
    double PCM2 = std::sqrt(std::pow(ECM2,2.) - std::pow(massEnergies[2],2.));
    double PCM3 = std::sqrt(std::pow(ECM3,2.) - std::pow(massEnergies[3],2.));
    
    G4ThreeVector Lab3Momentum0(0,0,std::sqrt(std::pow(TBeam,2.) + 2 * TBeam * massEnergies[0]));
    G4ThreeVector Lab3Momentum1(0,0,0);
    
    G4LorentzVector Lab4Momentum0(Lab3Momentum0,massEnergies[0] + TBeam);
    G4LorentzVector Lab4Momentum1(Lab3Momentum1,massEnergies[1]);
    
    double BetaCM = (Lab4Momentum0+Lab4Momentum1).beta();
    //std::cout << "BetaCM = " << BetaCM << std::endl;
    
    G4LorentzVector CoM4Momentum0 = Lab4Momentum0;
    CoM4Momentum0.boost(G4ThreeVector(0,0,-BetaCM));
    
    G4LorentzVector CoM4Momentum1 = Lab4Momentum1;
    CoM4Momentum1.boost(G4ThreeVector(0,0,-BetaCM));
    
    G4LorentzVector CoM4Momentum2 = G4LorentzVector(PCM2 * std::sin(ThetaAlphaCM * M_PI/180.) * std::cos(PhiAlphaCM * M_PI/180.),
                                                  PCM2 * std::sin(ThetaAlphaCM * M_PI/180.) * std::sin(PhiAlphaCM * M_PI/180.),
                                                  PCM2 * std::cos(ThetaAlphaCM * M_PI/180.),
                                                  ECM2);
    G4LorentzVector CoM4Momentum3 = f4MomentumCentreOfMass - CoM4Momentum2;
    
    G4LorentzVector Lab4Momentum2 = CoM4Momentum2;
    Lab4Momentum2.boost(0,0,BetaCM);
    G4LorentzVector Lab4Momentum3 = CoM4Momentum3;
    Lab4Momentum3.boost(0,0,BetaCM);
    
    //if(VerboseFlag)Lab4Momentum2.Print();
    //if(VerboseFlag)Lab4Momentum3.Print();
    
    if((Lab4Momentum0+Lab4Momentum1-Lab4Momentum2-Lab4Momentum3).mag()>1.e-6)
    {
        std::cout << "Problem with two-body kinematics" << std::endl;
//        (Lab4Momentum0+Lab4Momentum1-Lab4Momentum2-Lab4Momentum3).Print();
    }
    
    result.push_back(Lab4Momentum0);
    result.push_back(Lab4Momentum1);
    result.push_back(Lab4Momentum2);
    result.push_back(Lab4Momentum3);
    
    return result;

}

//std::vector<TLorentzVector> PrimaryGeneratorAction::CalculateBinaryDecayMomentumVectors(double *massEnergies, TLorentzVector parentMomentumVector, double thetaCM, double phiCM)
//{
//    //parentParticleKineticEnergy
//    
//    //Function to calculate the recoil TLorentzVector in the lab frame
//    std::vector<TLorentzVector> result;
//    
//    double s = std::pow(massEnergies[0],2.);
//    TLorentzVector f4MomentumCentreOfMass(0,0,0,std::sqrt(s));
//    
//    //    double BetaCM = (f4MomentumCentreOfMass).Beta();
//    
//    double ECM1 = (std::pow(massEnergies[0],2.) + std::pow(massEnergies[1],2.) - std::pow(massEnergies[2],2.))/(2*massEnergies[0]);
//    double ECM2 = (std::pow(massEnergies[0],2.) - std::pow(massEnergies[1],2.) + std::pow(massEnergies[2],2.))/(2*massEnergies[0]);
//    
//    std::cout << "ECM1: " << ECM1 << std::endl;
//    std::cout << "ECM2: " << ECM2 << std::endl;
//    
//    double PCM1 = std::sqrt(std::pow(ECM1,2.) - std::pow(massEnergies[1],2.));
//    double PCM2 = std::sqrt(std::pow(ECM2,2.) - std::pow(massEnergies[2],2.));
//    
//    //    TVector3 Lab3Momentum0(0,0,std::sqrt(std::pow(parentParticleKineticEnergy,2.) + 2*parentParticleKineticEnergy*massEnergies[0]));
//    //    TLorentzVector Lab4Momentum0(Lab3Momentum0,massEnergies[0] + parentParticleKineticEnergy);
//    
//    double BetaCM = (parentMomentumVector).Beta();
//    
//    //    TLorentzVector CoM4Momentum0 = Lab4Momentum0;
//    //    CoM4Momentum0.Boost(TVector3(0,0,-BetaCM));
//    
//    TLorentzVector CoM4Momentum1 = TLorentzVector(PCM1 * std::sin(thetaCM * M_PI/180.) * std::cos(phiCM * M_PI/180.),
//                                                  PCM1 * std::sin(thetaCM * M_PI/180.) * std::sin(phiCM * M_PI/180.),
//                                                  PCM1 * std::cos(thetaCM * M_PI/180.),
//                                                  ECM1);
//    TLorentzVector CoM4Momentum2 = f4MomentumCentreOfMass - CoM4Momentum1;
//    
//    //------------------------------------------------
//    TLorentzVector Lab4Momentum1 = CoM4Momentum1;
//    Lab4Momentum1.Boost(0,0,BetaCM);
//    TLorentzVector Lab4Momentum2 = CoM4Momentum2;
//    Lab4Momentum2.Boost(0,0,BetaCM);
//    
//    //------------------------------------------------
//    if((parentMomentumVector-Lab4Momentum1-Lab4Momentum2).Mag()>1.e-6)
//    {
//        std::cout << "Problem with two-body kinematics" << std::endl;
//        (parentMomentumVector-Lab4Momentum1-Lab4Momentum2).Print();
//    }
//    
//    //------------------------------------------------
//    result.push_back(parentMomentumVector);
//    result.push_back(Lab4Momentum1);
//    result.push_back(Lab4Momentum2);
//    
//    return result;
//    
//    //    //Function to calculate the recoil TLorentzVector in the lab frame
//    //    std::vector<TLorentzVector> result;
//    //
//    //    double s = std::pow(massEnergies[0],2.) + std::pow(massEnergies[1],2.) + 2 * massEnergies[1] * (TBeam + massEnergies[0]);
//    //
//    //    TLorentzVector f4MomentumCentreOfMass(0,0,0,std::sqrt(s));
//    //
//    //    double ECM0 = (s + std::pow(massEnergies[0],2.) - std::pow(massEnergies[1],2.))/(2*std::sqrt(s));
//    //    double ECM1 = (s + std::pow(massEnergies[1],2.) - std::pow(massEnergies[0],2.))/(2*std::sqrt(s));
//    //    double ECM2 = (s + std::pow(massEnergies[2],2.) - std::pow(massEnergies[3],2.))/(2*std::sqrt(s));
//    //    double ECM3 = (s + std::pow(massEnergies[3],2.) - std::pow(massEnergies[2],2.))/(2*std::sqrt(s));
//    //
//    //    //    std::cout << "ECM2: " << ECM2 << std::endl;
//    //    //    std::cout << "ECM3: " << ECM3 << std::endl;
//    //
//    //    double PCM0 = std::sqrt(std::pow(ECM0,2.) - std::pow(massEnergies[0],2.));
//    //    double PCM1 = std::sqrt(std::pow(ECM1,2.) - std::pow(massEnergies[1],2.));
//    //    double PCM2 = std::sqrt(std::pow(ECM2,2.) - std::pow(massEnergies[2],2.));
//    //    double PCM3 = std::sqrt(std::pow(ECM3,2.) - std::pow(massEnergies[3],2.));
//    //
//    //    TVector3 Lab3Momentum0(0,0,std::sqrt(std::pow(TBeam,2.) + 2 * TBeam * massEnergies[0]));
//    //    TVector3 Lab3Momentum1(0,0,0);
//    //
//    //    TLorentzVector Lab4Momentum0(Lab3Momentum0,massEnergies[0] + TBeam);
//    //    TLorentzVector Lab4Momentum1(Lab3Momentum1,massEnergies[1]);
//    //
//    //    double BetaCM = (Lab4Momentum0+Lab4Momentum1).Beta();
//    //    //std::cout << "BetaCM = " << BetaCM << std::endl;
//    //
//    //    TLorentzVector CoM4Momentum0 = Lab4Momentum0;
//    //    CoM4Momentum0.Boost(TVector3(0,0,-BetaCM));
//    //
//    //    TLorentzVector CoM4Momentum1 = Lab4Momentum1;
//    //    CoM4Momentum1.Boost(TVector3(0,0,-BetaCM));
//    //
//    //    TLorentzVector CoM4Momentum2 = TLorentzVector(PCM2 * std::sin(thetaCM * M_PI/180.) * std::cos(phiCM * M_PI/180.),
//    //                                                  PCM2 * std::sin(thetaCM * M_PI/180.) * std::sin(phiCM * M_PI/180.),
//    //                                                  PCM2 * std::cos(thetaCM * M_PI/180.),
//    //                                                  ECM2);
//    //    TLorentzVector CoM4Momentum3 = f4MomentumCentreOfMass - CoM4Momentum2;
//    //
//    //    TLorentzVector Lab4Momentum2 = CoM4Momentum2;
//    //    Lab4Momentum2.Boost(0,0,BetaCM);
//    //    TLorentzVector Lab4Momentum3 = CoM4Momentum3;
//    //    Lab4Momentum3.Boost(0,0,BetaCM);
//    //
//    //    //if(VerboseFlag)Lab4Momentum2.Print();
//    //    //if(VerboseFlag)Lab4Momentum3.Print();
//    //
//    //    if((Lab4Momentum0+Lab4Momentum1-Lab4Momentum2-Lab4Momentum3).Mag()>1.e-6)
//    //    {
//    //        std::cout << "Problem with two-body kinematics" << std::endl;
//    //        (Lab4Momentum0+Lab4Momentum1-Lab4Momentum2-Lab4Momentum3).Print();
//    //    }
//    //    
//    //    result.push_back(Lab4Momentum0);
//    //    result.push_back(Lab4Momentum1);
//    //    result.push_back(Lab4Momentum2);
//    //    result.push_back(Lab4Momentum3);
//    //    
//    //    return result;
//    
//}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

