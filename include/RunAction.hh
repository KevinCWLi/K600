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

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include <vector>

class G4Run;

/// Run action class
///
/// It accumulates statistic and computes dispersion of the energy deposit
/// and track lengths of charged particles with use of analysis tools:
/// H1D histograms are created in BeginOfRunAction() for the following
/// physics quantities:
/// - Edep in absorber
/// - Edep in gap
/// - Track length in absorber
/// - Track length in gap
/// The same values are also saved in the ntuple.
/// The histograms and ntuple are saved in the output file in a format
/// accoring to a selected technology in B4Analysis.hh.
///
/// In EndOfRunAction(), the accumulated statistic and computed
/// dispersion is printed.
///

class RunAction : public G4UserRunAction
{
public:
    RunAction();
    virtual ~RunAction();
    
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);
    
    //--------------------------------------------------------------------------------
    //      Initial Simulated Particles
    std::vector<double>   iSimulatedParticle_T; // MeV
    
    void SetInitialParticleKineticEnergy(std::vector<double> vec) {iSimulatedParticle_T = vec;};
    
    //--------------------------------------------------------------------------------
    //      Detectable energy

    std::vector<double>   detectableEnergy; // MeV
    
    void SetDetectableEnergy(std::vector<double> vec) {detectableEnergy = vec;};

    //--------------------------------------------------------------------------------
    //      CAKE
    std::vector<int>    CAKE_detNr;
    std::vector<int>    CAKE_ringNr;
    std::vector<int>    CAKE_sectorNr;
    std::vector<double> CAKE_hitRadius;
    
    void SetCAKE_detNr(std::vector<int> vec)     {CAKE_detNr = vec;};
    void SetCAKE_ringNr(std::vector<int> vec)    {CAKE_ringNr = vec;};
    void SetCAKE_sectorNr(std::vector<int> vec)  {CAKE_sectorNr = vec;};
    void SetCAKE_hitRadius(std::vector<double> vec)  {CAKE_hitRadius = vec;};

    //--------------------------------------------------------------------------------
    //      Reaction variables
    std::vector<int>    nReactionProducts;
    std::vector<int>    reaction_Z;
    std::vector<double> reaction_A;
    std::vector<double> reaction_P;
    std::vector<double> reaction_T;
    std::vector<double> reaction_Ex;
    
    void SetNReactionProducts(int n)
    {
        nReactionProducts.clear();
        nReactionProducts.push_back(n);
    };

    void SetReaction_Z(std::vector<int> vec)    {reaction_Z = vec;};
    void SetReaction_A(std::vector<double> vec) {reaction_A = vec;};
    void SetReaction_P(std::vector<double> vec) {reaction_P = vec;};
    void SetReaction_T(std::vector<double> vec) {reaction_T = vec;};
    void SetReaction_Ex(std::vector<double> vec) {reaction_Ex = vec;};

    std::vector<double> reaction_theta_LAB;
    std::vector<double> reaction_phi_LAB;
    std::vector<double> reaction_theta_reactionCOM;
    std::vector<double> reaction_phi_reactionCOM;
    
    void SetReaction_theta_LAB(std::vector<double> vec)         {reaction_theta_LAB = vec;};
    void SetReaction_phi_LAB(std::vector<double> vec)           {reaction_phi_LAB = vec;};
    void SetReaction_theta_reactionCOM(std::vector<double> vec) {reaction_theta_reactionCOM = vec;};
    void SetReaction_phi_reactionCOM(std::vector<double> vec)   {reaction_phi_reactionCOM = vec;};

    //--------------------------------------------------------------------------------
    //      Decay variables
    std::vector<int>    nDecayProducts;
    std::vector<int>    decay_Z;
    std::vector<double> decay_A;
    std::vector<double> decay_P;
    std::vector<double> decay_T;

    void SetNDecayProducts(int n)
    {
        nDecayProducts.clear();
        nDecayProducts.push_back(n);
    };

    void SetDecay_Z(std::vector<int> vec)       {decay_Z = vec;};
    void SetDecay_A(std::vector<double> vec)    {decay_A = vec;};
    void SetDecay_P(std::vector<double> vec)    {decay_P = vec;};
    void SetDecay_T(std::vector<double> vec)    {decay_T = vec;};

    std::vector<double> decay_theta_LAB;
    std::vector<double> decay_phi_LAB;
    std::vector<double> decay_theta_reactionCOM;
    std::vector<double> decay_phi_reactionCOM;
    std::vector<double> decay_theta_recoilCOM;
    std::vector<double> decay_phi_recoilCOM;

    void SetDecay_theta_LAB(std::vector<double> vec)            {decay_theta_LAB = vec;};
    void SetDecay_phi_LAB(std::vector<double> vec)              {decay_phi_LAB = vec;};
    void SetDecay_theta_recoilCOM(std::vector<double> vec)      {decay_theta_recoilCOM = vec;};
    void SetDecay_phi_recoilCOM(std::vector<double> vec)        {decay_phi_recoilCOM = vec;};
    void SetDecay_theta_reactionCOM(std::vector<double> vec)    {decay_theta_reactionCOM = vec;};
    void SetDecay_phi_reactionCOM(std::vector<double> vec)      {decay_phi_reactionCOM = vec;};
    
    
    //--------------------------------------------------------------------------------
    //      CLOVER detectors
    std::vector<int>    CLOVER_eventFold;
    std::vector<int>    CLOVER_iD;
    std::vector<int>    CLOVER_nCrystalsTriggered;
    std::vector<double> CLOVER_energyPerCrystal;
    std::vector<double> CLOVER_energy;
    std::vector<double> CLOVER_initialEnergy;
    std::vector<double> CLOVER_initialEnergy_COM;
    std::vector<double> CLOVER_detectorTheta;
    std::vector<double> CLOVER_detectorPhi;
    std::vector<int>    CLOVER_CrystalReflectionIndex;
    std::vector<double> CLOVER_initialInteractionTheta;
    std::vector<double> CLOVER_initialInteractionPhi;
    std::vector<double> CLOVER_initialParticleTheta;
    std::vector<double> CLOVER_initialParticlePhi;

    std::vector<int>    CLOVER_BGOCrystalsTriggered;

    void SetCLOVER_EventFold(int n)
    {
        CLOVER_eventFold.clear();
        CLOVER_eventFold.push_back(n);
    };
    
    void SetCLOVER_IDs(std::vector<int> vec) {CLOVER_iD = vec;};
    void SetCLOVER_NCrystalsTriggered(std::vector<int> vec) {CLOVER_nCrystalsTriggered = vec;};
    void SetCLOVER_EnergiesPerCrystal(std::vector<double> vec) {CLOVER_energyPerCrystal = vec;};
    void SetCLOVER_Energies(std::vector<double> vec) {CLOVER_energy = vec;};
    void SetCLOVER_InitialEnergies(std::vector<double> vec) {CLOVER_initialEnergy = vec;};
    void SetCLOVER_InitialEnergiesCOM(std::vector<double> vec) {CLOVER_initialEnergy_COM = vec;};
    void SetCLOVER_DetectorThetas(std::vector<double> vec) {CLOVER_detectorTheta = vec;};
    void SetCLOVER_DetectorPhis(std::vector<double> vec) {CLOVER_detectorPhi = vec;};
    void SetCLOVER_CrystalReflectionIndices(std::vector<int> vec) {CLOVER_CrystalReflectionIndex = vec;};
    void SetCLOVER_InitialInteractionThetas(std::vector<double> vec) {CLOVER_initialInteractionTheta = vec;};
    void SetCLOVER_InitialInteractionPhis(std::vector<double> vec) {CLOVER_initialInteractionPhi = vec;};
    void SetCLOVER_InitialParticleThetas(std::vector<double> vec) {CLOVER_initialParticleTheta = vec;};
    void SetCLOVER_InitialParticlePhis(std::vector<double> vec) {CLOVER_initialParticlePhi = vec;};

    void SetCLOVER_BGOCrystalsTriggered(std::vector<int> vec) {CLOVER_BGOCrystalsTriggered = vec;};

    //--------------------------------------------------------------------------------
    //      LaBr3Ce detectors
    std::vector<int>    laBr3Ce_eventFold;
    std::vector<int>    laBr3Ce_iD;
    std::vector<double> laBr3Ce_energy;
    std::vector<double> laBr3Ce_detectorTheta;
    std::vector<double> laBr3Ce_detectorPhi;
    std::vector<double> laBr3Ce_theta;
    std::vector<double> laBr3Ce_phi;
    std::vector<double> laBr3Ce_xPos;
    std::vector<double> laBr3Ce_yPos;
    std::vector<double> laBr3Ce_zPos;
    
    void SetLaBr3Ce_EventFold(int n)
    {
        laBr3Ce_eventFold.clear();
        laBr3Ce_eventFold.push_back(n);
    };

    void SetLaBr3Ce_IDs(std::vector<int> vec) {laBr3Ce_iD = vec;};
    void SetLaBr3Ce_Energies(std::vector<double> vec) {laBr3Ce_energy = vec;};
    void SetLaBr3Ce_Thetas(std::vector<double> vec) {laBr3Ce_theta = vec;};
    void SetLaBr3Ce_Phis(std::vector<double> vec) {laBr3Ce_phi = vec;};
    void SetLaBr3Ce_DetectorThetas(std::vector<double> vec) {laBr3Ce_detectorTheta = vec;};
    void SetLaBr3Ce_DetectorPhis(std::vector<double> vec) {laBr3Ce_detectorPhi = vec;};
    void SetLaBr3Ce_xPos(std::vector<double> vec) {laBr3Ce_xPos = vec;};
    void SetLaBr3Ce_yPos(std::vector<double> vec) {laBr3Ce_yPos = vec;};
    void SetLaBr3Ce_zPos(std::vector<double> vec) {laBr3Ce_zPos = vec;};

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

