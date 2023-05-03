//
// Dung Phan, Kyle Klein 
//

#include "G4Run.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ProcessManager.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4VProcess.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ThreeVector.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4LogicalVolume.hh"
#include "G4TrackingManager.hh"
#include "G4Electron.hh"
#include "G4EventManager.hh"


#include "PETRunAction.hh"
#include "PETSteppingAction.hh"
#include "PETDetectorConstruction.hh"
#include "PETTrackerSD.hh"
#include "G4EventManager.hh"

#include <sstream>
#include <iostream>
#include <fstream>
#include <random>

PETSteppingAction::PETSteppingAction(PETDetectorConstruction * detector) : fDetector(detector) {
  DetectedCerenCounter = 0, fOpProcess = NULL;
}

PETSteppingAction::~PETSteppingAction() {}

void PETSteppingAction::UserSteppingAction(const G4Step * theStep) {
  G4Track* theTrack = theStep->GetTrack();

  G4int EventID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

  G4String CPName;

  G4double Secondary_Energy;
  G4double TotalSecondary_Energy;
  G4int gammaCounter=0;
  G4int TrackID;
  G4int ParentID;
  G4int Step_Number = theTrack->GetCurrentStepNumber();
  G4double Energy_gamma = theTrack->GetTotalEnergy();
  G4double hPlanckTimesc = 1239.842; //Planck's constant times the speed of light, in eV*nm


  //  if(theTrack->GetCreatorProcess()!=0){
  //  if(theTrack->GetParentID()==0){
    //if(theTrack->GetCreatorProcess()->GetProcessName()=="compt"){
    //      if(Step_Number == 2){
    //CPName = theTrack->GetCreatorProcess()->GetProcessName();
    //G4double Energy_Loss = theStep->GetDeltaEnergy();
    //G4double Energy_Gamma = theTrack->GetTotalEnergy();
    //Secondary_Energy = theTrack->GetKineticEnergy()/keV;
    //ParentID = theTrack->GetParentID();
    //G4double x  = theTrack->GetPosition().x();//GetVertexPosition().x();
    //G4double y  = theTrack->GetPosition().y();//GetVertexPosition().y();
    //G4double z  = theTrack->GetPosition().z();//GetVertexPosition().z();
    //if(Energy_Loss!=0){
    //std::cout<<Step_Number<< "\t" <<Energy_Loss<< "\t" << Energy_Gamma <<std::endl;
    //std::ofstream myfile("gamma_energy.txt", std::ios_base::app);
    //myfile1 << Step_Number  << "," << Energy_Loss << "," << Energy_Gamma << "," << x << "," << y << "," << z <<std::endl;
    //myfile1.close();
    //}
    //}
  //}
  //    }


  G4StepPoint* thePrePoint  = theStep->GetPreStepPoint();
  G4StepPoint* thePostPoint = theStep->GetPostStepPoint();
  G4StepStatus Step_Status = thePostPoint->GetStepStatus();
  G4String procName = thePostPoint->GetProcessDefinedStep()->GetProcessName();

  G4VPhysicalVolume* thePrePV  = thePrePoint->GetPhysicalVolume();
  G4VPhysicalVolume* thePostPV = thePostPoint->GetPhysicalVolume();



  G4String thePrePVname  = " ";
  G4String thePostPVname = " ";

  if (thePostPV) {
    thePrePVname  = thePrePV->GetName();
    thePostPVname = thePostPV->GetName();
  }

  G4TouchableHistory* theTouchable = (G4TouchableHistory*)(thePostPoint->GetTouchable());
  G4int copyNumber = theTouchable->GetCopyNumber();



  // if (theTrack->GetParentID() > 0){
  //   std::cout<<thePostPVname<<std::endl;
  // }

  // Retrieve the status of the photon
  G4OpBoundaryProcessStatus theStatus = Undefined;
  G4ProcessManager* OpManager = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();

  if (OpManager) {
    G4int MAXofPostStepLoops = OpManager->GetPostStepProcessVector()->entries();
    G4ProcessVector* fPostStepDoItVector = OpManager->GetPostStepProcessVector(typeDoIt);

    for ( G4int i=0; i<MAXofPostStepLoops; i++) {
      G4VProcess* fCurrentProcess = (*fPostStepDoItVector)[i];
      fOpProcess = dynamic_cast<G4OpBoundaryProcess*>(fCurrentProcess);
      if (fOpProcess) { theStatus = fOpProcess->GetStatus(); break;}
    }
  }


  //  {
    //      switch (theStatus)
  //	{
	  // Detected by a detector
	  //case Detection:
	  //	case 12:
	  // Check if the photon hits the detector and process the hit if it does
    //if(theTrack->GetDefinition()==G4Electron::ElectronDefinition())
      //{
        //if (thePostPVname == "crystal") {
        //// theTrack = G4Step::GetPostStepPoint();
        //const G4ThreeVector EPos = theTrack->GetPosition();
        //G4double ETime = theTrack->GetGlobalTime()/CLHEP::ns;
          //std::ofstream myfile2("Electron_Time_Tolerance.txt", std::ios_base::app);
          //myfile2 << EventID << "," << ETime << "," << EPos << std::endl;
          //myfile2.close();
          //}

      //}

 // Important according to Kyle
  if(theTrack->GetDefinition()==G4Gamma::GammaDefinition())
    {
      if (thePostPVname == "crystal") {
//       theTrack = G4Step::GetPostStepPoint();
      const G4ThreeVector GammaPos = theTrack->GetPosition();
//      G4double GammaX  = theTrack->GetPosition().x();
//      G4double GammaY = theTrack->GetPosition().y();
//      G4double GammaZ = theTrack->GetPosition().z();
      G4double GammaTime = theTrack->GetGlobalTime()/CLHEP::ns;
      G4double GammaEnergy = theTrack->GetTotalEnergy() / CLHEP::keV;
        std::ofstream myfile3("GammaTracker.txt", std::ios_base::app);
        myfile3 << EventID << "," << copyNumber << "," << GammaPos << "," << thePostPVname << std::endl;
        myfile3.close();
        }

      }

  if(theTrack->GetDefinition()==G4OpticalPhoton::OpticalPhotonDefinition()) //specify optical photons
    {
		//G4int TrackID = theTrack->GetTrackID();
		//G4double PhotTime = theTrack->GetGlobalTime()/CLHEP::ns;
        //G4double x  = theTrack->GetPosition().x();
        //G4double y  = theTrack->GetPosition().y();
        //G4double z  = theTrack->GetPosition().z();
        //G4double Energy_photon = theTrack->GetTotalEnergy() / eV;
        //G4double Wavelength_photon = hPlanckTimesc / Energy_photon;
        //Step_Number = theTrack->GetCurrentStepNumber();
        //std::ofstream myfile10("DETAILEDTRACKING3x315.txt", std::ios_base::app);
        //myfile10 << EventID << "," << TrackID << "," << thePostPVname << "," << x << "," << y << "," << z << "," << PhotTime << "," << Wavelength_photon  << "," << Step_Number << std::endl;
        //myfile10.close();
        
		
		
	  //if((theTrack->GetTrackStatus() == fStopAndKill) and (thePostPVname != "detector")){
	        //G4double PhotTime = theTrack->GetGlobalTime()/CLHEP::ns;
	        //G4double x  = theTrack->GetPosition().x();
	        //G4double y  = theTrack->GetPosition().y();
	        //G4double z  = theTrack->GetPosition().z();
	        //TrackID = theTrack->GetTrackID();
	        //G4double Energy_photon = theTrack->GetTotalEnergy() / eV;
	        //G4double Wavelength_photon = hPlanckTimesc / Energy_photon;
	        //Step_Number = theTrack->GetCurrentStepNumber();
	        //std::ofstream myfile4("KILLED3x315.txt", std::ios_base::app);
	        //myfile4 << EventID << "," << x << "," << y << "," << z << "," << PhotTime << "," << Wavelength_photon  << "," << Step_Number << "," << procName << std::endl;
	        //myfile4.close();
      //}
      //if((theTrack->GetTrackStatus() == fStopAndKill) and (thePostPVname == "detector")){
	        //G4double PhotTime = theTrack->GetGlobalTime()/CLHEP::ns;
	        //G4double x  = theTrack->GetPosition().x();
	        //G4double y  = theTrack->GetPosition().y();
	        //G4double z  = theTrack->GetPosition().z();
	        //TrackID = theTrack->GetTrackID();
	        //G4double Energy_photon = theTrack->GetTotalEnergy() / eV;
	        //G4double Wavelength_photon = hPlanckTimesc / Energy_photon;
	        //Step_Number = theTrack->GetCurrentStepNumber();
	        //std::ofstream myfile5("DETECTORLOSS3x315.txt", std::ios_base::app);
	        //myfile5 << EventID << "," << x << "," << y << "," << z << "," << PhotTime << "," << Wavelength_photon  << "," << Step_Number << "," << procName << std::endl;
	        //myfile5.close();
      ////}
	  if(theTrack->GetCurrentStepNumber() == 1){
	        G4double PhotTime = theTrack->GetGlobalTime()/CLHEP::ns;
	        G4double x  = thePrePoint->GetPosition().x();
	        G4double y  = thePrePoint->GetPosition().y();
          G4String Masterprocesstype = theTrack->GetCreatorProcess()->GetProcessName();
	        G4double z  = thePrePoint->GetPosition().z();
	        TrackID = theTrack->GetTrackID();
	        G4String process = theTrack->GetCreatorProcess()->GetProcessName();
	        G4double Energy_photon = theTrack->GetTotalEnergy() / eV;
	        G4double Wavelength_photon = hPlanckTimesc / Energy_photon;
          std::ofstream myfile6("MasterPhotonFile.txt", std::ios_base::app);
          myfile6 << EventID << "," <<  Energy_photon << "," << Masterprocesstype << std::endl;
          myfile6.close();
      }
     		
      if ((thePostPVname == "detector") & (theTrack->GetTrackStatus() != fStopAndKill)) {
        double sipmResolution = 100; //Resolution of the SiPMs, in ps FWHM
		// random device class instance, source of 'true' randomness for initializing random seed
		std::random_device rd;
		// Mersenne twister PRNG, initialized with seed from previous random device instance
		std::mt19937 gen(rd());
		// instance of class std::normal_distribution with specific mean and stddev
		std::normal_distribution<float> d(sipmResolution/1000, sipmResolution/1000/2.355);
        
        G4double PhotTime = theTrack->GetGlobalTime()/CLHEP::ns;
        G4String processtype = theTrack->GetCreatorProcess()->GetProcessName();
        G4double photX  = theTrack->GetPosition().x();
        G4double photY = theTrack->GetPosition().y();
        G4double photZ = theTrack->GetPosition().z();
        G4double Energy_photon = theTrack->GetTotalEnergy() / eV;
        G4double Wavelength_photon = hPlanckTimesc / Energy_photon;
        TrackID = theTrack->GetTrackID();
        

       //if(theTrack->GetCreatorProcess()->GetProcessName()=="Cerenkov"){
        //DetectedCerenCounter++;
        //G4double Time1 = theTrack->GetGlobalTime()/CLHEP::ns;
        ////std::ofstream myfile7("CerenData.txt", std::ios_base::app);
        ////myfile7 << EventID <<  "," << DetectedCerenCounter << "," << Time1 << "," << copyNumber << std::endl;
      //}
    //  if( > 1){
      //  std::ofstream myfile8("COMPT_PHOT.txt", std::ios_base::app);
   	  // myfile8 << EventID << "," << copyNumber << "," << x << "," << PhotTime << std::endl;
   	  // myfile8.close();
    //  }
	   //std::cout<<copyNumber<<std::endl;\


          std::ofstream myfile9("photo.txt", std::ios_base::app);
          myfile9 << EventID << "," << copyNumber << "," << PhotTime + d(gen) << "," << processtype << "," << Wavelength_photon << std::endl;
          myfile9.close();
     G4SDManager* SDman = G4SDManager::GetSDMpointer();
     G4String SDname="TrackerChamberSD";
     PETTrackerSD* mppcSD = (PETTrackerSD*)SDman->FindSensitiveDetector(SDname);
     if (mppcSD) mppcSD->ProcessHits_constStep(theStep,NULL);
	// Stop Tracking when it hits the detector's surface
	   theTrack->SetTrackStatus(fStopAndKill);
	   return;
      }
      //	      break;
      //	      default: break;
      //	    }
    }


  }
  void PETSteppingAction::PrepareNewEvent()
  {DetectedCerenCounter = 0;}
