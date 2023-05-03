//Kyle Klein 
#include "PETStackingAction.hh"
#include "PETRunAction.hh"

#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Track.hh"
#include "G4RunManager.hh"
#include "PETUserEventInformation.hh"




PETStackingAction::PETStackingAction(PETRunAction *runAction)
  : fRunAction(runAction), gammaCounter(0), CerenCounter(0)
{;}

PETStackingAction::~PETStackingAction()
{;}

/*PETUserEventInformation* eventInformation=
  (PETUserEventInformation*)G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetUserInformation();*/


G4ClassificationOfNewTrack
PETStackingAction::ClassifyNewTrack(const G4Track * aTrack)
{


  if(aTrack->GetParentID()==0)return fUrgent;


  G4String name = aTrack->GetDefinition()->GetParticleName();
  G4int EventID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
  G4String CPName  = aTrack->GetCreatorProcess()->GetProcessName();
  G4int TrackID = aTrack->GetTrackID();

  if(aTrack->GetParentID()>0){
    if(name=="e-"){
		G4String process = aTrack->GetCreatorProcess()->GetProcessName();
		G4double Sec_energy = aTrack->GetKineticEnergy() / CLHEP::keV;
		G4double Sec_x  = aTrack->GetPosition().x();
		G4double Sec_y  = aTrack->GetPosition().y();
		G4double Sec_z  = aTrack->GetPosition().z();
		G4double Time = aTrack->GetGlobalTime()/CLHEP::ns;
		//	std::ofstream myfile("Sec_depth.txt", std::ios_base::app);
		//myfile << Sec_x <<std::endl;
		//myfile.close();
		// std::cout<<name<<"\t" <<  Sec_energy << std::endl;
		//fRunAction->FillSecondaryTree(Sec_x, Sec_y, Sec_z, Sec_energy);
	
	  std::ofstream myfile8("Edata15.txt", std::ios_base::app);
	  myfile8 << CPName << "," << EventID << "," << Sec_energy << "," << Time << "," << Sec_x << "," << Sec_y << "," << Sec_z << "," << TrackID << std::endl;
	  myfile8.close();

    }
  }
  //return fUrgent;

  //if(aTrack->GetDefinition()==G4OpticalPhoton::OpticalPhotonDefinition())
    //{ // particle is optical photon
      //G4double PhotTime1 = aTrack->GetGlobalTime()/CLHEP::ns;
      //if(aTrack->GetParentID()>0){


	 //// particle is secondary
	  //gammaCounter++;
    //if(aTrack->GetCreatorProcess()->GetProcessName()=="Scintillation"){
    ////eventInformation->IncPhotonCount_Scint();
    //}
    //else if(aTrack->GetCreatorProcess()->GetProcessName()=="Cerenkov"){
    ////eventInformation->IncPhotonCount_Ceren();
    ////CerenCounter++;
  ////  std::cout << Time1 << std::endl;
  //}
//}



  //}


    //std::ofstream myfile1("gammacount_15_40x3_pix_time_g5.txt", std::ios_base::app);
    //myfile1 << gammaCounter <<std::endl;
    //myfile1.close();
    ////fRunAction->FillGammaTree(gammaCounter);
  ////return fKill;
  return fUrgent;
}

void PETStackingAction::NewStage()
{

/*std::ofstream myfile2("CerenCount.txt", std::ios_base::app);
  myfile2 << CerenCounter <<std::endl;
  myfile2.close();*/


  //G4cout << "Number of gamma produces in this event : "
	// << gammaCounter << G4endl;

}

void PETStackingAction::PrepareNewEvent()
{ gammaCounter = 0;
  CerenCounter = 0;}
