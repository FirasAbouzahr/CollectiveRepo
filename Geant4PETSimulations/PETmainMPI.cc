#include "G4RunManager.hh"
#include "G4MTRunManager.hh"


#include "G4UImanager.hh"
#include "G4String.hh"

#include "PETDetectorConstruction.hh"
#include "PETPhysicsList.hh"
#include "PETActionInitialization.hh"


#include "QGSP_BIC_HP.hh"


#include "G4VisExecutive.hh"



#include "G4UIExecutive.hh"


// for parallelization for HPC
#include "G4MPImanager.hh"
#include "G4MPIsession.hh"


int main(int argc, char ** argv) {
    
  G4MPImanager* g4MPI= new G4MPImanager(argc,argv);
  G4MPIsession* session= g4MPI-> GetMPIsession();
  G4RunManager* runManager= new G4RunManager();
    
  //G4RunManager * runManager = new G4RunManager;
  
  PETDetectorConstruction* detector = new PETDetectorConstruction();


  runManager->SetUserInitialization(detector);

  runManager->SetUserInitialization(new PETPhysicsList());

  runManager->SetUserInitialization(new PETActionInitialization(detector));





    G4VisManager * visManager = new G4VisExecutive;
    visManager->Initialize();


    runManager->Initialize();

  // get the pointer to the UI manager and set verbosities
  G4UImanager * UImanager = G4UImanager::GetUIpointer();

  if (argc == 1) {
    UImanager->ApplyCommand("/control/execute geom.in");

    G4UIExecutive * ui = new G4UIExecutive(argc, argv);

    UImanager->ApplyCommand("/control/execute vis.mac");

    if (ui->IsGUI()) {
      UImanager->ApplyCommand("/control/execute gui.mac");
    }
    ui->SessionStart();
    delete ui;

  } else {
    G4String command = "/control/execute ";
    G4String filename = argv[1];
    UImanager->ApplyCommand(command + filename);
  }


  delete visManager;


  delete runManager;
  delete g4MPI;

  return 0;
}
