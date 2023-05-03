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


int main(int argc, char ** argv) {
	
	#ifdef G4MULTITHREADED
		  G4MTRunManager* runManager = new G4MTRunManager;
		  runManager->SetNumberOfThreads(68);
	#else
		  G4RunManager* runManager = new G4RunManager;
	#endif

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

  return 0;
}
