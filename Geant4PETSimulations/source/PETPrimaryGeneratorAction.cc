//Kyle Klein
#include "PETPrimaryGeneratorAction.hh"

PETPrimaryGeneratorAction::PETPrimaryGeneratorAction() : G4VUserPrimaryGeneratorAction(), fParticleGun(0) {
  G4int n_particle = 1;

  fParticleGun = new G4ParticleGun(n_particle);
  G4ParticleTable * particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition *particle = particleTable->FindParticle(particleName = "gamma");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleEnergy(0.511 * MeV);

  tParticleGun = new G4ParticleGun(n_particle);
  tParticleGun->SetParticleDefinition(particle);
  tParticleGun->SetParticleMomentumDirection(G4ThreeVector(-1, 0, 0.0));
  tParticleGun->SetParticleEnergy(0.511 * MeV);

}

PETPrimaryGeneratorAction::~PETPrimaryGeneratorAction() {
  delete fParticleGun;
}

void PETPrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent) {

  //Shoot straight into the middle of one crystal
  G4double mid_z = 324.317 * mm;
  G4double x0 = 0 * mm;
  G4double y0 = 0 * mm;
    
 // 20 cm long line source centered in middle of the cylinder
  G4double z0 = mid_z + G4UniformRand()*2*100 - 100;
//  G4double xMomentumDirection = 1;
//  G4double yMomentumDirection = 0;
//  G4double zMomentumDirection = 0;
    
  G4double theta = CLHEP::twopi * G4UniformRand();
  G4double phi = CLHEP::pi *acos(1-2* G4UniformRand());
  G4double xMomentumDirection = sin(phi)*cos(theta);
  G4double yMomentumDirection = sin(phi)*sin(theta);
  G4double zMomentumDirection = cos(phi);
    
  fParticleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(xMomentumDirection, yMomentumDirection, zMomentumDirection)); 
  fParticleGun->GeneratePrimaryVertex(anEvent);
  tParticleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));
  fParticleGun->SetParticleTime(1000*ns);
  tParticleGun->SetParticleMomentumDirection(G4ThreeVector(-xMomentumDirection, -yMomentumDirection, -zMomentumDirection)); 
  tParticleGun->GeneratePrimaryVertex(anEvent);

  ////Shoot straight from the origin
  //G4double x0 = 0 * mm;
  //G4double y0 = 0 * mm;
  //G4double z0 = (0) * mm;
  //G4double xMomentumDirection = 1;
  //G4double yMomentumDirection = 0;
  //G4double zMomentumDirection = 0;
  //fParticleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));
  //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(xMomentumDirection, yMomentumDirection, zMomentumDirection)); 
  //fParticleGun->GeneratePrimaryVertex(anEvent);

  ////Randomize over all solid angle, from a line source, and shoot back-to-back gammas, one much later than the other for easy differentiation:
  //G4double x0 = 100 * mm;
  ////G4double y0 = 0 * mm;
  //G4double y0 = 0 * mm;
  ////G4double z0 = (300*G4UniformRand()-150-164.4/2) * mm;
  //G4double z0 = (0) * mm;
  ////G4double z0 = (-210.6/2 + 25.795/2 + 1.6) * mm;
  //G4double theta = CLHEP::twopi * G4UniformRand();
  //G4double phi = CLHEP::pi *acos(1-2* G4UniformRand());
  //G4double xMomentumDirection = sin(phi)*cos(theta);
  //G4double yMomentumDirection = sin(phi)*sin(theta);
  //G4double zMomentumDirection = cos(phi);
  ////G4double zMomentumDirection = 0;
  //fParticleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));
  //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(xMomentumDirection, yMomentumDirection, zMomentumDirection)); 
  //fParticleGun->SetParticleTime(0*ns);
  //fParticleGun->GeneratePrimaryVertex(anEvent);
  //fParticleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));
  //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(-xMomentumDirection, -yMomentumDirection, -zMomentumDirection));
  //fParticleGun->SetParticleTime(1000*ns);
  //fParticleGun->GeneratePrimaryVertex(anEvent);
  
  ////Randomize over one 8x8 array, including its outer gaps:
  //G4double x0 = 0 * mm;
  //G4double y0 = 0.0 * mm;
  //G4double z0 = 0.0 * mm;
  //G4double theta = 0.0714*G4UniformRand();
  //G4double phi = 0.0714 * G4UniformRand();
  //G4double xMomentumDirection = 1;
  //G4double yMomentumDirection = tan(phi);
  //G4double zMomentumDirection = tan(theta);
  //fParticleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));
  //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(xMomentumDirection, yMomentumDirection, zMomentumDirection)); 
  //fParticleGun->GeneratePrimaryVertex(anEvent);
  
  ////Randomize over central 4x4 crystals, not including their outer gaps:
  //G4double xMomentumDirection = 1;
  //G4double yMomentumDirection = ((G4UniformRand() *  (2 * 6.3 / 163.91)) - (6.3 / 163.91));
  //G4double zMomentumDirection = ((G4UniformRand() *  (2 * 6.3 / 163.91)) - (6.3 / 163.91));
  
  //fParticleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));
  //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(xMomentumDirection, yMomentumDirection, zMomentumDirection)); 
  //fParticleGun->GeneratePrimaryVertex(anEvent);
  ////Include if you want back-to-back gammas
  //tParticleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));
  //tParticleGun->SetParticleMomentumDirection(G4ThreeVector(-xMomentumDirection, -yMomentumDirection, -zMomentumDirection));
  //tParticleGun->GeneratePrimaryVertex(anEvent);
  
  ////Shoot back-to-back gammas radially from center with random 2*pi angle
  //G4double x0 = 0 * mm;
  //G4double y0 = 0 * mm;
  //G4double z0 = 1.6 * mm;
  //G4double theta = CLHEP::twopi*G4UniformRand();
  //G4double xMomentumDirection = cos(theta);
  //G4double yMomentumDirection = sin(theta);
  //G4double zMomentumDirection = 0;
  //fParticleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));
  //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(xMomentumDirection, yMomentumDirection, zMomentumDirection)); 
  //fParticleGun->GeneratePrimaryVertex(anEvent);
  //tParticleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));
  //tParticleGun->SetParticleMomentumDirection(G4ThreeVector(-xMomentumDirection, -yMomentumDirection, -zMomentumDirection));
  //tParticleGun->SetParticleTime(1000*ns);
  //tParticleGun->GeneratePrimaryVertex(anEvent);

}
