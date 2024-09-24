#include "G4HadReentrentException.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTable.hh"

#include "SimG4Core/CustomPhysics/interface/FullModelHadronicProcess.h"
#include "SimG4Core/CustomPhysics/interface/G4ProcessHelper.h"
#include "SimG4Core/CustomPhysics/interface/Decay3Body.h"
#include "SimG4Core/CustomPhysics/interface/CustomPDGParser.h"
#include "SimG4Core/CustomPhysics/interface/CustomParticle.h"

using namespace CLHEP;

FullModelHadronicProcess::FullModelHadronicProcess(G4ProcessHelper* aHelper, const G4String& processName)
    : G4VDiscreteProcess(processName), theHelper(aHelper) {}

FullModelHadronicProcess::~FullModelHadronicProcess() {}

G4bool FullModelHadronicProcess::IsApplicable(const G4ParticleDefinition& aP) {
  return theHelper->ApplicabilityTester(aP);
}

G4double FullModelHadronicProcess::GetMicroscopicCrossSection(const G4DynamicParticle* aParticle,
                                                              const G4Element* anElement,
                                                              G4double aTemp) {
  //Get the cross section for this particle/element combination from the ProcessHelper
  G4double InclXsec = theHelper->GetInclusiveCrossSection(aParticle, anElement);
  //  G4cout<<"Returned cross section from helper was: "<<InclXsec/millibarn<<" millibarn"<<G4endl;
  return InclXsec;
}

G4double FullModelHadronicProcess::GetMeanFreePath(const G4Track& aTrack, G4double, G4ForceCondition*) {
  G4Material* aMaterial = aTrack.GetMaterial();
  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
  G4double sigma = 0.0;

  G4int nElements = aMaterial->GetNumberOfElements();

  const G4double* theAtomicNumDensityVector = aMaterial->GetAtomicNumDensityVector();
  G4double aTemp = aMaterial->GetTemperature();

  for (G4int i = 0; i < nElements; ++i) {
    G4double xSection = GetMicroscopicCrossSection(aParticle, (*aMaterial->GetElementVector())[i], aTemp);
    sigma += theAtomicNumDensityVector[i] * xSection;
  }
  G4double res = DBL_MAX;
  if (sigma > 0.0) {
    res = 1. / sigma;
  }
  return res;
}

G4VParticleChange* FullModelHadronicProcess::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep) {
  const G4TouchableHandle& thisTouchable(aTrack.GetTouchableHandle());

  // Define the R-hadron as a CustomParticle named CustomIncident. Declare other variables of use later in the script.
  aParticleChange.Initialize(aTrack);
  const G4DynamicParticle* IncidentRhadron = aTrack.GetDynamicParticle();
  CustomParticle* CustomIncident = static_cast<CustomParticle*>(IncidentRhadron->GetDefinition());
  const G4ThreeVector& aPosition = aTrack.GetPosition();
  const G4int theIncidentPDG = IncidentRhadron->GetDefinition()->GetPDGEncoding();
  G4ParticleTable* theParticleTable = G4ParticleTable::GetParticleTable();
  std::vector<G4ParticleDefinition*> theParticleDefinitions;
  G4bool IncidentSurvives = false;
  G4bool TargetSurvives = false;
  G4Nucleus targetNucleus(aTrack.GetMaterial());
  G4ParticleDefinition* outgoingRhadron = nullptr;
  G4ParticleDefinition* outgoingCloud = nullptr;
  G4ParticleDefinition* outgoingTarget = nullptr;
  G4ThreeVector p_0 = IncidentRhadron->GetMomentum();
  G4double e_kin_0 = IncidentRhadron->GetKineticEnergy();

  // Declare the quark cloud as a G4DynamicParticle
  G4DynamicParticle* cloudParticle = new G4DynamicParticle();
  cloudParticle->SetDefinition(CustomIncident->GetCloud());
  if (cloudParticle->GetDefinition() == nullptr) {
    G4cout << "FullModelHadronicProcess::PostStepDoIt  Definition of particle cloud not available!!" << G4endl;
  }

  // Define the gluino and quark cloud G4LorentzVector (momentum, total energy) based on the momentum of the R-hadron and the ratio of the masses
  double scale = cloudParticle->GetDefinition()->GetPDGMass() / IncidentRhadron->GetDefinition()->GetPDGMass();
  G4LorentzVector cloudMomentum(IncidentRhadron->GetMomentum() * scale, cloudParticle->GetTotalEnergy());
  G4LorentzVector gluinoMomentum(IncidentRhadron->GetMomentum() * (1. - scale), CustomIncident->GetTotalEnergy() - cloudParticle->GetTotalEnergy());

  //These two for getting CMS transforms later (histogramming purposes...) NEED TO FIGURE OUT WHAT THIS DOES
  G4LorentzVector FullRhadron4Momentum = IncidentRhadron->Get4Momentum();
  const G4LorentzVector& Cloud4Momentum = cloudMomentum;
  //////////////////////////////////////////////////////////////////////////////////////////////////////////

  //Set the momentum of the quark cloud
  cloudParticle->Set4Momentum(cloudMomentum);

  // Declare the kinetic energies of the R-hadron and the quark cloud. Calculate the kinetic energy of the target nucleus.
  double incidentRHadronKineticEnergy = IncidentRhadron->GetKineticEnergy();
  G4double cloudKineticEnergy = cloudParticle->GetKineticEnergy();
  G4double targetNucleusKineticEnergy = targetNucleus.Cinema(cloudKineticEnergy);
  cloudKineticEnergy += targetNucleusKineticEnergy;
  targetNucleusKineticEnergy = targetNucleus.EvaporationEffects(cloudKineticEnergy); // calculate black track energies
  cloudKineticEnergy -= targetNucleusKineticEnergy;

  // If the R-hadron kinetic energy is less than 0.1 MeV, or the cloud kinetic energy is less than or equal to 0, stop the track but keep it alive. This should be very rare.
  if (cloudKineticEnergy + gluinoMomentum.e() - gluinoMomentum.m() <= 0.1 * MeV || cloudKineticEnergy <= 0.) {
    G4cout << "Kinetic energy is sick" << G4endl;
    G4cout << "Full R-hadron: " << (cloudKineticEnergy + gluinoMomentum.e() - gluinoMomentum.m()) / MeV << " MeV" << G4endl;
    G4cout << "Quark system: " << cloudKineticEnergy / MeV << " MeV" << G4endl;
    aParticleChange.ProposeTrackStatus(fStopButAlive);  // AR_NEWCODE_IMPORT
    return &aParticleChange;
  }
  cloudParticle->SetKineticEnergy(cloudKineticEnergy);

  //Get the final state particles
  G4ParticleDefinition* aTarget;
  ReactionProduct reactionProduct = theHelper->GetFinalState(aTrack, aTarget);
  G4int NumberOfSecondaries = reactionProduct.size();

  //Getting CMS transforms. Boosting is done at histogram filling
  G4LorentzVector Target4Momentum(0., 0., 0., aTarget->GetTotalEnergy());
  G4LorentzVector psum_full = FullRhadron4Momentum + Target4Momentum;
  G4LorentzVector psum_cloud = Cloud4Momentum + Target4Momentum;
  G4ThreeVector trafo_full_cms = (-1) * psum_full.boostVector();
  G4ThreeVector trafo_cloud_cms = (-1) * psum_cloud.boostVector();

  //Process outgoing particles from reactions
  for (ReactionProduct::iterator it = reactionProduct.begin(); it != reactionProduct.end(); ++it) {
    G4ParticleDefinition* tempDef = theParticleTable->FindParticle(*it);
    CustomParticle* tempCust = dynamic_cast<CustomParticle*>(tempDef);
    if (tempDef == aTarget)
      TargetSurvives = true;

    if (tempDef->GetParticleType()=="rhadron") {
      outgoingRhadron = tempDef;
      outgoingCloud = tempCust->GetCloud();
      if (outgoingCloud == nullptr) {
        G4cout << "FullModelHadronicProcess::PostStepDoIt  Definition of outgoing particle cloud not available!" << G4endl;
      }
    }

    if (tempDef == G4Proton::Proton() || tempDef == G4Neutron::Neutron()) outgoingTarget = tempDef;
    if (tempCust == nullptr && reactionProduct.size() == 2) outgoingTarget = tempDef;
    if (tempDef->GetPDGEncoding() == theIncidentPDG) {
      IncidentSurvives = true;
    } else {
      theParticleDefinitions.push_back(tempDef);
    }
  }

  //If no reaction occured, set the outgoingTarget to the incoming particle
  if (outgoingTarget == nullptr)
    outgoingTarget = theParticleTable->FindParticle(reactionProduct[1]);

  //If the incident particle survives, decrement the number of secondaries
  if (IncidentSurvives)
    NumberOfSecondaries--;
  aParticleChange.SetNumberOfSecondaries(NumberOfSecondaries);

  //Calculate the Lorentz rotation of the cloud particle to the lab frame
  G4HadProjectile* originalIncident = new G4HadProjectile(*cloudParticle);
  G4LorentzRotation cloudParticleToLabFrameRotation = originalIncident->GetTrafoToLab();

  //Create the current and target particles with proper momenta and kinetic energy
  G4DynamicParticle* originalTarget = new G4DynamicParticle;
  originalTarget->SetDefinition(aTarget);
  G4ReactionProduct targetParticle(aTarget);
  G4ReactionProduct currentParticle(const_cast<G4ParticleDefinition*>(originalIncident->GetDefinition()));
  currentParticle.SetMomentum(originalIncident->Get4Momentum().vect());
  currentParticle.SetKineticEnergy(originalIncident->GetKineticEnergy());
  G4ReactionProduct modifiedOriginal = currentParticle; // modifiedOriginal will have Fermi motion and evaporative effects included

  //Set the hemisphere of the current and target particles. Initialize an empty vector for the secondary particles
  currentParticle.SetSide(1);  // incident always goes in forward hemisphere
  targetParticle.SetSide(-1);  // target always goes in backward hemisphere
  G4bool quasiElastic = false;
  if (reactionProduct.size() == 2)
    quasiElastic = true;
  G4FastVector<G4ReactionProduct, MYGHADLISTSIZE> vec;  // vec will contain the secondary particles
  G4int vecLen = 0;
  vec.Initialize(0);

  //Fill the vector with the secondary particles. Here secondary particle is defined as any particle that is not the incident or target particle.
  for (G4int i = 0; i != NumberOfSecondaries; i++) {
    if (theParticleDefinitions[i] != aTarget && theParticleDefinitions[i] != originalIncident->GetDefinition() &&
        theParticleDefinitions[i] != outgoingRhadron && theParticleDefinitions[i] != outgoingTarget) {
      G4ReactionProduct* pa = new G4ReactionProduct;
      pa->SetDefinition(theParticleDefinitions[i]);
      (G4UniformRand() < 0.5) ? pa->SetSide(-1) : pa->SetSide(1); //Here we randomly determine the hemisphere of the secondary particle
      vec.SetElement(vecLen++, pa);
    }
  }

  //Update the current and target particles based on wether or not they survive the reaction
  if (!IncidentSurvives) {
    currentParticle.SetDefinitionAndUpdateE(outgoingCloud);
    modifiedOriginal.SetDefinition(outgoingCloud);
  }
  if (!TargetSurvives)
    targetParticle.SetDefinitionAndUpdateE(outgoingTarget);

  CalculateMomenta(vec,
                   vecLen,
                   originalIncident,
                   originalTarget,
                   modifiedOriginal,
                   targetNucleus,
                   currentParticle,
                   targetParticle,
                   !IncidentSurvives,
                   !TargetSurvives,
                   quasiElastic);

 //STOPPED HERE//////////////

  G4String cPname = currentParticle.GetDefinition()->GetParticleName();

  aParticleChange.SetNumberOfSecondaries(vecLen + NumberOfSecondaries);
  G4double e_kin = 0;
  G4LorentzVector cloud_p4_new;  //Cloud 4-momentum in lab after collision

  cloud_p4_new.setVectM(currentParticle.GetMomentum(), currentParticle.GetMass());
  cloud_p4_new *= cloudParticleToLabFrameRotation;

  G4LorentzVector cloud_p4_old_full = Cloud4Momentum;  //quark system in CMS BEFORE collision
  cloud_p4_old_full.boost(trafo_full_cms);
  G4LorentzVector cloud_p4_old_cloud = Cloud4Momentum;  //quark system in cloud CMS BEFORE collision
  cloud_p4_old_cloud.boost(trafo_cloud_cms);
  G4LorentzVector cloud_p4_full = cloud_p4_new;  //quark system in CMS AFTER collision
  cloud_p4_full.boost(trafo_full_cms);
  G4LorentzVector cloud_p4_cloud = cloud_p4_new;  //quark system in cloud CMS AFTER collision
  cloud_p4_cloud.boost(trafo_cloud_cms);

  G4LorentzVector p_g_cms = gluinoMomentum;  //gluino in CMS BEFORE collision
  p_g_cms.boost(trafo_full_cms);

  G4LorentzVector p4_new = cloud_p4_new + gluinoMomentum;
  //  G4cout<<"P4-diff: "<<(p4_new-cloud_p4_new-gluinoMomentum)/GeV<<", magnitude: "
  // <<(p4_new-cloud_p4_new-gluinoMomentum).m()/MeV<<" MeV" <<G4endl;

  G4ThreeVector p_new = p4_new.vect();

  aParticleChange.ProposeLocalEnergyDeposit((p4_new - cloud_p4_new - gluinoMomentum).m());

  if (!IncidentSurvives) {
    G4DynamicParticle* p0 = new G4DynamicParticle;
    p0->SetDefinition(outgoingRhadron);
    p0->SetMomentum(p_new);

    // May need to run SetDefinitionAndUpdateE here...
    G4Track* Track0 = new G4Track(p0, aTrack.GetGlobalTime(), aPosition);
    Track0->SetTouchableHandle(thisTouchable);
    aParticleChange.AddSecondary(Track0);
    /*
      G4cout<<"Adding a particle "<<p0->GetDefinition()->GetParticleName()<<G4endl;
      G4cout<<"with momentum: "<<p0->GetMomentum()/GeV<<" GeV"<<G4endl;
      G4cout<<"and kinetic energy: "<<p0->GetKineticEnergy()/GeV<<" GeV"<<G4endl;
      */
    if (p0->GetKineticEnergy() > e_kin_0) {
      G4cout << "ALAAAAARM!!! (incident changed from " << IncidentRhadron->GetDefinition()->GetParticleName() << " to "
             << p0->GetDefinition()->GetParticleName() << ")" << G4endl;
      G4cout << "Energy loss: " << (e_kin_0 - p0->GetKineticEnergy()) / GeV << " GeV (should be positive)" << G4endl;
      //Turns out problem is only in 2 -> 3 (Won't fix 2 -> 2 energy deposition)
      if (reactionProduct.size() != 3)
        G4cout << "DOUBLE ALAAAAARM!!!" << G4endl;
    } /*else {
	G4cout<<"NO ALAAAAARM!!!"<<G4endl;
	}*/
    if (std::abs(p0->GetKineticEnergy() - e_kin_0) > 100 * GeV) {
      G4cout << "Diff. too big" << G4endl;
    }
    aParticleChange.ProposeTrackStatus(fStopAndKill);
  } else {
    G4double p = p_new.mag();
    if (p > DBL_MIN)
      aParticleChange.ProposeMomentumDirection(p_new.x() / p, p_new.y() / p, p_new.z() / p);
    else
      aParticleChange.ProposeMomentumDirection(1.0, 0.0, 0.0);
  }

  //    return G4VDiscreteProcess::PostStepDoIt( aTrack, aStep);
  if (targetParticle.GetMass() > 0.0)  // targetParticle can be eliminated in TwoBody
  {
    G4DynamicParticle* p1 = new G4DynamicParticle;
    p1->SetDefinition(targetParticle.GetDefinition());
    //G4cout<<"Target secondary: "<<targetParticle.GetDefinition()->GetParticleName()<<G4endl;
    G4ThreeVector momentum = targetParticle.GetMomentum();
    momentum = momentum.rotate(cache, what);
    p1->SetMomentum(momentum);
    p1->SetMomentum((cloudParticleToLabFrameRotation * p1->Get4Momentum()).vect());
    G4Track* Track1 = new G4Track(p1, aTrack.GetGlobalTime(), aPosition);
    Track1->SetTouchableHandle(thisTouchable);
    aParticleChange.AddSecondary(Track1);
  }
  G4DynamicParticle* pa;
  /*
    G4cout<<"vecLen: "<<vecLen<<G4endl;
    G4cout<<"#sec's: "<<aParticleChange.GetNumberOfSecondaries()<<G4endl;
  */

  for (int i = 0; i < vecLen; ++i) {
    pa = new G4DynamicParticle();
    pa->SetDefinition(vec[i]->GetDefinition());
    pa->SetMomentum(vec[i]->GetMomentum());
    pa->Set4Momentum(cloudParticleToLabFrameRotation * (pa->Get4Momentum()));
    G4ThreeVector pvec = pa->GetMomentum();
    G4Track* Trackn = new G4Track(pa, aTrack.GetGlobalTime(), aPosition);
    Trackn->SetTouchableHandle(thisTouchable);
    aParticleChange.AddSecondary(Trackn);

    delete vec[i];
  }

  // Histogram filling
  const G4DynamicParticle* theRhadron = FindRhadron(&aParticleChange);

  if (theRhadron != nullptr || IncidentSurvives) {
    double E_new;
    if (IncidentSurvives) {
      E_new = e_kin;
    } else {
      E_new = theRhadron->GetKineticEnergy();
      if (CustomPDGParser::s_isRMeson(theRhadron->GetDefinition()->GetPDGEncoding()) !=
              CustomPDGParser::s_isRMeson(theIncidentPDG) ||
          CustomPDGParser::s_isMesonino(theRhadron->GetDefinition()->GetPDGEncoding()) !=
              CustomPDGParser::s_isMesonino(theIncidentPDG)) {
        G4cout << "Rm: " << CustomPDGParser::s_isRMeson(theRhadron->GetDefinition()->GetPDGEncoding())
               << " vs: " << CustomPDGParser::s_isRMeson(theIncidentPDG) << G4endl;
        G4cout << "Sm: " << CustomPDGParser::s_isMesonino(theRhadron->GetDefinition()->GetPDGEncoding())
               << " vs: " << CustomPDGParser::s_isMesonino(theIncidentPDG) << G4endl;
      }
    }

    //Calculating relevant scattering angles.
    G4LorentzVector p4_old_full = FullRhadron4Momentum;  //R-hadron in CMS BEFORE collision
    p4_old_full.boost(trafo_full_cms);
    G4LorentzVector p4_old_cloud = FullRhadron4Momentum;  //R-hadron in cloud CMS BEFORE collision
    p4_old_cloud.boost(trafo_cloud_cms);
    G4LorentzVector p4_full = p4_new;  //R-hadron in CMS AFTER collision
    //      G4cout<<p4_full.v()/GeV<<G4endl;
    p4_full = p4_full.boost(trafo_full_cms);
    // G4cout<<p4_full.m()<<" / "<<(cloud_p4_new+gluinoMomentum).boost(trafo_full_cms).m()<<G4endl;
    G4LorentzVector p4_cloud = p4_new;  //R-hadron in cloud CMS AFTER collision
    p4_cloud.boost(trafo_cloud_cms);

    G4double AbsDeltaE = incidentRHadronKineticEnergy - E_new;
    //      G4cout <<"Energy loss: "<<AbsDeltaE/GeV<<G4endl;
    if (AbsDeltaE > 10 * GeV) {
      G4cout << "Energy loss larger than 10 GeV..." << G4endl;
      G4cout << "Incident RHadron kinetic energy: " << E_incidentRHadronKineticEnergy0 / GeV << " GeV" << G4endl;
      G4cout << "E_new: " << E_new / GeV << " GeV" << G4endl;
      G4cout << "Gamma: " << IncidentRhadron->GetTotalEnergy() / IncidentRhadron->GetDefinition()->GetPDGMass()
             << G4endl;
      G4cout << "x: " << aPosition.x() / cm << " cm" << G4endl;
    }
  }
  delete originalIncident;
  delete originalTarget;
  //  aParticleChange.DumpInfo();
  //  G4cout << "Exiting FullModelHadronicProcess::PostStepDoIt"<<G4endl;

  //clear interaction length
  ClearNumberOfInteractionLengthLeft();

  return &aParticleChange;
}

void FullModelHadronicProcess::CalculateMomenta(
    G4FastVector<G4ReactionProduct, MYGHADLISTSIZE>& vec,
    G4int& vecLen,
    const G4HadProjectile* originalIncident,  // the original incident particle
    const G4DynamicParticle* originalTarget,
    G4ReactionProduct& modifiedOriginal,  // Fermi motion and evap. effects included
    G4Nucleus& targetNucleus,
    G4ReactionProduct& currentParticle,
    G4ReactionProduct& targetParticle,
    G4bool& incidentHasChanged,
    G4bool& targetHasChanged,
    G4bool quasiElastic) {
  FullModelReactionDynamics theReactionDynamics;

  cache = 0;
  what = originalIncident->Get4Momentum().vect();

  if (quasiElastic) {
    //      G4cout<<"We are calling TwoBody..."<<G4endl;
    theReactionDynamics.TwoBody(
        vec, vecLen, modifiedOriginal, originalTarget, currentParticle, targetParticle, targetNucleus, targetHasChanged);

    return;
  }

  //If ProduceStrangeParticlePairs is commented out, let's cut this one as well
  G4ReactionProduct leadingStrangeParticle;
  G4bool leadFlag = MarkLeadingStrangeParticle(currentParticle, targetParticle, leadingStrangeParticle);

  //
  // Note: the number of secondaries can be reduced in GenerateXandPt and TwoCluster
  //
  G4bool finishedGenXPt = false;
  G4bool annihilation = false;
  if (originalIncident->GetDefinition()->GetPDGEncoding() < 0 && currentParticle.GetMass() == 0.0 &&
      targetParticle.GetMass() == 0.0) {
    // original was an anti-particle and annihilation has taken place
    annihilation = true;
    G4double ekcor = 1.0;
    G4double cloudKineticEnergy = originalIncident->GetKineticEnergy();
    G4double ekOrg = cloudKineticEnergy;

    const G4double tarmas = originalTarget->GetDefinition()->GetPDGMass();
    if (cloudKineticEnergy > 1.0 * GeV)
      ekcor = 1. / (cloudKineticEnergy / GeV);
    const G4double atomicWeight = G4double(targetNucleus.GetN_asInt());
    cloudKineticEnergy = 2 * tarmas + cloudKineticEnergy * (1. + ekcor / atomicWeight);

    G4double targetNucleusKineticEnergy = targetNucleus.Cinema(cloudKineticEnergy);
    ekOrg += targetNucleusKineticEnergy;
    modifiedOriginal.SetKineticEnergy(ekOrg);
  }

  const G4double twsup[] = {1.0, 0.7, 0.5, 0.3, 0.2, 0.1};
  G4double rand1 = G4UniformRand();
  G4double rand2 = G4UniformRand();
  if ((annihilation || (vecLen >= 6) || (modifiedOriginal.GetKineticEnergy() / GeV >= 1.0)) &&
      (((originalIncident->GetDefinition() == G4KaonPlus::KaonPlus()) ||
        (originalIncident->GetDefinition() == G4KaonMinus::KaonMinus()) ||
        (originalIncident->GetDefinition() == G4KaonZeroLong::KaonZeroLong()) ||
        (originalIncident->GetDefinition() == G4KaonZeroShort::KaonZeroShort())) &&
       ((rand1 < 0.5) || (rand2 > twsup[vecLen]))))
    finishedGenXPt = theReactionDynamics.GenerateXandPt(vec,
                                                        vecLen,
                                                        modifiedOriginal,
                                                        originalIncident,
                                                        currentParticle,
                                                        targetParticle,
                                                        targetNucleus,
                                                        incidentHasChanged,
                                                        targetHasChanged,
                                                        leadFlag,
                                                        leadingStrangeParticle);
  if (finishedGenXPt) {
    Rotate(vec, vecLen);
    return;
  }

  G4bool finishedTwoClu = false;
  if (modifiedOriginal.GetTotalMomentum() / MeV < 1.0) {
    for (G4int i = 0; i < vecLen; i++)
      delete vec[i];
    vecLen = 0;
  } else {
    theReactionDynamics.SuppressChargedPions(vec,
                                             vecLen,
                                             modifiedOriginal,
                                             currentParticle,
                                             targetParticle,
                                             targetNucleus,
                                             incidentHasChanged,
                                             targetHasChanged);

    try {
      finishedTwoClu = theReactionDynamics.TwoCluster(vec,
                                                      vecLen,
                                                      modifiedOriginal,
                                                      originalIncident,
                                                      currentParticle,
                                                      targetParticle,
                                                      targetNucleus,
                                                      incidentHasChanged,
                                                      targetHasChanged,
                                                      leadFlag,
                                                      leadingStrangeParticle);
    } catch (G4HadReentrentException& aC) {
      aC.Report(G4cout);
      throw G4HadReentrentException(__FILE__, __LINE__, "Failing to calculate momenta");
    }
  }
  if (finishedTwoClu) {
    Rotate(vec, vecLen);
    return;
  }

  //
  // PNBlackTrackEnergy is the kinetic energy available for
  //   proton/neutron black track particles [was enp(1) in fortran code]
  // DTABlackTrackEnergy is the kinetic energy available for
  //   deuteron/triton/alpha particles      [was enp(3) in fortran code]
  // the atomic weight of the target nucleus is >= 1.5            AND
  //   neither the incident nor the target particles have changed  AND
  //     there is no kinetic energy available for either proton/neutron
  //     or for deuteron/triton/alpha black track particles
  // For diffraction scattering on heavy nuclei use elastic routines instead

  theReactionDynamics.TwoBody(
      vec, vecLen, modifiedOriginal, originalTarget, currentParticle, targetParticle, targetNucleus, targetHasChanged);
}

G4bool FullModelHadronicProcess::MarkLeadingStrangeParticle(const G4ReactionProduct& currentParticle,
                                                            const G4ReactionProduct& targetParticle,
                                                            G4ReactionProduct& leadParticle) {
  // the following was in GenerateXandPt and TwoCluster
  // add a parameter to the GenerateXandPt function telling it about the strange particle
  //
  // assumes that the original particle was a strange particle
  //
  G4bool lead = false;
  if ((currentParticle.GetMass() >= G4KaonPlus::KaonPlus()->GetPDGMass()) &&
      (currentParticle.GetDefinition() != G4Proton::Proton()) &&
      (currentParticle.GetDefinition() != G4Neutron::Neutron())) {
    lead = true;
    leadParticle = currentParticle;  //  set lead to the incident particle
  } else if ((targetParticle.GetMass() >= G4KaonPlus::KaonPlus()->GetPDGMass()) &&
             (targetParticle.GetDefinition() != G4Proton::Proton()) &&
             (targetParticle.GetDefinition() != G4Neutron::Neutron())) {
    lead = true;
    leadParticle = targetParticle;  //   set lead to the target particle
  }
  return lead;
}

void FullModelHadronicProcess::Rotate(G4FastVector<G4ReactionProduct, MYGHADLISTSIZE>& vec, G4int& vecLen) {
  G4double rotation = 2. * pi * G4UniformRand();
  cache = rotation;
  G4int i;
  for (i = 0; i < vecLen; ++i) {
    G4ThreeVector momentum = vec[i]->GetMomentum();
    momentum = momentum.rotate(rotation, what);
    vec[i]->SetMomentum(momentum);
  }
}

const G4DynamicParticle* FullModelHadronicProcess::FindRhadron(G4ParticleChange* aParticleChange) {
  G4int nsec = aParticleChange->GetNumberOfSecondaries();
  if (nsec == 0)
    return nullptr;
  int i = 0;
  G4bool found = false;
  while (i != nsec && !found) {
    //    G4cout<<"Checking "<<aParticleChange->GetSecondary(i)->GetDynamicParticle()->GetDefinition()->GetParticleName()<<G4endl;
    //    if (aParticleChange->GetSecondary(i)->GetDynamicParticle()->GetDefinition()->GetParticleType()=="rhadron") found = true;
    if (dynamic_cast<CustomParticle*>(aParticleChange->GetSecondary(i)->GetDynamicParticle()->GetDefinition()) !=
        nullptr)
      found = true;
    i++;
  }
  i--;
  if (found)
    return aParticleChange->GetSecondary(i)->GetDynamicParticle();
  return nullptr;
}
