// -*- C++ -*-
//
// Package:    SpikedRHadronAnalyzer
// Class:      SpikedRHadronAnalyzer
//
/**\class SpikedRHadronAnalyzer src/SpikedRHadronAnalyzer.cc

 Description: [Reads SIM ROOTs with gluino samples and outputs histograms for: calorimiter hit location, 4 momenta of R-Hadrons 1 and 2, ...]

 Implementation:
     [Github repository: https://github.com/ctdax/SUSYBSMAnalysis]
*/
//
// Original Author:  Colby Thompson
//         Created:  Thu, 22 Feb 2024 20:45:33 GMT
//
//


//System include files
#include <memory>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>

//Triggers and Handles
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

//Tracker
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "SimDataFormats/Track/interface/CoreSimTrack.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/CoreSimVertex.h"

//ECAL
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

//HCAL
#include "Geometry/Records/interface/HcalGeometryRecord.h"
#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "Geometry/CaloTopology/interface/HcalTopology.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/HcalCommonData/interface/HcalDDDRecConstants.h"
#include "Geometry/HcalCommonData/interface/HcalHitRelabeller.h"
#include "Geometry/Records/interface/HcalRecNumberingRecord.h"
#include "CondFormats/HcalObjects/interface/HcalRespCorrs.h"
#include "CondFormats/DataRecord/interface/HcalRespCorrsRcd.h"

//Muon Chamber
#include "Geometry/Records/interface/MuonGeometryRcd.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/DTGeometry/interface/DTLayer.h"
#include "Geometry/DTGeometry/interface/DTChamber.h"
#include "Geometry/DTGeometry/interface/DTSuperLayer.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"

//Calculations
#include "DataFormats/Common/interface/HLTPathStatus.h"
#include "DataFormats/Common/interface/HLTenums.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackReco/interface/DeDxHitInfo.h"

//ROOT
#include "TFile.h"
#include "TH1.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TString.h"
#include "TCanvas.h"
#include "TEfficiency.h"
#include "TStyle.h"
#include "TObject.h"
#include "TDirectory.h"
#include "TMatrix.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TGraph.h"

//HSCP specific packages
#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCParticle.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/MuonSegment.h"
//#include "SUSYBSMAnalysis/Analyzer/interface/CommonFunction.h"
//#include "SUSYBSMAnalysis/Analyzer/interface/DeDxUtility.h"
//#include "SUSYBSMAnalysis/HSCP/interface/HSCPHelpers.h"

//FWCORE
#define FWCORE

// Thresholds for candidate preselection
float globalMaxEta_ = 1.0;
float globalMinPt_ = 55.0;
unsigned int globalMinNOPH_ = 2;
float globalMinFOVH_ = 0.8;
unsigned int globalMinNOM_ = 10;
float globalMaxChi2_ = 5.0;
float globalMaxEoP_ = 0.3; //unsure
float globalMaxDZ_ = 0.1;
float globalMaxDXY_ = 0.02;
float globalMaxTIsol_, globalMiniRelIsoAll_; //unsure
float globalMinIh_ = 3.47; //unsure
float trackProbQCut_;
unsigned int minMuStations_;
float globalMinIs_ = 0.0;
float globalMinTOF_;
float GlobalMinNDOF = 8;            // cut on number of     DegreeOfFreedom used for muon TOF measurement
float GlobalMinNDOFDT = 6;          // cut on number of DT  DegreeOfFreedom used for muon TOF measurement
float GlobalMinNDOFCSC = 6;         // cut on number of CSC DegreeOfFreedom used for muon TOF measurement
float GlobalMaxTOFErr = 0.15;       //0.07;   // cut on error on muon TOF measurement
bool useClusterCleaning = true; //unsure

int evtcount = 0;
int evtcount00 = 0;
int evtcount01 = 0;
int evtcount01n = 0;
int evtcount01c = 0;
int evtcount11 = 0;
int badiqstat = 0;
int passevtcount0 = 0;
int passevtcount1 = 0;
int passevtcount00 = 0;
int passevtcount01 = 0;
int passevtcount01c = 0;
int passevtcount01n = 0;
int passevtcount11 = 0;

int ngencheta09 = 0;
int nsimmatch[30] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; 
int nsimnomatch[30] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; 
int nsimnomatchqsame[30] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; 
int nsimnomatchqdiff[30] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; 
double layerr[30] = {0, 3., 7., 11., 16., 26., 33.5, 42., 50., 61., 69., 78., 87., 97., 108., 440., 525., 630., 740., 0., 0., 0.};

struct simhitinfo {
  int ilayer;
  int pdgID;
//  PSimHit *hit;
};

using namespace edm;
using namespace reco;
using namespace std;
using namespace __gnu_cxx;
using namespace trigger;

class TupleMaker;
class MCWeight;

#define DRMATCH 0.002
#define DPTFRACMIN -0.5
#define DPTFRACMAX 1.0

#define TIBEVT1 2
#define TIBEVT2 3
#define TIBEVT3 5
#define TRKEVT1 11
#define TRKEVT2 13
#define TRKEVT3 15
#define BPIXEVT1 32
#define BPIXEVT2 34
#define BPIXEVT3 27

#define NTRLAYERS 19
#define BPIX0 0
#define TIB0 4
#define TOB0 8
#define DT0 14

#define TYPEGLUE 9 
#define TYPEMESON 111
#define TYPEBARYON 1111
#define TYPELEPTON 11

#define ETACUT 0.9

class SpikedRHadronAnalyzer : public edm::EDAnalyzer {
public:
  explicit SpikedRHadronAnalyzer (const edm::ParameterSet&);
  ~SpikedRHadronAnalyzer();


private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

  float deltaphi(float phi1, float phi2);
  float deltaphiabs(float phi1, float phi2);
  float deltaRta(float phi1, float eta1, float phi2, float eta2);
  int determineCharge(int pid);
  int determineType(int pid);
  bool passHSCPPreselection(int typeMode, susybsm::HSCParticle hscpCan, reco::VertexCollection vertColl, reco::VertexCollection inclusiveSecondaryVertices, reco::PFCandidateCollection pf, reco::DeDxData *dedxSObj, reco::DeDxData *dedxMObj, bool passedCutsArray[]);
  bool comparePID(pair<int, int> p1, pair<int, int> p2); 

  edm::EDGetTokenT<vector<reco::GenParticle>> genParticlesToken_;

  const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> tGeomEsToken_;

  // Tracker hits
  edm::EDGetTokenT<edm::SimTrackContainer> edmSimTrackContainerToken_;
  edm::EDGetTokenT<edm::SimVertexContainer> edmSimVertexContainerToken_;
  edm::EDGetTokenT<edm::PSimHitContainer> edmPSimHitContainer_siTIBLow_Token_;
  edm::EDGetTokenT<edm::PSimHitContainer> edmPSimHitContainer_siTIBHigh_Token_;
  edm::EDGetTokenT<edm::PSimHitContainer> edmPSimHitContainer_siTOBLow_Token_;
  edm::EDGetTokenT<edm::PSimHitContainer> edmPSimHitContainer_siTOBHigh_Token_;
  edm::EDGetTokenT<edm::PSimHitContainer> edmPSimHitContainer_siTIDLow_Token_;
  edm::EDGetTokenT<edm::PSimHitContainer> edmPSimHitContainer_siTIDHigh_Token_;
  edm::EDGetTokenT<edm::PSimHitContainer> edmPSimHitContainer_siTECLow_Token_;
  edm::EDGetTokenT <edm::PSimHitContainer> edmPSimHitContainer_siTECHigh_Token_;
  edm::EDGetTokenT <edm::PSimHitContainer> edmPSimHitContainer_pxlBrlLow_Token_;
  edm::EDGetTokenT <edm::PSimHitContainer> edmPSimHitContainer_pxlBrlHigh_Token_;
  edm::EDGetTokenT <edm::PSimHitContainer> edmPSimHitContainer_pxlFwdLow_Token_;
  edm::EDGetTokenT <edm::PSimHitContainer> edmPSimHitContainer_pxlFwdHigh_Token_;

  // ECAL hits
  edm::EDGetTokenT <edm::PCaloHitContainer> edmCaloHitContainer_EcalHitsEB_Token_;
  edm::EDGetTokenT <edm::PCaloHitContainer> edmCaloHitContainer_EcalHitsEE_Token_;
  edm::EDGetTokenT <edm::PCaloHitContainer> edmCaloHitContainer_EcalHitsES_Token_;

  // HCAL hits
  const HcalDDDRecConstants *hcons_         ;
  const CaloGeometry        *geometry_      ;
  const CaloTopology        *topology_      ;
  const HcalTopology        *hcalTopology_  ;
  const HcalRespCorrs       *respCorrs_     ;
  HcalRespCorrs             *respCorrs      ;

  const edm::ESGetToken< HcalDDDRecConstants, HcalRecNumberingRecord > tok_HRNDC_        ;
  const edm::ESGetToken< CaloGeometry,        CaloGeometryRecord     > tok_geom_         ;
  const edm::ESGetToken< CaloTopology,        CaloTopologyRecord     > tok_caloTopology_ ;
  const edm::ESGetToken< HcalTopology,        HcalRecNumberingRecord > tok_hcalTopology_ ;
  const edm::ESGetToken< HcalRespCorrs,       HcalRespCorrsRcd       > tok_resp_         ;
  edm::EDGetTokenT <edm::PCaloHitContainer> edmCaloHitContainer_HcalHits_Token_;


  // Muon hits
  edm::EDGetTokenT <edm::PSimHitContainer> edmPSimHitContainer_muonCSC_Token_;
  edm::EDGetTokenT <edm::PSimHitContainer> edmPSimHitContainer_muonDT_Token_;
  edm::EDGetTokenT <edm::PSimHitContainer> edmPSimHitContainer_muonRPC_Token_;
  edm::EDGetTokenT <edm::PSimHitContainer> edmPSimHitContainer_muonGEM_Token_;

  // Miscellaneous
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<trigger::TriggerEvent> triggerSummary_;
  edm::EDGetTokenT<std::vector<reco::PFMET>> pfmet_;
  edm::EDGetTokenT<std::vector<reco::CaloMET>> calomet_;
  edm::EDGetTokenT<VertexCollection> vertexToken_;
  edm::EDGetTokenT<VertexCollection> inclusiveSecondaryvertexToken_;
  edm::EDGetTokenT<PFCandidateCollection> pfToken_;
  edm::EDGetTokenT<std::vector<reco::Track>> gentrkToken_;
  edm::EDGetTokenT<vector<susybsm::HSCParticle>> hscpCandToken_;
  edm::EDGetTokenT<reco::DeDxHitInfoAss> dedxToken_;

  //declare variables here
#define TLEN 1625
  string trig_name[TLEN] = {};
  int trig_pass[TLEN]    = {};
  int trig_total[TLEN]   = {};
  int trig_pass_or = 0;
  int trig_pass_and = 0;
  //int trig_total_special = 0;

   std::vector<std::pair<int, int>> genrhadcounts;

  ///double met_values[];

  int typeMode_;
  int sampleType_;
  string period_;
  int isData;
  bool useTemplateLayer_;
  bool skipPixel_;
  string dEdxTemplate_;
  bool  enableDeDxCalibration_;
  string dEdxCalibration_;

  //dedxGainCorrector trackerCorrector;
  //TH3F* dEdxTemplates = nullptr;
  float dEdxSF_0_, dEdxSF_1_;
  float dEdxSF[2] = {dEdxSF_0_, dEdxSF_1_};

  std::ofstream csv;
  TFile *outputFile_;
  TH1F *ECalHits_energy;
  TH1F *ECalHits_eta;
  TH1F *ECalHits_phi;
  TH2F *ECalHits_2DEtaPhi;
  TH2F *ECalHits1000_2DEtaPhi;
  TH2F *ECalHits_2DXY;
  TH2F *EBHits_RZ;
  TH2F *ECalHits1000_2DXY;
  TH2F *EBHits1000_RZ;
  TH1F *ECalHits_x;
  TH1F *ECalHits_y;
  TH1F *ECalHits_z;
  TH1F *HCalHits_energy;
  TH1F *HCalHits_eta;
  TH1F *HCalHits_phi;
  TH2F *HCalHits_2DEtaPhi;
  TH1F *HCalHits_x;
  TH1F *HCalHits_y;
  TH1F *HCalHits_z;
  TH1F *RHadron1_pdgId;
  TH1F *RHadron1_mass;
  TH1F *RHadron1_px;
  TH1F *RHadron1_py;
  TH1F *RHadron1_pz;
  TH2F *RHadron1_pxpy;
  TH1F *RHadron2_pdgId;
  TH1F *RHadron2_mass;
  TH1F *RHadron2_px;
  TH1F *RHadron2_py;
  TH1F *RHadron2_pz;
  TH2F *RHadron2_pxpy;
  TH1F *RHadron_ECalSubDetector;
  TH2F *SimVertex_2DXY;
  TH2F *SimVertex_2DRZ;

};

std::string vectorToString(const std::vector<int>& vec, const std::string& delimiter = "/") {
    std::ostringstream oss;
    for (size_t i = 0; i < vec.size(); ++i) {
        oss << vec[i];
        if (i != vec.size() - 1) {
            oss << delimiter;
        }
    }
    return oss.str();
}

int SpikedRHadronAnalyzer::determineCharge(int pdgID) {
  int charge = -99;
  int pid = abs(pdgID);
  if (pid==1000993) {charge = 0;} // ~g_glueball
  else if (pid==1009113) {charge = 0;} // ~g_rho0
  else if (pid==1009223) {charge = 0;} // ~g_omega
  else if (pid==1009213) {charge = 1;} // ~g_rho+
  else if (pid==1009313) {charge = 0;} // ~g_K*0
  else if (pid==1009323) {charge = 1;} // ~g_K*+
  else if (pid==1009333) {charge = 0;} // ~g_phi
  else if (pid==1091114) {charge = -1;} // ~g_Delta-
  else if (pid==1092114) {charge = 0;} // ~g_Delta0
  else if (pid==1092214) {charge = 1;} // ~g_Delta+
  else if (pid==1092224) {charge = 2;} // ~g_Delta++
  else if (pid==1093114) {charge = -1;} // ~g_Sigma*-
  else if (pid==1093214) {charge = 0;} // ~g_Sigma*0
  else if (pid==1093224) {charge = 1;} // ~g_Sigma*+
  else if (pid==1093314) {charge = -1;} // ~g_Xi*-
  else if (pid==1093324) {charge = 0;} // ~g_Xi*0
  else if (pid==1093334) {charge = -1;} // ~g_Omega-
  else if (pid==1000015) {charge = -1;} // ~stau-
  else if (pid==1000612) {charge = +1;} // ~T+
  else if (pid==1000622 ) {charge = 0;} // ~T0 
  else if (pid==1000632 ) {charge = +1;} // ~T_s+
  else if (pid==1000642 ) {charge = 0;} // ~T_c0
  else if (pid==1000652 ) {charge = +1;} // ~T_b+
  else if (pid==1006113 ) {charge = 0;} // ~T_dd10
  else if (pid==1006211 ) {charge = +1;} // ~T_ud0+
  else if (pid==1006213 ) {charge = +1;} // ~T_ud1+
  else if (pid==1006223 ) {charge = +2;} // ~T_uu1++ 
  else if (pid==1006311 ) {charge = 0;} // ~T_sd00  
  else if (pid==1006313 ) {charge = 0;} // ~T_sd10  
  else if (pid==1006321 ) {charge = +1;} // ~T_su0+  
  else if (pid==1006323 ) {charge = +1;} // ~T_su1+  
  else if (pid==1006333 ) {charge = 0;} // ~T_ss10  
/* 1000021   1500.0     # ~g
      1000993   1500.700   # ~g_glueball
      1009213   1500.650   # ~g_rho+
      1009313   1500.825   # ~g_K*0
      1009323   1500.825   # ~g_K*+
      1009113   1500.650   # ~g_rho0
      1009223   1500.650   # ~g_omega
      1009333   1501.800   # ~g_phi
      1091114   1500.975   # ~g_Delta-
      1092114   1500.975   # ~g_Delta0
      1092214   1500.975   # ~g_Delta+
      1092224   1500.975   # ~g_Delta++
      1093114   1501.150   # ~g_Sigma*-
      1093214   1501.150   # ~g_Sigma*0
      1093224   1501.150   # ~g_Sigma*+
      1093314   1501.300   # ~g_Xi*-
      1093324   1501.300   # ~g_Xi*0
      1093334   1501.600   # ~g_Omega- 
  stop
        1000006 800.000   # ~t_1
	1000612 800.325   # ~T+  
	1000622 800.325   # ~T0  
	1000632 800.500   # ~T_s+
	1000642 801.500   # ~T_c0
	1000652 804.800   # ~T_b+
	1006113 800.650   # ~T_dd10
	1006211 800.650   # ~T_ud0+
	1006213 800.650   # ~T_ud1+
	1006223 800.650   # ~T_uu1++ 
	1006311 800.825   # ~T_sd00  
	1006313 800.825   # ~T_sd10  
	1006321 800.825   # ~T_su0+  
	1006323 800.825   # ~T_su1+  
	1006333 801.000   # ~T_ss10  
        -1000006 800.000   # ~t_1bar
	-1000612 800.325   # ~Tbar-  
	-1000622 800.325   # ~Tbar0  
	-1000632 800.500   # ~Tbar_s-
	-1000642 801.500   # ~Tbar_c0
	-1000652 804.800   # ~Tbar_b-
	-1006113 800.650   # ~Tbar_dd10
	-1006211 800.650   # ~Tbar_ud0-
	-1006213 800.650   # ~Tbar_ud1-
	-1006223 800.650   # ~Tbar_uu1-- 
	-1006311 800.825   # ~Tbar_sd00  
	-1006313 800.825   # ~Tbar_sd10  
	-1006321 800.825   # ~Tbar_su0-  
	-1006323 800.825   # ~Tbar_su1-  
	-1006333 801.000   # ~Tbar_ss10  

*/
  //else {cout << "Unknown PID = " << pid << endl;}
  if ((charge>-98)&&(pdgID<0)) charge *= -1;
  return(charge);
}

int SpikedRHadronAnalyzer::determineType(int pdgID) {
     int pidabs = abs(pdgID);
     if ((pidabs>1000900)&&(pidabs<1000999)) {
       return(TYPEGLUE);
     } else if ((pidabs>1000999)&&(pidabs<1009999)) {
       return(TYPEMESON);
     } else if ((pidabs>1009999)&&(pidabs<1099999)) {
       return(TYPEBARYON);
     } else if ((pidabs>1000009)&&(pidabs<1000020)) {
       return(TYPELEPTON);
     }
  return(-1);
}

// delta(phi1-phi2) - return value from -pi to +pi
float SpikedRHadronAnalyzer::deltaphi(float phi1, float phi2) {
  float tdelta = phi1 - phi2;
  while (tdelta<(-1.0*M_PI)) {tdelta += (2.0*M_PI);}
  while (tdelta>(1.0*M_PI)) {tdelta -= (2.0*M_PI);}
  return(tdelta);
}

// |delta(phi1-phi2)| - return value from 0 to +pi
float SpikedRHadronAnalyzer::deltaphiabs(float phi1, float phi2) {
  float tdphi = deltaphi(phi1,phi2);
  return(fabs(tdphi));
}

float SpikedRHadronAnalyzer::deltaRta(float phi1, float eta1, float phi2, float eta2) {
  float diffphi1 = deltaphiabs(phi1,phi2);
  float diffeta1 = eta1 - eta2;
  float tmpdeltaR1 = sqrt(diffphi1*diffphi1 + diffeta1*diffeta1);
  return(tmpdeltaR1);
}

// preselection
bool SpikedRHadronAnalyzer::passHSCPPreselection(int typeMode, susybsm::HSCParticle hscpCan, reco::VertexCollection vertexColl, reco::VertexCollection inclusiveSecondaryVertices, reco::PFCandidateCollection pf, reco::DeDxData *dedxSObj, reco::DeDxData *dedxMObj, bool passedCutsArray[]) {

  reco::TrackRef track = hscpCan.trackRef();
  if (track.isNull()) return(false);

    // find closest vertex
  int closestZGoodVertex = -1;
  int goodVerts = 0;
  float dzMin = 10000;
    // Loop on the vertices in the event
  for (unsigned int i = 0; i < vertexColl.size(); i++) {
    if (vertexColl[i].isFake() || fabs(vertexColl[i].z()) > 24 || vertexColl[i].position().rho() > 2 || vertexColl[i].ndof() <= 4) continue;  //only consider good vertex
    goodVerts++;
    
    if (fabs(track->dz(vertexColl[i].position())) < fabs(dzMin)) {
      dzMin = fabs(track->dz(vertexColl[i].position()));
      closestZGoodVertex = i;
    }
  } // End loop on the vertices in the event
  
  if (closestZGoodVertex < 0) {
    closestZGoodVertex = 0;
  }
  
  // Impact paramters dz and dxy
  float dz = track->dz(vertexColl[closestZGoodVertex].position());
  float dxy = track->dxy(vertexColl[closestZGoodVertex].position());

  bool isMaterialTrack = false;

  // Loop on the secondary vertices to find tracks that original from the pixel layers
  // i.e. are due to NI
  for (unsigned int i = 0; i < inclusiveSecondaryVertices.size(); i++) {
    if (inclusiveSecondaryVertices[i].isFake()) {
      continue;
    }
    auto rho = inclusiveSecondaryVertices[i].position().rho();
    if ( (( 2.80-0.075 ) < rho && rho < ( 3.10+0.075 )) || (( 6.60-0.075 ) < rho && rho < ( 7.00+0.075 ))
        || (( 10.9-0.075 ) < rho && rho < ( 10.9+0.075 )) || (( 16.0-0.075 ) < rho && rho < ( 16.0+0.075 )) ) {
      for( const auto& rf_track : inclusiveSecondaryVertices[i].refittedTracks() ) {
        const reco::Track& origTrk = *( inclusiveSecondaryVertices[i].originalTrack( rf_track ));
        if( track->pt() == origTrk.pt() ){
          isMaterialTrack = true;
          break;
        } 
      } 
    } else {
      continue;
    } 
  } 

    // Save PF informations and isolation
    float pfIsolation_DZ_ = 0.1;

    float track_PFIso005_sumCharHadPt = 0, track_PFIso005_sumNeutHadPt = 0, track_PFIso005_sumPhotonPt = 0, track_PFIso005_sumPUPt = 0;
    float track_PFIso01_sumCharHadPt = 0, track_PFIso01_sumNeutHadPt = 0, track_PFIso01_sumPhotonPt = 0, track_PFIso01_sumPUPt = 0;
    float track_PFIso03_sumCharHadPt = 0, track_PFIso03_sumNeutHadPt = 0, track_PFIso03_sumPhotonPt = 0, track_PFIso03_sumPUPt = 0;
    float track_PFIso05_sumCharHadPt = 0, track_PFIso05_sumNeutHadPt = 0, track_PFIso05_sumPhotonPt = 0, track_PFIso05_sumPUPt = 0;

    float track_PFMiniIso_sumCharHadPt = 0, track_PFMiniIso_sumNeutHadPt = 0, track_PFMiniIso_sumPhotonPt = 0, track_PFMiniIso_sumPUPt = 0, track_PFMiniIso_sumMuonPt = 0;
    float pf_energy=0;

    float RMin = 9999.;
    unsigned int idx_pf_RMin = 9999;

      for (unsigned int i = 0; i < pf.size(); i++){
          const reco::PFCandidate pfCand = pf[i];
          float dr = deltaR(pfCand.eta(),pfCand.phi(),track->eta(),track->phi());
          if(dr < RMin){
              RMin = dr;
              idx_pf_RMin = i;
           }
       }//end loop PFCandidates

      // https://urldefense.proofpoint.com/v2/url?u=https-3A__github.com_cms-2Dsw_cmssw_blob_72d0fc00976da53d1fb745eb7f37b2a4ad965d7e_&d=DwIGAg&c=gRgGjJ3BkIsb5y6s49QqsA&r=iYf-W5o_XDmJ-u42pi_qsGfw1LmudygYMQ_Az9XJWjc&m=oBFJTDv2gfKrlPynZ9fSMGuSOoJZcNBR8epoz7KKL7P7LrQKFQ2Uhpkq3B3sODWH&s=Y8xRnPaFLOTBqfLjih3mXkPcBp2SP-_OPNc8fWzBgOg&e= 
      // PhysicsTools/PatAlgos/plugins/PATIsolatedTrackProducer.cc#L555
//      for(unsigned int i=0;i<pf->size();i++){
      for (unsigned int i=0;i<pf.size();i++){
        const reco::PFCandidate* pfCand = &(pf)[i];

        pf_energy = pfCand->ecalEnergy() + pfCand->hcalEnergy();
        bool pf_isPhotonForIdx = pfCand->translatePdgIdToType(pfCand->pdgId()) == reco::PFCandidate::ParticleType::gamma;
        bool pf_isChHadronForIdx = pfCand->translatePdgIdToType(pfCand->pdgId()) == reco::PFCandidate::ParticleType::h;
        bool pf_isNeutHadronForIdx = pfCand->translatePdgIdToType(pfCand->pdgId()) == reco::PFCandidate::ParticleType::h0;
        bool pf_isMuonForIdx = pfCand->translatePdgIdToType(pfCand->pdgId()) == reco::PFCandidate::ParticleType::mu;

        if(i == idx_pf_RMin) continue; //don't count itself
        float dr = deltaR(pfCand->eta(),pfCand->phi(),track->eta(),track->phi());
        bool fromPV = (fabs(dz) < pfIsolation_DZ_);
        int id = std::abs(pfCand->pdgId());
        float pt = pfCand->p4().pt();
        if(dr<0.05){
            // charged cands from PV get added to trackIso
          if(id == 211 && fromPV) track_PFIso005_sumCharHadPt+=pt;
            // charged cands not from PV get added to pileup iso
          else if(id == 211) track_PFIso005_sumPUPt+=pt;
            // neutral hadron iso
          if(id == 130) track_PFIso005_sumNeutHadPt+=pt;
            // photon iso
          if(id == 22) track_PFIso005_sumPhotonPt+=pt;
        }if(dr<0.1){
            if(id == 211 && fromPV) track_PFIso01_sumCharHadPt+=pt;
            else if(id == 211) track_PFIso01_sumPUPt+=pt;
            if(id == 130) track_PFIso01_sumNeutHadPt+=pt;
            if(id == 22) track_PFIso01_sumPhotonPt+=pt;
        }if(dr<0.3){
            if(id == 211 && fromPV) track_PFIso03_sumCharHadPt+=pt;
            else if(id == 211) track_PFIso03_sumPUPt+=pt;
            if(id == 130) track_PFIso03_sumNeutHadPt+=pt;
            if(id == 22) track_PFIso03_sumPhotonPt+=pt;
        }if(dr<0.5){
            if(id == 211 && fromPV) track_PFIso05_sumCharHadPt+=pt;
            else if(id == 211) track_PFIso05_sumPUPt+=pt;
            if(id == 130) track_PFIso05_sumNeutHadPt+=pt;
            if(id == 22) track_PFIso05_sumPhotonPt+=pt;
        }

        float drForMiniIso = 0.0;
        if (track->pt() < 50 ) {
          drForMiniIso = 0.2;
        } else if (track->pt() < 200) {
          drForMiniIso = 10/track->pt();
        } else {
          drForMiniIso = 0.05;
        }
        if (dr<drForMiniIso) {
            // charged cands from PV get added to trackIso
          if(pf_isChHadronForIdx && fromPV) track_PFMiniIso_sumCharHadPt+=pt;
            // charged cands not from PV get added to pileup iso
          else if(pf_isChHadronForIdx) track_PFMiniIso_sumPUPt+=pt;
            // neutral hadron iso
          if(pf_isNeutHadronForIdx) track_PFMiniIso_sumNeutHadPt+=pt;
            // photon iso
          if(pf_isPhotonForIdx) track_PFMiniIso_sumPhotonPt+=pt;
            // muon iso
          if(pf_isMuonForIdx) track_PFMiniIso_sumMuonPt+=pt;
        }
      }//end loop PFCandidates
//    }
            

  // Calculate PF mini relative isolation
  float miniRelIsoAll = (track_PFMiniIso_sumCharHadPt + std::max(0.0, track_PFMiniIso_sumNeutHadPt + track_PFMiniIso_sumPhotonPt - 0.5* track_PFMiniIso_sumPUPt))/track->pt();

  float EoP = pf_energy / track->p();

// Number of DeDx hits
  unsigned int numDeDxHits = (dedxSObj) ? dedxSObj->numberOfMeasurements() : 0;

  float Ih = (dedxMObj) ?  dedxMObj->dEdx() : 0.0;

  // No cut, i.e. events after trigger
  passedCutsArray[0]  = true;
  // Check if eta is inside the max eta cut
  passedCutsArray[1]  = (fabs(track->eta()) < globalMaxEta_) ? true : false;
  // Cut on minimum track pT
  passedCutsArray[2]  = (track->pt() > globalMinPt_) ? true : false;
  // Check the number of pixel hits
  passedCutsArray[3]  = (typeMode != 3 && fabs(track->hitPattern().numberOfValidPixelHits()) > globalMinNOPH_) ? true : false;
  // Check the min fraction of valid hits
  passedCutsArray[4]  = (typeMode != 3 && track->validFraction() > globalMinFOVH_) ? true : false;
  // Cut for the number of dEdx hits
  passedCutsArray[5]  = (numDeDxHits >= globalMinNOM_)  ? true : false;
  // Select only high purity tracks
  passedCutsArray[6]  = (typeMode != 3 && track->quality(reco::TrackBase::highPurity)) ? true : false;
  // Cut on the chi2 / ndof
  passedCutsArray[7] = (typeMode != 3 && track->chi2() / track->ndof() < globalMaxChi2_) ? true : false;

  // Cut on the energy over momenta
  passedCutsArray[8] = (EoP < globalMaxEoP_) ? true : false;

 // Cut on the impact parameter
 // for typeMode_ 5 dz is supposed to come from the beamspot, TODO
  passedCutsArray[9] = (  (typeMode != 5 && fabs(dz) < globalMaxDZ_)
                       || (typeMode == 5 && fabs(dz) < 4)) ? true : false;
 // for typeMode_ 5 dxy is supposed to come from the beamspot, TODO
  passedCutsArray[10] = (  (typeMode != 5 && fabs(dxy) < globalMaxDXY_)
                         || (typeMode == 5 && fabs(dxy) < 4)) ? true : false;

 // Cut on the tracker based isolation
  passedCutsArray[12] = (!isMaterialTrack) ? true : false;

 // Cut on the PF based mini-isolation
  passedCutsArray[13] = ( miniRelIsoAll < globalMiniRelIsoAll_) ? true : false;
 // Cut on the PF electron ID

 // Cut on min Ih (or max for fractionally charged)
  passedCutsArray[15] = (  (typeMode != 5 &&  Ih > globalMinIh_)
                         || (typeMode == 5 && Ih < globalMinIh_)) ? true : false;
 
  return(true);
}

bool SpikedRHadronAnalyzer::comparePID(pair<int, int> p1, pair<int, int> p2) {
  return (p1.first < p2.first);
}

//constructor
SpikedRHadronAnalyzer::SpikedRHadronAnalyzer(const edm::ParameterSet& iConfig) 
 :  tGeomEsToken_(esConsumes<TrackerGeometry, TrackerDigiGeometryRecord>())


{
  // Tracker
  edmSimTrackContainerToken_ = consumes<edm::SimTrackContainer>(iConfig.getParameter<edm::InputTag>("G4TrkSrc"));
  edmSimVertexContainerToken_ = consumes<edm::SimVertexContainer>(iConfig.getParameter<edm::InputTag>("G4VtxSrc"));
  edmPSimHitContainer_siTIBLow_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("TrackerHitsTIBLowTof"));
  edmPSimHitContainer_siTIBHigh_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("TrackerHitsTIBHighTof"));
  edmPSimHitContainer_siTOBLow_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("TrackerHitsTOBLowTof"));
  edmPSimHitContainer_siTOBHigh_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("TrackerHitsTOBHighTof"));
  edmPSimHitContainer_siTIDLow_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("TrackerHitsTIDLowTof"));
  edmPSimHitContainer_siTIDHigh_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("TrackerHitsTIDHighTof"));
  edmPSimHitContainer_siTECLow_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("TrackerHitsTECLowTof"));
  edmPSimHitContainer_siTECHigh_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("TrackerHitsTECHighTof"));
  edmPSimHitContainer_pxlBrlLow_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("TrackerHitsPixelBarrelLowTof"));
  edmPSimHitContainer_pxlBrlHigh_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("TrackerHitsPixelBarrelHighTof"));
  edmPSimHitContainer_pxlFwdLow_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("TrackerHitsPixelEndcapLowTof"));
  edmPSimHitContainer_pxlFwdHigh_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("TrackerHitsPixelEndcapHighTof"));

  // Calorimiter
  edmCaloHitContainer_EcalHitsEB_Token_ = consumes<edm::PCaloHitContainer>(iConfig.getParameter<edm::InputTag>("EcalHitsEB"));
  edmCaloHitContainer_EcalHitsEE_Token_ = consumes<edm::PCaloHitContainer>(iConfig.getParameter<edm::InputTag>("EcalHitsEE"));
  edmCaloHitContainer_EcalHitsES_Token_ = consumes<edm::PCaloHitContainer>(iConfig.getParameter<edm::InputTag>("EcalHitsES"));
  edmCaloHitContainer_HcalHits_Token_ = consumes<edm::PCaloHitContainer>(iConfig.getParameter<edm::InputTag>("HcalHits"));

  // Muon Chamber
  edmPSimHitContainer_muonCSC_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("MuonCSCHits"));
  edmPSimHitContainer_muonDT_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("MuonDTHits"));
  edmPSimHitContainer_muonRPC_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("MuonRPCHits"));
  edmPSimHitContainer_muonGEM_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("MuonGEMHits"));

  // Miscellaneous
  triggerBits_ = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"));
  triggerSummary_ = consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag>("trig_sum"));
  pfmet_ = consumes<std::vector<reco::PFMET>>(iConfig.getParameter<edm::InputTag>("pfmet_reco"));
  calomet_ = consumes<std::vector<reco::CaloMET>>(iConfig.getParameter<edm::InputTag>("calomet_reco"));
  vertexToken_ = mayConsume<VertexCollection>(iConfig.getParameter<edm::InputTag>("vertex_reco"));
  inclusiveSecondaryvertexToken_ = mayConsume<VertexCollection>(iConfig.getParameter<edm::InputTag>("secondaryvertex_reco"));
  genParticlesToken_ = mayConsume<vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("gen_info"));
  pfToken_ = mayConsume<PFCandidateCollection>(iConfig.getParameter<edm::InputTag>("pf_reco"));
  gentrkToken_ = mayConsume<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("gentrk_reco"));
  hscpCandToken_ = mayConsume<vector<susybsm::HSCParticle>>(iConfig.getParameter<edm::InputTag>("hscp_cand"));
  dedxToken_ = consumes<reco::DeDxHitInfoAss>(iConfig.getParameter<edm::InputTag>("DedxCollection"));

  typeMode_ = iConfig.getUntrackedParameter<int>("TypeMode");
  sampleType_ = iConfig.getUntrackedParameter<int>("SampleType");
  period_ = iConfig.getUntrackedParameter<string>("Period");
  useTemplateLayer_ = iConfig.getUntrackedParameter<bool>("UseTemplateLayer");
  skipPixel_ = iConfig.getUntrackedParameter<bool>("SkipPixel");
  //dEdxTemplate_ = iConfig.getUntrackedParameter<string>("DeDxTemplate");
  //enableDeDxCalibration_ = iConfig.getUntrackedParameter<bool>("EnableDeDxCalibration");
  //dEdxCalibration_ = iConfig.getUntrackedParameter<string>("DeDxCalibration");
  dEdxSF_0_ = iConfig.getUntrackedParameter<double>("DeDxSF_0");
  dEdxSF_1_ = iConfig.getUntrackedParameter<double>("DeDxSF_1");

  isData = (sampleType_ == 0);
  dEdxSF[0] = dEdxSF_0_;
  dEdxSF[1] = dEdxSF_1_;


  // Output File
  outputFile_ = new TFile("/uscms/home/cthompso/nobackup/CMSSW_10_6_30/src/SUSYBSMAnalysis/HSCP/test/Gluinos1800GeV.root", "RECREATE");  

  // Create csv for energy spike R-hadron analysis
  csv.open ("/uscms/home/cthompso/nobackup/CMSSW_10_6_30/src/SUSYBSMAnalysis/HSCP/test/blank.csv");
  csv << "Event,Calo Hit ID,Calo Hit Energy [GeV],Process Type,Calo Hit Parent Particle,Outgoing Particles From Vertex,Parent To Vertex\n";

  // Declare ROOT histograms
  ECalHits_energy = new TH1F("ECalHits_energy","ECalHits_energy",100,0.,8000.);
  ECalHits_energy->GetXaxis()->SetTitle("[GeV]");
  ECalHits_eta = new TH1F("ECalHits_eta","ECalHits_eta",50,-4.,4.);
  ECalHits_phi = new TH1F("ECalHits_phi","ECalHits_phi",50,-3.5,3.5);
  ECalHits_2DEtaPhi = new TH2F("ECalHits_2DEtaPhi","ECalHits_2DEtaPhi",50,-4.,4.,50,-3.5,3.5);
  ECalHits_2DEtaPhi->GetXaxis()->SetTitle("#eta");
  ECalHits_2DEtaPhi->GetYaxis()->SetTitle("#phi");
  ECalHits1000_2DEtaPhi = new TH2F("ECalHits > 1000GeV #eta #phi","ECalHits > 1000GeV #eta #phi",50,-4.,4.,50,-3.5,3.5);
  ECalHits1000_2DEtaPhi->GetXaxis()->SetTitle("#eta");
  ECalHits1000_2DEtaPhi->GetYaxis()->SetTitle("#phi");
  ECalHits_2DXY = new TH2F("ECalHits_2DXY","ECalHits_2DXY",50,-200.,200.,50,-200.,200.);
  ECalHits_2DXY->GetXaxis()->SetTitle("[cm]");
  ECalHits_2DXY->GetYaxis()->SetTitle("[cm]");
  EBHits_RZ = new TH2F("EBHits RZ","EBHits RZ",100,-300.,300.,20,0.,200.);
  EBHits_RZ->GetYaxis()->SetTitle("R [cm]");
  EBHits_RZ->GetXaxis()->SetTitle("Z [cm]");
  ECalHits1000_2DXY = new TH2F("ECalHits > 1000GeV XY","ECalHits > 1000GeV XY",50,-200.,200.,50,-200.,200.);
  ECalHits1000_2DXY->GetXaxis()->SetTitle("[cm]");
  ECalHits1000_2DXY->GetYaxis()->SetTitle("[cm]");
  EBHits1000_RZ = new TH2F("EBHits > 1000GeV RZ","EBHits > 1000GeV RZ",100,-300.,300.,20,0.,200.);
  EBHits1000_RZ->GetYaxis()->SetTitle("R [cm]");
  EBHits1000_RZ->GetXaxis()->SetTitle("Z [cm]");
  ECalHits_x = new TH1F("ECalHits_x","ECalHits_x",50,-200.,200.);
  ECalHits_x->GetXaxis()->SetTitle("[cm]");
  ECalHits_y = new TH1F("ECalHits_y","ECalHits_y",50,-200.,200.);
  ECalHits_y->GetXaxis()->SetTitle("[cm]");
  ECalHits_z = new TH1F("ECalHits_z","ECalHits_z",50,-200.,200.);
  ECalHits_z->GetXaxis()->SetTitle("[cm]");

  HCalHits_energy = new TH1F("HCalHits_energy","HCalHits_energy",100,0.,8000.);
  HCalHits_energy->GetXaxis()->SetTitle("[GeV]");
  HCalHits_eta = new TH1F("HCalHits_eta","HCalHits_eta",50,-4.,4.);
  HCalHits_phi = new TH1F("HCalHits_phi","HCalHits_phi",50,-3.5,3.5);
  HCalHits_2DEtaPhi = new TH2F("HCalHits_2DEtaPhi","HCalHits_2DEtaPhi",50,-4.,4.,50,-3.5,3.5);
  HCalHits_2DEtaPhi->GetXaxis()->SetTitle("#eta");
  HCalHits_2DEtaPhi->GetYaxis()->SetTitle("#phi");
  HCalHits_x = new TH1F("HCalHits_x","HCalHits_x",50,-1000.,1000.);
  HCalHits_x->GetXaxis()->SetTitle("[cm]");
  HCalHits_y = new TH1F("HCalHits_y","HCalHits_y",50,-1000.,1000.);
  HCalHits_y->GetXaxis()->SetTitle("[cm]");
  HCalHits_z = new TH1F("HCalHits_z","HCalHits_z",50,-1000.,1000.);
  HCalHits_z->GetXaxis()->SetTitle("[cm]");

  RHadron1_pdgId = new TH1F("RHadron1_pdgId", "RHadron1_pdgId", 100, 1000599., 1100000.);
  RHadron1_mass = new TH1F("RHadron1_mass", "RHadron1_mass", 100, 0., 3000.);
  RHadron1_mass->GetXaxis()->SetTitle("[GeV]");
  RHadron1_px = new TH1F("RHadron1_px","RHadron1_px",100,-10000.,10000.);
  RHadron1_px->GetXaxis()->SetTitle("[GeV]");
  RHadron1_py = new TH1F("RHadron1_py","RHadron1_py",100,-10000.,10000.);
  RHadron1_py->GetXaxis()->SetTitle("[GeV]");
  RHadron1_pz = new TH1F("RHadron1_pz","RHadron1_pz",100,-10000.,10000.);
  RHadron1_pz->GetXaxis()->SetTitle("[GeV]");
  RHadron1_pxpy = new TH2F("RHadron1_pxpy","RHadron1_pxpy",100,-3000.,3000.,100,-3000.,3000.);
  RHadron1_pxpy->GetXaxis()->SetTitle("[GeV]");
  RHadron1_pxpy->GetYaxis()->SetTitle("[GeV]");

  RHadron2_pdgId = new TH1F("RHadron2_pdgId", "RHadron2_pdgId", 100, 1000599., 1100000.);
  RHadron2_mass = new TH1F("RHadron2_mass", "RHadron2_mass", 100, 0., 3000.);
  RHadron2_mass->GetXaxis()->SetTitle("[GeV]");
  RHadron2_px = new TH1F("RHadron2_px","RHadron2_px",100,-10000.,10000.);
  RHadron2_px->GetXaxis()->SetTitle("[GeV]");
  RHadron2_py = new TH1F("RHadron2_py","RHadron2_py",100,-10000.,10000.);
  RHadron2_py->GetXaxis()->SetTitle("[GeV]");
  RHadron2_pz = new TH1F("RHadron2_pz","RHadron2_pz",100,-10000.,10000.);
  RHadron2_pz->GetXaxis()->SetTitle("[GeV]");
  RHadron2_pxpy = new TH2F("RHadron2_pxpy","RHadron2_pxpy",100,-3000.,3000.,100,-3000.,3000.);
  RHadron2_pxpy->GetXaxis()->SetTitle("[GeV]");
  RHadron2_pxpy->GetYaxis()->SetTitle("[GeV]");

  RHadron_ECalSubDetector = new TH1F("RHadron_ECalSubDetector","RHadron_ECalSubDetector",2,0.,2.);
  RHadron_ECalSubDetector->GetXaxis()->SetBinLabel(1,"Barrel");
  RHadron_ECalSubDetector->GetXaxis()->SetBinLabel(2,"Endcap");

  SimVertex_2DXY = new TH2F("SimVertex_2DXY","SimVertex_2DXY",50,-200.,200.,50,-200.,200.);
  SimVertex_2DXY->GetXaxis()->SetTitle("[cm]");
  SimVertex_2DXY->GetYaxis()->SetTitle("[cm]");
  SimVertex_2DRZ = new TH2F("SimVertex_2DRZ","SimVertex_2DRZ",100,-300.,300.,100,0.,1000.);
  SimVertex_2DRZ->GetXaxis()->SetTitle("[cm]");
  SimVertex_2DRZ->GetYaxis()->SetTitle("[cm]");
  
  evtcount = 0;
  evtcount00 = 0;
  evtcount01 = 0;
  evtcount01n = 0;
  evtcount01c = 0;
  evtcount11 = 0;
  badiqstat = 0;
  passevtcount0 = 0;
  passevtcount1 = 0;
  passevtcount00 = 0;
  passevtcount01 = 0;
  passevtcount01n = 0;
  passevtcount01c = 0;
  passevtcount11 = 0;

  ngencheta09 = 0;

}


//destructor
SpikedRHadronAnalyzer::~SpikedRHadronAnalyzer() {
  // Write histograms and file
  outputFile_->cd();

  ECalHits_energy->Write();
  ECalHits_eta->Write();
  ECalHits_phi->Write();
  ECalHits_2DEtaPhi->Write();
  ECalHits1000_2DEtaPhi->Write();
  ECalHits_2DXY->Write(); 
  ECalHits1000_2DXY->Write();
  EBHits_RZ->Write();
  EBHits1000_RZ->Write();
  ECalHits_x->Write();
  ECalHits_y->Write();
  ECalHits_z->Write();

  HCalHits_energy->Write();
  HCalHits_eta->Write();
  HCalHits_phi->Write();
  HCalHits_2DEtaPhi->Write();
  HCalHits_x->Write();
  HCalHits_y->Write();
  HCalHits_z->Write();

  RHadron1_pdgId->Write();
  RHadron1_mass->Write();
  RHadron1_px->Write();
  RHadron1_py->Write();
  RHadron1_pz->Write();
  RHadron1_pxpy->Write();

  RHadron2_pdgId->Write();
  RHadron2_mass->Write();
  RHadron2_px->Write();
  RHadron2_py->Write();
  RHadron2_pz->Write();
  RHadron2_pxpy->Write();

  RHadron_ECalSubDetector->Write();

  SimVertex_2DXY->Write();
  SimVertex_2DRZ->Write();

  std::cout << "saving MET trigger histograms..." << std::endl;
  outputFile_->Write();
  outputFile_->Close();
  csv.close();
}

void SpikedRHadronAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

   evtcount++;

   // Generator particles
   edm::Handle<vector<reco::GenParticle>> genColl;
   iEvent.getByToken(genParticlesToken_, genColl);
//   std::cout << " Gen particles: " << std::endl;

   vector<simhitinfo> linkedhits1;
   vector<simhitinfo> linkedhits2;

// simhit code inspired by
//   Validation/TrackerHits/src/TrackerHitAnalyzer.cc

   // PSimHits
   /////////////////////////////////
   // get Silicon TIB information
   //////////////////////////////////
   // extract TIB low container
   edm::Handle<edm::PSimHitContainer> SiTIBLowContainer;
   iEvent.getByToken(edmPSimHitContainer_siTIBLow_Token_, SiTIBLowContainer);
   if (!SiTIBLowContainer.isValid()) {
     edm::LogError("TrackerHitProducer::analyze") << "Unable to find TrackerHitsTIBLowTof in event!";
     return;
   }
   //////////////////////////////////
   // extract TIB high container
   edm::Handle<edm::PSimHitContainer> SiTIBHighContainer;
   iEvent.getByToken(edmPSimHitContainer_siTIBHigh_Token_, SiTIBHighContainer);
   if (!SiTIBHighContainer.isValid()) {
     edm::LogError("TrackerHitProducer::analyze") << "Unable to find TrackerHitsTIBHighTof in event!";
     return;
   }

   ///////////////////////////////////
   // get Silicon TOB information
   //////////////////////////////////
   // extract TOB low container
   edm::Handle<edm::PSimHitContainer> SiTOBLowContainer;
   iEvent.getByToken(edmPSimHitContainer_siTOBLow_Token_, SiTOBLowContainer);
   if (!SiTOBLowContainer.isValid()) {
     edm::LogError("TrackerHitProducer::analyze") << "Unable to find TrackerHitsTOBLowTof in event!";
     return;
   }
   //////////////////////////////////
   // extract TOB high container
   edm::Handle<edm::PSimHitContainer> SiTOBHighContainer;
   iEvent.getByToken(edmPSimHitContainer_siTOBHigh_Token_, SiTOBHighContainer);
   if (!SiTOBHighContainer.isValid()) {
     edm::LogError("TrackerHitProducer::analyze") << "Unable to find TrackerHitsTOBHighTof in event!";
     return;
   }

   //
   /////////////////////////////////////
   // get Silicon TID information
   //////////////////////////////////
   // extract TID low container
   edm::Handle<edm::PSimHitContainer> SiTIDLowContainer;
   iEvent.getByToken(edmPSimHitContainer_siTIDLow_Token_, SiTIDLowContainer);
   if (!SiTIDLowContainer.isValid()) {
     edm::LogError("TrackerHitProducer::analyze") << "Unable to find TrackerHitsTIDLowTof in event!";
     return;
   }
   //////////////////////////////////
   // extract TID high container
   edm::Handle<edm::PSimHitContainer> SiTIDHighContainer;
   iEvent.getByToken(edmPSimHitContainer_siTIDHigh_Token_, SiTIDHighContainer);
   if (!SiTIDHighContainer.isValid()) {
     edm::LogError("TrackerHitProducer::analyze") << "Unable to find TrackerHitsTIDHighTof in event!";
     return;
   }

   ///////////////////////////////////
   // get Silicon TEC information
   //////////////////////////////////
   // extract TEC low container
   edm::Handle<edm::PSimHitContainer> SiTECLowContainer;
   iEvent.getByToken(edmPSimHitContainer_siTECLow_Token_, SiTECLowContainer);
   if (!SiTECLowContainer.isValid()) {
     edm::LogError("TrackerHitProducer::analyze") << "Unable to find TrackerHitsTECLowTof in event!";
     return;
   }
   //////////////////////////////////
   // extract TEC high container
   edm::Handle<edm::PSimHitContainer> SiTECHighContainer;
   iEvent.getByToken(edmPSimHitContainer_siTECHigh_Token_, SiTECHighContainer);
   if (!SiTECHighContainer.isValid()) {
     edm::LogError("TrackerHitProducer::analyze") << "Unable to find TrackerHitsTECHighTof in event!";
     return;
   }

   /////////////////////////////////
   // get Pixel Barrel information
   ////////////////////////////////
   // extract low container
   edm::Handle<edm::PSimHitContainer> PxlBrlLowContainer;
   iEvent.getByToken(edmPSimHitContainer_pxlBrlLow_Token_, PxlBrlLowContainer);
   if (!PxlBrlLowContainer.isValid()) {
     edm::LogError("TrackerHitAnalyzer::analyze") << "Unable to find TrackerHitsPixelBarrelLowTof in event!";
     return;
   }
   // extract high container
   edm::Handle<edm::PSimHitContainer> PxlBrlHighContainer;
   iEvent.getByToken(edmPSimHitContainer_pxlBrlHigh_Token_, PxlBrlHighContainer);
   if (!PxlBrlHighContainer.isValid()) {
     edm::LogError("TrackerHitAnalyzer::analyze") << "Unable to find TrackerHitsPixelBarrelHighTof in event!";
     return;
   }
   /////////////////////////////////
   // get Pixel Forward information
   ////////////////////////////////
   // extract low container
   edm::Handle<edm::PSimHitContainer> PxlFwdLowContainer;
   iEvent.getByToken(edmPSimHitContainer_pxlFwdLow_Token_, PxlFwdLowContainer);
   if (!PxlFwdLowContainer.isValid()) {
     edm::LogError("TrackerHitAnalyzer::analyze") << "Unable to find TrackerHitsPixelEndcapLowTof in event!";
     return;
   }
   // extract high container
   edm::Handle<edm::PSimHitContainer> PxlFwdHighContainer;
   iEvent.getByToken(edmPSimHitContainer_pxlFwdHigh_Token_, PxlFwdHighContainer);
   if (!PxlFwdHighContainer.isValid()) {
     edm::LogError("TrackerHitAnalyzer::analyze") << "Unable to find TrackerHitsPixelEndcapHighTof in event!";
     return;
   }
    ////////////////////////////////
    // Get calorimiter information//
    ////////////////////////////////
    //EcalEB
    edm::Handle<edm::PCaloHitContainer> EcalEBContainer;
    iEvent.getByToken(edmCaloHitContainer_EcalHitsEB_Token_, EcalEBContainer);
    if (!EcalEBContainer.isValid()) {
        edm::LogError("CaloHitProducer::analyze") << "Unable to find EcalHitsEB in event!";
        return;
   }
    //EcalEE
    edm::Handle<edm::PCaloHitContainer> EcalEEContainer;
    iEvent.getByToken(edmCaloHitContainer_EcalHitsEE_Token_, EcalEEContainer);
    if (!EcalEEContainer.isValid()) {
        edm::LogError("CaloHitProducer::analyze") << "Unable to find EcalHitsEE in event!";
        return;
   }
    //EcalES
    edm::Handle<edm::PCaloHitContainer> EcalESContainer;
    iEvent.getByToken(edmCaloHitContainer_EcalHitsES_Token_, EcalESContainer);
    if (!EcalESContainer.isValid()) {
        edm::LogError("CaloHitProducer::analyze") << "Unable to find EcalHitsES in event!";
        return;
   }
    //Hcal
    edm::Handle<edm::PCaloHitContainer> HcalContainer;
    iEvent.getByToken(edmCaloHitContainer_HcalHits_Token_, HcalContainer);
    if (!HcalContainer.isValid()) {
        edm::LogError("CaloHitProducer::analyze") << "Unable to find HcalHits in event!";
        return;
   }

   /////////////////////////////////
   // get muon information
   ////////////////////////////////
   // extract DT container
   edm::Handle<edm::PSimHitContainer> MuonDTContainer;
   iEvent.getByToken(edmPSimHitContainer_muonDT_Token_, MuonDTContainer);
   if (!MuonDTContainer.isValid()) {
     edm::LogError("TrackerHitAnalyzer::analyze") << "Unable to find MuonDT simhits in event!";
     return;
   }

   // extract CSC container
   edm::Handle<edm::PSimHitContainer> MuonCSCContainer;
   iEvent.getByToken(edmPSimHitContainer_muonCSC_Token_, MuonCSCContainer);
   if (!MuonCSCContainer.isValid()) {
     edm::LogError("TrackerHitAnalyzer::analyze") << "Unable to find MuonCSC simhits in event!";
     return;
   }

    // extract RPC container
    edm::Handle<edm::PSimHitContainer> MuonRPCContainer;
    iEvent.getByToken(edmPSimHitContainer_muonRPC_Token_, MuonRPCContainer);
    if (!MuonRPCContainer.isValid()) {
        edm::LogError("TrackerHitAnalyzer::analyze") << "Unable to find MuonRPC simhits in event!";
        return;
    }

    // extract GEM container
    edm::Handle<edm::PSimHitContainer> MuonGEMContainer;
    iEvent.getByToken(edmPSimHitContainer_muonGEM_Token_, MuonGEMContainer);
    if (!MuonGEMContainer.isValid()) {
        edm::LogError("TrackerHitAnalyzer::analyze") << "Unable to find MuonGEM simhits in event!";
        return;
    }

    // Get G4SimTracks
    edm::Handle<edm::SimTrackContainer> G4TrkContainer;
    iEvent.getByToken(edmSimTrackContainerToken_, G4TrkContainer);
    if (!G4TrkContainer.isValid()) {
      edm::LogError("TrackerHitAnalyzer::analyze") << "Unable to find SimTrack in event!";
      return;
   }

    // Get G4SimVertices
    edm::Handle<edm::SimVertexContainer> G4VtxContainer;
    iEvent.getByToken(edmSimVertexContainerToken_, G4VtxContainer);
    //iEvent.getByLabel("g4SimHits",G4VtxContainer);
    if (!G4VtxContainer.isValid()) {
      edm::LogError("TrackerHitAnalyzer::analyze") << "Unable to find SimVertex in event!";
      return;
    }

   // Import tracker, calorimiter, and muon geometry
  edm::ESHandle<TrackerGeometry> tkGeometry;
  iSetup.get<TrackerDigiGeometryRecord>().get(tkGeometry);

  edm::ESHandle<CaloGeometry> caloGeometry;
  iSetup.get<CaloGeometryRecord>().get(caloGeometry); 

  edm::ESHandle<DTGeometry> DTGeometry;
  iSetup.get<MuonGeometryRecord>().get(DTGeometry);

  edm::ESHandle<CSCGeometry> CSCGeometry;
  iSetup.get<MuonGeometryRecord>().get(CSCGeometry);

  edm::ESHandle<RPCGeometry> RPCGeometry;
  iSetup.get<MuonGeometryRecord>().get(RPCGeometry);

  edm::ESHandle<GEMGeometry> GEMGeometry;
  iSetup.get<MuonGeometryRecord>().get(GEMGeometry);

   // Find R-hadrons
   const reco::GenParticle *genrhad1=0;
   const reco::GenParticle *genrhad2=0;

   for (unsigned int ig=0;ig<genColl->size();ig++) {
     const reco::GenParticle & part = (*genColl)[ig];
     if ((abs(part.pdgId())>1000600)&&(abs(part.pdgId())<1100000)&&(part.status()==1)) {
   //cout << "SUSY part = " << part.pdgId() << endl;
       if (genrhad1==0) {
         genrhad1 = &(*genColl)[ig];
       } else if (genrhad2==0) {
         genrhad2 = &(*genColl)[ig];
       } else {
         std::cout << "WARNING - found more than two R-hadrons" << std::endl;
       }
         // find and increment count
       bool found = false;

       for (std::vector<std::pair<int, int>>::iterator it = genrhadcounts.begin(); (it != genrhadcounts.end())&&(!found); ++it) {
         if (part.pdgId()==(*it).first) {
           (*it).second++;
           found = true;
         }  
       }
       if (!found) {
         std::pair<int,int> pair1(part.pdgId(),1);
         genrhadcounts.push_back(pair1);
         //cout << " NEW ID = " << part.pdgId() << " " << endl;
       }
     }
   }

  //Plot Rhadron momenta, mass, and pdgId
  RHadron1_pdgId->Fill(genrhad1->pdgId());
  RHadron1_mass->Fill(genrhad1->mass());
  RHadron1_px->Fill(genrhad1->p4().px());
  RHadron1_py->Fill(genrhad1->p4().py());
  RHadron1_pz->Fill(genrhad1->p4().pz());
  RHadron1_pxpy->Fill(genrhad1->p4().px(),genrhad1->p4().py());

  RHadron2_pdgId->Fill(genrhad2->pdgId());
  RHadron2_mass->Fill(genrhad2->mass());
  RHadron2_px->Fill(genrhad2->p4().px());
  RHadron2_py->Fill(genrhad2->p4().py());
  RHadron2_pz->Fill(genrhad2->p4().pz());
  RHadron2_pxpy->Fill(genrhad2->p4().px(),genrhad2->p4().py());

  //Plot which subdetector the R-hadrons hit in the ECAL based on eta
  if (abs(genrhad1->eta())<1.479) {
    RHadron_ECalSubDetector->Fill(0);
  } else if ((abs(genrhad1->eta())>=1.479)&&(abs(genrhad1->eta())<3.0)) {
    RHadron_ECalSubDetector->Fill(1.5);
  }

  if (abs(genrhad2->eta())<1.479) {
    RHadron_ECalSubDetector->Fill(0);
  } else if ((abs(genrhad2->eta())>=1.479)&&(abs(genrhad2->eta())<3.0)) {
    RHadron_ECalSubDetector->Fill(1.5);
  }


  //Check the tracker for R-hadron hits
  edm::PSimHitContainer::const_iterator itHit;

  int nhitstib[20] = {0,0,0,0,0,0,0,0,0,0};
  int nhitstib1[20] = {0,0,0,0,0,0,0,0,0,0};
  int nhitstib2[20] = {0,0,0,0,0,0,0,0,0,0};

  //TIB
  for (itHit = SiTIBLowContainer->begin(); itHit != SiTIBLowContainer->end(); ++itHit) {
    DetId detid = DetId(itHit->detUnitId());
    const GeomDetUnit *det = (const GeomDetUnit *)tkGeometry->idToDetUnit(detid);
    GlobalPoint gpos = det->toGlobal(itHit->localPosition());
    float r = sqrt(gpos.x()*gpos.x()+gpos.y()*gpos.y());

    if ((abs(itHit->particleType())>1000600)&&(abs(itHit->particleType())<1100000)) {
      float dphi1 = genrhad1->phi() - gpos.phi();
      float dphi2 = genrhad2->phi() - gpos.phi();
      float deta1 = genrhad1->eta() - gpos.eta();
      float deta2 = genrhad2->eta() - gpos.eta();
      float dR1 = sqrt(dphi1*dphi1+deta1*deta1);
      float dR2 = sqrt(dphi2*dphi2+deta2*deta2);
      
      //const reco::GenParticle *genrhadclosest=0;
      int iclosest = -1;
      float dRclosest = 999999.;
      if (dR1<dR2) {
        dRclosest = dR1;
        iclosest = 1;
        //genrhadclosest = genrhad1;
      } else {
        dRclosest = dR2;
        iclosest = 2;
        //genrhadclosest = genrhad2;
      }
      
      int ilayer = -1;
      if (dRclosest<1.0) {
        if ((r>23)&&(r<29)) { // first layer
          ilayer = TIB0 + 1;
        } else if ((r>31)&&(r<36)) { // second layer
          ilayer = TIB0 + 2;
        } else if ((r>39)&&(r<45)) { // third layer
          ilayer = TIB0 + 3;
        } else if ((r>47)&&(r<53)) { // fourth layer
          ilayer = TIB0 + 4;
        }
      }
      if (ilayer>-1) {
        nhitstib[ilayer]++;
        if (dR1<0.4) nhitstib1[ilayer]++;
        if (dR2<0.4) nhitstib2[ilayer]++;
        if (iclosest==1) {
          linkedhits1.push_back({ilayer,itHit->particleType()});
        } else if (iclosest==2) {
          linkedhits2.push_back({ilayer,itHit->particleType()});
        }
      }
    }
  }

  //////////////////////////////////
  // Grab sim hits in the tracker //
  //////////////////////////////////
/*
  edm::PSimHitContainer::const_iterator simHit;
  for (simHit = SiTIBLowContainer->begin(); simHit != SiTIBLowContainer->end(); ++simHit) {
    DetId detid = DetId(simHit->detUnitId());
    const GeomDetUnit *det = (const GeomDetUnit *)tkGeometry->idToDetUnit(detid);
    GlobalPoint gpos = det->toGlobal(simHit->localPosition());
    csv << evtcount << "," << genrhad1->p4().px() << "," << genrhad1->p4().py() << "," << genrhad2->p4().px() << "," << genrhad2->p4().py() << "," << simHit->particleType() << "," << "TrackerTIBLow" << "," << gpos.x() << "," << gpos.y() << "," << gpos.z() << "," << sqrt(gpos.x()*gpos.x()+gpos.y()*gpos.y()) << "," << simHit->energyLoss() << "\n";
  }

  for (simHit = SiTIBHighContainer->begin(); simHit != SiTIBHighContainer->end(); ++simHit) {
    DetId detid = DetId(simHit->detUnitId());
    const GeomDetUnit *det = (const GeomDetUnit *)tkGeometry->idToDetUnit(detid);
    GlobalPoint gpos = det->toGlobal(simHit->localPosition());
    csv << evtcount << "," << genrhad1->p4().px() << "," << genrhad1->p4().py() << "," << genrhad2->p4().px() << "," << genrhad2->p4().py() << "," << simHit->particleType() << "," << "TrackerTIBHigh" << "," << gpos.x() << "," << gpos.y() << "," << gpos.z() << "," << sqrt(gpos.x()*gpos.x()+gpos.y()*gpos.y()) << "," << simHit->energyLoss() << "\n";
  }

  for (simHit = SiTOBLowContainer->begin(); simHit != SiTOBLowContainer->end(); ++simHit) {
    DetId detid = DetId(simHit->detUnitId());
    const GeomDetUnit *det = (const GeomDetUnit *)tkGeometry->idToDetUnit(detid);
    GlobalPoint gpos = det->toGlobal(simHit->localPosition());
    csv << evtcount << "," << genrhad1->p4().px() << "," << genrhad1->p4().py() << "," << genrhad2->p4().px() << "," << genrhad2->p4().py() << "," << simHit->particleType() << "," << "TrackerTOBLow" << "," << gpos.x() << "," << gpos.y() << "," << gpos.z() << "," << sqrt(gpos.x()*gpos.x()+gpos.y()*gpos.y()) << "," << simHit->energyLoss() << "\n";
  }

  for (simHit = SiTOBHighContainer->begin(); simHit != SiTOBHighContainer->end(); ++simHit) {
    DetId detid = DetId(simHit->detUnitId());
    const GeomDetUnit *det = (const GeomDetUnit *)tkGeometry->idToDetUnit(detid);
    GlobalPoint gpos = det->toGlobal(simHit->localPosition());
    csv << evtcount << "," << genrhad1->p4().px() << "," << genrhad1->p4().py() << "," << genrhad2->p4().px() << "," << genrhad2->p4().py() << "," << simHit->particleType() << "," << "TrackerTOBHigh" << "," << gpos.x() << "," << gpos.y() << "," << gpos.z() << "," << sqrt(gpos.x()*gpos.x()+gpos.y()*gpos.y()) << "," << simHit->energyLoss() << "\n";
  }

  for (simHit = SiTIDLowContainer->begin(); simHit != SiTIDLowContainer->end(); ++simHit) {
    DetId detid = DetId(simHit->detUnitId());
    const GeomDetUnit *det = (const GeomDetUnit *)tkGeometry->idToDetUnit(detid);
    GlobalPoint gpos = det->toGlobal(simHit->localPosition());
    csv << evtcount << "," << genrhad1->p4().px() << "," << genrhad1->p4().py() << "," << genrhad2->p4().px() << "," << genrhad2->p4().py() << "," << simHit->particleType() << "," << "TrackerTIDLow" << "," << gpos.x() << "," << gpos.y() << "," << gpos.z() << "," << sqrt(gpos.x()*gpos.x()+gpos.y()*gpos.y()) << "," << simHit->energyLoss() << "\n";
  }

  for (simHit = SiTIDHighContainer->begin(); simHit != SiTIDHighContainer->end(); ++simHit) {
    DetId detid = DetId(simHit->detUnitId());
    const GeomDetUnit *det = (const GeomDetUnit *)tkGeometry->idToDetUnit(detid);
    GlobalPoint gpos = det->toGlobal(simHit->localPosition());
    csv << evtcount << "," << genrhad1->p4().px() << "," << genrhad1->p4().py() << "," << genrhad2->p4().px() << "," << genrhad2->p4().py() << "," << simHit->particleType() << "," << "TrackerTIDHigh" << "," << gpos.x() << "," << gpos.y() << "," << gpos.z() << "," << sqrt(gpos.x()*gpos.x()+gpos.y()*gpos.y()) << "," << simHit->energyLoss() << "\n";
  }

  for (simHit = SiTECLowContainer->begin(); simHit != SiTECLowContainer->end(); ++simHit) {
    DetId detid = DetId(simHit->detUnitId());
    const GeomDetUnit *det = (const GeomDetUnit *)tkGeometry->idToDetUnit(detid);
    GlobalPoint gpos = det->toGlobal(simHit->localPosition());
    csv << evtcount << "," << genrhad1->p4().px() << "," << genrhad1->p4().py() << "," << genrhad2->p4().px() << "," << genrhad2->p4().py() << "," << simHit->particleType() << "," << "TrackerTECLow" << "," << gpos.x() << "," << gpos.y() << "," << gpos.z() << "," << sqrt(gpos.x()*gpos.x()+gpos.y()*gpos.y()) << "," << simHit->energyLoss() << "\n";
  }

  for (simHit = SiTECHighContainer->begin(); simHit != SiTECHighContainer->end(); ++simHit) {
    DetId detid = DetId(simHit->detUnitId());
    const GeomDetUnit *det = (const GeomDetUnit *)tkGeometry->idToDetUnit(detid);
    GlobalPoint gpos = det->toGlobal(simHit->localPosition());
    csv << evtcount << "," << genrhad1->p4().px() << "," << genrhad1->p4().py() << "," << genrhad2->p4().px() << "," << genrhad2->p4().py() << "," << simHit->particleType() << "," << "TrackerTECHigh" << "," << gpos.x() << "," << gpos.y() << "," << gpos.z() << "," << sqrt(gpos.x()*gpos.x()+gpos.y()*gpos.y()) << "," << simHit->energyLoss() << "\n";
  }

  for (simHit = PxlBrlLowContainer->begin(); simHit != PxlBrlLowContainer->end(); ++simHit) {
    DetId detid = DetId(simHit->detUnitId());
    const GeomDetUnit *det = (const GeomDetUnit *)tkGeometry->idToDetUnit(detid);
    GlobalPoint gpos = det->toGlobal(simHit->localPosition());
    csv << evtcount << "," << genrhad1->p4().px() << "," << genrhad1->p4().py() << "," << genrhad2->p4().px() << "," << genrhad2->p4().py() << "," << simHit->particleType() << "," << "TrackerPxlBrlLow" << "," << gpos.x() << "," << gpos.y() << "," << gpos.z() << "," << sqrt(gpos.x()*gpos.x()+gpos.y()*gpos.y()) << "," << simHit->energyLoss() << "\n";
  }

  for (simHit = PxlBrlHighContainer->begin(); simHit != PxlBrlHighContainer->end(); ++simHit) {
    DetId detid = DetId(simHit->detUnitId());
    const GeomDetUnit *det = (const GeomDetUnit *)tkGeometry->idToDetUnit(detid);
    GlobalPoint gpos = det->toGlobal(simHit->localPosition());
    csv << evtcount << "," << genrhad1->p4().px() << "," << genrhad1->p4().py() << "," << genrhad2->p4().px() << "," << genrhad2->p4().py() << "," << simHit->particleType() << "," << "TrackerPxlBrlHigh" << "," << gpos.x() << "," << gpos.y() << "," << gpos.z() << "," << sqrt(gpos.x()*gpos.x()+gpos.y()*gpos.y()) << "," << simHit->energyLoss() << "\n";
  }

  for (simHit = PxlFwdLowContainer->begin(); simHit != PxlFwdLowContainer->end(); ++simHit) {
    DetId detid = DetId(simHit->detUnitId());
    const GeomDetUnit *det = (const GeomDetUnit *)tkGeometry->idToDetUnit(detid);
    GlobalPoint gpos = det->toGlobal(simHit->localPosition());
    csv << evtcount << "," << genrhad1->p4().px() << "," << genrhad1->p4().py() << "," << genrhad2->p4().px() << "," << genrhad2->p4().py() << "," << simHit->particleType() << "," << "TrackerPxlFwdLow" << "," << gpos.x() << "," << gpos.y() << "," << gpos.z() << "," << sqrt(gpos.x()*gpos.x()+gpos.y()*gpos.y()) << "," << simHit->energyLoss() << "\n";
  }

  for (simHit = PxlFwdHighContainer->begin(); simHit != PxlFwdHighContainer->end(); ++simHit) {
    DetId detid = DetId(simHit->detUnitId());
    const GeomDetUnit *det = (const GeomDetUnit *)tkGeometry->idToDetUnit(detid);
    GlobalPoint gpos = det->toGlobal(simHit->localPosition());
    csv << evtcount << "," << genrhad1->p4().px() << "," << genrhad1->p4().py() << "," << genrhad2->p4().px() << "," << genrhad2->p4().py() << "," << simHit->particleType() << "," << "TrackerPxlFwdHigh" << "," << gpos.x() << "," << gpos.y() << "," << gpos.z() << "," << sqrt(gpos.x()*gpos.x()+gpos.y()*gpos.y()) << "," << simHit->energyLoss() << "\n";
  }
  */
  ///////////////////////////////////////
  // Grab sim hits in the calorimiters //
  ///////////////////////////////////////

  edm::PCaloHitContainer::const_iterator caloHit;
  edm::SimTrackContainer::const_iterator simTrack;
  edm::SimVertexContainer::const_iterator simVertex;

  // EB
  for (caloHit = EcalEBContainer->begin(); caloHit != EcalEBContainer->end(); ++caloHit) {

    unsigned trackId = caloHit->geantTrackId();
    std::vector<int> myVector;
    for (simVertex = G4VtxContainer->begin(); simVertex != G4VtxContainer->end(); ++simVertex) {
      unsigned parentId = simVertex->parentIndex();

      if (parentId == trackId) {
        int processType = simVertex->processType();
        unsigned vertexId = simVertex->vertexId();

        if (caloHit->energy() >= 1000.) {
          int parentToVertex = 0;
          for (simTrack = G4TrkContainer->begin(); simTrack != G4TrkContainer->end(); ++simTrack) {
            unsigned vertid = simTrack->vertIndex();

            if (simTrack->trackId() == trackId) {
              parentToVertex = simTrack->type();
            }

            if (vertid == vertexId){
              myVector.push_back(simTrack->type());
            }
          }
          std::string vectorAsString = vectorToString(myVector);
          std::cout << "Sim Vertex ID = " << vertexId << ", Calohit Energy = " << caloHit->energy() << ", process type = " << simVertex->processType() << ", parent particle = " << parentToVertex << ", daughters = " << vectorAsString << std::endl;
        }
      }
    }
  }
/*
    if (processType == 0) continue;

    std::vector<int> myVector;
    for (simTrack = G4TrkContainer->begin(); simTrack != G4TrkContainer->end(); ++simTrack) {
      unsigned vertid = simTrack->vertIndex();
      if (simTrack->trackId() == trackId) {
        parentToVertex = simTrack->type();
      }
      if (vertid == vertexId){
        myVector.push_back(simTrack->type());
      }
    }
    std::string vectorAsString = vectorToString(myVector);
    std::cout << "Calohit Energy = " << caloHit->energy() << ", process type = " << processType << "parent particle = " << parentToVertex << "daughters = " << vectorAsString << std::endl;
  }
    /*
    // Declare variables
    unsigned trackId = caloHit->geantTrackId();
    unsigned vertexId = 0;
    int particleThatCausedHit = 0;
    int parentToVertex = 0;
    int processType = 0;

    // Grab the track that caused the calo hit
    for (simTrack = G4TrkContainer->begin(); simTrack != G4TrkContainer->end(); ++simTrack) {
      if (simTrack->trackId() == trackId) {
        vertexId = simTrack->vertIndex();
        particleThatCausedHit = simTrack->type();
        break;
      }
    }

    // If the track that caused the calo hit could not be found, continue to the next calo hit
    if (particleThatCausedHit == 0) continue;

    // Grab the vertex that the track came from, get the process type and the parents of the vertex
    for (simVertex = G4VtxContainer->begin(); simVertex != G4VtxContainer->end(); ++simVertex) {
      if (simVertex->vertexId() == vertexId) {
        processType = simVertex->processType();
        unsigned parentIndex = simVertex->parentIndex();
        for (simTrack = G4TrkContainer->begin(); simTrack != G4TrkContainer->end(); ++simTrack) {
          if (simTrack->trackId() == parentIndex) {
            parentToVertex = simTrack->type();
            break;
          }
        }
        break;
      }
    }

    // Grab the other daughters of the vertex, dump all of the information to a csv
    std::vector<int> myVector;
    for (simTrack = G4TrkContainer->begin(); simTrack != G4TrkContainer->end(); ++simTrack) {
      unsigned vertid = simTrack->vertIndex();
      if (vertid == vertexId){
        myVector.push_back(simTrack->type());
      }
    }
    std::string vectorAsString = vectorToString(myVector);
    csv << evtcount << "," << caloHitCounter << "," << caloHit->energy() << "," << processType << "," << particleThatCausedHit << "," << vectorAsString << "," << parentToVertex << "\n";
    caloHitCounter += 1;
  }
/*
    // Fill calohit energy histogram
    ECalHits_energy->Fill(caloHit->energy());

    // Fill calohit location histograms
    ECalHits_eta->Fill(gpos.eta());
    ECalHits_phi->Fill(gpos.phi());
    ECalHits_2DEtaPhi->Fill(gpos.eta(),gpos.phi());
    ECalHits_x->Fill(gpos.x());
    ECalHits_y->Fill(gpos.y());
    ECalHits_z->Fill(gpos.z());
    ECalHits_2DXY->Fill(gpos.x(),gpos.y());

    //Fill histograms after energy cut
    if (caloHit->energy()>=1000.) {
      ECalHits1000_2DEtaPhi->Fill(gpos.eta(),gpos.phi());
      ECalHits1000_2DXY->Fill(gpos.x(),gpos.y());
      EBHits1000_RZ->Fill(gpos.z(),sqrt(gpos.x()*gpos.x()+gpos.y()*gpos.y()));
    }


    //Fill EB specific histograms
    EBHits_RZ->Fill(gpos.z(),sqrt(gpos.x()*gpos.x()+gpos.y()*gpos.y()));

    // Fill csv for energy spike R-hadron analysis
    csv << evtcount << "," << genrhad1->p4().px() << "," << genrhad1->p4().py() << "," << genrhad2->p4().px() << "," << genrhad2->p4().py() << "," << 0 << "," << "EB" << "," << gpos.x() << "," << gpos.y() << "," << gpos.z() << "," << sqrt(gpos.x()*gpos.x()+gpos.y()*gpos.y()) << "," << caloHit->energy() << "\n";
  }
  /*
  // Grab calorimiter hits in ES
  for (caloHit = EcalESContainer->begin(); caloHit != EcalESContainer->end(); ++caloHit) {
    //DetId detid = DetId(caloHit->id());
    //GlobalPoint gpos = caloGeometry->getPosition(detid);

    //Get the simulated vertex position
    GlobalPoint gpos = GlobalPoint();
    unsigned trackId = caloHit->geantTrackId();
    for (simVertex = G4VtxContainer->begin(); simVertex != G4VtxContainer->end(); ++simVertex) {
      unsigned parentId = simVertex->parentIndex();
      if (parentId == trackId) {
        gpos = GlobalPoint(simVertex->position().x(),simVertex->position().y(),simVertex->position().z());
        break;
      }
    }

    //END OF WARNING//

    // Fill calohit energy histogram
    ECalHits_energy->Fill(caloHit->energy());

    // Fill calohit location histograms
    ECalHits_eta->Fill(gpos.eta());
    ECalHits_phi->Fill(gpos.phi());
    ECalHits_2DEtaPhi->Fill(gpos.eta(),gpos.phi());
    ECalHits_x->Fill(gpos.x());
    ECalHits_y->Fill(gpos.y());
    ECalHits_z->Fill(gpos.z());
    ECalHits_2DXY->Fill(gpos.x(),gpos.y());

    //Fill histograms after energy cut
    if (caloHit->energy()>=1000.) {
      ECalHits1000_2DEtaPhi->Fill(gpos.eta(),gpos.phi());
      ECalHits1000_2DXY->Fill(gpos.x(),gpos.y());
    }

    // Fill csv for energy spike R-hadron analysis
    csv << evtcount << "," << genrhad1->p4().px() << "," << genrhad1->p4().py() << "," << genrhad2->p4().px() << "," << genrhad2->p4().py() << "," << 0 << "," << "ES" << "," << gpos.x() << "," << gpos.y() << "," << gpos.z() << "," << sqrt(gpos.x()*gpos.x()+gpos.y()*gpos.y()) << "," << caloHit->energy() << "\n";
  }
  
  // Grab calorimiter hits in EE
  for (caloHit = EcalEEContainer->begin(); caloHit != EcalEEContainer->end(); ++caloHit) {
    //DetId detid = DetId(caloHit->id());
    //GlobalPoint gpos = caloGeometry->getPosition(detid);

    //Get the simulated vertex position
    GlobalPoint gpos = GlobalPoint();
    unsigned trackId = caloHit->geantTrackId();
    for (simVertex = G4VtxContainer->begin(); simVertex != G4VtxContainer->end(); ++simVertex) {
      unsigned parentId = simVertex->parentIndex();
      if (parentId == trackId) {
        gpos = GlobalPoint(simVertex->position().x(),simVertex->position().y(),simVertex->position().z());
        break;
      }
    }

    // Fill calohit energy histogram
    ECalHits_energy->Fill(caloHit->energy());

    // Fill calohit location histograms
    ECalHits_eta->Fill(gpos.eta());
    ECalHits_phi->Fill(gpos.phi());
    ECalHits_2DEtaPhi->Fill(gpos.eta(),gpos.phi());
    ECalHits_x->Fill(gpos.x());
    ECalHits_y->Fill(gpos.y());
    ECalHits_z->Fill(gpos.z());
    ECalHits_2DXY->Fill(gpos.x(),gpos.y());

    //Fill histograms after energy cut
    if (caloHit->energy()>=1000.) {
      ECalHits1000_2DEtaPhi->Fill(gpos.eta(),gpos.phi());
      ECalHits1000_2DXY->Fill(gpos.x(),gpos.y());
    }

    // Fill csv for energy spike R-hadron analysis
    csv << evtcount << "," << genrhad1->p4().px() << "," << genrhad1->p4().py() << "," << genrhad2->p4().px() << "," << genrhad2->p4().py() << "," << 0 << "," << "EE" << "," << gpos.x() << "," << gpos.y() << "," << gpos.z() << "," << sqrt(gpos.x()*gpos.x()+gpos.y()*gpos.y()) << "," << caloHit->energy() << "\n";
  }

  /*
  // Grab calorimiter hits in HCal
  //test::test(const edm::ParameterSet& iConfig)  
  //: tok_HRNDC_(esConsumes<HcalDDDRecConstants, HcalRecNumberingRecord, edm::Transition::BeginRun>()),
  //  tok_geom_(esConsumes<CaloGeometry, CaloGeometryRecord, edm::Transition::BeginRun>()), 
  //  tok_caloTopology_(esConsumes<CaloTopology, CaloTopologyRecord, edm::Transition::BeginRun>()), 
  //  tok_hcalTopology_(esConsumes<HcalTopology, HcalRecNumberingRecord, edm::Transition::BeginRun>()), 
  //  tok_resp_(esConsumes<HcalRespCorrs, HcalRespCorrsRcd, edm::Transition::BeginRun>()) {
  //for (caloHit = HcalContainer->begin(); caloHit != HcalContainer->end(); caloHit++) {
  //  HcalDetId hid;
  //  hid = HcalHitRelabeller::relabel(caloHit->id(), hcons_);
  //  GlobalPoint gpos = caloGeometry->getPosition(hid);
  //  csv << evtcount << "," << genrhad1->p4().px() << "," << genrhad1->p4().py() << "," << genrhad2->p4().px() << "," << genrhad2->p4().py() << "," << 0 << "," << "HCAL" << "," << gpos.x() << "," << gpos.y() << "," << gpos.z() << "," << sqrt(gpos.x()*gpos.x()+gpos.y()*gpos.y()) << "," << caloHit->energy() << "\n";
  //}

  //////////////////////////////////////
  // Grab sim hits in the muon system //
  //////////////////////////////////////

  for (simHit = MuonDTContainer->begin(); simHit != MuonDTContainer->end(); ++simHit) {
    DetId detid = DetId(simHit->detUnitId());
    const GeomDetUnit *det = (const GeomDetUnit *)DTGeometry->idToDetUnit(detid);
    GlobalPoint gpos = det->toGlobal(simHit->localPosition());
    csv << evtcount << "," << genrhad1->p4().px() << "," << genrhad1->p4().py() << "," << genrhad2->p4().px() << "," << genrhad2->p4().py() << "," << simHit->particleType() << "," << "MuonDT" << "," << gpos.x() << "," << gpos.y() << "," << gpos.z() << "," << sqrt(gpos.x()*gpos.x()+gpos.y()*gpos.y()) << "," << simHit->energyLoss() << "\n";
  }

  for (simHit = MuonCSCContainer->begin(); simHit != MuonCSCContainer->end(); ++simHit) {
    DetId detid = DetId(simHit->detUnitId());
    const GeomDetUnit *det = (const GeomDetUnit *)CSCGeometry->idToDetUnit(detid);
    GlobalPoint gpos = det->toGlobal(simHit->localPosition());
    csv << evtcount << "," << genrhad1->p4().px() << "," << genrhad1->p4().py() << "," << genrhad2->p4().px() << "," << genrhad2->p4().py() << "," << simHit->particleType() << "," << "MuonCSC" << "," << gpos.x() << "," << gpos.y() << "," << gpos.z() << "," << sqrt(gpos.x()*gpos.x()+gpos.y()*gpos.y()) << "," << simHit->energyLoss() << "\n";
  }

  for (simHit = MuonRPCContainer->begin(); simHit != MuonRPCContainer->end(); ++simHit) {
    DetId detid = DetId(simHit->detUnitId());
    const GeomDetUnit *det = (const GeomDetUnit *)RPCGeometry->idToDetUnit(detid);
    GlobalPoint gpos = det->toGlobal(simHit->localPosition());
    csv << evtcount << "," << genrhad1->p4().px() << "," << genrhad1->p4().py() << "," << genrhad2->p4().px() << "," << genrhad2->p4().py() << "," << simHit->particleType() << "," << "MuonRPC" << "," << gpos.x() << "," << gpos.y() << "," << gpos.z() << "," << sqrt(gpos.x()*gpos.x()+gpos.y()*gpos.y()) << "," << simHit->energyLoss() << "\n";
  }

  for (simHit = MuonGEMContainer->begin(); simHit != MuonGEMContainer->end(); ++simHit) {
    DetId detid = DetId(simHit->detUnitId());
    const GeomDetUnit *det = (const GeomDetUnit *)GEMGeometry->idToDetUnit(detid);
    GlobalPoint gpos = det->toGlobal(simHit->localPosition());
    csv << evtcount << "," << genrhad1->p4().px() << "," << genrhad1->p4().py() << "," << genrhad2->p4().px() << "," << genrhad2->p4().py() << "," << simHit->particleType() << "," << "MuonGEM" << "," << gpos.x() << "," << gpos.y() << "," << gpos.z() << "," << sqrt(gpos.x()*gpos.x()+gpos.y()*gpos.y()) << "," << simHit->energyLoss() << "\n";
  }
  

 ///////////////////////////////////
 // Plot all sim Vertex positions //
 ///////////////////////////////////
  for (simVertex = G4VtxContainer->begin(); simVertex != G4VtxContainer->end(); ++simVertex) {
    float r = sqrt(simVertex->position().x()*simVertex->position().x()+simVertex->position().y()*simVertex->position().y());
    SimVertex_2DXY->Fill(simVertex->position().x(),simVertex->position().y());
    SimVertex_2DRZ->Fill(simVertex->position().z(),r);
  }
  */
}
//define this as a plug-in
DEFINE_FWK_MODULE(SpikedRHadronAnalyzer);
