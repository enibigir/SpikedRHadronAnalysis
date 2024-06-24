import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag.globaltag = "106X_mcRun3_2021_realistic_v3"

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        ### Currently a ROOT file that contains gluinos, given by Todd
        #'file:/uscms_data/d2/tadams/hscp/fall22a/CMSSW_10_6_30/src/EXO-RunIISummer20UL18GENSIM-00038-stop.root'
        'file:data/EXO-RunIISummer20UL18GENSIM-00010-v3.root'
        #'file:/uscms_data/d2/tadams/hscp/fall22a/CMSSW_10_6_30/src/EXO-RunIISummer20UL18GENSIM-00073-stau.root'
    )
)

process.demo = cms.EDAnalyzer("SpikedRHadronAnalyzer",

#    gen_info = cms.InputTag("genParticles","","DIGI2RAW"),
    gen_info = cms.InputTag("genParticles","","SIM"),

    G4TrkSrc = cms.InputTag("g4SimHits"),
    G4VtxSrc = cms.InputTag("g4SimHits"),

    # Tracker
    TrackerHitsPixelBarrelLowTof = cms.InputTag("g4SimHits","TrackerHitsPixelBarrelLowTof"),
    TrackerHitsPixelBarrelHighTof = cms.InputTag("g4SimHits","TrackerHitsPixelBarrelHighTof"),
    TrackerHitsPixelEndcapLowTof = cms.InputTag("g4SimHits","TrackerHitsPixelEndcapLowTof"),
    TrackerHitsPixelEndcapHighTof = cms.InputTag("g4SimHits","TrackerHitsPixelEndcapHighTof"),
    TrackerHitsTIBLowTof = cms.InputTag("g4SimHits","TrackerHitsTIBLowTof"),
    TrackerHitsTIBHighTof = cms.InputTag("g4SimHits","TrackerHitsTIBHighTof"),
    TrackerHitsTOBLowTof = cms.InputTag("g4SimHits","TrackerHitsTOBLowTof"),
    TrackerHitsTOBHighTof = cms.InputTag("g4SimHits","TrackerHitsTOBHighTof"),
    TrackerHitsTECLowTof = cms.InputTag("g4SimHits","TrackerHitsTECLowTof"),
    TrackerHitsTECHighTof = cms.InputTag("g4SimHits","TrackerHitsTECHighTof"),
    TrackerHitsTIDLowTof = cms.InputTag("g4SimHits","TrackerHitsTIDLowTof"),
    TrackerHitsTIDHighTof = cms.InputTag("g4SimHits","TrackerHitsTIDHighTof"),

    # Calorimiter
    EcalHitsEB = cms.InputTag("g4SimHits","EcalHitsEB"),
    EcalHitsEE = cms.InputTag("g4SimHits","EcalHitsEE"),
    EcalHitsES = cms.InputTag("g4SimHits","EcalHitsES"),
    HcalHits = cms.InputTag("g4SimHits","HcalHits"),

    # Muon
    MuonCSCHits = cms.InputTag("g4SimHits","MuonCSCHits"),
    MuonDTHits = cms.InputTag("g4SimHits","MuonDTHits"),
    MuonRPCHits = cms.InputTag("g4SimHits","MuonRPCHits"),
    MuonGEMHits = cms.InputTag("g4SimHits","MuonGEMHits"),

    bits = cms.InputTag("TriggerResults","","HLT"),
#    objects = cms.InputTag("selectedPatTrigger"),
#    prescales = cms.InputTag("patTrigger"),
    trig_sum = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
    pfmet_reco = cms.InputTag("pfMet","","RECO"),
    calomet_reco = cms.InputTag("caloMet","","RECO"),
    dedx_info = cms.InputTag("dedxHitInfo","","RECO"),
    pf_reco = cms.InputTag("particleFlow"),
    vertex_reco = cms.InputTag("offlinePrimaryVertices"),
    secondaryvertex_reco = cms.InputTag("inclusiveSecondaryVertices"),
#did not work    vertex_reco = cms.InputTag("inclusiveSecondaryVertices"),
    gentrk_reco = cms.InputTag("generalTracks","","RECO"),
    hscp_cand = cms.InputTag("HSCParticleProducer","","HSCPAnalysis"),
    DedxCollection = cms.InputTag("dedxHitInfo","","RECO"),
#    TypeMode = cms.InputTag("dedxHitInfo","","RECO"),
    TypeMode = cms.untracked.int32(0),
    SampleType = cms.untracked.int32(0),
    Period = cms.untracked.string("2018"),
    UseTemplateLayer = cms.untracked.bool(False),
    SkipPixel = cms.untracked.bool(False),
    DeDxTemplate = cms.untracked.string("data/dEdxTemplate_harm2_SO_in_noC_CCC_MG_2017B.root"),
    EnableDeDxCalibration = cms.untracked.bool(False),
    DeDxCalibration = cms.untracked.string("Data13TeVGains_v2.root"),
    DeDxSF_0 = cms.untracked.double(1.0),
    DeDxSF_1 = cms.untracked.double(1.0325),


)

#process.demo.TypeMode=0
#process.analyzer.SampleType=0
#process.analyzer.saveTree=0 #all saved
#process.analyzer.saveGenTree=0
#process.analyzer.DeDxCalibration="Data13TeVGains_v2.root"
#process.analyzer.DeDxTemplate="dEdxTemplate_harm2_SO_in_noC_CCC_MG_2017B.root"

process.p = cms.Path(process.demo)
