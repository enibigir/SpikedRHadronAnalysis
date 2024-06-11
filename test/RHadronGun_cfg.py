# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: Configuration/GenProduction/python/EXO-RunIISummer20UL18GENSIM-00010-fragment.py --python_filename EXO-RunIISummer20UL18GENSIM-00010_1_cfg.py --eventcontent RAWSIM --customise SimG4Core/CustomPhysics/Exotica_HSCP_SIM_cfi.customise,Configuration/DataProcessing/Utils.addMonitoring --datatier GEN-SIM --fileout file:EXO-RunIISummer20UL18GENSIM-00010.root --conditions 106X_upgrade2018_realistic_v4 --beamspot Realistic25ns13TeVEarly2018Collision --step GEN,SIM --geometry DB:Extended --era Run2_2018 --no_exec --mc
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.Eras.Era_Run2_2018_cff import Run2_2018

options = VarParsing('python')
options.parseArguments()
process = cms.Process('SIM',Run2_2018)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic25ns13TeVEarly2018Collision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('Configuration/GenProduction/python/EXO-RunIISummer20UL18GENSIM-00010-fragment.py nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    ),
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(1),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(20971520),
    fileName = cms.untracked.string(options.outputFile),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.XMLFromDBSource.label = cms.string("Extended")
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_upgrade2018_realistic_v4', '')


process.generator = cms.EDFilter("FlatRandomEGunProducer",
    PGunParameters = cms.PSet(
		MaxEta = cms.vdouble(0.0, 0.0),
		MinEta = cms.vdouble(0.0, 0.0),
		MinP = cms.untracked.vdouble(0.0, 0.0),
		MinPt = cms.untracked.vdouble(50.0, 50.0),
		PartID = cms.vint32(1000993),
		ParticleCharge = cms.untracked.int32(0),
    ),
    SLHAFileForPythia8 = cms.string('Configuration/Generator/data/HSCP_gluino_1800_SLHA.spc'),
    comEnergy = cms.double(13000.0),
    crossSection = cms.untracked.double(-1),
    hscpFlavor = cms.untracked.string('gluino'),
    massPoint = cms.untracked.int32(1800),
    maxEventsToPrint = cms.untracked.int32(0),
    particleFile = cms.untracked.string('Configuration/Generator/data/particles_gluino_1800_GeV.txt'),
    pdtFile = cms.FileInPath('Configuration/Generator/data/hscppythiapdtgluino1800.tbl'),
    processFile = cms.untracked.string(options.inputFiles[0]),
    pythiaHepMCVerbosity = cms.untracked.bool(True),
    pythiaPylistVerbosity = cms.untracked.int32(10),
    slhaFile = cms.untracked.string('Configuration/Generator/data/HSCP_gluino_1800_SLHA.spc'),
    useregge = cms.bool(False)
)


process.ProductionFilterSequence = cms.Sequence(process.generator)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)