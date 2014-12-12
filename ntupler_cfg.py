from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('python')

options.register('outfilename',
                 "Run2012A_HT750_NTuple.root",
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Name of output file"
)


import FWCore.ParameterSet.Config as cms

process = cms.Process("Run2012ANtuple")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",

    fileNames = cms.untracked.vstring(
        '/store/data/Run2012A/HT/AOD/22Jan2013-v1/20000/0014C105-AB81-E211-9122-003048678B8E.root'
    )
)

process.TFileService = cms.Service("TFileService",
        fileName = cms.string(options.outfilename)
)

process.Run2012ANtuple = cms.EDAnalyzer('Ntupler',
        src = cms.InputTag("hltCaloJetIDPassed"),
        srcAK5 = cms.InputTag("ak5PFJets"),
        srcCA8 = cms.InputTag("ak7PFJets"),
        srcCSVBTag = cms.InputTag("combinedSecondarVertexBJetTags"),
)

process.p = cms.Path(process.Run2012ANtuple)
