import FWCore.ParameterSet.Config as cms

from PhysicsTools.NanoAOD.common_cff import Var
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
from RecoBTag.ONNXRuntime.pfParticleNetAK4_cff import _pfParticleNetAK4JetTagsAll
from PhysicsTools.NanoAOD.custom_jme_cff import AddParticleNetAK4Scores
from PhysicsTools.NanoAOD.common_cff import *
import CommonTools.PileupAlgos.Puppi_cff


def customize_xcone(process, isData=False):
    process.edTask = cms.Task()
    process.MessageLogger.cerr.FwkReport.reportEvery = 1

    # #################### #
    #         CHS          #
    # #################### #

    process.chs = cms.EDFilter("CandPtrSelector",
        src=cms.InputTag("packedPFCandidates"),
        cut=cms.string("fromPV(0) > 0")
    )
    process.edTask.add(process.chs)

    usePseudoXCone = cms.bool(True)
    process.xconeCHS = cms.EDProducer("XConeProducer",
      src=cms.InputTag("chs"),
      usePseudoXCone=usePseudoXCone,  # use PseudoXCone (faster) or XCone
      NJets = cms.uint32(2),          # number of fatjets
      RJets = cms.double(1.2),        # cone radius of fatjets
      BetaJets = cms.double(2.0),     # conical mesure (beta = 2.0 is XCone default)
      NSubJets = cms.uint32(3),       # number of subjets in each fatjet
      RSubJets = cms.double(0.4),     # cone radius of subjetSrc
      BetaSubJets = cms.double(2.0)   # conical mesure for subjets
    )
    process.edTask.add(getattr(process,"xconeCHS"))

    process.XConeJetCHSTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
      src = cms.InputTag("xconeCHS"),
      cut = cms.string(" pt > 30"), #probably already applied in miniaod
      name = cms.string("xconeCHS"),
      # doc  = cms.string(" "),  #slimmedJetsAK8, i.e. ak8 fat jets for boosted analysis
      singleton = cms.bool(False), # the number of entries is variable
      extension = cms.bool(False), # this is the main table for the jets
      variables = cms.PSet(P4Vars,
          area = Var("jetArea()", float, doc="jet catchment area, for JECs",precision=10),
          nConstituents = Var("numberOfDaughters()","uint8",doc="Number of particles in the jet"),
          rawFactor = Var("1.",float,doc="1 - Factor to get back to raw pT",precision=6),
          chHEF = Var("chargedHadronEnergyFraction()", float, doc="charged Hadron Energy Fraction", precision= 6),
          neHEF = Var("neutralHadronEnergyFraction()", float, doc="neutral Hadron Energy Fraction", precision= 6),
          chEmEF = Var("chargedEmEnergyFraction()", float, doc="charged Electromagnetic Energy Fraction", precision= 6),
          neEmEF = Var("neutralEmEnergyFraction()", float, doc="neutral Electromagnetic Energy Fraction", precision= 6),
          muEF = Var("muonEnergyFraction()", float, doc="muon Energy Fraction", precision= 6),
          btagDeepB = Var("?(bDiscriminator('pfDeepCSVJetTags:probb')+bDiscriminator('pfDeepCSVJetTags:probbb'))>=0?bDiscriminator('pfDeepCSVJetTags:probb')+bDiscriminator('pfDeepCSVJetTags:probbb'):-1",float,doc="DeepCSV b+bb tag discriminator",precision=10),
          btagDeepFlavB = Var("bDiscriminator('pfDeepFlavourJetTags:probb')+bDiscriminator('pfDeepFlavourJetTags:probbb')+bDiscriminator('pfDeepFlavourJetTags:problepb')",float,doc="DeepJet b+bb+lepb tag discriminator",precision=10),
          subJetIdx1 = Var("?nSubjetCollections()>0 && subjets().size()>0?subjets()[0].key():-1", int, doc="index of first subjet"),
          subJetIdx2 = Var("?nSubjetCollections()>0 && subjets().size()>1?subjets()[1].key():-1", int, doc="index of second subjet"),
          subJetIdx3 = Var("?nSubjetCollections()>0 && subjets().size()>2?subjets()[2].key():-1", int, doc="index of third subjet"),
          tau1 = Var("userFloat('tau1')",float, doc="Nsubjettiness (1 axis)",precision=10),
          tau2 = Var("userFloat('tau2')",float, doc="Nsubjettiness (2 axis)",precision=10),
          tau3 = Var("userFloat('tau3')",float, doc="Nsubjettiness (3 axis)",precision=10),
          tau4 = Var("userFloat('tau4')",float, doc="Nsubjettiness (4 axis)",precision=10),
          # msoftdrop = Var("groomedMass('softdropmass')",float, doc="Corrected soft drop mass with CHS",precision=10),
      ),
    )

    process.XConeSubJetCHSTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
      src = cms.InputTag("xconeCHS","SubJets"),
      cut = cms.string(""), #probably already applied in miniaod
      name = cms.string("XConeCHSsubjet"),
      # doc  = cms.string("slimmedJetsAK8, i.e. ak8 fat jets for boosted analysis"),
      singleton = cms.bool(False), # the number of entries is variable
      extension = cms.bool(False), # this is the main table for the jets
      variables = cms.PSet(P4Vars,
          area = Var("jetArea()", float, doc="jet catchment area, for JECs",precision=10),
          nConstituents = Var("numberOfDaughters()","uint8",doc="Number of particles in the jet"),
          rawFactor = Var("1.",float,doc="1 - Factor to get back to raw pT",precision=6),
          chHEF = Var("chargedHadronEnergyFraction()", float, doc="charged Hadron Energy Fraction", precision= 6),
          neHEF = Var("neutralHadronEnergyFraction()", float, doc="neutral Hadron Energy Fraction", precision= 6),
          chEmEF = Var("chargedEmEnergyFraction()", float, doc="charged Electromagnetic Energy Fraction", precision= 6),
          neEmEF = Var("neutralEmEnergyFraction()", float, doc="neutral Electromagnetic Energy Fraction", precision= 6),
          tau1 = Var("userFloat('tau1')",float, doc="Nsubjettiness (1 axis)",precision=10),
          tau2 = Var("userFloat('tau2')",float, doc="Nsubjettiness (2 axis)",precision=10),
          tau3 = Var("userFloat('tau3')",float, doc="Nsubjettiness (3 axis)",precision=10),
          tau4 = Var("userFloat('tau4')",float, doc="Nsubjettiness (4 axis)",precision=10),
          muEF = Var("muonEnergyFraction()", float, doc="muon Energy Fraction", precision= 6),
      )
    )

    process.XConeTopJetCHSTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
      src = cms.InputTag("xconeCHS","TopJets"),
      cut = cms.string(""), #probably already applied in miniaod
      name = cms.string("xconeCHStopjet"),
      # doc  = cms.string("slimmedJetsAK8, i.e. ak8 fat jets for boosted analysis"),
      singleton = cms.bool(False), # the number of entries is variable
      extension = cms.bool(False), # this is the main table for the jets
      variables = cms.PSet(P4Vars,
          area = Var("jetArea()", float, doc="jet catchment area, for JECs",precision=10),
          nConstituents = Var("numberOfDaughters()","uint8",doc="Number of particles in the jet"),
          rawFactor = Var("1.",float,doc="1 - Factor to get back to raw pT",precision=6),
          chHEF = Var("chargedHadronEnergyFraction()", float, doc="charged Hadron Energy Fraction", precision= 6),
          neHEF = Var("neutralHadronEnergyFraction()", float, doc="neutral Hadron Energy Fraction", precision= 6),
          chEmEF = Var("chargedEmEnergyFraction()", float, doc="charged Electromagnetic Energy Fraction", precision= 6),
          neEmEF = Var("neutralEmEnergyFraction()", float, doc="neutral Electromagnetic Energy Fraction", precision= 6),
          tau1 = Var("userFloat('tau1')",float, doc="Nsubjettiness (1 axis)",precision=10),
          tau2 = Var("userFloat('tau2')",float, doc="Nsubjettiness (2 axis)",precision=10),
          tau3 = Var("userFloat('tau3')",float, doc="Nsubjettiness (3 axis)",precision=10),
          tau4 = Var("userFloat('tau4')",float, doc="Nsubjettiness (4 axis)",precision=10),
          muEF = Var("muonEnergyFraction()", float, doc="muon Energy Fraction", precision= 6),
      )
    )

    # #################### #
    #        Puppi         #
    # #################### #

    process.load('CommonTools/PileupAlgos/Puppi_cff')
    process.puppi.candName = cms.InputTag('packedPFCandidates')
    process.puppi.vertexName = cms.InputTag('offlineSlimmedPrimaryVertices')
    process.puppi.clonePackedCands = cms.bool(True)
    process.puppi.useExistingWeights = cms.bool(True)
    process.edTask.add(process.puppi)

    usePseudoXCone = cms.bool(True)
    process.xconePuppi = cms.EDProducer("XConeProducer",
      src=cms.InputTag("puppi"),
      usePseudoXCone=usePseudoXCone,  # use PseudoXCone (faster) or XCone
      NJets = cms.uint32(2),          # number of fatjets
      RJets = cms.double(1.2),        # cone radius of fatjets
      BetaJets = cms.double(2.0),     # conical mesure (beta = 2.0 is XCone default)
      NSubJets = cms.uint32(3),       # number of subjets in each fatjet
      RSubJets = cms.double(0.4),     # cone radius of subjetSrc
      BetaSubJets = cms.double(2.0)   # conical mesure for subjets
    )
    process.edTask.add(getattr(process,"xconePuppi"))

    process.XConeJetPuppiTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
      src = cms.InputTag("xconePuppi"),
      cut = cms.string(" pt > 30"), #probably already applied in miniaod
      name = cms.string("xconePuppi"),
      # doc  = cms.string(" "),  #slimmedJetsAK8, i.e. ak8 fat jets for boosted analysis
      singleton = cms.bool(False), # the number of entries is variable
      extension = cms.bool(False), # this is the main table for the jets
      variables = cms.PSet(P4Vars,
          area = Var("jetArea()", float, doc="jet catchment area, for JECs",precision=10),
          nConstituents = Var("numberOfDaughters()","uint8",doc="Number of particles in the jet"),
          rawFactor = Var("1.",float,doc="1 - Factor to get back to raw pT",precision=6),
          chHEF = Var("chargedHadronEnergyFraction()", float, doc="charged Hadron Energy Fraction", precision= 6),
          neHEF = Var("neutralHadronEnergyFraction()", float, doc="neutral Hadron Energy Fraction", precision= 6),
          chEmEF = Var("chargedEmEnergyFraction()", float, doc="charged Electromagnetic Energy Fraction", precision= 6),
          neEmEF = Var("neutralEmEnergyFraction()", float, doc="neutral Electromagnetic Energy Fraction", precision= 6),
          muEF = Var("muonEnergyFraction()", float, doc="muon Energy Fraction", precision= 6),
          btagDeepB = Var("?(bDiscriminator('pfDeepCSVJetTags:probb')+bDiscriminator('pfDeepCSVJetTags:probbb'))>=0?bDiscriminator('pfDeepCSVJetTags:probb')+bDiscriminator('pfDeepCSVJetTags:probbb'):-1",float,doc="DeepCSV b+bb tag discriminator",precision=10),
          btagDeepFlavB = Var("bDiscriminator('pfDeepFlavourJetTags:probb')+bDiscriminator('pfDeepFlavourJetTags:probbb')+bDiscriminator('pfDeepFlavourJetTags:problepb')",float,doc="DeepJet b+bb+lepb tag discriminator",precision=10),
          subJetIdx1 = Var("?nSubjetCollections()>0 && subjets().size()>0?subjets()[0].key():-1", int, doc="index of first subjet"),
          subJetIdx2 = Var("?nSubjetCollections()>0 && subjets().size()>1?subjets()[1].key():-1", int, doc="index of second subjet"),
          subJetIdx3 = Var("?nSubjetCollections()>0 && subjets().size()>2?subjets()[2].key():-1", int, doc="index of third subjet"),
          tau1 = Var("userFloat('tau1')",float, doc="Nsubjettiness (1 axis)",precision=10),
          tau2 = Var("userFloat('tau2')",float, doc="Nsubjettiness (2 axis)",precision=10),
          tau3 = Var("userFloat('tau3')",float, doc="Nsubjettiness (3 axis)",precision=10),
          tau4 = Var("userFloat('tau4')",float, doc="Nsubjettiness (4 axis)",precision=10),
          # msoftdrop = Var("groomedMass('softdropmass')",float, doc="Corrected soft drop mass with Puppi",precision=10),
      ),
    )

    process.XConeSubJetPuppiTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
      src = cms.InputTag("xconePuppi","SubJets"),
      cut = cms.string(""), #probably already applied in miniaod
      name = cms.string("xconePuppisubjet"),
      # doc  = cms.string("slimmedJetsAK8, i.e. ak8 fat jets for boosted analysis"),
      singleton = cms.bool(False), # the number of entries is variable
      extension = cms.bool(False), # this is the main table for the jets
      variables = cms.PSet(P4Vars,
          area = Var("jetArea()", float, doc="jet catchment area, for JECs",precision=10),
          nConstituents = Var("numberOfDaughters()","uint8",doc="Number of particles in the jet"),
          rawFactor = Var("1.",float,doc="1 - Factor to get back to raw pT",precision=6),
          chHEF = Var("chargedHadronEnergyFraction()", float, doc="charged Hadron Energy Fraction", precision= 6),
          neHEF = Var("neutralHadronEnergyFraction()", float, doc="neutral Hadron Energy Fraction", precision= 6),
          chEmEF = Var("chargedEmEnergyFraction()", float, doc="charged Electromagnetic Energy Fraction", precision= 6),
          neEmEF = Var("neutralEmEnergyFraction()", float, doc="neutral Electromagnetic Energy Fraction", precision= 6),
          tau1 = Var("userFloat('tau1')",float, doc="Nsubjettiness (1 axis)",precision=10),
          tau2 = Var("userFloat('tau2')",float, doc="Nsubjettiness (2 axis)",precision=10),
          tau3 = Var("userFloat('tau3')",float, doc="Nsubjettiness (3 axis)",precision=10),
          tau4 = Var("userFloat('tau4')",float, doc="Nsubjettiness (4 axis)",precision=10),
          muEF = Var("muonEnergyFraction()", float, doc="muon Energy Fraction", precision= 6),
      )
    )

    process.XConeTopJetPuppiTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
      src = cms.InputTag("xconePuppi","TopJets"),
      cut = cms.string(""), #probably already applied in miniaod
      name = cms.string("XConePuppitopjet"),
      # doc  = cms.string("slimmedJetsAK8, i.e. ak8 fat jets for boosted analysis"),
      singleton = cms.bool(False), # the number of entries is variable
      extension = cms.bool(False), # this is the main table for the jets
      variables = cms.PSet(P4Vars,
          area = Var("jetArea()", float, doc="jet catchment area, for JECs",precision=10),
          nConstituents = Var("numberOfDaughters()","uint8",doc="Number of particles in the jet"),
          rawFactor = Var("1.",float,doc="1 - Factor to get back to raw pT",precision=6),
          chHEF = Var("chargedHadronEnergyFraction()", float, doc="charged Hadron Energy Fraction", precision= 6),
          neHEF = Var("neutralHadronEnergyFraction()", float, doc="neutral Hadron Energy Fraction", precision= 6),
          chEmEF = Var("chargedEmEnergyFraction()", float, doc="charged Electromagnetic Energy Fraction", precision= 6),
          neEmEF = Var("neutralEmEnergyFraction()", float, doc="neutral Electromagnetic Energy Fraction", precision= 6),
          tau1 = Var("userFloat('tau1')",float, doc="Nsubjettiness (1 axis)",precision=10),
          tau2 = Var("userFloat('tau2')",float, doc="Nsubjettiness (2 axis)",precision=10),
          tau3 = Var("userFloat('tau3')",float, doc="Nsubjettiness (3 axis)",precision=10),
          tau4 = Var("userFloat('tau4')",float, doc="Nsubjettiness (4 axis)",precision=10),
          muEF = Var("muonEnergyFraction()", float, doc="muon Energy Fraction", precision= 6),
      )
    )

    # #################### #
    #         GEN          #
    # #################### #

    process.packedGenParticlesForJetsNoNu = cms.EDFilter("CandPtrSelector",
        src=cms.InputTag("packedGenParticles"),
        cut=cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16")
    )
    process.edTask.add(process.packedGenParticlesForJetsNoNu)

    usePseudoXCone = cms.bool(True)
    process.xconeGEN = cms.EDProducer("GenXConeProducer",
        src=cms.InputTag("packedGenParticlesForJetsNoNu"),
        usePseudoXCone=usePseudoXCone,  # use PseudoXCone (faster) or XCone
        NJets = cms.uint32(2),          # number of fatjets
        RJets = cms.double(1.2),        # cone radius of fatjets
        BetaJets = cms.double(2.0),     # conical mesure (beta = 2.0 is XCone default)
        NSubJets = cms.uint32(3),       # number of subjets in each fatjet
        RSubJets = cms.double(0.4),     # cone radius of subjetSrc
        BetaSubJets = cms.double(2.0),  # conical mesure for subjets
        doLeptonSpecific = cms.bool(False),
    )
    process.edTask.add(process.xconeGEN)

    process.XConeJetGENTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
      src = cms.InputTag("xconeGEN"),
      cut = cms.string(" pt > 30"), #probably already applied in miniaod
      name = cms.string("xconeGEN"),
      # doc  = cms.string(" "),  #slimmedJetsAK8, i.e. ak8 fat jets for boosted analysis
      singleton = cms.bool(False), # the number of entries is variable
      extension = cms.bool(False), # this is the main table for the jets
      variables = cms.PSet(P4Vars,
          nConstituents = Var("numberOfDaughters()","uint8",doc="Number of particles in the jet"),
          subJetIdx1 = Var("?nSubjetCollections()>0 && subjets().size()>0?subjets()[0].key():-1", int, doc="index of first subjet"),
          subJetIdx2 = Var("?nSubjetCollections()>0 && subjets().size()>1?subjets()[1].key():-1", int, doc="index of second subjet"),
          subJetIdx3 = Var("?nSubjetCollections()>0 && subjets().size()>2?subjets()[2].key():-1", int, doc="index of third subjet"),
          tau1 = Var("userFloat('tau1')",float, doc="Nsubjettiness (1 axis)",precision=10),
          tau2 = Var("userFloat('tau2')",float, doc="Nsubjettiness (2 axis)",precision=10),
          tau3 = Var("userFloat('tau3')",float, doc="Nsubjettiness (3 axis)",precision=10),
          tau4 = Var("userFloat('tau4')",float, doc="Nsubjettiness (4 axis)",precision=10),
          # msoftdrop = Var("groomedMass('softdropmass')",float, doc="Corrected soft drop mass with CHS",precision=10),
      ),
    )

    process.XConeSubJetGENTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
      src = cms.InputTag("xconeGEN","SubJets"),
      cut = cms.string(""), #probably already applied in miniaod
      name = cms.string("xconeGENsubjet"),
      # doc  = cms.string("slimmedJetsAK8, i.e. ak8 fat jets for boosted analysis"),
      singleton = cms.bool(False), # the number of entries is variable
      extension = cms.bool(False), # this is the main table for the jets
      variables = cms.PSet(P4Vars,
          nConstituents = Var("numberOfDaughters()","uint8",doc="Number of particles in the jet"),
          tau1 = Var("userFloat('tau1')",float, doc="Nsubjettiness (1 axis)",precision=10),
          tau2 = Var("userFloat('tau2')",float, doc="Nsubjettiness (2 axis)",precision=10),
          tau3 = Var("userFloat('tau3')",float, doc="Nsubjettiness (3 axis)",precision=10),
          tau4 = Var("userFloat('tau4')",float, doc="Nsubjettiness (4 axis)",precision=10),
      )
    )

    # #################### #
    #         Write        #
    # #################### #

    XConeCHSSequence = cms.Sequence(process.chs + process.xconeCHS + process.XConeJetCHSTable + process.XConeSubJetCHSTable + process.XConeTopJetCHSTable)
    XConePuppiSequence = cms.Sequence(process.puppi + process.xconePuppi + process.XConeJetPuppiTable + process.XConeSubJetPuppiTable + process.XConeTopJetPuppiTable)
    XConeGENSequence = cms.Sequence(process.packedGenParticlesForJetsNoNu + process.xconeGEN + process.XConeJetGENTable + process.XConeSubJetGENTable)

    process.edTask.add(getattr(process,"updatedPatJetsAK8WithDeepInfo"))
    process.edTask.add(getattr(process,"selectedUpdatedPatJetsAK8WithDeepInfo"))
    process.edTask.add(getattr(process,"updatedPatJetsTransientCorrectedAK8WithDeepInfo"))
    process.edTask.add(getattr(process,"patJetCorrFactorsTransientCorrectedAK8WithDeepInfo"))
    process.edTask.add(getattr(process,"pfParticleNetTagInfosAK8WithDeepInfo"))
    process.edTask.add(getattr(process,"pfParticleNetMassRegressionJetTagsAK8WithDeepInfo"))
    process.edTask.add(getattr(process,"patJetCorrFactorsAK8WithDeepInfo"))

    if isData:
      process.edPath = cms.Path(XConeCHSSequence + XConePuppiSequence, process.edTask)
      process.schedule = cms.Schedule(process.edPath,process.nanoAOD_step,process.endjob_step, process.NANOAODoutput_step)
    else:
      process.edPath = cms.Path(XConeCHSSequence + XConePuppiSequence + XConeGENSequence, process.edTask)
      process.schedule = cms.Schedule(process.edPath,process.nanoAOD_step,process.endjob_step,process.NANOAODSIMoutput_step)

    return process
