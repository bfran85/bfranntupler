// system include files
#include <memory>
#include <vector>
#include <math.h>
#include <map>
#include <string>
#include <sstream>
#include <vector>
#include <tuple>
#include <cstdlib>
#include <fstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TTree.h"
#include "TFolder.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"


using namespace edm;
using namespace std;
                                                        /////////////////////
                                                        //class declaration//
                                                        ////////////////////

class Ntupler : public edm::EDAnalyzer {
        public:
                explicit Ntupler(const edm::ParameterSet&);
                ~Ntupler();

        private:

                virtual void beginJob() ;
                virtual void ntuplize(const edm::Event&, const edm::EventSetup&);
                virtual void endJob() ;
                int counter;
                edm::Service<TFileService> fs;  //TFile service for Tree generation
                TTree *newtree;                 //Pointer for tree

                int nJets_AK5-PF;
                int nJets_AK7-PF;
                vector<float> ak5jet_px;
                vector<float> ak5jet_py;
                vector<float> ak5jet_pz;
                vector<float> ak5jet_e;
                vector<float> ak5jet_eta;
                vector<float> ak5jet_phi;
                vector<float> ak5jet_pt;
                vector<float> ak7jet_px;
                vector<float> ak7jet_py;
                vector<float> ak7jet_pz;
                vector<float> ak7jet_e;
                vector<float> ak7jet_eta;
                vector<float> ak7jet_phi;
                vector<float> ak5jet_pt;
                vector<float> ak5jet_CSVbdisc;
                vector<float> ak7jet_CSVbdisc;

        // ----------member data ---------------------------
        edm::InputTag src_;
        edm::InputTag src2_;
        edm::InputTag srcAK7_;
        edm::InputTag srcAK5_;


/// constructors and destructor
///
Ntupler::Ntupler(const edm::ParameterSet& iConfig) {
        src_( iConfig.getParameter<edm::InputTag>( "src" ) ),                   // Obtain input
         src2_( iConfig.getParameter<edm::InputTag>( "srcGen" ) )                        // Obtain input
        srcAK7_( iConfig.getParameter<edm::InputTag>("srcAK7") ),
        srcAK5_( iConfig.getParameter<edm::InputTag>("srcAK5") )

        // Declare new tree
        newtree = fs->make<TTree>("DataSetTree","Analysis Tree for Data Set");

        // Create branches for new tree
//      newtree->Branch("nJets_AK5-PF",&nJets_AK5-PF,"nJets_AK5-PF/I");
//      newtree->Branch("nJets_AK7-PF",&nJets_AK7-PF,"nJets_AK7-PF/I");
        newtree->Branch("ak5jet_px",&ak5jet_px);
        newtree->Branch("ak5jet_py",&ak5jet_py);
        newtree->Branch("ak5jet_pz",&ak5jet_pz);
        newtree->Branch("ak5jet_e",&ak5jet_e);
        newtree->Branch("ak5jet_eta",&ak5jet_eta);
        newtree->Branch("ak5jet_phi",&ak5jet_phi);
        newtree->Branch("ak5jet_pt",&ak5jet_pt);
        newtree->Branch("ak7jet_px",&ak7jet_px);
        newtree->Branch("ak7jet_py",&ak7jet_py);
        newtree->Branch("ak7jet_pz",&ak7jet_pz);
        newtree->Branch("ak7jet_e",&ak7jet_e);
        newtree->Branch("ak7jet_eta",&ak7jet_eta);
        newtree->Branch("ak7jet_phi",&ak7jet_phi);
        newtree->Branch("ak7jet_pt",&ak7jet_pt);
//      newtree->Branch("ak5jet_CSVbdisc",&ak5jet_CSVbdisc);
//      newtree->Branch("ak7jet_CSVbdisc",&ak7jet_CSVbdisc);
//      newtree->Branch("ak5jet_DINKObdisc[nJets_AK5-PF]",&ak5jet_DINKObdisc);
//      newtree->Branch("ak7jet_DINKObdisc[nJets_AK7-PF]",&ak7jet_DINKObdisc);

}

Ntupler::~Ntupler()
{

        // do anything here that needs to be done at deconstruction time
        // (e.g. close files, deallocate resources, etc.)

}

// --------- method called for each event ----------
void
Ntupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
        using namespace std;
        using namespace edm;
        using namespace reco;

        // Way to call GenParticle
//      edm::Handle<std::vector<reco::GenParticle>> particles;
//      iEvent.getByLabel( src2_ , particles );
//      const std::vector<reco::GenParticle> & p = *particles;

        // Way to call CALOJETS
//      edm::Handle<std::vector<reco::CaloJet>> jets;
//      iEvent.getByLabel( src_ , jets );
//      const std::vector<reco::CaloJet> & jet = *jets;

        // Way to call ak5PFJETS
        edm::Handle<std::vector<reco::PFJet>> jet5;
        iEvent.getByLabel(srcAK5_ , jet5);
        const std::vector<reco::PFJet> & j5 = *jet5;     // AODSIM

        // Way to call ak7PFJETS
        edm::Handle<std::vector<reco::PFJet>> jet7;
        iEvent.getByLabel(srcAK7_ , jet7);
        const std::vector<reco::PFJet> & j7 = *jet7;

        for (int iJet = 0; iJet < (int)jets->size(); iJet++)
        {   reco::PFJet j5 = jet5->at(iJet);
            reco::PFJet j7 = jet7->at(iJet);

            ak5jet_px.push_back( jet5.px() );
            ak5jet_py.push_back( jet5.py() );
            ak5jet_pz.push_back( jet5.pz() );
            ak5jet_e.push_back( jet5.e() );
            ak5jet_eta.push_back( jet5.eta() );
            ak5jet_phi.push_back( jet5.phi() );
            ak5jet_pt.push_back( jet5.pt() );

            ak7jet_px.push_back( jet7.px() );
            ak7jet_py.push_back( jet7.py() );
            ak7jet_pz.push_back( jet7.pz() );
            ak7jet_e.push_back( jet7.e() );
            ak7jet_eta.push_back( jet7.eta() );
            ak7jet_phi.push_back( jet7.phi() );
            ak7jet_pt.push_back( jet7.pt() );


        }

        newtree->Fill();

        // clear the vectors for next event
        ak5jet_px.clear();
        ak5jet_py.clear();
        ak5jet_pz.clear();
        ak5jet_e.clear();
        ak5jet_eta.clear();
        ak5jet_phi.clear();
        ak5jet_pt.clear();

        ak7jet_px.clear();
        ak7jet_py.clear();
        ak7jet_pz.clear();
        ak7jet_e.clear();
        ak7jet_eta.clear();
        ak7jet_phi.clear();
        ak7jet_pt.clear();


}

void
Ntupler::beginJob()
{
}

void
Ntupler::endJob()
{
}
