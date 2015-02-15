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
#include "DataFormats/BTauReco/interface/JetTag.h"
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
            virtual void analyze(const edm::Event&, const edm::EventSetup&);
            virtual void endJob() ;
            int counter;
            edm::Service<TFileService> fs;  //TFile service for Tree generation
            TTree *newtree;                 //Pointer for tree

            int nJets_AK5Calo;
            vector<float> ak5jet_px;
            vector<float> ak5jet_py;
            vector<float> ak5jet_pz;
            vector<float> ak5jet_e;
            vector<float> ak5jet_eta;
            vector<float> ak5jet_phi;
            vector<float> ak5jet_pt;
            vector<float> ak5jet_CSVbdisc;

            // ----------member data ---------------------------
            edm::InputTag srcCal_;
            edm::InputTag srcAK7_;
            edm::InputTag srcAK5_;
            edm::InputTag srcCSVBTag_;
};


/// constructors and destructor
///
Ntupler::Ntupler(const edm::ParameterSet& iConfig) :
        srcCal_( iConfig.getParameter<edm::InputTag>( "srcCal" ) ),                   // Obtain input
        srcAK7_( iConfig.getParameter<edm::InputTag>("srcAK7") ),
        srcAK5_( iConfig.getParameter<edm::InputTag>("srcAK5") ),
        srcCSVBTag_( iConfig.getParameter<edm::InputTag>("srcCSVBTag") )
        {

        // Declare new tree
        newtree = fs->make<TTree>("DataSetTree","Analysis Tree for Data Set");

        // Create branches for new tree
        newtree->Branch("nJets_AK5Calo",&nJets_AK5PF,"nJets_AK5Calo/I");
        newtree->Branch("ak5jet_px[nJets_AK5Calo]",&ak5jet_px,"ak5jet_px[nJets_AK5Calo]/F");
        newtree->Branch("ak5jet_py[nJets_AK5Calo]",&ak5jet_py,"ak5jet_py[nJets_AK5Calo]/F");
        newtree->Branch("ak5jet_pz[nJets_AK5Calo]",&ak5jet_pz,"ak5jet_pz[nJets_AK5Calo]/F");
        newtree->Branch("ak5jet_e[nJets_AK5Calo]",&ak5jet_e,"ak5jet_e[nJets_AK5Calo]/F");
        newtree->Branch("ak5jet_eta[nJets_AK5Calo]",&ak5jet_eta,"ak5jet_eta[nJets_AK5Calo]/F");
        newtree->Branch("ak5jet_phi[nJets_AK5Calo]",&ak5jet_phi,"ak5jet_phi[nJets_AK5Calo]/F");
        newtree->Branch("ak5jet_pt[nJets_AK5Calo]",&ak5jet_pt,"ak5jet_pt[nJets_AK5Calo]/F");
        newtree->Branch("ak5jet_CSVbdisc[nJets_AK5Calo]",&ak5jet_CSVbdisc,"ak5jet_CSVbdisc[nJets_AK5Calo]/F");

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

        // Way to call CALOJETS
        edm::Handle<std::vector<reco::CaloJet>> cjet;
        iEvent.getByLabel( src_ , cjet );

        // Way to call ak5PFJETS
//      edm::Handle<std::vector<reco::PFJet>> jet5;
//      iEvent.getByLabel(srcAK5_ , jet5);

        // Way to call ak7PFJETS
//      edm::Handle<std::vector<reco::PFJet>> jet7;
//      iEvent.getByLabel(srcAK7_ , jet7);

        // Way to clal BTag Collection    
        edm::Handle<reco::JetTagCollection> bTagHandle;
        iEvent.getByLabel(srcCSVBTag_, bTagHandle);
        const reco::JetTagCollection & bTags = *( bTagHandle.product() );

//      # Loop over gets and study b tag info.
//      for (int i = 0, i != bTags.size(), ++i) {
//          cout<<" Jet "<< i
//               <<" has b tag discriminator = "<<bTags[i].second
//               << " and jet Pt = "<<bTags[i].first->pt()<<endl;

        for (int iJet = 0; iJet < (int)cjet->size(); iJet++)
        {   reco::CaloJet cj = cjet->at(iJet);
            ak5jet_px.push_back( cj.px() );
            ak5jet_py.push_back( cj.py() );
            ak5jet_pz.push_back( cj.pz() );
            ak5jet_e.push_back( cj.energy() );
            ak5jet_eta.push_back( cj.eta() );
            ak5jet_phi.push_back( cj.phi() );
            ak5jet_pt.push_back( cj.pt() );
            ak5jet_CSVbdisc.push_back( bTags[iJet].second );  
        }

        nJets_AK5Calo = cjets->size();
        newtree->Fill();

        // clear the vectors for next event
        ak5jet_px.clear();
        ak5jet_py.clear();
        ak5jet_pz.clear();
        ak5jet_e.clear();
        ak5jet_eta.clear();
        ak5jet_phi.clear();
        ak5jet_pt.clear();
        ak5jet_CSVbdisc.clear();

}

void
Ntupler::beginJob()
{
}

void
Ntupler::endJob()
{
}

DEFINE_FWK_MODULE(Ntupler);
