//TTre contains collection of GENIE particles 
//Because gst convertion doesn't work for nucleon events
//Note this is not the same particle list from GEANT
//A. Higuera
//ahiguera@central.uh.edu


#ifndef GENIEdumper_module
#define GENIEdumper_module
// LArSoft includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nutools/EventGeneratorBase/evgenbase.h"
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SummaryData/RunData.h"
// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"


// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"

#define MAX_TRACKS 1000

using namespace std;

//========================================================================

namespace DUNE{

class GENIE : public art::EDAnalyzer {
public:

    explicit GENIE(fhicl::ParameterSet const& pset);
    virtual ~GENIE();

    void beginJob();
    void endJob();
    void beginRun(const art::Run& run);
    void analyze(const art::Event& evt);

    void reconfigure(fhicl::ParameterSet const& pset);

    void Process(const art::Event& evt);
    void reset();

private:

    // the parameters we'll read from the .fcl
    std::string fGenieGenModuleLabel;
 
    TTree *fEventTree;
        
    // Event 
    int Event;
    int Run;
    int SubRun;

    //MC truth
    double MC_vertex[4];
    int    MC_npart;  
    int    MC_id[MAX_TRACKS];  
    int    MC_pdg[MAX_TRACKS]; 
    int    MC_mother[MAX_TRACKS];  
    int    MC_statusCode[MAX_TRACKS];
    double MC_startXYZT[MAX_TRACKS][4]; 
    double MC_endXYZT[MAX_TRACKS][4];  
    double MC_startMomentum[MAX_TRACKS][4]; 
    double MC_endMomentum[MAX_TRACKS][4];  
 

 }; // class GENIE


//========================================================================
GENIE::GENIE(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet)
{
    reconfigure(parameterSet);
}
//========================================================================
GENIE::~GENIE(){
  //destructor
}
//========================================================================
void GENIE::reconfigure(fhicl::ParameterSet const& p){

    fGenieGenModuleLabel = p.get<std::string>("GenieGenModuleLabel");
}
//========================================================================
void GENIE::beginJob(){
  cout<<"job begin..."<<endl;

  art::ServiceHandle<art::TFileService> tfs;
  fEventTree = tfs->make<TTree>("Event", "Event Tree from Sim & Reco");
  fEventTree->Branch("eventNo", &Event);
  fEventTree->Branch("runNo", &Run);
  fEventTree->Branch("subRunNo", &SubRun);
  fEventTree->Branch("mc_npart", &MC_npart);  // number of particles 
  fEventTree->Branch("mc_vertex", &MC_vertex, "mc_vertex[4]/D");
  fEventTree->Branch("mc_id", &MC_id, "mc_id[mc_npart]/I");  
  fEventTree->Branch("mc_pdg", &MC_pdg, "mc_pdg[mc_npart]/I"); 
  fEventTree->Branch("mc_statusCode", &MC_statusCode, "mc_statusCode[mc_npart]/I"); 
  fEventTree->Branch("mc_mother", &MC_mother, "mc_mother[mc_npart]/I"); 
  fEventTree->Branch("mc_startXYZT", &MC_startXYZT, "mc_startXYZT[mc_npart][4]/D");  
  fEventTree->Branch("mc_endXYZT", &MC_endXYZT, "mc_endXYZT[mc_npart][4]/D"); 
  fEventTree->Branch("mc_startMomentum", &MC_startMomentum, "mc_startMomentum[mc_npart][4]/D");  
  fEventTree->Branch("mc_endMomentum", &MC_endMomentum, "mc_endMomentum[mc_npart][4]/D"); 

 
}
//========================================================================
void GENIE::endJob(){
     
}
//========================================================================
void GENIE::beginRun(const art::Run& /*run*/){
  mf::LogInfo("GENIE")<<"begin run..."<<endl;
}
//========================================================================
void GENIE::analyze( const art::Event& event ){
    if (event.isRealData()) return;
    reset();

    Event  = event.id().event(); 
    Run    = event.run();
    SubRun = event.subRun();
    Process(event);
    fEventTree->Fill();
}
//========================================================================
void GENIE::Process( const art::Event& event){

    //GENIE particles 
    art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
    std::vector<art::Ptr<simb::MCTruth> > mclist;
    if(event.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
    art::fill_ptr_vector(mclist, mctruthListHandle);

    art::Ptr<simb::MCTruth> mctruth = mclist[0];
    MC_npart = mctruth->NParticles();
    for( int i =0; i<MC_npart; ++i ){ 
       simb::MCParticle particle = mctruth->GetParticle(i);
       MC_id[i] = particle.TrackId();
       MC_pdg[i] = particle.PdgCode();
       MC_mother[i] = particle.Mother();
       MC_statusCode[i] =particle.StatusCode();
       const TLorentzVector& positionStart = particle.Position(0);
       const TLorentzVector& positionEnd   = particle.EndPosition();
       const TLorentzVector& momentumStart = particle.Momentum(0);
       const TLorentzVector& momentumEnd   = particle.EndMomentum();
       //!Save the true vertex as the vertex using primaries ...hmmm do you have another suggestion?
       if( particle.Mother() == 0 ) positionStart.GetXYZT(MC_vertex); 
       positionStart.GetXYZT(MC_startXYZT[i]);
       positionEnd.GetXYZT(MC_endXYZT[i]);
       momentumStart.GetXYZT(MC_startMomentum[i]);
       momentumEnd.GetXYZT(MC_endMomentum[i]);
      
    }
}
//========================================================================
void GENIE::reset(){

   MC_npart =-999; 
   for(int i = 0; i<4; ++i){
      MC_vertex[i] = -999.0;
   }
  
   for(int i=0; i<MAX_TRACKS; ++i) {
       MC_id[i] = -999;
       MC_pdg[i] = -999;
       MC_mother[i] = -999;
       MC_statusCode[i] = -999;
       for(int j=0; j<4; ++j) {
          MC_startXYZT[i][j]      = -999.0;
          MC_endXYZT[i][j]        = -999.0;
          MC_startMomentum[i][j]  = -999.0;
          MC_endMomentum[i][j]    = -999.0;
       }
    }
}
//========================================================================
DEFINE_ART_MODULE(GENIE)

} 

#endif // GENIE
