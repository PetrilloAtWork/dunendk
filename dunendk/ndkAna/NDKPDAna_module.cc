//Module analyzer to look into Photon detector system information
//Ana TTree 
//ahiguera@central.uh.edu


#ifndef NDKPDAna_Module
#define NDKPDAna_Module

// LArSoft includes
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"  
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larsim/MCCheater/BackTracker.h"
#include "larsim/MCCheater/PhotonBackTracker.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/OpHit.h"

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
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"


// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"

#define MAX_TRACKS 20000
#define MAX_FLASHES 2000

using namespace std;

//========================================================================

namespace DUNE{

class NDKPDAna : public art::EDAnalyzer {
public:

    explicit NDKPDAna(fhicl::ParameterSet const& pset);
    virtual ~NDKPDAna();

    void beginJob();
    void endJob();
    void beginRun(const art::Run& run);
    void analyze(const art::Event& evt);

    void reconfigure(fhicl::ParameterSet const& pset);

    void Process(const art::Event& evt, bool &isFiducial);
    bool insideFV(double vertex[4]);
    void reset();

private:

    // the parameters we'll read from the .fcl
    std::string fMCTruthModuleLabel;
    std::string fOpHitModuleLabel;
    std::string fOpFlashModuleLabel;
    bool	fSaveMCTree; 
 
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
    double MC_t0[MAX_TRACKS]; 
    int    MC_mother[MAX_TRACKS];  
    double MC_startXYZT[MAX_TRACKS][4]; 
    double MC_endXYZT[MAX_TRACKS][4];  
    double MC_startMomentum[MAX_TRACKS][4]; 
    double MC_endMomentum[MAX_TRACKS][4];  

    int    n_ophits;
    double ophit_pe[MAX_FLASHES];
    double ophit_time[MAX_FLASHES];
    double ophit_abstime[MAX_FLASHES];
    double ophit_width[MAX_FLASHES];
    int    ophit_frame[MAX_FLASHES];
    int    ophit_ch[MAX_FLASHES];
    //int    ophit_pdg[MAX_FLASHES];
    int    ophit_npart[MAX_FLASHES];
    //double ophit_XYZ[MAX_FLASHES][3];

    int    n_ophits_k;
    double ophit_pe_k[MAX_FLASHES];
    double ophit_time_k[MAX_FLASHES];
    double ophit_abstime_k[MAX_FLASHES];
    double ophit_width_k[MAX_FLASHES];
    int    ophit_frame_k[MAX_FLASHES];
    int    ophit_ch_k[MAX_FLASHES];
    //double ophit_XYZ_k[MAX_FLASHES][3];

    int    n_ophits_mu;
    double ophit_pe_mu[MAX_FLASHES];
    double ophit_time_mu[MAX_FLASHES];
    double ophit_abstime_mu[MAX_FLASHES];
    double ophit_width_mu[MAX_FLASHES];
    int    ophit_frame_mu[MAX_FLASHES];
    int    ophit_ch_mu[MAX_FLASHES];
    //double ophit_XYZ_mu[MAX_FLASHES][3];

    int    n_ophits_e;
    double ophit_pe_e[MAX_FLASHES];
    double ophit_time_e[MAX_FLASHES];
    double ophit_abstime_e[MAX_FLASHES];
    double ophit_width_e[MAX_FLASHES];
    int    ophit_frame_e[MAX_FLASHES];
    int    ophit_ch_e[MAX_FLASHES];
    //double ophit_XYZ_e[MAX_FLASHES][3];

    int    n_flashes;
    double flash_time[MAX_FLASHES];
    double flash_pe[MAX_FLASHES];
    double flash_ycenter[MAX_FLASHES];
    double flash_zcenter[MAX_FLASHES];
    double flash_ywidth[MAX_FLASHES];
    double flash_zwidth[MAX_FLASHES];
    double flash_timewidth[MAX_FLASHES];
    double flash_abstime[MAX_FLASHES];
    int    flash_frame[MAX_FLASHES];

    double fFidVolCutX;
    double fFidVolCutY;
    double fFidVolCutZ;

    double fFidVolXmin;
    double fFidVolXmax;
    double fFidVolYmin;
    double fFidVolYmax;
    double fFidVolZmin;
    double fFidVolZmax;

    detinfo::DetectorProperties const *detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    detinfo::DetectorClocks const *ts = lar::providerFrom<detinfo::DetectorClocksService>();
    double XDriftVelocity = detprop->DriftVelocity()*1e-3; //cm/ns
    double WindowSize     = detprop->NumberTimeSamples() * ts->TPCClock().TickPeriod() * 1e3;
    art::ServiceHandle<geo::Geometry> geom;

 }; // class NDKPDAna


//========================================================================
NDKPDAna::NDKPDAna(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet)
{
    reconfigure(parameterSet);
}
//========================================================================
NDKPDAna::~NDKPDAna(){
  //destructor
}
//========================================================================
void NDKPDAna::reconfigure(fhicl::ParameterSet const& p){

    fMCTruthModuleLabel  = p.get<std::string>("MCTruthModuleLabel");
    fOpFlashModuleLabel  = p.get<std::string>("OpFlashModuleLabel");
    fOpHitModuleLabel    = p.get<std::string>("OpHitModuleLabel");
    fFidVolCutX          = p.get<double>("FidVolCutX");
    fFidVolCutY          = p.get<double>("FidVolCutY");
    fFidVolCutZ          = p.get<double>("FidVolCutZ");
}
//========================================================================
void NDKPDAna::beginJob(){
  cout<<"job begin..."<<endl;
  // Get geometry.
  auto const* geo = lar::providerFrom<geo::Geometry>();
  // Define histogram boundaries (cm).
  // For now only draw cryostat=0.
  double minx = 1e9;
  double maxx = -1e9;
  double miny = 1e9;
  double maxy = -1e9;
  double minz = 1e9;
  double maxz = -1e9;
  for (size_t i = 0; i<geo->NTPC(); ++i){
    double local[3] = {0.,0.,0.};
    double world[3] = {0.,0.,0.};
    const geo::TPCGeo &tpc = geo->TPC(i);
    tpc.LocalToWorld(local,world);
    if (minx>world[0]-geo->DetHalfWidth(i))
      minx = world[0]-geo->DetHalfWidth(i);
    if (maxx<world[0]+geo->DetHalfWidth(i))
      maxx = world[0]+geo->DetHalfWidth(i);
    if (miny>world[1]-geo->DetHalfHeight(i))
      miny = world[1]-geo->DetHalfHeight(i);
    if (maxy<world[1]+geo->DetHalfHeight(i))
      maxy = world[1]+geo->DetHalfHeight(i);
    if (minz>world[2]-geo->DetLength(i)/2.)
      minz = world[2]-geo->DetLength(i)/2.;
    if (maxz<world[2]+geo->DetLength(i)/2.)
      maxz = world[2]+geo->DetLength(i)/2.;
  }

  fFidVolXmin = minx + fFidVolCutX;
  fFidVolXmax = maxx - fFidVolCutX;
  fFidVolYmin = miny + fFidVolCutY;
  fFidVolYmax = maxy - fFidVolCutY;
  fFidVolZmin = minz + fFidVolCutZ;
  fFidVolZmax = maxz - fFidVolCutZ;

  cout<<"Fiducial volume:"<<"\n"
	   <<fFidVolXmin<<"\t< x <\t"<<fFidVolXmax<<"\n"
	   <<fFidVolYmin<<"\t< y <\t"<<fFidVolYmax<<"\n"
	   <<fFidVolZmin<<"\t< z <\t"<<fFidVolZmax<<"\n";

  
  art::ServiceHandle<art::TFileService> tfs;
  fEventTree = tfs->make<TTree>("Event", "Event Tree from Sim & Reco");
  
  fEventTree->Branch("eventNo", &Event);
  fEventTree->Branch("runNo", &Run);
  fEventTree->Branch("subRunNo", &SubRun);
  fEventTree->Branch("mc_vertex", &MC_vertex, "mc_vertex[4]/D");
  fEventTree->Branch("mc_npart", &MC_npart);  // number of particles 
  fEventTree->Branch("mc_id", &MC_id, "mc_id[mc_npart]/I");  
  fEventTree->Branch("mc_t0", &MC_t0, "mc_t0[mc_npart]/D");  
  fEventTree->Branch("mc_pdg", &MC_pdg, "mc_pdg[mc_npart]/I"); 
  fEventTree->Branch("mc_mother", &MC_mother, "mc_mother[mc_npart]/I"); 
  fEventTree->Branch("mc_startXYZT", &MC_startXYZT, "mc_startXYZT[mc_npart][4]/D");  
  fEventTree->Branch("mc_endXYZT", &MC_endXYZT, "mc_endXYZT[mc_npart][4]/D"); 
  fEventTree->Branch("mc_startMomentum", &MC_startMomentum, "mc_startMomentum[mc_npart][4]/D");  
  fEventTree->Branch("mc_endMomentum", &MC_endMomentum, "mc_endMomentum[mc_npart][4]/D"); 

  fEventTree->Branch("n_ophits", &n_ophits);
  fEventTree->Branch("ophit_time", &ophit_time,"ophit_time[n_ophits]/D");
  fEventTree->Branch("ophit_pe", &ophit_pe,"ophit_pe[n_ophits]/D");
  //fEventTree->Branch("ophit_XYZ", &ophit_XYZ,"ophit_XYZ[n_ophits][3]/D");
  fEventTree->Branch("ophit_width", &ophit_width,"ophit_width[n_ophits]/D");
  fEventTree->Branch("ophit_abstime", &ophit_abstime, "ophit_abstime[n_ophits]/D");
  fEventTree->Branch("ophit_frame", &ophit_frame, "ophit_frame[n_ophits]/I");
  fEventTree->Branch("ophit_ch", &ophit_ch, "ophit_ch[n_ophits]/I");
  //fEventTree->Branch("ophit_pdg", &ophit_pdg, "ophit_pdg[n_ophits]/I");
  fEventTree->Branch("ophit_npart", &ophit_npart, "ophit_npart[n_ophits]/I");

  
  fEventTree->Branch("n_ophits_k", &n_ophits_k);
  fEventTree->Branch("ophit_time_k", &ophit_time_k,"ophit_time_k[n_ophits_k]/D");
  fEventTree->Branch("ophit_pe_k", &ophit_pe_k,"ophit_pe_k[n_ophits_k]/D");
  //fEventTree->Branch("ophit_XYZ_k", &ophit_XYZ_k,"ophit_XYZ_k[n_ophits_k][3]/D");
  fEventTree->Branch("ophit_width_k", &ophit_width_k,"ophit_width_k[n_ophits_k]/D");
  fEventTree->Branch("ophit_abstime_k", &ophit_abstime_k, "ophit_abstime_k[n_ophits_k]/D");
  fEventTree->Branch("ophit_frame_k", &ophit_frame_k, "ophit_frame_k[n_ophits_k]/I");
  fEventTree->Branch("ophit_ch_k", &ophit_ch_k, "ophit_ch_k[n_ophits_k]/I");

  fEventTree->Branch("n_ophits_mu", &n_ophits_mu);
  fEventTree->Branch("ophit_time_mu", &ophit_time_mu,"ophit_time_mu[n_ophits_mu]/D");
  fEventTree->Branch("ophit_pe_mu", &ophit_pe_mu,"ophit_pe_mu[n_ophits_mu]/D");
  //fEventTree->Branch("ophit_XYZ_mu", &ophit_XYZ_mu,"ophit_XYZ_mu[n_ophits_mu][3]/D");
  fEventTree->Branch("ophit_width_mu", &ophit_width_mu,"ophit_width_mu[n_ophits_mu]/D");
  fEventTree->Branch("ophit_abstime_mu", &ophit_abstime_mu, "ophit_abstime_mu[n_ophits_mu]/D");
  fEventTree->Branch("ophit_frame_mu", &ophit_frame_mu, "ophit_frame_mu[n_ophits_mu]/I");
  fEventTree->Branch("ophit_ch_mu", &ophit_ch_mu, "ophit_ch_mu[n_ophits_mu]/I");
  
  fEventTree->Branch("n_ophits_e", &n_ophits_e);
  fEventTree->Branch("ophit_time_e", &ophit_time_e,"ophit_time_e[n_ophits_e]/D");
  fEventTree->Branch("ophit_pe_e", &ophit_pe_e,"ophit_pe_e[n_ophits_e]/D");
  //fEventTree->Branch("ophit_XYZ_e", &ophit_XYZ_e,"ophit_XYZ_e[n_ophits_e][3]/D");
  fEventTree->Branch("ophit_width_e", &ophit_width_e,"ophit_width_e[n_ophits_e]/D");
  fEventTree->Branch("ophit_abstime_e", &ophit_abstime_e, "ophit_abstime_e[n_ophits_e]/D");
  fEventTree->Branch("ophit_frame_e", &ophit_frame_e, "ophit_frame_e[n_ophits_e]/I");
  fEventTree->Branch("ophit_ch_e", &ophit_ch_e, "ophit_ch_e[n_ophits_e]/I");
  
  fEventTree->Branch("n_flashes", &n_flashes);
  fEventTree->Branch("flash_time", &flash_time,"flash_time[n_flashes]/D");
  fEventTree->Branch("flash_pe", &flash_pe,"flash_pe[n_flashes]/D");
  fEventTree->Branch("flash_ycenter", &flash_ycenter,"flash_ycenter[n_flashes]/D");
  fEventTree->Branch("flash_zcenter", &flash_zcenter,"flash_zcenter[n_flashes]/D");
  fEventTree->Branch("flash_ywidth", &flash_ywidth,"flash_ywidth[n_flashes]/D");
  fEventTree->Branch("flash_zwidth", &flash_zwidth,"flash_zwidth[n_flashes]/D");
  fEventTree->Branch("flash_timewidth", &flash_timewidth, "flash_timewidth[n_flashes]/D");
  fEventTree->Branch("flash_abstime", &flash_abstime, "flash_abstime[n_flashes]/D");
  fEventTree->Branch("flash_frame", &flash_frame, "flash_frame[n_flashes]/I");


}
//========================================================================
void NDKPDAna::endJob(){
}
//========================================================================
void NDKPDAna::beginRun(const art::Run& /*run*/){
  mf::LogInfo("NDKPDAna")<<"begin run..."<<endl;
}
//========================================================================
void NDKPDAna::analyze( const art::Event& event ){
    if (event.isRealData()) return;
    reset();

    Event  = event.id().event(); 
    Run    = event.run();
    SubRun = event.subRun();
    bool isFiducial = false;
    Process(event, isFiducial);
    if(isFiducial) fEventTree->Fill();
    
}
//========================================================================
void NDKPDAna::Process( const art::Event& event, bool &isFiducial){
    
    art::ServiceHandle<cheat::BackTracker> bt;
    const sim::ParticleList& plist = bt->ParticleList();
    simb::MCParticle *particle=0;
    int i=0; // particle index
    MC_npart = plist.size();
    std::vector<int> kaon_id;
    std::vector<int> muon_id;
    std::vector<int> michel_id;
    int tmp_muon_id =-999;
    int tmp_kaon_id =-999;
    if( MC_npart > MAX_TRACKS ) return;
    for( sim::ParticleList::const_iterator ipar = plist.begin(); ipar!=plist.end(); ++ipar){
       particle = ipar->second;
       MC_t0[i] = particle->T();
       MC_id[i] = particle->TrackId();
       MC_pdg[i] = particle->PdgCode();
       MC_mother[i] = particle->Mother();
       //cout<<MC_pdg[i]<<" "<<MC_id[i]<<" "<<particle->Mother()<<endl;
       if( particle->PdgCode() == 321 && particle->Mother() == 0 )  { tmp_kaon_id = particle->TrackId(); kaon_id.push_back(particle->TrackId()); } 
       else if ( particle->PdgCode() == -13 && particle->Mother() == tmp_kaon_id ){tmp_muon_id = particle->TrackId(); muon_id.push_back(particle->TrackId());}
       else if ( particle->PdgCode() == -11 && particle->Mother() == tmp_muon_id ) michel_id.push_back(particle->TrackId());
       
       const TLorentzVector& positionStart = particle->Position(0);
       const TLorentzVector& positionEnd   = particle->EndPosition();
       const TLorentzVector& momentumStart = particle->Momentum(0);
       const TLorentzVector& momentumEnd   = particle->EndMomentum();
       //!Save the true vertex as the vertex using primaries ...hmmm do you have another suggestion?
       if( particle->Mother() == 0 ) positionStart.GetXYZT(MC_vertex); 
       positionStart.GetXYZT(MC_startXYZT[i]);
       positionEnd.GetXYZT(MC_endXYZT[i]);
       momentumStart.GetXYZT(MC_startMomentum[i]);
       momentumEnd.GetXYZT(MC_endMomentum[i]);
       
       ++i; //paticle index
       
    }
     
    //cout<<MC_t0[0]<<" "<<kaon_id<<" "<<muon_id<<" "<<michel_id<<endl; 
    //isFiducial =insideFV( MC_vertex );
    //cout<<"fiducial "<<isFiducial<<endl;
    //if( !isFiducial ) return;
    isFiducial = true;
    /* 
    //========================================================================
    //========================================================================
    // Reco  stuff
    //========================================================================
    //========================================================================
    //PDS info...
    art::ServiceHandle<cheat::PhotonBackTracker> Photon_bt;
    art::Handle< std::vector< recob::OpHit > > OpHitHandle;
    if( !event.getByLabel(fOpHitModuleLabel, OpHitHandle)) return;
    std::vector< art::Ptr< recob::OpHit >> opHit_list;
    art::fill_ptr_vector(opHit_list, OpHitHandle);
    n_ophits = opHit_list.size();
    //int wtrkID=0;  //with track ID
    //int wotrkID =0;// with out track ID (noise)?
    for(int i = 0; i < n_ophits; ++i){
       art::Ptr<recob::OpHit> op_hit = opHit_list[i]; 
       ophit_time[i]= op_hit->PeakTime();
       ophit_abstime[i] = op_hit->PeakTimeAbs();
       ophit_pe[i] = op_hit->PE();
       ophit_ch[i] = op_hit->OpChannel();
       ophit_frame[i] = op_hit->Frame();
       ophit_width[i] = op_hit->Width();

       std::vector<sim::TrackSDP> trkSDPs = Photon_bt->OpHitToEveSDPs(op_hit);
       //cout<<"TrackIDs associated to this ophit: "<<trkSDPs.size()<<endl; 
       ophit_npart[i] = (int)trkSDPs.size();
       //bool pass_MCpart = false;
       //int tmp_pdg =-999;
       if( trkSDPs.size() != 0 ){ 
         //wtrkID ++;
         for(const sim::TrackSDP trkSDP : trkSDPs){
            //cout<<trkSDP.trackID<<endl;
            if( trkSDP.trackID != 0){//Why ??
              const simb::MCParticle* particle = Photon_bt->TrackIDToParticle(trkSDP.trackID); 
              //cout<<"Which track ID: "<<trkSDP.trackID<<" w/pdg: "<<particle->PdgCode()<<" particle T: "<<particle->T()<<endl;
              if( particle->PdgCode() == 321 )
              cout<<"Which track ID: "<<trkSDP.trackID<<" w/pdg: "<<particle->PdgCode()<<" particle T: "<<particle->T()<<endl;
              //tmp_pdg = particle->PdgCode();
              //pass_MCpart = true;
            }
         }
       }
       / *
       else wotrkID ++; 
       ophit_pdg[i] = tmp_pdg; //temporally 
       //cout<<"peak time "<<op_hit->PeakTime()<<" PE "<<op_hit->PE()<<" CH "<<op_hit->OpChannel()<<" Frame "<<op_hit->Frame()<<endl;
       std::vector<double> tmp_xyz;
       if( pass_MCpart ){ 
         tmp_xyz = Photon_bt->OpHitToXYZ(op_hit);
         ophit_XYZ[i][0]=tmp_xyz.at(0);
         ophit_XYZ[i][1]=tmp_xyz.at(1);
         ophit_XYZ[i][2]=tmp_xyz.at(2);
       }
        * /
    }
    cout<<"How many reco optical hits: "<<n_ophits<<endl;
    //cout<<"these many optical hits have a track ID associated: "<<wtrkID<<endl;
    //cout<<"these many optical hits do not have a track ID associated (noise): "<<wotrkID<<endl;
 
    //============================================================================== 
    //============================================================================== 
    int k=0; 
    cout<<kaon_id.size()<<" "<<muon_id.size()<<" "<<michel_id.size()<<endl; 
    cout<<kaon_id[0]<<" "<<muon_id[0]<<" "<<michel_id[0]<<endl; 
    if( kaon_id.size() != 0){
      const std::vector<std::vector<art::Ptr<recob::OpHit>>> opHits_k = Photon_bt->TrackIDsToOpHits(opHit_list,kaon_id); 
      n_ophits_k = opHits_k[0].size(); //there is only one kaon per event
      cout<<"How many ophit from the kaon: "<<opHits_k[0].size()<<endl;
      std::vector<art::Ptr<recob::OpHit>>::const_iterator it_k_ophit;// = opHits_k[0].begin();
      for( it_k_ophit = opHits_k[0].begin(); it_k_ophit !=  opHits_k[0].end(); ++ it_k_ophit  ){
         ophit_time_k[k]= (*it_k_ophit)->PeakTime();
         ophit_abstime_k[k] = (*it_k_ophit)->PeakTimeAbs();
         ophit_pe_k[k] = (*it_k_ophit)->PE();
         ophit_ch_k[k] = (*it_k_ophit)->OpChannel();
         ophit_frame_k[k] = (*it_k_ophit)->Frame();
         ophit_width_k[k] = (*it_k_ophit)->Width();
         k++;

       }
    }
    //============================================================================== 
    //============================================================================== 
    k =0; //reset index
    if( muon_id.size() != 0){
      const std::vector<std::vector<art::Ptr<recob::OpHit>>> opHits_mu = Photon_bt->TrackIDsToOpHits(opHit_list,muon_id); 
      n_ophits_mu = opHits_mu[0].size(); //there is only one muon per event
      cout<<"How many ophit from the muon: "<<opHits_mu[0].size()<<endl;
      std::vector<art::Ptr<recob::OpHit>>::const_iterator it_mu_ophit;
      for( it_mu_ophit = opHits_mu[0].begin(); it_mu_ophit !=  opHits_mu[0].end(); ++ it_mu_ophit  ){
         ophit_time_mu[k]= (*it_mu_ophit)->PeakTime();
         ophit_abstime_mu[k] = (*it_mu_ophit)->PeakTimeAbs();
         ophit_pe_mu[k] = (*it_mu_ophit)->PE();
         ophit_ch_mu[k] = (*it_mu_ophit)->OpChannel();
         ophit_frame_mu[k] = (*it_mu_ophit)->Frame();
         ophit_width_mu[k] = (*it_mu_ophit)->Width();
         k++;
      }
    }
    //============================================================================== 
    //============================================================================== 
    k=0;
    if( michel_id.size()!=0){
      const std::vector<std::vector<art::Ptr<recob::OpHit>>> opHits_e = Photon_bt->TrackIDsToOpHits(opHit_list,michel_id); 
      n_ophits_e = opHits_e[0].size(); //there is only one michel per event
      cout<<"How many ophit from the michel: "<<opHits_e[0].size()<<endl;
      std::vector<art::Ptr<recob::OpHit>>::const_iterator it_e_ophit;// = opHits_k[0].begin();
      for( it_e_ophit = opHits_e[0].begin(); it_e_ophit !=  opHits_e[0].end(); ++ it_e_ophit  ){
         ophit_time_e[k]= (*it_e_ophit)->PeakTime();
         ophit_abstime_e[k] = (*it_e_ophit)->PeakTimeAbs();
         ophit_pe_e[k] = (*it_e_ophit)->PE();
         ophit_ch_e[k] = (*it_e_ophit)->OpChannel();
         ophit_frame_e[k] = (*it_e_ophit)->Frame();
         ophit_width_e[k] = (*it_e_ophit)->Width();
         k++;
      }
    }

    //========================================================================
    //========================================================================
    art::Handle< std::vector<recob::OpFlash> > flashListHandle;
    std::vector<art::Ptr<recob::OpFlash> > flashlist;
    if(!event.getByLabel(fOpFlashModuleLabel,flashListHandle)) return;
    art::fill_ptr_vector(flashlist, flashListHandle); 
    
    art::FindManyP<recob::OpHit> flash_ophits(flashListHandle, event, fOpFlashModuleLabel);
    n_flashes = flashlist.size();
    cout<<"How many flashes "<<n_flashes<<endl;
    for(int i=0; i<n_flashes && i <MAX_FLASHES; ++i ){
       flash_time[i] = flashlist[i]->Time();
       flash_pe[i] = flashlist[i]->TotalPE();
       flash_ycenter[i] = flashlist[i]->YCenter();
       flash_zcenter[i] = flashlist[i]->ZCenter();
       flash_ywidth[i] = flashlist[i]->YWidth();
       flash_zwidth[i] = flashlist[i]->ZWidth();
       flash_timewidth[i] = flashlist[i]->TimeWidth();
       flash_abstime[i] = flashlist[i]->AbsTime();
       flash_frame[i] = flashlist[i]->Frame();
    }
    / *
    std::vector<art::Ptr<recob::OpHit>> flash1_Hits;  //one flash
    std::vector<art::Ptr<recob::OpHit>> flash2_Hits;  //two flash
    if( n_flashes == 1 ){
       std::vector<art::Ptr<recob::OpHit>> flash1_Hits = flash_ophits.at(0);  
       cout<<"one "<<flash1_Hits.size()<<endl;
      //for( size_t j =0; j<flash1_Hits.size(); ++j){
      //   cout<<"time "<<flash1_Hits[j]->PeakTime()<<endl;
      //}
    }
    if( n_flashes == 2 ){
       std::vector<art::Ptr<recob::OpHit>> flash1_Hits = flash_ophits.at(0);  
       std::vector<art::Ptr<recob::OpHit>> flash2_Hits = flash_ophits.at(1);  
       cout<<"one "<<flash1_Hits.size()<<" 2 "<<flash2_Hits.size()<<endl;
       for( size_t j =0; j<flash2_Hits.size(); ++j){
         cout<<"time "<<flash1_Hits[j]->PeakTime()<<" other "<<flash2_Hits[j]->PeakTime()<<endl;
       }

    }
    * /
    */
} 
//========================================================================
bool NDKPDAna::insideFV( double vertex[4]){ 

     double x = vertex[0];
     double y = vertex[1];
     double z = vertex[2];

     if (x>fFidVolXmin && x<fFidVolXmax&&
	 y>fFidVolYmin && y<fFidVolYmax&&
	 z>fFidVolZmin && z<fFidVolZmax)
       return true;
     else
       return false;
}
//========================================================================
//========================================================================
void NDKPDAna::reset(){

   MC_npart =0; 
   for(int i = 0; i<4; ++i){
      MC_vertex[i] = -999.0;
   }
  
   for(int i=0; i<MAX_TRACKS; ++i) {
       MC_id[i] = -999;
       MC_pdg[i] = -999;
       MC_mother[i] = -999;
       MC_t0[i] =-999.0;
       for(int j=0; j<4; ++j) {
          MC_startXYZT[i][j]      = -999.0;
          MC_endXYZT[i][j]        = -999.0;
          MC_startMomentum[i][j]  = -999.0;
          MC_endMomentum[i][j]    = -999.0;
       }
    }
    n_ophits_k =0;
    n_ophits_mu =0;
    n_ophits_e =0;
    n_flashes =0;
    for(int i=0; i<MAX_FLASHES; ++i){
       flash_time[i] =-999.0;
       flash_pe[i] =-999.0;
       flash_ycenter[i] =-999.0;
       flash_zcenter[i] =-999.0;
       flash_ywidth[i] =-999.0;
       flash_zwidth[i] =-999.0;
       flash_timewidth[i] = -999.0;
       flash_abstime[i] =-999.0;
       flash_frame[i] =-999;
   
       ophit_time[i] =-999.0;
       ophit_pe[i] =-999.0;
       ophit_abstime[i] =-999.0;
       ophit_frame[i] =-999;
       ophit_ch[i] = -999;
       //ophit_pdg[i] = -999;
       ophit_width[i] = -999.;
       ophit_npart[i] = -999;

       ophit_time_k[i] =-999.0;
       ophit_pe_k[i] =-999.0;
       ophit_abstime_k[i] =-999.0;
       ophit_frame_k[i] =-999;
       ophit_ch_k[i] = -999;
       ophit_width_k[i] = -999.;
 
       ophit_time_mu[i] =-999.0;
       ophit_pe_mu[i] =-999.0;
       ophit_abstime_mu[i] =-999.0;
       ophit_frame_mu[i] =-999;
       ophit_ch_mu[i] = -999;
       ophit_width_mu[i] = -999.;
 
       ophit_time_e[i] =-999.0;
       ophit_pe_e[i] =-999.0;
       ophit_abstime_e[i] =-999.0;
       ophit_frame_e[i] =-999;
       ophit_ch_e[i] = -999;
       ophit_width_e[i] = -999.;
 
    }
}
//========================================================================
DEFINE_ART_MODULE(NDKPDAna)

} 

#endif // NDKPDAna_module
