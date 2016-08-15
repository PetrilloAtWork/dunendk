#ifndef NDKAna_Module
#define NDKAna_Module

// LArSoft includes
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardata/RecoBase/Track.h"
#include "SimulationBase/MCParticle.h"
#include "SimulationBase/MCTruth.h"
#include "larcore/SimpleTypesAndConstants/geo_types.h"
#include "larsim/MCCheater/BackTracker.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/RecoBase/Vertex.h"
#include "lardata/RawData/ExternalTrigger.h"
#include "lardata/AnalysisBase/Calorimetry.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardata/RecoBase/OpFlash.h"
#include "lardata/RecoBase/Shower.h"
#include "lardata/AnalysisBase/T0.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/FindManyP.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Persistency/Common/Ptr.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"

#define MAX_TRACKS 1000
#define MAX_FLASHES 1000
#define MAX_SHOWERS 1000
#define MAX_HITS 200000

using namespace std;

//========================================================================

namespace DUNE{

class NDKAna : public art::EDAnalyzer {
public:

    explicit NDKAna(fhicl::ParameterSet const& pset);
    virtual ~NDKAna();

    void beginJob();
    void endJob();
    void beginRun(const art::Run& run);
    void analyze(const art::Event& evt);

    void reconfigure(fhicl::ParameterSet const& pset);

    void Process(const art::Event& evt, bool &isFiducial);
    void truthMatcher( std::vector<art::Ptr<recob::Hit>> track_hits, const simb::MCParticle *&MCparticle, double &Efrac);
    double truthLength( const simb::MCParticle *MCparticle );
    void calPID( std::vector<const anab::Calorimetry*> cal, double &dEdx, double &range, double &res_range, double &PIDA );
    bool insideFV(double vertex[4]);
    void reset();

private:

    // the parameters we'll read from the .fcl
    std::string fMCTruthModuleLabel;
    std::string fTrackModuleLabel;
    std::string fOpFlashModuleLabel;
    std::string fShowerModuleLabel;
    std::string fHitsModuleLabel;
    double 	fPIDA_endPoint;
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
    int    MC_mother[MAX_TRACKS];  
    double MC_startXYZT[MAX_TRACKS][4]; 
    double MC_endXYZT[MAX_TRACKS][4];  
    double MC_startMomentum[MAX_TRACKS][4]; 
    double MC_endMomentum[MAX_TRACKS][4];  
    double MC_truthlength[MAX_TRACKS];
    std::vector<string> MC_process;  //particle 

    int    n_vertices;
    double vertex[MAX_TRACKS][4];
    int    n_recoTracks;
    int    track_isContained[MAX_TRACKS];
    int    n_vtxFrom_track[MAX_TRACKS];
    double track_vtx[MAX_TRACKS][4];
    double track_end[MAX_TRACKS][4];
    double track_length[MAX_TRACKS]; 
    double track_PIDA[MAX_TRACKS];
    double track_Prange[MAX_TRACKS];
    double track_Efrac[MAX_TRACKS];
    int    track_mcID[MAX_TRACKS];
    int    track_mcPDG[MAX_TRACKS];


    int    n_recoShowers;
    double sh_direction_X[MAX_SHOWERS];
    double sh_direction_Y[MAX_SHOWERS];
    double sh_direction_Z[MAX_SHOWERS];
    double sh_start_X[MAX_SHOWERS];
    double sh_start_Y[MAX_SHOWERS];
    double sh_start_Z[MAX_SHOWERS];
    double sh_energy[MAX_SHOWERS][3];
    double sh_MIPenergy[MAX_SHOWERS][3];
    double sh_dEdx[MAX_SHOWERS][3];
    int    sh_bestplane[MAX_SHOWERS];
    double sh_length[MAX_SHOWERS];

   
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

    int     n_hits;                  
    int     hit_channel[MAX_HITS]; 
    double  hit_peakT[MAX_HITS];  
    double  hit_charge[MAX_HITS];
    double  hit_wireID[MAX_HITS];    
    double  hit_plane[MAX_HITS];    
    double  hit_tpc[MAX_HITS];    
    int     hit_key[MAX_HITS];
    

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

 }; // class NDKAna


//========================================================================
NDKAna::NDKAna(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet)
{
    reconfigure(parameterSet);
}
//========================================================================
NDKAna::~NDKAna(){
  //destructor
}
//========================================================================
void NDKAna::reconfigure(fhicl::ParameterSet const& p){

    fMCTruthModuleLabel  = p.get<std::string>("MCTruthModuleLabel");
    fTrackModuleLabel    = p.get<std::string>("TrackModuleLabel");
    fOpFlashModuleLabel  = p.get<std::string>("OpFlashModuleLabel");
    fShowerModuleLabel	 = p.get<std::string>("ShowerModuleLabel");
    fHitsModuleLabel     = p.get< std::string >("HitsModuleLabel");
    fSaveMCTree		 = p.get<bool>("SaveMCTree");
    fPIDA_endPoint	 = p.get<double>("PIDA_endPoint");
    fFidVolCutX          = p.get<double>("FidVolCutX");
    fFidVolCutY          = p.get<double>("FidVolCutY");
    fFidVolCutZ          = p.get<double>("FidVolCutZ");
}
//========================================================================
void NDKAna::beginJob(){
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
  if( fSaveMCTree ){
    fEventTree = tfs->make<TTree>("Event", "Event Tree from Sim & Reco");
    fEventTree->Branch("eventNo", &Event);
    fEventTree->Branch("runNo", &Run);
    fEventTree->Branch("subRunNo", &SubRun);
    fEventTree->Branch("mc_vertex", &MC_vertex, "mc_vertex[4]/D");
    fEventTree->Branch("mc_npart", &MC_npart);  // number of particles 
    fEventTree->Branch("mc_id", &MC_id, "mc_id[mc_npart]/I");  
    fEventTree->Branch("mc_pdg", &MC_pdg, "mc_pdg[mc_npart]/I"); 
    fEventTree->Branch("mc_mother", &MC_mother, "mc_mother[mc_npart]/I"); 
    fEventTree->Branch("mc_startXYZT", &MC_startXYZT, "mc_startXYZT[mc_npart][4]/D");  
    fEventTree->Branch("mc_endXYZT", &MC_endXYZT, "mc_endXYZT[mc_npart][4]/D"); 
    fEventTree->Branch("mc_startMomentum", &MC_startMomentum, "mc_startMomentum[mc_npart][4]/D");  
    fEventTree->Branch("mc_endMomentum", &MC_endMomentum, "mc_endMomentum[mc_npart][4]/D"); 
    fEventTree->Branch("mc_truthlength", &MC_truthlength, "mc_truthlength[mc_npart]/D"); 
    fEventTree->Branch("mc_process", &MC_process); 

    fEventTree->Branch("n_vertices", &n_vertices);
    fEventTree->Branch("vertex", &vertex,"vertex[n_vertices][4]/D");
    fEventTree->Branch("n_reco_tracks", &n_recoTracks);
    fEventTree->Branch("track_vtx", &track_vtx,"track_vtx[n_reco_tracks][4]/D");
    fEventTree->Branch("track_end", &track_end,"track_end[n_reco_tracks][4]/D");
    fEventTree->Branch("n_vtxFrom_track", &n_vtxFrom_track, "n_vtxFrom_track[n_reco_tracks]/I");
    fEventTree->Branch("track_isContained", &track_isContained,"track_isContained[n_reco_tracks]/I");
    fEventTree->Branch("track_length", &track_length,"track_length[n_reco_tracks]/D");
    fEventTree->Branch("track_PIDA", &track_PIDA,"track_PIDA[n_reco_tracks]/D");
    fEventTree->Branch("track_Prange", &track_Prange,"track_Prange[n_reco_tracks]/D");
    fEventTree->Branch("track_Efrac", &track_Efrac,"track_Efrac[n_reco_tracks]/D");
    fEventTree->Branch("track_mcID", &track_mcID,"track_mcID[n_reco_tracks]/I");
    fEventTree->Branch("track_mcPDG", &track_mcPDG,"track_mcPDG[n_reco_tracks]/I");

    fEventTree->Branch("n_showers", &n_recoShowers);
    fEventTree->Branch("sh_direction_X", &sh_direction_X, "sh_direction_X[n_showers]/D");
    fEventTree->Branch("sh_direction_Y", &sh_direction_Y, "sh_direction_Y[n_showers]/D");
    fEventTree->Branch("sh_direction_Z", &sh_direction_Z, "sh_direction_Z[n_showers]/D");
    fEventTree->Branch("sh_start_X", &sh_start_X, "sh_start_X[n_showers]/D");
    fEventTree->Branch("sh_start_Y", &sh_start_Y, "sh_start_Y[n_showers]/D");
    fEventTree->Branch("sh_start_Z", &sh_start_Z, "sh_start_Z[n_showers]/D");
    fEventTree->Branch("sh_energy", &sh_energy, "sh_energy[n_showers][3]/D");
    fEventTree->Branch("sh_MIPenergy", &sh_MIPenergy, "sh_MIPenergy[n_showers][3]/D");
    fEventTree->Branch("sh_dEdx", &sh_dEdx, "sh_dEdx[n_showers][3]/D");
    fEventTree->Branch("sh_bestplane", &sh_bestplane, "sh_bestplane[n_showers]/I");
    fEventTree->Branch("sh_length", &sh_length, "sh_length[n_showers]/D");

    fEventTree->Branch("n_hits", &n_hits);  //number of hits
    fEventTree->Branch("hit_channel", &hit_channel, "hit_channel[n_hits]/I");  // channel ID
    fEventTree->Branch("hit_peakT", &hit_peakT, "hit_peakT[n_hits]/D");  // peak time
    fEventTree->Branch("hit_charge", &hit_charge, "hit_charge[n_hits]/D");  // charge (area)
    fEventTree->Branch("hit_wireID", &hit_wireID, "hit_wireID[n_hits]/D");  // charge (area)
    fEventTree->Branch("hit_plane", &hit_plane, "hit_plane[n_hits]/D");  // charge (area)
    fEventTree->Branch("hit_tpc", &hit_tpc, "hit_tpc[n_hits]/D");  // charge (area)
    fEventTree->Branch("hit_key", &hit_key, "hit_key[n_hits]/I");  // charge (area)


    fEventTree->Branch("n_flashes", &n_flashes);
    fEventTree->Branch("flash_time", &flash_time,"flash_time[n_flashes]/D");
    fEventTree->Branch("flash_pe", &flash_pe,"flash_pe[n_flashes]/D");
    fEventTree->Branch("flash_ycenter", &flash_ycenter,"flash_ycenter[n_flashes]/D");
    fEventTree->Branch("flash_zcenter", &flash_zcenter,"flash_zcenter[n_flashes]/D");
    fEventTree->Branch("flash_ywidth", &flash_ywidth,"flash_ywidth[n_flashes]/D");
    fEventTree->Branch("flash_zwidth", &flash_zwidth,"flash_zwidth[n_flashes]/D");
    fEventTree->Branch("flash_timewidth", &flash_timewidth, "flash_timewidth[n_flashes]/D");

  }


}
//========================================================================
void NDKAna::endJob(){
     
}
//========================================================================
void NDKAna::beginRun(const art::Run& /*run*/){
  mf::LogInfo("NDKAna")<<"begin run..."<<endl;
}
//========================================================================
void NDKAna::analyze( const art::Event& event ){
    if (event.isRealData()) return;
    reset();

    Event  = event.id().event(); 
    Run    = event.run();
    SubRun = event.subRun();
    bool isFiducial = false;
    Process(event, isFiducial);
    if( fSaveMCTree ){
      if(isFiducial) fEventTree->Fill();
    }
}
//========================================================================
void NDKAna::Process( const art::Event& event, bool &isFiducial){

    art::ServiceHandle<cheat::BackTracker> bt;
    const sim::ParticleList& plist = bt->ParticleList();
    simb::MCParticle *particle=0;
    int i=0; // particle index
    MC_npart = plist.size();
    if( MC_npart > MAX_TRACKS ) return;
    for( sim::ParticleList::const_iterator ipar = plist.begin(); ipar!=plist.end(); ++ipar){
       particle = ipar->second;
       //cout<<particle->PdgCode()<<" "<<particle->Mother()<<" "<<particle->TrackId()<<endl;
       MC_id[i] = particle->TrackId();
       MC_pdg[i] = particle->PdgCode();
       MC_mother[i] = particle->Mother();
       MC_process.push_back(particle->Process());
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
       MC_truthlength[i] = truthLength(particle);
       ++i; //paticle index
    }

    isFiducial =insideFV( MC_vertex );
    if( !isFiducial ) return;
 
    //========================================================================
    //========================================================================
    // Reco  stuff
    //========================================================================
    //========================================================================
    //this is a track base analysis so it must be at least one track
    art::Handle< std::vector<recob::Track> > trackListHandle;
    if(! event.getByLabel(fTrackModuleLabel, trackListHandle)) return;
    std::vector<art::Ptr<recob::Track>> tracklist;
    art::fill_ptr_vector(tracklist, trackListHandle);
    n_recoTracks = tracklist.size();
    if( n_recoTracks > MAX_TRACKS ) return;
    //cout<<"Found these many reco tracks "<<n_recoTracks<<endl;
    
    art::Handle< std::vector<recob::Vertex> > vtxListHandle;
    std::vector<art::Ptr<recob::Vertex>> vtxlist;
    if(!event.getByLabel(fTrackModuleLabel, vtxListHandle)) return;
    art::fill_ptr_vector(vtxlist, vtxListHandle);
    n_vertices = vtxlist.size();

    //cout<<"Found these many vertices "<<n_vertices<<endl; 
     for( int i =0; i<n_vertices; ++i){
       double tmp_vtx[3] ={};
       vtxlist[i]->XYZ(tmp_vtx);
       for( int j=0; j<3; ++j) vertex[i][j]=tmp_vtx[j];
    } 

    //art::FindManyP<recob::Vertex> kink_vertices(trackListHandle, event, art::InputTag("pmtrack","kink"));
    art::FindManyP<recob::Hit> track_hits(trackListHandle, event, fTrackModuleLabel);
    art::FindManyP<recob::Vertex> vtxFromtrack(trackListHandle, event, fTrackModuleLabel);
    //cal dE/dx
    std::string calo_ModuleLabel = fTrackModuleLabel;
    calo_ModuleLabel +="calo";
    art::FindMany<anab::Calorimetry>  reco_cal(trackListHandle, event, calo_ModuleLabel);
    //T0
    art::FindMany<anab::T0>track_T0  (trackListHandle, event, fTrackModuleLabel);
    trkf::TrackMomentumCalculator trackP;
    for(int i=0; i<n_recoTracks; ++i) {
       //std::vector<art::Ptr<recob::Vertex>> kink_vertex_vec = kink_vertices.at(i);
       //cout<<"Kink tracks "<<kink_vertex_vec.size()<<endl;
       //cout<<"This track "<<i<<endl;
       std::vector<const anab::T0*> T0 = track_T0.at(i);
       //cout<<"size T0 "<<T0.size()<<endl; 
       std::vector<art::Ptr<recob::Vertex>> trk_vtxlist = vtxFromtrack.at(i);
       n_vtxFrom_track[i] = trk_vtxlist.size(); 
       art::Ptr<recob::Track> track = tracklist[i];
       track_length[i] = track->Length();
       const TVector3 tmp_track_vtx = track->Vertex();
       const TVector3 tmp_track_end = track->End();
       track_vtx[i][0] =tmp_track_vtx[0];
       track_vtx[i][1] =tmp_track_vtx[1];
       track_vtx[i][2] =tmp_track_vtx[2];
       track_vtx[i][3] = -999.0;

       track_end[i][0] =tmp_track_end[0];
       track_end[i][1] =tmp_track_end[1];
       track_end[i][2] =tmp_track_end[2];
       track_end[i][3] = -999.0;

       track_Prange[i] = trackP.GetTrackMomentum(track_length[i],13); 
       //cout<<"track vertex "<<tmp_track_vtx[0]<<" "<<tmp_track_vtx[1]<<" "<<tmp_track_vtx[2]<<endl;
       double trk_end[4] ={tmp_track_end[0],tmp_track_end[1],tmp_track_end[2],-999};
       bool track_isInside = insideFV( trk_end );
       //check if the track ends within the FV
       if( track_isInside ) track_isContained[i] =1;
       else track_isContained[i] =0;
       //bool track_isFiducial = insideFV( track_vtx );
       //if( !track_isFiducial ) continue;
       std::vector<art::Ptr<recob::Hit>> all_trackHits = track_hits.at(i);  
       //call dEdx
       std::vector<const anab::Calorimetry*> trk_cal = reco_cal.at(i);
       double PID_dEdx, range, res_range, PIDA;
       calPID(trk_cal, PID_dEdx, range, res_range, PIDA); 
       track_PIDA[i] = PIDA;
       double tmpEfrac = 0;
       const simb::MCParticle *particle;
       truthMatcher( all_trackHits, particle, tmpEfrac );
       if(!particle) continue;
       track_mcID[i] = particle->TrackId();
       track_mcPDG[i] = particle->PdgCode();
       track_Efrac[i] = tmpEfrac; 
       //if( track_PIDA[i] == 0 )
       //cout<<track_Prange[i]<<" "<<track_PIDA[i]<<" "<<track_mcPDG[i]<<" "<<PID_dEdx<<" "<<range<<" res "<<res_range<<endl;
       
    }

    /*
    //Save all hits in the event
    art::Handle< std::vector<recob::Hit> > hitListHandle;
    if (! event.getByLabel(fHitsModuleLabel, hitListHandle)) return;
    std::vector<art::Ptr<recob::Hit> > hitlist;
    art::fill_ptr_vector(hitlist, hitListHandle);
    n_hits = hitlist.size();
     
    for(int i=0; i<n_hits; i++) {
       art::Ptr<recob::Hit> hit = hitlist[i];
       hit_channel[i] = hit->Channel();
       hit_charge[i] = hit->Integral();
       hit_peakT[i] = hit->PeakTime();        
       hit_wireID[i] = hit->WireID().Wire;
       hit_tpc[i] = hit->WireID().TPC;
       hit_plane[i] = hit->WireID().Plane;
    }
    */ 

    //Showers... for background rejeciton?
    art::Handle<std::vector<recob::Shower>> showerHandle;
    if(!event.getByLabel(fShowerModuleLabel,showerHandle)) return;
    std::vector<art::Ptr<recob::Shower>> showerlist;
    art::fill_ptr_vector(showerlist, showerHandle); 
    n_recoShowers= showerlist.size();
    //cout<<"Found this many showers "<<n_recoShowers<<endl; 
    for(int i=0; i<n_recoShowers && i< MAX_SHOWERS; ++i){
       art::Ptr<recob::Shower> shower = showerlist[i];
       sh_direction_X[i] = shower->Direction().X();  
       sh_direction_Y[i] = shower->Direction().Y();  
       sh_direction_Z[i] = shower->Direction().Z();  
       sh_start_X[i] = shower->ShowerStart().X();
       sh_start_Y[i] = shower->ShowerStart().Y();
       sh_start_Z[i] = shower->ShowerStart().Z();
       sh_bestplane[i] = shower->best_plane();
       sh_length[i] = shower->Length();
       for( size_t j =0; j<shower->Energy().size(); j ++) sh_energy[i][j] = shower->Energy()[j];
       for( size_t j =0; j<shower->MIPEnergy().size(); j++) sh_MIPenergy[i][j] = shower->MIPEnergy()[j];
       for( size_t j =0; j<shower->dEdx().size(); j++) sh_dEdx[i][j] = shower->dEdx()[j];
    }
    /*
    //PDS info... this may be useful for background rejection
    art::Handle< std::vector<recob::OpFlash> > flashListHandle;
    std::vector<art::Ptr<recob::OpFlash> > flashlist;
    if(event.getByLabel(fOpFlashModuleLabel,flashListHandle))
    art::fill_ptr_vector(flashlist, flashListHandle); 
    
    n_flashes = flashlist.size();
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
    */

}
//========================================================================
void NDKAna::truthMatcher( std::vector<art::Ptr<recob::Hit>> track_hits, const simb::MCParticle *&MCparticle, double &Efrac){

    //cout<<"truthMatcher..."<<endl;
    art::ServiceHandle<cheat::BackTracker> bt;
    std::map<int,double> trkID_E;
    for(size_t j = 0; j < track_hits.size(); ++j){
       art::Ptr<recob::Hit> hit = track_hits[j];
       std::vector<sim::TrackIDE> TrackIDs = bt->HitToTrackID(hit);
       for(size_t k = 0; k < TrackIDs.size(); ++k){
          trkID_E[TrackIDs[k].trackID] += TrackIDs[k].energy;
       }            
    }
    double max_E = -999.0;
    double total_E = 0.0;
    int TrackID = -999;
    double partial_E =0.0; // amount of energy deposited by the particle that deposited more energy... tomato potato... blabla
    //!if the collection of hits have more than one particle associate save the particle w/ the highest energy deposition 
    //!since we are looking for muons/pions/protons this should be enough 
    if( !trkID_E.size() ) {
      MCparticle = 0;
      return; //Ghost track???
    }
    for(std::map<int,double>::iterator ii = trkID_E.begin(); ii!=trkID_E.end(); ++ii){
       total_E += ii->second;
       if((ii->second)>max_E){
         partial_E = ii->second;
         max_E = ii->second;
         TrackID = ii->first;
       }
    } 

    MCparticle = bt->TrackIDToParticle(TrackID);
    Efrac = partial_E/total_E;
    //std::cout<<"total "<<total_E<<" frac "<<Efrac<<" which particle "<<MCparticle->PdgCode()<<std::endl;
}

//========================================================================
void NDKAna::calPID ( std::vector<const anab::Calorimetry*> cal, double &dEdx, double &range, double &res_range, double &PIDA ){
  dEdx  = 0.0;
  range =0.0;
  res_range =0.0;
  PIDA = 0.0;
  int UsedHits = 0;
  //cout<<"size "<<cal.size()<<endl;
  for( int Plane=0; Plane < (int)cal.size(); ++Plane ) { // Loop through planes
     //================================================
     //================================================
     //PIDA
     if( cal[Plane]->dEdx().size() == 0 )
     //cout<<"plane "<<Plane<<" "<<(int)cal[Plane]->dEdx().size()<<endl;
     for( int PlaneHit=0; PlaneHit < (int)cal[Plane]->dEdx().size(); ++PlaneHit ) { // loop through hits on each plane
        //cout<<" calo "<<cal[Plane]->ResidualRange()[PlaneHit]<<" "<<cal[Plane]->dEdx()[PlaneHit]<<endl;
        if( cal[Plane]->ResidualRange()[PlaneHit] < fPIDA_endPoint ) { // Only want PIDA for last x cm
          dEdx += cal[Plane]->dEdx()[PlaneHit]; 
          PIDA += cal[Plane]->dEdx()[PlaneHit]* pow(cal[Plane]->ResidualRange()[PlaneHit], 0.41 ); 
          ++UsedHits;
          range +=cal[Plane]->Range();
          res_range += cal[Plane]->ResidualRange()[PlaneHit];
        } // If ResRange < x cm
     }// Loop over hits on each plane
     //===============================================
     //===============================================
  }// Loop over planes
 
  if ( UsedHits != 0 ){ // If had any hits, work out PIDA and calculate
    PIDA = PIDA /UsedHits;
    dEdx = dEdx / UsedHits;
    range = range/ UsedHits;
    res_range = res_range/UsedHits; 
  }  
} // CalcPIDA

//========================================================================
double NDKAna::truthLength( const simb::MCParticle *MCparticle ){
   //calculate the truth length considering only the part that is inside the TPC
   //Base on a peace of code from dune/TrackingAna/TrackingEfficiency_module.cc

   if( !MCparticle ) return -999.0;
   int numberTrajectoryPoints = MCparticle->NumberTrajectoryPoints();
   double TPCLengthHits[numberTrajectoryPoints];
   int FirstHit=0, LastHit=0;
   double TPCLength = 0.0;
   bool BeenInVolume = false;

   for(int MCHit=0; MCHit < numberTrajectoryPoints; ++MCHit) {
      const TLorentzVector& tmpPosition= MCparticle->Position(MCHit);
      double const tmpPosArray[]={tmpPosition[0],tmpPosition[1],tmpPosition[2]};
      if (MCHit!=0) TPCLengthHits[MCHit] = sqrt( pow( (MCparticle->Vx(MCHit-1)-MCparticle->Vx(MCHit)),2)+ pow( (MCparticle->Vy(MCHit-1)-MCparticle->Vy(MCHit)),2)+ pow( (MCparticle->Vz(MCHit-1)-MCparticle->Vz(MCHit)),2));
      geo::TPCID tpcid = geom->FindTPCAtPosition(tmpPosArray);
      if(tpcid.isValid) {
        // -- Check if hit is within drift window...
        geo::CryostatGeo const& cryo = geom->Cryostat(tpcid.Cryostat);
        geo::TPCGeo      const& tpc  = cryo.TPC(tpcid.TPC);
        double XPlanePosition      = tpc.PlaneLocation(0)[0];
        double DriftTimeCorrection = fabs( tmpPosition[0] - XPlanePosition ) / XDriftVelocity;
        double TimeAtPlane         = MCparticle->T() + DriftTimeCorrection; 
        if( TimeAtPlane < detprop->TriggerOffset() || TimeAtPlane > detprop->TriggerOffset() + WindowSize ) continue;
        LastHit = MCHit;
        if( !BeenInVolume ) {
	  BeenInVolume = true;
          FirstHit = MCHit;
	}
      }		
   }
   for (int Hit = FirstHit+1; Hit <= LastHit; ++Hit ) TPCLength += TPCLengthHits[Hit];
   return TPCLength;
}
//========================================================================
bool NDKAna::insideFV( double vertex[4]){ 

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
void NDKAna::reset(){

   MC_npart =-999; 
   n_recoTracks =-999;
   n_vertices = -999;
   for(int i = 0; i<4; ++i){
      MC_vertex[i] = -999.0;
   }
  
   for(int i=0; i<MAX_TRACKS; ++i) {
       MC_id[i] = -999;
       MC_pdg[i] = -999;
       MC_mother[i] = -999;
       MC_truthlength[i] = -999.0;
       track_PIDA[i] = -999.0;
       track_Prange[i] =-999.0;
       track_Efrac[i] = -999.0; 
       track_mcID[i] = -999;
       track_length[i] = -999;
       track_isContained[i]= -999;
       track_mcPDG[i] =-999;
       for(int j=0; j<4; ++j) {
          MC_startXYZT[i][j]      = -999.0;
          MC_endXYZT[i][j]        = -999.0;
          MC_startMomentum[i][j]  = -999.0;
          MC_endMomentum[i][j]    = -999.0;
	  track_vtx[i][j] 	  = -999.0;
	  track_end[i][j] 	  = -999.0;
       }
    }
    n_recoShowers = -999;
    for(int i=0; i<MAX_SHOWERS; ++i){
      sh_direction_X[i] = -999.0;
      sh_direction_Y[i] = -999.0;
      sh_direction_Z[i] = -999.0;
      sh_start_X[i] = -999.0;
      sh_start_Y[i] = -999.0;
      sh_start_Z[i] = -999.0;
      sh_bestplane[i] = -999.0;
      sh_length[i] = -999.0;
      for( int j=0; j<3; j++){
         sh_energy[i][j] = -999.0;
         sh_MIPenergy[i][j] = -999.0;
         sh_dEdx[i][j] = -999.0;
      }
    } 
    /*
    n_flashes =-999;
    for(int i=0; i<MAX_FLASHES; ++i){
       flash_time[i] =-999.0;
       flash_pe[i] =-999.0;
       flash_ycenter[i] =-999.0;
       flash_zcenter[i] =-999.0;
       flash_ywidth[i] =-999.0;
       flash_zwidth[i] =-999.0;
       flash_timewidth[i] = -999.0;
       flash_abstime[i] =-999;
       flash_frame[i] =-999;
    }
    n_hits =-999;
    for(int i=0; i<MAX_HITS; ++i) {
       hit_channel[i] = 0;
       hit_peakT[i] = 0;
       hit_charge[i] = 0;
       hit_wireID[i] = 0;
       hit_plane[i] = 0;
       hit_tpc[i] =0;
       hit_key[i] =0;
    }
    */
}
//========================================================================
DEFINE_ART_MODULE(NDKAna)

} 

#endif // NDKAna_module
