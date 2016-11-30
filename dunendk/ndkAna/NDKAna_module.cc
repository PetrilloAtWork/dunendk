//Module analyzer
//Ana module for nucleon decay and atmospheric background analysis
//Ana TTTre contains MC truth info and Reconstruction 
//ahiguera@central.uh.edu


#ifndef NDKAna_Module
#define NDKAna_Module

// LArSoft includes
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"  
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larsim/MCCheater/BackTracker.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/AnalysisBase/T0.h"

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

#define MAX_TRACKS 200
#define MAX_FLASHES 100
#define MAX_SHOWERS 100

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
    void truthMatcher( std::vector<art::Ptr<recob::Hit>> all_hits, std::vector<art::Ptr<recob::Hit>> track_hits, const simb::MCParticle *&MCparticle, double &Efrac, double &Ecomplet);
    double myPrange( double track_range );
    double truthLength( const simb::MCParticle *MCparticle );
    void cal( std::vector<const anab::Calorimetry*> cal, double &dEdx, double &range, double &res_range, double &PIDA, double &KE );
    bool insideFV(double vertex[4]);
    void reset();

private:

    // the parameters we'll read from the .fcl
    std::string fMCTruthModuleLabel;
    std::string fTrackModuleLabel;
    std::string fOpFlashModuleLabel;
    std::string fShowerModuleLabel;
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
    double MC_Prange[MAX_TRACKS];
 

    int    n_vertices;
    double vertex[MAX_TRACKS][4];
    int    n_recoTracks;
    int    track_isContained[MAX_TRACKS];
    double track_vtx[MAX_TRACKS][4];
    double track_end[MAX_TRACKS][4];
    double track_length[MAX_TRACKS]; 
    double track_dir_vtx[MAX_TRACKS][4]; 
    double track_PIDA[MAX_TRACKS];
    double track_KE[MAX_TRACKS];
    double track_Prange[MAX_TRACKS];
    double track_myPrange[MAX_TRACKS];
    double track_Efrac[MAX_TRACKS];
    double track_complet[MAX_TRACKS];
    int    track_mcID[MAX_TRACKS];
    int    track_mcPDG[MAX_TRACKS];

    int    n_NDKrecoTracks;
    double track_NDKvtx[MAX_TRACKS][4];
    double track_NDKend[MAX_TRACKS][4];
    int    track_NDK_isContained[MAX_TRACKS];
    double track_NDK_length[MAX_TRACKS]; 
    double track_NDK_dir_vtx[MAX_TRACKS][4]; 
    double track_NDK_PIDA[MAX_TRACKS];
    double track_NDK_KE[MAX_TRACKS];
    double track_NDK_Prange[MAX_TRACKS];
    double track_NDK_myPrange[MAX_TRACKS];
    double track_NDK_Efrac[MAX_TRACKS];
    double track_NDK_complet[MAX_TRACKS];
    int    track_NDK_mcID[MAX_TRACKS];
    int    track_NDK_mcPDG[MAX_TRACKS];




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
  fEventTree->Branch("mc_Prange", &MC_Prange, "mc_Prange[mc_npart]/D"); 
  fEventTree->Branch("mc_truthlength", &MC_truthlength, "mc_truthlength[mc_npart]/D"); 

  fEventTree->Branch("n_vertices", &n_vertices);
  fEventTree->Branch("vertex", &vertex,"vertex[n_vertices][4]/D");
  fEventTree->Branch("n_reco_tracks", &n_recoTracks);
  fEventTree->Branch("track_vtx", &track_vtx,"track_vtx[n_reco_tracks][4]/D");
  fEventTree->Branch("track_end", &track_end,"track_end[n_reco_tracks][4]/D");
  fEventTree->Branch("track_isContained", &track_isContained,"track_isContained[n_reco_tracks]/I");
  fEventTree->Branch("track_length", &track_length,"track_length[n_reco_tracks]/D");
  fEventTree->Branch("track_PIDA", &track_PIDA,"track_PIDA[n_reco_tracks]/D");
  fEventTree->Branch("track_KE", &track_KE,"track_KE[n_reco_tracks]/D");
  fEventTree->Branch("track_Prange", &track_Prange,"track_Prange[n_reco_tracks]/D");
  fEventTree->Branch("track_myPrange", &track_myPrange,"track_myPrange[n_reco_tracks]/D");
  fEventTree->Branch("track_complet", &track_complet,"track_complet[n_reco_tracks]/D");
  fEventTree->Branch("track_Efrac", &track_Efrac,"track_Efrac[n_reco_tracks]/D");
  fEventTree->Branch("track_mcID", &track_mcID,"track_mcID[n_reco_tracks]/I");
  fEventTree->Branch("track_mcPDG", &track_mcPDG,"track_mcPDG[n_reco_tracks]/I");

  fEventTree->Branch("n_NDKreco_tracks", &n_NDKrecoTracks);
  fEventTree->Branch("track_NDKvtx", &track_NDKvtx,"track_NDKvtx[n_NDKreco_tracks][4]/D");
  fEventTree->Branch("track_NDKend", &track_NDKend,"track_NDKend[n_NDKreco_tracks][4]/D");
  fEventTree->Branch("track_NDK_isContained", &track_NDK_isContained,"track_NDK_isContained[n_NDKreco_tracks]/I");
  fEventTree->Branch("track_NDK_length", &track_NDK_length,"track_NDK_length[n_NDKreco_tracks]/D");
  fEventTree->Branch("track_NDK_PIDA", &track_NDK_PIDA,"track_NDK_PIDA[n_NDKreco_tracks]/D");
  fEventTree->Branch("track_NDK_KE", &track_NDK_KE,"track_NDK_KE[n_NDKreco_tracks]/D");
  fEventTree->Branch("track_NDK_Prange", &track_NDK_Prange,"track_NDK_Prange[n_NDKreco_tracks]/D");
  fEventTree->Branch("track_NDK_myPrange", &track_NDK_myPrange,"track_NDK_myPrange[n_NDKreco_tracks]/D");
  fEventTree->Branch("track_NDK_Efrac", &track_NDK_Efrac,"track_NDK_Efrac[n_NDKreco_tracks]/D");
  fEventTree->Branch("track_NDK_complet", &track_NDK_complet,"track_NDK_complet[n_NDKreco_tracks]/D");
  fEventTree->Branch("track_NDK_mcID", &track_NDK_mcID,"track_NDK_mcID[n_NDKreco_tracks]/I");
  fEventTree->Branch("track_NDK_mcPDG", &track_NDK_mcPDG,"track_NDK_mcPDG[n_NDKreco_tracks]/I");

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

/*
  fEventTree->Branch("n_flashes", &n_flashes);
  fEventTree->Branch("flash_time", &flash_time,"flash_time[n_flashes]/D");
  fEventTree->Branch("flash_pe", &flash_pe,"flash_pe[n_flashes]/D");
  fEventTree->Branch("flash_ycenter", &flash_ycenter,"flash_ycenter[n_flashes]/D");
  fEventTree->Branch("flash_zcenter", &flash_zcenter,"flash_zcenter[n_flashes]/D");
  fEventTree->Branch("flash_ywidth", &flash_ywidth,"flash_ywidth[n_flashes]/D");
  fEventTree->Branch("flash_zwidth", &flash_zwidth,"flash_zwidth[n_flashes]/D");
  fEventTree->Branch("flash_timewidth", &flash_timewidth, "flash_timewidth[n_flashes]/D");
  */ 
  

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
    if(isFiducial) fEventTree->Fill();
    
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
       MC_id[i] = particle->TrackId();
       MC_pdg[i] = particle->PdgCode();
       MC_mother[i] = particle->Mother();
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
       MC_Prange[i] = myPrange(MC_truthlength[i]);
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
    std::vector<art::Ptr<recob::Track>> tracklist;
    if( event.getByLabel(fTrackModuleLabel, trackListHandle)); 
      art::fill_ptr_vector(tracklist, trackListHandle);

    n_recoTracks = tracklist.size();
    if( n_recoTracks > MAX_TRACKS || n_recoTracks == 0) return;
    
    art::Handle< std::vector<recob::Vertex> > vtxListHandle;
    std::vector<art::Ptr<recob::Vertex>> vtxlist;
    if(!event.getByLabel(fTrackModuleLabel, vtxListHandle)) return; 
    art::fill_ptr_vector(vtxlist, vtxListHandle);
    n_vertices = vtxlist.size();

    for( int i =0; i<n_vertices; ++i){
       double tmp_vtx[3] ={-999.0,-999.0,-999.0};
       vtxlist[i]->XYZ(tmp_vtx);
       for( int j=0; j<3; ++j) vertex[i][j]=tmp_vtx[j];
    } 
    art::FindManyP<recob::Hit> track_hits(trackListHandle, event, fTrackModuleLabel);
    std::string calo_ModuleLabel = fTrackModuleLabel;
    calo_ModuleLabel +="calo";
    art::FindMany<anab::Calorimetry>  reco_cal(trackListHandle, event, calo_ModuleLabel);
    trkf::TrackMomentumCalculator trackP;

    std::vector<art::Ptr<recob::Hit>> tmp_all_trackHits = track_hits.at(0);  
    std::vector<art::Ptr<recob::Hit>> all_hits;
    art::Handle<std::vector<recob::Hit>> hithandle;
    if(event.get(tmp_all_trackHits[0].id(), hithandle))  art::fill_ptr_vector(all_hits, hithandle);

    for(int i=0; i<n_recoTracks; ++i) {
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
       track_myPrange[i] = myPrange( track_length[i]);
       double trk_end[4] ={tmp_track_end[0],tmp_track_end[1],tmp_track_end[2],-999};
       bool track_isInside = insideFV( trk_end );
       //check if the track ends within the FV
       if( track_isInside ) track_isContained[i] =1;
       else track_isContained[i] =0;
       std::vector<art::Ptr<recob::Hit>> all_trackHits = track_hits.at(i);  
       //call dEdx
       std::vector<const anab::Calorimetry*> trk_cal = reco_cal.at(i);
       double PID_dEdx, range, res_range;// PIDA;
       double PIDA =-999.0;
       double KE =-999.0;
       cal(trk_cal, PID_dEdx, range, res_range, PIDA, KE); 
       track_PIDA[i] = PIDA;
       track_KE[i] = KE;
       double tmpEfrac = 0;
       const simb::MCParticle *particle;
       double tmpComplet = 0;
       truthMatcher( all_hits,  all_trackHits, particle, tmpEfrac, tmpComplet );
       if(!particle) continue;
       track_mcID[i] = particle->TrackId();
       track_mcPDG[i] = particle->PdgCode();
       track_Efrac[i] = tmpEfrac; 
       track_complet[i] = tmpComplet;
    }

    
    //========================================================================
    // NDK Reco  stuff
    //========================================================================
    //This part is for new track reconstruction in addition to the standard one
    //It is still under development, see https://indico.fnal.gov/getFile.py/access?contribId=4&resId=0&materialId=slides&confId=13197 for more details 
    art::Handle< std::vector<recob::Track> > trackListNDKHandle;
    std::vector<art::Ptr<recob::Track>> trackNDKlist;
    if(! event.getByLabel("pmtrackNDK", trackListNDKHandle)) return;
    art::fill_ptr_vector(trackNDKlist, trackListNDKHandle);
    n_NDKrecoTracks  = trackNDKlist.size();

    if( n_NDKrecoTracks > MAX_TRACKS || n_NDKrecoTracks == 0) return;

    art::FindManyP<recob::Hit> trackNDK_hits(trackListNDKHandle, event,"pmtrackNDK" );
    art::FindMany<anab::Calorimetry>  reco_NDKcal(trackListNDKHandle, event, "pmtrackNDKcalo");

    std::vector<art::Ptr<recob::Hit>> tmp_all_NDKtrackHits = trackNDK_hits.at(0);  
    std::vector<art::Ptr<recob::Hit>> all_NDKhits;
    art::Handle<std::vector<recob::Hit>> NDKhithandle;
    if(event.get(tmp_all_NDKtrackHits[0].id(), NDKhithandle))  art::fill_ptr_vector(all_NDKhits, NDKhithandle);


    for(int i=0; i<n_NDKrecoTracks; ++i) {
       art::Ptr<recob::Track> track = trackNDKlist[i];
       track_NDK_length[i] = track->Length();
       const TVector3 tmp_track_vtx = track->Vertex();
       const TVector3 tmp_track_end = track->End();
       track_NDKvtx[i][0] =tmp_track_vtx[0];
       track_NDKvtx[i][1] =tmp_track_vtx[1];
       track_NDKvtx[i][2] =tmp_track_vtx[2];
       track_NDKvtx[i][3] = -999.0;

       track_NDKend[i][0] =tmp_track_end[0];
       track_NDKend[i][1] =tmp_track_end[1];
       track_NDKend[i][2] =tmp_track_end[2];
       track_NDKend[i][3] = -999.0;

       track_NDK_Prange[i] = trackP.GetTrackMomentum(track_NDK_length[i],13); 
       track_NDK_myPrange[i] = myPrange( track_NDK_length[i]);
       double trk_end[4] ={tmp_track_end[0],tmp_track_end[1],tmp_track_end[2],-999};
       bool track_isInside = insideFV( trk_end );
       //check if the track ends within the FV
       if( track_isInside ) track_NDK_isContained[i] =1;
       else track_NDK_isContained[i] =0;
       std::vector<art::Ptr<recob::Hit>> all_trackHits = trackNDK_hits.at(i);  
       //call dEdx
       std::vector<const anab::Calorimetry*> trk_cal = reco_NDKcal.at(i);
       double PID_dEdx, range, res_range;// PIDA;
       double PIDA =-999.0;
       double KE =-999.0;
       cal(trk_cal, PID_dEdx, range, res_range, PIDA, KE); 
       track_NDK_PIDA[i] = PIDA;
       track_NDK_KE[i] = KE;
       double tmpEfrac = 0.0;
       double tmpComplet =0.0;
       const simb::MCParticle *particle;
       truthMatcher( all_NDKhits,  all_trackHits, particle, tmpEfrac, tmpComplet );
       if(!particle) continue;
       track_NDK_mcID[i] = particle->TrackId();
       track_NDK_mcPDG[i] = particle->PdgCode();
       track_NDK_Efrac[i] = tmpEfrac; 
       track_NDK_complet[i] = tmpComplet;
 
    }

    
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
void NDKAna::truthMatcher( std::vector<art::Ptr<recob::Hit>> all_hits, std::vector<art::Ptr<recob::Hit>> track_hits, const simb::MCParticle *&MCparticle, double &Efrac, double &Ecomplet){

    //std::cout<<"truthMatcher..."<<std::endl;
    art::ServiceHandle<cheat::BackTracker> bt;
    std::map<int,double> trkID_E;
    for(size_t j = 0; j < track_hits.size(); ++j){
       art::Ptr<recob::Hit> hit = track_hits[j];
       std::vector<sim::TrackIDE> TrackIDs = bt->HitToTrackID(hit);
       for(size_t k = 0; k < TrackIDs.size(); k++){
          trkID_E[TrackIDs[k].trackID] += TrackIDs[k].energy;
       }            
    }
    double E_em =0.0;
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
         if( TrackID < 0 ) E_em += ii->second;
       }
    } 
    MCparticle = bt->TrackIDToParticle(TrackID);

    //In the current simulation, we do not save EM Shower daughters in GEANT. But we do save the energy deposition in TrackIDEs. If the energy deposition is from a particle that is the daughter of 
    //an EM particle, the negative of the parent track ID is saved in TrackIDE for the daughter particle
    //we don't want to track gammas or any other EM activity 
    if( TrackID < 0 ) return;

    //Efrac = (partial_E+E_em)/total_E;
    Efrac = (partial_E)/total_E;

    //completeness
    double totenergy =0;
    for(size_t k = 0; k < all_hits.size(); ++k){
       art::Ptr<recob::Hit> hit = all_hits[k];
       std::vector<sim::TrackIDE> TrackIDs = bt->HitToTrackID(hit);
       for(size_t l = 0; l < TrackIDs.size(); ++l){
          if(TrackIDs[l].trackID==TrackID) totenergy += TrackIDs[l].energy;
       }
    } 
    Ecomplet = partial_E/totenergy;
}

//========================================================================
void NDKAna::cal ( std::vector<const anab::Calorimetry*> cal, double &dEdx, double &range, double &res_range, double &PIDA, double &KE ){
  dEdx  = 0.0;
  range =0.0;
  res_range =0.0;
  PIDA = 0.0;

  int UsedHits = 0;
  //cout<<"size "<<cal.size()<<endl;
  int best_plane =-1;
  int plane0 = (int)cal[0]->dEdx().size(); 
  int plane1 = (int)cal[1]->dEdx().size(); 
  int plane2 = (int)cal[2]->dEdx().size(); 

  if((plane0 >= plane1) && (plane0 >= plane2)) {
    best_plane = 0;
  }
  else if ((plane1 >= plane0) && (plane1 >= plane2)) {
    best_plane = 1;
  }
  else {
    best_plane = 2;
  }
  //cout<<"best "<<best_plane<<endl;
  KE = cal[best_plane]->KineticEnergy();
  for( int PlaneHit=0; PlaneHit < (int)cal[best_plane]->dEdx().size(); ++PlaneHit ) { // loop through hits on each plane
     if( cal[best_plane]->ResidualRange()[PlaneHit] < fPIDA_endPoint ) { // Only want PIDA for last x cm
       dEdx += cal[best_plane]->dEdx()[PlaneHit]; 
       PIDA += cal[best_plane]->dEdx()[PlaneHit]* pow(cal[best_plane]->ResidualRange()[PlaneHit], 0.42 ); 
       range +=cal[best_plane]->Range();
       res_range += cal[best_plane]->ResidualRange()[PlaneHit];
       ++UsedHits;
          //cout<<"used hits "<<UsedHits<<endl;
     } // If ResRange < x cm
  }// Loop over hits on each plane
  
  if ( UsedHits != 0 ){ // If had any hits, work out PIDA and calculate
    PIDA =PIDA /UsedHits;
    dEdx = dEdx / UsedHits;
    range = range/ UsedHits;
    res_range = res_range/UsedHits; 
    //cout<<PIDA<<endl;
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
double NDKAna::myPrange( double track_range ){

  double KE =0.0; 
  double Prange = 0.0;
  double M = 105.658;
  if(track_range>0 && track_range<=200)
    KE = 10.658+ (3.71181*track_range) + (-0.0302898*track_range*track_range) +
         (0.000228501*track_range*track_range*track_range)+
         (-5.86201E-7*track_range*track_range*track_range*track_range);
  else if (track_range>200 && track_range<=1.91E4)                         
    KE = 30.6157+ (2.11089*track_range) + (0.000255267*track_range*track_range)+
         (-5.19254E-8*track_range*track_range*track_range)+
         (6.39454E-12*track_range*track_range*track_range*track_range)+
         (-3.99392E-16*track_range*track_range*track_range*track_range*track_range)+
         (9.73174E-21*track_range*track_range*track_range*track_range*track_range*track_range); 

  return  Prange = sqrt((KE*KE)+(2*M*KE));
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

   MC_npart =0; 
   n_recoTracks =0;
   n_NDKrecoTracks =0;
   n_vertices = 0;
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
       track_myPrange[i] =-999.0;
       track_Efrac[i] = -999.0; 
       track_complet[i] =-999.0;
       track_mcID[i] = -999;
       track_length[i] = -999;
       track_isContained[i]= -999;
       track_mcPDG[i] =-999;
       track_KE[i] =-999.0;

       track_NDK_PIDA[i] = -999.0;
       track_NDK_Prange[i] =-999.0;
       track_NDK_myPrange[i] =-999.0;
       track_NDK_Efrac[i] = -999.0; 
       track_NDK_complet[i] = -999.0;
       track_NDK_mcID[i] = -999;
       track_NDK_length[i] = -999;
       track_NDK_isContained[i]= -999;
       track_NDK_mcPDG[i] =-999;
       track_NDK_KE[i] =-999.0;
       for(int j=0; j<4; ++j) {
          MC_startXYZT[i][j]      = -999.0;
          MC_endXYZT[i][j]        = -999.0;
          MC_startMomentum[i][j]  = -999.0;
          MC_endMomentum[i][j]    = -999.0;
	  track_vtx[i][j] 	  = -999.0;
	  track_end[i][j] 	  = -999.0;
	  track_NDKvtx[i][j] 	  = -999.0;
	  track_NDKend[i][j] 	  = -999.0;
      
       }
    }
    n_recoShowers = 0;
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
    n_flashes =0;
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
    */
}
//========================================================================
DEFINE_ART_MODULE(NDKAna)

} 

#endif // NDKAna_module
