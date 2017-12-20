//Module analyzer
//Ana module for nucleon decay and atmospheric analysis
//Ana TTree contains MC truth and reconstruction info
//ahiguera@central.uh.edu


#ifndef NDKAna_Module
#define NDKAna_Module

// LArSoft includes
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Cluster.h"  
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larsim/MCCheater/BackTracker.h"
#include "larsim/MCCheater/PhotonBackTracker.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardata/ArtDataHelper/MVAReader.h"
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
#define MAX_FLASHES 20000
#define MAX_SHOWERS 1000
#define MVA_LENGTH 4
#define MAX_HITS 50000
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
    void PIDAcal( std::vector<const anab::Calorimetry*> cal, std::vector<double> &PIDA );
    void Process(const art::Event& evt, bool &isFiducial);
    void truthMatcher( std::vector<art::Ptr<recob::Hit>> all_hits, std::vector<art::Ptr<recob::Hit>> track_hits, const simb::MCParticle *&MCparticle, double &Efrac, double &Ecomplet);
    double truthLength( const simb::MCParticle *MCparticle );
    bool insideFV(double vertex[4]);
    void reset();

private:

    // the parameters we'll read from the .fcl
    std::string fHitModuleLabel;
    std::string fTrackModuleLabel;
    std::string fOpFlashModuleLabel;
    std::string fShowerModuleLabel;
    std::string fPointIDModuleLabel;
    std::string fNNetModuleLabel;
    std::string fPFParticleModuleLabel;
    double      fExponentConstant;
    double 	fPIDA_endPoint;
    double	fMaxPIDAValue;
    double	fMinPIDAValue;
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
    int    n_decayVtx;
    double decayVtx[MAX_TRACKS][3];
    int    n_recoTracks;
    int    track_isContained[MAX_TRACKS];
    double track_vtxDir[MAX_TRACKS][3];
    double track_vtx[MAX_TRACKS][4];
    double track_end[MAX_TRACKS][4];
    double track_length[MAX_TRACKS]; 
    double track_dir_vtx[MAX_TRACKS][4]; 
    double track_PIDA_bp[MAX_TRACKS];
    double track_PIDA[MAX_TRACKS][3];
    int	   track_PID_pdg[MAX_TRACKS][3];
    int	   track_PID_pdg_bp[MAX_TRACKS];
    double track_KE[MAX_TRACKS][3];
    double track_KE_bp[MAX_TRACKS];
    double track_Prange[MAX_TRACKS];
    double track_Efrac[MAX_TRACKS];
    double track_complet[MAX_TRACKS];
    int    track_mcID[MAX_TRACKS];
    int    track_mcPDG[MAX_TRACKS];
    int    n_cal_points;
    double track_dQ_dx[10000];
    double track_dE_dx[10000];
    double track_range[10000];
   
    int     n_recoHits;
    int     hit_channel[MAX_HITS];
    int     hit_tpc[MAX_HITS];
    int     hit_plane[MAX_HITS];
    int     hit_wire[MAX_HITS];
    double  hit_peakT[MAX_HITS];
    double  hit_charge[MAX_HITS];
    double  hit_ph[MAX_HITS];
    double  hit_startT[MAX_HITS];
    double  hit_endT[MAX_HITS];
    double  hit_rms[MAX_HITS];
    double  hit_electrons[MAX_HITS];
    
  
    double Em_ch;
    double Em_e;
    double trk_e;
    double Emichel_e;

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
    double flash_PE_ndk[MAX_FLASHES];
    double flash_PE_Ar39[MAX_FLASHES];
 
    double fFidVolCutX;
    double fFidVolCutY;
    double fFidVolCutZ;

    double fFidVolXmin;
    double fFidVolXmax;
    double fFidVolYmin;
    double fFidVolYmax;
    double fFidVolZmin;
    double fFidVolZmax;

    double fPidValue;
    unsigned int    fView;

    bool fSaveHits;

    detinfo::DetectorProperties const *detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    detinfo::DetectorClocks const *ts = lar::providerFrom<detinfo::DetectorClocksService>();
    double XDriftVelocity = detprop->DriftVelocity()*1e-3; //cm/ns
    double WindowSize     = detprop->NumberTimeSamples() * ts->TPCClock().TickPeriod() * 1e3;
    art::ServiceHandle<geo::Geometry> geom;
    calo::CalorimetryAlg fCalorimetryAlg;
 }; // class NDKAna


//========================================================================

NDKAna::NDKAna(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet), fCalorimetryAlg(parameterSet.get< fhicl::ParameterSet >("CalorimetryAlg"))
    
{
    reconfigure(parameterSet);
}
//========================================================================
NDKAna::~NDKAna(){
  //destructor
}
//========================================================================
void NDKAna::reconfigure(fhicl::ParameterSet const& p)
{

    fHitModuleLabel      = p.get<std::string>("HitModuleLabel");
    fTrackModuleLabel    = p.get<std::string>("TrackModuleLabel");
    fOpFlashModuleLabel  = p.get<std::string>("OpFlashModuleLabel");
    fShowerModuleLabel	 = p.get<std::string>("ShowerModuleLabel");
    fPointIDModuleLabel  = p.get<std::string>("PointIDModuleLabel");
    fNNetModuleLabel     = p.get<std::string>("NNetModuleLabel");
    fSaveHits		 = p.get<bool>("SaveHits");
    fExponentConstant 	 = p.get<double>("ExponentConstant");
    fMaxPIDAValue 	 = p.get<double>("MaxPIDAValue");
    fMinPIDAValue 	 = p.get<double>("MinPIDAValue");
    fPIDA_endPoint	 = p.get<double>("PIDA_endPoint");
    fView                = p.get<double>("View");
    fPidValue		 = p.get<double>("PidValue");
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
  fEventTree->Branch("n_decayVtx", &n_decayVtx);
  fEventTree->Branch("decayVtx", &decayVtx,"decayVtx[n_decayVtx][3]/D");
  fEventTree->Branch("track_vtx", &track_vtx,"track_vtx[n_reco_tracks][4]/D");
  fEventTree->Branch("track_vtxDir", &track_vtxDir,"track_vtxDir[n_reco_tracks][3]/D");
  fEventTree->Branch("track_end", &track_end,"track_end[n_reco_tracks][4]/D");
  fEventTree->Branch("track_isContained", &track_isContained,"track_isContained[n_reco_tracks]/I");
  fEventTree->Branch("track_length", &track_length,"track_length[n_reco_tracks]/D");
  fEventTree->Branch("track_PIDA_bp", &track_PIDA_bp,"track_PIDA_bp[n_reco_tracks]/D");
  fEventTree->Branch("track_PIDA", &track_PIDA,"track_PIDA[n_reco_tracks][3]/D");
  fEventTree->Branch("track_PID_pdg", &track_PID_pdg,"track_PID_pdg[n_reco_tracks][3]/I");
  fEventTree->Branch("track_PID_pdg_bp", &track_PID_pdg_bp,"track_PID_pdg_bp[n_reco_tracks]/I");
  fEventTree->Branch("track_KE", &track_KE,"track_KE[n_reco_tracks][3]/D");
  fEventTree->Branch("track_KE_bp", &track_KE_bp,"track_KE_bp[n_reco_tracks]/D");
  fEventTree->Branch("track_Prange", &track_Prange,"track_Prange[n_reco_tracks]/D");
  fEventTree->Branch("n_cal_points", &n_cal_points);
  fEventTree->Branch("track_dQ_dx", &track_dQ_dx,"track_dQ_dx[n_cal_points]/D");
  fEventTree->Branch("track_dE_dx", &track_dE_dx,"track_dE_dx[n_cal_points]/D");
  fEventTree->Branch("track_range", &track_range,"track_range[n_cal_points]/D");
  fEventTree->Branch("track_complet", &track_complet,"track_complet[n_reco_tracks]/D");
  fEventTree->Branch("track_Efrac", &track_Efrac,"track_Efrac[n_reco_tracks]/D");
  fEventTree->Branch("track_mcID", &track_mcID,"track_mcID[n_reco_tracks]/I");
  fEventTree->Branch("track_mcPDG", &track_mcPDG,"track_mcPDG[n_reco_tracks]/I");

  fEventTree->Branch("Em_ch", &Em_ch);
  fEventTree->Branch("Em_e", &Em_e);
  fEventTree->Branch("trk_e", &trk_e);
  fEventTree->Branch("Emichel_e", &Emichel_e);
  if (fSaveHits){ 
  fEventTree->Branch("n_recoHits", &n_recoHits);
  fEventTree->Branch("hit_channel", &hit_channel, "hit_channel[n_recoHits]/I");
  fEventTree->Branch("hit_tpc", &hit_tpc, "hit_tpc[n_recoHits]/I");
  fEventTree->Branch("hit_plane", &hit_plane, "hit_plane[n_recoHits]/I");
  fEventTree->Branch("hit_wire", &hit_wire, "hit_wire[n_recoHits]/I");
  fEventTree->Branch("hit_peakT", &hit_peakT, "hit_peakT[n_recoHits]/D");
  fEventTree->Branch("hit_charge", &hit_charge, "hit_charge[n_recoHits]/D");
  fEventTree->Branch("hit_ph", &hit_ph, "hit_ph[n_recoHits]/D");
  fEventTree->Branch("hit_endT", &hit_endT, "hit_endT[n_recoHits]/D");
  fEventTree->Branch("hit_rms", &hit_rms, "hit_rms[n_recoHits]/D");
  fEventTree->Branch("hit_electrons", &hit_electrons, "hit_electrons[n_recoHits]/D");
  }

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

  fEventTree->Branch("flash_PE_ndk", &flash_PE_ndk, "flash_PE_ndk[n_flashes]/D");
  fEventTree->Branch("flash_PE_Ar39", &flash_PE_Ar39, "flash_PE_Ar39[n_flashes]/D");

 

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
    //reset();

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
       //MC_Prange[i] = myPrange(MC_truthlength[i]);
      
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
    if( event.getByLabel(fTrackModuleLabel, trackListHandle)) 
      art::fill_ptr_vector(tracklist, trackListHandle);
    n_recoTracks = tracklist.size();
    if( n_recoTracks > MAX_TRACKS || n_recoTracks == 0) return;
    
    art::Handle< std::vector<recob::Vertex> > vtxListHandle;
    std::vector<art::Ptr<recob::Vertex>> vtxlist;
    if(event.getByLabel("pmtrack", vtxListHandle)) 
      art::fill_ptr_vector(vtxlist, vtxListHandle);

    n_vertices = vtxlist.size();
    if( n_vertices != 0 )
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

    std::string PID_ModuleLabel = fTrackModuleLabel;
    PID_ModuleLabel +="pid";
    art::FindMany<anab::ParticleID> reco_PID(trackListHandle, event, PID_ModuleLabel); 

    art::Handle<std::vector<recob::Hit>> HitHandle;
    std::vector<art::Ptr<recob::Hit>> all_hits;
    if(event.getByLabel(fHitModuleLabel,HitHandle)) 
      art::fill_ptr_vector(all_hits, HitHandle); 
    
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
   
       const TVector3 tmp_vtx_dir = track->VertexDirection();
       track_vtxDir[i][0] = tmp_vtx_dir[0]; 
       track_vtxDir[i][1] = tmp_vtx_dir[1]; 
       track_vtxDir[i][2] = tmp_vtx_dir[2]; 

       track_Prange[i] = trackP.GetTrackMomentum(track_length[i],13); 
       double trk_end[4] ={tmp_track_end[0],tmp_track_end[1],tmp_track_end[2],-999};
       bool track_isInside = insideFV( trk_end );
       //check if the track ends within the FV
       if( track_isInside ) track_isContained[i] =1;
       else track_isContained[i] =0;
       //calculate PID 
       std::vector<const anab::Calorimetry*> trk_cal = reco_cal.at(i);
       std::vector<const anab::ParticleID*> trk_pid = reco_PID.at(i);

       int plane0 =   trk_pid[0]->Ndf();
       int plane1 =   trk_pid[1]->Ndf();
       int plane2 =   trk_pid[2]->Ndf();
       int best_plane =-1;
       int most_ndf = std::max({plane0, plane1, plane2});
       for( size_t p =0; p<3; p++) if( most_ndf == trk_pid[p]->Ndf() ) best_plane = p;

       track_PID_pdg[i][0] = trk_pid[0]->Pdg();
       track_PID_pdg[i][1] = trk_pid[1]->Pdg();
       track_PID_pdg[i][2] = trk_pid[2]->Pdg();
       track_KE[i][0] = trk_cal[0]->KineticEnergy();
       track_KE[i][1] = trk_cal[1]->KineticEnergy();
       track_KE[i][2] = trk_cal[2]->KineticEnergy();

       std::vector<double> PIDAval;
       std::vector<double> chi2;
       PIDAcal( trk_cal, PIDAval); 
       int idx =0;
       for(auto const& val : PIDAval){
          track_PIDA[i][idx] = val;
          idx ++;
       }
       track_PIDA_bp[i] = track_PIDA[i][best_plane]; 
       track_PID_pdg_bp[i] = trk_pid[best_plane]->Pdg();
       track_KE_bp[i] = trk_cal[best_plane]->KineticEnergy();

      
       //save dE/dx & dQ/dx
       n_cal_points = trk_cal[best_plane]->dEdx().size();
       for( unsigned i =0; i<trk_cal[best_plane]->dEdx().size(); ++i ) {
          track_dQ_dx[i] = trk_cal[best_plane]->dQdx()[i];
          track_dE_dx[i] = trk_cal[best_plane]->dEdx()[i];
          track_range[i] = trk_cal[best_plane]->ResidualRange()[i];
       }
       //truth matcher
       double tmpEfrac = 0;
       const simb::MCParticle *particle;
       double tmpComplet = 0;
       std::vector<art::Ptr<recob::Hit>> all_trackHits = track_hits.at(i);  
       truthMatcher( all_hits,  all_trackHits, particle, tmpEfrac, tmpComplet );
       if(!particle) continue;
       track_mcID[i] = particle->TrackId();
       track_mcPDG[i] = particle->PdgCode();
       track_Efrac[i] = tmpEfrac; 
       track_complet[i] = tmpComplet;
    }
    //CNN dacayID vertex 
    art::Handle<std::vector<recob::Vertex>> dcy_vtxHandle;
    std::vector<art::Ptr<recob::Vertex>> dcy_vtxlist;
    if(event.getByLabel(fPointIDModuleLabel,dcy_vtxHandle)) 
      art::fill_ptr_vector(dcy_vtxlist, dcy_vtxHandle); 
    n_decayVtx= dcy_vtxlist.size();
   
    //art::FindManyP<recob::Track> decay_tracklist(dcy_vtxHandle, event, fPointIDModuleLabel);
    if( n_decayVtx !=0 )
    for( int i=0; i< n_decayVtx; ++i){
       double tmp_vtx[3] ={-999.0,-999.0,-999.0};
       dcy_vtxlist[i]->XYZ(tmp_vtx);
       for( int j=0; j<3; ++j) decayVtx[i][j]=tmp_vtx[j]; 
       //std::vector<art::Ptr<recob::Track>>  decay_track = decay_tracklist.at(i);
       //cout<<"how many tracks? "<<decay_track.size()<<endl;
    }
    //CNN Em-trk hits
    //Imported code from PointIdEffTest_module.cc 
    //to included hit and CNN output info into the analysis module
    //to quantify how much EM activity we have in order to reduce background 
    anab::MVAReader<recob::Hit, MVA_LENGTH> hitResults(event, fNNetModuleLabel);                     // hit-by-hit outpus just to be dumped to file for debugging
    auto cluResults = anab::MVAReader<recob::Cluster, MVA_LENGTH>::create(event, fNNetModuleLabel);  // outputs for clusters recovered in not-throwing way 
    int trk_idx = hitResults.getIndex("track");
    int em_idx = hitResults.getIndex("em");
    //int michel_idx = hitResults.getIndex("michel");
    Em_ch =0.0;
    Em_e =0.0;
    trk_e =0.0;
    if(cluResults){
      const art::FindManyP<recob::Hit> hitsFromClusters(cluResults->dataHandle(), event, cluResults->dataTag());
      for(size_t c = 0; c < cluResults->size(); ++c){
	 const recob::Cluster & clu = cluResults->item(c);
         if(clu.Plane().Plane != fView) continue;
    	 const std::vector< art::Ptr<recob::Hit> > & hits = hitsFromClusters.at(c);
         std::vector< anab::FeatureVector<MVA_LENGTH> >  hit_outs = hitResults.outputs();
         //EM hits
         for(auto const & h : hits){
            auto const & vout = hit_outs[h.key()];
            double p_trk_or_sh = vout[trk_idx] + vout[em_idx];
	    if(p_trk_or_sh > 0){ 
              double PidValue = vout[trk_idx] / p_trk_or_sh;
              if( PidValue < fPidValue ){
                Em_ch += h->SummedADC()* fCalorimetryAlg.LifetimeCorrection( h->PeakTime() );
                Em_e  += fCalorimetryAlg.ElectronsFromADCArea( h->Integral(), h->WireID().Plane) * fCalorimetryAlg.LifetimeCorrection( h->PeakTime() );
                //Michel hits
                //if( vout[michel_idx] > 0.1 ) //temporay
                //Emichel_e  += fCalorimetryAlg.ElectronsFromADCArea( h->Integral(), h->WireID().Plane) * fCalorimetryAlg.LifetimeCorrection( h->PeakTime() );
              }
              else  trk_e += fCalorimetryAlg.ElectronsFromADCArea( h->Integral(), h->WireID().Plane) * fCalorimetryAlg.LifetimeCorrection( h->PeakTime() ); 
	   }           
         }
      }
    }
    //Showers... for background rejection
    art::Handle<std::vector<recob::Shower>> showerHandle;
    std::vector<art::Ptr<recob::Shower>> showerlist;
    if(event.getByLabel(fShowerModuleLabel,showerHandle)) 
      art::fill_ptr_vector(showerlist, showerHandle); 
    n_recoShowers= showerlist.size();
    if( n_recoShowers != 0 )
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

    if( fSaveHits ){
      //Hits
      n_recoHits= all_hits.size();
      if( n_recoHits != 0 )
      for(int i = 0; i < n_recoHits && i < MAX_HITS ; ++i){//loop over hits
         hit_channel[i] = all_hits[i]->Channel();
         hit_tpc[i]   = all_hits[i]->WireID().TPC;
         hit_plane[i]   = all_hits[i]->WireID().Plane;
         hit_wire[i]    = all_hits[i]->WireID().Wire;
         hit_peakT[i]   = all_hits[i]->PeakTime();
         hit_charge[i]  = all_hits[i]->Integral();
         hit_ph[i]  = all_hits[i]->PeakAmplitude();
         hit_startT[i] = all_hits[i]->PeakTimeMinusRMS();
         hit_endT[i] = all_hits[i]->PeakTimePlusRMS();
         hit_rms[i] = all_hits[i]->RMS();
         hit_electrons[i]  = fCalorimetryAlg.ElectronsFromADCArea( all_hits[i]->Integral(), all_hits[i]->WireID().Plane) * fCalorimetryAlg.LifetimeCorrection( all_hits[i]->PeakTime() );
      }
    }

    /*
    //PDS info... this may be useful for background rejection

    art::ServiceHandle<cheat::PhotonBackTracker> Photon_bt;
    art::Handle< std::vector<recob::OpFlash> > flashListHandle;
    std::vector<art::Ptr<recob::OpFlash> > flashlist;
    if(event.getByLabel(fOpFlashModuleLabel,flashListHandle))
      art::fill_ptr_vector(flashlist, flashListHandle); 
    //cout<<"flashes "<<flashlist.size()<<endl; 
    n_flashes = flashlist.size();
    if( n_flashes != 0 ){
      art::FindManyP<recob::OpHit> flash_ophits(flashListHandle, event, fOpFlashModuleLabel);

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

         std::vector<art::Ptr<recob::OpHit>> thisflash_Hits = flash_ophits.at(i);

         double Pe_Ar39 =0.0;
         double Pe_ndk =0.0;
         for( size_t j =0; j<thisflash_Hits.size(); ++j){
            std::vector<sim::TrackSDP> trkSDPs = Photon_bt->OpHitToEveSDPs(thisflash_Hits[j]);
            if( trkSDPs.size() != 0 ){ 
              for(const sim::TrackSDP trkSDP : trkSDPs){
                 if( trkSDP.trackID != 0){//Why ??
                   double E_frac = trkSDP.energyFrac;  //fraction of OpHit energy from the particle with this trackID
                   const simb::MCParticle* particle = Photon_bt->TrackIDToParticle(trkSDP.trackID); 
                   //cout<<"e frac "<<E_frac<<" "<<particle->PdgCode()<<" "<<thisflash_Hits[j]->PE()<<endl;
                   if( particle->PdgCode() == 11 &&  particle->Mother() == 0 ){ 
                     Pe_Ar39 += thisflash_Hits[j]->PE()*E_frac;
                   }
                   else{
                     Pe_ndk +=  thisflash_Hits[j]->PE()*E_frac; 
                   }
                 }
              }
            }
         }
       
         flash_PE_ndk[i]= Pe_ndk;
         flash_PE_Ar39[i]= Pe_Ar39;
      }
    }
     */     
    //cout<<"pass all "<<isFiducial<<endl; 
}
//========================================================================
void NDKAna::PIDAcal( std::vector<const anab::Calorimetry*> cal, std::vector<double> &PIDA){
  std::vector<double> pida_vec_v0;
  std::vector<double> pida_vec_v1;
  std::vector<double> pida_vec_v2;

  int used_points[3] ={0,0,0};
  double tmp_pida[3];
  for( unsigned j =0; j<3; ++j){
     for( unsigned i =0; i<cal[j]->dEdx().size(); ++i ) { // loop through hits on each plane
        if( cal[j]->ResidualRange()[i] < fPIDA_endPoint ) { // Only want PIDA for last x cm
          tmp_pida[j] = cal[j]->dEdx()[i]* pow(cal[j]->ResidualRange()[i], fExponentConstant ); 
          if(fMinPIDAValue > tmp_pida[j] || tmp_pida[j] > fMaxPIDAValue) continue;
          if( j ==  0 )pida_vec_v0.push_back(tmp_pida[j]);
          if( j ==  1 )pida_vec_v1.push_back(tmp_pida[j]);
          if( j ==  2 )pida_vec_v2.push_back(tmp_pida[j]);
          used_points[j] ++;
        } // If ResRange < x cm
     }// Loop over hits on each plane
  }

  
  //for each view calculate PIDA median value
  std::sort(pida_vec_v0.begin(), pida_vec_v0.end() );
  int size_v0 = pida_vec_v0.size();
  double median_v0 = -999;
  if( size_v0 > 0 ) median_v0 = size_v0 % 2 ? pida_vec_v0[size_v0 / 2] : (pida_vec_v0[size_v0 / 2 - 1] + pida_vec_v0[size_v0 / 2]) / 2;

  std::sort(pida_vec_v1.begin(), pida_vec_v1.end() );
  int size_v1 = pida_vec_v1.size();
  double median_v1 = -999;
  if( size_v1 > 0 ) median_v1 = size_v1 % 2 ? pida_vec_v1[size_v1 / 2] : (pida_vec_v1[size_v1 / 2 - 1] + pida_vec_v1[size_v1 / 2]) / 2;

  std::sort(pida_vec_v2.begin(), pida_vec_v2.end() );
  int size_v2 = pida_vec_v2.size();
  double median_v2 = -999;
  if( size_v2 > 0 ) median_v2 = size_v2 % 2 ? pida_vec_v2[size_v2 / 2] : (pida_vec_v2[size_v2 / 2 - 1] + pida_vec_v2[size_v2 / 2]) / 2;

  PIDA.push_back(median_v0);
  PIDA.push_back(median_v1);
  PIDA.push_back(median_v2);

}
//========================================================================
void NDKAna::truthMatcher( std::vector<art::Ptr<recob::Hit>> all_hits, std::vector<art::Ptr<recob::Hit>> track_hits, const simb::MCParticle *&MCparticle, double &Efrac, double &Ecomplet){

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
double NDKAna::truthLength( const simb::MCParticle *MCparticle ){
   //calculate the truth length considering only the part that is inside the TPC
   //Base on a peace of code from dune/TrackingAna/TrackingEfficiency_module.cc

   if( !MCparticle ) return -999.0;
   double TPCLength = 0.0;
   int numberTrajectoryPoints = MCparticle->NumberTrajectoryPoints();
   std::vector<double> TPCLengthHits(numberTrajectoryPoints, 0);
   int FirstHit=0, LastHit=0;
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
   n_vertices = 0;
   n_decayVtx =0;
   for(int i = 0; i<4; ++i){
      MC_vertex[i] = -999.0;
   }
  
   for(int i=0; i<MAX_TRACKS; ++i) {
       MC_id[i] = -999;
       MC_pdg[i] = -999;
       MC_mother[i] = -999;
       MC_truthlength[i] = -999.0;
       //track_PIDA[i] = -999.0;
       track_Prange[i] =-999.0;
       track_Efrac[i] = -999.0; 
       track_complet[i] =-999.0;
       track_mcID[i] = -999;
       track_length[i] = -999;
       track_isContained[i]= -999;
       track_mcPDG[i] =-999;
       //track_KE[i] =-999.0;
       //track_PID_pdg[i] = -999;
       for(int j=0; j<4; ++j) {
          MC_startXYZT[i][j]      = -999.0;
          MC_endXYZT[i][j]        = -999.0;
          MC_startMomentum[i][j]  = -999.0;
          MC_endMomentum[i][j]    = -999.0;
	  track_vtx[i][j] 	  = -999.0;
	  track_end[i][j] 	  = -999.0;
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
       flash_abstime[i] =-999.0;
       flash_frame[i] =-999;
       flash_PE_ndk[i] = -999;  
       flash_PE_Ar39[i] = -999;  
    }
    */ 
}
//========================================================================
DEFINE_ART_MODULE(NDKAna)

} 

#endif // NDKAna_module
