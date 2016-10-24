


// LArSoft includes
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardataobj/RecoBase/Track.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larreco/RecoAlg/DBScanAlg.h"
#include "larreco/ClusterFinder/ClusterCreator.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "larreco/RecoAlg/ClusterParamsImportWrapper.h"
#include "larreco/RecoAlg/ClusterRecoUtil/StandardClusterParamsAlg.h"
#include "larreco/RecoAlg/MergeClusterAlg.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
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


using namespace std;


namespace cluster{
//========================================================================

class NDKReco : public art::EDProducer {
public:

    explicit NDKReco(fhicl::ParameterSet const& pset);
    virtual ~NDKReco();

    void produce(art::Event& evt);
    void reconfigure(fhicl::ParameterSet const& pset);


private:

    // the parameters we'll read from the .fcl
    std::string fTrackModuleLabel;
    
    DBScanAlg fDBScan;
}; 
}
namespace cluster{
//========================================================================
NDKReco::NDKReco(fhicl::ParameterSet const& parameterSet)
    : fDBScan(parameterSet.get< fhicl::ParameterSet >("DBScanAlg")) 
{
    reconfigure(parameterSet);
    produces<std::vector<recob::Cluster>>(); 
    produces<art::Assns<recob::Cluster,recob::Hit>>();
}
//========================================================================
NDKReco::~NDKReco(){
}
//========================================================================
void NDKReco::reconfigure(fhicl::ParameterSet const& p){

    fTrackModuleLabel    = p.get<std::string>("TrackModuleLabel");
    fDBScan.reconfigure(p.get< fhicl::ParameterSet >("DBScanAlg"));
}
//========================================================================
void NDKReco::produce(art::Event& event ){

    
    art::Handle< std::vector<recob::Track> > trackListHandle;
    if(! event.getByLabel(fTrackModuleLabel, trackListHandle)) return;
    std::vector<art::Ptr<recob::Track>> tracklist;
    art::fill_ptr_vector(tracklist, trackListHandle);

    //get a collection of clusters   
    std::unique_ptr<std::vector<recob::Cluster> > ccol(new std::vector<recob::Cluster>);
    std::unique_ptr< art::Assns<recob::Cluster, recob::Hit> > assn(new art::Assns<recob::Cluster, recob::Hit>);
    //Compute the cluster characteristics
    ClusterParamsImportWrapper<StandardClusterParamsAlg> ClusterParamAlgo;


    art::FindManyP<recob::Hit> track_hits(trackListHandle, event, fTrackModuleLabel);
    art::Ptr<recob::Track> mu_track_candidate;     
    std::vector<art::Ptr<recob::Hit>> mu_trackHits;  
    double tmp_length =-999.0;
    //the longest track should be the muon track
    for(size_t i=0; i<tracklist.size(); ++i) {
       art::Ptr<recob::Track> track = tracklist[i];
       if( track->Length() > tmp_length){
         tmp_length = track->Length();
         mu_track_candidate = tracklist[i];
         mu_trackHits = track_hits.at(i);  
       }
    }
    cout<<"Muon Track "<<mu_track_candidate->Length()<<" w/Hits "<<mu_trackHits.size()<<endl;
    std::vector<double> startTrk;
    std::vector<double> endTrk;
    mu_track_candidate->Extent(startTrk,endTrk);

    art::FindManyP<recob::SpacePoint> space_pts(mu_trackHits, event, fTrackModuleLabel);

    std::vector<art::Ptr<recob::Hit>> hits_after_cleaning;
    for(size_t i=0; i<mu_trackHits.size(); ++i){
       std::vector<art::Ptr<recob::SpacePoint>> hit_spts = space_pts.at(i);
       double dis_start = sqrt(pow(hit_spts[0]->XYZ()[0]-startTrk[0],2)+ pow(hit_spts[0]->XYZ()[1]-startTrk[1],2)+pow(hit_spts[0]->XYZ()[2]-startTrk[2],2));
       double dis_end = sqrt(pow(hit_spts[0]->XYZ()[0]-endTrk[0],2)+ pow(hit_spts[0]->XYZ()[1]-endTrk[1],2)+pow(hit_spts[0]->XYZ()[2]-endTrk[2],2));
       art::Ptr<recob::Hit> hit = mu_trackHits[i];
       if( dis_start < 10.0 ){
         if( hit->Integral() > 250.0 ) hits_after_cleaning.push_back(hit);
       }
       else if( dis_end < 10.0){
         if( hit->Integral() > 250.0 ) hits_after_cleaning.push_back(hit);
       }
    }
  
    cout<<"After track cleaning there are these many hits "<<hits_after_cleaning.size()<<endl; 

   // check if there are any bad channels to consider  
   lariov::ChannelStatusProvider const& channelStatus = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
   lariov::ChannelStatusProvider::ChannelSet_t const BadChannels = channelStatus.BadChannels();
   //create clusters
   fDBScan.InitScan(hits_after_cleaning, BadChannels);
   fDBScan.run_cluster();

   
   for(size_t i = 0; i < fDBScan.fclusters.size(); ++i){
      art::PtrVector<recob::Hit> clusterHits;
      for(size_t j = 0; j < fDBScan.fpointId_to_clusterId.size(); ++j){	  
         if(fDBScan.fpointId_to_clusterId[j]==i){
           clusterHits.push_back(hits_after_cleaning[j]);
         }
      }
      if(clusterHits.size()>0){
        const geo::WireID& wireID = clusterHits.front()->WireID();
        unsigned int sw = wireID.Wire;
        unsigned int ew = clusterHits.back()->WireID().Wire;
        // Put cluster hits in the algorithm
        ClusterParamAlgo.ImportHits(clusterHits); 
        
        ClusterCreator cluster(
          ClusterParamAlgo,				// algo
          float(sw),					// start_wire
          0.,						// sigma_start_wire
          clusterHits.front()->PeakTime(),		// start_tick
          clusterHits.front()->SigmaPeakTime(),		// sigma_start_tick
          float(ew),					// end_wire
          0.,						// sigma_end_wire,
          clusterHits.back()->PeakTime(),		// end_tick
          clusterHits.back()->SigmaPeakTime(),		// sigma_end_tick
          ccol->size(),					// ID
          clusterHits.front()->View(),			// view
          clusterHits.front()->WireID().planeID(),	// plane
          recob::Cluster::Sentry			// sentry
        );
        ccol->emplace_back(cluster.move());
        
        
        // associate the hits to this cluster
        util::CreateAssn(*this, event, *(ccol.get()), clusterHits, *(assn.get()));
        clusterHits.clear();
      }       
   }
   
   event.put(std::move(ccol));
   event.put(std::move(assn));
   return;

} 

//========================================================================
}
namespace cluster{
  DEFINE_ART_MODULE(NDKReco)

}//cluster 

