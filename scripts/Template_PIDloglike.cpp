//===================================================
// Template for log-likelihood PID
// Based on C. Marshall's idea 
// see: https://indico.fnal.gov/event/18097/contribution/0/material/slides/0.pdf
// Skol...
//===================================================
#include <iostream>
#include <algorithm>
#include "include/NDKAna.h"
#include "include/include.h"
using namespace std;
int Template_PIDloglike(  const char *filename, const char *outfile){

  TFile *f_data = new TFile(filename,"READ");

  TFile *f = new TFile(outfile,"RECREATE");
  TH2D *h_dEdx = new TH2D("h_dEdx",";cal point;dE/dx",50,0,50,50,0,50); 
  TH2D *h_dEdx_rev = new TH2D("h_dEdx_rev",";cal point;dE/dx",50,0,50,50,0,50); 
 
  //========================================================
  TTree *tree =(TTree*)f_data->Get("NDKAna/Event");
  NDKAna *signal = new NDKAna(tree);
  int n_entries = tree->GetEntries();
  cout<<"..this many entries "<<n_entries<<endl;
  

  for( size_t i=0; i<n_entries; ++i){
     signal->GetEntry(i);

     //==========================================================================================
     // RECO INFO 
     //==========================================================================================

     if( signal->n_reco_tracks < 2 ) continue;

     //================================================
     // longest and shortest tracks
     //================================================
     float length_tmp = -999.0;
     int longest_trk_idx =-1; 
     float short_tmp = 999999;
     int shortest_trk_idx =-1;
     int longest_vtx_ID =-999;
     for( size_t ii=0; ii<signal->n_reco_tracks; ++ii){
        if( signal->track_length[ii] > length_tmp ){
          length_tmp = signal->track_length[ii];
          longest_trk_idx = ii;
        } 
        if( signal->track_length[ii] < short_tmp){
           short_tmp = signal->track_length[ii]; 
           shortest_trk_idx = ii;
        }
     }
     if( shortest_trk_idx  == longest_trk_idx ) continue; //this should never happened though 

     //longest track must be a muon
     //short track cannot be a pion or lepton 
     if( abs(signal->track_mcPDG[longest_trk_idx]) != 13 || signal->track_mcPDG[shortest_trk_idx] < 321 ) continue;  

     TVector3 sh_trk_vtx(signal->track_vtx[shortest_trk_idx][0],signal->track_vtx[shortest_trk_idx][1],signal->track_vtx[shortest_trk_idx][2]);
     TVector3 sh_trk_end(signal->track_end[shortest_trk_idx][0],signal->track_end[shortest_trk_idx][1],signal->track_end[shortest_trk_idx][2]);
     TVector3 lo_trk_vtx(signal->track_vtx[longest_trk_idx][0],signal->track_vtx[longest_trk_idx][1],signal->track_vtx[longest_trk_idx][2]);
     TVector3 lo_trk_end(signal->track_end[longest_trk_idx][0],signal->track_end[longest_trk_idx][1],signal->track_end[longest_trk_idx][2]);

     double min_start = std::min((sh_trk_vtx-lo_trk_vtx).Mag(), (sh_trk_end-lo_trk_vtx).Mag() );    
     double min_end = std::min((sh_trk_vtx-lo_trk_end).Mag(), (sh_trk_end-lo_trk_end).Mag() );    

     int sh_trk_dir =0; 
     int lo_trk_dir =0; 
     //check muon and short track directions & common vertex?
     // dir = 1 correct 
     // dir =-1 reverse dir
     if( min_start < min_end ){ 
       lo_trk_dir =1;  
       if( (sh_trk_vtx-lo_trk_vtx).Mag() < (sh_trk_end - lo_trk_vtx).Mag() ) sh_trk_dir = 1;
       else sh_trk_dir = -1;
     } 
     else{  
       lo_trk_dir =-1;
       if( (sh_trk_vtx-lo_trk_end).Mag() < (sh_trk_end - lo_trk_end).Mag() ) sh_trk_dir = 1;
       else sh_trk_dir = -1;
     }

     //follow Marshall's convention pg 9,10
     for( int j=0; j<signal->n_cal_points[shortest_trk_idx]; ++j){
        if( sh_trk_dir == 1 ){
          h_dEdx->Fill(signal->n_cal_points[shortest_trk_idx]-1-j,signal->track_dE_dx[shortest_trk_idx][j]); 
          h_dEdx_rev->Fill(j,signal->track_dE_dx[shortest_trk_idx][j]); 
        }
        else if( sh_trk_dir == -1){
          h_dEdx_rev->Fill(signal->n_cal_points[shortest_trk_idx]-1-j ,signal->track_dE_dx[shortest_trk_idx][j]); 
          h_dEdx->Fill(j,signal->track_dE_dx[shortest_trk_idx][j]); 
          
        }  
     } 
 
  }/// All entries end
  f->Write();
  f->Close();
  return 0;

}
//======================================================================
int main( int argc, char *argv[] ){
  cout<<"**************************** "<<endl;
  cout<<"*    WELCOME TO JAMAICA    * "<<endl;
  cout<<"*      HAVE A NICE DAY     * "<<endl;
  cout<<"**************************** "<<endl;

  if( argc == 1){
    cout<<"Enter a ROOT file to process and output.ROOT file name ..."<<endl;
    return 0;
  }
  if( argc == 3) return Template_PIDloglike(argv[1],argv[2]);

  else return 0;
}

