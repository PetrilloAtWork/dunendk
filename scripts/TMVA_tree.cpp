#include <iostream>
#include "include/NDKAna.h"
#include "include/include.h"
using namespace std;

int TMVA( const char *filename, const char *outfile ){

  TFile *f_temple = new TFile("../PIDA-loglikelihood/bkgd.root");
  TH2D *h_dEdx_p =  (TH2D*)f_temple->Get("h_dEdx"); 
  TH2D *h_dEdx_rev_p =  (TH2D*)f_temple->Get("h_dEdx_rev"); 

  TFile *f_kaon = new TFile("../PIDA-loglikelihood/signal.root");
  TH2D *h_dEdx_k =  (TH2D*)f_kaon->Get("h_dEdx"); 
  TH2D *h_dEdx_rev_k =  (TH2D*)f_kaon->Get("h_dEdx_rev"); 


  TFile *f_data = new TFile(filename,"READ");

  TFile *f = new TFile(outfile,"RECREATE");

  cout<<"Lets dance..."<<std::endl;

  TTree *tmva_tree = new TTree("TMVA_tree","MVA tree");

  float n_tracks;
  float P_range;
  float PIDA_longest;
  float PIDA_2longest;
  float PIDA_shortest;
  float rev_PIDA;
  float pdg_longest;
  float pdg_shortest;
  float trk_longest;
  float trk_shortest;
  float n_showers;
  float etotal_shower;
  float Em_e;
  float Ev_ccqe;
  float trk_e;
  float n_vtx;
  float KE;
  float frac_trkE;
  float trk_sh_KE;
  float PID_loglike;
  float cnn_score;
  float wgt;
  float mc_K_KE;
  int ev_n;
  int sub_run;
  int run;

  tmva_tree->Branch("ev_n",&ev_n);
  tmva_tree->Branch("sub_run",&sub_run);
  tmva_tree->Branch("run",&run); 
  tmva_tree->Branch("mc_K_KE",&mc_K_KE); 
  tmva_tree->Branch("n_tracks",&n_tracks); 
  tmva_tree->Branch("n_showers",&n_showers); 
  tmva_tree->Branch("n_vtx",&n_vtx); 
  tmva_tree->Branch("trk_e",&trk_e); 
  tmva_tree->Branch("EM_e",&Em_e); 
  tmva_tree->Branch("PIDA_longest",&PIDA_longest); 
  tmva_tree->Branch("PIDA_shortest",&PIDA_shortest); 
  tmva_tree->Branch("PID_loglike",&PID_loglike); 
  tmva_tree->Branch("trk_longest",&trk_longest); 
  tmva_tree->Branch("trk_shortest",&trk_shortest); 
  tmva_tree->Branch("P_range",&P_range); 
  tmva_tree->Branch("KE",&KE); 
  tmva_tree->Branch("frac_trkE",&frac_trkE); 
  tmva_tree->Branch("cnn_score",&cnn_score); 

  //========================================================
  float saved =0;
  float n_ine=0;
  float no_ine=0;
  double kaon_mass = 0.4936;
  //========================================================
  TTree *tree =(TTree*)f_data->Get("NDKAna/Event");
  NDKAna *signal = new NDKAna(tree);
  float n_entries = tree->GetEntries();
  cout<<"..this many entries "<<n_entries<<endl;
  for( size_t i=0; i<n_entries; ++i){
     
     signal->GetEntry(i);
     //======= MC truth
     bool inelastic = false;
     bool noninteraction = false;
     for( size_t it =0; it<signal->mcgenie_npart; ++it){
       if( signal->mcgenie_pdg[it] == 321 ){
         mc_K_KE = signal->mcgenie_startMomentum[it][3]-kaon_mass;
         if( signal->mcgenie_fate[it] ==4 || signal->mcgenie_fate[it] == 5) 
           inelastic = true;
         else if( signal->mcgenie_fate[it] == 1 )
           noninteraction = true;
       }
     }
     if(inelastic) n_ine ++;
     else no_ine ++;

     //======= Reco 
     if( signal->n_reco_tracks < 2 ) continue;

     run = signal->runNo;
     sub_run = signal->subRunNo;
     ev_n = signal->eventNo;

     n_showers = signal->n_showers;
     Em_e= signal->Em_e*1000.0;  //EM vis E
     trk_e= signal->trk_e*1000.0;//trk-like vis E

     KE = Em_e + trk_e;  //sum of all visible energy
     n_vtx = signal->n_vertices;
 
     n_tracks = signal->n_reco_tracks;
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
 
     double like_p=1;
     double like_k=1;
     double loglike =1; 
     double like_p_rev=1;
     double like_k_rev=1;
     double loglike_rev =1; 

     for( int k=0; k<signal->n_cal_points[shortest_trk_idx]; k++){
          if( k > 49 ) break; //only 50 bins in the template
          double dEdx =signal->track_dE_dx[shortest_trk_idx][k]; 
          int bin_p =h_dEdx_p->GetYaxis()->FindBin(dEdx);
          int bin_k =h_dEdx_k->GetYaxis()->FindBin(dEdx);
          double tmp_like_p =h_dEdx_p->GetBinContent(k+1,bin_p)/h_dEdx_p->ProjectionY("",k+1,k+1)->Integral(1,h_dEdx_p->GetYaxis()->GetNbins()); 
          double tmp_like_k =h_dEdx_k->GetBinContent(k+1,bin_k)/h_dEdx_k->ProjectionY("",k+1,k+1)->Integral(1,h_dEdx_k->GetYaxis()->GetNbins()); 
              
          if( tmp_like_p == 0 ) tmp_like_p = 1e-3;
          if( tmp_like_k == 0 ) tmp_like_k = 1e-3;
          like_p *= tmp_like_p;
          like_k *= tmp_like_k;

          //using reverse temple
          int rev_idx = signal->n_cal_points[shortest_trk_idx]-1-k;
          int bin_p_rev =h_dEdx_rev_p->GetYaxis()->FindBin(dEdx);
          int bin_k_rev =h_dEdx_rev_k->GetYaxis()->FindBin(dEdx);
          double tmp_like_p_rev =h_dEdx_rev_p->GetBinContent(rev_idx+1,bin_p_rev)/h_dEdx_rev_p->ProjectionY("",rev_idx+1,rev_idx+1)->Integral(1,h_dEdx_rev_p->GetYaxis()->GetNbins()); 
          double tmp_like_k_rev =h_dEdx_rev_k->GetBinContent(rev_idx+1,bin_k_rev)/h_dEdx_rev_k->ProjectionY("",rev_idx+1,rev_idx+1)->Integral(1,h_dEdx_rev_k->GetYaxis()->GetNbins()); 
              
          if( tmp_like_p_rev == 0 ) tmp_like_p_rev = 1e-3;
          if( tmp_like_k_rev == 0 ) tmp_like_k_rev = 1e-3;
          like_p_rev *= tmp_like_p_rev;
          like_k_rev *= tmp_like_k_rev;
       
     }
     double like_ratio = like_p/like_k;
     double like_ratio_rev = like_p_rev/like_k_rev;
     loglike = TMath::Log(like_ratio);
     loglike_rev = TMath::Log(like_ratio_rev);
     double d_like_ratio = (like_p+like_p_rev)/(like_k+like_k_rev);
     PID_loglike = TMath::Log(d_like_ratio);
 
     PIDA_longest = signal->track_PIDA[longest_trk_idx][signal->track_bestplane[longest_trk_idx]];
     if( PIDA_longest <0 ) PIDA_longest=signal->track_PIDA[longest_trk_idx][0];
     if( PIDA_longest <0 ) PIDA_longest=signal->track_PIDA[longest_trk_idx][1];
     if( PIDA_longest <0 ) PIDA_longest=signal->track_PIDA[longest_trk_idx][2];
     if( PIDA_longest <0 ) PIDA_longest=0.0;
     if( isnan(PIDA_longest) ) PIDA_longest=0.0;
     
     PIDA_shortest = signal->track_PIDA[shortest_trk_idx][signal->track_bestplane[shortest_trk_idx]];
     if( PIDA_shortest <0 ) PIDA_shortest=signal->track_PIDA[shortest_trk_idx][0];
     if( PIDA_shortest <0 ) PIDA_shortest=signal->track_PIDA[shortest_trk_idx][1];
     if( PIDA_shortest <0 ) PIDA_shortest=signal->track_PIDA[shortest_trk_idx][2];
     if( PIDA_shortest <0 ) PIDA_shortest=0.0;
     if( isnan(PIDA_shortest) ) PIDA_shortest=0.0;

     P_range = signal->track_Prange[longest_trk_idx];
     P_range *= 1000.0;
     trk_longest = signal->track_length[longest_trk_idx];
     trk_shortest = signal->track_length[shortest_trk_idx];

     frac_trkE = trk_e/KE;
     if( isnan(frac_trkE) ) frac_trkE=0.0;
     if( isnan(PID_loglike) || isinf(PID_loglike) ) PID_loglike=0.0;
     if( trk_longest > 100) continue;

     saved ++;
     tmva_tree->Fill();
  }/// All entries end
  cout<<"All "<<n_entries<<endl; 
  cout<<"Passed "<<saved<<" eff "<<double(saved)/n_entries<<endl;
  cout<<"non ine "<<no_ine<<endl;
  cout<<"All ine "<<n_ine<<endl;
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
    cout<<"Enter a ROOT file to process..."<<endl;
    return 0;
  }
  else if( argc == 2) return TMVA(argv[1],"TMVA_tree.root");
  else if( argc == 3) return TMVA(argv[1],argv[2]);

  else return 0;
}



