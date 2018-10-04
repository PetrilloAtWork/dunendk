//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Apr 23 17:37:10 2018 by ROOT version 6.13/01
// from TTree Event/Event Tree from Sim & Reco
// found on file: ../Data_reco_v2/ndk_test.root
//////////////////////////////////////////////////////////

#ifndef NDKAna_h
#define NDKAna_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
using namespace std;
class NDKAna {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           eventNo;
   Int_t           runNo;
   Int_t           subRunNo;
   Double_t        MC_Ev;
   Int_t           MC_cc;
   Double_t        MC_Q2;
   Int_t           MC_nuPDG;
   Double_t        MC_hit_nucleon;
   Int_t           mcgenie_npart;
   Int_t           mcgenie_id[110];   //[mcgenie_npart]
   Int_t           mcgenie_fate[110];   //[mcgenie_npart]
   Int_t           mcgenie_statusCode[110];   //[mcgenie_npart]
   Int_t           mcgenie_pdg[110];   //[mcgenie_npart]
   Int_t           mcgenie_mother[110];   //[mcgenie_npart]
   Double_t        mcgenie_startMomentum[110][4];   //[mcgenie_npart]
   Double_t        mcgenie_endMomentum[110][4];   //[mcgenie_npart]
   Double_t        mc_vertex[4];
   Int_t           mc_npart;
   Int_t           mc_id[2770];   //[mc_npart]
   Int_t           mc_pdg[2770];   //[mc_npart]
   Int_t           mc_statusCode[2770];   //[mc_npart]
   Int_t           mc_mother[2770];   //[mc_npart]
   Double_t        mc_startXYZT[2770][4];   //[mc_npart]
   Double_t        mc_endXYZT[2770][4];   //[mc_npart]
   Double_t        mc_startMomentum[2770][4];   //[mc_npart]
   Double_t        mc_endMomentum[2770][4];   //[mc_npart]
   Double_t        mc_Prange[2770];   //[mc_npart]
   Double_t        mc_truthlength[2770];   //[mc_npart]
   vector<string>  *mc_process;
   vector<string>  *mc_Endprocess;
   Int_t           n_vertices;
   Double_t        vertex[60][4];   //[n_vertices]
   Int_t           n_reco_tracks;
   Int_t           n_decayVtx;
   Double_t        decayVtx[10][3];   //[n_decayVtx]
   Double_t        track_vtx[60][4];   //[n_reco_tracks]
   Double_t        track_vtxDir[60][3];   //[n_reco_tracks]
   Double_t        track_end[60][4];   //[n_reco_tracks]
   Int_t           track_isContained[60];   //[n_reco_tracks]
   Double_t        track_length[60];   //[n_reco_tracks]
   Double_t        track_PIDA[60][3];   //[n_reco_tracks]
   Int_t           track_PID_pdg[60][3];   //[n_reco_tracks]
   Double_t        track_KE[60][3];   //[n_reco_tracks]
   Double_t        track_Prange[60];   //[n_reco_tracks]
   Int_t           track_bestplane[60];   //[n_reco_tracks]
   Int_t           n_cal_points[60];
   Double_t        track_dQ_dx[60][500];   //[n_cal_points]
   Double_t        track_dE_dx[60][500];   //[n_cal_points]
   Double_t        track_range[60][500];   //[n_cal_points]
   Double_t        track_pitch[60][500];   //[n_cal_points]
   Double_t        track_complet[60];   //[n_reco_tracks]
   Double_t        track_Efrac[60];   //[n_reco_tracks]
   Int_t           track_mcID[60];   //[n_reco_tracks]
   Int_t           track_mcPDG[60];   //[n_reco_tracks]
   Double_t        Em_ch;
   Double_t        Em_e;
   Double_t        trk_e;
   Double_t        Emichel_e;
   Int_t           n_showers;
   Double_t        sh_direction_X[80];   //[n_showers]
   Double_t        sh_direction_Y[80];   //[n_showers]
   Double_t        sh_direction_Z[80];   //[n_showers]
   Double_t        sh_start_X[80];   //[n_showers]
   Double_t        sh_start_Y[80];   //[n_showers]
   Double_t        sh_start_Z[80];   //[n_showers]
   Double_t        sh_energy[80][3];   //[n_showers]
   Double_t        sh_MIPenergy[80][3];   //[n_showers]
   Double_t        sh_dEdx[80][3];   //[n_showers]
   Int_t           sh_bestplane[80];   //[n_showers]
   Double_t        sh_length[80];   //[n_showers]
   Int_t           n_flashes;
   Double_t        flash_time[1];   //[n_flashes]
   Double_t        flash_pe[1];   //[n_flashes]
   Double_t        flash_ycenter[1];   //[n_flashes]
   Double_t        flash_zcenter[1];   //[n_flashes]
   Double_t        flash_ywidth[1];   //[n_flashes]
   Double_t        flash_zwidth[1];   //[n_flashes]
   Double_t        flash_timewidth[1];   //[n_flashes]
   Double_t        flash_abstime[1];   //[n_flashes]
   Int_t           flash_frame[1];   //[n_flashes]
   Double_t        flash_PE_ndk[1];   //[n_flashes]
   Double_t        flash_PE_Ar39[1];   //[n_flashes]

   // List of branches
   TBranch        *b_eventNo;   //!
   TBranch        *b_runNo;   //!
   TBranch        *b_subRunNo;   //!
   TBranch        *b_MC_Ev;   //!
   TBranch        *b_MC_cc;   //!
   TBranch        *b_MC_Q2;   //!
   TBranch        *b_MC_nuPDG;   //!
   TBranch        *b_MC_hit_nucleon;   //!
   TBranch        *b_mcgenie_npart;   //!
   TBranch        *b_mcgenie_id;   //!
   TBranch        *b_mcgenie_fate;   //!
   TBranch        *b_mcgenie_statusCode;   //!
   TBranch        *b_mcgenie_pdg;   //!
   TBranch        *b_mcgenie_mother;   //!
   TBranch        *b_mcgenie_startMomentum;   //!
   TBranch        *b_mcgenie_endMomentum;   //!
   TBranch        *b_mc_vertex;   //!
   TBranch        *b_mc_npart;   //!
   TBranch        *b_mc_id;   //!
   TBranch        *b_mc_pdg;   //!
   TBranch        *b_mc_statusCode;   //!
   TBranch        *b_mc_mother;   //!
   TBranch        *b_mc_startXYZT;   //!
   TBranch        *b_mc_endXYZT;   //!
   TBranch        *b_mc_startMomentum;   //!
   TBranch        *b_mc_endMomentum;   //!
   TBranch        *b_mc_Prange;   //!
   TBranch        *b_mc_truthlength;   //!
   TBranch        *b_mc_process;   //!
   TBranch        *b_mc_Endprocess;   //!
   TBranch        *b_n_vertices;   //!
   TBranch        *b_vertex;   //!
   TBranch        *b_n_reco_tracks;   //!
   TBranch        *b_n_decayVtx;   //!
   TBranch        *b_decayVtx;   //!
   TBranch        *b_track_vtx;   //!
   TBranch        *b_track_vtxDir;   //!
   TBranch        *b_track_end;   //!
   TBranch        *b_track_isContained;   //!
   TBranch        *b_track_length;   //!
   TBranch        *b_track_PIDA;   //!
   TBranch        *b_track_PID_pdg;   //!
   TBranch        *b_track_KE;   //!
   TBranch        *b_track_Prange;   //!
   TBranch        *b_track_bestplane;   //!
   TBranch        *b_n_cal_points;   //!
   TBranch        *b_track_dQ_dx;   //!
   TBranch        *b_track_dE_dx;   //!
   TBranch        *b_track_range;   //!
   TBranch        *b_track_pitch;   //!
   TBranch        *b_track_complet;   //!
   TBranch        *b_track_Efrac;   //!
   TBranch        *b_track_mcID;   //!
   TBranch        *b_track_mcPDG;   //!
   TBranch        *b_Em_ch;   //!
   TBranch        *b_Em_e;   //!
   TBranch        *b_trk_e;   //!
   TBranch        *b_Emichel_e;   //!
   TBranch        *b_n_showers;   //!
   TBranch        *b_sh_direction_X;   //!
   TBranch        *b_sh_direction_Y;   //!
   TBranch        *b_sh_direction_Z;   //!
   TBranch        *b_sh_start_X;   //!
   TBranch        *b_sh_start_Y;   //!
   TBranch        *b_sh_start_Z;   //!
   TBranch        *b_sh_energy;   //!
   TBranch        *b_sh_MIPenergy;   //!
   TBranch        *b_sh_dEdx;   //!
   TBranch        *b_sh_bestplane;   //!
   TBranch        *b_sh_length;   //!
   TBranch        *b_n_flashes;   //!
   TBranch        *b_flash_time;   //!
   TBranch        *b_flash_pe;   //!
   TBranch        *b_flash_ycenter;   //!
   TBranch        *b_flash_zcenter;   //!
   TBranch        *b_flash_ywidth;   //!
   TBranch        *b_flash_zwidth;   //!
   TBranch        *b_flash_timewidth;   //!
   TBranch        *b_flash_abstime;   //!
   TBranch        *b_flash_frame;   //!
   TBranch        *b_flash_PE_ndk;   //!
   TBranch        *b_flash_PE_Ar39;   //!

   NDKAna(TTree *tree=0);
   virtual ~NDKAna();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef NDKAna_cxx
NDKAna::NDKAna(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../Data_reco_v2/ndk_test.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../Data_reco_v2/ndk_test.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("../Data_reco_v2/ndk_test.root:/NDKAna");
      dir->GetObject("Event",tree);

   }
   Init(tree);
}

NDKAna::~NDKAna()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t NDKAna::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t NDKAna::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void NDKAna::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   mc_process = 0;
   mc_Endprocess = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("eventNo", &eventNo, &b_eventNo);
   fChain->SetBranchAddress("runNo", &runNo, &b_runNo);
   fChain->SetBranchAddress("subRunNo", &subRunNo, &b_subRunNo);
   fChain->SetBranchAddress("MC_Ev", &MC_Ev, &b_MC_Ev);
   fChain->SetBranchAddress("MC_cc", &MC_cc, &b_MC_cc);
   fChain->SetBranchAddress("MC_Q2", &MC_Q2, &b_MC_Q2);
   fChain->SetBranchAddress("MC_nuPDG", &MC_nuPDG, &b_MC_nuPDG);
   fChain->SetBranchAddress("MC_hit_nucleon", &MC_hit_nucleon, &b_MC_hit_nucleon);
   fChain->SetBranchAddress("mcgenie_npart", &mcgenie_npart, &b_mcgenie_npart);
   fChain->SetBranchAddress("mcgenie_id", mcgenie_id, &b_mcgenie_id);
   fChain->SetBranchAddress("mcgenie_fate", mcgenie_fate, &b_mcgenie_fate);
   fChain->SetBranchAddress("mcgenie_statusCode", mcgenie_statusCode, &b_mcgenie_statusCode);
   fChain->SetBranchAddress("mcgenie_pdg", mcgenie_pdg, &b_mcgenie_pdg);
   fChain->SetBranchAddress("mcgenie_mother", mcgenie_mother, &b_mcgenie_mother);
   fChain->SetBranchAddress("mcgenie_startMomentum", mcgenie_startMomentum, &b_mcgenie_startMomentum);
   fChain->SetBranchAddress("mcgenie_endMomentum", mcgenie_endMomentum, &b_mcgenie_endMomentum);
   fChain->SetBranchAddress("mc_vertex", mc_vertex, &b_mc_vertex);
   fChain->SetBranchAddress("mc_npart", &mc_npart, &b_mc_npart);
   fChain->SetBranchAddress("mc_id", mc_id, &b_mc_id);
   fChain->SetBranchAddress("mc_pdg", mc_pdg, &b_mc_pdg);
   fChain->SetBranchAddress("mc_statusCode", mc_statusCode, &b_mc_statusCode);
   fChain->SetBranchAddress("mc_mother", mc_mother, &b_mc_mother);
   fChain->SetBranchAddress("mc_startXYZT", mc_startXYZT, &b_mc_startXYZT);
   fChain->SetBranchAddress("mc_endXYZT", mc_endXYZT, &b_mc_endXYZT);
   fChain->SetBranchAddress("mc_startMomentum", mc_startMomentum, &b_mc_startMomentum);
   fChain->SetBranchAddress("mc_endMomentum", mc_endMomentum, &b_mc_endMomentum);
   fChain->SetBranchAddress("mc_Prange", mc_Prange, &b_mc_Prange);
   fChain->SetBranchAddress("mc_truthlength", mc_truthlength, &b_mc_truthlength);
   fChain->SetBranchAddress("mc_process", &mc_process, &b_mc_process);
   fChain->SetBranchAddress("mc_Endprocess", &mc_Endprocess, &b_mc_Endprocess);
   fChain->SetBranchAddress("n_vertices", &n_vertices, &b_n_vertices);
   fChain->SetBranchAddress("vertex", vertex, &b_vertex);
   fChain->SetBranchAddress("n_reco_tracks", &n_reco_tracks, &b_n_reco_tracks);
   fChain->SetBranchAddress("n_decayVtx", &n_decayVtx, &b_n_decayVtx);
   fChain->SetBranchAddress("decayVtx", &decayVtx, &b_decayVtx);
   fChain->SetBranchAddress("track_vtx", track_vtx, &b_track_vtx);
   fChain->SetBranchAddress("track_vtxDir", track_vtxDir, &b_track_vtxDir);
   fChain->SetBranchAddress("track_end", track_end, &b_track_end);
   fChain->SetBranchAddress("track_isContained", track_isContained, &b_track_isContained);
   fChain->SetBranchAddress("track_length", track_length, &b_track_length);
   fChain->SetBranchAddress("track_PIDA", track_PIDA, &b_track_PIDA);
   fChain->SetBranchAddress("track_PID_pdg", track_PID_pdg, &b_track_PID_pdg);
   fChain->SetBranchAddress("track_KE", track_KE, &b_track_KE);
   fChain->SetBranchAddress("track_Prange", track_Prange, &b_track_Prange);
   fChain->SetBranchAddress("track_bestplane", track_bestplane, &b_track_bestplane);
   fChain->SetBranchAddress("n_cal_points", &n_cal_points, &b_n_cal_points);
   fChain->SetBranchAddress("track_dQ_dx", track_dQ_dx, &b_track_dQ_dx);
   fChain->SetBranchAddress("track_dE_dx", track_dE_dx, &b_track_dE_dx);
   fChain->SetBranchAddress("track_range", track_range, &b_track_range);
   fChain->SetBranchAddress("track_pitch", track_pitch, &b_track_pitch);
   fChain->SetBranchAddress("track_complet", track_complet, &b_track_complet);
   fChain->SetBranchAddress("track_Efrac", track_Efrac, &b_track_Efrac);
   fChain->SetBranchAddress("track_mcID", track_mcID, &b_track_mcID);
   fChain->SetBranchAddress("track_mcPDG", track_mcPDG, &b_track_mcPDG);
   fChain->SetBranchAddress("Em_ch", &Em_ch, &b_Em_ch);
   fChain->SetBranchAddress("Em_e", &Em_e, &b_Em_e);
   fChain->SetBranchAddress("trk_e", &trk_e, &b_trk_e);
   fChain->SetBranchAddress("Emichel_e", &Emichel_e, &b_Emichel_e);
   fChain->SetBranchAddress("n_showers", &n_showers, &b_n_showers);
   fChain->SetBranchAddress("sh_direction_X", sh_direction_X, &b_sh_direction_X);
   fChain->SetBranchAddress("sh_direction_Y", sh_direction_Y, &b_sh_direction_Y);
   fChain->SetBranchAddress("sh_direction_Z", sh_direction_Z, &b_sh_direction_Z);
   fChain->SetBranchAddress("sh_start_X", sh_start_X, &b_sh_start_X);
   fChain->SetBranchAddress("sh_start_Y", sh_start_Y, &b_sh_start_Y);
   fChain->SetBranchAddress("sh_start_Z", sh_start_Z, &b_sh_start_Z);
   fChain->SetBranchAddress("sh_energy", sh_energy, &b_sh_energy);
   fChain->SetBranchAddress("sh_MIPenergy", sh_MIPenergy, &b_sh_MIPenergy);
   fChain->SetBranchAddress("sh_dEdx", sh_dEdx, &b_sh_dEdx);
   fChain->SetBranchAddress("sh_bestplane", sh_bestplane, &b_sh_bestplane);
   fChain->SetBranchAddress("sh_length", sh_length, &b_sh_length);
   fChain->SetBranchAddress("n_flashes", &n_flashes, &b_n_flashes);
   fChain->SetBranchAddress("flash_time", &flash_time, &b_flash_time);
   fChain->SetBranchAddress("flash_pe", &flash_pe, &b_flash_pe);
   fChain->SetBranchAddress("flash_ycenter", &flash_ycenter, &b_flash_ycenter);
   fChain->SetBranchAddress("flash_zcenter", &flash_zcenter, &b_flash_zcenter);
   fChain->SetBranchAddress("flash_ywidth", &flash_ywidth, &b_flash_ywidth);
   fChain->SetBranchAddress("flash_zwidth", &flash_zwidth, &b_flash_zwidth);
   fChain->SetBranchAddress("flash_timewidth", &flash_timewidth, &b_flash_timewidth);
   fChain->SetBranchAddress("flash_abstime", &flash_abstime, &b_flash_abstime);
   fChain->SetBranchAddress("flash_frame", &flash_frame, &b_flash_frame);
   fChain->SetBranchAddress("flash_PE_ndk", &flash_PE_ndk, &b_flash_PE_ndk);
   fChain->SetBranchAddress("flash_PE_Ar39", &flash_PE_Ar39, &b_flash_PE_Ar39);
   Notify();
}

Bool_t NDKAna::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void NDKAna::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t NDKAna::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef NDKAna_cxx
