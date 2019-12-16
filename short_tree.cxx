#include "Riostream.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TClasTool.h"
#include "TIdentificator.h"
#include "TMath.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TBenchmark.h"

// Variables used

bool g_Debug = false;
TString gVarList_particles =
  "TargType:Q2:Nu:Xb:W:SectorEl:ThetaPQ:PhiPQ:Zh:Pt:W2p:Xf:T:P:T4:deltaZ:evnt:pid:ThetaLab:PhiLab:Px:Py:Pz";
TString gVarList_selectron =
  "Q2:W:Nu:vxec:vyec:vzec:vxe:vye:vze:Pex:Pey:Pez:evnt";

void SetParticleVars(Float_t * particle_vars, TIdentificator * t, Int_t k, Int_t i, bool sim, TString category);
void SetElectronVars(Float_t  * e_vars, TIdentificator * t, Int_t k, bool sim);

void Clear_variables(Int_t & TargType,
		     Float_t & Q2,
		     Float_t & Nu,
		     Float_t & Xb,
		     Float_t & W,
		     Int_t & SectorEl,
		     std::vector<Float_t> & Zh,
		     std::vector<Float_t> & Pt,
		     std::vector<Float_t> & W2p,
		     std::vector<Float_t> & Xf,
		     std::vector<Float_t> & T,
		     std::vector<Float_t> & P,
		     std::vector<Float_t> & deltaZ,
		     std::vector<Float_t> & Px,
		     std::vector<Float_t> & Py,
		     std::vector<Float_t> & Pz,
		     std::vector<Int_t>   & pid,
		     Int_t & evnt,
		     std::vector<Float_t> & ThetaPQ,
		     std::vector<Float_t> & PhiPQ,
		     std::vector<Float_t> & Theta,
		     std::vector<Float_t> & Phi
		     )
{
  TargType = -1;
  Q2 = 0;
  Nu = 0;
  Xb = 0;
  W  = 0;
  evnt= 0;
  SectorEl = -1;
  Zh.clear();
  Pt.clear();
  W2p.clear();
  Xf.clear();
  T.clear();
  P.clear();
  deltaZ.clear();
  Px.clear();
  Py.clear();
  Pz.clear();
  pid.clear();
  ThetaPQ.clear();
  PhiPQ.clear();
  Theta.clear();
  Phi.clear();
}
  

int main(int argc, char **argv)
{

  TBenchmark bench;
  bool verbose = true;

  bench.Start("bench");
  
  //////////////////////
  // FILES
  bool simul_key = 0;
  
  TClasTool *input = new TClasTool();
  input->InitDSTReader("ROOTDSTR");

  if(argc == 2) {
    cout << "one arguments: target" << endl;
    char File[200];
    ifstream in("dataFiles.txt", ios::in);
    if (!in) {
      cerr << "File Not Opened!\n";
      cerr << "To use simulFiles.txt execute with any argument e.g. ./short_tuple asdf" << endl;
      exit(1);
    }
    while (in >> File) {
      input->Add(File);
    }
    in.close();
  } else {
    simul_key = 1;
    char File[200];
    ifstream in("simulFiles.txt", ios::in);
    if (!in) {
      cerr << "File Not Opened!" << endl;
      exit(1);
    }
    while (in >> File) {
      input->Add(File);
    }
    in.close();
  }

  //////////////////////
  //////////////////////
  ////// FILES /////////
  //////////////////////
  //////////////////////

  
  //////////////////////
  // VARIABLES and NTUPLES
  TDatabasePDG pdg;
    
  ///// TTrees implementation
  // OUTPUT FILE
  
  std::string out_filename_tree = "local/CFFTree_";
  
  // DATA
  std::string target_name = string(argv[1]);
  std::string data_simulation = target_name;
  if ( !simul_key ){
    data_simulation = data_simulation+"_data";
  } else {
    data_simulation = data_simulation+"_simul";
  }
    
  // General
  out_filename_tree.append(data_simulation);
  std::string out_extension = ".root";
  out_filename_tree.append(out_extension);
  TFile out_tree(out_filename_tree.c_str(), "RECREATE", data_simulation.c_str());

  /////// TTREE
  // Data
  
  TTree tree_data("tree_data","data of reconstructed particles");
  TTree tree_thrown("tree_thrown","tree of all thrown particles");
  TTree tree_accept("tree_accept","tree of all reconstructed particles");
  
  // VECTORS REWRITE
  // "TargType:Q2:Nu:Xb:W:SectorEl:ThetaPQ:PhiPQ:Zh:Pt:W2p:Xf:T:P:T4:deltaZ:evnt:pid:ThetaLab:PhiLab:Px:Py:Pz";
  // "Q2:W:Nu:vxec:vyec:vzec:vxe:vye:vze:Pex:Pey:Pez:evnt";
  Int_t TargType;
  Float_t Q2;
  Float_t Nu;
  Float_t Xb;
  Float_t W;
  Int_t SectorEl;
  std::vector<Float_t> Zh;
  std::vector<Float_t> Pt;
  std::vector<Float_t> W2p;
  std::vector<Float_t> Xf;
  std::vector<Float_t> T;
  std::vector<Float_t> P;
  std::vector<Float_t> deltaZ;
  std::vector<Float_t> Px;
  std::vector<Float_t> Py;
  std::vector<Float_t> Pz;
  std::vector<Int_t>   pid;
  Int_t evnt;
  std::vector<Float_t> ThetaPQ;
  std::vector<Float_t> PhiPQ;
  std::vector<Float_t> Theta;
  std::vector<Float_t> Phi;

  std::vector<TTree *> trees;
  
  trees.push_back(&tree_data);
  trees.push_back(&tree_thrown);
  trees.push_back(&tree_accept);

  for ( auto tree: trees ) {
    tree->Branch("TargType",&TargType);
    tree->Branch("Q2",&Q2);
    tree->Branch("Nu",&Nu);
    tree->Branch("Xb",&Xb);
    tree->Branch("W",&W);
    tree->Branch("SectorEl",&SectorEl);
    tree->Branch("Zh",&Zh);
    tree->Branch("Pt",&Pt);
    tree->Branch("W2p",&W2p);
    tree->Branch("Xf",&Xf);
    tree->Branch("T",&T);
    tree->Branch("P",&P);
    tree->Branch("deltaZ",&deltaZ);
    tree->Branch("Px",&Px);
    tree->Branch("Py",&Py);
    tree->Branch("Pz",&Pz);
    tree->Branch("pid",&pid);
    tree->Branch("evnt",&evnt);
    tree->Branch("ThetaPQ",&ThetaPQ);
    tree->Branch("PhiPQ",&PhiPQ);
    tree->Branch("Theta",&Theta);
    tree->Branch("Phi",&Phi);
  }
  
  TIdentificator *t = new TIdentificator(input);
  
  Long_t nEntries = (Long_t) input->GetEntries();
  
  TNtuple *e_recons = new TNtuple("e_rec","Reconstructed Electrons","Q2:W:Nu:vxec:vyec:vzec:vxe:vye:vze:Pex:Pey:Pez:evnt");
  Float_t e_vars[e_recons->GetNvar()];
  
  TNtuple *e_thrown = 0;
  if( simul_key == 1 ) {
    e_thrown = new TNtuple("e_thrown","thrown Electrons","Q2:W:Nu:vxec:vyec:vzec:vxe:vye:vze:Pex:Pey:Pez:evnt");
  }
  
  // VARIABLES and NTUPLES
  //////////////////////  
  

  cout.width(4);
  input->Next();
  bool desync = false;
  Int_t GSIMrows = 0;
  ProcInfo_t procinfo;
  Int_t EventCounter = 0;
  for (Int_t k = 0; k < nEntries; k++) {

    gSystem->GetProcInfo(&procinfo);
    Double_t totmem = (Double_t) (procinfo.fMemResident);
    
    Int_t nRows = input->GetNRows("EVNT");
    
    Clear_variables(TargType,Q2,Nu, Xb, W, SectorEl, Zh, Pt, W2p, Xf, T, P, deltaZ, Px, Py, Pz, pid, evnt, ThetaPQ, PhiPQ, Theta, Phi);
      
    if ( simul_key == 1 ) {
      GSIMrows = input->GetNRows("GSIM");
      if ( GSIMrows == 0 ){
	input->Next();
	continue;
      }
    }
    //const char * tt = "C";
    const char * tt = target_name.c_str();
    
    //////////////////////////
    // FIRST PARTICLE ELECTRON
    //////////////////////////
    if ( t->GetCategorization(0,tt) != "electron" ) {
      input->Next();
      continue;
    }
    
    if ( simul_key  && t->Id(0,1) != 3 ) {
      input->Next();
      continue;
    }

    ////////////////////////
    // TargType 1 or 2
    ////////////////////////
    if ( t->ElecVertTarg() != 1 &&  t->ElecVertTarg() != 2 ) {
      input->Next();
      continue;
    }

    if ( simul_key && t->ElecVertTarg(1) != 1 &&  t->ElecVertTarg(1) != 2 ) {
      input->Next();
      continue;
    }
    
    ////////////////////////
    // PION+ FILTER
    bool PionEvent      = false;
    bool PionMinusEvent = false;
    bool PhotonEvent    = false;
    bool ParticleSelection = false;
    for ( Int_t i = 1; i < nRows; i++ ) {
      TString categ = t->GetCategorization(i,tt);
      if ( categ == "high energy pion +"
	   || categ == "low energy pion +"
	   ) {
	PionEvent = true;
      }
      if ( categ == "pi-" ) {
	PionMinusEvent = true;
      }
      if ( categ == "photon"
	   || categ == "gamma"
	   ) {
	PhotonEvent = true;
      }
      if ( PionEvent
	   || PionMinusEvent
	   ) {
	ParticleSelection = true;
	break;
      }
    }
    
    if ( simul_key == 1 ) {
      GSIMrows = input->GetNRows("GSIM");
      for ( Int_t i = 1; i < GSIMrows; i++ ) {
	if ( t->Id(i,1) == 8) {
	  PionEvent = true;
	}
	if ( t->Id(i,1) == 9 ) {
	  PionMinusEvent = true;
	}
	if ( t->Id(i,1) == 1 ) {
	  PhotonEvent = true;
	}
	if ( PionEvent
	     || PionMinusEvent
	     ) {
	  ParticleSelection = true;
	  break;
	}
      }
    }
    
    if ( ParticleSelection == false ) {
      if ( verbose == true && EventCounter%100 == 0 ) {
	cout << std::right <<  std::setw(12) << float(k+1)/nEntries*100.0 << "% mem: " << totmem << "\r";
	cout.flush();
      }
      input->Next();
      continue;
    }
    // PION+ FILTER
    /////////////////////////

    ///////////////////////////////
    ///////////////////////////////
    // R E C O N S T R U C T E D //
    ///////////////////////////////
    ///////////////////////////////
    
    if ( nRows == 0 || nRows == 1 ) {
      //eventsize=0;
      for ( Int_t l = 0; l < 13; l++ ) {
	e_vars[l] = 0;
      }

      pid.emplace_back(0);
      TargType = 0;
      Q2 = 0;
      Nu = 0;
      Xb = 0;
      W  = 0;
      SectorEl = 0;
      Zh.emplace_back(0);
      Pt.emplace_back(0);
      W2p.emplace_back(0);
      Xf.emplace_back(0);
      T.emplace_back(0);
      P.emplace_back(0);
      deltaZ.emplace_back(0);
      Px.emplace_back(0);
      Py.emplace_back(0);
      Pz.emplace_back(0);
      evnt = k;
      ThetaPQ.emplace_back(0);
      PhiPQ.emplace_back(0);
      Theta.emplace_back(0);
      Phi.emplace_back(0);
      
      e_recons->Fill(e_vars);

      if (t -> Id(0,1)==3 ) {
	if ( simul_key == 1 )
	  tree_accept.Fill();
	if ( simul_key == 0 )
	  tree_data.Fill();
      }
    }
    
    if( nRows > 1 && (t->GetCategorization(0,tt)) == "electron" ) {
      
      // flag -x
      if ( simul_key == 1 && t->Id(0,1)!=3 ) {
	input->Next();
	continue;
      }
      
      // ELECTRON 
      // set variables
      SetElectronVars(e_vars, t, k, 0);
      // fill ntuple
      e_recons->Fill(e_vars);
      EventCounter++;
      // RECONSTRUCTED PARTICLES
      for (Int_t i = 1; i < nRows; i++) {
	TString category = t->GetCategorization(i,tt);
	// set variables		  
	//SetParticleVars(particle_vars,t,k,i,0,category);
	if ( category == "pi-" ) {
	  pid.emplace_back(-211);
	} else if ( category == "high energy pion +" 
		    || category == "low energy pion +"
		    ) { 
	  pid.emplace_back(211);
	} else if ( category == "photon" ) {
	  pid.emplace_back(21);
	} else if ( category == "gamma" ) {
	  pid.emplace_back(22);
	} else {
	  continue;
	}
	TargType = t->ElecVertTarg();
	Q2 = t->Q2();
	Nu = t->Nu();
	Xb = t->Xb();
	W  = t->W() ;
	SectorEl = t->Sector(0);
	Zh.emplace_back(t->Zh(i));
	Pt.emplace_back(TMath::Sqrt(t->Pt2(i)));
	W2p.emplace_back(t->Mx2(i));
	Xf.emplace_back(t->Xf(i));
	T.emplace_back(t->T(i));
	P.emplace_back(t->Momentum(i));
	deltaZ.emplace_back( (t->Z(i)) - (t->Z(0)) );
	Px.emplace_back(t->Px(i));
	Py.emplace_back(t->Py(i));
	Pz.emplace_back(t->Pz(i));
	evnt = k;
	ThetaPQ.emplace_back(t->ThetaPQ(i));
	PhiPQ.emplace_back(t->PhiPQ(i));
	Theta.emplace_back(t->ThetaLab(i));
	Phi.emplace_back(t->PhiLab(i));
      }
	
      if ( simul_key == 1 )
	tree_accept.Fill();
      if ( simul_key == 0 )
	tree_data.Fill();

    } else if ( nRows > 1 && (t->GetCategorization(0,tt)) != "electron" ) {
      // also throwing away if first particle not electron
      input->Next();
      continue;
      /*
      for ( Int_t l = 0; l < 13; l++ )
	e_vars[l] = 0;
      e_recons->Fill(e_vars);
      */
    }

    
    if ( simul_key == 0 ) {
      input->Next();
      continue;
    }
    
    //////////////////////////////
    //////////////////////////////
    //////// T H R O W N /////////
    //////////////////////////////
    //////////////////////////////


    Clear_variables(TargType,Q2,Nu, Xb, W, SectorEl, Zh, Pt, W2p, Xf, T, P, deltaZ, Px, Py, Pz, pid, evnt, ThetaPQ, PhiPQ, Theta, Phi);

    // very old comment:
    //&& t -> Id(0,1)==3 /*&& t -> Q2(1) > 1. && t -> W(1) > 2. && t -> Nu(1) / 5.015 < 0.85*/ )
    if( simul_key == 1 ) { 
      
      if (t -> Id(0,1)==3 ) {

	// THROWN ELECTRONS
	// set variables
	SetElectronVars(e_vars, t, k, 1);
	// fill ntuple
	e_thrown->Fill(e_vars);
	// THROWN PARTICLES (for simulations)
	/*
	if ( GSIMrows == 0 ) {
	  //eventsize=0;
	  tree_thrown.Fill();
	}
	*/
	for( Int_t i=1; i < GSIMrows; i++ ) {
	  if ( t->Id(i,1) == 8 ) {
	    pid.emplace_back(211);
	  } else if ( t->Id(i,1) == 9 ) {
	    pid.emplace_back(-211);
	  } else if ( t->Id(i,1) == 1 ) {	    
	    pid.emplace_back(22);
	  } else {
	    continue;
	  }
	  TargType = t->ElecVertTarg(1);
	  Q2 = t->Q2(1);
	  Nu = t->Nu(1);
	  Xb = t->Xb(1);
	  W  = t->W(1) ;
	  SectorEl = t->Sector(0,1);
	  Zh.emplace_back(t->Zh(i,1));
	  Pt.emplace_back( TMath::Sqrt(t->Pt2(i,1)) );
	  W2p.emplace_back(t->Mx2(i,1));
	  Xf.emplace_back(t->Xf(i,1));
	  T.emplace_back(t->T(i,1));
	  P.emplace_back(t->Momentum(i,1));
	  deltaZ.emplace_back( (t->Z(i,1)) - (t->Z(0,1)) );
	  Px.emplace_back(t->Px(i,1));
	  Py.emplace_back(t->Py(i,1));
	  Pz.emplace_back(t->Pz(i,1));
	  evnt = k;
	  ThetaPQ.emplace_back(t->ThetaPQ(i,1));
	  PhiPQ.emplace_back(t->PhiPQ(i,1));
	  Theta.emplace_back(t->ThetaLab(i,1));
	  Phi.emplace_back(t->PhiLab(i,1));
	}
	tree_thrown.Fill();
      } else if ( t->Id(0,1) != 3 ) {
	// FLAG -x
	input->Next(); 
	continue;
      }
    }
    
    input->Next();
    
    // DESYNC check 
    if ( tree_accept.GetEntries() != tree_thrown.GetEntries() && desync == false ) {
      cout << "k = " << k << endl;
      cout << "t->GetCat = " << t->GetCategorization(0,tt) << endl;
      cout << "t->Id(0,1) = " << t->Id(0,1) << endl;
      cout << "unsync starting " << k << endl;
      cout << "ac = " << tree_accept.GetEntries() << " th = " << tree_thrown.GetEntries()  << endl;
      desync = true;
      break;
    }

  }
  // END of FILES READING
  TTree version("version",VERSION);
  version.Write();
  
  if ( simul_key == 0 ) {
    tree_data.Write();
    e_recons->Write();
  }
  
  if ( simul_key == 1 ) {
    tree_thrown.Write();
    tree_accept.Write();
    e_thrown->Write();
    e_recons->Write();
  }
    
  //out_tree.Write();
  out_tree.Close();
  cout << "Done." << endl;
  bench.Show("bench");
  return 0;
}


void SetElectronVars(Float_t  * e_vars, TIdentificator * t, Int_t k, bool sim)
{
  TVector3 *vert;
  if ( sim == 0 ) {
    // variables reminder
    //  0:1: 2:   3:   4:   5:  6:  7:  8:  9: 10: 11:   12
    // Q2:W:Nu:vxec:vyec:vzec:vxe:vye:vze:Pex:Pey:Pez:event
    e_vars[0] = t -> Q2();
    e_vars[1] = t -> W();
    e_vars[2] = t -> Nu();
    vert = t->GetCorrectedVert();
    Float_t vxec=vert->X(); 
    Float_t vyec=vert->Y(); 
    Float_t vzec=vert->Z(); 
    e_vars[3] = vxec; 
    e_vars[4] = vyec; 
    e_vars[5] = vzec;
    e_vars[6] = t->X(0);
    e_vars[7] = t->Y(0);
    e_vars[8] = t->Z(0);
    e_vars[9] = t -> Px(0);
    e_vars[10] = t -> Py(0);
    e_vars[11] = t -> Pz(0);
    e_vars[12] = k;
  } else if ( sim == 1 ) {
    e_vars[0] = t -> Q2(1);
    e_vars[1] = t -> W(1);
    e_vars[2] = t -> Nu(1);
    e_vars[3] = 0;
    e_vars[4] = 0;
    e_vars[5] = 0;
    e_vars[6] = t -> X(0,1);
    e_vars[7] = t -> Y(0,1);
    e_vars[8] = t -> Z(0,1);
    e_vars[9] = t -> Px(0,1);
    e_vars[10] = t -> Py(0,1);
    e_vars[11] = t -> Pz(0,1);
    e_vars[12] = k;
  }
}

// legacy code:
void SetParticleVars(Float_t * particle_vars, TIdentificator * t, Int_t k, Int_t i, bool sim, TString category)
{
  Int_t f = 0;
  if ( sim == 0 ) {
    particle_vars[f] = t -> ElecVertTarg(); f++;
    particle_vars[f] = t -> Q2(); f++;
    particle_vars[f] = t -> Nu(); f++;
    particle_vars[f] = t -> Xb(); f++;
    particle_vars[f] = t -> W(); f++;
    particle_vars[f] = t -> Sector(0); f++;
    particle_vars[f] = t -> ThetaPQ(i); f++;
    particle_vars[f] = t -> PhiPQ(i); f++;
    particle_vars[f] = t -> Zh(i); f++;
    particle_vars[f] = TMath::Sqrt(t -> Pt2(i)); f++;
    particle_vars[f] = t -> Mx2(i); f++;
    particle_vars[f] = t -> Xf(i); f++;
    particle_vars[f] = t -> T(i); f++;
    particle_vars[f] = t -> Momentum(i); f++;
    particle_vars[f] = t -> TimeCorr4(0.139570,i); f++;
    particle_vars[f] = (t -> Z(i)) - (t -> Z(0)); f++;
    particle_vars[f] = k; f++;	      
    particle_vars[f] = ((category == "gamma")?22:
			((category == "pi-")?-211:
			 (( category == "high energy pion +" || category == "low energy pion +")?211:
			  ((category == "s_electron")?11:-11)
			  )
			 )
			); f++;
    particle_vars[f] = t -> ThetaLab(i); f++;
    particle_vars[f] = t -> PhiLab(i); f++;
    particle_vars[f] = t -> Px(i); f++;
    particle_vars[f] = t -> Py(i); f++;
    particle_vars[f] = t -> Pz(i); f++;
  } else if ( sim == 1 ) {
    Int_t f = 0;
    particle_vars[f] = t -> ElecVertTarg(1); f++;
    particle_vars[f] = t -> Q2(1); f++;
    particle_vars[f] = t -> Nu(1); f++;
    particle_vars[f] = t -> Xb(1); f++;
    particle_vars[f] = t -> W(1); f++;
    particle_vars[f] = t -> Sector(0,1); f++;
    particle_vars[f] = t -> ThetaPQ(i,1); f++;
    particle_vars[f] = t -> PhiPQ(i,1); f++;
    particle_vars[f] = t -> Zh(i,1); f++;
    particle_vars[f] = TMath::Sqrt(t -> Pt2(i,1)); f++;
    particle_vars[f] = t -> Mx2(i,1); f++;
    particle_vars[f] = t -> Xf(i,1); f++;
    particle_vars[f] = t -> T(i,1); f++;
    particle_vars[f] = t -> Momentum(i,1); f++;
    particle_vars[f] = 0; f++;//t -> TimeCorr4(0.139570,i);
    particle_vars[f] = (t -> Z(i,1)) - (t -> Z(0,1)); f++;
    particle_vars[f] = k; f++;
    particle_vars[f] = t -> Id(i,1); f++;
    particle_vars[f] = t -> ThetaLab(i,1); f++;
    particle_vars[f] = t -> PhiLab(i,1); f++;
    particle_vars[f] = t -> Px(i,1); f++;
    particle_vars[f] = t -> Py(i,1); f++;
    particle_vars[f] = t -> Pz(i,1); f++;
  }
}
