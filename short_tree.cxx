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

void SetParticleVars(Float_t * particle_vars, TIdentificator * t, Int_t k, Int_t i, bool sim, TString category);
void SetElectronVars(Float_t  * e_vars, TIdentificator * t, Int_t k, bool sim);

int main(int argc, char **argv)
{
  TBenchmark bench;
  bool verbose = false;

  bench.Start("bench");
  
  //////////////////////
  // FILES
  bool simul_key = 0;
  
  TClasTool *input = new TClasTool();
  input->InitDSTReader("ROOTDSTR");
  
  if(argc == 1)
    {    
      char File[200];
      ifstream in("dataFiles.txt", ios::in);
      if (!in)
	{
	  cerr << "File Not Opened!\n";
	  cerr << "To use simulFiles.txt execute with any argument e.g. ./short_tuple asdf" << endl;
	  exit(1);
	}
      while (in >> File)
	{
	  input->Add(File);
	}
      in.close();
    }
  else
    {
      simul_key = 1;
      char File[200];
      ifstream in("simulFiles.txt", ios::in);
      if (!in)
	{
	  cerr << "File Not Opened!" << endl;
	  exit(1);
	}
      while (in >> File)
	{
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
  
  //const char* NtupleName;
  
  TString     VarList = "TargType:Q2:Nu:Xb:W:SectorEl:ThetaPQ:PhiPQ:Zh:Pt:W2p:Xf:T:P:T4:deltaZ:evnt:pid";
  Int_t Nvar = VarList.CountChar(':')+1;
  Float_t *particle_vars = new Float_t[Nvar];
  Int_t evntpos = 16;

  ///// TTrees implementation
  // OUTPUT FILE
  
  std::string out_filename_tree = "local/CFFTree_";
  
  // DATA
  std::string target_name2 = "data";
  // SIMULATION
  if ( simul_key == 1 )
    {
      target_name2 = string(argv[1]);
    }
  
  // General
  out_filename_tree.append(target_name2);
  std::string out_extension = ".root";
  out_filename_tree.append(out_extension);
  TFile out_tree(out_filename_tree.c_str(), "RECREATE", target_name2.c_str());

  /////// TTREE
  // Data
  
  TTree tree_data("tree_data","data of reconstructed particles");
  TTree tree_thrown("tree_thrown","tree of all thrown particles");
  TTree tree_accept("tree_accept","tree of all reconstructed particles");
  
  
  // point the tree to the right places
  Float_t tree_variables[Nvar][50]; // max of 50 particles per event
  TString leafname;
  TString leaflist;
  Ssiz_t from = 0;
  Int_t eventsize = 0;

  tree_data.Branch("eventsize",&eventsize,"eventsize/I");
  tree_thrown.Branch("eventsize",&eventsize,"eventsize/I");
  tree_accept.Branch("eventsize",&eventsize,"eventsize/I");
  
  for ( Int_t i = 0; i < Nvar; ++i )
    {
      VarList.Tokenize(leafname,from,":");
      leaflist = leafname+"[eventsize]/F"; // F is for float
      tree_data.Branch(leafname.Data(), &tree_variables[i][1], leaflist.Data());
      tree_thrown.Branch(leafname.Data(), &tree_variables[i][1], leaflist.Data()); 
      tree_accept.Branch(leafname.Data(), &tree_variables[i][1], leaflist.Data());
    }
  

  
  TVector3 *vert;
  TIdentificator *t = new TIdentificator(input);
  
  Long_t nEntries = (Long_t) input->GetEntries();
  
 
  TNtuple *e_recons = new TNtuple("e_rec","Reconstructed Electrons","Q2:W:Nu:vxec:vyec:vzec:vxe:vye:vze:Pex:Pey:Pez:evnt");
  Float_t e_vars[e_recons->GetNvar()];
  TNtuple *e_thrown = 0;
  if( simul_key == 1 )
    {
      e_thrown = new TNtuple("e_thrown","thrown Electrons","Q2:W:Nu:vxec:vyec:vzec:vxe:vye:vze:Pex:Pey:Pez:evnt");
    }
  
  // VARIABLES and NTUPLES
  //////////////////////  
  

  cout.width(4);
  input->Next();
  bool desync = false;
  Int_t GSIMrows, to_fill = -8;
  ProcInfo_t procinfo;
  for (Int_t k = 0; k < nEntries; k++) 
    {
      
      gSystem->GetProcInfo(&procinfo);
      Double_t totmem = (Double_t) (procinfo.fMemResident);
      
      Int_t nRows = input->GetNRows("EVNT");


      if ( simul_key == 1 )
	{
	  GSIMrows = input->GetNRows("GSIM");
	  to_fill = std::abs(GSIMrows-nRows); // in order to fill the difference
	  //cout << "PREVIO: " << GSIMrows << " " << nRows << " " << to_fill << endl;
	}
      const char * tt = "C";
      
      ////////////////////////
      // PION+ FILTER
      bool PionEvent = false;
      for ( Int_t i = 1; i < nRows; i++ )
	{
	  TString categ = t->GetCategorization(i,tt);
	  if ( t->GetCategorization(0,tt) != "electron" )  // electron == trigger
	    {
	      break;
	    }
	  if ( categ == "high energy pion +" || categ == "low energy pion +" )
	    {
	      PionEvent = true;
	      break;
	    }
	}

      if ( simul_key == 1 )
	{
	  GSIMrows = input->GetNRows("GSIM");
	  for ( Int_t i = 1; i < GSIMrows; i++ )
	    {
	      if ( t->Id(i,1) == 8)
		{
		  PionEvent = true;
		  break;
		}
	    }
	}
      
      if ( PionEvent == false )
	{
	  if ( verbose == true )
	    {
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
      
      if( nRows > 0 && (t->GetCategorization(0,tt)) == "electron" ) 
	{
	  
	  // flag -x
	  if ( simul_key == 1 && t->Id(0,1)!=3 )
	    {
	      input->Next();
	      continue;
	    }

	  // ELECTRON 
	  // set variables
	  SetElectronVars(e_vars, t, k, 0);
	  // fill ntuple
	  e_recons->Fill(e_vars);

	  if ( nRows == 1 )
	    {
	      eventsize=0;
	      if ( simul_key == 1 )
		tree_accept.Fill();
	      if ( simul_key == 0 )
		tree_data.Fill();
	    }
	  else if ( nRows != 1 )
	    {
	      for (Int_t i = 1; i < nRows; i++) 
		{
		  TString category = t->GetCategorization(i,tt);
		  // RECONSTRUCTED PARTICLES
		  // set variables		  
		  SetParticleVars(particle_vars,t,k,i,0,category);
		  for ( Int_t w = 0; w < Nvar; ++w )
		    {
		      tree_variables[w][i] = particle_vars[w];
		    }
		}
	      eventsize=nRows-1;
	      if ( simul_key == 1 )
		tree_accept.Fill();
	      if ( simul_key == 0 )
		tree_data.Fill();
	    }
	} // end:  if(nRows>0 && (t->GetCategorization(0,tt)) == "electron")
      
      else // if(nRows>0 && (t->GetCategorization(0,tt)) != "electron")
	{
	  if ( simul_key == 1  && t->Id(0,1)!= 3)
	    {
	      input->Next();
	      continue;
	    }
	  for ( Int_t l = 0; l < 13; l++ )
	    e_vars[l] = 0;
	  e_recons->Fill(e_vars);

	  eventsize=0;
	  if ( simul_key == 1 )
	    tree_accept.Fill();
	  if ( simul_key == 0 )
	    tree_data.Fill();
	  if ( simul_key == 1)
	    {
	      for ( Int_t i = 1; i < std::max(GSIMrows,nRows); i++ )
		{
		  for ( Int_t ll = 0; ll < Nvar; ll++ )
		    particle_vars[ll] = 0;
		  particle_vars[evntpos] = k;
		}
	    }
	}
      
      if ( simul_key == 0 )
	{
	  input->Next();
	  continue;
	}
      
      //////////////////////////////
      //////////////////////////////
      //////// T H R O W N /////////
      //////////////////////////////
      //////////////////////////////

      if( simul_key == 1 )//&& t -> Id(0,1)==3 /*&& t -> Q2(1) > 1. && t -> W(1) > 2. && t -> Nu(1) / 5.015 < 0.85*/ )
	{
	  if (t -> Id(0,1)==3 )
	    {
	      // THROWN ELECTRONS
	      // set variables
	      SetElectronVars(e_vars, t, k, 1);
	      // fill ntuple
	      e_thrown->Fill(e_vars);

	      // THROWN PARTICLES (for simulations)
	      if ( GSIMrows == 0 )
		{
		  eventsize=0;
		  tree_thrown.Fill();
		} 
	      for( Int_t i=1; i < GSIMrows; i++ )
		{
		  // set variables
		  SetParticleVars(particle_vars, t, k, i, 1, "");
		  for ( Int_t w = 0; w < Nvar; ++w )
		    {
		      tree_variables[w][i] = particle_vars[w];
		    }
		}
	      eventsize = GSIMrows-1;
	      tree_thrown.Fill();
	    }

	  // FLAG -x
	  else if ( t->Id(0,1) != 3 )
	    {
	      input->Next(); 
	      continue;
	    }
	}

      input->Next();

      // DESYNC check 
      if ( tree_accept.GetEntries() != tree_thrown.GetEntries() && desync == false )
	{
	  cout << "k = " << k << endl;
	  cout << "t->Id(0,1) = " << t->Id(0,1) << endl;
	  cout << "unsync starting " << k << endl;
	  cout << "ac = " << tree_accept.GetEntries() << " th = " << tree_thrown.GetEntries()  << endl;
	  desync = true;
	  break;
	}

    }
  // END of MAIN

  
  TTree version("version",VERSION);
  version.Write();
  
  if ( simul_key == 0 )
    {
      tree_data.Write();
      e_recons->Write();
    }

  if ( simul_key == 1 )
    {
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
  if ( sim == 0 )
    {
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
    }
  else if ( sim == 1 )
    {
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


void SetParticleVars(Float_t * particle_vars, TIdentificator * t, Int_t k, Int_t i, bool sim, TString category)
{
  Int_t f = 0;
  if ( sim == 0 )
    {
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
    }

  else if ( sim == 1 )
    {
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
    }
}
