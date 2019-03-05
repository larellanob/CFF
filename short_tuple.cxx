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


int main(int argc, char **argv)
{

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

  // OUTPUT FILE
  std::string out_filename = "local/CFF_";
  std::string target_name = string(argv[1]);
  std::string out_extension = ".root";
  out_filename.append(target_name);
  out_filename.append(out_extension);
  target_name.append(" simulation: events which include pi+");
  
  // FILES
  //////////////////////


  //////////////////////
  // VARIABLES and NTUPLES
  TDatabasePDG pdg;
  
  const char* NtupleName;
  
  TString     VarList = "TargType:Q2:Nu:Xb:W:SectorEl:ThetaPQ:PhiPQ:Zh:Pt:W2p:Xf:T:P:T4:deltaZ:evnt:pid";
  Int_t Nvar = VarList.CountChar(':')+1;
  Float_t *particle_vars = new Float_t[Nvar];
  Int_t evntpos = 16;
  
  
  TVector3 *vert;
  TIdentificator *t = new TIdentificator(input);
  
  Long_t nEntries = (Long_t) input->GetEntries();

  TFile *output;
  if(simul_key == 0)
    {
      NtupleName = "ntuple_data";
      output = new TFile("local/prune_data_test.root", "RECREATE", "Data of particles");
    }
  else
    { 
      NtupleName = "ntuple_accept";
      output = new TFile(out_filename.c_str(), "RECREATE", target_name.c_str());
    }

  TNtuple *e_recons = new TNtuple("e_rec","Reconstructed Electrons","Q2:W:Nu:vxec:vyec:vzec:vxe:vye:vze:Pex:Pey:Pez:evnt");
  Float_t e_vars[e_recons->GetNvar()];

  TNtuple *ntuple_accept = new TNtuple(NtupleName,"reconstructed particles",VarList);
  TNtuple *ntuple_thrown = 0;
  TNtuple *e_thrown=0;
    
  
  if(simul_key == 1)
    {
      ntuple_thrown = new TNtuple("ntuple_thrown","thrown particles",VarList);
      e_thrown = new TNtuple("e_thrown","thrown Electrons","Q2:W:Nu:vxec:vyec:vzec:vxe:vye:vze:Pex:Pey:Pez:evnt");
    }
  // VARIABLES and NTUPLES
  //////////////////////  
  

  cout.width(4);
  input->Next();

  for (Int_t k = 0; k < nEntries; k++) 
    {
      Int_t nRows = input->GetNRows("EVNT");
      //Int_t nRows = input->GetNRows("GSIM");
      
      Int_t GSIMrows = input->GetNRows("GSIM");
      Int_t to_fill = std::abs(GSIMrows-nRows); // in order to fill the difference
      const char * tt = "C";
      
      ////////////////////////
      // PION+ FILTER
      bool PionEvent = false;
      for ( Int_t i=1; i < nRows; i++ )
	{
	  TString categ = t->GetCategorization(i,tt);
	  if ( t->GetCategorization(0,tt) != "electron" )  // POR QUE???
	    break;
	  if ( categ == "high energy pion +" || categ == "low energy pion +" )
	    {
	      PionEvent = true;
	      break;
	    }
	}

      for ( Int_t i = 1; i < GSIMrows; i++ )
	{
	  if ( t->Id(i,1) == 8)
	    {
	      PionEvent = true;
	      break;
	    }
	}
      
      if ( PionEvent == false )
	{
	  cout<<std::right<<float(k+1)/nEntries*100<<"%\r";
	  cout.flush();
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
      
      // the second condition makes the code only work for simulation files 
      if(nRows>0 && (t->GetCategorization(0,tt)) == "electron" && t -> Id(0,1)==3 )
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
	  
	  e_recons->Fill(e_vars);

	  // fill with zeroes if needed
	  if ( nRows == 1 )
	    {
	      for ( Int_t l = 1; l < GSIMrows; l++ )
		{
		  for ( Int_t ll = 0; ll < Nvar; ll++)
		    particle_vars[ll] = 0;
		  particle_vars[evntpos] = k;
		  ntuple_accept->Fill(particle_vars);
		}
	    }
	  for (Int_t i = 1; i < nRows; i++) 
	    {
	      
	      TString category = t->GetCategorization(i,tt);

	      // other possible filters
	      // if (category == "gamma" || category == "pi-" || category == "high energy pion +" || category == "low energy pion +" || category == "s_electron" || category == "positron") 
	      
	      // NTUPLE_ACCEPT FILLING
	      Int_t f = 0;
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
	      ntuple_accept->Fill(particle_vars);


	      if ( i  == nRows-1 && GSIMrows > nRows )
		{
		  for ( Int_t l = 0; l < to_fill; l++ )
		    {
		      for ( Int_t ll = 0; ll < Nvar; ll++)
			particle_vars[ll] = 0;
		      particle_vars[evntpos] = k;
		      ntuple_accept->Fill(particle_vars);
		    }
		}

	      
	    } // end: for (Int_t i = 1; i < nRows; i++)
	} // end:  if(nRows>0 && (t->GetCategorization(0,tt)) == "electron")
      
      else // if(nRows>0 && (t->GetCategorization(0,tt)) != "electron")
	{
	  for ( Int_t l = 0; l < 13; l++ )
	    e_vars[l] = 0;
	  e_recons->Fill(e_vars);

	  for ( Int_t i = 1; i < std::max(GSIMrows,nRows); i++ )
	    {
	      for ( Int_t ll = 0; ll < Nvar; ll++ )
		particle_vars[ll] = 0;
	      particle_vars[evntpos] = k;
	      //particle_vars[Nvar] = k;
	      ntuple_accept->Fill(particle_vars);
	    }
	}
      
      
      
      //////////////////////////////
      //////////////////////////////
      //////// T H R O W N /////////
      //////////////////////////////
      //////////////////////////////
      
      
      
      // CLASEVENT->Show(150) electron no es el primero
      // y no lo toma
      if( simul_key == 1 )//&& t -> Id(0,1)==3 /*&& t -> Q2(1) > 1. && t -> W(1) > 2. && t -> Nu(1) / 5.015 < 0.85*/ )
	{
	  if (t -> Id(0,1)==3 )
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
	      
	      e_thrown->Fill(e_vars);

	      for( Int_t i=1; i < GSIMrows; i++ )
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
		  ntuple_thrown->Fill(particle_vars);
		}

	      if ( nRows > GSIMrows )
		{
		  for ( Int_t i = 0; i < to_fill; i++ )
		    {
		      for ( Int_t ll = 0; ll < Nvar; ll++ )
			particle_vars[ll] = 0;
		      particle_vars[evntpos] = k;
		      ntuple_thrown->Fill(particle_vars);
		    }
		}
	    }
	  
	  
	  else
	    {
	      for ( Int_t i = 0; i < 13; i++ )
		{
		  e_vars[i] = 0;
		}
	      for ( Int_t i = 1; i < GSIMrows; i++ )
		{
		  for ( Int_t j = 0;  j< Nvar; j++ )
		    {
		      if ( j != Nvar-2 )
			particle_vars[j] = 0;
		      else if ( j == Nvar-2 )
			particle_vars[j] = k;
		    }
		  ntuple_thrown->Fill(particle_vars);
		}
	      e_thrown->Fill(e_vars);
	      
	      if ( nRows > GSIMrows )
		{
		  for ( Int_t i = 0; i < to_fill; i++ )
		    {
		      for ( Int_t ll = 0; ll < Nvar; ll++ )
			particle_vars[ll] = 0;
		      particle_vars[evntpos] = k;
		      ntuple_thrown->Fill(particle_vars);
		    }
		}
	    }


	} // if( simul_key == 1 )
      else
	{
	  for ( Int_t i = 0; i < 13; i++ )
	    e_vars[i] = 0;
	  e_thrown->Fill(e_vars);
	  for ( Int_t i = 1; i < std::max(GSIMrows,nRows); i++ )
	    {
	      for ( Int_t ll = 0; ll < Nvar; ll++ )
		particle_vars[ll] = 0;
	      particle_vars[evntpos] = k;
	      ntuple_thrown->Fill(particle_vars);
	    } 
	}
      
      cout<<std::right<<float(k+1)/nEntries*100<<"%\r";
      cout.flush();
      input->Next();
    } // for (Int_t k = 0; k < nEntries; k++) 
  output->Write();
  output->Close();
  cout << "Done." << endl;
  return 0;
}
