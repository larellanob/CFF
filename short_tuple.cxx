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
  
  if(argc == 1) {    
    char File[200];
    ifstream in("dataFiles.txt", ios::in);
    if (!in) {
        cerr << "File Not Opened!" << endl;
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
  // FILES
  //////////////////////


  //////////////////////
  // VARIABLES and NTUPLES
  TDatabasePDG pdg;
  
  const char* NtupleName;
  
  TString     VarList = "TargType:Q2:Nu:Xb:W:SectorEl:ThetaPQ:PhiPQ:Zh:Pt:W2p:Xf:T:P:T4:deltaZ:evnt:pid";
  Int_t Nvar = VarList.CountChar(':')+1;
  Float_t *vars = new Float_t[Nvar];
  Int_t evntpos = 16;
  
  
  TVector3 *vert;
  TIdentificator *t = new TIdentificator(input);
  
  Long_t nEntries = (Long_t) input->GetEntries();

  TFile *output;
  if(simul_key == 0) {
    NtupleName = "ntuple_data";
    output = new TFile("local/prune_data_test.root", "RECREATE", "Data of particles");
  } else { 
    NtupleName = "particles_recons";
    output = new TFile("local/prune_simul.root", "RECREATE", "simul_Fe");
  }

  TNtuple *e_recons = new TNtuple("e_rec","Reconstructed Electrons","Q2:W:Nu:vxec:vyec:vzec:vxe:vye:vze:Pex:Pey:Pez:evnt");
  Float_t DataElec[e_recons->GetNvar()];

  TNtuple *particles_recons = new TNtuple(NtupleName,"reconstructed particles",VarList);
  TNtuple *particles_thrown = 0;
  TNtuple *e_thrown=0;
    
  
  if(simul_key == 1) {
    particles_thrown = new TNtuple("particles_thrown","thrown particles",VarList);
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
      bool RareTH = false;      // rare thrown event
      bool RareREC = false;     // rare reconstructed event
      for ( Int_t i=1; i < nRows; i++ )
	{
	  TString categ = t->GetCategorization(i,tt);
	  if ( t->GetCategorization(0,tt) != "electron" )  // POR QUE???
	    break;
	  if ( categ == "high energy pion +" || categ == "low energy pion +" )
	    {
	      PionEvent = true;
	      RareREC = true;
	      //cout << "reconstr. true at " << k << endl;
	      break;
	    }
	  
	}

      for ( Int_t i = 1; i < GSIMrows; i++ )
	{
	  if ( t->Id(i,1) == 8)
	    {
	      PionEvent = true;
	      RareTH = true;
	      //cout << "thrown true at " << k << endl;
	      break;
	    }
	}

      if ( RareREC == true && RareTH == false )
	{
	  //cout << "rare event at " << k << endl;
	  //cout << PionEvent << endl;
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
            
      

      
      
      if(nRows>0 && (t->GetCategorization(0,tt)) == "electron")  
	{
	  // variables reminder
	  // Q2:W:Nu:vxec:vyec:vzec:vxe:vye:vze:Pex:Pey:Pez:event
	  DataElec[0] = t -> Q2();
	  DataElec[1] = t -> W();
	  DataElec[2] = t -> Nu();
	  vert = t->GetCorrectedVert();
	  Float_t vxec=vert->X(); 
	  Float_t vyec=vert->Y(); 
	  Float_t vzec=vert->Z(); 
	  DataElec[3] = vxec; 
	  DataElec[4] = vyec; 
	  DataElec[5] = vzec;
	  DataElec[6] = t->X(0);
	  DataElec[7] = t->Y(0);
	  DataElec[8] = t->Z(0);
	  DataElec[9] = t -> Px(0);
	  DataElec[10] = t -> Py(0);
	  DataElec[11] = t -> Pz(0);
	  DataElec[12] = k;
	  
	  e_recons->Fill(DataElec);

	  // fill with zeroes if needed
	  if ( nRows == 1 )
	    {
	      for ( Int_t l = 1; l < GSIMrows; l++ )
		{
		  for ( Int_t ll = 0; ll < Nvar; ll++)
		    vars[ll] = 0;
		  vars[evntpos] = k;
		  particles_recons->Fill(vars);
		}
	    }
	  for (Int_t i = 1; i < nRows; i++) 
	    {
	      
	      TString category = t->GetCategorization(i,tt);
	      
	      // if (category == "gamma" || category == "pi-" || category == "high energy pion +" || category == "low energy pion +" || category == "s_electron" || category == "positron") 
	      
	      // NTUPLE_ACCEPT FILLING
	      Int_t f = 0;
	      vars[f] = t -> ElecVertTarg(); f++;
	      vars[f] = t -> Q2(); f++;
	      vars[f] = t -> Nu(); f++;
	      vars[f] = t -> Xb(); f++;
	      vars[f] = t -> W(); f++;
	      vars[f] = t -> Sector(0); f++;
	      vars[f] = t -> ThetaPQ(i); f++;
	      vars[f] = t -> PhiPQ(i); f++;
	      vars[f] = t -> Zh(i); f++;
	      vars[f] = TMath::Sqrt(t -> Pt2(i)); f++;
	      vars[f] = t -> Mx2(i); f++;
	      vars[f] = t -> Xf(i); f++;
	      vars[f] = t -> T(i); f++;
	      vars[f] = t -> Momentum(i); f++;
	      vars[f] = t -> TimeCorr4(0.139570,i); f++;
	      vars[f] = (t -> Z(i)) - (t -> Z(0)); f++;
	      vars[f] = k; f++;	      
	      vars[f] = ((category == "gamma")?22:
			  ((category == "pi-")?-211:
			   (( category == "high energy pion +" || category == "low energy pion +")?211:
			    ((category == "s_electron")?11:-11)
			    )
			   )
			 ); f++;
	      particles_recons->Fill(vars);


	      if ( i  == nRows-1 && GSIMrows > nRows )
		{
		  for ( Int_t l = 0; l < to_fill; l++ )
		    {
		      for ( Int_t ll = 0; ll < Nvar; ll++)
			vars[ll] = 0;
		      vars[evntpos] = k;
		      particles_recons->Fill(vars);
		    }
		}

	      
	    } // end: for (Int_t i = 1; i < nRows; i++)
	} // end:  if(nRows>0 && (t->GetCategorization(0,tt)) == "electron")
      
      else // if(nRows>0 && (t->GetCategorization(0,tt)) != "electron")
	{
	  for ( Int_t l = 0; l < 13; l++ )
	    DataElec[l] = 0;
	  e_recons->Fill(DataElec);

	  for ( Int_t i = 1; i < std::max(GSIMrows,nRows); i++ )
	    {
	      for ( Int_t ll = 0; ll < Nvar; ll++ )
		vars[ll] = 0;
	      vars[evntpos] = k;
	      //vars[Nvar] = k;
	      particles_recons->Fill(vars);
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
	      DataElec[0] = t -> Q2(1);
	      DataElec[1] = t -> W(1);
	      DataElec[2] = t -> Nu(1);
	      DataElec[3] = t -> Z(0,1);
	      DataElec[4] = t -> Px(0,1);
	      DataElec[5] = t -> Py(0,1);
	      DataElec[6] = t -> Pz(0,1);
	      DataElec[12] = k;
	      
	      e_thrown->Fill(DataElec);
	    }
	  else
	    {
	      DataElec[0] = 0;
	      DataElec[1] = 0;
	      DataElec[2] = 0;
	      DataElec[3] = 0;
	      DataElec[4] = 0;
	      DataElec[5] = 0;
	      DataElec[6] = 0;
	      DataElec[7] = 0;
	      DataElec[12] = 0;
	      
	      e_thrown->Fill(DataElec);
	    }
	  
	  for( Int_t i=1; i < GSIMrows; i++ )
	    {
	      Int_t f = 0;
	      vars[f] = t -> ElecVertTarg(1); f++;
	      vars[f] = t -> Q2(1); f++;
	      vars[f] = t -> Nu(1); f++;
	      vars[f] = t -> Xb(1); f++;
	      vars[f] = t -> W(1); f++;
	      vars[f] = t -> Sector(0,1); f++;
	      vars[f] = t -> ThetaPQ(i,1); f++;
	      vars[f] = t -> PhiPQ(i,1); f++;
	      vars[f] = t -> Zh(i,1); f++;
	      vars[f] = TMath::Sqrt(t -> Pt2(i,1)); f++;
	      vars[f] = t -> Mx2(i,1); f++;
	      vars[f] = t -> Xf(i,1); f++;
	      vars[f] = t -> T(i,1); f++;
	      vars[f] = t -> Momentum(i,1); f++;
	      vars[f] = 0; f++;//t -> TimeCorr4(0.139570,i);
	      vars[f] = (t -> Z(i,1)) - (t -> Z(0,1)); f++;
	      vars[f] = k; f++;
	      vars[f] = t -> Id(i,1); f++;
	      particles_thrown->Fill(vars);
	    }
	  
	  
	  if ( nRows > GSIMrows )
	    {
	      for ( Int_t i = 0; i < to_fill; i++ )
		{
		  for ( Int_t ll = 0; ll < Nvar; ll++ )
		    vars[ll] = 0;
		  vars[evntpos] = k;
		  particles_thrown->Fill(vars);
		}
	    }
	} // if( simul_key == 1 )
      else
	{
	  for ( Int_t i = 0; i < 13; i++ )
	    DataElec[i] = 0;
	  e_thrown->Fill(DataElec);
	  for ( Int_t i = 1; i < std::max(GSIMrows,nRows); i++ )
	    {
	      for ( Int_t ll = 0; ll < Nvar; ll++ )
		vars[ll] = 0;
	      vars[evntpos] = k;
	      particles_thrown->Fill(vars);
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
