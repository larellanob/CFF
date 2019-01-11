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
  //gROOT->Reset(); a

  bool simul_key = 0;
  
  TClasTool *input = new TClasTool();
  input->InitDSTReader("ROOTDSTR");
  
  if(argc == 1) {    
    char File[200];
    //system("ls -l *.root > dataFiles.txt");
    ifstream in("dataFiles.txt", ios::in);
    if (!in) {
        cerr << "File Not Opened!" << endl;
        exit(1);
    }
    while (in >> File) {
        input->Add(File);
    }
    in.close();
    //system("rm dataFiles.txt");
  } else {
    simul_key = 1;
    char File[200];
    //system("ls -1 *.root > simulFiles.txt");
    ifstream in("simulFiles.txt", ios::in);
    if (!in) {
        cerr << "File Not Opened!" << endl;
        exit(1);
    }
    while (in >> File) {
        input->Add(File);
    }
    in.close();
    //system("rm simulFiles.txt");
  }
  TDatabasePDG pdg;
  //Double_t kMe =pdg.GetParticle(11)->Mass();
  const char* NtupleName;
  //Float_t is_pion=0;
  TString     VarList = "TargType:Q2:Nu:Xb:W:SectorEl:ThetaPQ:PhiPQ:Zh:Pt:W2p:Xf:T:P:T4:deltaZ:evnt:pid";
  Int_t Nvar = VarList.CountChar(':')+1;
  Float_t *vars = new Float_t[Nvar];
  Int_t evntpos = 16;
  //Float_t vars2[Nvar];
  // ORIGINAL FULL VARIABLE LIST
  //"TargType:Q2:Nu:Xb:W:SectorEl:ThetaPQ:PhiPQ:Zh:Pt:W2p:Xf:T:P:T4:deltaZ:E:Ee:Pe:Ect:Sct:Ecr:Scr:evnt:Px:Py:Pz:Xe:Ye:Ze:Xec:Yec:Zec:TEc:ECX:ECY:ECZ:Pex:Pey:Pez:Ein:Eout:Eine:Eoute:pid:Betta:vxh:vyh:vzh:is_pion:k";
  
  TVector3 *vert;
  TIdentificator *t = new TIdentificator(input);
  
  Long_t nEntries = (Long_t) input->GetEntries();

  TFile *output;
  if(simul_key == 0) {
    NtupleName = "ntuple_data";
    output = new TFile("local/prune_data_test.root", "RECREATE", "Data of particles");
  } else { 
    NtupleName = "ntuple_accept";
    output = new TFile("local/prune_simul.root", "RECREATE", "simul_Fe");
  }

  TNtuple *tElec = new TNtuple("e_rec","Reconstructed Electrons","Q2:W:Nu:vxec:vyec:vzec:vxe:vye:vze:Pex:Pey:Pez:evnt");
  Float_t DataElec[tElec->GetNvar()];

  TNtuple *ntuple = new TNtuple(NtupleName,"reconstructed particles",VarList);
  TNtuple *ntuple_thrown = 0;
  TNtuple *e_thrown=0;
  //TTree *tree_thrown  = new TTree("T","particle production events");
  //Int_t particle_event;
  //tree_thrown->Branch("branch_thrown1",&particle_event,"particle_event/I");
  //Float_t *ntuple_branch;
  //tree_thrown->Branch("branch_thrown",&ntuple_branch,VarList);
  
  
  if(simul_key == 1) {
    ntuple_thrown = new TNtuple("ntuple_thrown","thrown particles",VarList);
    e_thrown = new TNtuple("e_thrown","thrown Electrons","Q2:W:Nu:vxec:vyec:vzec:vxe:vye:vze:Pex:Pey:Pez:evnt");
  }
  
  
//  TH1F *ht = new TH1F("ht","tdiff",1000,-15,15); 
  cout.width(4);
  input->Next();

  for (Int_t k = 0; k < nEntries; k++) 
    // for (Int_t k = 0; k < 100; k++) 
    {

      Int_t nRows = input->GetNRows("EVNT");
      //Int_t nRows = input->GetNRows("GSIM");
      
      Int_t GSIMrows = input->GetNRows("GSIM");
      Int_t to_fill = std::abs(GSIMrows-nRows); // need to fill the difference
      const char * tt = "C";

      //cout << k << endl;
      
      ////////////////////////
      // PION+ FILTER
      bool PionEvent = false;
      bool RareTH = false;
      bool RareREC = false;
      for ( Int_t i=1; i < nRows; i++ )
	{
	  TString categ = t->GetCategorization(i,tt);
	  if ( t->GetCategorization(0,tt) != "electron" )
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
            
      

      
      //if(nRows>0 && (t->GetCategorization(0,tt)) == "electron" && t -> Q2() > 1. && t -> W() > 2. && t -> Nu() / 5.015 < 0.85)
      if(nRows>0 && (t->GetCategorization(0,tt)) == "electron")  
	{
	  //Q2:W:Nu:vxec:vyec:vzec:vxe:vye:vze:Pex:Pey:Pez:event
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
	  
	  tElec->Fill(DataElec);
	  //std::cout<<"event: "<<input->GetCurrentEvent()<<std::endl;
	  //std::cout<<"got electron data"<<std::endl;
	  //Int_t NmbPion = 0;
	  if ( nRows == 1 )
	    {
	      for ( Int_t l = 1; l < GSIMrows; l++ )
		{
		  for ( Int_t ll = 0; ll < Nvar; ll++)
		    vars[ll] = 0;
		  vars[evntpos] = k;
		  ntuple->Fill(vars);
		}
	    }
	  for (Int_t i = 1; i < nRows; i++) 
	    {
	      
	      TString category = t->GetCategorization(i,tt);
	      //if (category == "high energy pion +" || category == "low energy pion +")
	      //is_pion=1;
	      //particle_event = nRows;
	      //std::cout<<"\tnr: "<<i<<std::endl; 
	      
	      //	      if (category == "gamma" || category == "pi-" || category == "high energy pion +" || category == "low energy pion +" || category == "s_electron" || category == "positron") 
	      //{

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
	      //vars[f] = TMath::Max(t->Etot(i),t->Ein(i)+t->Eout(i)); f++;
	      //vars[f] = TMath::Max(t->Etot(0),t->Ein(0)+t->Eout(0)); f++;
	      //vars[f] = t->Momentum(0); f++;
	      //vars[f] = t->TimeEC(0); f++;
	      //vars[f] = t->TimeSC(0); f++;
	      //vars[f] = t->PathEC(0); f++;
	      //vars[f] = t->PathSC(0); f++;
	      vars[f] = k; f++;
	      //vars[f] = t->Px(i); f++;
	      //vars[f] = t->Py(i); f++;
	      //vars[f] = t->Pz(i); f++;
	      //vars[f] = t->X(0); f++;
	      //vars[f] = t->Y(0); f++;
	      //vars[f] = t->Z(0); f++;
	      vert = t->GetCorrectedVert();
	      //vars[f] = vert->X();  f++;
	      //vars[f] = vert->Y();  f++;
	      //vars[f] = vert->Z();  f++;
	      //vars[f] = t->TimeEC(i); f++;
	      //vars[f] = t->XEC(i); f++;
	      //vars[f] = t->YEC(i); f++;
	      //vars[f] = t->ZEC(i); f++;
	      //vars[f] = t->Px(0); f++;
	      //vars[f] = t->Py(0); f++;
	      //vars[f] = t->Pz(0); f++;
	      
	      //vars[f] = t->Ein(i); f++;
	      //vars[f] = t->Eout(i); f++;
	      //vars[f] = t->Ein(0); f++;
	      //vars[f] = t->Eout(0); f++;
	      
	      vars[f] = ((category == "gamma")?22:
			  ((category == "pi-")?-211:
			   (( category == "high energy pion +" || category == "low energy pion +")?211:
			    ((category == "s_electron")?11:-11)
			    )
			   )
			 ); f++;
	      //vars[f] = t->Betta(i); f++;
	      //vars[f] = t->X(i); f++;
	      //vars[f] = t->Y(i); f++;
	      //vars[f] = t->Z(i); f++;
	      //vars[f] = is_pion; //f++;
	      //vars[f] = k;
	      ntuple->Fill(vars);
	      //is_pion = 0;
	      if ( i  == nRows-1 && GSIMrows > nRows )
		{
		  //cout << i << " " << nRows << endl;
		  //cout << "row : " << k << " " ;
		  //cout << "filling : " << GSIMrows  << " - " << nRows;
		  //cout << " = " << to_fill << endl;
		  for ( Int_t l = 0; l < to_fill; l++ )
		    {
		      for ( Int_t ll = 0; ll < Nvar; ll++)
			vars[ll] = 0;
		      vars[evntpos] = k;
		      ntuple->Fill(vars);
		    }
		}
	    }
	}
      else 
	{
	  for ( Int_t l = 0; l < 13; l++ )
	    DataElec[l] = 0;
	  tElec->Fill(DataElec);

	  for ( Int_t i = 1; i < std::max(GSIMrows,nRows); i++ )
	    {
	      for ( Int_t ll = 0; ll < Nvar; ll++ )
		vars[ll] = 0;
	      vars[evntpos] = k;
	      //vars[Nvar] = k;
	      ntuple->Fill(vars);
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
	  //      std::cout<<"got electron gsim"<<std::endl;
	  //Int_t NmbPion = 0;
	  //particle_event = input->GetNRows("GSIM");
	  //cout << particle_event << endl;
	  //ntuple_branch = new Float_t[particle_event];
	  for( Int_t i=1; i < GSIMrows; i++ )
	    {
	      //if(t -> Id(i,1)==1 || t -> Id(i,1)==9 || t -> Id(i,1)==8 ) //gamma: 1/22, pi0,+,-: 7/111,8/211,9 (Geant3/pdg)
	      {
		//if ( t-> Id(i,1) == 8 )
		//  {
		//    is_pion = 1;
		//  }
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
		//vars[f] = t->Momentum(i,1); f++;//TMath::Max(t->Etot(i),t->Ein(i)+t->Eout(i));;
		//vars[f] = TMath::Sqrt(t->Momentum(0,1)*t->Momentum(0,1)+kMe*kMe);  f++;//TMath::Max(t->Etot(0),t->Ein(0)+t->Eout(0));
		//vars[f] =t->Momentum(0,1); f++;
		//vars[f] = 0; f++;//t->TimeEC(0);
		//vars[f] = 0; f++;//t->TimeSC(0);
		//vars[f] = 0; f++;//t->PathEC(0);
		//vars[f] = 0; f++;//t->PathSC(0);
		vars[f] = k; f++;
		//vars[f] = t->Px(i,1); f++;
		//vars[f] = t->Py(i,1); f++;
		//vars[f] = t->Pz(i,1); f++;
		//vars[f] = t->X(0,1); f++;
		//vars[f] = t->Y(0,1); f++;
		//vars[f] = t->Z(0,1); f++;
		//vert = t->GetCorrectedVert();
		//vars[f] = t->X(0,1); f++;//vert->X(); 
		//vars[f] = t->Y(0,1); f++;//vert->Y(); 
		//vars[f] = t->Z(0,1); f++;//vert->Z(); 
		//vars[f] = 0; f++;//t->TimeEC(i);
		//vars[f] = 0; f++;//t->XEC(i);
		//vars[f] = 0; f++;//t->YEC(i);
		//vars[f] = 0; f++;//t->ZEC(i);
		//vars[f] = t->Px(0,1); f++;
		//vars[f] = t->Py(0,1); f++;
		//vars[f] = t->Pz(0,1); f++;
		
		//vars[f] = 0; f++;
		//vars[f] = 0; f++;
		//vars[f] = 0; f++;
		//vars[f] = 0; f++;
		vars[f] = t -> Id(i,1); f++;
		//vars[f] = t->Betta(i,1); f++;
		//vars[f] = t->X(i,1); f++;
		//vars[f] = t->Y(i,1); f++;
		//vars[f] = t->Z(i,1); f++;
		//vars[f] = is_pion;// f++;
		//vars[f] = k;
		ntuple_thrown->Fill(vars);
		//is_pion = 0;
		//ntuple_branch = ntuple_thrown->GetArgs();
		//ntuple_branch[i].Fill(vars);;
		//tree_thrown->Fill();
	      }	
	    }
	  //delete ntuple_branch;
	  if ( nRows > GSIMrows )
	    {
	      for ( Int_t i = 0; i < to_fill; i++ )
		{
		  for ( Int_t ll = 0; ll < Nvar; ll++ )
		    vars[ll] = 0;
		  vars[evntpos] = k;
		  ntuple_thrown->Fill(vars);
		}
	    }
	}
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
	      ntuple_thrown->Fill(vars);
	    }
	  
	  
	  // if there is pair creation between thrown and accepted
	  // if ( nRows > GSIMrows )
	  //   {
	  //     //if (k == 1858 )
	  //     //cout << "im (also?) on line 374!" << endl;
	  //     for ( Int_t i = 0; i < to_fill; i++ )
	  // 	{
	  // 	  for ( Int_t ll = 0; ll < Nvar; ll++ )
	  // 	    vars[ll] = 0;
	  // 	  vars[50] = k;
	  // 	  ntuple_thrown->Fill(vars);
	  // 	}
	  //   }
	}
      
      cout<<std::right<<float(k+1)/nEntries*100<<"%\r";
      cout.flush();
      //cout<<std::right<<float(k+1)/nEntries*100<<"%\n";
      input->Next();
    }
  
  output->Write();
  output->Close();
  cout << "Done." << endl;
  return 0;
}
