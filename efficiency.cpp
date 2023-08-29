#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TSystemDirectory.h"
#include <TStyle.h>
#include "MUonERealDataEvent.h"

using namespace std;
int main(int argc, char* argv[]){

//      TFile *inputfile = new TFile("dataReconstruction.root");
  if (argc != 1) {
    //    cerr << "Usage : "<< argv[0] << endl;
    cerr << "argc = " << argc << endl;
    cerr << "argv[0] = "<< argv[0] << endl;
    cerr << "argv[1] = "<< argv[1] << endl;
    //    exit(100);
  }

 bool _debug_ = true;

const char *dirname=argv[1];
TString path_s(dirname);

const char *ext=".root";
TChain chain("cbmsim"); // chain event trees
TSystemDirectory dir(dirname, dirname);
   TList *files = dir.GetListOfFiles();
   if (files) {
int i=0;
      TSystemFile *file;
      TString fname;
      TIter next(files);
      while ((file=(TSystemFile*)next()) and i<15) { if(i<10){i++; continue;}
         fname = path_s+file->GetName();
         if (!file->IsDirectory() && fname.EndsWith(ext)) {

          cout << fname.Data() << endl;
        chain.Add(fname.Data());
        i++;
         }
      }
   }
Long64_t nevtot = chain.GetEntries();
  cerr << "Total number of input events  = " << nevtot << endl;

  // Set it to kTRUE if you do not run interactively
  gROOT->SetBatch(kTRUE); 
  //gROOT->SetBatch(kFALSE); 

	TRint* app = new TRint("Root Application", &argc, argv);
        MUonERealDataEvent* event = new MUonERealDataEvent();

	chain.SetBranchAddress("Event", &event);

	TH1F* histo_nstubs_BX=new TH1F("histo_nstubs_BX", "nstubs per BX",50,0,50);

	std::vector<TH1F*> histo_nstubs_per_mod_BX(6);
	histo_nstubs_per_mod_BX.at(0)=new TH1F("histo_nstubs_mod0_BX", "nstubs per BX mod 0",40,0,40);
        histo_nstubs_per_mod_BX.at(1)=new TH1F("histo_nstubs_mod1_BX", "nstubs per BX mod 1",40,0,40);
        histo_nstubs_per_mod_BX.at(2)=new TH1F("histo_nstubs_mod2_BX", "nstubs per BX mod 2",40,0,40);
        histo_nstubs_per_mod_BX.at(3)=new TH1F("histo_nstubs_mod3_BX", "nstubs per BX mod 3",40,0,40);
        histo_nstubs_per_mod_BX.at(4)=new TH1F("histo_nstubs_mod4_BX", "nstubs per BX mod 4",40,0,40);
        histo_nstubs_per_mod_BX.at(5)=new TH1F("histo_nstubs_mod5_BX", "nstubs per BX mod 5",40,0,40);

	std::vector<TH2F*> corr_mod;
	  corr_mod.push_back(new TH2F("corr_01", "corr mod0 vs mod1", 10000,-5,5,10000,-5,5));
	  corr_mod.push_back(new TH2F("corr_45", "corr mod4 vs mod5", 10000,-5,5,10000,-5,5));
	  corr_mod.push_back(new TH2F("corr_04", "corr mod0 vs mod4", 10000,-1,1,10000,-1,1));
	  corr_mod.push_back(new TH2F("corr_15", "corr mod1 vs mod5", 10000,-1,1,10000,-1,1));

	TH1F* modY=new TH1F("Ymod1", "Y without X", 10000,-5,5);
        TH1F* modX=new TH1F("Xmod0", "X without Y", 10000,-5,5);
        TH1F* modwithY=new TH1F("Ymod1w", "Y with X", 10000,-5,5);
        TH1F* modwithX=new TH1F("Xmod0w", "X with Y", 10000,-5,5);
        TH1F* mod4=new TH1F("X4", "X mod 4 no X0", 10000,-5,5);


double ynox,yx,ynox_all,y_all,y=0;

for(Long64_t i = 0; i < chain.GetEntries(); i++) {
		chain.GetEntry(i);
		if(i%1000 == 0) cout<<"Entry "<<i<<endl;


	vector<MUonERealDataStub> stub= event->stubs();
        vector<int> nstubs_per_mod_BX(10,0);
	vector<vector<MUonERealDataStub>> stubs_mod(40);

	for(int i = 0; i<stub.size(); i++)
	{
	 int link= stub.at(i).link();
	 nstubs_per_mod_BX.at(link)+=1;
	 stubs_mod.at(link).push_back(stub.at(i));
	}

	if(nstubs_per_mod_BX.at(1)==1) y++;
        if(nstubs_per_mod_BX.at(0)==1 and nstubs_per_mod_BX.at(1)==1) yx++;
	if(nstubs_per_mod_BX.at(0)==0 and nstubs_per_mod_BX.at(1)==1) ynox++;
	if(nstubs_per_mod_BX.at(1)==1 and nstubs_per_mod_BX.at(4)==1 and nstubs_per_mod_BX.at(5)==1) y_all++;
        if(nstubs_per_mod_BX.at(0)==0 and nstubs_per_mod_BX.at(1)==1 and nstubs_per_mod_BX.at(4)==1 and nstubs_per_mod_BX.at(5)==1) 
	{ynox_all++;modY->Fill(stubs_mod.at(1).at(0).absY());mod4->Fill(stubs_mod.at(4).at(0).absX());
	 corr_mod.at(1)->Fill(stubs_mod.at(4).at(0).absX(),stubs_mod.at(5).at(0).absY());}
	if(nstubs_per_mod_BX.at(0)==1 and nstubs_per_mod_BX.at(1)==0 and nstubs_per_mod_BX.at(4)==1 and nstubs_per_mod_BX.at(5)==1) modX->Fill(stubs_mod.at(0).at(0).absX());

	if(nstubs_per_mod_BX.at(0)==1 and nstubs_per_mod_BX.at(1)==1 and nstubs_per_mod_BX.at(4)==1 and nstubs_per_mod_BX.at(5)==1){modwithY->Fill(stubs_mod.at(1).at(0).absY());
																	modwithX->Fill(stubs_mod.at(0).at(0).absX());}

	if(nstubs_per_mod_BX.at(0)==1 and nstubs_per_mod_BX.at(1)==1 and nstubs_per_mod_BX.at(4)==1 and nstubs_per_mod_BX.at(5)==1)
	{
	corr_mod.at(0)->Fill(stubs_mod.at(0).at(0).absX(),stubs_mod.at(4).at(0).absX());
	}

} //end of general for


cout <<"numero di eventi Y0 senza X0 " << ynox << " su " << y << " eventi y " << endl;
cout <<"numero di eventi Y0 con Y1 e X1 senza X0 " << ynox_all << " su " << y_all << " eventi y_all " << endl;
cout <<"numero di eventi Y0 con X0 " << yx << " su " << y << " eventi y " << endl;


TCanvas n2 ("X_Y", "X and Y per events with 6 hits", 700,500);
  n2.Divide(2,2);
  n2.cd(1);
  corr_mod.at(1)->Draw("COLZ");
  n2.cd(2);
  mod4->Draw();
  n2.cd(3);
modY->SetLineColor(kRed);
  modwithY->Draw();
  modY->Draw("same");
  n2.cd(4);
modX->SetLineColor(kRed);
  modwithX->Draw();
  modX->Draw("same");
  n2.SaveAs("impact_coo_6hits.pdf");
  n2.SaveAs("impact_coo_6hits.root");

delete event;
for(int d=0; d<4; d++){delete histo_nstubs_per_mod_BX[d];}
delete histo_nstubs_BX;

}


