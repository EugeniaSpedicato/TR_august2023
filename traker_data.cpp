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
      while ((file=(TSystemFile*)next())) {
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

int n_mod=6;


	TRint* app = new TRint("Root Application", &argc, argv);

    std::vector<uint16_t>* stub_Bx = new std::vector<uint16_t>;
    std::vector<uint16_t>* stub_Link = new std::vector<uint16_t>;
    std::vector<uint32_t>* stub_SuperID = new std::vector<uint32_t>;
    std::vector<float>* stub_Bend = new std::vector<float>;
    std::vector<float>* stub_BendCode = new std::vector<float>;
    std::vector<float>* stub_LocalX = new std::vector<float>;
    std::vector<float>* stub_LocalY = new std::vector<float>;
    std::vector<float>* stub_Z = new std::vector<float>;

        chain.SetBranchAddress("Z", &stub_Z);
        chain.SetBranchAddress("Bend", &stub_Bend);
        chain.SetBranchAddress("BendCode", &stub_BendCode);
        chain.SetBranchAddress("Bx", &stub_Bx);
        chain.SetBranchAddress("Link", &stub_Link);
        chain.SetBranchAddress("SuperID", &stub_SuperID);
        chain.SetBranchAddress("LocalX", &stub_LocalX);
        chain.SetBranchAddress("LocalY", &stub_LocalY);

	std::vector<uint16_t> Nstubs4_Bx;
        std::vector<uint16_t> Nstubs4_SuperID;
        std::vector<uint16_t> Nstubs5_Bx;
        std::vector<uint16_t> Nstubs5_SuperID;
        std::vector<uint16_t> Nstubs6_Bx;
        std::vector<uint16_t> Nstubs6_SuperID;

for(Long64_t i = 0; i < chain.GetEntries(); i++) {
		chain.GetEntry(i);
		if(i%1000 == 0) cout<<"Entry "<<i<<endl;


        std::vector<int> nstubs_per_mod_Bx(6,0);
	std::vector<int> module; module.reserve(40);
	int n_stubs=stub_Link->size();

	int Bx = stub_Bx->at(0);
	int SuperID = stub_SuperID->at(0);

	for(int j = 0; j<n_stubs; j++){nstubs_per_mod_Bx.at(stub_Link->at(j))+=1; module.push_back(stub_Link->at(j));}

	if(n_stubs==4){ sort(module.begin(), module.end());
			if(adjacent_find(module.begin(), module.end())==module.end()){
											Nstubs4_Bx.push_back(Bx);
											Nstubs4_SuperID.push_back(SuperID);
			    }
		      }

        if(n_stubs==5){ sort(module.begin(), module.end());
                       	     if(adjacent_find(module.begin(), module.end())==module.end()){

                                                                                        Nstubs5_Bx.push_back(Bx);
                                                                                        Nstubs5_SuperID.push_back(SuperID);
                            }
		      }

	if(n_stubs==6){ sort(module.begin(), module.end());
                       	     if(adjacent_find(module.begin(), module.end())==module.end()){
                                                                                        Nstubs6_Bx.push_back(Bx);
                                                                                        Nstubs6_SuperID.push_back(SuperID);
                            }
                      }


} //end of general for


}
