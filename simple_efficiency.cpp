//////////////////////////////////////////////////////////////////////////////////
//
//This code evaluatea simple efficiency for n_modules in n_stations, so pleas
//BEFORE using the code set n_mod to 6 or 12 and n_stat to 1 or 2.
//It will give some numbers for efficiency without traking or alignment and
//gives some idea on if there's some asynchronized module
//			TO RUN THE CODE:
//   g++ simple_efficiency.cpp -o _re `root-config --cflags --glibs` -std=c++14
//./_re /eos/experiment/mu-e/staging/daq/2023/decoded/commissioning/NAME_OF_THE_RUN
//
//////////////////////////////////////////////////////////////////////////////////


#include <string>
#include <iostream>
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLine.h"
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
      while ((file=(TSystemFile*)next()) and i<5) { if(i<0){i++; continue;}
//      while ((file=(TSystemFile*)next()) ) {
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





//THIS NEED TO BE SET BEFORE STARTING
int const n_mod=12;
int const n_stat=2;

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
        //chain.SetBranchAddress("BendCode", &stub_Bend);
        chain.SetBranchAddress("Bx", &stub_Bx);
        chain.SetBranchAddress("Link", &stub_Link);
        chain.SetBranchAddress("SuperID", &stub_SuperID);
        chain.SetBranchAddress("LocalX", &stub_LocalX);
        chain.SetBranchAddress("LocalY", &stub_LocalY);

        TH1F* histo_1stub= new TH1F("ones","when #stub=1, which module is?",6,0,6);
        TH1F* histo_nstubs_BX=new TH1F("histo_nstubs_BX", "nstubs per BX",50,0,50);
	//TH1F* histo_1stub= new TH1F("ones","when #stub=1, which module is?",6,0,6);

        std::vector<TH1D*> histo_nstubs_per_mod_BX(n_mod);
        for(int m=0; m<n_mod; m++)
        {	string name="h_stubs_permod_BX"+to_string(m);
                string title="nstubs per BX mod"+to_string(m);
	 histo_nstubs_per_mod_BX.at(m)=new TH1D(name.c_str(),title.c_str(),40,0,40);
	}

        std::vector<TH1D*> h_hist(n_mod);

	for(int m=0; m<n_mod; m++)
	{      	string name="h_hist"+to_string(m);
		string title="#stubs_mod==1; localX of module"+to_string(m);
		h_hist.at(m)=new TH1D(name.c_str(),title.c_str(),1000,-5.,5.);
	}

        std::vector<TH1D*> h_lost(n_stat);
        for(int s=0; s<n_stat; s++)
        {	string name="h_lost"+to_string(s);
                string title="When I have 5 stubs, which module is missing? station "+to_string(s);
	        h_lost.at(s)=new TH1D(name.c_str(),title.c_str(),6,0,6.);
	}

	std::vector<TH1D*> h_bx_prec(n_mod);
        std::vector<TH1D*> h_bx_next(n_mod);
        for(int m=0; m<n_mod; m++)
        {	string name1="h_bx_next"+to_string(m);
                string title1="Bx with stub (event->bx()-Bx) module_link"+to_string(m);
                h_bx_next.at(m)=new TH1D(name1.c_str(),title1.c_str(),20,-10,10);
		string name2="h_bx_prec"+to_string(m);
                string title2="Bx with stub (event->bx()-Bx) module_link"+to_string(m);
                h_bx_prec.at(m)=new TH1D(name2.c_str(),title2.c_str(),20,-10,10);

        }

TH1D* h_BX = new TH1D("BX","BX value per event", 20,3550.,3570.);

double c_mod=0.;

double stubsix_st0=0;
double stubsix_st1=0;

double stub5_st0[6]={0.};
double stub5_st1[6]={0.};

double mod_st0_given_st1[6]={0.};
double mod_st1_given_st0[6]={0.};


for(Long64_t i = 0; i < 300000; i++) {//300000; i++) {//chain.GetEntries(); i++) {
                chain.GetEntry(i);
                if(i%1000 == 0) cout<<"Entry "<<i<<endl;

// From strip number to cm in local reference frame



	              	auto local_coo=[&](int mod,int stub,std::array<std::vector<float>,n_mod> locX_strip_units){
				 int a=1;
				 if(mod % 2==0) a=-1;
				 return a*( + (locX_strip_units.at(mod).at(stub) - 1016./2.)*0.009 +0.009/2.);};


// X and Y rotation in local UV frame
                auto newU=[](double angle, auto X, auto Y){return cos(angle)*X - sin(angle)*Y;};
                auto newV=[](double angle, auto X, auto Y){return sin(angle)*X + cos(angle)*Y;};


int found_st0[6]={0};
int found_st1[6]={0};

        std::vector<uint16_t> nstubs_per_mod_BX_st0(6,0);
        std::vector<uint16_t> nstubs_per_mod_BX_st1(6,0);
        std::vector<int> module_st0; module_st0.reserve(40);
        std::vector<int> module_st1; module_st1.reserve(40);
        std::array<std::vector<float>,n_mod> localX;
        int n_stubs_tot=stub_Link->size();
	std::array<int,n_stat> n_stubs;


	for(int s=0;s<n_stat;s++){n_stubs.at(s)=0.;}

	histo_nstubs_BX->Fill(n_stubs_tot);

	if(stub_Link->size()==1) histo_1stub->Fill(stub_Link->at(0));

cout << "stub at module : " << stub_Link->size() << endl;

        for(int j = 0; j<stub_Link->size(); j++)
        {
		cout << "link " << stub_Link->at(j) << " ,";
		int link= stub_Link->at(j);

         localX.at(link).push_back(stub_LocalX->at(j));

         if(link<6){nstubs_per_mod_BX_st0.at(link)+=1;
		    module_st0.push_back(link);
		    n_stubs.at(0)+=1;}
	 else{nstubs_per_mod_BX_st1.at(link-6)+=1;
		cout << " link-6 " << link-6 << endl;
	      module_st1.push_back(link);
	      n_stubs.at(1)+=1;}
        }

	cout << " " << endl;

//-------STATION 0--------

int six_st0=0;

   	for(int m=0; m<n_mod/n_stat; m++){
	cout << "station 0, module " << m << "  n_stub: " << nstubs_per_mod_BX_st0.at(m) << " " << endl;
		if(nstubs_per_mod_BX_st0.at(m)==1) h_hist[m]->Fill(local_coo(m,0,localX));
	}
//efficiency station 0

        cout << "n_stubs staton 0: " << n_stubs.at(0) <<endl;

        sort(module_st0.begin(), module_st0.end());
if(adjacent_find(module_st0.begin(), module_st0.end())==module_st0.end())
 {
     for(int m=0; m<n_mod/n_stat; m++){
	nstubs_per_mod_BX_st0.at(m)==1? found_st0[m]=1 : found_st0[m]=0;
	cout << "module " << m << ", found_st0 = " << found_st0[m] << endl;
     }
if(n_stubs.at(0)==5){auto it=find(begin(found_st0), end(found_st0),0);
		     auto pos=std::distance(begin(found_st0), it);
		     h_lost.at(0)->Fill(pos);
		     cout << "manca modulo " << std::distance(begin(found_st0), it) << endl;

int Bx=stub_Bx->at(0);
int found=0;
chain.GetEntry(i-1);
        vector<int> nstubs_per_mod_BX1(6,0);
        vector<vector<float>> stubs_mod1(40);
        double ns1=0.;
        for(int i = 0; i<stub_Link->size(); i++)
        {
         int link= stub_Link->at(i);
         if(link<6){ns1++; nstubs_per_mod_BX1.at(link)+=1;
	            stubs_mod1.at(link).push_back(stub_LocalX->at(i));}
        }
if(nstubs_per_mod_BX1.at(pos)==1 and ns1==1){
h_bx_prec[pos]->Fill(stub_Bx->at(0)-Bx);
	}
chain.GetEntry(i+1);
        vector<int> nstubs_per_mod_BX2(6,0);
        vector<vector<float>> stubs_mod2(40);
	double ns2=0.;
        for(int i = 0; i<stub_Link->size(); i++)
        {
         int link= stub_Link->at(i);
         if(link<6){ns2++; nstubs_per_mod_BX2.at(link)+=1;
         stubs_mod2.at(link).push_back(stub_LocalX->at(i));}
        }
if(nstubs_per_mod_BX2.at(pos)==1 and ns2==1){
h_bx_next[pos]->Fill(stub_Bx->at(0)-Bx);
        }
chain.GetEntry(i);
}
else if(n_stubs.at(0)==6){stubsix_st0++; six_st0=1;}

     for(int m=0; m<6; m++){
	int count=0;
      for(int z=0; z<6; z++){
		//counting if in the modules different from module m (of which I want the efficiency) I have 1 hit
		 if(z!=m and nstubs_per_mod_BX_st0.at(z)==1) count++;
		}
	 if(count==5) stub5_st0[m]++;
	}

} //end of station 0 efficiency


//-------STATION 1--------


 if(n_stat==2){

int six_st1=0;

   	for(int m=0; m<n_mod/n_stat; m++){
        cout << "station 1, module " << m << "  n_stub: " << nstubs_per_mod_BX_st1.at(m) << endl;
                //if(nstubs_per_mod_BX_st1.at(m)==1) h_hist[m+6]->Fill(local_coo(m+6,0,localX));
        }

//efficiency station 1

        cout << "n_stubs staton 1: " << n_stubs.at(1) <<endl;

        sort(module_st1.begin(), module_st1.end());
if(adjacent_find(module_st1.begin(), module_st1.end())==module_st1.end())
 {
     for(int m=0; m<6; m++){
        nstubs_per_mod_BX_st1.at(m)==1? found_st1[m]=1 : found_st1[m]=0;
        cout << "module " << m << ", found_st1 = " << found_st1[m] << endl;
     }
if(n_stubs.at(1)==5){auto it=find(begin(found_st1), end(found_st1),0);
                     auto pos=std::distance(begin(found_st1), it);
                     h_lost.at(1)->Fill(pos);
                     cout << "manca modulo " << std::distance(begin(found_st1), it) << endl;
 int Bx=stub_Bx->at(0);
 int found=0;
//search roughly a stub in other BX
chain.GetEntry(i-1);
        vector<int> nstubs_per_mod_BX1(6,0);
        vector<vector<float>> stubs_mod1(40);
        double ns1=0.;
         for(int i = 0; i<stub_Link->size(); i++)
         {
          int link= stub_Link->at(i);
         if(link>5){ns1++; nstubs_per_mod_BX1.at(link-6)+=1;
          stubs_mod1.at(link-6).push_back(stub_LocalX->at(i));}
         }
	if(nstubs_per_mod_BX1.at(pos)==1 and ns1==1){
	h_bx_prec[pos+6]->Fill(stub_Bx->at(0)-Bx);
        }
chain.GetEntry(i+1);
        vector<int> nstubs_per_mod_BX2(6,0);
        vector<vector<float>> stubs_mod2(40);
	double ns2=0.;
         for(int i = 0; i<stub_Link->size(); i++)
         {
          int link= stub_Link->at(i);
          if(link>5){ns2++; nstubs_per_mod_BX2.at(link-6)+=1;
          stubs_mod2.at(link-6).push_back(stub_LocalX->at(i));}
         }
	if(nstubs_per_mod_BX2.at(pos)==1 and ns2==1){
	h_bx_next[pos+6]->Fill(stub_Bx->at(0)-Bx);
        }
chain.GetEntry(i);
 }
else if(n_stubs.at(1)==6){stubsix_st1++; six_st1=1;}

     for(int m=0; m<6; m++){
        int count=0;
      for(int z=0; z<6; z++){
                //counting if in the modules different from module m (of which I want the efficiency) I have 1 hit
                 if(z!=m and nstubs_per_mod_BX_st1.at(z)==1) count++;
                }
         if(count==5) stub5_st1[m]++;
        }

      } //end of station 1 efficiency




// EVALUATE EFF. OF STATION 1 WHEN STATION 0 HAS 6 STUBS

if(six_st0==1){

     for(int m=0; m<6; m++){
        if(nstubs_per_mod_BX_st1.at(m)==1) mod_st1_given_st0[m]++;
	}

}


// EVALUATE EFF. OF STATION 0 WHEN STATION 1 HAS 6 STUBS

if(six_st1==1){

     for(int m=0; m<6; m++){
        if(nstubs_per_mod_BX_st0.at(m)==1) mod_st0_given_st1[m]++;
        }

}

if(n_stubs.at(0)>=4 and n_stubs.at(1)>=4) c_mod++;



   }//end of n_stat==2 if

 cout << " fine " << endl;
nstubs_per_mod_BX_st0.clear();
nstubs_per_mod_BX_st1.clear();
module_st0.clear();
module_st1.clear();
for(int l=0; l<localX.size();l++) localX.at(l).clear();
} //end of general for

	for(int m=0; m<6; m++){cout << "Efficiency module " << m << " station 0 is " << (stubsix_st0/stub5_st0[m]) << endl;}
        if(n_stat==2){
		for(int m=0; m<6; m++){cout << "Efficiency module " << m << " station 1 is " << (stubsix_st1/stub5_st1[m]) << endl;}
                for(int m=0; m<6; m++){cout << "Given 6 stubs in station 0, efficiency module " << m << " of station 1 is " << (mod_st1_given_st0[m]/stubsix_st0) << endl;}
                for(int m=0; m<6; m++){cout << "Given 6 stubs in station 1, efficiency module " << m << " of station 0 is " << (mod_st0_given_st1[m]/stubsix_st1) << endl;}
	}

cout << "risp st0 " << c_mod/stubsix_st0 <<endl;
cout << "risp st1 " << c_mod/stubsix_st1 <<endl;


TCanvas n ("canvas_stubs_BX","stubs per BX", 700, 500);
  histo_1stub->Draw();
  n.SaveAs("histo_1.pdf");

TCanvas n1("n1","n1",1000,1000);
histo_nstubs_BX->Draw();
n1.SaveAs("nstubstot_perBX.pdf");

TCanvas n3("n3","n3",1000,1000);
n3.Divide(2,n_mod/2);
for(int m=0; m<n_mod; m++){
n3.cd(m+1);
h_hist[m]->Draw();}
n3.SaveAs("local.pdf");

TCanvas n2("n2","n2",1000,1000);
n2.Divide(1,n_stat);
for(int s=0; s<n_stat; s++){
n2.cd(s+1);
h_lost[s]->Draw();}
n2.SaveAs("lost.pdf");

TCanvas n4("n4","n4",1000,1000);
n4.Divide(2,n_mod/2);
for(int m=0; m<n_mod; m++){
n4.cd(m+1);
h_bx_prec[m]->SetLineColor(kRed);
if(h_bx_next[m]->GetEntries()>h_bx_prec[m]->GetEntries()){
h_bx_next[m]->Draw();
h_bx_prec[m]->Draw("same");}
else{
h_bx_prec[m]->Draw();
h_bx_next[m]->Draw("same");
 }
}
n4.SaveAs("deltaBX.pdf");

TCanvas n5("n5","n5",1000,1000);
h_BX->Draw();
n5.SaveAs("BX_check.pdf");

}
