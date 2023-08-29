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
TChain chain("Cereal"); // chain event trees
TSystemDirectory dir(dirname, dirname);
   TList *files = dir.GetListOfFiles();
   if (files) {
int i=0;
      TSystemFile *file;
      TString fname;
      TIter next(files);
      while ((file=(TSystemFile*)next()) and i<11) { if(i<10){i++; continue;}
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

    std::vector<uint16_t>* stub_Bx = new std::vector<uint16_t>;
    std::vector<uint16_t>* stub_Link = new std::vector<uint16_t>;
    std::vector<uint32_t>* stub_SuperID = new std::vector<uint32_t>;
    std::vector<float>* stub_Bend = new std::vector<float>;
    std::vector<float>* stub_LocalX = new std::vector<float>;
    std::vector<float>* stub_LocalY = new std::vector<float>;
    std::vector<float>* stub_Z = new std::vector<float>;

        chain.SetBranchAddress("Z", &stub_Z);
        chain.SetBranchAddress("Bend", &stub_Bend);
        chain.SetBranchAddress("Bx", &stub_Bx);
        chain.SetBranchAddress("Link", &stub_Link);
        chain.SetBranchAddress("SuperID", &stub_SuperID);
        chain.SetBranchAddress("LocalX", &stub_LocalX);
        chain.SetBranchAddress("LocalY", &stub_LocalY);


        TH1F* histo_nstubs_BX=new TH1F("histo_nstubs_BX", "nstubs per BX",50,0,50);

        std::vector<TH1F*> histo_nstubs_per_mod_BX(6);
        histo_nstubs_per_mod_BX.at(0)=new TH1F("histo_nstubs_mod0_BX", "nstubs per BX mod 0",40,0,40);
        histo_nstubs_per_mod_BX.at(1)=new TH1F("histo_nstubs_mod1_BX", "nstubs per BX mod 1",40,0,40);
        histo_nstubs_per_mod_BX.at(2)=new TH1F("histo_nstubs_mod2_BX", "nstubs per BX mod 2",40,0,40);
        histo_nstubs_per_mod_BX.at(3)=new TH1F("histo_nstubs_mod3_BX", "nstubs per BX mod 3",40,0,40);
        histo_nstubs_per_mod_BX.at(4)=new TH1F("histo_nstubs_mod4_BX", "nstubs per BX mod 4",40,0,40);
        histo_nstubs_per_mod_BX.at(5)=new TH1F("histo_nstubs_mod5_BX", "nstubs per BX mod 5",40,0,40);

        TH1D* h_StartX=new TH1D("h_StartX","X mod0",120,-6,6);
        TH1D* h_StartY=new TH1D("h_StartY","Y mod1",120,-6,6);
        TH1D* h_divX=new TH1D("h_divX", "h_divX",300,-0.015,0.015);
        TH1D* h_divY=new TH1D("h_divY", "h_divY",300,-0.015,0.015); 

        std::vector<TH2F*> corr_mod;
          corr_mod.push_back(new TH2F("corr_u3", "expected UV when there are 2 stubs", 100,-5,5,100,-5,5));

double stub2=0;
double stub3=0;
double stub4=0;


//double Xoffsets[6]={-0.070300968495475138, 0.043476838127823117, 0.081218458321752882, 0.068428399849415442, -0.079017060072807038, -0.051013420927518687};
//double Yoffsets[6]={-0.17741353189702608, 0.0079311725718771414, -0.080868219919505421, 0.068374266162434874, -0.05630130387604667, 0.0046751032833262096};

double Xoffsets[6]={0.,0.,0.,0.,0.,0.};
double Yoffsets[6]={0.,0.,0.,0.,0.,0.};

//double tilt[6]={0.233,-0.233,0.,0.,0.233,-0.233};
//double posZ[6]={18.0218, 21.8693, 55.4700-0.2, 56.7305-0.2, 89.9218, 93.7693};

double tilt[6]={0.,0.,0.,0.,0.,0.};
double posZ[6]={0.,0.,0.,0.,0.,0.};

for(Long64_t i = 0; i < chain.GetEntries(); i++) {
                chain.GetEntry(i);
                if(i%1000 == 0) cout<<"Entry "<<i<<endl;


        vector<uint16_t> nstubs_per_mod_BX(6,0);
        vector<vector<float>> stubs_mod0(40);

// U and V from strip number to cm in local reference frame
                auto absU_=[&](int mod,int stub,vector<vector<float>> stubs_mod){ double posStripU = stubs_mod.at(mod).at(stub)/2. -1;
                                      return (- (posStripU - 1016/2.)*0.009 -0.009/2.);};
                auto absV_=[&](int mod,int stub,vector<vector<float>> stubs_mod){ double posStripV = stubs_mod.at(mod).at(stub)/2. -1;
                                      return ( + (posStripV - 1016/2.)*0.009 +0.009/2.);};
// X and Y from strip number to cm in local reference frame
                auto absX_=[&](int mod,int stub,vector<vector<float>> stubs_mod){ double posStripU = stubs_mod.at(mod).at(stub)/2. -1;
                                      return (cos(0.233)*(- (posStripU - 1016/2.)*0.009 -0.009/2.));};
                auto absY_=[&](int mod,int stub,vector<vector<float>> stubs_mod){ double posStripV = stubs_mod.at(mod).at(stub)/2. -1;
                                      return (cos(-0.233)*( + (posStripV - 1016/2.)*0.009 +0.009/2.));};
// Z position accounting for tilt
                auto absZ04_=[&](int mod,auto pos, double off){return - (pos/cos(0.233)+off)*sin(tilt[mod]) + posZ[mod] ;};
                auto absZ15_=[&](int mod,auto pos, double off){return - (pos/cos(-0.233)+off)*sin(tilt[mod]) + posZ[mod] ;};

// U and V rotation in global XY frame
                auto newX=[](double angle, auto U, auto V){return cos(angle)*U + sin(angle)*V;};
                auto newY=[](double angle, auto U, auto V){return -sin(angle)*U + cos(angle)*V;};

// X and Y rotation in local UV frame
                auto newU=[](double angle, auto X, auto Y){return cos(angle)*X - sin(angle)*Y;};
                auto newV=[](double angle, auto X, auto Y){return sin(angle)*X + cos(angle)*Y;};

                auto l_mXZ = [&](int posX0, int posX4){return (posX4-posX0)/(absZ04_(4,posX4,0.)-absZ04_(0,posX0,0.));};
                auto l_mYZ = [&](int posY1, int posY5){return (posY5-posY1)/(absZ15_(5,posY5,0.)-absZ15_(1,posY1,0.));};


        for(int j = 0; j<stub_Link->size(); j++)
        {
         int link= stub_Link->at(j);
         nstubs_per_mod_BX.at(link)+=1;
         stubs_mod0.at(link).push_back(stub_LocalX->at(j));
        }

histo_nstubs_BX->Fill(stub_Link->size());
for(int m=0; m<6; m++){ histo_nstubs_per_mod_BX.at(m)->Fill(nstubs_per_mod_BX.at(m));}

if(nstubs_per_mod_BX.at(0)!=0 and nstubs_per_mod_BX.at(1)!=0 and nstubs_per_mod_BX.at(4)!=0 and nstubs_per_mod_BX.at(5)!=0){
double pos_X0=absX_(0,0,stubs_mod0)+Xoffsets[0];
double pos_Y1=absY_(1,0,stubs_mod0)+Yoffsets[1];
corr_mod.at(0)->Fill(pos_X0,pos_Y1);

double pos_X4=absX_(4,0,stubs_mod0)+Xoffsets[4];
double pos_Y5=absY_(5,0,stubs_mod0)+Yoffsets[5];

h_StartX->Fill(pos_X0);
h_StartY->Fill(pos_Y1);
double mXZ=(pos_X4-pos_X0)/(89.9218-18.0218);
double mYZ=(pos_Y1-pos_Y5)/(93.7693-21.8693);
h_divX->Fill(mXZ);
h_divY->Fill(mYZ);
}

stubs_mod0.clear();

} //end of general for

TCanvas n0 ("canvas_stubs_BX","stubs per BX", 700, 500);
n0.Divide(1,1);
 n0.cd(1);
  histo_nstubs_BX->Draw();
  n0.SaveAs("histo_nstubs_BX.pdf");

TCanvas n1 ("canvas_stubs_BX_mod","stubs per mod per BX", 700, 500);
  n1.Divide(2,3);
  n1.cd(1);
  histo_nstubs_per_mod_BX.at(0)->Draw();
        for(int c=1; c<6; c++)
        {n1.cd(c+1);
         histo_nstubs_per_mod_BX.at(c)->Draw(); 
        }
  n1.SaveAs("histo_nstubs_per_mod_BX.pdf");

TCanvas n2("n2","n2",1000,1000);
n2.Divide(2,3);
n2.cd(1);
h_StartX->Draw( );
n2.cd(2);
h_StartY->Draw( );
n2.cd(3);
h_divX->Draw( );
n2.cd(4);
h_divY->Draw( );
n2.cd(5);
corr_mod.at(0)->Draw("COLZ");
n2.SaveAs("XY.pdf");
}

