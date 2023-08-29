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
      while ((file=(TSystemFile*)next()) and i<12) { if(i<10){i++; continue;}
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


        TH1F* histo_nstubs_BX=new TH1F("histo_nstubs_BX", "nstubs per BX",50,0,50);

        std::vector<TH1F*> histo_nstubs_per_mod_BX(6);
        histo_nstubs_per_mod_BX.at(0)=new TH1F("histo_nstubs_mod0_BX", "nstubs per BX mod 0",40,0,40);
        histo_nstubs_per_mod_BX.at(1)=new TH1F("histo_nstubs_mod1_BX", "nstubs per BX mod 1",40,0,40);
        histo_nstubs_per_mod_BX.at(2)=new TH1F("histo_nstubs_mod2_BX", "nstubs per BX mod 2",40,0,40);
        histo_nstubs_per_mod_BX.at(3)=new TH1F("histo_nstubs_mod3_BX", "nstubs per BX mod 3",40,0,40);
        histo_nstubs_per_mod_BX.at(4)=new TH1F("histo_nstubs_mod4_BX", "nstubs per BX mod 4",40,0,40);
        histo_nstubs_per_mod_BX.at(5)=new TH1F("histo_nstubs_mod5_BX", "nstubs per BX mod 5",40,0,40);

        std::vector<TH2F*> corr_mod;
          corr_mod.push_back(new TH2F("corr_u1", "expected UV plane when there is a stub in U plane", 1000,-5,5,1000,-5,5));
          corr_mod.push_back(new TH2F("corr_v1", "expected UV plane when there is a stub in V plane", 1000,-5,5,1000,-5,5));
          corr_mod.push_back(new TH2F("corr_u2", "expected UV plane when there is NOT a stub in U plane", 1000,-5,5,1000,-5,5));
          corr_mod.push_back(new TH2F("corr_v2", "expected UV plane when there is NOT a stub in V plane", 1000,-5,5,1000,-5,5));
          corr_mod.push_back(new TH2F("corr_u3", "expected UV when there are 2 stubs", 1000,-5,5,1000,-5,5));
          corr_mod.push_back(new TH2F("corr_v3", "expected UV when there are NOT 2 stubs", 1000,-5,5,1000,-5,5));

        TH1F* h_posU_expected=new TH1F("exp_u", "expected U (cm)", 1111,-5,5);
        TH1F* h_posU=new TH1F("h_u", " U position (cm)", 1111,-5,5);
        TH1F* h_posU_expected_no=new TH1F("exp_v", "expected V (cm) no hit", 1111,-5,5);
        TH1F* h_posV_expected=new TH1F("exp_v", "expected V (cm)", 1111,-5,5);
        TH1F* h_posV_expected_no=new TH1F("exp_v", "expected V (cm) no hit", 1111,-5,5);
        TH1F* h_posV=new TH1F("h_v", " V position (cm)", 1111,-5,5); 

        TH1F* residual_X=new TH1F("residual_V1", "residual_U Bx-1",222,-0.1,0.1);
        TH1F* residual_Y=new TH1F("residual_V2", "residual_U Bx+1",222,-0.1,0.1);
        TH1F* residual_Z=new TH1F("residual_V3", "residual_U Bx",222,-0.1,0.1);

        TH1F* residual_X_v=new TH1F("residual_V1", "residual_V Bx-1",222,-0.1,0.1);
        TH1F* residual_Y_v=new TH1F("residual_V2", "residual_V Bx+1",222,-0.1,0.1);
        TH1F* residual_Z_v=new TH1F("residual_V3", "residual_V Bx",222,-0.1,0.1);

	TH1F* h_bx0=new TH1F("bx0","Bx with stub (event->bx()-Bx) modU",30,-15,15);
        TH1F* h_bx1=new TH1F("bx1","Bx with stub (event->bx()-Bx) modU",30,-15,15);

        TH1F* h_bx0_v=new TH1F("bx0_v","Bx with stub (event->bx()-Bx) modV",30,-15,15);
        TH1F* h_bx1_v=new TH1F("bx1_v","Bx with stub (event->bx()-Bx) modV",30,-15,15);

double stubU1=0; double stubU2=0; double stubU3=0;
double stubV1=0; double stubV2=0; double stubV3=0;
double stub4=0;
int found=0; double found_v=0; double two=0;

double Xoffsets[6]={-0.070300968495475138, 0.043476838127823117, 0.081218458321752882, 0.068428399849415442, -0.079017060072807038, -0.051013420927518687};
double Yoffsets[6]={-0.17741353189702608, 0.0079311725718771414, -0.080868219919505421, 0.068374266162434874, -0.05630130387604667, 0.0046751032833262096};

//double Xoffsets[6]={0.,0.,0.,0.,0.,0.};
//double Yoffsets[6]={0.,0.,0.,0.,0.,0.};

double tilt[6]={0.233,-0.233,0.,0.,0.233,-0.233};
double posZ[6]={18.0218, 21.8693, 55.4700-0.2, 56.7305-0.2, 89.9218, 93.7693};


for(Long64_t i = 0; i < chain.GetEntries(); i++) {
		chain.GetEntry(i);
		if(i%1000 == 0) cout<<"Entry "<<i<<endl;


        vector<uint16_t> nstubs_per_mod_BX(6,0);
        vector<vector<float>> stubs_mod0(40); 


// U and V from strip number to cm in local reference frame
                auto absU_=[&](int mod,int stub,vector<vector<float>> stubs_mod){ double posStripU = stubs_mod.at(mod).at(stub);
                                      return (- (posStripU - 1016/2.)*0.009 -0.009/2.);};
                auto absV_=[&](int mod,int stub,vector<vector<float>> stubs_mod){ double posStripV = stubs_mod.at(mod).at(stub);
                                      return ( + (posStripV - 1016/2.)*0.009 +0.009/2.);};
// X and Y from strip number to cm in local reference frame
                auto absX_=[&](int mod,int stub,vector<vector<float>> stubs_mod){ double posStripU = stubs_mod.at(mod).at(stub);
                                      return (cos(0.233)*(- (posStripU - 1016/2.)*0.009 -0.009/2.));};
                auto absY_=[&](int mod,int stub,vector<vector<float>> stubs_mod){ double posStripV = stubs_mod.at(mod).at(stub);
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

// m and q parameters for XZ and YZ track
                auto l_mXZ = [&](int stub_mod0, int stub_mod4,vector<vector<float>> stubs_mod){return ((absX_(4,stub_mod4,stubs_mod)+Xoffsets[4])-(absX_(0,stub_mod0,stubs_mod)+Xoffsets[0]))/(absZ04_(4,absX_(4,stub_mod4,stubs_mod),Xoffsets[4])-absZ04_(0,absX_(0,stub_mod0,stubs_mod),Xoffsets[0]));};
                auto l_mYZ = [&](int stub_mod1, int stub_mod5,vector<vector<float>> stubs_mod){return ((absY_(5,stub_mod5,stubs_mod)+Yoffsets[5])-(absY_(1,stub_mod1,stubs_mod)+Yoffsets[1]))/(absZ15_(5,absY_(5,stub_mod5,stubs_mod),Yoffsets[5])-absZ15_(1,absY_(1,stub_mod1,stubs_mod),Yoffsets[1]));};
                auto l_qXZ = [&](int stub_mod0, int stub_mod4,vector<vector<float>> stubs_mod){return (absZ04_(4,absX_(4,stub_mod4,stubs_mod),Xoffsets[4])*(absX_(0,stub_mod0,stubs_mod)+Xoffsets[0]) - absZ04_(0,absX_(0,stub_mod0,stubs_mod),Xoffsets[0])*(absX_(4,stub_mod4,stubs_mod)+Xoffsets[4]))/(absZ04_(4,absX_(4,stub_mod4,stubs_mod),Xoffsets[4])-absZ04_(0,absX_(0,stub_mod0,stubs_mod),Xoffsets[0]));};
                auto l_qYZ = [&](int stub_mod1, int stub_mod5,vector<vector<float>> stubs_mod){return (absZ15_(5,absY_(5,stub_mod5,stubs_mod),Yoffsets[5])*(absY_(1,stub_mod1,stubs_mod)+Yoffsets[1]) - absZ15_(1,absY_(1,stub_mod1,stubs_mod),Yoffsets[1])*(absY_(5,stub_mod5,stubs_mod)+Yoffsets[5]))/(absZ15_(5,absY_(5,stub_mod5,stubs_mod),Yoffsets[5])-absZ15_(1,absY_(1,stub_mod1,stubs_mod),Yoffsets[1]));};

// X or Y position on the track at a given Z
                auto pos_on_track = [](double q, double m, double z){return (q + m*z);};



        for(int i = 0; i<stub_Link->size(); i++)
        {
         int link= stub_Link->at(i);
         nstubs_per_mod_BX.at(link)+=1;
	 stubs_mod0.at(link).push_back(stub_LocalX->at(i));
        }

histo_nstubs_BX->Fill(stub_Link->size());


        if(nstubs_per_mod_BX.at(0)==1 and nstubs_per_mod_BX.at(1)==1 and nstubs_per_mod_BX.at(4)==1 and nstubs_per_mod_BX.at(5)==1)
        {

	  double diffXZ=l_mXZ(0,0,stubs_mod0);
	  double diffYZ=l_mYZ(0,0,stubs_mod0);
	  double qx = l_qXZ(0,0,stubs_mod0);
	  double qy = l_qYZ(0,0,stubs_mod0);

                 double X_ru = pos_on_track(qx,diffXZ,55.47-0.2);
                 double X_rv = pos_on_track(qx,diffXZ,56.7305-0.2);
                 double Y_ru = pos_on_track(qy,diffYZ,55.47-0.2);
                 double Y_rv = pos_on_track(qy,diffYZ,56.7305-0.2);


                 double posU=newU(0.7854,X_ru,Y_ru); // expected U
                 double posV=newV(0.7854,X_rv,Y_rv); // expected V


                h_posU_expected->Fill(posU);
                h_posV_expected->Fill(posV);


if(abs(posU)<4 and abs(posV)<4){ // Restricting to 64cm^2 the sensors area where I search for stubs
                stub4++;

  //MODULE U

//If i don't find anything on the current BX, I search the stub in previous and next BX
          if(nstubs_per_mod_BX.at(2)==0){

int Bx=stub_Bx->at(0);

//--------Previous BX----------

chain.GetEntry(i-1);
        vector<int> nstubs_per_mod_BX1(6,0);
        vector<vector<float>> stubs_mod1(40);

        for(int i = 0; i<stub_Link->size(); i++)
        {
         int link= stub_Link->at(i);
         nstubs_per_mod_BX1.at(link)+=1;
         stubs_mod1.at(link).push_back(stub_LocalX->at(i));
        }

double res_U1;
	  if(nstubs_per_mod_BX1.at(2)>0){
for(int s=0; s<nstubs_per_mod_BX1.at(2); s++){

//evaluate the position in the local reference frame and then calculate the residual res_U1 with expected position posU
                 double pos_off_U1=absU_(2,s,stubs_mod1)+newU(0.7854,Xoffsets[2],Yoffsets[2]);
                 res_U1 = pos_off_U1 - posU;

//If the residual is compatible with the residual's distribution, then the stub is found
if(res_U1>-0.02 and res_U1<0.015){
		 found=1;two=1;  residual_X->Fill(res_U1);
			}
		}
	}
if(two==1){stubU1++;h_bx0->Fill(stub_Bx->at(0)-Bx); }
two=0;

//--------Next BX--------
chain.GetEntry(i+1);

        vector<int> nstubs_per_mod_BX2(6,0);
        vector<vector<float>> stubs_mod2(40);

        for(int i = 0; i<stub_Link->size(); i++)
        {
         int link= stub_Link->at(i);
         nstubs_per_mod_BX2.at(link)+=1;
         stubs_mod2.at(link).push_back(stub_LocalX->at(i));
        }

 double res_U2;
          if(nstubs_per_mod_BX2.at(2)>0){
for(int s=0; s<nstubs_per_mod_BX2.at(2); s++){

//evaluate the position in the local reference frame and then calculate the residual res_U1 with expected position posU
                 double pos_off_U2=absU_(2,s,stubs_mod2)+newU(0.7854,Xoffsets[2],Yoffsets[2]);
                  res_U2 = pos_off_U2 - posU;

//If the residual is compatible with the residual's distribution, then the stub is found
if(res_U2>-0.02 and res_U2<0.015){
		found=1;two=1;
			}
		}
              }
if(two==1){stubU2++; h_bx1->Fill(stub_Bx->at(0)-Bx);        residual_Y->Fill(res_U2);}
two=0;
chain.GetEntry(i);

//If I don't find any stub in any BX, then I look where it is expected to be
          if(found==0){corr_mod.at(5)->Fill(posU,posV);}

stubs_mod1.clear();
stubs_mod2.clear();
nstubs_per_mod_BX1.clear();
nstubs_per_mod_BX2.clear();
	}

//current BX
else if(nstubs_per_mod_BX.at(2)>0){

for(int s=0; s<nstubs_per_mod_BX.at(2); s++){

//If the residual is compatible with the residual's distribution, then the stub is found
                 double pos_off_U3=absU_(2,s,stubs_mod0)+newU(0.785,Xoffsets[2],Yoffsets[2]);
                 double res_U3 = pos_off_U3 - posU;
                 residual_Z->Fill(res_U3);

//If I don't find any stub in any BX, then I look where it is expected to be
if(res_U3>-0.02 and res_U3<0.015){
		 found=1;two=1;
                  }
		}
	}
if(two==1)stubU3++;
two=0;




    //MODULE V
         if(nstubs_per_mod_BX.at(3)==0){

int Bx_v=stub_Bx->at(0);

//Previous BX
chain.GetEntry(i-1);
        vector<int> nstubs_per_mod_BX1_v(6,0);
        vector<vector<float>> stubs_mod1_v(40);

        for(int i = 0; i<stub_Link->size(); i++)
        {
         int link= stub_Link->at(i);
         nstubs_per_mod_BX1_v.at(link)+=1;
         stubs_mod1_v.at(link).push_back(stub_LocalX->at(i));
        }

          if(nstubs_per_mod_BX1_v.at(3)>0){
for(int s=0; s<nstubs_per_mod_BX1_v.at(3); s++){

//If the residual is compatible with the residual's distribution, then the stub is found
                 double pos_off_V1=absV_(3,s,stubs_mod1_v)+newV(0.7854,Xoffsets[3],Yoffsets[3]);
                 double res_V1 = pos_off_V1 - posV;
                 residual_X_v->Fill(res_V1);
//If I don't find any stub in any BX, then I look where it is expected to be
if(res_V1>-0.05 and res_V1<0.015){
                 found_v=1; two=1;
                        }
                }
	}
if(two==1){stubV1++; h_bx0_v->Fill(stub_Bx->at(0)-Bx_v);}
two=0;


//next BX
chain.GetEntry(i+1);

        vector<int> nstubs_per_mod_BX2_v(6,0);
        vector<vector<float>> stubs_mod2_v(40);

        for(int i = 0; i<stub_Link->size(); i++)
        {
         int link= stub_Link->at(i);
         nstubs_per_mod_BX2_v.at(link)+=1;
         stubs_mod2_v.at(link).push_back(stub_LocalX->at(i));
        }

          if(nstubs_per_mod_BX2_v.at(3)>0){

for(int s=0; s<nstubs_per_mod_BX2_v.at(3); s++){

//If the residual is compatible with the residual's distribution, then the stub is found
                 double pos_off_V2=absV_(3,s,stubs_mod2_v)+newV(0.7854,Xoffsets[3],Yoffsets[3]);
                 double res_V2 = pos_off_V2 - posV;
        residual_Y_v->Fill(res_V2);

//If I don't find any stub in any BX, then I look where it is expected to be
if(res_V2>-0.05 and res_V2<0.015){
                found_v=1; two=1;
                        }
		}
              }
if(two==1){stubV2++; h_bx1_v->Fill(stub_Bx->at(0)-Bx_v);}
two=0;

chain.GetEntry(i);

stubs_mod1_v.clear();
stubs_mod2_v.clear();
nstubs_per_mod_BX1_v.clear();
nstubs_per_mod_BX2_v.clear();
	}

//Current BX
          if(nstubs_per_mod_BX.at(3)>0){
for(int s=0; s<nstubs_per_mod_BX.at(3); s++){

//If the residual is compatible with the residual's distribution, then the stub is found
                 double pos_off_V3=absV_(3,s,stubs_mod0)+newV(0.7854,Xoffsets[3],Yoffsets[3]);
                 double res_V3 = pos_off_V3 - posV;
                 residual_Z_v->Fill(res_V3);

//If I don't find any stub in any BX, then I look where it is expected to be
if(res_V3>-0.05 and res_V3<0.015){
                 found_v=1;two=1;
                  }
                }
	}
if(two==1)stubV3++;
two=0;

if(found==0){corr_mod.at(2)->Fill(posU,posV);}
if(found_v==0){corr_mod.at(3)->Fill(posU,posV);}
if(found==0 and found_v==0){corr_mod.at(5)->Fill(posU,posV);}


}
}
stubs_mod0.clear();
nstubs_per_mod_BX.clear();
found=0;
found_v=0;

nstubs_per_mod_BX.clear();

} //end of general for



TCanvas n ("canvas_stubs_BX","stubs per BX", 700, 500);
  histo_nstubs_BX->Draw();
  n.SaveAs("histo_nstubs_BX.pdf");

TCanvas m ("canvas_stubs_BX_mod","stubs per mod per BX", 700, 500);
  m.Divide(2,3);
  m.cd(1);
  histo_nstubs_per_mod_BX.at(0)->Draw();
        for(int c=1; c<6; c++)
        {m.cd(c+1);
         histo_nstubs_per_mod_BX.at(c)->Draw(); 
        }
  m.SaveAs("histo_nstubs_per_mod_BX.pdf");


double ratioU=stubU1/stub4;
double ratioV=stubU2/stub4;
double ratio5=stubU3/stub4;

double ratioV1=stubV1/stub4;
double ratioV2=stubV2/stub4;
double ratioV3=stubV3/stub4;

cout << "events with stub in U when 0145 BX-1 have stubs: " << stubU1 << " over " << stub4 << ", so ratio is " << ratioU << endl;
cout << "events with stub in U when 0145 BX+1 have stubs: " << stubU2 << " over " << stub4 << ", so ratio is " << ratioV << endl;
cout << "events with stub in U when 0145 BX have stubs: " << stubU3 << " over " << stub4 << ", so ratio is " << ratio5 << endl;

cout << "events with stub in V when 0145 BX-1 have stubs: " << stubV1 << " over " << stub4 << ", so ratio is " << ratioV1 << endl;
cout << "events with stub in V when 0145 BX+1 have stubs: " << stubV2 << " over " << stub4 << ", so ratio is " << ratioV2 << endl;
cout << "events with stub in V when 0145 BX have stubs: " << stubV3 << " over " << stub4 << ", so ratio is " << ratioV3 << endl;

TCanvas n4("res","res", 700,700);
n4.Divide(1,3);
n4.cd(1);
gPad->SetLogy();
residual_X->Fit("gaus");
residual_X->Draw();
gStyle->SetOptFit(1011);
n4.cd(2);
gPad->SetLogy();
residual_Y->Fit("gaus");
residual_Y->Draw();
gStyle->SetOptFit(1011);
n4.cd(3);
gPad->SetLogy();
residual_Z->Fit("gaus");
residual_Z->Draw();
gStyle->SetOptFit(1011);
n4.SaveAs("BX_resU.pdf");

TCanvas n0("res_v","res_v", 700,700);
n0.Divide(1,3);
n0.cd(1);
gPad->SetLogy();
residual_X_v->Fit("gaus");
residual_X_v->Draw();
gStyle->SetOptFit(1011);
n0.cd(2);
gPad->SetLogy();
residual_Y_v->Fit("gaus");
residual_Y_v->Draw();
gStyle->SetOptFit(1011);
n0.cd(3);
gPad->SetLogy();
residual_Z_v->Fit("gaus");
residual_Z_v->Draw();
gStyle->SetOptFit(1011);
n0.SaveAs("BX_resV.pdf");


TCanvas n1("bxc","bxc", 700,700);
n1.Divide(1,2);
n1.cd(1);
h_bx0->SetLineColor(kRed);
h_bx0->Draw();
h_bx1->Draw("same");
n1.cd(2);
h_bx0_v->SetLineColor(kRed);
h_bx0_v->Draw();
h_bx1_v->Draw("same");
n1.SaveAs("bx.pdf");

TCanvas n2("pos","pos", 700,700);
n2.Divide(1,3);
n2.cd(1);
corr_mod.at(2)->Draw("COLZ");
n2.cd(2);
corr_mod.at(3)->Draw("COLZ");
n2.cd(3);
corr_mod.at(5)->Draw("COLZ");
n2.SaveAs("pos.pdf");


}

