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
        TH1F* histo_nstubs_BX0145=new TH1F("histo_nstubs_BX0145", "nstubs per BX just in mod. 0145",50,0,50);
	TH1F* histo_1stub= new TH1F("ones","when #stub=1, which module is?",6,0,6);

        std::vector<TH1F*> histo_nstubs_per_mod_BX(6);
        histo_nstubs_per_mod_BX.at(0)=new TH1F("histo_nstubs_mod0_BX", "nstubs per BX mod 0",40,0,40);
        histo_nstubs_per_mod_BX.at(1)=new TH1F("histo_nstubs_mod1_BX", "nstubs per BX mod 1",40,0,40);
        histo_nstubs_per_mod_BX.at(2)=new TH1F("histo_nstubs_mod2_BX", "nstubs per BX mod 2",40,0,40);
        histo_nstubs_per_mod_BX.at(3)=new TH1F("histo_nstubs_mod3_BX", "nstubs per BX mod 3",40,0,40);
        histo_nstubs_per_mod_BX.at(4)=new TH1F("histo_nstubs_mod4_BX", "nstubs per BX mod 4",40,0,40);
        histo_nstubs_per_mod_BX.at(5)=new TH1F("histo_nstubs_mod5_BX", "nstubs per BX mod 5",40,0,40);

	TH1D* h_hist0=new TH1D("h_hist0","localX of module0 when #stubs==1",1250,0,2500);
        TH1D* h_hist1=new TH1D("h_hist1","localX of module1 when #stubs==1",1250,0,2500);
        TH1D* h_hist2=new TH1D("h_hist2","localX of module2 when #stubs==1",1250,0,2500);
        TH1D* h_hist3=new TH1D("h_hist3","localX of module3 when #stubs==1",1250,0,2500);
        TH1D* h_hist4=new TH1D("h_hist4","localX of module4 when #stubs==1",1250,0,2500);
        TH1D* h_hist5=new TH1D("h_hist5","localX of module5 when #stubs==1",1250,0,2500);

        TH1D* h_hist00=new TH1D("h_hist00","localX of module0 when #stubs>0",1250,0,2500);
        TH1D* h_hist10=new TH1D("h_hist10","localX of module1 when #stubs>0",1250,0,2500);
        TH1D* h_hist20=new TH1D("h_hist20","localX of module2 when #stubs>0",1250,0,2500);
        TH1D* h_hist30=new TH1D("h_hist30","localX of module3 when #stubs>0",1250,0,2500);
        TH1D* h_hist40=new TH1D("h_hist40","localX of module4 when #stubs>0",1250,0,2500);
        TH1D* h_hist50=new TH1D("h_hist50","localX of module5 when #stubs>0",1250,0,2500);

        TH1D* h_histSID=new TH1D("h_histSID","localX of module3 when #stubs>0",63284,3136716,3200000);
        std::vector<TH2F*> corr_mod;
          corr_mod.push_back(new TH2F("corr_u3", "expected UV when there are 2 stubs", 100,-5,5,100,-5,5));

double stub2=0;
double stub3=0;
double stub4=0;


double Xoffsets[6]={-0.070300968495475138, 0.043476838127823117, 0.081218458321752882, 0.068428399849415442, -0.079017060072807038, -0.051013420927518687};
double Yoffsets[6]={-0.17741353189702608, 0.0079311725718771414, -0.080868219919505421, 0.068374266162434874, -0.05630130387604667, 0.0046751032833262096};

//double Xoffsets[6]={0.,0.,0.,0.,0.,0.};
//double Yoffsets[6]={0.,0.,0.,0.,0.,0.};

double tilt[6]={0.233,-0.233,0.,0.,0.233,-0.233};
double posZ[6]={18.0218, 21.8693, 55.4700-0.2, 56.7305-0.2, 89.9218, 93.7693};


for(Long64_t i = 0; i < 10000; i++) { //chain.GetEntries(); i++) {
                chain.GetEntry(i);
                if(i%1000 == 0) cout<<"Entry "<<i<<endl;


        vector<uint16_t> nstubs_per_mod_BX(6,0);
        vector<vector<float>> stubs_mod0(40);
        vector<vector<uint32_t>> stubs_supID(40);

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


int no=0;
        for(int j = 0; j<stub_Link->size(); j++)
        {
         int link= stub_Link->at(j);
         nstubs_per_mod_BX.at(link)+=1;
	 if(link!=2 and link!=3) no+=1;
         stubs_mod0.at(link).push_back(stub_LocalX->at(j));
         stubs_supID.at(link).push_back(stub_SuperID->at(j));
	cout << "stub_Link : " << link << ", size " << stub_Link->size() << endl;
        }

        if(nstubs_per_mod_BX.at(0)==1)  h_hist0->Fill(stubs_mod0.at(0).at(0));
        if(nstubs_per_mod_BX.at(1)==1)  h_hist1->Fill(stubs_mod0.at(1).at(0));
        if(nstubs_per_mod_BX.at(2)==1 and nstubs_per_mod_BX.at(0)==0 and nstubs_per_mod_BX.at(1)==0 and nstubs_per_mod_BX.at(3)==0 and nstubs_per_mod_BX.at(4)==0 and nstubs_per_mod_BX.at(5)==0)  h_hist2->Fill(stubs_mod0.at(2).at(0));
        if(nstubs_per_mod_BX.at(3)==1)  h_hist3->Fill(stubs_mod0.at(3).at(0));
        if(nstubs_per_mod_BX.at(4)==1)  h_hist4->Fill(stubs_mod0.at(4).at(0));
        if(nstubs_per_mod_BX.at(5)==1)  h_hist5->Fill(stubs_mod0.at(5).at(0));

        if(nstubs_per_mod_BX.at(0)!=0)  h_hist00->Fill(stubs_mod0.at(0).at(0));
        if(nstubs_per_mod_BX.at(1)!=0)  h_hist10->Fill(stubs_mod0.at(1).at(0));
        if(nstubs_per_mod_BX.at(2)!=0)  h_hist20->Fill(stubs_mod0.at(2).at(0));
        if(nstubs_per_mod_BX.at(3)!=0 and stubs_mod0.at(3).at(0)>=2022 and stubs_mod0.at(3).at(0)<=2024) 
		{h_hist30->Fill(stubs_mod0.at(3).at(0)); h_histSID->Fill(stubs_supID.at(3).at(0));}
        if(nstubs_per_mod_BX.at(4)!=0)  h_hist40->Fill(stubs_mod0.at(4).at(0));
        if(nstubs_per_mod_BX.at(5)!=0)  h_hist50->Fill(stubs_mod0.at(5).at(0));

histo_nstubs_BX0145->Fill(no);
histo_nstubs_BX->Fill(stub_Link->size());
//         if(nstubs_per_mod_BX.at(3)==0 and nstubs_per_mod_BX.at(2)==0) histo_nstubs_BX0145->Fill(stub_Link->size());

histo_nstubs_BX0145->Fill(no);
if(stub_Link->size()==1) histo_1stub->Fill(stub_Link->at(0));

for(int m=0; m<6; m++){ histo_nstubs_per_mod_BX.at(m)->Fill(nstubs_per_mod_BX.at(m));}

        if(nstubs_per_mod_BX.at(0)==1 and nstubs_per_mod_BX.at(1)==1 and nstubs_per_mod_BX.at(4)==1 and nstubs_per_mod_BX.at(5)==1)
	{
//	 if(nstubs_per_mod_BX.at(3)==0 and nstubs_per_mod_BX.at(2)==0)
double pos_X=absX_(0,0,stubs_mod0)+Xoffsets[0];
double pos_Y=absY_(1,0,stubs_mod0)+Yoffsets[1];
corr_mod.at(0)->Fill(pos_X,pos_Y);
stub4++;
double pos_off_U=-999;
double pos_off_V=-999;
	 if(nstubs_per_mod_BX.at(2)==1) {pos_off_U=absU_(2,0,stubs_mod0);//+newU(0.7854-0.11530836310634285,Xoffsets[2],Yoffsets[2]);
					 stub2++;}
	 if(nstubs_per_mod_BX.at(3)==1) {pos_off_V=absV_(3,0,stubs_mod0);//+newV(0.7854-0.18331295053815419,Xoffsets[3],Yoffsets[3]);
					 stub3++;}
	 //if(nstubs_per_mod_BX.at(3)==1 and nstubs_per_mod_BX.at(2)==1) corr_mod.at(0)->Fill(pos_off_U,pos_off_V);
	}

stubs_supID.clear();
stubs_mod0.clear();

} //end of general for

double ratio2=stub2/stub4;
double ratio3=stub3/stub4;

cout << "ratio2: " << ratio2 << "stub2, stub4  " << stub2<< ", " << stub4 <<endl;
cout << "ratio3: " << ratio3 << endl;

TCanvas n ("canvas_stubs_BX","stubs per BX", 700, 500);
  histo_1stub->Draw();
  n.SaveAs("histo_1.pdf");

TCanvas n0 ("canvas_stubs_BX","stubs per BX", 700, 500);
n0.Divide(1,2);
 n0.cd(1);
  histo_nstubs_BX->Draw();
 n0.cd(2);
  histo_nstubs_BX0145->Draw();
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

TCanvas n2("xy","xy", 700,700);
corr_mod.at(0)->Draw("COLZ");
n2.SaveAs("UV.pdf");

TCanvas n3("n3","n3",1000,1000);
n3.Divide(2,3);
n3.cd(1);
h_hist0->Draw( );
n3.cd(2);
h_hist1->Draw( );
n3.cd(3);
h_hist2->Draw( );
n3.cd(4);
h_hist3->Draw( );
n3.cd(5);
h_hist4->Draw( );
n3.cd(6);
h_hist5->Draw( );
n3.SaveAs("local.pdf");

TCanvas n4("n4","n4",1000,1000);
n4.Divide(2,3);
n4.cd(1);
h_hist00->Draw( );
n4.cd(2);
h_hist10->Draw( );
n4.cd(3);
h_hist20->Draw( );
n4.cd(4);
h_hist30->Draw( );
n4.cd(5);
h_hist40->Draw( );
n4.cd(6);
h_hist50->Draw( );
n4.SaveAs("local0.pdf");
}
