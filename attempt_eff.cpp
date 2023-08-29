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

        TH1F* residual_X=new TH1F("residual_U", "residual_U",111,-0.05,0.05);
        TH1F* residual_Y=new TH1F("residual_V", "residual_V",111,-0.05,0.05);

TH1F* h_stub5 = new TH1F("stub5","modules where stubs are when #stubs==5",6,0,6);
TH1F* h_stub6 = new TH1F("stub6","number of events in U and V when mod 0,1,4,5 have 1 stub",6,0,6);


double stub6,stub5,stub4=0;
double stubU,stubV=0;



double Xoffsets[6]={-0.070300968495475138, 0.043476838127823117, 0.081218458321752882, 0.068428399849415442, -0.079017060072807038, -0.051013420927518687};
double Yoffsets[6]={-0.17741353189702608, 0.0079311725718771414, -0.080868219919505421, 0.068374266162434874, -0.05630130387604667, 0.0046751032833262096};

//double Xoffsets[6]={0.,0.,0.,0.,0.,0.};
//double Yoffsets[6]={0.,0.,0.,0.,0.,0.};

double tilt[6]={0.233,-0.233,0.,0.,0.233,-0.233};
double posZ[6]={18.0218, 21.8693, 55.4700-0.2, 56.7305-0.2, 89.9218, 93.7693};

for(Long64_t i = 0; i < 10000; i++) {//chain.GetEntries(); i++) {
                chain.GetEntry(i);
                if(i%1000 == 0) cout<<"Entry "<<i<<endl;


        vector<uint16_t> nstubs_per_mod_BX(6,0);
        vector<vector<float>> stubs_mod(40);
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


        for(int j = 0; j<stub_Link->size(); j++)
        {
         int link= stub_Link->at(j);
         nstubs_per_mod_BX.at(link)+=1;
         stubs_mod.at(link).push_back(stub_LocalX->at(j));
	cout << stub_LocalY->at(j) << endl;
        }

histo_nstubs_BX->Fill(stub_Link->size());



        for(int link=0; link<6; link++)
        {histo_nstubs_per_mod_BX.at(link)->Fill(nstubs_per_mod_BX.at(link));}

// requiring 1 hit in 0145 
	if(nstubs_per_mod_BX.at(0)==1 and nstubs_per_mod_BX.at(1)==1 and nstubs_per_mod_BX.at(4)==1 and nstubs_per_mod_BX.at(5)==1)
	{
        double diffXZ=l_mXZ(0,0,stubs_mod);
          double diffYZ=l_mYZ(0,0,stubs_mod);
          double qx = l_qXZ(0,0,stubs_mod);
          double qy = l_qYZ(0,0,stubs_mod);

                 double X_ru = pos_on_track(qx,diffXZ,55.47-0.2);
                 double X_rv = pos_on_track(qx,diffXZ,56.7305-0.2);
                 double Y_ru = pos_on_track(qy,diffYZ,55.47-0.2);
                 double Y_rv = pos_on_track(qy,diffYZ,56.7305-0.2);


                 double posU=newU(0.7854,X_ru,Y_ru); // expected U
                 double posV=newV(0.7854,X_rv,Y_rv); // expected V


                h_posU_expected->Fill(posU);
                h_posV_expected->Fill(posV);

//		if(abs(Y_ru)<3 and abs(Y_rv)<3 and abs(X_rv)<3 and abs(X_ru)<3){

                //h_posU_expected->Fill(posU);
                //h_posV_expected->Fill(posV);


		stub4++;

          if(nstubs_per_mod_BX.at(2)==0 and nstubs_per_mod_BX.at(3)==0){corr_mod.at(5)->Fill(posU,posV);}

	  if(nstubs_per_mod_BX.at(2)==1){

                 double pos_off_U=absU_(2,0,stubs_mod)+newU(0.785398,Xoffsets[2],Yoffsets[2]);
                 double res_U = pos_off_U - posU; cout << "res_U " << res_U << endl;
                if(res_U>-0.02 and res_U<0.01){
		corr_mod.at(0)->Fill(posU,posV);
                //h_posU->Fill(pos_off_U);
		h_posU->Fill(posU);
		 stubU++;
                 residual_X->Fill(res_U);
		} //else{corr_mod.at(2)->Fill(posU,posV); h_posU_expected_no->Fill(posU);}
	}
	else if(nstubs_per_mod_BX.at(2)==0){corr_mod.at(2)->Fill(posU,posV); h_posU_expected_no->Fill(posU);}


          if(nstubs_per_mod_BX.at(3)==1) {
                double pos_off_V=absV_(3,0,stubs_mod)+newV(45331295053815419,Xoffsets[3],Yoffsets[3]);
                 double res_V = pos_off_V - posV;
		if(res_V>-0.04 and res_V<0.01){
          	if(nstubs_per_mod_BX.at(2)==1) {corr_mod.at(4)->Fill(posU,posV);}
		corr_mod.at(1)->Fill(posU,posV);
                //h_posV->Fill(pos_off_V);
		h_posV->Fill(posV);
		 stubV++;
		 residual_Y->Fill(res_V);
		} //else{corr_mod.at(3)->Fill(posU,posV); h_posV_expected_no->Fill(posV);}
	}
		else if(nstubs_per_mod_BX.at(3)==0){corr_mod.at(3)->Fill(posU,posV); h_posV_expected_no->Fill(posV);}


//	}
}
stubs_mod.clear();
nstubs_per_mod_BX.clear();

} //end of general for
double ratioU=stubU/stub4;
double ratioV=stubV/stub4;
double ratio5=stub5/nevtot;
double ratio6=stub6/nevtot;
cout << "events with stub in U when 0145 have stubs: " << stubU << " over " << stub4 << ", so ratio is " << ratioU << endl;
cout << "events with stub in V when 0145 have stubs: " << stubV << " over " << stub4 << ", so ratio is " << ratioV << endl;

h_posU->Sumw2();
h_posU_expected->Sumw2();
h_posU_expected_no->Sumw2();
h_posV->Sumw2();
h_posV_expected_no->Sumw2();
h_posV_expected->Sumw2();

TH1F *h_ratioU = (TH1F*)h_posU->Clone("h_ratioU");
h_ratioU->Divide(h_posU,h_posU_expected);

TH1F *h_ratioV = (TH1F*)h_posV->Clone("h_ratioV");
h_ratioV->Divide(h_posV,h_posV_expected);

TCanvas n ("canvas_stubs_BX","stubs per BX", 700, 700);
  histo_nstubs_BX->Draw();
  n.SaveAs("histo_nstubs_BX.pdf");

TCanvas n1 ("canvas_stubs_BX_mod","stubs per mod per BX", 700, 700);
  n1.Divide(2,3);
n1.cd(1);
         histo_nstubs_per_mod_BX.at(0)->Draw(); 
        for(int c=1; c<6; c++)
        {n1.cd(c+1);
         histo_nstubs_per_mod_BX.at(c)->Draw(); 
        }
  n1.SaveAs("histo_nstubs_per_mod_BX.pdf");

TCanvas n2("n2","n2",700,700);
n2.Divide(1,2);
n2.cd(1);
gPad->SetLogy();
h_stub6->Draw();
h_stub6->Draw("TEXT SAME");
n2.cd(2);
h_stub5->Draw();
h_stub5->Draw("TEXT SAME");
n2.SaveAs("h_stub56.pdf");


TCanvas n4("res","res", 700,700);
n4.Divide(1,2);
n4.cd(1);
residual_X->Fit("gaus");
residual_X->Draw();
gStyle->SetOptFit(1011);
n4.cd(2);
residual_Y->Fit("gaus");
residual_Y->Draw();
gStyle->SetOptFit(1011);
n4.SaveAs("sing_res.pdf");
n4.SaveAs("sing_res.root");

TLine *line1 = new TLine(-5,1,5,1); 
line1->SetLineColor(kRed);

TCanvas n5("pos","pos", 700,700);
n5.BuildLegend();
n5.Divide(2,3);
n5.cd(1);
gPad->SetLogy();
//corr_mod[0]->Draw("COLZ");
h_posU->SetLineColor(kRed);
h_posU_expected->Draw();
h_posU_expected_no->SetLineColor(kOrange);
h_posU_expected_no->Draw("same");
h_posU->Draw("same");
n5.cd(2);
gPad->SetLogy();
//corr_mod[1]->Draw("COLZ");
h_posV->SetLineColor(kRed);
h_posV_expected_no->SetLineColor(kOrange);
h_posV_expected->Draw();
h_posV_expected_no->Draw("same");
h_posV->Draw("same");
n5.cd(3);
//h_ratioU->SetAxisRange(2, 0, "Y");
h_ratioU->Draw();
line1->Draw("same");
n5.cd(4);
//h_ratioV->SetAxisRange(2, 0, "Y");
h_ratioV->Draw();
line1->Draw("same");
n5.cd(5);
h_posU->SetLineColor(kRed);
h_posU_expected_no->SetLineColor(kOrange);
h_posU_expected->Draw();
h_posU_expected_no->Draw("same");
h_posU->Draw("same");
n5.cd(6);
//corr_mod[1]->Draw("COLZ");
h_posV->SetLineColor(kRed);
h_posV_expected->Draw();
h_posV_expected_no->SetLineColor(kOrange);
h_posV_expected_no->Draw("same");
h_posV->Draw("same");
n5.SaveAs("posUV.pdf");
n5.SaveAs("posUV.root");

TCanvas n6("xy","xy", 700,700);
n6.Divide(2,3);
n6.cd(1);
corr_mod[0]->Draw("COLZ");
n6.cd(2);
corr_mod[1]->Draw("COLZ");
n6.cd(3);
corr_mod[2]->Draw("COLZ");
n6.cd(4);
corr_mod[3]->Draw("COLZ");
n6.cd(5);
corr_mod[4]->Draw("COLZ");
n6.cd(6);
corr_mod[5]->Draw("COLZ");
n6.SaveAs("xy.pdf");

}

