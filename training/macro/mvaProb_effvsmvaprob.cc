#include "TROOT.h"
#include "TKey.h"
#include "TFile.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1F.h"
#include "TObjArray.h"
#include "THStack.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TCut.h"
#include "TF1.h"
#include "TStyle.h"
#include <iostream>

#include <vector>
#include <fstream>
#include <sstream>

#include <cstring>
#include <string>
#include <TROOT.h>

// Double_t myfunction(Double_t *x, Double_t *par)
// {
//    Float_t xx =x[0];
//    // Double_t f = TMath::Abs(par[0]*sin(par[1]*xx)/xx);
//    Double_t f = par[0]+par[1]−par[2]+par[3]−par[4]+par[1]*xx+par[2]*xx*xx+par[3]*xx*xx*xx+par[4]*xx*xx*xx*xx ;
//    return f;
// }
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

void mvaProb_effvsmvaprob(){

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
  // gStyle->SetOptStat(0);
  // gStyle->SetOptFit(0);
  // gROOT->ForceStyle();

  int nBinMVAprov = 50;
  float minMVAprob = -1;
  float maxMVAprob = 1;
  TFile *File;

  //File = TFile::Open("Vtx942ProbTrainingTree_GluGluHToGG_M120_13TeV_amcatnloFXFX_pythia8.root");
  File = TFile::Open("GluGluHToGG_M120_13TeV_amcatnloFXFX_pythia8_Training2.root");
  TDirectory *dir = (TDirectory*)File->Get("commissioning");

  float VtxProb = -99.;
  float dZtrue = -99.;
  int dipho_index = -99;
  float dipho_pt = -99.;
  float evWeight =1;
  float nVtx;
  float NConv;

  TTree *t_tree = (TTree*)dir->Get("diphoTree");

  int nevent = t_tree->GetEntries();

  Double_t vpTBins[] = {0,3,5,7,10,13,17,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160};
  Int_t nbins_pT = sizeof(vpTBins)/sizeof(Double_t) - 1;


  Double_t vpTBins2[] = {0,5,10,15,20,30,40,50,75,100,150};
  Int_t nbins_pT2 = sizeof(vpTBins)/sizeof(Double_t) - 1;

  TH1F *hGood[3];
  TH1F *hTot[3];
  TH1F *hGoodPt[3];
  TH1F *hTotPt[3];

  t_tree->SetBranchAddress("vtxProb",&VtxProb);
  t_tree->SetBranchAddress("NConv",&NConv);
  t_tree->SetBranchAddress("evWeight",&evWeight);
  t_tree->SetBranchAddress("DZtrue",&dZtrue);
  t_tree->SetBranchAddress("dipho_index",&dipho_index);
  t_tree->SetBranchAddress("pt",&dipho_pt);
  //  t_tree->SetBranchAddress("NVert",&nVtx);


  TString hName[3];
  TString hNameGood[3];
  TString hNamePt[3];
  TString hNameGoodPt[3];

  for (int i=0;i<3;i++) hName[i]="hTot";
  for (int i=0;i<3;i++) hNameGood[i]="hGood";
  for (int i=0;i<3;i++) hNamePt[i]="hTotPt";
  for (int i=0;i<3;i++) hNameGoodPt[i]="hGoodPt";
  hName[1]+="NoConv";
  hName[2]+="Conv";
  hNameGood[1]+="NoConv";
  hNameGood[2]+="Conv";
  hNamePt[1]+="NoConv";
  hNamePt[2]+="Conv";
  hNameGoodPt[1]+="NoConv";
  hNameGoodPt[2]+="Conv";

  for (int i=0;i<3;i++){
    
    hGood[i] = new TH1F(hNameGood[i].Data(),hNameGood[i].Data(),nBinMVAprov,minMVAprob,maxMVAprob);
    hTot[i] = new TH1F(hName[i].Data(),hName[i].Data(),nBinMVAprov,minMVAprob,maxMVAprob);
    hGoodPt[i] = new TH1F(hNameGoodPt[i].Data(),hNameGoodPt[i].Data(),nbins_pT,vpTBins);
    hTotPt[i] = new TH1F(hNamePt[i].Data(),hNamePt[i].Data(),nbins_pT,vpTBins);

  }
  
  TCut cutGood[3];
  cutGood[0]="(abs(DZtrue)<1)*evWeight";
  cutGood[1]="(abs(DZtrue<1) && NConv==0)*evWeight";
  cutGood[2]="(abs(DZtrue<1) && NConv>0)*evWeight";
  
  TCut cut[3];
  cut[0]="evWeight"; 
  cut[1]="(NConv==0)*evWeight";
  cut[2]="(NConv>0)*evWeight";


  TCanvas *test = new TCanvas("test","test");
  for (int i=0;i<3;i++){    
    cout<<hNameGood[i]<<" "<<hName[i]<< endl;
    cout<<hNameGoodPt[i]<<" "<<hNamePt[i]<< endl;
    t_tree->Draw("vtxProb>>"+hNameGood[i],cutGood[i]);
    t_tree->Draw("vtxProb>>"+hName[i],cut[i]);
    t_tree->Draw("pt>>"+hNameGoodPt[i],cutGood[i]);
    t_tree->Draw("pt>>"+hNamePt[i],cut[i]);
    if(i==1){
      hTotPt[i]->Draw();
      hGoodPt[i]->Draw("same");
    }
  }

  TEfficiency* eff[3]; 
  TEfficiency* effPt[3]; 
  for (int i=0;i<3;i++){  
    eff[i]=new TEfficiency(*hGood[i],*hTot[i]);
    effPt[i]=new TEfficiency(*hGoodPt[i],*hTotPt[i]);
  }

  TF1* myfunc[3];
  TString fname[3];


  for (int i=0;i<3;i++){
    fname[i]="fit";
    fname[i]+=i;
    // myfunc[i] = new TF1(fname[i].Data(), "1+[0]−[1]+[2]−[3]+[0]*x+[1]*x*x+[2]*x*x*x+[3]*x*x*x*x",-1.,1.);
    // myfunc[i] = new TF1(fname[i].Data(), "1+[1]−[2]+[3]−[4]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x",-1.,1.);
    // myfunc[i] = new TF1("func", myfunction,-1.,1.);
    // myfunc[i]->SetParameters(0,1);
    //myfunc[i]
    // myfunc[i] = new TF1("func", "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x",-1,1);
    // myfunc[i]->SetParameters(0, myfunc[i]->Eval(-1.) );
    myfunc[i] = new TF1(fname[i].Data(), "1+[0]-[1]+[2]-[3]+[0]*x+[1]*x*x+[2]*x*x*x+[3]*x*x*x*x",-1,1);
    // myfunc[i] = new TF1(fname[i].Data(), "[0]-[1]*exp(log([0]/[1])*x)",-1,0.8);
    // myfunc[i] = new TF1(fname[i].Data(), "1-[0]*exp([1])+[0]*exp(-[1]*x)",-1,0.8);
   
    // if(i!=2) myfunc[i] = new TF1(fname[i].Data(), "[0]+[1]*exp(x*[2])",-1,0.8); //works no constraint
    // else myfunc[i] = new TF1(fname[i].Data(), "[0]+[1]*exp(x*[2])",-1,0.6); //works no constraint
    // if(i==2)  myfunc[i]->SetRange(-1,0.6);
    //myfunc[i] = new TF1(fname[i].Data(), "1-[1]*exp(-[2])+[1]*exp(x*[2])",-1,0.8); 
    //myfunc[i] = new TF1(fname[i].Data(), "[0]-[1]*exp(log([1]/([0]-1))*x)",-1,0.8); 
    
   
    myfunc[i]->SetLineColor(i+1);
    //myfunc[i] = new TF1(fname[i].Data(), "1.-[0]*exp([1])-[2]*exp([3])+[0]*exp(-x*[1])+[2]*exp(-x*[3])",-1,0.8);
    if(i>0)eff[i]->Fit(myfunc[i],"E");
  }


  Double_t xl1=.25, yl1=0.3, xl2=xl1+.3, yl2=yl1+.175;
  TLegend leg(xl1,yl1,xl2,yl2);
  vector<TString> leg_Caption;  
  leg_Caption.push_back("all");
  leg_Caption.push_back("nConv=0");
  leg_Caption.push_back("nConv>0");

  for (int i=0;i<3;i++)
    leg.AddEntry(eff[i],leg_Caption.at(i),"lem");
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextFont(42);

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TLatex lat;
  //lat.SetTextColor(6);
  //lat.SetTextSize(0.1);
  char buf[500]; 

  TCanvas *can = new TCanvas("can","can");
  can->cd();
  // gStyle->SetOptFit(0);
  // gROOT->ForceStyle();
  //gStyle->SetFitFormat("%3.4g"); 

  //p1 is 0 in fit parameters
  double p1,p1err,p2,p2err,p3,p3err,p4,p4err;
  char buff[500]; 

  for (int i=0;i<3;i++){  
    eff[i]->SetLineColor(i+1); 
    if(i==0) eff[i]->Draw();
    else eff[i]->Draw("same");
    //gPad->Update();
    eff[i]->SetTitle(";MVAprob;Fraction of |z_{selected} - z_{true}| < 10 mm");
    //if(i==0) eff[i]->Draw();
    //else eff[i]->Draw("same");
    myfunc[i]->SetLineColor(i+1);
    
    if(i==0) {std::cout << "Category & p_1 & p_2 & p_3 & p_4  \\\\" << std::endl;}
    
    if(i>0){
      //Print parameters and errors
      p1 = myfunc[i]->GetParameter(0);
      p1err = myfunc[i]->GetParError(0);
      p2 = myfunc[i]->GetParameter(1);
      p2err = myfunc[i]->GetParError(1);
      p3 = myfunc[i]->GetParameter(2);
      p3err = myfunc[i]->GetParError(2);
      p4 = myfunc[i]->GetParameter(3);
      p4err = myfunc[i]->GetParError(3);
    
        if(i==1) {
	        std::cout << "non converted & " << std::setprecision(4) << p1 << " \\pm " << std::setprecision(4) << p1err << " & " << std::setprecision(4) << p2 << " \\pm " << std::setprecision(4) << p2err << 
	        " & " << std::setprecision(4) << p3 << " \\pm " << std::setprecision(4) << p3err << " & " << std::setprecision(4) << p4 << " \\pm " << std::setprecision(4) << p4err <<  " \\\\ " << std::endl;
        }
        if(i==2) {
	        std::cout << std::setprecision(4) << "converted & " << p1 << " \\pm " << p1err << " & " << p2 << " \\pm " << p2err << 
	        " & " << p3 << " \\pm " << p3err << " & " << p4 << " \\pm " << p4err <<  " \\\\ " << std::endl;
        }

      myfunc[i]->Draw("same");
      cout<<" FUNC :"<<i<<" at -1:"<< myfunc[i]->Eval(-1)<< endl;
    }
    //gPad->Update();
  }

  sprintf(buf,"f(x) = p_{0} + p_{1} x + p_{2} x + p_{3} x + p_{4} x ");
  lat.DrawLatexNDC(0.5,0.85,buf);
  sprintf(buf,"p_{0} = 1 + p_{1} - p_{2} + p_{3} - p_{4} ");
  lat.DrawLatexNDC(0.65,0.75,buf);
  

  eff[2]->Draw("same");
  leg.Draw("same");
  //gPad->Update();
  // can->SaveAs("efficienciesVersusMVAProbWithFits.pdf","pdf");
  //can->SaveAs("efficienciesVersusMVAProbWithFits.png","png");
  can->Print("efficienciesVersusMVAProbWithFits.pdf");
  can->SaveAs("efficienciesVersusMVAProbWithFits.root",".root");



}

