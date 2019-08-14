#include <iostream>
#include <fstream> 
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <map>
#include <algorithm>

#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TLatex.h"

#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "boost/bind.hpp"
#include "boost/function.hpp"
#include "boost/format.hpp"
#include "boost/algorithm/string.hpp"
#include <boost/assign/list_of.hpp>
#include <boost/range/algorithm_ext/for_each.hpp>

#include "TDRStyle.h"

int compareMadridNew(){//main

  SetTdrStyle();
  gStyle->SetPadRightMargin(0.05);

  std::string finName1 = "L1ntuple_VBFH_200PU.root";
  TFile *fin1 = TFile::Open(finName1.c_str());
  if (!fin1) {
    std::cout << " Input file " << finName1 << " not found." << std::endl;
    return 1;
  }
  else {
    std::cout << " File found: " << fin1->GetName() << std::endl;
  }
  fin1->cd("l1PhaseIITree");

  TTree *tree1 = (TTree*)gDirectory->Get("L1PhaseIITree");

  if (!tree1) return 1;

  std::string finName2 = "VBFH125_L1NtuplePhaseII_160.root";
  TFile *fin2 = TFile::Open(finName2.c_str());
  if (!fin2) {
    std::cout << " Input file " << finName2 << " not found." << std::endl;
    return 1;
  }
  else {
    std::cout << " File found: " << fin2->GetName() << std::endl;
  }
  fin2->cd("l1PhaseIITree");

  TTree *tree2 = (TTree*)gDirectory->Get("L1PhaseIITree");

  if (!tree2) return 1;
  
  
  const unsigned nV = 8;
  std::string varName[nV] = {
    "nPuppiJets",
    "puppiJetEt[0]",
    "puppiJetEta[0]",
    "puppiJetEt[1]",
    "puppiJetEta[1]",
    "puppiJetEt[2]",
    "puppiJetEta[2]",
    "puppiMETEt"
  };
  std::string varLabel[nV] = {
    "nJets",
    "Jet1pt",
    "Jet1eta",
    "Jet2pt",
    "Jet2eta",
    "Jet3pt",
    "Jet3eta",
    "MET"
  };

  //0 = left, 1 = middle, 2 = right
  const unsigned legPlace[nV] = {2,2,1,2,1,2,1,2};
  
  const unsigned nBins = 50;
  const double binMin[nV] = {0,0,-5,0,-5,0,-5,0};
  const double binMax[nV] = {50,400,5,400,5,400,5,400};
  std::string axisTitles[nV] = {
    ";nJets;Events",
    ";p_{T}^{j1} (GeV);Events",
    ";#eta^{j1};Events",
    ";p_{T}^{j2} (GeV);Events",
    ";#eta^{j2};Events",
    ";p_{T}^{j3} (GeV);Events",
    ";#eta^{j3};Events",
    ";MET (GeV);Events"
  };

  std::string baseline[nV] = {
    "",
    "nPuppiJets>0",
    "nPuppiJets>0",
    "nPuppiJets>1",
    "nPuppiJets>1",
    "nPuppiJets>2",
    "nPuppiJets>2",
    ""
  };

  TCanvas *myc[nV];
  gStyle->SetOptStat(0);

  TFile *outFile = TFile::Open("HistosCompaMadrid.root","RECREATE");
  outFile->cd();
  
  TH1F *h1[nV];
  TH1F *h2[nV];
  
  std::ostringstream label;

  for (unsigned iV(0); iV<nV; ++iV){
    label.str("");
    label << "myc_" << varLabel[iV];
    myc[iV] = new TCanvas(label.str().c_str(),label.str().c_str(),1);
  }
  
  for (unsigned iV(0); iV<nV; ++iV){

    //double maxY = 0;
    label.str("");
    label << "h1_" << varLabel[iV];
    h1[iV] = new TH1F(label.str().c_str(),axisTitles[iV].c_str(),nBins,binMin[iV],binMax[iV]);
    h1[iV]->SetLineColor(1);
    h1[iV]->SetLineWidth(2);

    myc[iV]->cd();
    tree1->Draw((varName[iV]+">>"+label.str()).c_str(),baseline[iV].c_str());

    label.str("");
    label << "h2_" << varLabel[iV];
    h2[iV] = new TH1F(label.str().c_str(),axisTitles[iV].c_str(),nBins,binMin[iV],binMax[iV]);
    h2[iV]->SetLineColor(2);
    h2[iV]->SetLineWidth(2);

    tree2->Draw((varName[iV]+">>"+label.str()).c_str(),baseline[iV].c_str());

    h2[iV]->DrawNormalized();
    h1[iV]->DrawNormalized("same");
    //maxY = h1[iV]->GetMaximum();
    //std::cout << maxY << std::endl;
    //if (h2[iV]->GetMaximum()>maxY) maxY = h2[iV]->GetMaximum();

    bool leftLeg = legPlace[iV]==0;
    bool midLeg = legPlace[iV]==1;
    TLegend *leg = new TLegend(leftLeg?0.17:midLeg?0.41:0.7,0.7,leftLeg?0.41:midLeg?0.67:0.94,0.94);
    leg->AddEntry(h1[iV],"L1 Madrid","L");
    leg->AddEntry(h2[iV],"L1 New","L");
    leg->Draw("same");  
    myc[iV]->Update();
    myc[iV]->Print(("PLOTS/CompaMadridNew_"+varLabel[iV]+".pdf").c_str());
    
    
  }//loop on vars
  
  outFile->Write();
  
  return 0;
  
}//main
