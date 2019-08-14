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

int compareDelphesL1(){//main
  SetTdrStyle();
  gStyle->SetPadRightMargin(0.05);

  std::string finName1 = "L1BitStudyRecluster.root";
  TFile *fin1 = TFile::Open(finName1.c_str());
  if (!fin1) {
    std::cout << " Input file " << finName1 << " not found." << std::endl;
    return 1;
  }
  else {
    std::cout << " File found: " << fin1->GetName() << std::endl;
  }
  fin1->cd();

  TTree *tree1 = (TTree*)gDirectory->Get("LightTree");

  if (!tree1) return 1;

  std::string finName2 = "VBFH_200PU.root";
  TFile *fin2 = TFile::Open(finName2.c_str());
  if (!fin2) {
    std::cout << " Input file " << finName2 << " not found." << std::endl;
    return 1;
  }
  else {
    std::cout << " File found: " << fin2->GetName() << std::endl;
  }
  fin2->cd();

  TTree *tree2 = (TTree*)gDirectory->Get("LightTree");

  if (!tree2) return 1;
  
  
  const unsigned nV = 7;
  std::string varName[nV] = {
    "GenJet1_pt",
    "GenJet2_pt",
    "GenJet1_eta",
    "GenJet2_eta",
    "Gen_Mjj",
    "Gen_detajj",
    "Gen_dphijj"
  };
  std::string varLabel[nV] = {
    "Jet1pt",
    "Jet2pt",
    "Jet1eta",
    "Jet2eta",
    "Mjj",
    "detajj",
    "dphijj"
  };

  //0 = left, 1 = middle, 2 = right
  const unsigned legPlace[nV] = {2,2,1,1,2,0,2};
  
  const unsigned nBins = 100;
  const double binMin[nV] = {30,30,-5,-5,500,0,0};
  const double binMax[nV] = {400,300,5,5,5000,10,3.1416};
  std::string axisTitles[nV] = {
    ";p_{T}^{j1} (GeV);Events",
    ";p_{T}^{j2} (GeV);Events",
    ";#eta^{j1};Events",
    ";#eta^{j2};Events",
    ";M_{jj} (GeV);Events",
    ";#Delta#eta_{jj};Events" ,
    ";#Delta#phi_{jj};Events"
  };

  std::string baseline = "Gen_Mjj>1000 && Gen_detajj>4 && GenJet1_pt>70 && GenJet2_pt>40 && GenJet1_eta*GenJet2_eta<0 && Gen_dphijj<2";


  TCanvas *myc[nV];
  gStyle->SetOptStat(0);

  TFile *outFile = TFile::Open("HistosCompa.root","RECREATE");
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
    label.str("");
    label << "h1_" << varLabel[iV];
    h1[iV] = new TH1F(label.str().c_str(),axisTitles[iV].c_str(),nBins,binMin[iV],binMax[iV]);
    h1[iV]->SetLineColor(1);
    h1[iV]->SetLineWidth(2);

    myc[iV]->cd();
    tree1->Draw((varName[iV]+">>"+label.str()).c_str(),baseline.c_str());

    label.str("");
    label << "h2_" << varLabel[iV];
    h2[iV] = new TH1F(label.str().c_str(),axisTitles[iV].c_str(),nBins,binMin[iV],binMax[iV]);
    h2[iV]->SetLineColor(2);
    h2[iV]->SetLineWidth(2);

    tree2->Draw((varName[iV]+">>"+label.str()).c_str(),baseline.c_str());

    h1[iV]->DrawNormalized();
    h2[iV]->DrawNormalized("same");

    bool leftLeg = legPlace[iV]==0;
    bool midLeg = legPlace[iV]==1;
    TLegend *leg = new TLegend(leftLeg?0.17:midLeg?0.41:0.7,0.7,leftLeg?0.41:midLeg?0.67:0.94,0.94);
    leg->AddEntry(h1[iV],"L1","L");
    leg->AddEntry(h2[iV],"Delphes","L");
    leg->Draw("same");  
    myc[iV]->Update();
    myc[iV]->Print(("PLOTS/Compa_"+varLabel[iV]+".pdf").c_str());
    
    
  }//loop on vars
  
  outFile->Write();
  
  return 0;
  
}//main
