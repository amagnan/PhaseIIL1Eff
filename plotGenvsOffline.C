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

int plotGenvsOffline(){//main

  SetTdrStyle();
  
  std::string finName = "VBFH_200PU.root";
  TFile *fin = TFile::Open(finName.c_str());
  if (!fin) {
    std::cout << " Input file " << finName << " not found." << std::endl;
    return 1;
  }
  else {
    std::cout << " File found: " << fin->GetName() << std::endl;
  }
  fin->cd();

  TTree *tree = (TTree*)gDirectory->Get("LightTree");

  if (!tree) return 1;
  
  const unsigned nV = 12;
  std::string varName[nV] = {
    "Jet1_pt:GenJet1_pt",
    "Jet2_pt:GenJet2_pt",
    "Jet1_eta:GenJet1_eta",
    "Jet2_eta:GenJet2_eta",
    "Mjj:Gen_Mjj",
    "detajj:Gen_detajj",
    "dphijj:Gen_dphijj",
    "(Jet1_pt-GenJet1_pt)/GenJet1_pt",
    "(Jet2_pt-GenJet2_pt)/GenJet2_pt",
    "(Mjj-Gen_Mjj)/Gen_Mjj",
    "(detajj-Gen_detajj)/Gen_detajj",
    "(dphijj-Gen_dphijj)/Gen_dphijj"
  };
  std::string varLabel[nV] = {
    "Jet1pt",
    "Jet2pt",
    "Jet1eta",
    "Jet2eta",
    "Mjj",
    "detajj",
    "dphijj",
    "ptdiff1",
    "ptdiff2",
    "Mjjdiff",
    "detajjdiff",
    "dphijjdiff"
  };
  
  const bool is1D[nV] = {0,0,0,0,0,0,0,1,1,1,1,1};
  
  const unsigned nBins = 100;
  const double binMin[nV] = {30,30,-5,-5,500,0,0,-0.6,-0.6,-0.6,-0.5,-0.5};
  const double binMax[nV] = {400,300,5,5,5000,10,3.1416,1,1,1,0.5,0.5};
  std::string axisTitles[nV] = {
     ";p_{T}^{genj1} (GeV);p_{T}^{j1} (GeV);Events",
     ";p_{T}^{genj2} (GeV);p_{T}^{j2} (GeV);Events",
     ";#eta^{genj1};#eta^{j1};Events",
     ";#eta^{genj2};#eta^{j2};Events",
     ";M_{jj}^{gen} (GeV);M_{jj} (GeV);Events",
     ";#Delta#eta_{jj}^{gen};#Delta#eta_{jj};Events" ,
     ";#Delta#phi_{jj}^{gen};#Delta#phi_{jj};Events",
     ";#Delta p_{T}^{j1} /p_{T}^{genj1} ;Events",
     ";#Delta p_{T}^{j2} /p_{T}^{genj2} ;Events",
     ";#Delta M_{jj} /M_{jj}^{gen} ;Events",
     ";#Delta #Delta#eta_{jj} /#Delta#eta_{jj}^{gen} ;Events",
     ";#Delta #Delta#phi_{jj} /#Delta#phi_{jj}^{gen} ;Events"
  };
    
  std::string baseline = "Gen_Mjj>1000 && Gen_detajj>4 && Gen_dphijj<2 && GenJet1_pt>70 && GenJet2_pt>40 && GenJet1_eta*GenJet2_eta<0";
  //std::string baseline = "Mjj>800 && detajj>1 && dphijj<2.5 && Jet1_pt>30 && Jet2_pt>30 && Jet1_eta*Jet2_eta<0";
  std::string vetocuts = " && nmediumbjets==0 && ntaus==0 && nlooseEle==0 && nlooseMu==0 && ntightGamma==0";
  std::string genmatch = " && TMath::Abs(Jet1_eta-GenJet1_eta)<0.5 && TMath::Abs(Jet2_eta-GenJet2_eta)<0.5 && GenJet1_pt>0 && GenJet2_pt > 0 && Gen_Mjj > 0";

  
  TCanvas *myc[nV];
  gStyle->SetOptStat("eMR");

  TLatex lat;
  lat.SetTextSize(0.035);
  lat.SetTextColor(1);
 
  TFile *outFile = TFile::Open("HistosOfflinevsGen.root","RECREATE");
  outFile->cd();
  
  TH2F *h2D[nV];
  TH2F *h2Dvetos[nV];
  TH1F *h1D[nV];
  TH1F *h1Dvetos[nV];

  std::ostringstream label;

  for (unsigned iV(0); iV<nV; ++iV){
    label.str("");
    label << "myc_" << varLabel[iV];
    myc[iV] = new TCanvas(label.str().c_str(),label.str().c_str(),1);
    myc[iV]->cd();
  }
  
  for (unsigned iV(0); iV<nV; ++iV){
    label.str("");
    if (!is1D[iV]) {
      label << "h2D_" << varLabel[iV];
      h2D[iV] = new TH2F(label.str().c_str(),axisTitles[iV].c_str(),nBins,binMin[iV],binMax[iV],nBins,binMin[iV],binMax[iV]);
      h2D[iV]->SetStats(0);
    }
    else {
      label << "h1D_" << varLabel[iV];
      h1D[iV] = new TH1F(label.str().c_str(),axisTitles[iV].c_str(),nBins,binMin[iV],binMax[iV]);
      h1D[iV]->SetLineColor(1);
      h1D[iV]->SetLineWidth(2);
    }
    
    myc[iV]->cd();
    gPad->SetGridx(1);
    gPad->SetGridy(1);

    tree->Draw((varName[iV]+">>"+label.str()).c_str(),(baseline+genmatch).c_str(),is1D[iV]?"":"colz");

    lat.DrawLatexNDC(0.1,0.97,"Delphes 200PU");

    myc[iV]->Update();
    myc[iV]->Print(("PLOTS/GenReco_"+varLabel[iV]+".pdf").c_str());
    
    /*    label.str("");
    if (!is1D[iV]) {
      label << "h2Dvetos_" << varLabel[iV];
      h2Dvetos[iV] = new TH2F(label.str().c_str(),axisTitles[iV].c_str(),nBins,binMin[iV],binMax[iV],nBins,binMin[iV],binMax[iV]);
      h2Dvetos[iV]->SetStats(0);
    }
    else {
      label << "h1Dvetos_" << varLabel[iV];
      h1Dvetos[iV] = new TH1F(label.str().c_str(),axisTitles[iV].c_str(),nBins,binMin[iV],binMax[iV]);
      h1Dvetos[iV]->SetLineColor(1);
      h1Dvetos[iV]->SetLineWidth(2);
    }
    
    tree->Draw((varName[iV]+">>"+label.str()).c_str(),(baseline+vetocuts+genmatch).c_str(),is1D[iV]?"":"colz");

    
    if (!is1D[iV]) std::cout << varLabel[iV] << " " << h2D[iV]->GetEntries() << " " << h2Dvetos[iV]->GetEntries() << std::endl;
    else std::cout << varLabel[iV] << " " << h1D[iV]->GetEntries() << " " << h1Dvetos[iV]->GetEntries() << std::endl;

    lat.DrawLatexNDC(0.1,0.95,"Delphes 200PU");

    myc[iV]->Update();
    myc[iV]->Print(("PLOTS/GenReco_"+varLabel[iV]+"_vetos.pdf").c_str());
    */
    
  }//loop on vars

  outFile->Write();

  return 0;
  
}//main
