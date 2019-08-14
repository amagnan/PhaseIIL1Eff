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

int plotTrigEff(){//main

  const bool dovsL1 = false;
  const bool doMatching = false;
  
  SetTdrStyle();
  gStyle->SetPadRightMargin(0.05);

  std::string finName = "L1BitStudyRecluster.root";
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
  
  const unsigned nV = 10;
  std::string varName[nV] = {
    "GenJet1_pt",
    "GenJet2_pt",
    "GenJet1_eta",
    "GenJet2_eta",
    "Gen_Mjj",
    "GenJet1_l1pt",
    "GenJet2_l1pt",
    "Gen_detajj",
    "Gen_dphijj",
    "GenMet"
  };
  std::string l1varName[nV] = {
    "L1JetLead_pt",
    "L1JetSublead_pt",
    "L1JetLead_eta",
    "L1JetSublead_eta",
    "L1_Mjj",
    "L1_detajj",
    "L1_dphijj",
    "L1Met"
  };
  std::string varLabel[nV] = {
    "Jet1pt",
    "Jet2pt",
    "Jet1eta",
    "Jet2eta",
    "Mjj",
    "L1Jet1pt",
    "L1Jet2pt",
    "detajj",
    "dphijj",
    "met"
  };

  const unsigned nT = 5;
  std::string trigbit[nT] = {
    "passL1MET==1",
    "passNN==1",
    //"passL1DoublePFJet160==1",
    "passL1DoublePFJet75==1",
    //"passMjj620==1",
    "passMjj500==1"
    ,"passMjj500==1 && passL1MET==1"
  };
  std::string trigbitLabel[nT] = {
    "passL1MET",
    "passNN",
    //"passL1DoublePFJet160",
    "passL1DoublePFJet75",
    //"passMjj620",
    "passMjj500"
    ,"passMjj500L1MET"
  };

  int color[8] = {1,2,3,4,6,7,8,9};
  
  const unsigned nBins = 40;
  //const double binMin[nV] = {0,0,-5,-5,0,0,0,0};
  //const double binMax[nV] = {400,400,5,5,4000,10,3.1416,400};
  const double binMin[nV] = {0,0,-5,-5,0,0,0,0,0,0};
  const double binMax[nV] = {400,400,5,5,4000,400,400,10,3.1416,400};
  std::string axisTitles[nV] = {
    ";p_{T}^{j1} (GeV);Events",
    ";p_{T}^{j2} (GeV);Events",
    ";#eta^{j1};Events",
    ";#eta^{j2};Events",
    ";M_{jj} (GeV);Events",
    ";p_{T}^{L1j1} (GeV);Events",
    ";p_{T}^{L1j2} (GeV);Events",
    ";#Delta#eta_{jj};Events" ,
    ";#Delta#phi_{jj};Events",
    ";MET (GeV);Events"
  };

  std::string baseline = dovsL1 ? "L1JetLead_pt>50 && L1JetSublead_pt>40" : "GenMet>150 && Gen_Mjj>1000 && Gen_detajj>4 && Gen_dphijj < 2. && GenJet1_pt>70 && GenJet2_pt>40 && GenJet1_eta*GenJet2_eta<0 && TMath::Abs(GenJet1_eta) < 5. && TMath::Abs(GenJet2_eta) < 5.";
  if (doMatching && !dovsL1) baseline += " && GenJet1_l1dR<0.4 && GenJet2_l1dR<0.4";
  //std::string baseline = "L1JetLead_pt>50 && L1JetSublead_pt>40";
  //std::string baseline = "GenJet1_pt>70 && GenJet2_pt>40 && TMath::Abs(GenJet1_eta) < 5. && TMath::Abs(GenJet2_eta) < 5.";// && nL1Jets > 1 && GenJet1_l1dR<0.4 && GenJet2_l1dR<0.4";

  TCanvas *myc[nV];
  TCanvas *mycEff[nV];
  gStyle->SetOptStat(0);

  TFile *outFile = TFile::Open(dovsL1? "HistosTrigEffvsL1.root" : doMatching?"HistosTrigEffMatchvsGEN.root":"HistosTrigEffNoMatchvsGEN_metcut.root","RECREATE");
  outFile->cd();
  
  TH1F *h1Dpass[nT][nV];
  TH1F *h1Dtotal[nV];

  TGraphAsymmErrors *gr[nT][nV];
  
  std::ostringstream label;

  TLatex lat;
  lat.SetTextSize(0.035);
  
  for (unsigned iV(0); iV<nV; ++iV){
    label.str("");
    label << "myc_" << varLabel[iV];
    myc[iV] = new TCanvas(label.str().c_str(),label.str().c_str(),1);
    myc[iV]->Divide(2,nT%2==0?nT/2:nT/2+1);
    label.str("");
    label << "mycEff_" << varLabel[iV];
    mycEff[iV] = new TCanvas(label.str().c_str(),label.str().c_str(),1);
  }
  
  for (unsigned iV(0); iV<nV; ++iV){
    label.str("");
    label << "h1Dtotal_" << varLabel[iV];
    h1Dtotal[iV] = new TH1F(label.str().c_str(),axisTitles[iV].c_str(),nBins,binMin[iV],binMax[iV]);
    h1Dtotal[iV]->SetLineColor(1);
    h1Dtotal[iV]->SetLineWidth(2);

    mycEff[iV]->cd();
    std::string varSel = dovsL1?l1varName[iV] : varName[iV];
    tree->Draw((varSel+">>"+label.str()).c_str(),baseline.c_str());
  }
  
  for (unsigned iV(0); iV<nV; ++iV){
    for (unsigned iT(0); iT<nT; ++iT){

      myc[iV]->cd(iT+1);
      label.str("");
      label << "h1Dpass_" << varLabel[iV] << "_" << trigbitLabel[iT];
      h1Dpass[iT][iV] = new TH1F(label.str().c_str(),axisTitles[iV].c_str(),nBins,binMin[iV],binMax[iV]);
      h1Dpass[iT][iV]->SetLineColor(2);
      h1Dpass[iT][iV]->SetLineWidth(2);
      
      std::string varSel = dovsL1?l1varName[iV] : varName[iV];
      tree->Draw((varSel+">>"+label.str()).c_str(),(baseline+" && "+trigbit[iT]).c_str());
      
      myc[iV]->cd(iT+1);
      gPad->SetGridx(1);
      gPad->SetGridy(1);
      h1Dtotal[iV]->Draw();
      h1Dpass[iT][iV]->Draw("same");

      lat.SetTextColor(1);
      lat.DrawLatexNDC(0.2,0.96,trigbitLabel[iT].c_str());
      
      gr[iT][iV] = new TGraphAsymmErrors(h1Dpass[iT][iV],h1Dtotal[iV],"b");
      mycEff[iV]->cd();
      gPad->SetGridx(1);
      gPad->SetGridy(1);
      gr[iT][iV]->SetMarkerStyle(20+iT);
      gr[iT][iV]->SetMarkerColor(color[iT]);
      gr[iT][iV]->SetLineColor(color[iT]);
      gr[iT][iV]->SetLineWidth(2);
      gr[iT][iV]->SetTitle(axisTitles[iV].c_str());
      gr[iT][iV]->GetYaxis()->SetTitle("pass/total");
      gr[iT][iV]->Draw(iT==0?"APE":"PEsame");
      gr[iT][iV]->GetYaxis()->SetRangeUser(0,1.1);

      lat.SetTextColor(color[iT]);
      lat.DrawLatexNDC(0.4,0.2+0.04*iT,trigbitLabel[iT].c_str());

      std::cout << varLabel[iV] << " " << trigbitLabel[iT] << " " << h1Dtotal[iV]->GetEntries() << " " << h1Dpass[iT][iV]->GetEntries() << std::endl;
    }//loop on trigs
    
    myc[iV]->Update();
    if (dovsL1) myc[iV]->Print(("PLOTS/TrigL1_"+varLabel[iV]+".pdf").c_str());
    else {
      if (doMatching) myc[iV]->Print(("PLOTS/TrigMatch_"+varLabel[iV]+".pdf").c_str());
      else myc[iV]->Print(("PLOTS/TrigNoMatch_metcut_"+varLabel[iV]+".pdf").c_str());
    }
    mycEff[iV]->cd();
    
    mycEff[iV]->Update();
    if (dovsL1) mycEff[iV]->Print(("PLOTS/TrigEffL1_"+varLabel[iV]+".pdf").c_str());
    else {
      if (doMatching) mycEff[iV]->Print(("PLOTS/TrigEffMatch_"+varLabel[iV]+".pdf").c_str());
      else mycEff[iV]->Print(("PLOTS/TrigEffNoMatch_metcut_"+varLabel[iV]+".pdf").c_str());
    }
    

  }//loop on vars

  outFile->Write();

  return 0;
  
}//main
