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

int plotJetPurity(){//main

  const bool doTightVBF = true;
  
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
  
  const unsigned nV = 4;
  std::string varName[nV] = {
    "GenJet1_pt",
    "GenJet2_pt",
    "GenJet1_eta",
    "GenJet2_eta"
  };

  std::string varLabel[nV] = {
    "Jet1pt",
    "Jet2pt",
    "Jet1eta",
    "Jet2eta"
  };

  const unsigned nBins = doTightVBF?40:80;
  const double binMin[nV] = {0,0,-5,-5};
  const double binMax[nV] = {400,400,5,5};
  std::string axisTitles[nV] = {
    ";p_{T}^{j1} (GeV);Events",
    ";p_{T}^{j2} (GeV);Events",
    ";#eta^{j1};Events",
    ";#eta^{j2};Events"
  };

  std::string baseline = "";
  if (doTightVBF) baseline = "Gen_Mjj>1000 && Gen_detajj>4 && Gen_dphijj < 2 && GenJet1_pt>70 && GenJet2_pt>40 && GenJet1_eta*GenJet2_eta<0 && TMath::Abs(GenJet1_eta) < 5. && TMath::Abs(GenJet2_eta) < 5.";
  else baseline = "GenJet1_pt>0 && GenJet2_pt>0 && TMath::Abs(GenJet1_eta) < 5. && TMath::Abs(GenJet2_eta) < 5.";


  
  TCanvas *myc[nV];
  TCanvas *mycEff[nV];
  gStyle->SetOptStat(0);

  TFile *outFile = TFile::Open(doTightVBF?"HistosPurityTightVBFSel.root":"HistosPurity.root","RECREATE");
  outFile->cd();

  const unsigned nT = 4;
  std::string binsel[nV][nT];
  std::string l1sel[nV];
  for (unsigned iV(0); iV<nV; ++iV){
    if (iV%2==0) l1sel[iV] = "GenJet1_l1dR<0.4";
    else l1sel[iV] = "GenJet2_l1dR<0.4";
  }

  binsel[0][0] = "TMath::Abs(GenJet1_eta)<1";
  binsel[0][1] = "TMath::Abs(GenJet1_eta)>=1 && TMath::Abs(GenJet1_eta)<3";
  binsel[0][2] = "TMath::Abs(GenJet1_eta)>=3 && TMath::Abs(GenJet1_eta)<3.5";
  binsel[0][3] = "TMath::Abs(GenJet1_eta)>=3.5 && TMath::Abs(GenJet1_eta)<5";

  binsel[1][0] = "TMath::Abs(GenJet2_eta)<1";
  binsel[1][1] = "TMath::Abs(GenJet2_eta)>=1 && TMath::Abs(GenJet2_eta)<3";
  binsel[1][2] = "TMath::Abs(GenJet2_eta)>=3 && TMath::Abs(GenJet2_eta)<3.5";
  binsel[1][3] = "TMath::Abs(GenJet2_eta)>=3.5 && TMath::Abs(GenJet2_eta)<5";

  binsel[2][0] = "GenJet1_pt<50";
  binsel[2][1] = "GenJet1_pt>=50 && GenJet1_pt<100";
  binsel[2][2] = "GenJet1_pt>=100 && GenJet1_pt<150";
  binsel[2][3] = "GenJet1_pt>=150";

  binsel[3][0] = "GenJet2_pt<50";
  binsel[3][1] = "GenJet2_pt>=50 && GenJet2_pt<100";
  binsel[3][2] = "GenJet2_pt>=100 && GenJet2_pt<150";
  binsel[3][3] = "GenJet2_pt>=150";

  
  TH1F *h1Dpass[nT][nV];
  TH1F *h1Dtotal[nT][nV];

  TGraphAsymmErrors *gr[nT][nV];
  
  std::ostringstream label;

  TLatex lat;
  lat.SetTextSize(0.035);
  
  for (unsigned iV(0); iV<nV; ++iV){
    label.str("");
    label << "myc_" << varLabel[iV];
    myc[iV] = new TCanvas(label.str().c_str(),label.str().c_str(),1);
    myc[iV]->Divide(2,2);
    label.str("");
    label << "mycEff_" << varLabel[iV];
    mycEff[iV] = new TCanvas(label.str().c_str(),label.str().c_str(),1);
  }
  
  for (unsigned iV(0); iV<nV; ++iV){
    for (unsigned iT(0); iT<nT; ++iT){
      label.str("");
      label << "h1Dtotal_" << varLabel[iV] << "_" << iT;
      h1Dtotal[iT][iV] = new TH1F(label.str().c_str(),axisTitles[iV].c_str(),nBins,binMin[iV],binMax[iV]);
      h1Dtotal[iT][iV]->SetLineColor(1);
      h1Dtotal[iT][iV]->SetLineWidth(2);
      
      mycEff[iV]->cd();
      std::string varSel = varName[iV];
      tree->Draw((varSel+">>"+label.str()).c_str(),(baseline+" && "+binsel[iV][iT]).c_str());
    }
  }
  
  for (unsigned iV(0); iV<nV; ++iV){
    for (unsigned iT(0); iT<nT; ++iT){
      
      myc[iV]->cd(iT+1);
      label.str("");
      label << "h1Dpass_" << varLabel[iV] << "_" << iT;
      h1Dpass[iT][iV] = new TH1F(label.str().c_str(),axisTitles[iV].c_str(),nBins,binMin[iV],binMax[iV]);
      h1Dpass[iT][iV]->SetLineColor(2);
      h1Dpass[iT][iV]->SetLineWidth(2);
      
      std::string varSel = varName[iV];
      tree->Draw((varSel+">>"+label.str()).c_str(),(baseline+" && "+binsel[iV][iT]+" && "+l1sel[iV]).c_str());
      
      myc[iV]->cd(iT+1);
      gPad->SetGridx(1);
      gPad->SetGridy(1);
      h1Dtotal[iT][iV]->Draw();
      h1Dpass[iT][iV]->Draw("same");

      lat.SetTextColor(1);
      lat.DrawLatexNDC(0.1,0.96,binsel[iV][iT].c_str());
      
      gr[iT][iV] = new TGraphAsymmErrors(h1Dpass[iT][iV],h1Dtotal[iT][iV],"b");
      mycEff[iV]->cd();
      gPad->SetGridx(1);
      gPad->SetGridy(1);
      gr[iT][iV]->SetMarkerStyle(20+iT);
      gr[iT][iV]->SetMarkerColor(1+iT);
      gr[iT][iV]->SetLineColor(1+iT);
      gr[iT][iV]->SetLineWidth(2);
      gr[iT][iV]->SetTitle(axisTitles[iV].c_str());
      gr[iT][iV]->GetYaxis()->SetTitle("pass/total");
      gr[iT][iV]->Draw(iT==0?"APE":"PEsame");
      gr[iT][iV]->GetYaxis()->SetRangeUser(0,1.1);

      lat.SetTextColor(iT+1);
      lat.DrawLatexNDC(0.16,0.2+0.05*iT,binsel[iV][iT].c_str());

      std::cout << varLabel[iV] << " " << binsel[iV][iT] << " " << l1sel[iV] << " " << h1Dtotal[iT][iV]->GetEntries() << " " << h1Dpass[iT][iV]->GetEntries() << std::endl;
    }//loop on trigs
    
    myc[iV]->Update();
    if (doTightVBF) myc[iV]->Print(("PLOTS/TightVBFSel_GenL1_"+varLabel[iV]+".pdf").c_str());
    else myc[iV]->Print(("PLOTS/GenL1_"+varLabel[iV]+".pdf").c_str());
    mycEff[iV]->cd();
    
    mycEff[iV]->Update();
    if (doTightVBF) mycEff[iV]->Print(("PLOTS/TightVBFSel_GenL1MatchEff_"+varLabel[iV]+".pdf").c_str());
    else mycEff[iV]->Print(("PLOTS/GenL1MatchEff_"+varLabel[iV]+".pdf").c_str());
    

  }//loop on vars

  outFile->Write();

  return 0;
  
}//main
