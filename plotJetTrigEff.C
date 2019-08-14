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

int plotJetTrigEff(){//main

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
  
  const unsigned nV = 2;
  std::string varName[nV] = {
    "GenJet1_pt",
    "GenJet1_eta"
  };
  std::string varLabel[nV] = {
    "Jet1pt",
    "Jet1eta"
  };

  const unsigned nT = 3;
  std::string trigbit[nT] = {
    "passL1Jet35",
    "passL1Jet75",
    "passL1Jet160"
  };
  
  const unsigned nBins = 40;
  const double binMin[nV] = {0,-5};
  const double binMax[nV] = {400,5};
  std::string axisTitles[nV] = {
    ";p_{T}^{j1} (GeV);Events",
    ";#eta^{j1};Events"
  };

  std::string baseline = "GenJet1_l1dR<0.4";// && Gen_Mjj>1000 && Gen_detajj>4 && Gen_dphijj<2. && GenJet1_pt>70 && GenJet2_pt>40 && GenJet1_eta*GenJet2_eta<0";


  TCanvas *myc[nV];
  TCanvas *mycEff[nV];
  gStyle->SetOptStat(0);

  TFile *outFile = TFile::Open("HistosJetTrigEff.root","RECREATE");
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
    myc[iV]->Divide(2,2);
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
    tree->Draw((varName[iV]+">>"+label.str()).c_str(),baseline.c_str());
  }
  
  for (unsigned iV(0); iV<nV; ++iV){
    for (unsigned iT(0); iT<nT; ++iT){

      myc[iV]->cd(iT+1);
      label.str("");
      label << "h1Dpass_" << varLabel[iV] << "_" << trigbit[iT];
      h1Dpass[iT][iV] = new TH1F(label.str().c_str(),axisTitles[iV].c_str(),nBins,binMin[iV],binMax[iV]);
      h1Dpass[iT][iV]->SetLineColor(2);
      h1Dpass[iT][iV]->SetLineWidth(2);
      
      tree->Draw((varName[iV]+">>"+label.str()).c_str(),(baseline+" && "+trigbit[iT]+" == 1").c_str());
      
      myc[iV]->cd(iT+1);
      gPad->SetGridx(1);
      gPad->SetGridy(1);
      h1Dtotal[iV]->Draw();
      h1Dpass[iT][iV]->Draw("same");

      lat.SetTextColor(1);
      lat.DrawLatexNDC(0.2,0.96,trigbit[iT].c_str());
      
      gr[iT][iV] = new TGraphAsymmErrors(h1Dpass[iT][iV],h1Dtotal[iV],"b");
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
      lat.DrawLatexNDC(0.7,0.2+0.08*iT,trigbit[iT].c_str());

      std::cout << varLabel[iV] << " " << trigbit[iT] << " " << h1Dtotal[iV]->GetEntries() << " " << h1Dpass[iT][iV]->GetEntries() << std::endl;
    }//loop on trigs

    myc[iV]->cd();
    myc[iV]->Update();
    myc[iV]->Print(("PLOTS/JetTrig_"+varLabel[iV]+".pdf").c_str());

    mycEff[iV]->cd();
    if (iV==0){
      TLine *l1 = new TLine(35,0,35,1.1);
      l1->SetLineColor(1);
      l1->SetLineWidth(2);
      l1->Draw("same");
      TLine *l2 = new TLine(75,0,75,1.1);
      l2->SetLineColor(2);
      l2->SetLineWidth(2);
      l2->Draw("same");
      TLine *l3 = new TLine(160,0,160,1.1);
      l3->SetLineColor(3);
      l3->SetLineWidth(2);
      l3->Draw("same");
    }
    
    
    mycEff[iV]->Update();
    mycEff[iV]->Print(("PLOTS/JetTrigEff_"+varLabel[iV]+".pdf").c_str());
    

  }//loop on vars

  outFile->Write();

  return 0;
  
}//main
