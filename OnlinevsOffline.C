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


#include "TDRStyle.h"


double PuppiJetMin25OfflineEt(double offline, double Eta){ 
  if (fabs(Eta)<2.4) return (offline-5.08672)/1.68317; 
  else return (offline-12.7425)/1.77404;
};


double PuppiMETOfflineEt(const double & offline){
  return (offline-16.4021)/1.05826;
};

int OnlinevsOffline(){

  SetTdrStyle();
  gStyle->SetPadRightMargin(0.05);

  TCanvas *myc[3];
  gStyle->SetOptStat(0);

  std::ostringstream label;

  for (unsigned iV(0); iV<3; ++iV){
    label.str("");
    label << "myc_" << iV;
    myc[iV] = new TCanvas(label.str().c_str(),label.str().c_str(),1);
  }
  
  const unsigned nPt = 1000;

  TGraph *hJetEtBarrel = new TGraph();
  hJetEtBarrel->SetName("hJetEtBarrel");
  hJetEtBarrel->SetTitle(";offline E_{T} (GeV);online E_{T} (GeV)");

  TGraph *hJetEtEndcap = new TGraph();
  hJetEtEndcap->SetName("hJetEtEndcap");
  hJetEtEndcap->SetTitle(";offline E_{T} (GeV);online E_{T} (GeV)");

  TGraph *hMet = new TGraph();
  hMet->SetName("hMet");
  hMet->SetTitle(";offline MET (GeV); online MET (GeV)");

  for (unsigned iP(0); iP<nPt; ++iP){
    double step = 300./nPt;
    double et = iP*step;
    if (et<25) continue;
    hJetEtBarrel->SetPoint(iP,et,PuppiJetMin25OfflineEt(et,0));
    hJetEtEndcap->SetPoint(iP,et,PuppiJetMin25OfflineEt(et,3.0));
    hMet->SetPoint(iP,et,PuppiMETOfflineEt(et));
    //std::cout << iP << " " << et << " " << PuppiJetMin25OfflineEt(et,0) << " " << PuppiJetMin25OfflineEt(et,3.) << " " << PuppiMETOfflineEt(et) << std::endl;

  }

  TLatex lat;
  
  myc[0]->cd();
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  hJetEtBarrel->SetLineColor(1);
  hJetEtBarrel->SetLineColor(1);
  hJetEtBarrel->SetMarkerColor(1);
  hJetEtBarrel->SetMarkerStyle(6);
  hJetEtBarrel->Draw("AP");

  lat.DrawLatexNDC(0.2,0.8,"Jet |#eta| < 2.4");
  lat.DrawLatexNDC(0.2,0.9,"L1 = (offline-5.08672)/1.68317");

  myc[0]->Update();
  myc[0]->Print("PLOTS/L1vsOffline_barrelJet.pdf");
  
  myc[1]->cd();
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  hJetEtEndcap->SetMarkerStyle(6);
  hJetEtEndcap->Draw("AP");
  lat.DrawLatexNDC(0.2,0.8,"Jet |#eta| #geq 2.4");
  lat.DrawLatexNDC(0.2,0.9,"L1 = (offline-12.7425)/1.77404");

  myc[1]->Update();
  myc[1]->Print("PLOTS/L1vsOffline_endcapJet.pdf");

  myc[2]->cd();
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  hMet->SetMarkerStyle(6);
  hMet->Draw("AP");
  lat.DrawLatexNDC(0.2,0.9,"L1 = (offline-16.4021)/1.05826");

  myc[2]->Update();
  myc[2]->Print("PLOTS/L1vsOffline_met.pdf");

  
  return 0;
};
