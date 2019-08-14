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
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

//#include "boost/lexical_cast.hpp"
//#include "boost/program_options.hpp"
//#include "boost/bind.hpp"
//#include "boost/function.hpp"
//#include "boost/format.hpp"
//#include "boost/algorithm/string.hpp"
//#include <boost/assign/list_of.hpp>
//#include <boost/range/algorithm_ext/for_each.hpp>

#include "fastjet/ClusterSequence.hh"

using namespace fastjet;

bool CheckValue(ROOT::Internal::TTreeReaderValueBase& value) {
   if (value.GetSetupStatus() < 0) {
      std::cerr << "Error " << value.GetSetupStatus()
                << "setting up reader for " << value.GetBranchName() << '\n';
      return false;
   }
   return true;
}

struct PFJet{
  double Et;
  double Eta;
  double Phi;
  int Bx;
};

double pairInvMass(const PFJet & jet1, const PFJet & jet2){
  return sqrt(2.0*jet1.Et*jet2.Et*(cosh(jet1.Eta-jet2.Eta)-cos(jet1.Phi-jet2.Phi)));
};

bool PuppiJetMin25OfflineEtCut(double offline, double Et, double Eta){ 
  if (fabs(Eta)<2.4) return Et>(offline-5.08672)/1.68317; 
  else return Et>(offline-12.7425)/1.77404;
};

double PuppiMETOfflineEt(const double & offline){
  return (offline-16.4021)/1.05826;
};

bool L1_PFMet(const double Et, const double & value){
  return Et>PuppiMETOfflineEt(value);
};

bool L1_jet(const PFJet & jet, const double value){
  return PuppiJetMin25OfflineEtCut(value,jet.Et,jet.Eta) &&  jet.Et>25.0 && jet.Bx==0;
};

bool L1_pair(const PFJet & jet1, const PFJet & jet2, const double & value){
  return pairInvMass(jet1,jet2)>value;
};

void L1_jet2(const PFJet & jet, const PFJet & jet1, const double pt1, const double & pt2, const double & mjj, unsigned & counter){
  if ( L1_jet(jet,pt2) && L1_jet(jet1,pt1) && L1_pair(jet,jet1,mjj) ) counter++;
  return;
};



int main(){//main

  unsigned nEvtsToProc = 0;// 0 = all

  const bool pMakeJets = true;

  TFile *outFile = TFile::Open("L1BitStudyRecluster.root","RECREATE");
  outFile->cd();
  TTree *outtree = new TTree("LightTree","Analysis light tree");
  unsigned nGenJet = 0;
  double GenJet1_pt = 0;
  double GenJet2_pt = 0;
  double GenJet1_eta = 0;
  double GenJet2_eta = 0;
  double GenJet1_phi = 0;
  double GenJet2_phi = 0;
  double Gen_Mjj = 0;
  double Gen_detajj = 0;
  double Gen_dphijj = 0;
  float GenMet = 0;
  double GenJet1_l1dR = 10;
  double GenJet2_l1dR = 10;
  double GenJet1_l1pt = 0;
  double GenJet2_l1pt = 0;
  double GenJet1_l1eta = 0;
  double GenJet2_l1eta = 0;
  double GenJet1_l1phi = 0;
  double GenJet2_l1phi = 0;

  unsigned nL1Jets = 0;
  double L1JetLead_pt = 0;
  double L1JetSublead_pt = 0;
  double L1Jet1_pt = 0;
  double L1Jet2_pt = 0;
  double L1JetLead_eta = 0;
  double L1JetSublead_eta = 0;
  double L1Jet1_eta = 0;
  double L1Jet2_eta = 0;
  double L1JetLead_phi = 0;
  double L1JetSublead_phi = 0;
  double L1Jet1_phi = 0;
  double L1Jet2_phi = 0;
  double L1_Mjj = 0;
  double L1_detajj = 0;
  double L1_dphijj = 0;
  double L1Met = 0;

  bool passL1Jet35 = false;
  bool passL1Jet75 = false;
  bool passL1Jet160 = false;
  bool passMET100 = false;
  bool passMET150 = false;
  bool passNN = false;
  bool passL1DoublePFJet160 = false;
  bool passL1DoublePFJet75 = false;
  bool passMjj620 = false;
  bool passMjj500 = false;

  
  outtree->Branch("nGenJet",&nGenJet);
  outtree->Branch("GenJet1_pt",&GenJet1_pt);
  outtree->Branch("GenJet2_pt",&GenJet2_pt);
  outtree->Branch("GenJet1_eta",&GenJet1_eta);
  outtree->Branch("GenJet2_eta",&GenJet2_eta);
  outtree->Branch("GenJet1_phi",&GenJet1_phi);
  outtree->Branch("GenJet2_phi",&GenJet2_phi);
  outtree->Branch("Gen_Mjj",&Gen_Mjj);
  outtree->Branch("Gen_detajj",&Gen_detajj);
  outtree->Branch("Gen_dphijj",&Gen_dphijj);
  outtree->Branch("GenMet",&GenMet);
  outtree->Branch("GenJet1_l1dR",&GenJet1_l1dR);
  outtree->Branch("GenJet2_l1dR",&GenJet2_l1dR);
  outtree->Branch("GenJet1_l1pt",&GenJet1_l1pt);
  outtree->Branch("GenJet2_l1pt",&GenJet2_l1pt);
  outtree->Branch("GenJet1_l1eta",&GenJet1_l1eta);
  outtree->Branch("GenJet2_l1eta",&GenJet2_l1eta);
  outtree->Branch("GenJet1_l1phi",&GenJet1_l1phi);
  outtree->Branch("GenJet2_l1phi",&GenJet2_l1phi);

  outtree->Branch("nL1Jets",&nL1Jets);
  outtree->Branch("L1JetLead_pt",&L1JetLead_pt);
  outtree->Branch("L1JetSublead_pt",&L1JetSublead_pt);
  outtree->Branch("L1Jet1_pt",&L1Jet1_pt);
  outtree->Branch("L1Jet2_pt",&L1Jet2_pt);
  outtree->Branch("L1JetLead_eta",&L1JetLead_eta);
  outtree->Branch("L1JetSublead_eta",&L1JetSublead_eta);
  outtree->Branch("L1Jet1_eta",&L1Jet1_eta);
  outtree->Branch("L1Jet2_eta",&L1Jet2_eta);
  outtree->Branch("L1JetLead_phi",&L1JetLead_phi);
  outtree->Branch("L1JetSublead_phi",&L1JetSublead_phi);
  outtree->Branch("L1Jet1_phi",&L1Jet1_phi);
  outtree->Branch("L1Jet2_phi",&L1Jet2_phi);
  outtree->Branch("L1_Mjj",&L1_Mjj);
  outtree->Branch("L1_detajj",&L1_detajj);
  outtree->Branch("L1_dphijj",&L1_dphijj);
  outtree->Branch("L1Met",&L1Met);

  outtree->Branch("passL1Jet35",&passL1Jet35);
  outtree->Branch("passL1Jet75",&passL1Jet75);
  outtree->Branch("passL1Jet160",&passL1Jet160);
  outtree->Branch("passMET100",&passMET100);
  outtree->Branch("passMET150",&passMET150);
  outtree->Branch("passNN",&passNN);
  outtree->Branch("passL1DoublePFJet160",&passL1DoublePFJet160);
  outtree->Branch("passL1DoublePFJet75",&passL1DoublePFJet75);
  outtree->Branch("passMjj620",&passMjj620);
  outtree->Branch("passMjj500",&passMjj500);
  
  std::string finNameNN = "Hinv_rate_1.root";

  TFile *finNN = TFile::Open(finNameNN.c_str());
  if (!finNN) {
    std::cout << " Input file " << finNameNN << " not found." << std::endl;
    return 1;
  }
  else {
    std::cout << " File found: " << finNN->GetName() << std::endl;
  }
  finNN->cd();

  TTreeReader theNNReader("tree",gDirectory);
  TTreeReaderValue<float> NNval(theNNReader, "dense");

  std::string finName = pMakeJets ? "L1ntuple_VBFH_200PU.root" : "VBFH125_L1NtuplePhaseII_160.root";

  TFile *fin = TFile::Open(finName.c_str());
  if (!fin) {
    std::cout << " Input file " << finName << " not found." << std::endl;
    return 1;
  }
  else {
    std::cout << " File found: " << fin->GetName() << std::endl;
  }
  fin->cd("l1PhaseIITree");

  /*TTree *l1tree = (TTree*)gDirectory->Get("L1PhaseIITree");
  if (!l1tree){
    std::cout << " L1 tree not found !" << std::endl;
    return 1;
  }
  TTree *gentree = (TTree*)gDirectory->Get("L1GenTree");
  if (!gentree){
    std::cout << " Gen tree not found !" << std::endl;
    return 1;
    }*/


  TTreeReader theReader("L1PhaseIITree",gDirectory);
  TTreeReaderValue<UInt_t> nJets(theReader, "nPuppiJets");
  TTreeReaderValue<std::vector<double> > jetEt(theReader, "puppiJetEt");
  TTreeReaderValue<std::vector<double> > jetEta(theReader, "puppiJetEta");
  TTreeReaderValue<std::vector<double> > jetPhi(theReader, "puppiJetPhi");
  TTreeReaderValue<std::vector<int> > jetBx(theReader, "puppiJetBx");
  TTreeReaderValue<double> met(theReader, "puppiMETEt");
  TTreeReaderValue<double> metphi(theReader, "puppiMETPhi");
  
  fin->cd("genTree");
  TTreeReader theGenReader("L1GenTree",gDirectory);
  TTreeReaderValue<Int_t> nGenJets(theGenReader, "nJet");
  TTreeReaderValue<std::vector<float> > genjetPt(theGenReader, "jetPt");
  TTreeReaderValue<std::vector<float> > genjetEta(theGenReader, "jetEta");
  TTreeReaderValue<std::vector<float> > genjetPhi(theGenReader, "jetPhi");
  TTreeReaderValue<std::vector<float> > genjetMass(theGenReader, "jetM");
  TTreeReaderValue<Float_t> metTrue(theGenReader, "genMetTrue");

  std::vector<PseudoJet> lParticles;
  // choose a jet definition
  double R = 0.4;
  JetDefinition jet_def(antikt_algorithm, R);
  
  TTreeReaderValue<Int_t> nPart(theGenReader, "nPart");
  TTreeReaderValue<std::vector<int> > partId(theGenReader, "partId");
  TTreeReaderValue<std::vector<int> > partStat(theGenReader, "partStat");
  TTreeReaderValue<std::vector<float> > partPt(theGenReader, "partPt");
  TTreeReaderValue<std::vector<float> > partEta(theGenReader, "partEta");
  TTreeReaderValue<std::vector<float> > partPhi(theGenReader, "partPhi");
  TTreeReaderValue<std::vector<float> > partE(theGenReader, "partE");


  
  unsigned ievt = 0;
  unsigned nEvts = theReader.GetEntries(true);
  nEvtsToProc = (nEvtsToProc>0 && nEvtsToProc<nEvts) ? nEvtsToProc : nEvts;
  while(((pMakeJets && theNNReader.Next()) || !pMakeJets) && theReader.Next() && theGenReader.Next() && ievt < nEvtsToProc){
    //while(theReader.Next() && theGenReader.Next() && ievt < nEvtsToProc){
    if (ievt%1000==0) 
      std::cout << ".... Processing entry " << ievt << std::endl;

    //if (!metTrue.Get()) return 1;

    nGenJet = *nGenJets;

    if (nGenJet < 2) {
      ievt++;
      continue;
    }

    unsigned nJ = *nJets;
    nL1Jets = nJ;

    std::vector<PFJet> jetVec;
    jetVec.reserve(nJ);
    //std::cout << " - Event " << ievt << " nJets = " << nJ << " check " << (*jetEt).size() << std::endl;

    //std::cout << " -- met = " << *met << " ngenJets " << *nGenJets << std::endl;


    L1Met = *met;
    
    passMET100 = L1_PFMet(*met,100.);
    passMET150 = L1_PFMet(*met,150.);
    passNN = (*NNval)>0.85;

    //initialise variables...
    passL1DoublePFJet160 = false;
    passL1DoublePFJet75 = false;
    passMjj620 = false;
    passMjj500 = false;
    passL1Jet35 = false;
    passL1Jet75 = false;
    passL1Jet160 = false;
    
    GenJet1_l1dR = 10;
    GenJet2_l1dR = 10;
    GenJet1_l1pt = 0;
    GenJet2_l1pt = 0;
    GenJet1_l1eta = 0;
    GenJet2_l1eta = 0;
    GenJet1_l1phi = 0;
    GenJet2_l1phi = 0;
    L1JetLead_pt = 0;
    L1JetSublead_pt = 0;
    L1Jet1_pt = 0;
    L1Jet2_pt = 0;
    L1JetLead_eta = 0;
    L1JetSublead_eta = 0;
    L1Jet1_eta = 0;
    L1Jet2_eta = 0;
    L1JetLead_phi = 0;
    L1JetSublead_phi = 0;
    L1Jet1_phi = 0;
    L1Jet2_phi = 0;
    L1_Mjj = 0;
    L1_detajj = 0;
    L1_dphijj = 0;
    L1Met = 0;
    
    for (unsigned iJ(0); iJ<nJ; ++iJ){
      PFJet pfJet;
      pfJet.Et = (*jetEt)[iJ] ;
      pfJet.Eta = (*jetEta)[iJ] ;
      pfJet.Phi = (*jetPhi)[iJ] ;
      pfJet.Bx = (*jetBx)[iJ] ;

      if (iJ==0){
	passL1Jet35 = L1_jet(pfJet,35);
	passL1Jet75 = L1_jet(pfJet,75);
	passL1Jet160 = L1_jet(pfJet,160);
	L1JetLead_pt = pfJet.Et;
	L1JetLead_eta = pfJet.Eta;
	L1JetLead_phi = pfJet.Phi;
      }
      else if (iJ==1){
	L1JetSublead_pt = pfJet.Et;
	L1JetSublead_eta = pfJet.Eta;
	L1JetSublead_phi = pfJet.Phi;
      }
      
      //if (pass)
      jetVec.push_back(pfJet);
      //std::cout << " -- Jet " << iJ 
      //<< " pT = " << (*jetEt)[iJ] 
      //	<< " eta = " << (*jetEta)[iJ]
      //	<< " phi = " << (*jetPhi)[iJ]
      //	<< " bx = " << (*jetBx)[iJ]
      //	<< std::endl;
    }
    
    nJ = jetVec.size();
    unsigned jet1 = nJ;
    unsigned jet2 = nJ;
    double massSel = 0;
    for (unsigned iJ(0); iJ<nJ; ++iJ){      
      
      bool pass1 = L1_jet(jetVec[iJ],160);
      bool pass1met = L1_jet(jetVec[iJ],75);
      //if (!pass1 && !pass1met) continue;
      for (unsigned jJ(iJ+1); jJ<nJ; ++jJ){
	bool pass2 = L1_jet(jetVec[jJ],35);
	//if (!pass2) continue;
	double mass = pairInvMass(jetVec[iJ],jetVec[jJ]);
	//std::cout << " -- M" << iJ << jJ << " = " << mass << std::endl;
	passL1DoublePFJet160 = passL1DoublePFJet160 || (pass1 && pass2);
	passMjj620 = passMjj620 || (pass1 && pass2 && mass>620);
	passL1DoublePFJet75 = passL1DoublePFJet75 || (pass1met && pass2);
	passMjj500 = passMjj500 || (pass1met && pass2 && mass>500);

	if (mass>massSel){
	  massSel = mass;
	  jet1 = iJ;
	  jet2 = jJ;
	}
      }
    }
    
    //std::cout << " -- Selected highest mass pair with mass = " << massSel << " jet1 " << jet1 << " jet2 " << jet2 << std::endl;

    //unsigned n_pass1 = std::count_if(jetVec.begin(), jetVec.end(), bind(L1_jet,_1,160.0)) ;
    //jet1Vec.erase(std::remove_if(jet1Vec.begin(), jet1Vec.end(),  bind(L1_jet,_1, 160.)), jet1Vec.end());
    //boost::for_each(jet1Vec, jet2Vec, boost::bind(L1_jet2, _1, _2, 160.,35.,620,n_pass2));

    //std::cout << " - L1MET bit: " << passMET100 << std::endl
    //	      << " 2j triggerBit " << passL1DoublePFJetMassMin
    //	      << std::endl
    //	      << " 2j+met triggerBit " << passL1DoublePFJetMassMinPfMet
    //	      << std::endl;
    if (jet1 < nJ){
      L1Jet1_pt = jetVec[jet1].Et;
      L1Jet1_eta = jetVec[jet1].Eta;
      L1Jet1_phi = jetVec[jet1].Phi;
    }
    if (jet2 < nJ){
      L1Jet2_pt = jetVec[jet2].Et;
      L1Jet2_eta = jetVec[jet2].Eta;
      L1Jet2_phi = jetVec[jet2].Phi;
    }
    
    TLorentzVector l11;
    l11.SetPtEtaPhiM(L1Jet1_pt,L1Jet1_eta,L1Jet1_phi,0.135);
    TLorentzVector l12;
    l12.SetPtEtaPhiM(L1Jet2_pt,L1Jet2_eta,L1Jet2_phi,0.135);
    TLorentzVector l1Pair = l11 + l12;


    L1_Mjj = l1Pair.M();
    L1_detajj = fabs(L1Jet1_eta-L1Jet2_eta);
    L1_dphijj = fabs(l11.DeltaPhi(l12));

    TLorentzVector gen1;
    TLorentzVector gen2;

    if (pMakeJets){//pMakeJets
      lParticles.clear();
      
      for (unsigned iP(0); iP<(*nPart); ++iP){
	int status = (*partStat)[iP];
	if (status != 1) continue;
	int id = (*partId)[iP];
	if (abs(id)==12 || abs(id) == 14 || abs(id) == 16) continue;
	TLorentzVector lpart;
	lpart.SetPtEtaPhiE((*partPt)[iP],(*partEta)[iP],(*partPhi)[iP],(*partE)[iP]);
	
	lParticles.push_back( PseudoJet(lpart.Px(),lpart.Py(),lpart.Pz(),lpart.E()));
      }
	   
      // run the clustering, extract the jets
      ClusterSequence cs(lParticles, jet_def);
      std::vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
      
      // print the jets
      //std::cout <<   "-- evt " << ievt << ": found " << jets.size() << " Jets." << std::endl;
      //std::cout << " -- genjets size: " << nGenJet << std::endl;
      for (unsigned i = 0; i < jets.size(); i++) {
	const PseudoJet & lFastJet = jets[i];

	double phi = lFastJet.phi()>3.14159 ? lFastJet.phi()-2*3.14159 : lFastJet.phi();
	
	/*	std::cout << " -------- jet " << i << ": "
		  << lFastJet.perp() << " " 
		  << lFastJet.eta() << " " << phi
	  //<< " " << lFastJet.constituents().size()
		  << std::endl;
	// std::vector<PseudoJet> constituents = lFastJet.constituents();
	// for (unsigned j = 0; j < constituents.size(); j++) {
	//   std::cout << "    constituent " << j << "'s pt: " << constituents[j].perp()
	// 	      << std::endl;
	// }
	if (i<nGenJet){
	  std::cout << " -------- genjet " << i << ": "
		    << (*genjetPt)[i] << " " 
		    << (*genjetEta)[i] << " "
		    << (*genjetPhi)[i]
		    << std::endl;

		    }*/

	if (i==0) {
	  GenJet1_pt = lFastJet.perp();
	  GenJet1_eta = lFastJet.eta();
	  GenJet1_phi = phi;
	  gen1.SetPtEtaPhiE(GenJet1_pt,GenJet1_eta,GenJet1_phi,lFastJet.E());
	}
	else if (i==1){
	  GenJet2_pt = lFastJet.perp();
	  GenJet2_eta = lFastJet.eta();
	  GenJet2_phi = phi;
	  gen2.SetPtEtaPhiE(GenJet2_pt,GenJet2_eta,GenJet2_phi,lFastJet.E());
	}
	

      }//loop on fastjets
      
    }//pMakeJets
    else {
      GenJet1_pt = (*genjetPt)[0];
      GenJet2_pt = (*genjetPt)[1];
      GenJet1_eta = (*genjetEta)[0];
      GenJet2_eta = (*genjetEta)[1];
      GenJet1_phi = (*genjetPhi)[0];
      GenJet2_phi = (*genjetPhi)[1];
      
      gen1.SetPtEtaPhiM(GenJet1_pt,GenJet1_eta,GenJet1_phi,(*genjetMass)[0]);
      gen2.SetPtEtaPhiM(GenJet2_pt,GenJet2_eta,GenJet2_phi,(*genjetMass)[1]);
    }
    TLorentzVector genPair = gen1 + gen2;

    double dr1min = 10, dr2min=10;
    int idx1 = -1, idx2 = -1;
    for (unsigned iJ(0); iJ<nJ; ++iJ){      
      TLorentzVector lrec;
      lrec.SetPtEtaPhiM(jetVec[iJ].Et,jetVec[iJ].Eta,jetVec[iJ].Phi,0.135);
      double dR1 = gen1.DeltaR(lrec);
      double dR2 = gen2.DeltaR(lrec);
      if (dR1 < dr1min){
	dr1min = dR1;
	idx1 = iJ;
      }
      if (dR2 < dr2min){
	dr2min = dR2;
	idx2 = iJ;
      }
    }

    Gen_Mjj = genPair.M();
    Gen_detajj = fabs(GenJet1_eta-GenJet2_eta);
    Gen_dphijj = fabs(gen1.DeltaPhi(gen2));
    GenMet = *metTrue;

    if (idx1 >=0){
      GenJet1_l1dR = dr1min;
      GenJet1_l1pt = jetVec[idx1].Et;
      GenJet1_l1eta = jetVec[idx1].Eta;
      GenJet1_l1phi = jetVec[idx1].Phi;
    }
    if (idx2 >= 0){
      GenJet2_l1dR = dr2min;
      GenJet2_l1pt = jetVec[idx2].Et;
      GenJet2_l1eta = jetVec[idx2].Eta;
      GenJet2_l1phi = jetVec[idx2].Phi;
    }
    
    outtree->Fill();
    
    ievt++;
  }


  outFile->cd();
  outtree->Write();
  outFile->Write();
  
  std::cout << " -- File " << outFile->GetName() << " successfully written." << std::endl;

  return 0;

}//main
