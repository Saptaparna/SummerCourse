#include <TRandom.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <TGraphAsymmErrors.h>
#include <TVector2.h>
#include <TF1.h>
#include <iomanip>

using namespace std;

int makeFakeDataForDY(std::string outTree, std::string outfile)
{

 TRandom *r1 = new TRandom();
 TH2F *h_response=new TH2F("h_response", "Correlation Plot; Reconstructed Mass [GeV]; Generator Mass [GeV]", 100, 50.0, 150.0, 100, 50.0, 150.0); h_response->Sumw2();
 TH1F *h_InvariantMass_MuMu=new TH1F("h_InvariantMass_MuMu", "Di-muon invariant mass; m_{#mu#mu} [GeV]; Events/GeV", 200, 0, 200); h_InvariantMass_MuMu->Sumw2();
 TH1F *h_InvariantMass_MuMuGen = new TH1F("h_InvariantMass_MuMuGen", "Z generator particle mass; m_{Z} [GeV]; Events/GeV ", 100, 50, 150); h_InvariantMass_MuMuGen->Sumw2();

 TFile *outputFile;
 TTree *outputTree;
 double Mll_gen, Mll;
 
 std::string outputtreename=(outTree+".root").c_str();
 outputFile = new TFile((outputtreename).c_str(),"RECREATE");
 outputTree=new TTree("MLL_Tree", "MLL_Tree");
 outputTree->Branch("Mll", &Mll, "Mll/D");
 outputTree->Branch("Mll_gen", &Mll_gen, "Mll_gen/D");

 for(int i=0; i<100000; i++)
 {
   double genMass = r1->Uniform(50, 150);
   double recoMass = genMass + gRandom->Gaus(0, 0.1*genMass);
   h_InvariantMass_MuMu->Fill(recoMass);
   h_InvariantMass_MuMuGen->Fill(genMass);
   h_response->Fill(recoMass, genMass);
   Mll_gen = genMass;
   Mll = recoMass;
   outputTree->Fill();
 } 

 outputTree->Write();

 std::string histfilename=(outfile+".root").c_str();
 TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
 h_InvariantMass_MuMuGen->Write();
 h_InvariantMass_MuMu->Write();
 h_response->Write();
 tFile->Close();
 std::cout<<"Wrote output file "<<histfilename<<std::endl;

 return 0;

}
