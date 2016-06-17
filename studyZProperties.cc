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
#include <TMatrixD.h>

using namespace std;

TLorentzVector fillTLorentzVector(double pT, double eta, double phi, double M)
{
  TLorentzVector object_p4;
  object_p4.SetPtEtaPhiM(pT, eta, phi, M);
  return object_p4;
}

int makeDiffPlot()
{
  gROOT->SetStyle("Plain");
  TCanvas *c1 = new TCanvas("c1_fit1","The Fit Canvas",200,10,700,500);
  c1->SetGridx();
  c1->SetGridy();
  c1->GetFrame()->SetFillColor(21);
  c1->GetFrame()->SetBorderMode(-1);
  c1->GetFrame()->SetBorderSize(5);

  TFile* file = TFile::Open("output_LowPtSUSY_Tree_Delphes_ZLL_10_v2.root");
  TH1F *h_InvMass = (TH1F*) file->Get("h_InvariantMass_MuMu");
  h_InvMass->SetLineColor(kBlue);
  h_InvMass->SetLineWidth(2);
  h_InvMass->Rebin(900);
  h_InvMass->GetYaxis()->SetTitle("Events/GeV");
  h_InvMass->GetYaxis()->SetTitleOffset(1.4);
  //h_InvMass->GetXaxis()->SetRangeUser(60, 120);
  double integral = h_InvMass->Integral();
  //h_InvMass->Scale(1.0/integral);
  h_InvMass->SetMaximum(60000.0);
  h_InvMass->Draw("hist");  

  TH1F *h_InvMassGen = (TH1F*) file->Get("h_InvariantMass_MuMuGen");
  h_InvMassGen->SetLineColor(kRed);
  h_InvMassGen->SetLineWidth(2);
  h_InvMassGen->Rebin(900);
  //h_InvMassGen->Scale(1.0/integral);
  h_InvMassGen->Draw("same hist");

  c1->SaveAs("InvariantMassDifference.pdf");
  c1->SaveAs("InvariantMassDifference.png");

  return 0;
}

int computeResponseMatrix()
{

  TFile* file = TFile::Open("output_LowPtSUSY_Tree_Delphes_ZLL_10_v2.root");
  TH1F *h_InvMass = (TH1F*) file->Get("h_InvariantMass_MuMu");
  h_InvMass->Rebin(900);
  std::cout << "h_InvMass->GetNbinsX() = " << h_InvMass->GetNbinsX() << std::endl;
  int nbinsReco = h_InvMass->GetNbinsX();  

  TH1F *h_InvMassGen = (TH1F*) file->Get("h_InvariantMass_MuMuGen");
  h_InvMassGen->Rebin(900);
  std::cout << "h_InvMassGen->GetNbinsX() = " << h_InvMassGen->GetNbinsX() << std::endl; 
  int nbinsGen = h_InvMassGen->GetNbinsX();

  TMatrixD response(nbinsReco, nbinsGen);
  for(int i=0; i<=nbinsGen; i++)
  {
    int j = i+1;
    std::cout << "h_InvMass->GetBinContent(i) = " << h_InvMass->GetBinContent(i) << std::endl;   
    std::cout << "h_InvMassGen->GetBinContent(i) = " << h_InvMassGen->GetBinContent(i) << std::endl;
    double responseRow_11 = ((double) h_InvMass->GetBinContent(i))/((double) h_InvMassGen->GetBinContent(i));
    std::cout << "responseRow_" << i << j << " = " << responseRow_11 << std::endl;    
    double responseRow_12 = ((double) h_InvMass->GetBinContent(j))/((double) h_InvMassGen->GetBinContent(i));
    std::cout << "responseRow_" << i << j << " = " << responseRow_12 << std::endl;
    double responseRow_21 = ((double) h_InvMass->GetBinContent(i))/((double) h_InvMassGen->GetBinContent(j));
    std::cout << "responseRow_" << i << j << " = " << responseRow_21 << std::endl;
    double responseRow_22 = ((double) h_InvMass->GetBinContent(j))/((double) h_InvMassGen->GetBinContent(j));
    std::cout << "responseRow_" << i << j << " = " << responseRow_22 << std::endl;
  }
  return 0;
}
