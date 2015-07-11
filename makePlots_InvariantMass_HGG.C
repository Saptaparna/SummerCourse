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

using namespace std;

void makePlots_InvariantMass_HGG(){


  gROOT->SetStyle("Plain");
  TCanvas *c1 = new TCanvas("c1_fit1","The Fit Canvas",200,10,700,500);
  c1->SetGridx();
  c1->SetGridy();
  c1->GetFrame()->SetFillColor(21);
  c1->GetFrame()->SetBorderMode(-1);
  c1->GetFrame()->SetBorderSize(5); 
 
  TFile* file = TFile::Open("Output_LowPtSUSY_Tree_HGG.root");
  TH1F *h_InvMass = (TH1F*) file->Get("h_InvariantMass_PhPh");
  h_InvMass->SetLineColor(kBlue);
  h_InvMass->SetLineWidth(2);
  h_InvMass->Rebin(60);
  h_InvMass->GetYaxis()->SetTitle("Events/2.0 GeV");
  h_InvMass->GetYaxis()->SetTitleOffset(1.4);
  h_InvMass->GetXaxis()->SetRangeUser(100, 150); 
  TF1 *fSignal = new TF1("fSignal","gaus",100.0,140.0);
  fSignal->SetParLimits(0, 29000.0, 150000.0);
  fSignal->SetParLimits(1, 100.0, 140.0);
  fSignal->SetParLimits(2, 0.0, 1.0);
  h_InvMass->Draw("HIST");
  h_InvMass->Fit("fSignal", "", "", 10.0, 140.0);
  c1->SaveAs("h_InvariantMass_PhPh.pdf");
  c1->SaveAs("h_InvariantMass_PhPh.png");

} 
