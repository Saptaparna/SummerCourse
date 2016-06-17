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
#include <TVector3.h>
#include <TRandom.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TEllipse.h>

using std::string;

int makeBBplotDY()
{

  TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
  c1->SetLogx(); 
  
  //Double_t metricComb[16] = {0.995, 2.381, 4.687, 3.275, 5.481, 6.680, 7.593, 7.375, 7.769, 8.279, 8.336, 7.722, 7.099, 6.285, 2.568, 0.000};
  //Double_t binning[16] = {0.5, 1, 2, 2.5, 2.5, 3, 4, 5, 6, 7, 8, 9, 10, 25, 50, 100};
  Double_t metricComb[15] = {0.995, 2.381, 4.687, 5.481, 6.680, 7.593, 7.375, 7.769, 8.279, 8.336, 7.722, 7.099, 6.285, 2.568, 0.000};
  Double_t binning[15] = {0.5, 1, 2, 2.5, 3, 4, 5, 6, 7, 8, 9, 10, 25, 50, 100}; 
  Double_t metricComb_BB[15] = {0.995, 2.381, 4.687, 3.275, 6.680, 7.593, 7.375, 7.769, 8.279, 8.336, 7.722, 7.099, 6.285, 2.568, 0.000};
  Double_t binning_BB[15] = {0.5, 1, 2, 2.5, 3, 4, 5, 6, 7, 8, 9, 10, 25, 50, 100};

  TGraph *gr = new TGraph(15, binning, metricComb);
  gr->SetTitle("Combined Metric");
  gr->SetLineColor(kBlue);
  gr->SetLineWidth(3);
  gr->GetYaxis()->CenterTitle();
  gr->GetXaxis()->CenterTitle();
  gr->GetYaxis()->SetTitle("Combined Metric");
  gr->GetXaxis()->SetTitle("Binning");
  gr->Draw("AC*");
  
  TGraph *gr_bb = new TGraph(15, binning_BB, metricComb_BB);
  gr_bb->SetLineColor(kRed);
  gr_bb->SetLineWidth(2);
  gr_bb->SetMarkerStyle(20);
  gr_bb->SetMarkerColor(kBlack);
  gr_bb->SetMarkerSize(1.3);
  gr_bb->Draw("C*");

  c1->SaveAs("Combined_Metric.pdf");
  c1->SaveAs("Combined_Metric.png");
  return 0;
}
