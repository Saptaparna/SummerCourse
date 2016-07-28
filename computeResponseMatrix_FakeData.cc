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
#include <TUnfold.h>

using namespace std;

int computeResponseMatrix_FakeData(std::string infile, std::string outfile){

  std::string inputfilename=(infile+".root").c_str();
  TChain *tree=new TChain("MLL_Tree");
  tree->Add(inputfilename.c_str());
  std::cout<<"Opened input file "<<inputfilename<<std::endl;

  Double_t        Mll;
  Double_t        Mll_gen;

  tree->SetBranchAddress("Mll", &(Mll));
  tree->SetBranchAddress("Mll_gen", &(Mll_gen));

  TH2F *h_response=new TH2F("h_response", "Response Matrix; Reconstructed Mass [GeV]; Generated Mass [GeV]", 170, 30.0, 200.0, 170, 30.0, 200.0); h_response->Sumw2();
  //TH2F *h_response=new TH2F("h_response", "Response Matrix; Reconstructed Mass [GeV]; Generated Mass [GeV]", 34,  30.0, 200.0, 34, 30.0, 200.0); h_response->Sumw2();
  //TH2F *h_response=new TH2F("h_response", "Response Matrix; Reconstructed Mass [GeV]; Generated Mass [GeV]", 17,  30.0, 200.0, 17, 30.0, 200.0); h_response->Sumw2();
  TH1F *h_InvariantMass_MuMu=new TH1F("h_InvariantMass_MuMu", "Di-muon invariant mass; m_{#mu#mu} [GeV]; Events/GeV", 200, 0, 200); h_InvariantMass_MuMu->Sumw2();
  TH1F *h_InvariantMass_MuMuGen = new TH1F("h_InvariantMass_MuMuGen", "Z generator particle mass; m_{Z} [GeV]; Events/GeV ", 200, 0, 200); h_InvariantMass_MuMuGen->Sumw2();

  Float_t rebin_array[] = {31.7716571, 38.03552205, 41.20818079, 42.95655956, 45.43517219, 47.36059869, 49.15550893, 52.07252724, 55.66729689, 128.51950869, 132.30855742, 137.39460633, 142.17852438, 146.64154812, 152.19551123, 158.13295279, 164.09264608, 167.87511639, 172.84042825, 176.99907853, 184.54875878, 199.17566438};

  Int_t  nBins = sizeof(rebin_array)/sizeof(Float_t) - 1;
  std::cout << "nBins = " << nBins << std::endl;

  TH2F *h_responseBB = new TH2F("h_responseBB", "h_responseBB", nBins, rebin_array, nBins, rebin_array); h_responseBB->Sumw2();

  int nEvents=tree->GetEntries();
  std::cout << "nEvents= " << nEvents << std::endl;

  for (int i=0; i<nEvents ; ++i)
  {
    tree->GetEvent(i);
    double genMass = Mll_gen;
    double recoMass = Mll;    
    h_InvariantMass_MuMu->Fill(recoMass);
    h_InvariantMass_MuMuGen->Fill(genMass);
    h_response->Fill(recoMass, genMass);
    h_responseBB->Fill(recoMass, genMass);
  }//event loop closed

  TH2F *h_finalResponse=new TH2F("h_finalResponse", "Response Matrix; Reconstructed Mass [GeV]; Generator Mass [GeV]", 170, 30.0, 200.0, 170, 30.0, 200.0); h_finalResponse->Sumw2();
  //TH2F *h_finalResponse=new TH2F("h_finalResponse", "Response Matrix; Reconstructed Mass [GeV]; Generator Mass [GeV]", 34, 30.0, 200.0, 34, 30.0, 200.0); h_finalResponse->Sumw2();
  for(int k = 1; k<= h_response->GetNbinsX(); k++) //reco information
  { 
    double sumY = 0;
    for(int l = 1; l<= h_response->GetNbinsY(); l++)//gen information
    {
      sumY += h_response->GetBinContent(k, l);
    }
    
    for(int l = 1; l<= h_response->GetNbinsY(); l++)
    {
      if(sumY>0) h_finalResponse->SetBinContent(k, l, (h_response->GetBinContent(k, l)/sumY));
      else h_finalResponse->SetBinContent(k, l, 0);
      double fsrBinLow = h_response->GetYaxis()->FindBin(90); 
      //std::cout << "h_finalResponse->GetBinContent(k, l) = " << h_finalResponse->GetBinContent(k, l) << std::endl;
    }
  }

  double figureOfmerit = 1.0;
  for(int i = 1; i<= h_finalResponse->GetNbinsX(); i++) //without nested loop
  {
    double offDiagonal = h_finalResponse->GetBinContent(i, i+1);
    double diagonal_d = h_finalResponse->GetBinContent(i, i);
    double diagonal_r = h_finalResponse->GetBinContent(i+1, i+1); 
 
    double lowerDiagonal = h_finalResponse->GetBinContent(i, i-1);
    double lowerDiagonal_u = diagonal_d;
    double lowerDiagonal_l = h_finalResponse->GetBinContent(i-1, i-1);

    //if(offDiagonal == 0.0 or diagonal_d == 0.0 or diagonal_r == 0.0 or lowerDiagonal == 0.0 or lowerDiagonal_l == 0.0 or lowerDiagonal_u == 0.0) std::cout << "Some empty cells found" << std::endl;

    if(offDiagonal > 0.0 and diagonal_d > 0.0 and diagonal_r > 0.0)
    {
      figureOfmerit *= offDiagonal/(sqrt(diagonal_d*diagonal_r)); 
      std::cout << "offDiagonal = " << offDiagonal << " diagonal_d = " << diagonal_d << " diagonal_r = " << diagonal_r  << std::endl;
    } 
    if(lowerDiagonal > 0.0 and lowerDiagonal_l > 0.0 and lowerDiagonal_u > 0.0)
    {
      double figureOfmerit_lower = lowerDiagonal/(sqrt(lowerDiagonal_u*lowerDiagonal_l));
      figureOfmerit *= figureOfmerit_lower;
      std::cout << "lowerDiagonal = " << lowerDiagonal << " lowerDiagonal_u = " << lowerDiagonal_u << " lowerDiagonal_l = " << lowerDiagonal_l  << std::endl;
      //std::cout << "figureOfmerit_lower = " << figureOfmerit_lower << std::endl;
    }
    //std::cout << "figureOfmerit 1 = " << figureOfmerit << std::endl;
  }

  //std::cout << "figureOfmerit standard GeV binning = " << figureOfmerit/(h_finalResponse->GetNbinsX()) << std::endl;
  std::cout << "figureOfmerit standard GeV binning = " << figureOfmerit << std::endl;

  TH2F *h_finalResponseBB=new TH2F("h_finalResponseBB", "Response Matrix; Reconstructed Mass [GeV]; Generator Mass [GeV]", nBins, rebin_array, nBins, rebin_array); h_finalResponseBB->Sumw2();
  for(int k = 1; k<= h_responseBB->GetNbinsX(); k++) //reco information
  { 
    double sumY = 0;
    for(int l = 1; l<= h_responseBB->GetNbinsY(); l++)//gen information
    {
      sumY += h_responseBB->GetBinContent(k, l);
    }
    
    for(int l = 1; l<= h_responseBB->GetNbinsY(); l++)
    {
      if(sumY>0) h_finalResponseBB->SetBinContent(k, l, (h_responseBB->GetBinContent(k, l)/sumY));
      else h_finalResponseBB->SetBinContent(k, l, 0);
    }
  }

  double figureOfmeritBB = 1.0;
  for(int i = 1; i<= h_finalResponseBB->GetNbinsX(); i++) //without nested loop
  {
    double offDiagonal = h_finalResponseBB->GetBinContent(i, i+1);
    double diagonal_d = h_finalResponseBB->GetBinContent(i, i);
    double diagonal_r = h_finalResponseBB->GetBinContent(i+1, i+1);

    double lowerDiagonal = h_finalResponseBB->GetBinContent(i, i-1);
    double lowerDiagonal_u = diagonal_d;
    double lowerDiagonal_l = h_finalResponseBB->GetBinContent(i-1, i-1);

    if(offDiagonal > 0.0 and diagonal_d > 0.0 and diagonal_r > 0.0)
    {
      figureOfmeritBB *= offDiagonal/(sqrt(diagonal_d*diagonal_r));
    }
    if(lowerDiagonal > 0.0 and lowerDiagonal_l > 0.0 and lowerDiagonal_u > 0.0)
    {
      double figureOfmerit_lower = lowerDiagonal/(sqrt(lowerDiagonal_u*lowerDiagonal_l));
      figureOfmeritBB *= figureOfmerit_lower;
    }
    //if(offDiagonal == 0.0 or diagonal_d == 0.0 or diagonal_r == 0.0 or lowerDiagonal == 0.0 or lowerDiagonal_l == 0.0 or lowerDiagonal_u == 0.0) std::cout << "Some empty cells found" << std::endl;
    //std::cout << "figureOfmerit 2 = " << figureOfmeritBB << std::endl;    
  }

  //std::cout << "figureOfmerit Bayesian Blocks binning = " << figureOfmeritBB/(h_finalResponseBB->GetNbinsX()) << std::endl;
  std::cout << "figureOfmerit Bayesian Blocks binning = " << figureOfmeritBB << std::endl;

  std::string histfilename=(outfile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  h_response->Write();
  h_responseBB->Write();
  h_finalResponse->Write();
  h_finalResponseBB->Write();
  h_InvariantMass_MuMu->Write();
  h_InvariantMass_MuMuGen->Write();
  tFile->Close();
  std::cout<<"Wrote output file "<<histfilename<<std::endl;
  
  return 0;
}

int computeRelevantArea(std::string infile)
{
  std::string inputfilename=(infile+".root").c_str();
  TFile* file = TFile::Open((infile+".root").c_str());

  TH1F *h_finalResponse = (TH1F*) file->Get("h_finalResponse");
  std::cout << "Total area = " << h_finalResponse->Integral() << std::endl;
  double sumOfDiagonals = 0.0;
  for(int i = 1; i<= h_finalResponse->GetNbinsX(); i++) //without nested loop
  {
    sumOfDiagonals += h_finalResponse->GetBinContent(i, i); 
  }

  double sumOfAllCells = 0.0;
  for(int i = 1; i<= h_finalResponse->GetNbinsX(); i++) 
  {
    for(int j = 1; j<= h_finalResponse->GetNbinsY(); j++)
    { 
      sumOfAllCells += h_finalResponse->GetBinContent(i, j); 
    }
  }
  std::cout << "sumOfAllCells regular binning = " << sumOfAllCells << std::endl;
  std::cout << "sumOfDiagonals regular binning = " << sumOfDiagonals << std::endl;

  TH1F *h_finalResponseBB = (TH1F*) file->Get("h_finalResponseBB");
  std::cout << "Total area = " << h_finalResponseBB->Integral() << std::endl;
  double sumOfDiagonalsBB = 0.0;
  for(int i = 1; i<= h_finalResponseBB->GetNbinsX(); i++) //without nested loop
  {
    sumOfDiagonalsBB += h_finalResponseBB->GetBinContent(i, i);
  }

  double sumOfAllCellsBB = 0.0;
  for(int i = 1; i<= h_finalResponseBB->GetNbinsX(); i++)
  {
    for(int j = 1; j<= h_finalResponseBB->GetNbinsY(); j++)
    {
      sumOfAllCellsBB += h_finalResponseBB->GetBinContent(i, j);
    }
  }
  std::cout << "sumOfAllCells BB binning = " << sumOfAllCellsBB << std::endl;
  std::cout << "sumOfDiagonals BB binning = " << sumOfDiagonalsBB << std::endl;

  return 0;

}

int computeAreaForMichael()
{
  Float_t rebin_array[] = {31.7716571, 38.03552205, 41.20818079, 42.95655956, 45.43517219, 47.36059869, 49.15550893, 52.07252724, 55.66729689, 128.51950869, 132.30855742, 137.39460633, 142.17852438, 146.64154812, 152.19551123, 158.13295279, 164.09264608, 167.87511639, 172.84042825, 176.99907853, 184.54875878, 199.17566438};

  Int_t  nBins = sizeof(rebin_array)/sizeof(Float_t) - 1;

  double areaBB = 0.0;
  for(int i=0; i<nBins; i++)
  {
    double length = rebin_array[i+1]-rebin_array[i];
    //std::cout << "rebin_array[i+1] = " << rebin_array[i+1] << std::endl;
    //std::cout << "rebin_array[i] = " << rebin_array[i] << std::endl;
    std::cout << "length = " << length << std::endl;
    double areaPerCell = 1.0;
    areaPerCell = length*length; 
    std::cout << "areaPerCell = " << areaPerCell << std::endl;
    areaBB += areaPerCell;
  }

  double area10GeV = (double)10*10*17/(double)((150-50)*(200-30));
  double area5GeV = (double)5*5*34/(double)((150-50)*(200-30));
  double area1GeV = (double)170*1*1/(double)((150-50)*(200-30));
  double areaBB_allBins = areaBB/((150-50)*(rebin_array[nBins]-rebin_array[0]));

  std::cout << "area10GeV = " << 1.0/area10GeV << std::endl;
  std::cout << "area5GeV = " << 1.0/area5GeV << std::endl;
  std::cout << "area1GeV = " << 1.0/area1GeV << std::endl;
  std::cout << "areaBB = " << 1.0/areaBB_allBins << std::endl;

  return 0;
}
 
