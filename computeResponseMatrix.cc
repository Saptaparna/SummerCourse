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

TLorentzVector fillTLorentzVector(double pT, double eta, double phi, double M)
{
  TLorentzVector object_p4;
  object_p4.SetPtEtaPhiM(pT, eta, phi, M);
  return object_p4;
}

typedef struct
{
  double pT;
  double eta;
  double phi;
  double mass;
  double energy;
  int pdgID;
} GenParticleInfo;

typedef struct
{
  double pT;
  double eta;
  double phi;
  int charge;
} LeptonInfo;


bool sortLeptonsInDescendingpT(LeptonInfo lep1, LeptonInfo lep2)
{
  return (lep1.pT > lep2.pT);
}

int computeResponseMatrix(std::string infile, std::string outfile){

  std::string inputfilename=(infile+".root").c_str();
  TChain *tree=new TChain("LowPtSUSY_Tree");
  tree->Add(inputfilename.c_str());
  std::cout<<"Opened input file "<<inputfilename<<std::endl;

  std::vector<double>   *ph_pt;
  std::vector<double>   *ph_phi;
  vector<double>   *ph_eta;
  Int_t           nPhotons;
  vector<double>   *el_pt;
  vector<double>   *el_phi;
  vector<double>   *el_eta;
  vector<int>     *el_charge;
  Int_t           nElectrons;
  vector<double>   *mu_pt;
  vector<double>   *mu_phi;
  vector<double>   *mu_eta;
  vector<int>     *mu_charge;
  Int_t           nMuons;
  vector<double>   *jet_pt;
  vector<double>   *jet_phi;
  vector<double>   *jet_eta;
  vector<double>   *jet_mass;
  vector<int>     *jet_btag;
  Int_t           nJets;
  Float_t         MET;
  Float_t         MET_Phi;
  vector<double>   *GenParticle_PDGId;
  vector<double>   *GenParticle_Pt;
  vector<double>   *GenParticle_Phi;
  vector<double>   *GenParticle_Eta;
  vector<double>   *GenParticle_Mass;
  vector<double>   *GenParticle_Energy;
  Int_t            nGenParticles;

  ph_pt = 0;
  ph_phi = 0;
  ph_eta = 0;
  el_pt = 0;
  el_phi = 0;
  el_eta = 0;
  el_charge = 0;
  mu_pt = 0;
  mu_phi = 0;
  mu_eta = 0;
  mu_charge = 0;
  jet_pt = 0;
  jet_phi = 0;
  jet_eta = 0;
  jet_mass = 0;
  jet_btag = 0;
  GenParticle_PDGId = 0;
  GenParticle_Pt = 0;
  GenParticle_Phi = 0;
  GenParticle_Eta = 0;
  GenParticle_Mass = 0;
  GenParticle_Energy = 0;

  tree->SetBranchAddress("ph_pt", &(ph_pt));
  tree->SetBranchAddress("ph_phi", &(ph_phi));
  tree->SetBranchAddress("ph_eta", &(ph_eta));
  tree->SetBranchAddress("nPhotons", &(nPhotons));
  tree->SetBranchAddress("el_pt", &(el_pt));
  tree->SetBranchAddress("el_eta", &(el_eta));
  tree->SetBranchAddress("el_phi", &(el_phi));
  tree->SetBranchAddress("el_charge", &(el_charge));
  tree->SetBranchAddress("nElectrons", &(nElectrons));
  tree->SetBranchAddress("mu_pt", &(mu_pt));
  tree->SetBranchAddress("mu_eta", &(mu_eta));
  tree->SetBranchAddress("mu_phi", &(mu_phi));
  tree->SetBranchAddress("mu_charge", &(mu_charge));
  tree->SetBranchAddress("nMuons", &(nMuons));
  tree->SetBranchAddress("jet_pt", &(jet_pt));
  tree->SetBranchAddress("jet_phi", &(jet_phi));
  tree->SetBranchAddress("jet_eta", &(jet_eta));
  tree->SetBranchAddress("jet_mass", &(jet_mass));
  tree->SetBranchAddress("jet_btag", &(jet_btag));
  tree->SetBranchAddress("nJets", &(nJets));
  tree->SetBranchAddress("MET", &(MET));
  tree->SetBranchAddress("MET_Phi", &(MET_Phi));
  tree->SetBranchAddress("GenParticle_PDGId", &(GenParticle_PDGId));
  tree->SetBranchAddress("GenParticle_Pt", &(GenParticle_Pt));
  tree->SetBranchAddress("GenParticle_Phi", &(GenParticle_Phi));
  tree->SetBranchAddress("GenParticle_Eta", &(GenParticle_Eta));
  tree->SetBranchAddress("GenParticle_Mass", &(GenParticle_Mass));
  tree->SetBranchAddress("GenParticle_Energy", &(GenParticle_Energy));

  //TH2F *h_response=new TH2F("h_response", "Response Matrix; Reconstructed Mass [GeV]; Generated Mass [GeV]", 20,  50.0, 150.0, 20, 50.0, 150.0); h_response->Sumw2();
  TH2F *h_response=new TH2F("h_response", "Response Matrix; Reconstructed Mass [GeV]; Generated Mass [GeV]", 10, 50.0, 150.0, 10, 50.0, 150.0); h_response->Sumw2();
  TH1F *h_InvariantMass_MuMu=new TH1F("h_InvariantMass_MuMu", "Di-muon invariant mass; m_{#mu#mu} [GeV]; Events/GeV", 9000, 0, 300); h_InvariantMass_MuMu->Sumw2();
  TH1F *h_InvariantMass_MuMuGen = new TH1F("h_InvariantMass_MuMuGen", "Z generator particle mass; m_{Z} [GeV]; Events/GeV ", 9000, 0, 300); h_InvariantMass_MuMuGen->Sumw2();

  Float_t rebin_array[] = {50.0437002 ,   53.88551609,   60.76757885,   66.05113402, 69.82073386,   72.30532915,   75.81628342,   78.63675151, 80.47660605,   81.69053905,   82.32946915,   83.36438819, 84.26875909,   84.99499012,   85.51975447,   86.41110934, 86.99699361,   87.63302003,   88.29013696,   88.79892309, 89.5140167 ,   90.13442635,   92.01304163,   92.62524283, 93.04615559,   93.5202137 ,   94.02254874,   94.75193932, 95.29608935,   96.00046705,   96.56100526,   97.39413587, 98.58067156,  100.16170276,  102.20873394,  105.07039449, 107.40429574,  112.31184178,  121.2490133 ,  134.39090819, 149.84399272};

  Int_t  nBins = sizeof(rebin_array)/sizeof(Float_t) - 1;

  TH2F *h_responseBB = new TH2F("h_responseBB", "h_responseBB", nBins, rebin_array, nBins, rebin_array); h_responseBB->Sumw2();

  int nEvents=tree->GetEntries();
  //int nEvents=200000;
  std::cout << "nEvents= " << nEvents << std::endl;
  TMatrixD responseMatrix(h_response->GetNbinsX(), h_response->GetNbinsY());

  for (int i=0; i<nEvents ; ++i)
  {
     tree->GetEvent(i);
     std::vector<LeptonInfo> muons;
     for (unsigned int j=0; j<mu_pt->size(); ++j)
     {
       LeptonInfo muon;
       muon.pT=mu_pt->at(j);
       muon.eta=mu_eta->at(j);
       muon.phi=mu_phi->at(j);
       muon.charge=mu_charge->at(j);
       muons.push_back(muon);
     }

     std::vector<GenParticleInfo> genParticles;
     for (unsigned int j=0; j<GenParticle_PDGId->size(); j++)
     {
       GenParticleInfo genParticle;
       genParticle.pT = GenParticle_Pt->at(j);
       genParticle.eta = GenParticle_Eta->at(j);
       genParticle.phi = GenParticle_Phi->at(j);
       genParticle.mass = GenParticle_Mass->at(j);
       genParticle.energy = GenParticle_Energy->at(j);
       genParticle.pdgID = GenParticle_PDGId->at(j);
       if(genParticle.pdgID==23) genParticles.push_back(genParticle);
     }

    std::sort (muons.begin(), muons.end(), sortLeptonsInDescendingpT);
    TLorentzVector mu1_p4;
    TLorentzVector mu2_p4;
    mu1_p4 = fillTLorentzVector(0.0, 0.0, 0.0, 0.0);
    mu2_p4 = fillTLorentzVector(0.0, 0.0, 0.0, 0.0);

    if (muons.size() > 0) mu1_p4=fillTLorentzVector(muons.at(0).pT, muons.at(0).eta, muons.at(0).phi, 0.0);
    if (muons.size() > 1) mu2_p4=fillTLorentzVector(muons.at(1).pT, muons.at(1).eta, muons.at(1).phi, 0.0);
  
   

    if(genParticles.size() <=  0) continue;
    if(mu1_p4.Pt()>20.0)
    {
      if(mu2_p4.Pt()>20.0) 
      {
        double genMass = genParticles.at(0).mass;
        double recoMass = (mu1_p4+mu2_p4).M();
        h_InvariantMass_MuMu->Fill(recoMass);
        h_InvariantMass_MuMuGen->Fill(genMass);
        h_response->Fill(recoMass, genMass);
        h_responseBB->Fill(recoMass, genMass);
      }
    }

  }//event loop closed

  TH2F *h_finalResponse=new TH2F("h_finalResponse", "Response Matrix; Reconstructed Mass [GeV]; Generator Mass [GeV]", 10, 50.0, 150.0, 10, 50.0, 150.0); h_finalResponse->Sumw2();
  //TH2F *h_finalResponse=new TH2F("h_finalResponse", "Response Matrix; Reconstructed Mass [GeV]; Generator Mass [GeV]", 10, 50.0, 150.0, 10, 50.0, 150.0); h_finalResponse->Sumw2(); 
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
      //std::cout << "h_finalResponse->GetBinContent(k, l) = " << h_finalResponse->GetBinContent(k, l) << std::endl;
    }
  }

  //compute figure of merit
  double figureOfmerit = 0.0;
  for(int k = 1; k<= h_finalResponse->GetNbinsX(); k++) //reco information
  {
    for(int l = 1; l<= h_finalResponse->GetNbinsY(); l++)  
    {
      if(k==l){
        if(h_finalResponse->GetBinContent(l, l) > 0 and h_finalResponse->GetBinContent(k, k) > 0 and h_finalResponse->GetBinContent(k+1, l) > 0) std::cout << "figureOfmerit first diagonal = " << h_finalResponse->GetBinContent(k+1, l)/(sqrt(h_finalResponse->GetBinContent(l, l)*h_finalResponse->GetBinContent(k, k)))  << " (" << k << " , " << l << ") "                           
      <<  std::endl;
        if(h_finalResponse->GetBinContent(k, l) > 0 and h_finalResponse->GetBinContent(k-1, l) > 0 and h_finalResponse->GetBinContent(k, l+1) > 0) std::cout << "figureOfmerit second diagonal = " << h_finalResponse->GetBinContent(k, l)/(sqrt(h_finalResponse->GetBinContent(k-1, l)*h_finalResponse->GetBinContent(k, l+1)))  << " (" << k << " , " << l << ") "                           
      <<  std::endl;

      }
    }
  }

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
