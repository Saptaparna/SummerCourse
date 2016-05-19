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

int ReadPPZFiles_Delphes_BB(std::string infile, std::string outfile, std::string outTree){

  std::string inputfilename=(infile+".root").c_str();
  TChain *tree=new TChain("LowPtSUSY_Tree");
  tree->Add(inputfilename.c_str());
  std::cout<<"Opened input file "<<inputfilename<<std::endl;

  TFile *outputFile;
  TTree *outputTree;
  double Mll, Muon1_Pt, Muon2_Pt, Muon1_Eta, Muon2_Eta, Muon1_Phi, Muon2_Phi;

  //output tree for Brian
  std::string outputtreename=(outTree+".root").c_str();
  outputFile = new TFile((outputtreename).c_str(),"RECREATE");
  outputTree=new TTree("MLL_Tree", "MLL_Tree");
  outputTree->Branch("Mll", &Mll, "Mll/D");
  outputTree->Branch("Muon1_Pt", &Muon1_Pt, "Muon1_Pt/D");
  outputTree->Branch("Muon2_Pt", &Muon2_Pt, "Muon2_Pt/D");
  outputTree->Branch("Muon1_Eta", &Muon1_Eta, "Muon1_Eta/D");
  outputTree->Branch("Muon2_Eta", &Muon2_Eta, "Muon2_Eta/D");
  outputTree->Branch("Muon1_Phi", &Muon1_Phi, "Muon1_Phi/D");
  outputTree->Branch("Muon2_Phi", &Muon2_Phi, "Muon2_Phi/D");

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
  /*vector<double>   *GenParticle_PDGId;
  vector<double>   *GenParticle_Pt;
  vector<double>   *GenParticle_Phi;
  vector<double>   *GenParticle_Eta;
  vector<double>   *GenParticle_Mass;
  vector<double>   *GenParticle_Energy;
  Int_t            nGenParticles;
  */
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
  /*GenParticle_PDGId = 0;
  GenParticle_Pt = 0;
  GenParticle_Phi = 0;
  GenParticle_Eta = 0;
  GenParticle_Mass = 0;
  GenParticle_Energy = 0;
*/
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
  /*tree->SetBranchAddress("GenParticle_PDGId", &(GenParticle_PDGId));
  tree->SetBranchAddress("GenParticle_Pt", &(GenParticle_Pt));
  tree->SetBranchAddress("GenParticle_Phi", &(GenParticle_Phi));
  tree->SetBranchAddress("GenParticle_Eta", &(GenParticle_Eta));
  tree->SetBranchAddress("GenParticle_Mass", &(GenParticle_Mass));
  tree->SetBranchAddress("GenParticle_Energy", &(GenParticle_Energy));
*/
  TH1F *h_mu_pt_leading_MuMu=new TH1F("h_mu_pt_leading_MuMu", "Leading muon pT; pT [GeV]; Events/GeV", 1000, 0, 1000); h_mu_pt_leading_MuMu->Sumw2();
  TH1F *h_mu_pt_trailing_MuMu=new TH1F("h_mu_pt_trailing_MuMu", "Trailing muon pT; pT [GeV]; Events/GeV", 1000, 0, 1000); h_mu_pt_trailing_MuMu->Sumw2();
  TH1F *h_mu_eta_leading_MuMu=new TH1F("h_mu_eta_leading_MuMu", "Leading muon #eta ; #eta ; Events", 600, -3.0, 3.0); h_mu_eta_leading_MuMu->Sumw2();
  TH1F *h_mu_eta_trailing_MuMu=new TH1F("h_mu_eta_trailing_MuMu", "Trailing muon #eta; #eta ; Events", 600, -3.0, 3.0); h_mu_eta_trailing_MuMu->Sumw2();
  TH1F *h_mu_phi_leading_MuMu=new TH1F("h_mu_phi_leading_MuMu", "Leading muon #phi ; #phi ; Events", 800, -4.0, 4.0); h_mu_phi_leading_MuMu->Sumw2();
  TH1F *h_mu_phi_trailing_MuMu=new TH1F("h_mu_phi_trailing_MuMu", "Trailing muon #phi; #phi ; Events", 800, -4.0, 4.0); h_mu_phi_trailing_MuMu->Sumw2();
  TH1F *h_InvariantMass_MuMu=new TH1F("h_InvariantMass_MuMu", "Di-muon invariant mass; m_{#mu#mu} [GeV]; Events/GeV", 9000, 0, 300); h_InvariantMass_MuMu->Sumw2();
  TH1F *h_InvariantMass_MuMuGen = new TH1F("h_InvariantMass_MuMuGen", "Z generator particle mass; m_{Z} [GeV]; Events/GeV ", 9000, 0, 300); h_InvariantMass_MuMuGen->Sumw2();

  int nEvents=tree->GetEntries();
  std::cout << "nEvents= " << nEvents << std::endl;

  vector <double> checkDuplicates;
  checkDuplicates.clear();
  int nDup = 0;

  for (int i=0; i<nEvents ; ++i)
    {
     tree->GetEvent(i);

     Mll = 0.0;
     Muon1_Pt = 0.0;
     Muon2_Pt = 0.0;
     Muon1_Eta = 0.0;
     Muon2_Eta = 0.0;
     Muon1_Phi = 0.0;
     Muon2_Phi = 0.0;
 
    // filling the muon's properties into a vector of struct
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
/*
    std::vector<GenParticleInfo> genParticles;
    std::vector<GenParticleInfo> genParticlesLep;
    for (unsigned int j=0; j<GenParticle_PDGId->size(); j++)
    {
       GenParticleInfo genParticle;
       genParticle.pT = GenParticle_Pt->at(j);
       genParticle.eta = GenParticle_Eta->at(j);
       genParticle.phi = GenParticle_Phi->at(j);
       genParticle.mass = GenParticle_Mass->at(j);
       genParticle.energy = GenParticle_Energy->at(j);
       genParticle.pdgID = GenParticle_PDGId->at(j); 
       //cout << "GenParticle_PDGId->at(j) = " << GenParticle_PDGId->at(j) << endl;
       //cout << "genParticle.pdgID = " << genParticle.pdgID << endl;
       if(genParticle.pdgID==23) genParticles.push_back(genParticle);
       if(abs(genParticle.pdgID)==13) genParticlesLep.push_back(genParticle);
    }
    
    //cout << "genParticlesLep.size() = " << genParticlesLep.size() << endl;
    //cout << "genParticles.size() = " << genParticles.size() << endl;
    //if(genParticles.size() > 0) std::cout << "genParticles.at(0).pdgID = " << genParticles.at(0).pdgID << std::endl;
    //if(genParticles.size() > 0 and genParticles.at(0).pdgID==23) std::cout << "genParticles.at(0).mass true Z mass = " << genParticles.at(0).mass << std::endl;
    if(genParticlesLep.size() > 1) {
      TLorentzVector mu1Gen_p4;
      TLorentzVector mu2Gen_p4;
      mu1Gen_p4 = fillTLorentzVector(0.0, 0.0, 0.0, 0.0);
      mu2Gen_p4 = fillTLorentzVector(0.0, 0.0, 0.0, 0.0);
      mu1Gen_p4 = fillTLorentzVector(genParticlesLep.at(0).pT, genParticlesLep.at(0).eta, genParticlesLep.at(0).phi, 0.0);
      mu2Gen_p4 = fillTLorentzVector(genParticlesLep.at(1).pT, genParticlesLep.at(1).eta, genParticlesLep.at(1).phi, 0.0);
      //std::cout << "(mu1Gen_p4 + mu2Gen_p4).M() lepton mass = " << (mu1Gen_p4 + mu2Gen_p4).M() << std::endl;
    }
*/
     // Now sorting this vector of structs
    std::sort (muons.begin(), muons.end(), sortLeptonsInDescendingpT);
    TLorentzVector mu1_p4;
    TLorentzVector mu2_p4;
    mu1_p4 = fillTLorentzVector(0.0, 0.0, 0.0, 0.0);
    mu2_p4 = fillTLorentzVector(0.0, 0.0, 0.0, 0.0);

    if (muons.size() > 0) mu1_p4=fillTLorentzVector(muons.at(0).pT, muons.at(0).eta, muons.at(0).phi, 0.0);
    if (muons.size() > 1) mu2_p4=fillTLorentzVector(muons.at(1).pT, muons.at(1).eta, muons.at(1).phi, 0.0);

    //if(genParticles.size() <=  0) continue;

    if(mu1_p4.Pt()>20.0){
      if(mu2_p4.Pt()>20.0) {
        h_mu_phi_leading_MuMu->Fill(mu1_p4.Phi());
        h_mu_eta_leading_MuMu->Fill(mu1_p4.Eta());
        h_mu_pt_leading_MuMu->Fill(mu1_p4.Pt());
        h_mu_phi_trailing_MuMu->Fill(mu2_p4.Phi());
        h_mu_eta_trailing_MuMu->Fill(mu2_p4.Eta());
        h_mu_pt_trailing_MuMu->Fill(mu2_p4.Pt());
        h_InvariantMass_MuMu->Fill((mu1_p4+mu2_p4).M());
        //h_InvariantMass_MuMuGen->Fill(genParticles.at(0).mass);
        Muon1_Pt = mu1_p4.Pt();
        Muon2_Pt = mu2_p4.Pt();
        Muon1_Eta = mu1_p4.Eta();
        Muon2_Eta = mu2_p4.Eta();
        Muon1_Phi = mu1_p4.Phi();
        Muon2_Phi = mu2_p4.Phi();
        Mll = (mu1_p4+mu2_p4).M();
        if(Mll > 1000)
        {
          std::cout << "Muon1_Pt = " << Muon1_Pt << std::endl;
          std::cout << "Muon2_Pt = " << Muon2_Pt << std::endl;
        }
        double dupCheck = Mll;
        bool bDuplicate = false;
        for(unsigned int uid = 0; uid < checkDuplicates.size(); uid++){
          if (checkDuplicates.at(uid) == dupCheck){
            cout<<dupCheck<<endl;
            bDuplicate = true;
            nDup++;
            //break;
          }
        }
        //if (bDuplicate) continue;
        //else checkDuplicates.push_back(dupCheck);
       }//trailing muon "if"
     }//closing leading muon "if" statement 

    //if(checkDuplicates.size()>0) cout << "checkDuplicates.at(0) = " << checkDuplicates.at(0) << endl;

    outputTree->Fill();

  }//event loop closed

  std::cout << "Number of duplicate records found = " << nDup << std::endl;
  outputTree->Write();
  std::string histfilename=(outfile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  h_mu_pt_leading_MuMu->Write();
  h_mu_pt_trailing_MuMu->Write();
  h_mu_eta_leading_MuMu->Write();
  h_mu_eta_trailing_MuMu->Write();
  h_mu_phi_leading_MuMu->Write();
  h_mu_phi_trailing_MuMu->Write();
  h_InvariantMass_MuMu->Write();
  h_InvariantMass_MuMuGen->Write();
  tFile->Close();
  std::cout<<"Wrote output file "<<histfilename<<std::endl;

  return 0;

}
