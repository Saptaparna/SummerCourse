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

TLorentzVector fillTLorentzVector(double pT, double eta, double phi, double M)
{
  TLorentzVector object_p4;
  object_p4.SetPtEtaPhiM(pT, eta, phi, M);
  return object_p4;
}

typedef struct
{
  float pT;
  float eta;
  float phi;
  int charge;
} LeptonInfo;

typedef struct
{
  float pT;
  float eta;
  float phi;
} PhotonInfo;

typedef struct
{
  float pT;
  float eta;
  float phi;
  float mass;
  int btag;
} JetInfo;

bool sortLeptonsInDescendingpT(LeptonInfo lep1, LeptonInfo lep2)
{
  return (lep1.pT > lep2.pT);
}

bool sortJetsInDescendingpT(JetInfo jet1, JetInfo jet2)
{
  return (jet1.pT > jet2.pT);
}

bool sortJetVectorsInDescendingpT(TLorentzVector jet1, TLorentzVector jet2)
{
  return (jet1.Pt() > jet2.Pt());
}

bool sortPhotonsInDescendingpT(PhotonInfo pho1, PhotonInfo pho2)
{
  return (pho1.pT > pho2.pT);
}

int ReadPPZFiles_Delphes(std::string infile, std::string outfile, std::string channel){

  std::string inputfilename=(infile+".root").c_str();
  TChain *tree=new TChain("LowPtSUSY_Tree");
  tree->Add(inputfilename.c_str());
  std::cout<<"Opened input file "<<inputfilename<<std::endl;

  std::vector<float>   *ph_pt;
  std::vector<float>   *ph_phi;
  vector<float>   *ph_eta;
  Int_t           nPhotons;
  vector<float>   *el_pt;
  vector<float>   *el_phi;
  vector<float>   *el_eta;
  vector<int>     *el_charge;
  Int_t           nElectrons;
  vector<float>   *mu_pt;
  vector<float>   *mu_phi;
  vector<float>   *mu_eta;
  vector<int>     *mu_charge;
  Int_t           nMuons;
  vector<float>   *jet_pt;
  vector<float>   *jet_phi;
  vector<float>   *jet_eta;
  vector<float>   *jet_mass;
  vector<int>     *jet_btag;
  Int_t           nJets;
  Float_t         MET;
  Float_t         MET_Phi; 

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

  TH1F *h_mu_pt_leading_MuMu=new TH1F("h_mu_pt_leading_MuMu", "Leading muon pT; pT [GeV]; Events/GeV", 1000, 0, 1000); h_mu_pt_leading_MuMu->Sumw2();
  TH1F *h_mu_pt_trailing_MuMu=new TH1F("h_mu_pt_trailing_MuMu", "Trailing muon pT; pT [GeV]; Events/GeV", 1000, 0, 1000); h_mu_pt_trailing_MuMu->Sumw2();
  TH1F *h_el_pt_leading_ElEl=new TH1F("h_el_pt_leading_ElEl", "Leading electron pT; pT [GeV]; Events/GeV", 1000, 0, 1000); h_el_pt_leading_ElEl->Sumw2();
  TH1F *h_el_pt_trailing_ElEl=new TH1F("h_el_pt_trailing_ElEl", "Trailing electron pT; pT [GeV]; Events/GeV", 1000, 0, 1000); h_el_pt_trailing_ElEl->Sumw2();
  TH1F *h_mu_eta_leading_MuMu=new TH1F("h_mu_eta_leading_MuMu", "Leading muon #eta ; #eta ; Events", 600, -3.0, 3.0); h_mu_eta_leading_MuMu->Sumw2();
  TH1F *h_mu_eta_trailing_MuMu=new TH1F("h_mu_eta_trailing_MuMu", "Trailing muon #eta; #eta ; Events", 600, -3.0, 3.0); h_mu_eta_trailing_MuMu->Sumw2();
  TH1F *h_el_eta_leading_ElEl=new TH1F("h_el_eta_leading_ElEl", "Leading electron #eta; #eta ; Events", 600, -3.0, 3.0); h_el_eta_leading_ElEl->Sumw2();
  TH1F *h_el_eta_trailing_ElEl=new TH1F("h_el_eta_trailing_ElEl", "Trailing electron #eta; #eta ; Events", 600.0, -3.0, 3.0); h_el_eta_trailing_ElEl->Sumw2();
  TH1F *h_mu_phi_leading_MuMu=new TH1F("h_mu_phi_leading_MuMu", "Leading muon #phi ; #phi ; Events", 800, -4.0, 4.0); h_mu_phi_leading_MuMu->Sumw2();
  TH1F *h_mu_phi_trailing_MuMu=new TH1F("h_mu_phi_trailing_MuMu", "Trailing muon #phi; #phi ; Events", 800, -4.0, 4.0); h_mu_phi_trailing_MuMu->Sumw2();
  TH1F *h_el_phi_leading_ElEl=new TH1F("h_el_phi_leading_ElEl", "Leading electron #phi; #phi ; Events", 800, -4.0, 4.0); h_el_phi_leading_ElEl->Sumw2();
  TH1F *h_el_phi_trailing_ElEl=new TH1F("h_el_phi_trailing_ElEl", "Trailing electron #phi; #phi ; Events", 800.0, -4.0, 4.0); h_el_phi_trailing_ElEl->Sumw2();

  TH1F *h_mu_pt_leading_ElMu=new TH1F("h_mu_pt_leading_ElMu", "Leading muon pT; pT [GeV]; Events/GeV", 1000, 0, 1000); h_mu_pt_leading_ElMu->Sumw2();
  TH1F *h_mu_eta_leading_ElMu=new TH1F("h_mu_eta_leading_ElMu", "Leading muon #eta ; #eta ; Events", 600, -3.0, 3.0); h_mu_eta_leading_ElMu->Sumw2();
  TH1F *h_mu_phi_leading_ElMu=new TH1F("h_mu_phi_leading_ElMu", "Leading muon #phi ; #phi ; Events", 800, -4.0, 4.0); h_mu_phi_leading_ElMu->Sumw2();

  TH1F *h_el_pt_leading_ElMu=new TH1F("h_el_pt_leading_ElMu", "Leading electron pT; pT [GeV]; Events/GeV", 1000, 0, 1000); h_el_pt_leading_ElMu->Sumw2();
  TH1F *h_el_eta_leading_ElMu=new TH1F("h_el_eta_leading_ElMu", "Leading electron #eta ; #eta ; Events", 600, -3.0, 3.0); h_el_eta_leading_ElMu->Sumw2();
  TH1F *h_el_phi_leading_ElMu=new TH1F("h_el_phi_leading_ElMu", "Leading electron #phi ; #phi ; Events", 800, -4.0, 4.0); h_el_phi_leading_ElMu->Sumw2();

  TH1F *h_DeltaPhi_met_mu1_MuMu = new TH1F("h_DeltaPhi_met_mu1_MuMu", "#Delta #phi between MET and the leading muon; #Delta #phi(MET, leading muon); Events/GeV", 3500, -3.5, 3.5); h_DeltaPhi_met_mu1_MuMu->Sumw2();
  TH1F *h_DeltaPhi_met_mu2_MuMu = new TH1F("h_DeltaPhi_met_mu2_MuMu", "#Delta #phi between MET and the trailing muon; #Delta #phi(MET, trailing muon); Events/GeV", 3500, -3.5, 3.5); h_DeltaPhi_met_mu2_MuMu->Sumw2();
  TH1F *h_DeltaPhi_ph_mu1_MuMu = new TH1F("h_DeltaPhi_ph_mu1_MuMu", "#Delta #phi between the photon and the leading muon; #Delta #phi(#gamma, leading muon); Events/GeV", 3500, -3.5, 3.5); h_DeltaPhi_ph_mu1_MuMu->Sumw2();
  TH1F *h_DeltaPhi_ph_mu2_MuMu = new TH1F("h_DeltaPhi_ph_mu2_MuMu", "#Delta #phi between the photon and the trailing muon; #Delta #phi(#gamma, trailing muon); Events/GeV", 3500, -3.5, 3.5); h_DeltaPhi_ph_mu2_MuMu->Sumw2();

  TH1F *h_DeltaPhi_met_mu1_ElMu = new TH1F("h_DeltaPhi_met_mu1_ElMu", "#Delta #phi between MET and the leading muon; #Delta #phi(MET, leading muon); Events/GeV", 3500, -3.5, 3.5); h_DeltaPhi_met_mu1_ElMu->Sumw2();
  TH1F *h_DeltaPhi_met_el1_ElMu = new TH1F("h_DeltaPhi_met_el1_ElMu", "#Delta #phi between MET and the leading electron; #Delta #phi(MET, leading electron); Events/GeV", 3500, -3.5, 3.5); h_DeltaPhi_met_el1_ElMu->Sumw2();
  TH1F *h_DeltaPhi_ph_mu1_ElMu = new TH1F("h_DeltaPhi_ph_mu1_ElMu", "#Delta #phi between the photon and the leading muon; #Delta #phi(#gamma, leading muon); Events/GeV", 3500, -3.5, 3.5); h_DeltaPhi_ph_mu1_ElMu->Sumw2();
  TH1F *h_DeltaPhi_ph_el1_ElMu = new TH1F("h_DeltaPhi_ph_el1_ElMu", "#Delta #phi between the photon and the leading electron; #Delta #phi(#gamma, leading electron); Events/GeV", 3500, -3.5, 3.5); h_DeltaPhi_ph_el1_ElMu->Sumw2();

  TH1F *h_DeltaPhi_met_el1_ElEl = new TH1F("h_DeltaPhi_met_el1_ElEl", "#Delta #phi between MET and the leading electron; #Delta #phi(MET, leading electron); Events/GeV", 3500, -3.5, 3.5); h_DeltaPhi_met_el1_ElEl->Sumw2();
  TH1F *h_DeltaPhi_met_el2_ElEl = new TH1F("h_DeltaPhi_met_el2_ElEl", "#Delta #phi between MET and the trailing electron; #Delta #phi(MET, trailing electron); Events/GeV", 3500, -3.5, 3.5); h_DeltaPhi_met_el2_ElEl->Sumw2();
  TH1F *h_DeltaPhi_ph_el1_ElEl = new TH1F("h_DeltaPhi_ph_el1_ElEl", "#Delta #phi between the photon and the leading electron; #Delta #phi(#gamma, leading electron); Events/GeV", 3500, -3.5, 3.5); h_DeltaPhi_ph_el1_ElEl->Sumw2();
  TH1F *h_DeltaPhi_ph_el2_ElEl = new TH1F("h_DeltaPhi_ph_el2_ElEl", "#Delta #phi between the photon and the trailing electron; #Delta #phi(#gamma, trailing electron); Events/GeV", 3500, -3.5, 3.5); h_DeltaPhi_ph_el2_ElEl->Sumw2();


  TH1F *h_photon_pt =new TH1F("h_photon_pt", "Photon pT; pT [GeV]; Events/GeV", 1000, 0, 1000); h_photon_pt->Sumw2();
  TH1F *h_photon_eta =new TH1F("h_photon_eta", "Photon #eta; #eta ; Events", 600, -3.0, 3.0); h_photon_eta->Sumw2();
  TH1F *h_photon_phi =new TH1F("h_photon_phi", "Photon #phi; #phi ; Events", 800, -4.0, 4.0); h_photon_phi->Sumw2();
  TH1F *h_MET=new TH1F("h_MET", "Missing ET; MET [GeV]; Events/GeV", 600, 0, 600); h_MET->Sumw2();

  TH1F *h_HT = new TH1F("h_HT", "HT (scalar sum of jet pT); H_T [GeV]; Events/GeV", 5000, 0, 5000.0);h_HT->Sumw2();
  TH1F *h_HTb = new TH1F("h_HTb", "HTb (scalar sum of b-jet pT); H_Tb [GeV]; Events/GeV", 5000, 0, 5000.0);h_HTb->Sumw2();
  TH1F *h_jet_pt_leading=new TH1F("h_jet_pt_leading", "Leading jet pT; pT [GeV]; Events/GeV", 10000, 0, 1000); h_jet_pt_leading->Sumw2();
  TH1F *h_jet_pt_trailing=new TH1F("h_jet_pt_trailing", "Trailing jet pT; pT [GeV]; Events/GeV", 10000, 0, 1000); h_jet_pt_trailing->Sumw2();
  TH1F *h_jet_pt_3rd=new TH1F("h_jet_pt_3rd", "3rd jet pT; pT [GeV]; Events/GeV", 10000, 0, 1000); h_jet_pt_3rd->Sumw2();
  TH1F *h_jet_pt_4th=new TH1F("h_jet_pt_4th", "4th jet pT; pT [GeV]; Events/GeV", 10000, 0, 1000); h_jet_pt_4th->Sumw2();
  TH1F *h_jet_pt_5th=new TH1F("h_jet_pt_5th", "5th jet pT; pT [GeV]; Events/GeV", 10000, 0, 1000); h_jet_pt_5th->Sumw2();
  TH1F *h_jet_pt_6th=new TH1F("h_jet_pt_6th", "6th jet pT; pT [GeV]; Events/GeV", 10000, 0, 1000); h_jet_pt_6th->Sumw2();

  TH1F *h_jet_eta_leading=new TH1F("h_jet_eta_leading", "Leading jet #eta; #eta; Events", 600, -3.0, 3.0); h_jet_eta_leading->Sumw2();
  TH1F *h_jet_eta_trailing=new TH1F("h_jet_eta_trailing", "Trailing jet #eta; #eta; Events", 600, -3.0, 3.0); h_jet_eta_trailing->Sumw2();
  TH1F *h_jet_eta_3rd=new TH1F("h_jet_eta_3rd", "3rd jet #eta; #eta; Events/GeV", 600, -3.0, 3.0); h_jet_eta_3rd->Sumw2();
  TH1F *h_jet_eta_4th=new TH1F("h_jet_eta_4th", "4th jet #eta; #eta; Events/GeV", 600, -3.0, 3.0); h_jet_eta_4th->Sumw2();
  TH1F *h_jet_eta_5th=new TH1F("h_jet_eta_5th", "5th jet #eta; #eta; Events/GeV", 600, -3.0, 3.0); h_jet_eta_5th->Sumw2();
  TH1F *h_jet_eta_6th=new TH1F("h_jet_eta_6th", "6th jet #eta; #eta; Events/GeV", 600, -3.0, 3.0); h_jet_eta_6th->Sumw2();

  TH1F *h_jet_phi_leading=new TH1F("h_jet_phi_leading", "Leading jet #phi; #phi; Events/GeV", 800, -4.0, 4.0); h_jet_phi_leading->Sumw2();
  TH1F *h_jet_phi_trailing=new TH1F("h_jet_phi_trailing", "Trailing jet #phi; #phi; Events/GeV", 800, -4.0, 4.0); h_jet_phi_trailing->Sumw2();
  TH1F *h_jet_phi_3rd=new TH1F("h_jet_phi_3rd", "3rd jet #phi; #phi; Events/GeV", 800, -4.0, 4.0); h_jet_phi_3rd->Sumw2();
  TH1F *h_jet_phi_4th=new TH1F("h_jet_phi_4th", "4th jet #phi; #phi; Events/GeV", 800, -4.0, 4.0); h_jet_phi_4th->Sumw2();
  TH1F *h_jet_phi_5th=new TH1F("h_jet_phi_5th", "5th jet #phi; #phi; Events/GeV", 800, -4.0, 4.0); h_jet_phi_5th->Sumw2();
  TH1F *h_jet_phi_6th=new TH1F("h_jet_phi_6th", "6th jet #phi; #phi; Events/GeV", 800, -4.0, 4.0); h_jet_phi_6th->Sumw2();

  TH1F *h_jet_mass_leading=new TH1F("h_jet_mass_leading", "Leading jet Mass; Mass [GeV]; Events/GeV", 10000, 0, 1000); h_jet_mass_leading->Sumw2();
  TH1F *h_jet_mass_trailing=new TH1F("h_jet_mass_trailing", "Trailing jet Mass; Mass [GeV]; Events/GeV", 10000, 0, 1000); h_jet_mass_trailing->Sumw2();
  TH1F *h_jet_mass_3rd=new TH1F("h_jet_mass_3rd", "3rd jet Mass; Mass [GeV]; Events/GeV", 10000, 0, 1000); h_jet_mass_3rd->Sumw2();
  TH1F *h_jet_mass_4th=new TH1F("h_jet_mass_4th", "4th jet Mass; Mass [GeV]; Events/GeV", 10000, 0, 1000); h_jet_mass_4th->Sumw2();
  TH1F *h_jet_mass_5th=new TH1F("h_jet_mass_5th", "5th jet Mass; Mass [GeV]; Events/GeV", 10000, 0, 1000); h_jet_mass_5th->Sumw2();
  TH1F *h_jet_mass_6th=new TH1F("h_jet_mass_6th", "6th jet Mass; Mass [GeV]; Events/GeV", 10000, 0, 1000); h_jet_mass_6th->Sumw2();
  TH1F *h_nJets = new TH1F("h_nJets", "Number of Jets; Number of Jets; Events", 20, -0.5, 19.5);h_nJets->Sumw2();
  TH1F *h_nbJets = new TH1F("h_nbJets", "Number of b-Jets; Number of b-Jets; Events", 20, -0.5, 19.5);h_nbJets->Sumw2();
  TH2F *h_PhotonPt_MET = new TH2F("h_PhotonPt_MET", "Photon pT versus MET; Photon pT p[GeV]; MET [GeV]", 500, 0.0, 500.0, 500, 0.0, 500.0);h_PhotonPt_MET->Sumw2(); 
  TH1F *h_InvariantMass_MuMu=new TH1F("h_InvariantMass_MuMu", "Di-muon invariant mass; m_{#mu#mu} [GeV]; Events/GeV", 9000, 0, 300); h_InvariantMass_MuMu->Sumw2();
  TH1F *h_InvariantMass_MuMuPh=new TH1F("h_InvariantMass_MuMuPh", "Di-muon and photon invariant mass; m_{#mu#mu#gamma} [GeV]; Events/GeV", 9000, 0, 300); h_InvariantMass_MuMuPh->Sumw2();
  TH1F *h_InvariantMass_ElMu=new TH1F("h_InvariantMass_ElMu", "Electron-Muon invariant mass; m_{e#mu} [GeV]; Events/GeV", 9000, 0, 300); h_InvariantMass_ElMu->Sumw2();
  TH1F *h_InvariantMass_ElMuPh=new TH1F("h_InvariantMass_ElMuPh", "Electron-Muon and photon invariant mass; m_{#mu#mu#gamma} [GeV]; Events/GeV", 9000, 0, 300); h_InvariantMass_ElMuPh->Sumw2();

  TH1F *h_InvariantMass_ElEl=new TH1F("h_InvariantMass_ElEl", "Di-electron invariant mass; m_{ee} [GeV]; Events/GeV", 9000, 0, 300); h_InvariantMass_ElEl->Sumw2();
  TH1F *h_InvariantMass_ElElPh=new TH1F("h_InvariantMass_ElElPh", "Di-electron and photon invariant mass; m_{ee} [GeV]; Events/GeV", 9000, 0, 300); h_InvariantMass_ElElPh->Sumw2();
  TH1F *h_InvariantMass_PhPh=new TH1F("h_InvariantMass_PhPh", "Di-photon invariant mass; m_{#gamma#gamma} [GeV]; Events/GeV", 9000, 0, 300); h_InvariantMass_PhPh->Sumw2();
  TH1F *h_InvariantMass_HZZ=new TH1F("h_InvariantMass_HZZ", "Higgs to ZZ invariant mass; m_{#mu#mu#mu#mu} [GeV]; Events/GeV", 9000, 0, 300); h_InvariantMass_HZZ->Sumw2();

  int nEvents=tree->GetEntries();
  std::cout << "nEvents= " << nEvents << std::endl;
  for (int i=0; i<nEvents ; ++i)
    {
     tree->GetEvent(i);

     std::vector<PhotonInfo> photons;
     for (unsigned int j=0; j<ph_pt->size(); ++j)
       {
       PhotonInfo photon;
       photon.pT=ph_pt->at(j);
       photon.eta=ph_eta->at(j);
       photon.phi=ph_phi->at(j);
       photons.push_back(photon);
      }
     // Now sorting this vector of structs
     std::sort (photons.begin(), photons.end(), sortPhotonsInDescendingpT);

     TLorentzVector ph1_p4;
     TLorentzVector ph2_p4;
     ph1_p4 = fillTLorentzVector(0.0, 0.0, 0.0, 0.0);
     ph2_p4 = fillTLorentzVector(0.0, 0.0, 0.0, 0.0);

     //here working with the leading photon.
     if (photons.size() > 1) ph1_p4=fillTLorentzVector(photons.at(0).pT, photons.at(0).eta, photons.at(0).phi, 0.0);
     //also working with the sub-leading photon.
     if (photons.size() > 1) ph2_p4=fillTLorentzVector(photons.at(1).pT, photons.at(1).eta, photons.at(1).phi, 0.0);

     if(ph1_p4.Pt() > 0.0) h_PhotonPt_MET->Fill(ph1_p4.Pt(), MET); 

     std::vector<LeptonInfo> electrons;
     for (unsigned int j=0; j<el_pt->size(); ++j)
     {
        LeptonInfo electron;
        electron.pT=el_pt->at(j);
        electron.eta=el_eta->at(j);
        electron.phi=el_phi->at(j);
        electron.charge=el_charge->at(j);
        electrons.push_back(electron);
      }

     // Now sorting this vector of structs
     std::sort (electrons.begin(), electrons.end(), sortLeptonsInDescendingpT);
     TLorentzVector el1_p4;
     TLorentzVector el2_p4;
     TLorentzVector el3_p4;
     TLorentzVector el4_p4;
     el1_p4 = fillTLorentzVector(0.0, 0.0, 0.0, 0.0);
     el2_p4 = fillTLorentzVector(0.0, 0.0, 0.0, 0.0);
     el3_p4 = fillTLorentzVector(0.0, 0.0, 0.0, 0.0);
     el4_p4 = fillTLorentzVector(0.0, 0.0, 0.0, 0.0);
   
     if (electrons.size() > 0)
     {
       el1_p4=fillTLorentzVector(electrons.at(0).pT, electrons.at(0).eta, electrons.at(0).phi, 0.0);
     }//only execute if an electron exists.

     if (electrons.size() > 1)
     {
       el2_p4=fillTLorentzVector(electrons.at(1).pT, electrons.at(1).eta, electrons.at(1).phi, 0.0);
     }//fill only if the second electron exists.
    
     if (electrons.size() > 2)
     {
       el3_p4=fillTLorentzVector(electrons.at(2).pT, electrons.at(2).eta, electrons.at(2).phi, 0.0);
     }//fill only if the third electron exists.

    if (electrons.size() > 3)
     {
       el4_p4=fillTLorentzVector(electrons.at(3).pT, electrons.at(3).eta, electrons.at(3).phi, 0.0);
     }//fill only if the fourth electron exists.
    

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
     // Now sorting this vector of structs
    std::sort (muons.begin(), muons.end(), sortLeptonsInDescendingpT);
    TLorentzVector mu1_p4;
    TLorentzVector mu2_p4;
    TLorentzVector mu3_p4;
    TLorentzVector mu4_p4;
    mu1_p4 = fillTLorentzVector(0.0, 0.0, 0.0, 0.0);
    mu2_p4 = fillTLorentzVector(0.0, 0.0, 0.0, 0.0);
    mu3_p4 = fillTLorentzVector(0.0, 0.0, 0.0, 0.0);
    mu4_p4 = fillTLorentzVector(0.0, 0.0, 0.0, 0.0);

    if (muons.size() > 0) mu1_p4=fillTLorentzVector(muons.at(0).pT, muons.at(0).eta, muons.at(0).phi, 0.0);
    if (muons.size() > 1) mu2_p4=fillTLorentzVector(muons.at(1).pT, muons.at(1).eta, muons.at(1).phi, 0.0);
    if (muons.size() > 2) mu3_p4=fillTLorentzVector(muons.at(2).pT, muons.at(2).eta, muons.at(2).phi, 0.0);
    if (muons.size() > 3) mu4_p4=fillTLorentzVector(muons.at(3).pT, muons.at(3).eta, muons.at(3).phi, 0.0); 

    if (muons.size() > 3) h_InvariantMass_HZZ->Fill((mu1_p4+mu2_p4+mu3_p4+mu4_p4).M());

    int mumu_event = 0;
    if(mu1_p4.Pt()>20.0){
      if(mu2_p4.Pt()>20.0) {
        mumu_event = 1;
        TVector2 mu1_transverse;
        TVector2 mu2_transverse;
        TVector2 met_transverse;
        TVector2 ph_transverse;
        mu1_transverse.SetMagPhi(mu1_p4.Pt(), mu1_p4.Phi());
        mu1_transverse.SetMagPhi(mu2_p4.Pt(), mu2_p4.Phi());
        ph_transverse.SetMagPhi(ph1_p4.Pt(), ph1_p4.Phi());
        met_transverse.SetMagPhi(MET, MET_Phi);
        h_mu_phi_leading_MuMu->Fill(mu1_p4.Phi());
        h_mu_eta_leading_MuMu->Fill(mu1_p4.Eta());
        h_mu_pt_leading_MuMu->Fill(mu1_p4.Pt());
        h_mu_phi_trailing_MuMu->Fill(mu2_p4.Phi());
        h_mu_eta_trailing_MuMu->Fill(mu2_p4.Eta());
        h_mu_pt_trailing_MuMu->Fill(mu2_p4.Pt());
        h_InvariantMass_MuMu->Fill((mu1_p4+mu2_p4).M());
        if(ph1_p4.Pt() > 0.0) h_InvariantMass_MuMuPh->Fill((mu1_p4+mu2_p4+ph1_p4).M());
        h_DeltaPhi_met_mu1_MuMu->Fill(met_transverse.DeltaPhi(mu1_transverse));
        h_DeltaPhi_met_mu2_MuMu->Fill(met_transverse.DeltaPhi(mu2_transverse));
        h_DeltaPhi_ph_mu1_MuMu->Fill(ph_transverse.DeltaPhi(mu1_transverse));
        h_DeltaPhi_ph_mu2_MuMu->Fill(ph_transverse.DeltaPhi(mu2_transverse));
       }//trailing muon "if"
     }//closing leading muon "if" statement    

   if(photons.size() > 1){
    h_InvariantMass_PhPh->Fill((ph1_p4+ph2_p4).M()); 
   }


   int emu_event = 0;
   if(mu1_p4.Pt()>20.0){
      if(el1_p4.Pt()>20.0 and mumu_event==0) {
        emu_event = 1;
        TVector2 mu1_transverse;
        TVector2 el1_transverse;
        TVector2 met_transverse;
        TVector2 ph_transverse;
        mu1_transverse.SetMagPhi(mu1_p4.Pt(), mu1_p4.Phi());
        el1_transverse.SetMagPhi(el1_p4.Pt(), el1_p4.Phi());
        ph_transverse.SetMagPhi(ph1_p4.Pt(), ph1_p4.Phi());
        met_transverse.SetMagPhi(MET, MET_Phi);
        h_mu_phi_leading_ElMu->Fill(mu1_p4.Phi());
        h_mu_eta_leading_ElMu->Fill(mu1_p4.Eta());
        h_mu_pt_leading_ElMu->Fill(mu1_p4.Pt());
        h_el_phi_leading_ElMu->Fill(el1_p4.Phi());
        h_el_eta_leading_ElMu->Fill(el1_p4.Eta());
        h_el_pt_leading_ElMu->Fill(el1_p4.Pt());
        h_InvariantMass_ElMu->Fill((mu1_p4+el1_p4).M());
        if(ph1_p4.Pt() > 0.0) h_InvariantMass_ElMuPh->Fill((mu1_p4+el1_p4+ph1_p4).M()); 
        h_DeltaPhi_met_mu1_ElMu->Fill(met_transverse.DeltaPhi(mu1_transverse));
        h_DeltaPhi_met_el1_ElMu->Fill(met_transverse.DeltaPhi(el1_transverse));
        h_DeltaPhi_ph_mu1_ElMu->Fill(ph_transverse.DeltaPhi(mu1_transverse));
        h_DeltaPhi_ph_el1_ElMu->Fill(ph_transverse.DeltaPhi(el1_transverse));
     }
   }

   int elel_event = 0;
   if(el1_p4.Pt()>20.0){
      if(el2_p4.Pt()>20.0 and mumu_event==0 and emu_event==0) {
        elel_event = 1;
        TVector2 el1_transverse;
        TVector2 el2_transverse;
        TVector2 met_transverse;
        TVector2 ph_transverse;
        el1_transverse.SetMagPhi(el1_p4.Pt(), el1_p4.Phi());
        el1_transverse.SetMagPhi(el2_p4.Pt(), el2_p4.Phi());
        ph_transverse.SetMagPhi(ph1_p4.Pt(), ph1_p4.Phi());
        met_transverse.SetMagPhi(MET, MET_Phi);
        h_el_phi_leading_ElEl->Fill(el1_p4.Phi());
        h_el_eta_leading_ElEl->Fill(el1_p4.Eta());
        h_el_pt_leading_ElEl->Fill(el1_p4.Pt());
        h_el_phi_trailing_ElEl->Fill(el2_p4.Phi());
        h_el_eta_trailing_ElEl->Fill(el2_p4.Eta());
        h_el_pt_trailing_ElEl->Fill(el2_p4.Pt());
        h_InvariantMass_ElEl->Fill((el1_p4+el2_p4).M());
        if(ph1_p4.Pt() > 0.0) h_InvariantMass_ElElPh->Fill((el1_p4+el2_p4+ph1_p4).M());
        h_DeltaPhi_met_el1_ElEl->Fill(met_transverse.DeltaPhi(el1_transverse));
        h_DeltaPhi_met_el2_ElEl->Fill(met_transverse.DeltaPhi(el2_transverse));
        h_DeltaPhi_ph_el1_ElEl->Fill(ph_transverse.DeltaPhi(el1_transverse));
        h_DeltaPhi_ph_el2_ElEl->Fill(ph_transverse.DeltaPhi(el2_transverse));
       }//trailing elon "if"
     }//closing leading elon "if" statement

   std::vector<JetInfo> jets;
   for (unsigned int j=0; j<jet_pt->size(); ++j)
   {
     JetInfo jet;
     jet.pT = jet_pt->at(j);
     jet.eta = jet_eta->at(j);
     jet.phi = jet_phi->at(j);
     jet.mass = jet_mass->at(j);
     jet.btag = jet_btag->at(j);
     jets.push_back(jet);
   }

   // Now sorting this vector of structs
   std::sort (jets.begin(), jets.end(), sortJetsInDescendingpT);   
   
   double HT = 0.0; //jets are sorted. Don't care as far as HT is concerned.
   double HTb = 0.0; 
   vector<TLorentzVector> Jet_vector;
   Jet_vector.clear();
   vector<TLorentzVector> bJet_vector;
   bJet_vector.clear();
   for(unsigned int k=0; k<jets.size(); ++k)
     {

     TLorentzVector Jet;
     TLorentzVector bJet;
     if(fabs(jets.at(k).eta)<2.4 and jets.at(k).pT>25.0)
       {
       Jet.SetPtEtaPhiM(jets.at(k).pT, jets.at(k).eta, jets.at(k).phi, jets.at(k).mass);

       bool isGoodJet=true;
       for(unsigned int j=0; j<electrons.size(); ++j)
        {
        TLorentzVector Electron;
        Electron.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
        Electron.SetPtEtaPhiM(electrons.at(j).pT, electrons.at(j).eta, electrons.at(j).phi, 0.0);
        double DRjet_el = Jet.DeltaR(Electron);
        if(DRjet_el<0.5) isGoodJet=false;
      }

       for(unsigned int j=0; j<muons.size(); ++j)
         {
         TLorentzVector Muon;
         Muon.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
         Muon.SetPtEtaPhiM(muons.at(j).pT, muons.at(j).eta, muons.at(j).phi, 0.0);
         double DRjet_mu = Jet.DeltaR(Muon);
         if(DRjet_mu<0.5) isGoodJet=false;
      }  
       for(unsigned int l=0; l<photons.size(); ++l)
        {
        TLorentzVector Photon;
        Photon.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
        Photon.SetPtEtaPhiM(photons.at(l).pT, photons.at(l).eta, photons.at(l).phi, 0.0);
        double DRjet_ph = Jet.DeltaR(Photon);
        if(DRjet_ph<0.5) isGoodJet=false;
        }
      if(isGoodJet) Jet_vector.push_back(Jet); 
     }//close four vector if 
       if(fabs(jets.at(k).eta)<2.4 and jets.at(k).pT>25.0 and jets.at(k).btag==1)
       {
       bJet.SetPtEtaPhiM(jets.at(k).pT, jets.at(k).eta, jets.at(k).phi, jets.at(k).mass);

       bool isGoodbJet=true;
       for(unsigned int j=0; j<electrons.size(); ++j)
        {
        TLorentzVector Electron;
        Electron.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
        Electron.SetPtEtaPhiM(electrons.at(j).pT, electrons.at(j).eta, electrons.at(j).phi, 0.0);
        double DRjet_el = bJet.DeltaR(Electron);
        if(DRjet_el<0.5) isGoodbJet=false;
      }

       for(unsigned int j=0; j<muons.size(); ++j)
         {
         TLorentzVector Muon;
         Muon.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
         Muon.SetPtEtaPhiM(muons.at(j).pT, muons.at(j).eta, muons.at(j).phi, 0.0);
         double DRjet_mu = bJet.DeltaR(Muon);
         if(DRjet_mu<0.5) isGoodbJet=false;
      } 
       for(unsigned int l=0; l<photons.size(); ++l)
        {
        TLorentzVector Photon;
        Photon.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
        Photon.SetPtEtaPhiM(photons.at(l).pT, photons.at(l).eta, photons.at(l).phi, 0.0);
        double DRjet_ph = bJet.DeltaR(Photon);
        if(DRjet_ph<0.5) isGoodbJet=false;
        }

       if(isGoodbJet) Jet_vector.push_back(bJet);
     }//close four vector if
  }//close jet loop

  // Now sorting this vector of structs
  std::sort (Jet_vector.begin(), Jet_vector.end(), sortJetVectorsInDescendingpT);

  for(unsigned int m=0; m<Jet_vector.size(); m++)
    {
    HT += Jet_vector.at(m).Pt();
    }

  for(unsigned int m=0; m<bJet_vector.size(); m++)
    {
     HTb += bJet_vector.at(m).Pt();
    }

  if(channel=="MuMu" and mumu_event==1 and emu_event==0){

   if(ph1_p4.Pt() > 0.0){
       h_photon_pt->Fill(ph1_p4.Pt());
       h_photon_eta->Fill(ph1_p4.Eta());
       h_photon_phi->Fill(ph1_p4.Phi());
     }

   h_MET->Fill(MET);
 
  if(Jet_vector.size()>0 ) h_jet_pt_leading->Fill(Jet_vector.at(0).Pt());
  if(Jet_vector.size()>1 ) h_jet_pt_trailing->Fill(Jet_vector.at(1).Pt());
  if(Jet_vector.size()>2 ) h_jet_pt_3rd->Fill(Jet_vector.at(2).Pt());
  if(Jet_vector.size()>3 ) h_jet_pt_4th->Fill(Jet_vector.at(3).Pt());
  if(Jet_vector.size()>4 ) h_jet_pt_5th->Fill(Jet_vector.at(4).Pt());
  if(Jet_vector.size()>5 ) h_jet_pt_6th->Fill(Jet_vector.at(5).Pt());

  if(Jet_vector.size()>0 ) h_jet_eta_leading->Fill(Jet_vector.at(0).Eta());
  if(Jet_vector.size()>1 ) h_jet_eta_trailing->Fill(Jet_vector.at(1).Eta());
  if(Jet_vector.size()>2 ) h_jet_eta_3rd->Fill(Jet_vector.at(2).Eta());
  if(Jet_vector.size()>3 ) h_jet_eta_4th->Fill(Jet_vector.at(3).Eta());
  if(Jet_vector.size()>4 ) h_jet_eta_5th->Fill(Jet_vector.at(4).Eta());
  if(Jet_vector.size()>5 ) h_jet_eta_6th->Fill(Jet_vector.at(5).Eta());

  if(Jet_vector.size()>0 ) h_jet_phi_leading->Fill(Jet_vector.at(0).Phi());
  if(Jet_vector.size()>1 ) h_jet_phi_trailing->Fill(Jet_vector.at(1).Phi());
  if(Jet_vector.size()>2 ) h_jet_phi_3rd->Fill(Jet_vector.at(2).Phi());
  if(Jet_vector.size()>3 ) h_jet_phi_4th->Fill(Jet_vector.at(3).Phi());
  if(Jet_vector.size()>4 ) h_jet_phi_5th->Fill(Jet_vector.at(4).Phi());
  if(Jet_vector.size()>5 ) h_jet_phi_6th->Fill(Jet_vector.at(5).Phi());

  if(Jet_vector.size()>0 ) h_jet_mass_leading->Fill(Jet_vector.at(0).M());
  if(Jet_vector.size()>1 ) h_jet_mass_trailing->Fill(Jet_vector.at(1).M());
  if(Jet_vector.size()>2 ) h_jet_mass_3rd->Fill(Jet_vector.at(2).M());
  if(Jet_vector.size()>3 ) h_jet_mass_4th->Fill(Jet_vector.at(3).M());
  if(Jet_vector.size()>4 ) h_jet_mass_5th->Fill(Jet_vector.at(4).M());
  if(Jet_vector.size()>5 ) h_jet_mass_6th->Fill(Jet_vector.at(5).M());

  h_nJets->Fill(Jet_vector.size());
  h_nbJets->Fill(bJet_vector.size());
  h_HT->Fill(HT);
  h_HTb->Fill(HTb);
  }

  if(channel=="ElMu" and mumu_event==0 and emu_event==1){

  if(ph1_p4.Pt() > 0.0){
       h_photon_pt->Fill(ph1_p4.Pt());
       h_photon_eta->Fill(ph1_p4.Eta());
       h_photon_phi->Fill(ph1_p4.Phi());
     }

  h_MET->Fill(MET);
 
  if(Jet_vector.size()>0 ) h_jet_pt_leading->Fill(Jet_vector.at(0).Pt());
  if(Jet_vector.size()>1 ) h_jet_pt_trailing->Fill(Jet_vector.at(1).Pt());
  if(Jet_vector.size()>2 ) h_jet_pt_3rd->Fill(Jet_vector.at(2).Pt());
  if(Jet_vector.size()>3 ) h_jet_pt_4th->Fill(Jet_vector.at(3).Pt());
  if(Jet_vector.size()>4 ) h_jet_pt_5th->Fill(Jet_vector.at(4).Pt());
  if(Jet_vector.size()>5 ) h_jet_pt_6th->Fill(Jet_vector.at(5).Pt());

  if(Jet_vector.size()>0 ) h_jet_eta_leading->Fill(Jet_vector.at(0).Eta());
  if(Jet_vector.size()>1 ) h_jet_eta_trailing->Fill(Jet_vector.at(1).Eta());
  if(Jet_vector.size()>2 ) h_jet_eta_3rd->Fill(Jet_vector.at(2).Eta());
  if(Jet_vector.size()>3 ) h_jet_eta_4th->Fill(Jet_vector.at(3).Eta());
  if(Jet_vector.size()>4 ) h_jet_eta_5th->Fill(Jet_vector.at(4).Eta());
  if(Jet_vector.size()>5 ) h_jet_eta_6th->Fill(Jet_vector.at(5).Eta());

  if(Jet_vector.size()>0 ) h_jet_phi_leading->Fill(Jet_vector.at(0).Phi());
  if(Jet_vector.size()>1 ) h_jet_phi_trailing->Fill(Jet_vector.at(1).Phi());
  if(Jet_vector.size()>2 ) h_jet_phi_3rd->Fill(Jet_vector.at(2).Phi());
  if(Jet_vector.size()>3 ) h_jet_phi_4th->Fill(Jet_vector.at(3).Phi());
  if(Jet_vector.size()>4 ) h_jet_phi_5th->Fill(Jet_vector.at(4).Phi());
  if(Jet_vector.size()>5 ) h_jet_phi_6th->Fill(Jet_vector.at(5).Phi());

  if(Jet_vector.size()>0 ) h_jet_mass_leading->Fill(Jet_vector.at(0).M());
  if(Jet_vector.size()>1 ) h_jet_mass_trailing->Fill(Jet_vector.at(1).M());
  if(Jet_vector.size()>2 ) h_jet_mass_3rd->Fill(Jet_vector.at(2).M());
  if(Jet_vector.size()>3 ) h_jet_mass_4th->Fill(Jet_vector.at(3).M());
  if(Jet_vector.size()>4 ) h_jet_mass_5th->Fill(Jet_vector.at(4).M());
  if(Jet_vector.size()>5 ) h_jet_mass_6th->Fill(Jet_vector.at(5).M());

  h_nJets->Fill(Jet_vector.size());
  h_nbJets->Fill(bJet_vector.size());
  h_HT->Fill(HT);
  h_HTb->Fill(HTb);
  }

 if(channel=="ElEl" and mumu_event==0 and emu_event==0 and elel_event==1){

  if(ph1_p4.Pt() > 0.0){
       h_photon_pt->Fill(ph1_p4.Pt());
       h_photon_eta->Fill(ph1_p4.Eta());
       h_photon_phi->Fill(ph1_p4.Phi());
     }

  h_MET->Fill(MET);

  if(Jet_vector.size()>0 ) h_jet_pt_leading->Fill(Jet_vector.at(0).Pt());
  if(Jet_vector.size()>1 ) h_jet_pt_trailing->Fill(Jet_vector.at(1).Pt());
  if(Jet_vector.size()>2 ) h_jet_pt_3rd->Fill(Jet_vector.at(2).Pt());
  if(Jet_vector.size()>3 ) h_jet_pt_4th->Fill(Jet_vector.at(3).Pt());
  if(Jet_vector.size()>4 ) h_jet_pt_5th->Fill(Jet_vector.at(4).Pt());
  if(Jet_vector.size()>5 ) h_jet_pt_6th->Fill(Jet_vector.at(5).Pt());

  if(Jet_vector.size()>0 ) h_jet_eta_leading->Fill(Jet_vector.at(0).Eta());
  if(Jet_vector.size()>1 ) h_jet_eta_trailing->Fill(Jet_vector.at(1).Eta());
  if(Jet_vector.size()>2 ) h_jet_eta_3rd->Fill(Jet_vector.at(2).Eta());
  if(Jet_vector.size()>3 ) h_jet_eta_4th->Fill(Jet_vector.at(3).Eta());
  if(Jet_vector.size()>4 ) h_jet_eta_5th->Fill(Jet_vector.at(4).Eta());
  if(Jet_vector.size()>5 ) h_jet_eta_6th->Fill(Jet_vector.at(5).Eta());

  if(Jet_vector.size()>0 ) h_jet_phi_leading->Fill(Jet_vector.at(0).Phi());
  if(Jet_vector.size()>1 ) h_jet_phi_trailing->Fill(Jet_vector.at(1).Phi());
  if(Jet_vector.size()>2 ) h_jet_phi_3rd->Fill(Jet_vector.at(2).Phi());
  if(Jet_vector.size()>3 ) h_jet_phi_4th->Fill(Jet_vector.at(3).Phi());
  if(Jet_vector.size()>4 ) h_jet_phi_5th->Fill(Jet_vector.at(4).Phi());
  if(Jet_vector.size()>5 ) h_jet_phi_6th->Fill(Jet_vector.at(5).Phi());

  if(Jet_vector.size()>0 ) h_jet_mass_leading->Fill(Jet_vector.at(0).M());
  if(Jet_vector.size()>1 ) h_jet_mass_trailing->Fill(Jet_vector.at(1).M());
  if(Jet_vector.size()>2 ) h_jet_mass_3rd->Fill(Jet_vector.at(2).M());
  if(Jet_vector.size()>3 ) h_jet_mass_4th->Fill(Jet_vector.at(3).M());
  if(Jet_vector.size()>4 ) h_jet_mass_5th->Fill(Jet_vector.at(4).M());
  if(Jet_vector.size()>5 ) h_jet_mass_6th->Fill(Jet_vector.at(5).M());

  h_nJets->Fill(Jet_vector.size());
  h_HT->Fill(HT);
  }

  }//event loop closed

  std::string histfilename=(outfile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  h_nJets->Write();
  h_HTb->Write();
  h_nbJets->Write();
  h_jet_pt_leading->Write();
  h_jet_pt_trailing->Write();
  h_jet_pt_3rd->Write();
  h_jet_pt_4th->Write();
  h_jet_pt_5th->Write();
  h_jet_pt_6th->Write();

  h_jet_eta_leading->Write();
  h_jet_eta_trailing->Write();
  h_jet_eta_3rd->Write();
  h_jet_eta_4th->Write();
  h_jet_eta_5th->Write();
  h_jet_eta_6th->Write();

  h_jet_phi_leading->Write();
  h_jet_phi_trailing->Write();
  h_jet_phi_3rd->Write();
  h_jet_phi_4th->Write();
  h_jet_phi_5th->Write();
  h_jet_phi_6th->Write();

  h_jet_mass_leading->Write();
  h_jet_mass_trailing->Write();
  h_jet_mass_3rd->Write();
  h_jet_mass_4th->Write();
  h_jet_mass_5th->Write();
  h_jet_mass_6th->Write(); 

  h_HT->Write();
  h_MET->Write();

  h_mu_pt_leading_ElMu->Write();
  h_mu_eta_leading_ElMu->Write();
  h_mu_phi_leading_ElMu->Write();

  h_el_pt_leading_ElMu->Write();
  h_el_eta_leading_ElMu->Write();
  h_el_phi_leading_ElMu->Write();

  h_mu_pt_leading_MuMu->Write();
  h_mu_pt_trailing_MuMu->Write();
  h_mu_eta_leading_MuMu->Write();
  h_mu_eta_trailing_MuMu->Write();
  h_mu_phi_leading_MuMu->Write();
  h_mu_phi_trailing_MuMu->Write();

  h_DeltaPhi_met_mu1_MuMu->Write();
  h_DeltaPhi_met_mu2_MuMu->Write();
  h_DeltaPhi_ph_mu1_MuMu->Write();
  h_DeltaPhi_ph_mu2_MuMu->Write();

  h_DeltaPhi_met_mu1_ElMu->Write();
  h_DeltaPhi_met_el1_ElMu->Write();
  h_DeltaPhi_ph_mu1_ElMu->Write();
  h_DeltaPhi_ph_el1_ElMu->Write();

  h_photon_pt->Write();
  h_photon_eta->Write();
  h_photon_phi->Write();

  h_el_pt_leading_ElEl->Write();
  h_el_pt_trailing_ElEl->Write();
  h_el_eta_leading_ElEl->Write();
  h_el_eta_trailing_ElEl->Write();
  h_el_phi_leading_ElEl->Write();
  h_el_phi_trailing_ElEl->Write();

  h_DeltaPhi_met_el1_ElEl->Write();
  h_DeltaPhi_met_el2_ElEl->Write();
  h_DeltaPhi_ph_el1_ElEl->Write();
  h_DeltaPhi_ph_el2_ElEl->Write();


  h_PhotonPt_MET->Write();
  h_InvariantMass_MuMu->Write();
  h_InvariantMass_MuMuPh->Write();
  h_InvariantMass_ElMu->Write(); 
  h_InvariantMass_ElMuPh->Write();
  h_InvariantMass_ElEl->Write();
  h_InvariantMass_ElElPh->Write();
  h_InvariantMass_PhPh->Write();
  h_InvariantMass_HZZ->Write();
  tFile->Close();
  std::cout<<"Wrote output file "<<histfilename<<std::endl;

  return 0;

}
