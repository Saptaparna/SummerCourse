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

int ReadPPWFiles_Delphes(std::string infile, std::string outfile){

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

  TH1F *h_MET=new TH1F("h_MET", "Missing ET; MET [GeV]; Events/GeV", 600, 0, 600); h_MET->Sumw2();

  TH1F *h_HT = new TH1F("h_HT", "HT (scalar sum of jet pT); H_T [GeV]; Events/GeV", 5000, 0, 5000.0);h_HT->Sumw2();
  TH1F *h_HTb = new TH1F("h_HTb", "HTb (scalar sum of b-jet pT); H_Tb [GeV]; Events/GeV", 5000, 0, 5000.0);h_HTb->Sumw2();
  TH1F *h_jet_pt_leading=new TH1F("h_jet_pt_leading", "Leading jet pT; pT [GeV]; Events/GeV", 10000, 0, 1000); h_jet_pt_leading->Sumw2();
  TH1F *h_jet_pt_trailing=new TH1F("h_jet_pt_trailing", "Trailing jet pT; pT [GeV]; Events/GeV", 10000, 0, 1000); h_jet_pt_trailing->Sumw2();

  TH1F *h_jet_eta_leading=new TH1F("h_jet_eta_leading", "Leading jet #eta; #eta; Events", 600, -3.0, 3.0); h_jet_eta_leading->Sumw2();
  TH1F *h_jet_eta_trailing=new TH1F("h_jet_eta_trailing", "Trailing jet #eta; #eta; Events", 600, -3.0, 3.0); h_jet_eta_trailing->Sumw2();

  TH1F *h_jet_phi_leading=new TH1F("h_jet_phi_leading", "Leading jet #phi; #phi; Events/GeV", 800, -4.0, 4.0); h_jet_phi_leading->Sumw2();
  TH1F *h_jet_phi_trailing=new TH1F("h_jet_phi_trailing", "Trailing jet #phi; #phi; Events/GeV", 800, -4.0, 4.0); h_jet_phi_trailing->Sumw2();

  TH1F *h_jet_mass_leading=new TH1F("h_jet_mass_leading", "Leading jet Mass; Mass [GeV]; Events/GeV", 10000, 0, 1000); h_jet_mass_leading->Sumw2();
  TH1F *h_jet_mass_trailing=new TH1F("h_jet_mass_trailing", "Trailing jet Mass; Mass [GeV]; Events/GeV", 10000, 0, 1000); h_jet_mass_trailing->Sumw2();
  TH1F *h_nJets = new TH1F("h_nJets", "Number of Jets; Number of Jets; Events", 20, -0.5, 19.5);h_nJets->Sumw2();
  TH1F *h_nbJets = new TH1F("h_nbJets", "Number of b-Jets; Number of b-Jets; Events", 20, -0.5, 19.5);h_nbJets->Sumw2();
  
  TH2F *h_PhotonPt_MET = new TH2F("h_PhotonPt_MET", "Photon pT versus MET; Photon pT p[GeV]; MET [GeV]", 500, 0.0, 500.0, 500, 0.0, 500.0);h_PhotonPt_MET->Sumw2(); 

  TH1F *h_TransverseMass_Mu=new TH1F("h_TransverseMass_Mu", "Tranverse Mass with Muons; m_{T} [GeV]; Events/GeV", 9000, 0, 300); h_TransverseMass_Mu->Sumw2();
  TH1F *h_TransverseMass_El=new TH1F("h_TransverseMass_El", "Tranverse Mass with Electrons; m_{T} [GeV]; Events/GeV", 9000, 0, 300); h_TransverseMass_El->Sumw2();

  TH1F *h_mu_pt_leading_Mu=new TH1F("h_mu_pt_leading_Mu", "Leading muon pT; pT [GeV]; Events/GeV", 1000, 0, 1000); h_mu_pt_leading_Mu->Sumw2();
  TH1F *h_el_pt_leading_El=new TH1F("h_el_pt_leading_El", "Leading electron pT; pT [GeV]; Events/GeV", 1000, 0, 1000); h_el_pt_leading_El->Sumw2();
  TH1F *h_mu_eta_leading_Mu=new TH1F("h_mu_eta_leading_Mu", "Leading muon #eta ; #eta ; Events", 600, -3.0, 3.0); h_mu_eta_leading_Mu->Sumw2();
  TH1F *h_el_eta_leading_El=new TH1F("h_el_eta_leading_El", "Leading electron #eta; #eta ; Events", 600, -3.0, 3.0); h_el_eta_leading_El->Sumw2();
  TH1F *h_mu_phi_leading_Mu=new TH1F("h_mu_phi_leading_Mu", "Leading muon #phi ; #phi ; Events", 800, -4.0, 4.0); h_mu_phi_leading_Mu->Sumw2();
  TH1F *h_el_phi_leading_El=new TH1F("h_el_phi_leading_El", "Leading electron #phi; #phi ; Events", 800, -4.0, 4.0); h_el_phi_leading_El->Sumw2();

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
    el1_p4 = fillTLorentzVector(0.0, 0.0, 0.0, 0.0);

    if (electrons.size() > 0) el1_p4=fillTLorentzVector(electrons.at(0).pT, electrons.at(0).eta, electrons.at(0).phi, 0.0);
 
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
   mu1_p4 = fillTLorentzVector(0.0, 0.0, 0.0, 0.0);

   if (muons.size() > 0) mu1_p4=fillTLorentzVector(muons.at(0).pT, muons.at(0).eta, muons.at(0).phi, 0.0);
   
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

  h_MET->Fill(MET);
  TVector2 mu1_transverse;
  TVector2 met_transverse;
  TVector2 el1_transverse;
  
  met_transverse.SetMagPhi(MET, MET_Phi);
  
  if(mu1_p4.Pt() > 20.0) mu1_transverse.SetMagPhi(mu1_p4.Pt(), mu1_p4.Phi());
  if(mu1_p4.Pt() > 20.0){
    h_TransverseMass_Mu->Fill(TMath::Sqrt(2*mu1_p4.Pt()*MET*(1 - TMath::Cos(met_transverse.DeltaPhi(mu1_transverse)))));
    h_mu_pt_leading_Mu->Fill(mu1_p4.Pt());
    h_mu_eta_leading_Mu->Fill(mu1_p4.Eta());
    h_mu_phi_leading_Mu->Fill(mu1_p4.Phi());
  }  

  if(el1_p4.Pt() > 20.0) el1_transverse.SetMagPhi(el1_p4.Pt(), el1_p4.Phi());
  if(el1_p4.Pt() > 20.0){
    h_TransverseMass_El->Fill(TMath::Sqrt(2*el1_p4.Pt()*MET*(1 - TMath::Cos(met_transverse.DeltaPhi(el1_transverse)))));
    h_el_pt_leading_El->Fill(el1_p4.Pt());
    h_el_eta_leading_El->Fill(el1_p4.Eta());
    h_el_phi_leading_El->Fill(el1_p4.Phi());
  }

  if(Jet_vector.size()>0 ) h_jet_pt_leading->Fill(Jet_vector.at(0).Pt());
  if(Jet_vector.size()>1 ) h_jet_pt_trailing->Fill(Jet_vector.at(1).Pt());

  if(Jet_vector.size()>0 ) h_jet_eta_leading->Fill(Jet_vector.at(0).Eta());
  if(Jet_vector.size()>1 ) h_jet_eta_trailing->Fill(Jet_vector.at(1).Eta());

  if(Jet_vector.size()>0 ) h_jet_phi_leading->Fill(Jet_vector.at(0).Phi());
  if(Jet_vector.size()>1 ) h_jet_phi_trailing->Fill(Jet_vector.at(1).Phi());

  if(Jet_vector.size()>0 ) h_jet_mass_leading->Fill(Jet_vector.at(0).M());
  if(Jet_vector.size()>1 ) h_jet_mass_trailing->Fill(Jet_vector.at(1).M());

  h_nJets->Fill(Jet_vector.size());
  h_nbJets->Fill(bJet_vector.size());
  h_HT->Fill(HT);
  h_HTb->Fill(HTb);

  }//event loop closed

  std::string histfilename=(outfile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  h_nJets->Write();
  h_HTb->Write();
  h_nbJets->Write();
  h_jet_pt_leading->Write();
  h_jet_pt_trailing->Write();

  h_jet_eta_leading->Write();
  h_jet_eta_trailing->Write();

  h_jet_phi_leading->Write();
  h_jet_phi_trailing->Write();

  h_jet_mass_leading->Write();
  h_jet_mass_trailing->Write();

  h_HT->Write();
  h_MET->Write();

  h_PhotonPt_MET->Write();
  h_TransverseMass_Mu->Write();
  h_TransverseMass_El->Write();
  h_mu_pt_leading_Mu->Write();
  h_mu_phi_leading_Mu->Write();
  h_mu_eta_leading_Mu->Write();
  h_el_pt_leading_El->Write();
  h_el_phi_leading_El->Write();
  h_el_eta_leading_El->Write();
  tFile->Close();
  std::cout<<"Wrote output file "<<histfilename<<std::endl;

  return 0;

}
