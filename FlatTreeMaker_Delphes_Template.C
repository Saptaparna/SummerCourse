#define FlatTreeMaker_Delphes_cxx
//
// Root > T->Process("FlatTreeMaker_Delphes.C")
// Root > T->Process("FlatTreeMaker_Delphes.C","some options")
// Root > T->Process("FlatTreeMaker_Delphes.C+")
//

#include "FlatTreeMaker_Delphes.h"
#include <TH2.h>
#include <TStyle.h>
using namespace std;


void FlatTreeMaker_Delphes::Begin(TTree * /*tree*/)
{
   TString option = GetOption();
   outputFile = new TFile("LowPtSUSY_Tree_SUFFIX.root","RECREATE");
   outtree=new TTree("LowPtSUSY_Tree", "LowPtSUSY_Tree");
   outtree->Branch("ph_pt", &ph_pt);
   outtree->Branch("ph_phi", &ph_phi);
   outtree->Branch("ph_eta", &ph_eta);
   outtree->Branch("nPhotons", &nPhotons, "nPhotons/I");
   outtree->Branch("el_pt", &el_pt);
   outtree->Branch("el_phi", &el_phi);
   outtree->Branch("el_eta", &el_eta);
   outtree->Branch("el_charge", &el_charge);
   outtree->Branch("nElectrons", &nElectrons, "nElectrons/I");
   outtree->Branch("mu_pt", &mu_pt);
   outtree->Branch("mu_phi", &mu_phi);
   outtree->Branch("mu_eta", &mu_eta);
   outtree->Branch("mu_charge", &mu_charge);
   outtree->Branch("nMuons", &nMuons, "nMuons/I");
   outtree->Branch("jet_pt", &jet_pt);
   outtree->Branch("jet_phi", &jet_phi);
   outtree->Branch("jet_eta", &jet_eta);
   outtree->Branch("jet_mass", &jet_mass);
   outtree->Branch("nJets", &nJets, "nJets/I");
   outtree->Branch("jet_btag", &jet_btag); 
   outtree->Branch("MET", &MET, "MET/F");
   outtree->Branch("MET_Phi", &MET_Phi, "MET_Phi/F");

}

void FlatTreeMaker_Delphes::SlaveBegin(TTree * /*tree*/)
{
   TString option = GetOption();

}

Bool_t FlatTreeMaker_Delphes::Process(Long64_t entry)
{
  GetEntry(entry);

  if(entry % 1000 == 0) std::cout << "Processing event number: " << entry << std::endl;
  ph_pt.clear();
  ph_phi.clear();
  ph_eta.clear();
  nPhotons = -1;
  el_pt.clear();
  el_phi.clear();
  el_eta.clear();
  el_charge.clear();
  nElectrons = -1;
  mu_pt.clear();
  mu_phi.clear();
  mu_eta.clear();
  mu_charge.clear();
  nMuons = -1;
  jet_pt.clear();
  jet_phi.clear();
  jet_eta.clear();
  jet_mass.clear();
  jet_btag.clear();
  nJets = -1;
  MET = 0;
  MET_Phi = -99.0;

  for (Int_t j=0;j<Jet_size;j++){
    
    jet_pt.push_back(Jet_PT[j]);
    jet_eta.push_back(Jet_Eta[j]);
    jet_phi.push_back(Jet_Phi[j]);
    jet_mass.push_back(Jet_Mass[j]);
    jet_btag.push_back(Jet_BTag[j]);
  } 

  nJets = Jet_size;

  for (Int_t j=0;j<Electron_size;j++){

    el_pt.push_back(Electron_PT[j]);
    el_eta.push_back(Electron_Eta[j]);
    el_phi.push_back(Electron_Phi[j]);
    el_charge.push_back(Electron_Charge[j]);
  }

  nElectrons = Electron_size;

  for (Int_t j=0;j<Muon_size;j++){

    mu_pt.push_back(Muon_PT[j]);
    mu_eta.push_back(Muon_Eta[j]);
    mu_phi.push_back(Muon_Phi[j]);
    mu_charge.push_back(Muon_Charge[j]);
  }

  nMuons = Muon_size;

  for (Int_t j=0;j<Photon_size;j++){

    ph_pt.push_back(Photon_PT[j]);
    ph_eta.push_back(Photon_Eta[j]);
    ph_phi.push_back(Photon_Phi[j]);
  }

  nPhotons = Photon_size;


  MET = MissingET_MET[0];
  MET_Phi = MissingET_Phi[0];

  outtree->Fill();
  return kTRUE;
}

void FlatTreeMaker_Delphes::SlaveTerminate()
{

}

void FlatTreeMaker_Delphes::Terminate()
{
  outtree->Write();
}
