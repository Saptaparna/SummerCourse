
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <stdio.h>
#include "TChain.h"
#include "TROOT.h"
#include "TStopwatch.h" 

using namespace std;
        
    void run_Delphes_ZLL_8() {

    gROOT->LoadMacro("FlatTreeMaker_Delphes_Delphes_ZLL_8_C.so");
  
    TChain* fChain = new TChain("Delphes");
    ifstream sourceFiles("Delphes_ZLL_8.txt");
    char line[128];
    int  count = 0;
    cout<< "Adding files from Delphes_ZLL_8 to chain..."<< endl;
     while (sourceFiles >> line) {
        fChain->Add(line);
        ++count;
     }
    cout << count<<" files added!"<<endl;
    sourceFiles.close();
    TStopwatch timer;
    timer.Start();    
    fChain->Process("FlatTreeMaker_Delphes");

    cout << "\n\nDone!" << endl;
    cout << "CPU Time : " << timer.CpuTime() <<endl;
    cout << "RealTime : " << timer.RealTime() <<endl;                             
    cout <<"\n";
}
