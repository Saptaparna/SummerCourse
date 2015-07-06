#!/bin/csh                                                                                                                                                                         
cp FlatTreeMaker_Delphes_Template.C FlatTreeMaker_Delphes.C
#sed -i "s/SUFFIX/$1/g" FlatTreeMaker_Delphes.C

cat > Executable.C << +EOF

void Executable(){

gROOT->ProcessLine(".L FlatTreeMaker_Delphes.C++");

}

+EOF

root -l -b -q Executable.C
rm -f Executable.C
rm -f FlatTreeMaker_Delphes.C
mv FlatTreeMaker_Delphes_C.so FlatTreeMaker_Delphes_${1}_C.so

echo "FINISHED COMPILING"

cat > run_${1}.C <<EOF

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
        
    void run_${1}() {

    gROOT->LoadMacro("FlatTreeMaker_Delphes_${1}_C.so");
  
    TChain* fChain = new TChain("Delphes");
    ifstream sourceFiles("$1.txt");
    char line[128];
    int  count = 0;
    cout<< "Adding files from $1 to chain..."<< endl;
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
EOF
root -l -b -q run_$1.C
echo "WROTE A NEW RUN_{SOURCE}.C FILE"
