#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "THStack.h"
#include "TLegend.h"
#include "TF1.h"
#include "TF2.h"
#include "TLine.h"
#include "TMath.h"
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <cmath>
#include <algorithm>

void checkNEv(){
        //insert variables here
        int rnTrig;
        int rnMRD;
        int evNTrig;
        int evNMRD;
        int trigNum;
        int trigword, hasTank, hasMRD, tankMRDCoinc, noveto;

        //Open file and trees

	int runs = 10;
	int subruns = 1;
	//Loop through runs
	for(int rn = 0; rn < runs; rn++){
//		std::cout << "Looping through run " << std::to_string(rn) << std::endl;
	//Loop through run parts (sub runs)
	for(int srn = 0; srn < subruns; srn++){
//		std::cout << "Looping through subrun " << std::to_string(srn) << std::endl;
        //Open file and trees
        string file_path = "/pnfs/annie/persistent/users/jminock/v1_3_3_world_test/PhaseIITree_0." + std::to_string(rn) + "." + std::to_string(srn) + ".root";

	//check if files exist
	if(gSystem->AccessPathName(file_path.c_str())){
		std::cout << "WARNING: " << file_path << " does not exist. Skipping." << std::endl;
		continue;
	}
 
        TFile *f = new TFile(file_path.c_str(),"read");
        gSystem->Load("/exp/annie/app/users/jminock/ToolAnalysis/lib/libDataModel.so");
//        gInterpreter->GenerateDictionary("map<string,vector<double>>", "map;string;vector");
        TTree *tTrig = (TTree*)f->Get("phaseIITriggerTree");

        Long64_t nentriesTrig = tTrig->GetEntries();
        std::cout << "Number of events for Run " << rn << ": " << nentriesTrig << std::endl;

        delete f;
	} //end of subrun
	} //end of run
}
