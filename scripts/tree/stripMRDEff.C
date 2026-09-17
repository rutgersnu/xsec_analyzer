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

//Script that strips MRDEff because James made an oopsie and MRDEff needs further development
void stripMRDEff(){

	int runs = 5000;//5000;
	int subruns = 1;
	//Loop through runs
	for(int rn = 2500; rn < runs; rn++){
		std::cout << "Looping through run " << std::to_string(rn) << std::endl;
        //Open file and trees
        string file_path = "/pnfs/annie/persistent/users/jminock/v1_3_4_world_stv_ntuples/PhaseIITree_0." + std::to_string(rn) + ".0.root";
        string out_path = "/exp/annie/data/users/jminock/temp_add_branches/PhaseIITree_0." + std::to_string(rn) + ".0.root";

	//check if files exist
	if(gSystem->AccessPathName(file_path.c_str())){
		std::cout << "WARNING: " << file_path << " does not exist. Skipping." << std::endl;
		continue;
	}
        //Open file and trees
        TFile *f = new TFile(file_path.c_str(),"read");
//        gSystem->Load("/exp/annie/app/users/jminock/ToolAnalysis/lib/libDataModel.so");
//      gSystem->Load("/exp/annie/app/users/jminock/ToolAnalysis/lib/libDict.so");
//        gInterpreter->GenerateDictionary("map<string,vector<double>>", "map;string;vector");
        TTree *T = (TTree*)f->Get("phaseIITriggerTree");

	T->SetBranchStatus("*", 1);
	//Misc
	//T->SetBranchStatus("MRDEff", 0);
	T->SetBranchStatus("weight_MRDUnc", 0);
	//T->SetBranchStatus("weight_DirtUnc", 0);

        TFile *fnew = new TFile(out_path.c_str(),"recreate");
	auto Tnew = T->CloneTree();

	fnew->Write();
	delete f;
	delete fnew;
	}
}
