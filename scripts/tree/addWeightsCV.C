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

void addWeightsCV(){
        //weight branches
	vector<double>* fTunedCentralValue = new vector<double>();

	//sample branches
	vector<double>* sfTunedCentralValue = new vector<double>();

	int runs = 1;
	int subruns = 1;
	//Loop through runs
	for(int rn = 0; rn < runs; rn++){
		std::cout << "Looping through run " << std::to_string(rn) << std::endl;
	//Loop through run parts (sub runs)
	for(int srn = 0; srn < subruns; srn++){
//		std::cout << "Looping through subrun " << std::to_string(srn) << std::endl;
        //Open file and trees
        string weight_path = "/pnfs/annie/persistent/users/jminock/v1_3_4_world_stv_ntuples/PhaseIITree_numuMC_MRDY_stv_ntuple.root";
	string sample_path = "/exp/annie/data/users/jminock/temp_add_branches/PhaseIITree_jupyter_stv_ntuple.root";


	//check if files exist
	if(gSystem->AccessPathName(weight_path.c_str())){
		std::cout << "WARNING: " << weight_path << " does not exist. Skipping." << std::endl;
		continue;
	}
	if(gSystem->AccessPathName(sample_path.c_str())){
		std::cout << "WARNING: " << sample_path << " does not exist. Skipping." << std::endl;
		continue;
	}
        TFile *fweight = new TFile(weight_path.c_str(),"read");
        TFile *fsample = new TFile(sample_path.c_str(),"update");
//        gSystem->Load("/exp/annie/app/users/jminock/ToolAnalysis/lib/libDataModel.so");
//      gSystem->Load("/exp/annie/app/users/jminock/ToolAnalysis/lib/libDict.so");
//        gInterpreter->GenerateDictionary("map<string,vector<double>>", "map;string;vector");
        TTree *tWeight = (TTree*)fweight->Get("phaseIITriggerTree");
        TTree *tSample = (TTree*)fsample->Get("phaseIITriggerTree");

	//Weights
	tWeight->SetBranchAddress("weight_TunedCentralValue_UBGenie",&fTunedCentralValue);

	TBranch *WXTCV = tSample->Branch("weight_TunedCentralValue_UBGenie",&sfTunedCentralValue);

        double muon_m = 105.66;
        Long64_t nentriesTrig = tWeight->GetEntries();
        Long64_t nentriesTrig2 = tSample->GetEntries();
        //fill histograms
	for (Long64_t i = 0; i < nentriesTrig2; i++) {
//	for (Long64_t i = 0; i < 10; i++) {
		bool evflag = (i < nentriesTrig);
                if(evflag) tWeight->GetEntry(i);
		tSample->GetEntry(i);
		for(int j = 0; j < fTunedCentralValue->size(); j++){
			if(evflag) sfTunedCentralValue->push_back(fTunedCentralValue->at(j));
			else sfTunedCentralValue->push_back(0.);
		}
		WXTCV->Fill();

		sfTunedCentralValue->clear();
	}
	tWeight->ResetBranchAddresses();
	delete fweight;
	tSample->Write("",TObject::kOverwrite);
	tSample->ResetBranchAddresses();
	delete fsample;
	}//end of subrun loop
	}//end of run loop
}
