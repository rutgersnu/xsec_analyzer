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

//Script that strips excess weights, variables, and trees out of ANNIE MC files
//Used for numuMC files
void slim2(){

	int runs = 5000;//5000;
	int subruns = 1;
	//Loop through runs
	for(int rn = 3000; rn < runs; rn++){
		std::cout << "Looping through run " << std::to_string(rn) << std::endl;
        //Open file and trees
        string file_path = "/pnfs/annie/persistent/users/jminock/v1_3_4_world_ntuples/PhaseIITree_0." + std::to_string(rn) + ".0.root";
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

	T->SetBranchStatus("*", 0);
	//Weights
	T->SetBranchStatus("weight_All_UBGenie", 1);
	T->SetBranchStatus("weight_AxFFCCQEshape_UBGenie", 1);
	T->SetBranchStatus("weight_DecayAngMEC_UBGenie", 1);
	T->SetBranchStatus("weight_NormCCCOH_UBGenie", 1);
	T->SetBranchStatus("weight_NormNCCOH_UBGenie", 1);
	T->SetBranchStatus("weight_RPA_CCQE_UBGenie", 1);
	T->SetBranchStatus("weight_RootinoFix_UBGenie", 1);
	T->SetBranchStatus("weight_ThetaDelta2NRad_UBGenie", 1);
	T->SetBranchStatus("weight_Theta_Delta2Npi_UBGenie", 1);
	T->SetBranchStatus("weight_TunedCentralValue_UBGenie", 1);
	T->SetBranchStatus("weight_VecFFCCQEshape_UBGenie", 1);
	T->SetBranchStatus("weight_XSecShape_CCMEC_UBGenie", 1);
	T->SetBranchStatus("weight_flux_all", 1);
	//Truth
	T->SetBranchStatus("trueCC", 1);
	T->SetBranchStatus("trueNC", 1);
	T->SetBranchStatus("trueQEL", 1);
	T->SetBranchStatus("trueRES", 1);
	T->SetBranchStatus("trueMEC", 1);
	T->SetBranchStatus("trueCOH", 1);
	T->SetBranchStatus("trueDIS", 1);
	T->SetBranchStatus("trueFSLPdg", 1);
	T->SetBranchStatus("truePi0", 1);
	T->SetBranchStatus("truePiPlusCher", 1);
	T->SetBranchStatus("truePiMinusCher", 1);
	T->SetBranchStatus("trueFollowerParentPDG", 1);
	T->SetBranchStatus("trueKPlusCher", 1);
	T->SetBranchStatus("trueKMinusCher", 1);
	T->SetBranchStatus("trueFV", 1);
	T->SetBranchStatus("true_no_mesons", 1);
	T->SetBranchStatus("true_no_followers", 1);
	T->SetBranchStatus("trueAngle", 1);
	T->SetBranchStatus("trueCosTheta", 1);
	T->SetBranchStatus("trueMuonMomentum", 1);
	T->SetBranchStatus("trueNuIntxVtx_X", 1);
	T->SetBranchStatus("trueNuIntxVtx_Y", 1);
	T->SetBranchStatus("trueNuIntxVtx_Z", 1);
	T->SetBranchStatus("trueTrackLengthInMRD", 1);
	//Reco
	T->SetBranchStatus("trigword", 1);
	T->SetBranchStatus("HasTank", 1);
	T->SetBranchStatus("HasMRD", 1);
	T->SetBranchStatus("TankMRDCoinc", 1);
	T->SetBranchStatus("NoVeto", 1);
	T->SetBranchStatus("simpleRecoFlag", 1);
	T->SetBranchStatus("reco0pi_contained_in_MRD", 1);
	T->SetBranchStatus("recoInc_contained_in_MRD", 1);
	T->SetBranchStatus("recoFV", 1);
	T->SetBranchStatus("recoPE", 1);
	T->SetBranchStatus("reco_0pi", 1);
	T->SetBranchStatus("simpleRecoCosTheta", 1);
	T->SetBranchStatus("simpleRecoMomentumCor", 1);
	T->SetBranchStatus("simpleRecoTrackLengthInMRD", 1);
	T->SetBranchStatus("simpleRecoVtxX", 1);
	T->SetBranchStatus("simpleRecoVtxY", 1);
	T->SetBranchStatus("simpleRecoVtxZ", 1);
	T->SetBranchStatus("numMRDTracks", 1);
	T->SetBranchStatus("MRDStop", 1);
	T->SetBranchStatus("promptMuonTotalPE", 1);
	T->SetBranchStatus("Qij", 1);
	//Misc
	T->SetBranchStatus("MRDEff", 1);
	T->SetBranchStatus("DirtMu", 1);
	T->SetBranchStatus("weight_MRDUnc", 1);
	T->SetBranchStatus("weight_DirtUnc", 1);
	T->SetBranchStatus("mcEntryNumber", 1);

        TFile *fnew = new TFile(out_path.c_str(),"recreate");
	auto Tnew = T->CloneTree();

	fnew->Write();
	delete f;
	delete fnew;
	}
}
