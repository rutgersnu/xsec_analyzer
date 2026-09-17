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

void checkBranches(){
        //weight branches
        double mrd_eff;
	double dirt_mu;

        int rnTrig;
        int rnMRD;
        int evNTrig;
        int evNMRD;
        int trigNum;
        int trigword, hasTank, hasMRD, tankMRDCoinc, noveto;
        int nrings;
        double nuE, nuvtxx, nuvtxy, nuvtxz, nupx, nupy, nupz;
        double fslpx, fslpy, fslpz, mcfslE;
        double fslvtxx, fslvtxy, fslvtxz;
        int isCC, isQEL, isRES, isDIS, isCOH, isMEC, mcfslpdg;
        double mcvtxx, mcvtxy, mcvtxz, mcdirx, mcdiry, mcdirz, mcangle, mcmuonE, mctanktracklength, mcmrdtracklength;
        int hasPi0, hasPiP, hasPiM, hasPiPC, hasPiMC, hasKP, hasKM, hasKPC, hasKMC;
        int numMRDTracks;
        int mcentersmrd, mcexitsmrd, mcpenetratesmrd;

        bool mcfv;
	double mcp, mcct;
	double recop, recopc, recoPE;

        int simpleflag, simplefv;
        double simpleenergy, simplecostheta, simplept, simplemrdenergy, simplemrdtrack, simpletanktrack;
        double simplevtxx, simplevtxy, simplevtxz, simplestopvtxx, simplestopvtxy, simplestopvtxz;
        double simplemrdstartx, simplemrdstarty, simplemrdstartz, simplemrdstopx, simplemrdstopy, simplemrdstopz;
	double Qij, PE;
        double recovtxx, recovtxy, recovtxz;

        vector<double>* MRDTrackStartY = new vector<double>();
        vector<double>* hitT = new vector<double>();
        vector<double>* hitPE = new vector<double>();

	//sample branches
	double smrd_eff;
	double sdirt_mu;

        int snoveto;
        double snuvtxx, snuvtxy, snuvtxz;
        double smcfslE;
        int smcfslpdg;
        double smcangle, smcmuonE;

        double ssimpleenergy, ssimplecostheta;
        double ssimplevtxx, ssimplevtxy, ssimplevtxz;
	double sQij, sPE;

        vector<double>* sMRDTrackStartY = new vector<double>();
        vector<double>* shitT = new vector<double>();
        vector<double>* shitPE = new vector<double>();

	int runs = 1;
	int subruns = 1;
	//Loop through runs
	for(int rn = 0; rn < runs; rn++){
		std::cout << "Looping through run " << std::to_string(rn) << std::endl;
	//Loop through run parts (sub runs)
	for(int srn = 0; srn < subruns; srn++){
//		std::cout << "Looping through subrun " << std::to_string(srn) << std::endl;
        //Open file and trees
        string weight_path = "/pnfs/annie/persistent/users/jminock/v1_3_4_world_stv_ntuples/PhaseIITree_0.0.0.root";
	string sample_path = "/exp/annie/app/users/jminock/ToolAnalysis/PhaseIITree_0.0.0.root";


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
        TFile *fsample = new TFile(sample_path.c_str(),"read");
//        gSystem->Load("/exp/annie/app/users/jminock/ToolAnalysis/lib/libDataModel.so");
//      gSystem->Load("/exp/annie/app/users/jminock/ToolAnalysis/lib/libDict.so");
//        gInterpreter->GenerateDictionary("map<string,vector<double>>", "map;string;vector");
        TTree *tWeight = (TTree*)fweight->Get("phaseIITriggerTree");
        TTree *tSample = (TTree*)fsample->Get("phaseIITriggerTree");

	//Weights
        tWeight->SetBranchAddress("NoVeto",&noveto);
        tWeight->SetBranchAddress("MRDTrackStartY",&MRDTrackStartY);
        tWeight->SetBranchAddress("trueNuIntxVtx_X",&nuvtxx);
        tWeight->SetBranchAddress("trueNuIntxVtx_Y",&nuvtxy);
        tWeight->SetBranchAddress("trueNuIntxVtx_Z",&nuvtxz);
        tWeight->SetBranchAddress("trueFSLEnergy",&mcfslE);
        tWeight->SetBranchAddress("trueFSLPdg",&mcfslpdg);
        tWeight->SetBranchAddress("simpleRecoEnergy",&simpleenergy);
        tWeight->SetBranchAddress("simpleRecoVtxX",&simplevtxx);
        tWeight->SetBranchAddress("simpleRecoVtxY",&simplevtxy);
        tWeight->SetBranchAddress("simpleRecoVtxZ",&simplevtxz);
        tWeight->SetBranchAddress("simpleRecoCosTheta",&simplecostheta);
	tWeight->SetBranchAddress("Qij",&Qij);
	tWeight->SetBranchAddress("promptMuonTotalPE",&PE);
        tWeight->SetBranchAddress("trueMuonEnergy",&mcmuonE);
        tWeight->SetBranchAddress("trueAngle",&mcangle);
	tWeight->SetBranchAddress("MRDEff",&mrd_eff);
	tWeight->SetBranchAddress("DirtMu",&dirt_mu);
	tWeight->SetBranchAddress("hitT",&hitT);
	tWeight->SetBranchAddress("hitPE",&hitPE);
	tSample->SetBranchAddress("hitT",&shitT);
	tSample->SetBranchAddress("hitPE",&shitPE);
	tSample->SetBranchAddress("MRDEff",&smrd_eff);
	tSample->SetBranchAddress("DirtMu",&sdirt_mu);
        tSample->SetBranchAddress("NoVeto",&snoveto);
        tSample->SetBranchAddress("MRDTrackStartY",&sMRDTrackStartY);
        tSample->SetBranchAddress("trueNuIntxVtx_X",&snuvtxx);
        tSample->SetBranchAddress("trueNuIntxVtx_Y",&snuvtxy);
        tSample->SetBranchAddress("trueNuIntxVtx_Z",&snuvtxz);
        tSample->SetBranchAddress("trueFSLEnergy",&smcfslE);
        tSample->SetBranchAddress("trueFSLPdg",&smcfslpdg);
        tSample->SetBranchAddress("simpleRecoEnergy",&ssimpleenergy);
        tSample->SetBranchAddress("simpleRecoVtxX",&ssimplevtxx);
        tSample->SetBranchAddress("simpleRecoVtxY",&ssimplevtxy);
        tSample->SetBranchAddress("simpleRecoVtxZ",&ssimplevtxz);
        tSample->SetBranchAddress("simpleRecoCosTheta",&ssimplecostheta);
	tSample->SetBranchAddress("Qij",&sQij);
	tSample->SetBranchAddress("promptMuonTotalPE",&sPE);
        tSample->SetBranchAddress("trueMuonEnergy",&smcmuonE);
        tSample->SetBranchAddress("trueAngle",&smcangle);

        double muon_m = 105.66;
        Long64_t nentriesTrig = tWeight->GetEntries();
        Long64_t nentriesTrig2 = tSample->GetEntries();
	std::cout << nentriesTrig << " " << nentriesTrig2 << std::endl;
        //fill histograms
	for (Long64_t i = 0; i < nentriesTrig2; i++) {
//	for (Long64_t i = 0; i < 10; i++) {
		bool evflag = (i < nentriesTrig);
                if(evflag) tWeight->GetEntry(i);
		tSample->GetEntry(i);
		if(mrd_eff != smrd_eff) std::cout << "MRD "<<mrd_eff<<" " <<smrd_eff << std::endl;
		if(dirt_mu != sdirt_mu) std::cout << "DIRT " <<dirt_mu << " " <<sdirt_mu<< std::endl;
		if(noveto != snoveto) std::cout << "1 MISMATCH" << std::endl;
		if(nuvtxy != snuvtxy) std::cout << "nuvtxy "<<nuvtxy<< " " <<snuvtxy << std::endl;
//		if(mcfslE != smcfslE) std::cout << "MCFSLE " << mcfslE << " " << smcfslE << std::endl;
		if(mcfslpdg != smcfslpdg) std::cout << "MCFSLPDG " << mcfslpdg << " " << smcfslpdg << std::endl;
		if(simpleenergy != ssimpleenergy) std::cout << "SIMPLEENERGY " << simpleenergy << " " << ssimpleenergy << std::endl;
		if(simplevtxy != ssimplevtxy) std::cout << "SIMPLEVTX " << simplevtxy << " " << ssimplevtxy << std::endl;
		if(simplecostheta != ssimplecostheta) std::cout << "7 MISMATCH" << std::endl;
		if(Qij != sQij) std::cout << "QIJ "<< Qij << " " << sQij << std::endl;
		if(PE != sPE) std::cout << "PE " << PE << " " << sPE << std::endl;
/*		if(mcangle != smcangle) std::cout << "10 MISMATCH" << std::endl;
		for(int j = 0; j < MRDTrackStartY->size(); ++j){
			if(MRDTrackStartY->at(j) != sMRDTrackStartY->at(j)) std::cout << "MRD " << MRDTrackStartY->at(j) << " " << sMRDTrackStartY->at(j) <<std::endl;
		}
		for(int j = 0; j < hitT->size(); ++j){
			if(hitT->at(j) != shitT->at(j)) std::cout << "HITT " << hitT->at(j) << " " << shitT->at(j) <<std::endl;
		}
		for(int j = 0; j < hitPE->size(); ++j){
			if(hitPE->at(j) != shitPE->at(j)) std::cout << "HITPE " << hitPE->at(j) << " " << shitPE->at(j) <<std::endl;
		}
*/


	}
	tWeight->ResetBranchAddresses();
	delete fweight;
	tSample->ResetBranchAddresses();
	delete fsample;
	}//end of subrun loop
	}//end of run loop
}
