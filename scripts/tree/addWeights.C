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

void addWeights(){
        //weight branches
        int rnTrig;
        int evNTrig;

	vector<double>* All0_weight = new vector<double>();
	vector<double>* All1_weight = new vector<double>();
	vector<double>* All2_weight = new vector<double>();
	vector<double>* All3_weight = new vector<double>();
	vector<double>* All4_weight = new vector<double>();
	vector<double>* fAxFFCCQEshape = new vector<double>();
	vector<double>* fDecayAngMEC = new vector<double>();
	vector<double>* fNormCCCOH = new vector<double>();
	vector<double>* fNorm_NCCOH = new vector<double>();
	vector<double>* fRPA_CCQE = new vector<double>();
	vector<double>* fRootinoFix = new vector<double>();
	vector<double>* fThetaDelta2NRad = new vector<double>();
	vector<double>* fTheta_Delta2Npi = new vector<double>();
	vector<double>* fTunedCentralValue = new vector<double>();
	vector<double>* fVecFFCCQEshape = new vector<double>();
	vector<double>* fXSecShape_CCMEC = new vector<double>();
	vector<double>* flux_horncurrent = new vector<double>();
	vector<double>* flux_expskin = new vector<double>();
	vector<double>* flux_piplus = new vector<double>();
	vector<double>* flux_piminus = new vector<double>();
	vector<double>* flux_kplus = new vector<double>();
	vector<double>* flux_kminus = new vector<double>();
	vector<double>* flux_kzero = new vector<double>();
	vector<double>* flux_pionine = new vector<double>();
	vector<double>* flux_pionqe = new vector<double>();
	vector<double>* flux_piontot = new vector<double>();
	vector<double>* flux_nucine = new vector<double>();
	vector<double>* flux_nucqe = new vector<double>();
	vector<double>* flux_nuctot = new vector<double>();

	//sample branches
        int srnTrig;
        int sevNTrig;

	vector<double>* sAll0_weight = new vector<double>();
	vector<double>* sAll1_weight = new vector<double>();
	vector<double>* sAll2_weight = new vector<double>();
	vector<double>* sAll3_weight = new vector<double>();
	vector<double>* sAll4_weight = new vector<double>();
	vector<double>* sfAxFFCCQEshape = new vector<double>();
	vector<double>* sfDecayAngMEC = new vector<double>();
	vector<double>* sfNormCCCOH = new vector<double>();
	vector<double>* sfNorm_NCCOH = new vector<double>();
	vector<double>* sfRPA_CCQE = new vector<double>();
	vector<double>* sfRootinoFix = new vector<double>();
	vector<double>* sfThetaDelta2NRad = new vector<double>();
	vector<double>* sfTheta_Delta2Npi = new vector<double>();
	vector<double>* sfTunedCentralValue = new vector<double>();
	vector<double>* sfVecFFCCQEshape = new vector<double>();
	vector<double>* sfXSecShape_CCMEC = new vector<double>();
	vector<double>* sflux_horncurrent = new vector<double>();
	vector<double>* sflux_expskin = new vector<double>();
	vector<double>* sflux_piplus = new vector<double>();
	vector<double>* sflux_piminus = new vector<double>();
	vector<double>* sflux_kplus = new vector<double>();
	vector<double>* sflux_kminus = new vector<double>();
	vector<double>* sflux_kzero = new vector<double>();
	vector<double>* sflux_pionine = new vector<double>();
	vector<double>* sflux_pionqe = new vector<double>();
	vector<double>* sflux_piontot = new vector<double>();
	vector<double>* sflux_nucine = new vector<double>();
	vector<double>* sflux_nucqe = new vector<double>();
	vector<double>* sflux_nuctot = new vector<double>();

	int runs = 10;
	int subruns = 20;
	//Loop through runs
	for(int rn = 0; rn < runs; rn++){
		std::cout << "Looping through run " << std::to_string(rn) << std::endl;
	//Loop through run parts (sub runs)
	for(int srn = 0; srn < subruns; srn++){
//		std::cout << "Looping through subrun " << std::to_string(srn) << std::endl;
        //Open file and trees
        string weight_path = "/pnfs/annie/persistent/analysis/v1.3.0/MC/tank/PhaseIITree_0." + std::to_string(rn) + "." + std::to_string(srn) + ".root";
	string sample_path = "/exp/annie/data/users/jminock/temp_add_weights/PhaseIITree_0." + std::to_string(rn) + "." + std::to_string(srn) + ".root";


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
        gSystem->Load("/exp/annie/app/users/jminock/ToolAnalysis/lib/libDataModel.so");
//      gSystem->Load("/exp/annie/app/users/jminock/ToolAnalysis/lib/libDict.so");
        gInterpreter->GenerateDictionary("map<string,vector<double>>", "map;string;vector");
        TTree *tWeight = (TTree*)fweight->Get("phaseIITriggerTree");
        TTree *tSample = (TTree*)fsample->Get("phaseIITriggerTree");

        //Set branch addresses
        tWeight->SetBranchAddress("runNumber",&rnTrig);
        tWeight->SetBranchAddress("eventNumber",&evNTrig);
        tSample->SetBranchAddress("runNumber",&srnTrig);
        tSample->SetBranchAddress("eventNumber",&sevNTrig);

	//Weights
	tWeight->SetBranchAddress("weight_All0_UBGenie",&All0_weight);
	tWeight->SetBranchAddress("weight_All1_UBGenie",&All1_weight);
	tWeight->SetBranchAddress("weight_All2_UBGenie",&All2_weight);
	tWeight->SetBranchAddress("weight_All3_UBGenie",&All3_weight);
	tWeight->SetBranchAddress("weight_All4_UBGenie",&All4_weight);

	tWeight->SetBranchAddress("weight_AxFFCCQEshape_UBGenie",&fAxFFCCQEshape);
	tWeight->SetBranchAddress("weight_DecayAngMEC_UBGenie",&fDecayAngMEC);
	tWeight->SetBranchAddress("weight_NormCCCOH_UBGenie",&fNormCCCOH);
	tWeight->SetBranchAddress("weight_Norm_NCCOH_UBGenie",&fNorm_NCCOH);
	tWeight->SetBranchAddress("weight_RPA_CCQE_UBGenie",&fRPA_CCQE);
	tWeight->SetBranchAddress("weight_RootinoFix_UBGenie",&fRootinoFix);
	tWeight->SetBranchAddress("weight_ThetaDelta2NRad_UBGenie",&fThetaDelta2NRad);
	tWeight->SetBranchAddress("weight_Theta_Delta2Npi_UBGenie",&fTheta_Delta2Npi);
	tWeight->SetBranchAddress("weight_TunedCentralValue_UBGenie",&fTunedCentralValue);
	tWeight->SetBranchAddress("weight_VecFFCCQEshape_UBGenie",&fVecFFCCQEshape);
	tWeight->SetBranchAddress("weight_XSecShape_CCMEC_UBGenie",&fXSecShape_CCMEC); 

	tWeight->SetBranchAddress("weight_horncurrent_FluxUnisim",&flux_horncurrent);
	tWeight->SetBranchAddress("weight_expskin_FluxUnisim",&flux_expskin);
	tWeight->SetBranchAddress("weight_pioninexsec_FluxUnisim",&flux_pionine);
	tWeight->SetBranchAddress("weight_pionqexsec_FluxUnisim",&flux_pionqe);
	tWeight->SetBranchAddress("weight_piontotxsec_FluxUnisim",&flux_piontot);
	tWeight->SetBranchAddress("weight_nucleoninexsec_FluxUnisim",&flux_nucine);
	tWeight->SetBranchAddress("weight_nucleonqexsec_FluxUnisim",&flux_nucqe);
	tWeight->SetBranchAddress("weight_nucleontotxsec_FluxUnisim",&flux_nuctot);
	tWeight->SetBranchAddress("weight_piplus_PrimaryHadronSWCentralSplineVariation",&flux_piplus);
	tWeight->SetBranchAddress("weight_piminus_PrimaryHadronSWCentralSplineVariation",&flux_piminus);
	tWeight->SetBranchAddress("weight_kminus_PrimaryHadronNormalization",&flux_kminus);
	tWeight->SetBranchAddress("weight_kzero_PrimaryHadronSanfordWang",&flux_kzero);
	tWeight->SetBranchAddress("weight_kplus_PrimaryHadronFeynmanScaling",&flux_kplus);

	TBranch *WXA0 = tSample->Branch("weight_All0_UBGenie",&sAll0_weight);
	TBranch *WXA1 = tSample->Branch("weight_All1_UBGenie",&sAll1_weight);
	TBranch *WXA2 = tSample->Branch("weight_All2_UBGenie",&sAll2_weight);
	TBranch *WXA3 = tSample->Branch("weight_All3_UBGenie",&sAll3_weight);
	TBranch *WXA4 = tSample->Branch("weight_All4_UBGenie",&sAll4_weight);

	TBranch *WXAF = tSample->Branch("weight_AxFFCCQEshape_UBGenie",&sfAxFFCCQEshape);
	TBranch *WXDA = tSample->Branch("weight_DecayAngMEC_UBGenie",&sfDecayAngMEC);
	TBranch *WXNCC = tSample->Branch("weight_NormCCCOH_UBGenie",&sfNormCCCOH);
	TBranch *WXNNC = tSample->Branch("weight_NormNCCOH_UBGenie",&sfNorm_NCCOH);
	TBranch *WXRPA = tSample->Branch("weight_RPA_CCQE_UBGenie",&sfRPA_CCQE);
	TBranch *WXRF = tSample->Branch("weight_RootinoFix_UBGenie",&sfRootinoFix);
	TBranch *WXTDR = tSample->Branch("weight_ThetaDelta2NRad_UBGenie",&sfThetaDelta2NRad);
	TBranch *WXTDP = tSample->Branch("weight_Theta_Delta2Npi_UBGenie",&sfTheta_Delta2Npi);
	TBranch *WXTCV = tSample->Branch("weight_TunedCentralValue_UBGenie",&sfTunedCentralValue);
	TBranch *WXVF = tSample->Branch("weight_VecFFCCQEshape_UBGenie",&sfVecFFCCQEshape);
	TBranch *WXXS = tSample->Branch("weight_XSecShape_CCMEC_UBGenie",&sfXSecShape_CCMEC);

	TBranch *WFhc = tSample->Branch("weight_horncurrent_FluxUnisim",&sflux_horncurrent);
	TBranch *WFes = tSample->Branch("weight_expskin_FluxUnisim",&sflux_expskin);
	TBranch *WFpi = tSample->Branch("weight_pioninexsec_FluxUnisim",&sflux_pionine);
	TBranch *WFpq = tSample->Branch("weight_pionqexsec_FluxUnisim",&sflux_pionqe);
	TBranch *WFpt = tSample->Branch("weight_piontotxsec_FluxUnisim",&sflux_piontot);
	TBranch *WFni = tSample->Branch("weight_nucleoninexsec_FluxUnisim",&sflux_nucine);
	TBranch *WFnq = tSample->Branch("weight_nucleonqexsec_FluxUnisim",&sflux_nucqe);
	TBranch *WFnt = tSample->Branch("weight_nucleontotxsec_FluxUnisim",&sflux_nuctot);
	TBranch *WFpp = tSample->Branch("weight_piplus_PrimaryHadronSWCentralSplineVariation",&sflux_piplus);
	TBranch *WFpm = tSample->Branch("weight_piminus_PrimaryHadronSWCentralSplineVariation",&sflux_piminus);
	TBranch *WFkm = tSample->Branch("weight_kminus_PrimaryHadronNormalization",&sflux_kminus);
	TBranch *WFkz = tSample->Branch("weight_kzero_PrimaryHadronSanfordWang",&sflux_kzero);
	TBranch *WFkp = tSample->Branch("weight_kplus_PrimaryHadronFeynmanScaling",&sflux_kplus);

        double muon_m = 105.66;
        Long64_t nentriesTrig = tWeight->GetEntries();
        Long64_t nentriesTrig2 = tSample->GetEntries();
        //fill histograms
	for (Long64_t i = 0; i < nentriesTrig2; i++) {
//	for (Long64_t i = 0; i < 10; i++) {
		bool evflag = (i < nentriesTrig);
                if(evflag) tWeight->GetEntry(i);
		tSample->GetEntry(i);
		for(int j = 0; j < All0_weight->size(); j++){
			if(evflag){
				sAll0_weight->push_back(All0_weight->at(j));
				sAll1_weight->push_back(All1_weight->at(j));
				sAll2_weight->push_back(All2_weight->at(j));
				sAll3_weight->push_back(All3_weight->at(j));
				sAll4_weight->push_back(All4_weight->at(j));
			}
			else{
				sAll0_weight->push_back(0.);
				sAll1_weight->push_back(0.);
				sAll2_weight->push_back(0.);
				sAll3_weight->push_back(0.);
				sAll4_weight->push_back(0.);
			}
		}
		for(int j = 0; j < fAxFFCCQEshape->size(); j++){
			if(evflag) sfAxFFCCQEshape->push_back(fAxFFCCQEshape->at(j));
			else sfAxFFCCQEshape->push_back(0.);
		}
		for(int j = 0; j < fDecayAngMEC->size(); j++){
			if(evflag) sfDecayAngMEC->push_back(fDecayAngMEC->at(j));
			else sfDecayAngMEC->push_back(0.);
		}
		for(int j = 0; j < fNormCCCOH->size(); j++){
			if(evflag) sfNormCCCOH->push_back(fNormCCCOH->at(j));
			else sfNormCCCOH->push_back(0.);
		}
		for(int j = 0; j < fNorm_NCCOH->size(); j++){
			if(evflag) sfNorm_NCCOH->push_back(fNorm_NCCOH->at(j));
			else sfNorm_NCCOH->push_back(0.);
		}
		for(int j = 0; j < fRPA_CCQE->size(); j++){
			if(evflag) sfRPA_CCQE->push_back(fRPA_CCQE->at(j));
			else sfRPA_CCQE->push_back(0.);
		}
		for(int j = 0; j < fRootinoFix->size(); j++){
			if(evflag) sfRootinoFix->push_back(fRootinoFix->at(j));
			else sfRootinoFix->push_back(0.);
		}
		for(int j = 0; j < fThetaDelta2NRad->size(); j++){
			if(evflag) sfThetaDelta2NRad->push_back(fThetaDelta2NRad->at(j));
			else sfThetaDelta2NRad->push_back(0.);
		}
		for(int j = 0; j < fTheta_Delta2Npi->size(); j++){
			if(evflag) sfTheta_Delta2Npi->push_back(fTheta_Delta2Npi->at(j));
			else sfTheta_Delta2Npi->push_back(0.);
		}
		for(int j = 0; j < fTunedCentralValue->size(); j++){
			if(evflag) sfTunedCentralValue->push_back(fTunedCentralValue->at(j));
			else sfTunedCentralValue->push_back(0.);
		}
		for(int j = 0; j < fVecFFCCQEshape->size(); j++){
			if(evflag) sfVecFFCCQEshape->push_back(fVecFFCCQEshape->at(j));
			else sfVecFFCCQEshape->push_back(0.);
		}
		for(int j = 0; j < fXSecShape_CCMEC->size(); j++){
			if(evflag) sfXSecShape_CCMEC->push_back(fXSecShape_CCMEC->at(j));
			else sfXSecShape_CCMEC->push_back(0.);
		}
		for(int j = 0; j < flux_horncurrent->size(); j++){
			if(evflag){
				sflux_horncurrent->push_back(flux_horncurrent->at(j));
				sflux_expskin->push_back(flux_expskin->at(j));
				sflux_pionine->push_back(flux_pionine->at(j));
				sflux_pionqe->push_back(flux_pionqe->at(j));
				sflux_piontot->push_back(flux_piontot->at(j));
				sflux_nucine->push_back(flux_nucine->at(j));
				sflux_nucqe->push_back(flux_nucqe->at(j));
				sflux_nuctot->push_back(flux_nuctot->at(j));
				sflux_piplus->push_back(flux_piplus->at(j));
				sflux_piminus->push_back(flux_piminus->at(j));
				sflux_kminus->push_back(flux_kminus->at(j));
				sflux_kzero->push_back(flux_kzero->at(j));
				sflux_kplus->push_back(flux_kplus->at(j));
			}else{
				sflux_horncurrent->push_back(0.);
				sflux_expskin->push_back(0.);
				sflux_pionine->push_back(0.);
				sflux_pionqe->push_back(0.);
				sflux_piontot->push_back(0.);
				sflux_nucine->push_back(0.);
				sflux_nucqe->push_back(0.);
				sflux_nuctot->push_back(0.);
				sflux_piplus->push_back(0.);
				sflux_piminus->push_back(0.);
				sflux_kminus->push_back(0.);
				sflux_kzero->push_back(0.);
				sflux_kplus->push_back(0.);
			}
		}
		WXA0->Fill();
		WXA1->Fill();
		WXA2->Fill();
		WXA3->Fill();
		WXA4->Fill();
		WXAF->Fill();
		WXDA->Fill();
		WXNCC->Fill();
		WXNNC->Fill();
		WXRPA->Fill();
		WXRF->Fill();
		WXTDR->Fill();
		WXTDP->Fill();
		WXTCV->Fill();
		WXVF->Fill();
		WXXS->Fill();
		WFhc->Fill();
		WFes->Fill();
		WFpi->Fill();
		WFpq->Fill();
		WFpt->Fill();
		WFni->Fill();
		WFnq->Fill();
		WFnt->Fill();
		WFpp->Fill();
		WFpm->Fill();
		WFkm->Fill();
		WFkz->Fill();
		WFkp->Fill();

		sAll0_weight->clear();
		sAll1_weight->clear();
		sAll2_weight->clear();
		sAll3_weight->clear();
		sAll4_weight->clear();
		sfAxFFCCQEshape->clear();
		sfDecayAngMEC->clear();
		sfNormCCCOH->clear();
		sfNorm_NCCOH->clear();
		sfRPA_CCQE->clear();
		sfRootinoFix->clear();
		sfThetaDelta2NRad->clear();
		sfTheta_Delta2Npi->clear();
		sfTunedCentralValue->clear();
		sfVecFFCCQEshape->clear();
		sfXSecShape_CCMEC->clear();
		sflux_horncurrent->clear();
		sflux_expskin->clear();
		sflux_pionine->clear();
		sflux_pionqe->clear();
		sflux_piontot->clear();
		sflux_nucine->clear();
		sflux_nucqe->clear();
		sflux_nuctot->clear();
		sflux_piplus->clear();
		sflux_piminus->clear();
		sflux_kminus->clear();
		sflux_kzero->clear();
		sflux_kplus->clear();
	}
	tWeight->ResetBranchAddresses();
	delete fweight;
	tSample->Write("",TObject::kOverwrite);
	tSample->ResetBranchAddresses();
	delete fsample;
	}//end of subrun loop
	}//end of run loop
}
