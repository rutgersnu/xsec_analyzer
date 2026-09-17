#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "THStack.h"
#include "TLegend.h"
#include "TF1.h"
#include "TF2.h"
#include "TLine.h"
#include "TMath.h"
#include "TRandom3.h"
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

void breakCSV(string line, vector<double> &tokens){
	//change line into stream
	std::stringstream line_in(line);
	std::string temp_token;
	double token;
	while (line_in.good()){
		std::getline(line_in, temp_token, ','); //break up by comma
		std::stringstream to_token(temp_token);
		to_token >> token;          //turn into double
		tokens.push_back(token);    //save to array
	}
}

int findBin(double Y, int iter, vector<double> const &bins){
	if(iter == bins.size()) return 0;
	if(0.0001 > std::abs(bins[iter] - Y)) return iter;
	else return findBin(Y, iter+1, bins);
}

void DVShift(){
        //CV branches
        int trigword, hasTank, hasMRD, tankMRDCoinc, noveto;
        int isCC, mcfslpdg, nentry;

        bool mcfv;
	double mcp, mcct;
	bool mc_no_mesons, no_followers;
	bool recofv, recoMRDInc, recoMRD0pi, reco0pi;
	double recop, recopc, recoPE, recoE, tracklength, mctracklength, promptPE;
	double mcvtxx, mcvtxy, mcvtxz, recovtxx, recovtxy, recovtxz;

        int simpleflag;
        double simplecostheta;
	double mrd_eff;
	double dirt_muon;

        vector<double>* MRDTrackStartX = new vector<double>();
        vector<double>* MRDTrackStartY = new vector<double>();
	vector<bool>* MRDStop = new vector<bool>();

	//DV branches
        int dtrigword, dhasTank, dhasMRD, dtankMRDCoinc, dnoveto;
        int disCC, dmcfslpdg, dnentry;

        bool dmcfv;
	double dmcp, dmcct;
	bool dmc_no_mesons, dno_followers;
	bool drecofv, drecoMRDInc, drecoMRD0pi, dreco0pi;
	double drecop, drecopc, drecoPE, dtracklength, dmctracklength, dpromptPE;
	double dmcvtxx, dmcvtxy, dmcvtxz, drecovtxx, drecovtxy, drecovtxz;

        int dsimpleflag;
        double dsimplecostheta;
	double dmrd_eff;
	double ddirt_muon;

	//read in calibration file
	std::vector<double> bins_front;
	std::vector<double> bins_TL;
	std::vector<double> factorX;
	std::vector<double> factorY;
	std::vector<double> factorTL;
	std::vector<double> uncsX;
	std::vector<double> uncsY;
	std::vector<double> uncsTL;
	string bins_front_str = "";
	string bins_TL_str = "";
	string factorX_str = "";
	string uncsX_str = "";
	string factorY_str = "";
	string uncsY_str = "";
	string factorTL_str = "";
	string uncsTL_str = "";

	string mrd_cal_file = "/exp/annie/app/users/jminock/ANNIE_AuxFiles/MRDEffUncYlarge.txt";
	std::ifstream cal_file(mrd_cal_file.c_str(), ios::in);
	cal_file >> bins_front_str;
	cal_file >> factorY_str;
	cal_file >> uncsY_str;
	cal_file.close();

	breakCSV(bins_front_str, bins_front);
	breakCSV(factorY_str, factorY);
	breakCSV(uncsY_str, uncsY);

	TRandom3 rnd;
	rnd.SetSeed(1337);

	int runs = 60;
	//Loop through runs
	for(int rn = 1; rn <= runs; rn++){
		std::cout << "Looping through run " << std::to_string(rn) << std::endl;
	//Loop through run parts (sub runs)
	

        //Open file and trees
        TFile *fcv = new TFile("/pnfs/annie/persistent/users/jminock/v1_3_4_world_stv_DV_ntuples/PhaseIITree_CV_stv_ntuple.root","read");
	string file_name = "/exp/annie/data/users/jminock/temp_add_branches/PhaseIITree_DV"+std::to_string(rn)+"_stv_ntuple.root";
        TFile *fdv = new TFile(file_name.c_str(),"create");
//        gSystem->Load("/exp/annie/app/users/jminock/ToolAnalysis/lib/libDataModel.so");
//      gSystem->Load("/exp/annie/app/users/jminock/ToolAnalysis/lib/libDict.so");
//        gInterpreter->GenerateDictionary("map<string,vector<double>>", "map;string;vector");
        TTree *tCV = (TTree*)fcv->Get("phaseIITriggerTree");
        TTree *tDV = new TTree("phaseIITriggerTree", "phaseIITriggerTree");

	//check if calibration uncertainty file exists
	if(gSystem->AccessPathName(mrd_cal_file.c_str())){
		std::cout << "WARNING: " << mrd_cal_file << " does not exist. Stopping." << std::endl;
		return false;
	}

//	TRandom3 rnd;
//	rnd.SetSeed(rn);

        //Set branch addresses
        tCV->SetBranchAddress("trigword",&trigword);
        tCV->SetBranchAddress("HasTank",&hasTank);
        tCV->SetBranchAddress("HasMRD",&hasMRD);
        tCV->SetBranchAddress("TankMRDCoinc",&tankMRDCoinc);
        tCV->SetBranchAddress("NoVeto",&noveto);
        tCV->SetBranchAddress("trueCC",&isCC);
        tCV->SetBranchAddress("trueFSLPdg",&mcfslpdg);
        tCV->SetBranchAddress("trueFV",&mcfv);
        tCV->SetBranchAddress("trueMuonMomentum",&mcp);
        tCV->SetBranchAddress("trueCosTheta",&mcct);
        tCV->SetBranchAddress("trueNuIntxVtx_X",&mcvtxx);
        tCV->SetBranchAddress("trueNuIntxVtx_Y",&mcvtxy);
        tCV->SetBranchAddress("trueNuIntxVtx_Z",&mcvtxz);
        tCV->SetBranchAddress("trueTrackLengthInMRD",&mctracklength);
        tCV->SetBranchAddress("true_no_mesons",&mc_no_mesons);
        tCV->SetBranchAddress("true_no_followers",&no_followers);
        tCV->SetBranchAddress("recoInc_contained_in_MRD",&recoMRDInc);
        tCV->SetBranchAddress("reco0pi_contained_in_MRD",&recoMRD0pi);
        tCV->SetBranchAddress("recoFV",&recofv);
        tCV->SetBranchAddress("recoPE",&recoPE);
        tCV->SetBranchAddress("simpleRecoTrackLengthInMRD",&tracklength);
        tCV->SetBranchAddress("reco_0pi",&reco0pi);
        tCV->SetBranchAddress("MRDEff",&mrd_eff);
        tCV->SetBranchAddress("DirtMu",&dirt_muon);
        tCV->SetBranchAddress("mcEntryNumber",&nentry);
        tCV->SetBranchAddress("simpleRecoFlag",&simpleflag);
        tCV->SetBranchAddress("simpleRecoMomentumCor",&recopc);
        tCV->SetBranchAddress("simpleRecoEnergy",&recoE);
        tCV->SetBranchAddress("simpleRecoCosTheta",&simplecostheta);
        tCV->SetBranchAddress("simpleRecoVtxX",&recovtxx);
        tCV->SetBranchAddress("simpleRecoVtxY",&recovtxy);
        tCV->SetBranchAddress("simpleRecoVtxZ",&recovtxz);
        tCV->SetBranchAddress("promptMuonTotalPE",&promptPE);

//        tCV->SetBranchAddress("MRDTrackStartX",&MRDTrackStartX);
//        tCV->SetBranchAddress("MRDTrackStartY",&MRDTrackStartY);
//        tCV->SetBranchAddress("MRDStop",&MRDStop);

	//DV Branches
        TBranch *Dtw = tDV->Branch("trigword",&dtrigword);
        TBranch *Dhs = tDV->Branch("HasTank",&dhasTank);
        TBranch *Dhm = tDV->Branch("HasMRD",&dhasMRD);
        TBranch *Dtmc = tDV->Branch("TankMRDCoinc",&dtankMRDCoinc);
        TBranch *Dnv = tDV->Branch("NoVeto",&dnoveto);
        TBranch *Dtcc = tDV->Branch("trueCC",&disCC);
        TBranch *Dtpdg = tDV->Branch("trueFSLPdg",&dmcfslpdg);
        TBranch *Dtfv = tDV->Branch("trueFV",&dmcfv);
        TBranch *Dtmm = tDV->Branch("trueMuonMomentum",&dmcp);
        TBranch *Dtct = tDV->Branch("trueCosTheta",&dmcct);
        TBranch *Dtvx = tDV->Branch("trueNuIntxVtx_X",&dmcvtxx);
        TBranch *Dtvy = tDV->Branch("trueNuIntxVtx_Y",&dmcvtxy);
        TBranch *Dtvz = tDV->Branch("trueNuIntxVtx_Z",&dmcvtxz);
        TBranch *Dttl = tDV->Branch("trueTrackLengthInMRD",&dmctracklength);
        TBranch *Dtnm = tDV->Branch("true_no_mesons",&dmc_no_mesons);
        TBranch *Dtnf = tDV->Branch("true_no_followers",&dno_followers);
        TBranch *Dric = tDV->Branch("recoInc_contained_in_MRD",&drecoMRDInc);
        TBranch *Dr0c = tDV->Branch("reco0pi_contained_in_MRD",&drecoMRD0pi);
        TBranch *Drfv = tDV->Branch("recoFV",&drecofv);
        TBranch *Drpe = tDV->Branch("recoPE",&drecoPE);
	TBranch *Dtl  = tDV->Branch("simpleRecoTrackLengthInMRD",&dtracklength);
        TBranch *Dr0pi = tDV->Branch("reco_0pi",&dreco0pi);
        TBranch *Dmrd = tDV->Branch("MRDEff",&dmrd_eff);
        TBranch *Ddm = tDV->Branch("DirtMu",&ddirt_muon);
        TBranch *Dmen = tDV->Branch("mcEntryNumber",&dnentry);
        TBranch *Drf = tDV->Branch("simpleRecoFlag",&dsimpleflag);
        TBranch *Drmo = tDV->Branch("simpleRecoMomentumCor",&drecopc);
        TBranch *Drct = tDV->Branch("simpleRecoCosTheta",&dsimplecostheta);
        TBranch *Drvx = tDV->Branch("simpleRecoVtxX",&drecovtxx);
        TBranch *Drvy = tDV->Branch("simpleRecoVtxY",&drecovtxy);
        TBranch *Drvz = tDV->Branch("simpleRecoVtxZ",&drecovtxz);
        TBranch *Dpe = tDV->Branch("promptMuonTotalPE",&dpromptPE);

        double muon_m = 105.7;
	double threshold = 0.0001;
	double tank_unc = 16.78;//from Luis
	double mrd_unc = 31.05; //from Luis
	double unc_corr = 0.0;  //assume uncorrelated between tank and MRD
	double E_unc = std::sqrt(tank_unc*tank_unc + mrd_unc*mrd_unc + 2*unc_corr*tank_unc*mrd_unc);
        Long64_t nentriesCV = tCV->GetEntries();
        //fill histograms
	//std::cout << "# of events: " << nentriesCV << std::endl;
	for (Long64_t i = 0; i < nentriesCV; i++) {
                tCV->GetEntry(i);
//		if(i%200000 == 0) std::cout << i << std::endl;
		//Keep the same
	        dtrigword = trigword;
		dhasTank = hasTank;
		dhasMRD = hasMRD;
		dtankMRDCoinc = tankMRDCoinc;
		dnoveto = noveto;
        	disCC = isCC;
		dmcfslpdg = mcfslpdg;
		dnentry = nentry;
		dmcfv = mcfv;
		dmcp = mcp;
		dmcct = mcct;
		dmc_no_mesons = mc_no_mesons;
		dno_followers = no_followers;
		drecofv = recofv;
		drecoMRDInc = recoMRDInc;
		drecoMRD0pi = recoMRD0pi;
		dreco0pi = reco0pi;
		drecoPE = recoPE;
		dsimpleflag = simpleflag;
		dsimplecostheta = simplecostheta;
		dtracklength = tracklength;
		dmctracklength = mctracklength;
		dmcvtxx = mcvtxx;
		drecovtxx = recovtxx;
		dmcvtxy = mcvtxy;
		drecovtxy = recovtxy;
		dmcvtxz = mcvtxz;
		drecovtxz = recovtxz;
		dpromptPE = promptPE;


		dmrd_eff = mrd_eff;          //moving factor from DV to own weighted category
		ddirt_muon = dirt_muon;      //moving factor from DV to own weighted category
		//Apply shifts
		if(recoE < 0){//if invalid energy, skip
			drecopc = recopc;
		} else {
			double Emod = recoE*(1.+0.0589-0.0187); //from Luis
			double newE = recoE + rnd.Gaus(0.,E_unc);
			//convert KE to corrected momentum
			drecopc = std::sqrt((newE+muon_m)*(newE+muon_m) - muon_m*muon_m)*0.82 + 160.;
		}

		//MRD Eff
                //Establish MRD X/Y variable to use if there are multiple tracks
/*		int ntracks = MRDTrackStartY->size();
		if(ntracks <= 0) {
			dmrd_eff = 0.0;//No tracks so no uncertainty and no weight
		} else {
			double MRD_X = -9999;
			double MRD_Y = -9999;
			for(int j = 0; j < ntracks; j++){
				if(MRDStop->at(j)) MRD_X = MRDTrackStartX->at(j); //confirm its stopping track
				if(MRDStop->at(j)) MRD_Y = MRDTrackStartY->at(j);
			}
			//Establish which bin the event falls into
			int binX = findBin(MRD_X, 0, bins_front);
			int binY = findBin(MRD_Y, 0, bins_front);
			int binTL = findBin(tracklength*100, 0, bins_TL);
			//Assign uncertainty based on bin
			double mrd_X = factorX[binX];
			double mrd_Y = factorY[binY];
			double mrd_TL = factorTL[binTL];
			double mrd_uncX = uncsX[binX]; //Corresponding uncertainty for event
			double mrd_uncY = uncsY[binY];
			double mrd_uncTL = uncsTL[binTL];

			double mrd_unc = mrd_eff*std::sqrt((mrd_uncX/mrd_X)*(mrd_uncX/mrd_X) + (mrd_uncY/mrd_Y)*(mrd_uncY/mrd_Y) + (mrd_uncTL/mrd_TL)*(mrd_uncTL/mrd_TL));
			dmrd_eff = rnd.Gaus(mrd_eff,mrd_unc);
		}
*/
		//MRD Eff
                //Establish MRD Y variable to use if there are multiple tracks
/*		if(mrd_eff < threshold) {
			dmrd_eff = 0.0;//No tracks so no uncertainty and no weight
		} else {
			//Establish which bin the event falls into
			int binY = findBin(mrd_eff, 0, factorY);
			//Assign uncertainty based on bin
			double mrd_uncY = uncsY[binY];
			//dmrd_eff = rnd.Gaus(mrd_eff,mrd_uncY);
			dmrd_eff = rnd.Gaus(mrd_eff,std::abs(mrd_eff-1.));
			while(dmrd_eff < 0.) dmrd_eff = rnd.Gaus(mrd_eff,std::abs(mrd_eff-1.));//if negative, try again
		}
		//Dirt Muon Correction
		if(std::abs(dirt_muon - 1) < threshold){
			ddirt_muon = 1.; //not a dirt muon
		} else { //is a dirt muon
			//0.028 is calculated uncertainty from total # of events / POT
			ddirt_muon = rnd.Gaus(dirt_muon, 0.0028);
			while(ddirt_muon < 0.) ddirt_muon = rnd.Gaus(dirt_muon, 0.0028);//if negative, try again
		}
*/
		tDV->Fill();
	}
	tCV->ResetBranchAddresses();
	delete fcv;
	tDV->Write();
	tDV->ResetBranchAddresses();
	delete fdv;
	}
}
