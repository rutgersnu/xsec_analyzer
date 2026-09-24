/**
 * annie_stv_prep
 *
 * Prepare ANNIE PhaseIITrees for use in xsec_analyzer:
 *
 * * Compute complex selection flags/variables that define bins
 * * Drop some branches for space
 * * Merge weights into products
 * * Add weights for MRD and dirt corrections.
 *
 * ATM 2026/09/18, based on J. Minock's stvPrep.
 */

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <cstdlib>
#include "TFile.h"
#include "TTree.h"
#include "EventCategory.hh"

static const double muon_m = 105.7;  // MeV/c^2


bool FidVol(double x, double y, double z) {
  static const double radius   = 100.;  //cm
  static const double y_min    = -100.; //cm
  static const double y_max    = 100.;  //cm
  static const double z_center = 168.1; //cm
  static const double y_offset = 14.46; //cm

  return (
    y + y_offset > y_min &&
    y + y_offset < y_max &&
    radius > std::sqrt((z - z_center)*(z - z_center) + x*x)
  );
}


void breakCSV(std::string line, std::vector<double> &tokens){
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


int findBin(double Y, int iter, std::vector<double> const &bins){
  if(iter == (int)bins.size()) return 0;
  if(Y >= bins.at(iter)) return findBin(Y, iter+1, bins);
  else return iter;
}


void stvPrep(TString& infile, TString& outfile, TString* weightfile=nullptr){
  TFile f_in(infile.Data(), "read");
  assert(f_in.IsOpen() && !f_in.IsZombie());

  TTree* t = nullptr;
  f_in.GetObject("phaseIITriggerTree", t);
  assert(t && t->GetEntries());

  // Prune branches
  std::vector<const char*> branches = {
    "runNumber",
    "subrunNumber",
    "runType",
    "startTime",
    "eventNumber",
    "eventTimeTank",
    "eventTimeMRD",
    "nhits",
    "trigword",
    "HasTank",
    "HasMRD",
    "TankMRDCoinc",
    "NoVeto",
    "Extended",
    "beam_pot",
    "beam_ok",
    "filter",
    "numMRDTracks",
    "MRDEnergyLoss",
    "MRDEnergyLossError",
    "MRDTrackStartX",
    "MRDTrackStartY",
    "MRDTrackStartZ",
    "MRDTrackStopX",
    "MRDTrackStopY",
    "MRDTrackStopZ",
    "MRDSide",
    "MRDStop",
    "MRDThrough",
    "promptMuonTotalPE",
    "simpleRecoFlag",
    "simpleRecoEnergy",
    "simpleRecoVtxX",
    "simpleRecoVtxY",
    "simpleRecoVtxZ",
    "simpleRecoStopVtxX",
    "simpleRecoStopVtxY",
    "simpleRecoStopVtxZ",
    "simpleRecoCosTheta",
    "simpleRecoPt",
    "simpleRecoFV",
    "simpleRecoMrdEnergyLoss",
    "simpleRecoTrackLengthInMRD",
    "simpleRecoMRDStartX",
    "simpleRecoMRDStartY",
    "simpleRecoMRDStartZ",
    "simpleRecoMRDStopX",
    "simpleRecoMRDStopY",
    "simpleRecoMRDStopZ",
    "simpleRecoTrackLengthInTank",
    "triggerNumber",
    "mcEntryNumber",
    "trueVtxX",
    "trueVtxY",
    "trueVtxZ",
    "trueVtxTime",
    "trueDirX",
    "trueDirY",
    "trueDirZ",
    "trueAngle",
    "truePhi",
    "trueMuonEnergy",
    "truePrimaryPdg",
    "trueTrackLengthInWater",
    "trueTrackLengthInMRD",
    "trueMultiRing",
    "Pi0Count",
    "PiPlusCount",
    "PiMinusCount",
    "K0Count",
    "KPlusCount",
    "KMinusCount",
    "truePrimaryPdgs",
    "trueNeutrinoEnergy",
    "trueNuPDG",
    "trueNeutrinoMomentum_X",
    "trueNeutrinoMomentum_Y",
    "trueNeutrinoMomentum_Z",
    "trueNuIntxVtx_X",
    "trueNuIntxVtx_Y",
    "trueNuIntxVtx_Z",
    "trueNuIntxVtx_T",
    "trueFSLVtx_X",
    "trueFSLVtx_Y",
    "trueFSLVtx_Z",
    "trueFSLMomentum_X",
    "trueFSLMomentum_Y",
    "trueFSLMomentum_Z",
    "trueFSLTime",
    "trueFSLMass",
    "trueFSLPdg",
    "trueFSLEnergy",
    "trueQ2",
    "trueCC",
    "trueNC",
    "trueQEL",
    "trueRES",
    "trueDIS",
    "trueCOH",
    "trueMEC",
    "trueNeutrons",
    "trueProtons",
    "truePi0",
    "truePiPlus",
    "truePiPlusCher",
    "truePiMinus",
    "truePiMinusCher",
    "trueKPlus",
    "trueKPlusCher",
    "trueKMinus",
    "trueKMinusCher",
    "trueq0",
    "trueq3",
    "trueTargetZ",
    "weight_All0_UBGenie",
    "weight_All1_UBGenie",
    "weight_All2_UBGenie",
    "weight_All3_UBGenie",
    "weight_All4_UBGenie",
    "weight_All5_UBGenie",
    "weight_AxFFCCQEshape_UBGenie",
    "weight_DecayAngMEC_UBGenie",
    "weight_NormCCCOH_UBGenie",
    "weight_NormNCCOH_UBGenie",
    "weight_RPA_CCQE_UBGenie",
    "weight_RootinoFix_UBGenie",
    "weight_ThetaDelta2NRad_UBGenie",
    "weight_Theta_Delta2Npi_UBGenie",
    "weight_TunedCentralValue_UBGenie",
    "weight_VecFFCCQEshape_UBGenie",
    "weight_XSecShape_CCMEC_UBGenie",
    "weight_horncurrent_FluxUnisim",
    "weight_expskin_FluxUnisim",
    "weight_pioninexsec_FluxUnisim",
    "weight_pionqexsec_FluxUnisim",
    "weight_piontotxsec_FluxUnisim",
    "weight_nucleoninexsec_FluxUnisim",
    "weight_nucleonqexsec_FluxUnisim",
    "weight_nucleontotxsec_FluxUnisim",
    "weight_piplus_PrimaryHadronSWCentralSplineVariation",
    "weight_piminus_PrimaryHadronSWCentralSplineVariation",
    "weight_kminus_PrimaryHadronNormalization",
    "weight_kzero_PrimaryHadronSanfordWang",
    "weight_kplus_PrimaryHadronFeynmanScaling"
  };

  t->SetBranchStatus("*", 0);

  for (const char* k : branches) {
    t->SetBranchStatus(k, 1);
  }

  // Branches to read
  double nuvtxx, nuvtxy, nuvtxz;
  double mcangle, mcmuonE;
  int simpleflag, trueNuPDG;
  int trueCC, trueNC, trueQEL, trueRES, trueDIS, trueCOH, trueMEC;
  double simpleenergy, simplevtxx, simplevtxy, simplevtxz, PE;
  int hasPi0, hasPiP, hasPiM, hasPiPC, hasPiMC, hasKP, hasKM, hasKPC, hasKMC;
  int numMRDTracks;
  //double Qij;
  std::vector<double>* MRDTrackStartY = new std::vector<double>();
  std::vector<bool>* MRDStop = new std::vector<bool>();
  //std::vector<int>* mcFolPPDG = new std::vector<int>();

  t->SetBranchAddress("simpleRecoFlag",&simpleflag);
  t->SetBranchAddress("numMRDTracks",&numMRDTracks);
  t->SetBranchAddress("MRDTrackStartY",&MRDTrackStartY);
  t->SetBranchAddress("MRDStop",&MRDStop);
  t->SetBranchAddress("trueNuPDG",&trueNuPDG);
  t->SetBranchAddress("trueNuIntxVtx_X",&nuvtxx);
  t->SetBranchAddress("trueNuIntxVtx_Y",&nuvtxy);
  t->SetBranchAddress("trueNuIntxVtx_Z",&nuvtxz);
  t->SetBranchAddress("trueCC",&trueCC);
  t->SetBranchAddress("trueNC",&trueNC);
  t->SetBranchAddress("trueQEL",&trueQEL);
  t->SetBranchAddress("trueRES",&trueRES);
  t->SetBranchAddress("trueDIS",&trueDIS);
  t->SetBranchAddress("trueCOH",&trueCOH);
  t->SetBranchAddress("trueMEC",&trueMEC);
  t->SetBranchAddress("truePi0",&hasPi0);
  t->SetBranchAddress("truePiPlus",&hasPiP);
  t->SetBranchAddress("truePiMinus",&hasPiM);
  t->SetBranchAddress("truePiPlusCher",&hasPiPC);
  t->SetBranchAddress("truePiMinusCher",&hasPiMC);
  t->SetBranchAddress("trueKPlus",&hasKP);
  t->SetBranchAddress("trueKMinus",&hasKM);
  t->SetBranchAddress("trueKPlusCher",&hasKPC);
  t->SetBranchAddress("trueKMinusCher",&hasKMC);
  t->SetBranchAddress("simpleRecoFlag",&simpleflag);
  t->SetBranchAddress("simpleRecoEnergy",&simpleenergy);
  t->SetBranchAddress("simpleRecoVtxX",&simplevtxx);
  t->SetBranchAddress("simpleRecoVtxY",&simplevtxy);
  t->SetBranchAddress("simpleRecoVtxZ",&simplevtxz);
  t->SetBranchAddress("promptMuonTotalPE",&PE);
  t->SetBranchAddress("trueMuonEnergy",&mcmuonE);
  t->SetBranchAddress("trueAngle",&mcangle);
  //t->SetBranchAddress("Qij",&Qij);
  //t->SetBranchAddress("trueFollowerParentPDG",&mcFolPPDG);

  // Weights
  std::vector<double>* All0_weight = new std::vector<double>();
  std::vector<double>* All1_weight = new std::vector<double>();
  std::vector<double>* All2_weight = new std::vector<double>();
  std::vector<double>* All3_weight = new std::vector<double>();
  std::vector<double>* All4_weight = new std::vector<double>();
  std::vector<double>* flux_horncurrent = new std::vector<double>();
  std::vector<double>* flux_expskin = new std::vector<double>();
  std::vector<double>* flux_piplus = new std::vector<double>();
  std::vector<double>* flux_piminus = new std::vector<double>();
  std::vector<double>* flux_kplus = new std::vector<double>();
  std::vector<double>* flux_kminus = new std::vector<double>();
  std::vector<double>* flux_kzero = new std::vector<double>();
  std::vector<double>* flux_pionine = new std::vector<double>();
  std::vector<double>* flux_pionqe = new std::vector<double>();
  std::vector<double>* flux_piontot = new std::vector<double>();
  std::vector<double>* flux_nucine = new std::vector<double>();
  std::vector<double>* flux_nucqe = new std::vector<double>();
  std::vector<double>* flux_nuctot = new std::vector<double>();

  t->SetBranchAddress("weight_All0_UBGenie",&All0_weight);
  t->SetBranchAddress("weight_All1_UBGenie",&All1_weight);
  t->SetBranchAddress("weight_All2_UBGenie",&All2_weight);
  t->SetBranchAddress("weight_All3_UBGenie",&All3_weight);
  t->SetBranchAddress("weight_All4_UBGenie",&All4_weight);
  t->SetBranchAddress("weight_horncurrent_FluxUnisim",&flux_horncurrent);
  t->SetBranchAddress("weight_expskin_FluxUnisim",&flux_expskin);
  t->SetBranchAddress("weight_pioninexsec_FluxUnisim",&flux_pionine);
  t->SetBranchAddress("weight_pionqexsec_FluxUnisim",&flux_pionqe);
  t->SetBranchAddress("weight_piontotxsec_FluxUnisim",&flux_piontot);
  t->SetBranchAddress("weight_nucleoninexsec_FluxUnisim",&flux_nucine);
  t->SetBranchAddress("weight_nucleonqexsec_FluxUnisim",&flux_nucqe);
  t->SetBranchAddress("weight_nucleontotxsec_FluxUnisim",&flux_nuctot);
  t->SetBranchAddress("weight_piplus_PrimaryHadronSWCentralSplineVariation",&flux_piplus);
  t->SetBranchAddress("weight_piminus_PrimaryHadronSWCentralSplineVariation",&flux_piminus);
  t->SetBranchAddress("weight_kminus_PrimaryHadronNormalization",&flux_kminus);
  t->SetBranchAddress("weight_kzero_PrimaryHadronSanfordWang",&flux_kzero);
  t->SetBranchAddress("weight_kplus_PrimaryHadronFeynmanScaling",&flux_kplus);

  // Add branches to output
  TFile fo(outfile.Data(), "recreate");
  TTree* to = t->CloneTree(0);
  to->SetAutoSave(0);

  int category;
  bool mc_ccinc_signal;
  bool mcfv;
  bool recofv, recoMRDInc, recoMRD0pi, reco0pi;
  double recop, recopc, recoPE;
  double mcp, mcct;
  bool mc_no_mesons, no_followers;
  double mrd_eff, dirt_muon;
  std::vector<double>* All_weight = new std::vector<double>();
  std::vector<double>* flux_All = new std::vector<double>();

  TBranch* b_category = to->Branch("category", &category);
  TBranch* b_ccinc_signal = to->Branch("mc_ccinc_signal", &mc_ccinc_signal);
  TBranch* TFV = to->Branch("trueFV", &mcfv);
  TBranch* RFV = to->Branch("recoFV", &recofv);
  TBranch* RPE = to->Branch("recoPE", &recoPE);
  TBranch* TMuP = to->Branch("trueMuonMomentum", &mcp);
  TBranch* TCT = to->Branch("trueCosTheta", &mcct);
  TBranch* RMuP = to->Branch("simpleRecoMomentum", &recop);
  TBranch* RMuPC = to->Branch("simpleRecoMomentumCor", &recopc);
  TBranch* TnoPi = to->Branch("true_no_mesons", &mc_no_mesons);
  TBranch* TnoF = to->Branch("true_no_followers", &no_followers);
  TBranch* RinMRDInc = to->Branch("recoInc_contained_in_MRD", &recoMRDInc);
  TBranch* RinMRD0pi = to->Branch("reco0pi_contained_in_MRD", &recoMRD0pi);
  TBranch* R0pi = to->Branch("reco_0pi", &reco0pi);
  TBranch* MRDEff = to->Branch("MRDEff", &mrd_eff);
  TBranch* DirtMu = to->Branch("DirtMu", &dirt_muon);
  TBranch* WAll = to->Branch("weight_All_UBGenie", &All_weight);
  TBranch* WfAll = to->Branch("weight_flux_all", &flux_All);

  // MRD calibration inputs. TODO: Use a TH1D.
  const char* data_dir = getenv("STV_DATA_DIR");
  std::string mrd_cal_file = std::string(data_dir) + "MRDEffUncYlarge.txt";

  std::vector<double> bins_front, factorY;
  std::string bins_front_str, factorY_str;

  std::ifstream cal_file(mrd_cal_file.c_str(), std::ios::in);
  cal_file >> bins_front_str;
  cal_file >> factorY_str;
  cal_file.close();

  breakCSV(bins_front_str, bins_front);
  breakCSV(factorY_str, factorY);

  // Event loop
  for (Long64_t i=0; i<t->GetEntries(); i++) {
    t->GetEntry(i);

    // To save space, only keep reconstructed events
    if (simpleflag != 1 || std::isnan(simpleenergy)) {
      continue;
    }

    mcfv = FidVol(nuvtxx, nuvtxy, nuvtxz);
    mcp = std::sqrt(mcmuonE*mcmuonE - muon_m*muon_m);
    mcct = std::cos(mcangle*M_PI/180.);
    bool hasVisPi = (hasPi0) || (hasPiPC) || (hasPiMC);
    mc_no_mesons = !((hasVisPi) || (hasKPC) || (hasKMC));
    //no_followers = !(
    //  std::find(mcFolPPDG->begin(), mcFolPPDG->end(), 211) != mcFolPPDG->end() ||
    //  std::find(mcFolPPDG->begin(), mcFolPPDG->end(), -211) != mcFolPPDG->end());

    recofv = FidVol(simplevtxx*100, simplevtxy*100, simplevtxz*100);
    recoPE = (PE > 500) && (PE < 3000);
    double recototE = simpleenergy + muon_m;
    recop  = std::sqrt(recototE*recototE - muon_m*muon_m);
    recopc = std::sqrt(recototE*recototE - muon_m*muon_m)*0.82 + 160.;
    recoMRD0pi = numMRDTracks == 1 ? MRDStop->at(0) : false;
    recoMRDInc = numMRDTracks > 0 ? MRDStop->at(0) : false;
    //reco0pi = ((PE > 200.*Qij*Qij) && (PE < 2000.*std::cbrt(4.5-Qij)+1500.));

    // Event categories
    category = kUnknown;
    int abs_nu_pdg = std::abs(trueNuPDG);
    bool is_mc = (abs_nu_pdg == 12 || abs_nu_pdg == 14 || abs_nu_pdg == 16);

    if (is_mc) {
      if (!mcfv) {
        category = kOOFV;
      }
      else if (trueCC != 1) {
        category = kNC;
      }
      else if (trueNuPDG != 14) {
        if (trueNuPDG == 12 && trueCC == 1) category = kNuECC;
        else category = kOther;
      }
      else {
        mc_ccinc_signal = (
             (trueCC == 1)
          && (trueNuPDG == 14)
          && mcfv
          && (mcct > 0.8)
          && (mcp >= 600) && (mcp < 1200)
        );

        if (mc_ccinc_signal) {
          if (trueQEL == 1) category = kSignalCCQE;
          else if (trueMEC == 1) category = kSignalCCMEC;
          else if (trueRES == 1) category = kSignalCCRES;
          else category = kSignalOther;
        }
        else {
          category = kNuMuCCOther;
        }
      }
    }

    // MRD efficiency
    mrd_eff = 1.0;  // assign 1 for no MRD tracks
    for (size_t j=0; j<MRDTrackStartY->size(); j++) {
      if (MRDStop->at(j)) {
        double MRD_Y = MRDTrackStartY->at(j);
        int binY = findBin(MRD_Y, 0, bins_front);
        mrd_eff *= factorY[binY];
      }
    }

    // Dirt muon correction
    dirt_muon = nuvtxz < 0 ? 0.0697 : 1;

    // Weight vectors
    All_weight->insert(All_weight->end(), All0_weight->begin(), All0_weight->end());
    All_weight->insert(All_weight->end(), All1_weight->begin(), All1_weight->end());
    All_weight->insert(All_weight->end(), All2_weight->begin(), All2_weight->end());
    All_weight->insert(All_weight->end(), All3_weight->begin(), All3_weight->end());
    All_weight->insert(All_weight->end(), All4_weight->begin(), All4_weight->end());

    for (int j=0; j<1000; j++) {
      flux_All->push_back(
        flux_horncurrent->at(j) *
        flux_expskin->at(j) *
        flux_piplus->at(j) *
        flux_piminus->at(j) *
        flux_kplus->at(j) *
        flux_kminus->at(j) *
        flux_kzero->at(j) *
        flux_pionine->at(j) *
        flux_pionqe->at(j) *
        flux_piontot->at(j) *
        flux_nucine->at(j) *
        flux_nucqe->at(j) *
        flux_nuctot->at(j)
      );
    }

    to->Fill();

    All_weight->clear();
    flux_All->clear();
  }

  // Drop weights folded into "all" categories
  std::vector<const char*> weight_branches = {
    "weight_All0_UBGenie",
    "weight_All1_UBGenie",
    "weight_All2_UBGenie",
    "weight_All3_UBGenie",
    "weight_All4_UBGenie",
    "weight_All5_UBGenie",
    "weight_horncurrent_FluxUnisim",
    "weight_expskin_FluxUnisim",
    "weight_pioninexsec_FluxUnisim",
    "weight_pionqexsec_FluxUnisim",
    "weight_piontotxsec_FluxUnisim",
    "weight_nucleoninexsec_FluxUnisim",
    "weight_nucleonqexsec_FluxUnisim",
    "weight_nucleontotxsec_FluxUnisim",
    "weight_piplus_PrimaryHadronSWCentralSplineVariation",
    "weight_piminus_PrimaryHadronSWCentralSplineVariation",
    "weight_kminus_PrimaryHadronNormalization",
    "weight_kzero_PrimaryHadronSanfordWang",
    "weight_kplus_PrimaryHadronFeynmanScaling"
  };

  for (const char* k : weight_branches) {
    to->SetBranchStatus(k, 0);
  }

  TTree* to_prune = to->CloneTree();

  // Write and close
  to_prune->Write();
  fo.Write(0, TObject::kOverwrite);
  fo.Close();
  f_in.Close();
}


int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " INPUT.root [OUTPUT.root]" << std::endl;
    return 0;
  }

  // Input file
  TString infile(argv[1]);

  // Optional output file, defaults to X.root -> X.stv.root
  TString outfile(TString(basename(infile)).ReplaceAll(".root", ".stv.root"));

  if (argc >= 3) {
    outfile = argv[2];
  }

  stvPrep(infile, outfile);
}

