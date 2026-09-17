#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

void cc0pi_slim(TString filename) {

  TChain tc("phaseIITriggerTree");
  tc.Add(filename);
  assert(tc && tc.IsOpen());

  TObjArray* branches = tc.GetListOfBranches();
  for (int i=0; i<branches->GetEntries(); i++) {
    TString name = branches->At(i)->GetName();
    if(name.Contains("trueNuIntxVtx")) continue;
    if(name.Contains("trueAngle")) continue;
    if(name.Contains("trueMuonEnergy")) continue;
    if(name.Contains("trueCC")) continue;
    if(name.Contains("trueFSLPdg")) continue;
    if(name.Contains("trueTrackLengthInMRD")) continue;
    if(name.Contains("truePi0")) continue;
    if(name.Contains("truePiPlusCher")) continue;
    if(name.Contains("truePiMinusCher")) continue;
    if(name.Contains("trueKPlusCher")) continue;
    if(name.Contains("trueKMinusCher")) continue;
    if(name.Contains("simpleRecoVtx")) continue;
    if(name.Contains("simpleRecoCosTheta")) continue;
    if(name.Contains("simpleRecoEnergy")) continue;
    if(name.Contains("simpleRecoFlag")) continue;
    if(name.Contains("MRDStop")) continue;
    if(name.Contains("RCSRPred")) continue;
    if(name.Contains("numMRDTracks")) continue;
    if(name.Contains("trueq0")) continue;
    if(name.Contains("trueq3")) continue;
    std::cout << "DISABLE " << name << std::endl;
    tc.SetBranchStatus(name, 0);
//    if (name.Contains("weight_")) {
//      std::cout << "DISABLE " << name << std::endl;
//      tc.SetBranchStatus(name, 0);
//    }
  }

  //Create a new file + a clone of old tree header. Do not copy events
  TString n = "PhaseIITree_4mil_slim_ntuple.root";
  TFile* newfile = new TFile(n, "recreate");
  TTree* newtree = tc.CloneTree();

  newfile->Write();
  delete newfile;
}
