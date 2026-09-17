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

//Script that strips all weights (except TunedCentralValue) out of ANNIE MC files
//Used for non numuMC files that do not need weights
void stripNuWro(){
        //Open file and trees
        TFile *f = new TFile("/pnfs/annie/persistent/users/jminock/fake-data-nuwro/stv_ntuples/PhaseIITree_1mil_nuwro_ntuple.root","read");
        gSystem->Load("/exp/annie/app/users/jminock/ToolAnalysis/lib/libDataModel.so");
//      gSystem->Load("/exp/annie/app/users/jminock/ToolAnalysis/lib/libDict.so");
        gInterpreter->GenerateDictionary("map<string,vector<double>>", "map;string;vector");
        TTree *T = (TTree*)f->Get("phaseIITriggerTree");

	T->SetBranchStatus("trueCC", 0);

        TFile *fnew = new TFile("PhaseIITree_1mil_temp_nuwro_ntuple.root","recreate");
	auto Tnew = T->CloneTree();

	fnew->Write();
}
