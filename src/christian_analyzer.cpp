/**
 *This Code is to Anaylizer flat trees
 *
 * **/

#include "FlatTreeAnalyzer.h"




void FlatTreeAnalyzer::Loop(std::string inputTreeType) {

	//----------------------------------------//	

	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;

	double Units = 1E38; // so that the extracted cross-section is in 10^{-38} cm^{2}
	double A = 40; // so that we can have xsecs per nucleus Argon
    //double A = 16; //oxygen
    //double A = 18; //h20
    //double A = 0;
    //if(target_type=="Argon")A = 40;
    //else if(target_type=="Oxygen")A = 16;
    std::cout<<"Selected Target = " << A << std::endl;
    
	int NInte = 6; // Interaction processes: All, QE, MEC, RES, DIS, COH
	std::vector<TString> InteractionLabels = {"","_interaction_QE","_interaction_MEC","_interaction_RES","_interaction_DIS","_interaction_COH"};

   if(A==0){std::cout<<"Target failed to be set"<< std::endl; return;}

	//----------------------------------------//	

        // Output file

	TString FileNameAndPath = fOutputFile + ".root";
	TFile* file = new TFile(FileNameAndPath,"recreate");

	std::cout << std::endl << "------------------------------------------------" << std::endl << std::endl;
	std::cout << "File " << FileNameAndPath << " to be created" << std::endl << std::endl;
	
	//----------------------------------------//
    std::vector<double> PmuBinnEdges_wc{0.60, 0.74, 0.86, 1.00, 1.10, 1.2}; // Annie Binning 
    std::vector<double> Costheaedges{.8, 0.95, 1.0};

   std::vector<double> BinN = {1, 2, 3, 4, 5, 6, 7};
   int BinN_wc = 5; // 5 Pmu bins and 1 out of signal 
     std::map< double, std::vector<double> > MUON_2D_BIN_EDGES_WC = {
   // the 2D binning of inclusive but took the max limit of 2 GeV and not 2.5 GeV which its given  
   { .8, {PmuBinnEdges_wc}},
   { .95,{PmuBinnEdges_wc}},
    { 1.0, {} }
   };
   
   std::vector<double> Kp_Edges = {0.25, 0.28, 0.315, 0.35, 0.42, 0.525, 0.6, 0.8, 1.0};
   
   
   
//0 0 "(mc_is_cc0pi_wc_signal && mc_p3_mu.CosTheta() > 0.8 && mc_p3_mu.Mag() >= 0.600 && mc_p3_mu.Mag() < 0.740)"
//0 0 "(mc_is_cc0pi_wc_signal && mc_p3_mu.CosTheta() > 0.8 && mc_p3_mu.Mag() >= 0.740 && mc_p3_mu.Mag() < 0.860)"
//0 0 "(mc_is_cc0pi_wc_signal && mc_p3_mu.CosTheta() > 0.8 && mc_p3_mu.Mag() >= 0.860 && mc_p3_mu.Mag() < 1.000)"
//0 0 "(mc_is_cc0pi_wc_signal && mc_p3_mu.CosTheta() > 0.8 && mc_p3_mu.Mag() >= 1.000 && mc_p3_mu.Mag() < 1.100)"
//0 0 "(mc_is_cc0pi_wc_signal && mc_p3_mu.CosTheta() > 0.8 && mc_p3_mu.Mag() >= 1.100 && mc_p3_mu.Mag() < 1.200)"
//0 0 "(mc_is_cc0pi_wc_signal && (mc_p3_mu.CosTheta() <= 0.8 || mc_p3_mu.Mag() < 0.600 || mc_p3_mu.Mag() >= 1.200))"
   
   
   

   
   std::vector<double> theta_P_Edges = {-1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    int NBinsWC2D = 10;
    int NBinsWC_total = 45;
    
   UBTH2Poly* h_UBTH2Poly_WC = new UBTH2Poly("h_UBTH2Poly_WC", "h_UBTH2Poly_WC",  MUON_2D_BIN_EDGES_WC, false);
   
   TH1D* h_pmu_WC = new TH1D("h_pmu_WC" ,"h_pmu_WC", PmuBinnEdges_wc.size() - 1, PmuBinnEdges_wc.data());
   
   TH1D* h_BinN_WC = new TH1D("h_BinN_WC" ,"h_BinN_WC", BinN_wc, 1.0, BinN_wc+1);
   
   
   TH1D* h_costheta_WC = new TH1D("h_costheta_WC" ,"h_costheta_WC", Costheaedges.size() - 1, Costheaedges.data());
   TH1D* h_p_p_WC = new TH1D("h_p_p_WC" ,"h_p_p_WC", Kp_Edges.size() - 1, Kp_Edges.data());
   TH1D* h_costhetap_WC = new TH1D("h_costhetap_WC" ,"h_costhetap_WC", theta_P_Edges.size() - 1, theta_P_Edges.data());
   TH1D* h_2D_BinN_WC = new TH1D("h_2D_BinN_WC" ,"h_2D_BinN_WC", NBinsWC2D,1.,NBinsWC2D+1);
   
   TH1D* h_2D_BinN_WC_Total = new TH1D("h_2D_BinN_WC_total" ,"h_2D_BinN_WC_total", NBinsWC_total,1.,NBinsWC_total+1);
   TH1D* TrueMuon_PmuCosTheta_WCSection_binN[NInte];
   

	//TH1D* TrueMuonCosThetaPlot[NInte];
	
    TH1D* h_pmu_WC_interaction[NInte];
    TH1D* h_BinN_WC_interaction[NInte];
    
    TH1D* h_costheta_WC_interaction[NInte];
    TH1D* h_p_p_WC_interaction[NInte];
    TH1D* h_costhetap_WC_interaction[NInte];
	
	
  
	// Loop over the interaction processes
    // Better to fill by Bin Number 
    // Binning starts at 1 

	// Loop over the interaction processes


	for (int inte = 0; inte < NInte; inte++) {

	  //--------------------------------------------------//
	  //TrueMuonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueMuonCosThetaPlot",";cos(#theta_{#mu})",10,-1.,1.);
       h_BinN_WC_interaction[inte]= new TH1D("h_BinN_WC_interaction"+ InteractionLabels[inte] ,"h_BinN_WC_interaction"+ InteractionLabels[inte], BinN_wc,1.0, BinN_wc+1);
	  
	   //h_pmu_WC_interaction[inte]= new TH1D("h_pmu_WC"+ InteractionLabels[inte] ,"h_pmu_WC"+ InteractionLabels[inte], PmuBinnEdges_wc.size() - 1, PmuBinnEdges_wc.data());
       //h_costheta_WC_interaction[inte]= new TH1D("h_costheta_WC"+ InteractionLabels[inte] ,"h_costheta_WC"+ InteractionLabels[inte], Costheaedges.size() - 1, Costheaedges.data());
       //h_p_p_WC_interaction[inte]= new TH1D("h_p_p_WC"+ InteractionLabels[inte] ,"h_p_p_WC"+ InteractionLabels[inte], Kp_Edges.size() - 1, Kp_Edges.data());
       //h_costhetap_WC_interaction[inte]= new TH1D("h_costhetap_WC"+ InteractionLabels[inte] ,"h_costhetap_WC"+ InteractionLabels[inte], theta_P_Edges.size() - 1, theta_P_Edges.data());
	  
	  //--------------------------------------------------//:

	} // End of the loop over the interaction processes							

	//----------------------------------------//

	// Counters

	int CounterEventsPassedSelection = 0;
	int CounterQEEventsPassedSelection = 0;
	int CounterMECEventsPassedSelection = 0;
	int CounterRESEventsPassedSelection = 0;
	int CounterDISEventsPassedSelection = 0;
	int CounterCOHEventsPassedSelection = 0;	

	//----------------------------------------//
	
	// Loop over the events

	for (Long64_t jentry=0; jentry<nentries;jentry++) {

	  //----------------------------------------//	

	
	  Long64_t ientry = LoadTree(jentry);
	  if (ientry < 0) break; nb = fChain->GetEntry(jentry); nbytes += nb;
	  if (jentry%100000 == 0) std::cout << jentry/100000 << " 100 k " << std::setprecision(3) << double(jentry)/nentries*100. << " %"<< std::endl;

	  //----------------------------------------//	
		
	  //double weight = fScaleFactor*Units*A*Weight; //*Weight Gibuu wg
	  //weight=weight/250.0;
       //if(inputTreeType=="GiBUU_2023"){weight=weight/81.0;}
       //weight=weight/81.0;
       //double weight = fScaleFactor*Units*A; //*Weight
	  //weight=weight/150.0;
      
       double weight = fScaleFactor*Units*A; //*Weight
       
       	  //std::cout<<"inputTreeType="<< inputTreeType<< std::endl;
       
         if(inputTreeType=="GiBUU_2025_Argon"){
	  if (jentry%100000 == 0) std::cout<<"using GIBUU_2025"<< fScaleFactor<<std::endl;
	  weight = fScaleFactor*Units*A*Weight;
	  weight=weight/250.0;
	  }
      else if(inputTreeType=="GiBUU_2025_Oxygen"){
	  if (jentry%100000 == 0) std::cout<<"using GIBUU_2025"<< fScaleFactor<<std::endl;
	  weight = fScaleFactor*Units*A*Weight;
	  weight=weight/250.0;
	  }
	  
	  
	  
       
       
     // double weight =fScaleFactor*Units*A  /  POT_inter;

     //double weight =fScaleFactor * ( A ) *  1e38 / POT_inter ;
	  double weight_noUnits = fScaleFactor;// *Weight
	  //if(inputTreeType=="GiBUU_2023"){weight_noUnits=weight_noUnits/81.0;}
	  //weight_noUnits=weight_noUnits/81.0;
	  weight_noUnits=weight_noUnits;
	  //----------------------------------------//	
      //std::cout<<"fScaleFactor = "<< fScaleFactor << std::endl;
	  // Signal definition

	  if (PDGLep != 13) { continue; } // make sure that we have only a muon in the final state
	  if (cc != 1) { continue; } // make sure that we have only CC interactions		

	  int ProtonTagging = 0, ChargedPionTagging = 0, NeutralPionTagging = 0, MuonTagging = 0;
          int ElectronTagging = 0, PhotonTagging = 0;
	  //std::vector <int> ProtonID; ProtonID.clear();
	  std::vector <int> MuonID; MuonID.clear();		

	  // Example selection with CC0pi (units in GeV/c)
	  // Loop over final state particles
		double Pmu =0; 
		int indexProtonGreatestP= -99; 
		double Proton_greatest = .25; // has to greater than .25 to be a leading protoon 
		double p_angle = -99;
		bool GoodleadP =false; 
		
	  for (int i = 0; i < nfsp; i++) {
		
	    double pf = TMath::Sqrt( px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i]);

	    if (pdg[i] == 13) {
         
         // const float MUON_P_MIN_WC_MOM_CUT_jointcc0pi =.121; // GeV/c
         // adding ANNIE accpence 
	      MuonTagging ++;
	      MuonID.push_back(i);
          Pmu = pf;
	    }
   

	    if (fabs(pdg[i]) == 211 && pf > 0.160)  {
         //const float CHARGED_PI_WC_MOM_CUT_jointcc0pi =.161; // GeV/c increase the pion threhold for WC 

	      ChargedPionTagging ++;

	    }

	    if (pdg[i] == 111)  {

	      NeutralPionTagging ++;

	    }

	    if (fabs(pdg[i]) == 11)  {

	      ElectronTagging ++;

	    }

	    if (fabs(pdg[i]) == 22)  {

	      PhotonTagging ++;

	    }


         if (fabs(pdg[i]) == 2212)  {

        if(pf > Proton_greatest && pf < 1.0 ){
        GoodleadP=true;
        Proton_greatest = pf;
        indexProtonGreatestP = indexProtonGreatestP;
        p_angle = pz[i] / pf;
        }

	    }


	  } // End of the loop over the final state particles

	  // If the signal definition is not satisfied, continue

	  if (
	   ChargedPionTagging != 0 || 
		NeutralPionTagging != 0 || 
		MuonTagging !=1 )
		{ continue; }

	  //----------------------------------------//	
	  //-------apply Opimz ANNIE phase space---------------------------------//	
         if(CosLep <= 0.6 || Pmu < 0.600 || Pmu >= 1.200) continue;
	  // https://arxiv.org/pdf/2106.15809.pdf


	  CounterEventsPassedSelection++;
	
	  // Classify the events based on the interaction type

	  int genie_mode = -1.;
	  if (TMath::Abs(Mode) == 1) { CounterQEEventsPassedSelection++; genie_mode = 1; } // QE
	  else if (TMath::Abs(Mode) == 2) { CounterMECEventsPassedSelection++; genie_mode = 2; } // MEC
	  else if (
		   TMath::Abs(Mode) == 10 ||
		   TMath::Abs(Mode) == 11 || TMath::Abs(Mode) == 12 || TMath::Abs(Mode) == 13 ||
		   TMath::Abs(Mode) == 17 || TMath::Abs(Mode) == 22 || TMath::Abs(Mode) == 23
		   ) { CounterRESEventsPassedSelection++; genie_mode = 3; } // RES
	  else if (TMath::Abs(Mode) == 21 || TMath::Abs(Mode) == 26) { CounterDISEventsPassedSelection++; genie_mode = 4; } // DIS
	  else if (TMath::Abs(Mode) == 16) { CounterCOHEventsPassedSelection++; genie_mode = 5;} // COH
	  else { continue; }  

	  // Feb 8 2022: Only case that is not covered is 15 = diffractive

	  //----------------------------------------//

	  // filling in the histo regardless of interaction mode

	  //TrueMuonCosThetaPlot[0]->Fill(CosLep,weight);

     int inputBinnN = 99;
      
      /*
       // IF in overlow bin else not 
      if(CosLep <= 0.8 || Pmu < 0.600 || Pmu >= 1.200){
      
      inputBinnN = 6;
      
      }
      
      else{
      inputBinnN = h_pmu_WC->Fill(Pmu,weight);
      
      }
      */
    
     inputBinnN = h_pmu_WC->Fill(Pmu,weight);
      //std::cout<<"inputBinnN = "<< inputBinnN << " CosLep = "<<CosLep<< " Pmu = "<< Pmu<<  std::endl;

      h_BinN_WC->Fill(inputBinnN,weight); 
      h_BinN_WC_interaction[genie_mode]->Fill(inputBinnN,weight); 
	  //----------------------------------------//

	  // filling in the histo based on the interaction mode

	  //TrueMuonCosThetaPlot[genie_mode]->Fill(CosLep,weight);

      int WC2D_binN =h_UBTH2Poly_WC->Fill(CosLep,Pmu,weight);
      
      
      h_costheta_WC->Fill(CosLep,weight);
      h_2D_BinN_WC->Fill(WC2D_binN,weight);
      
      if(GoodleadP==true){
       h_p_p_WC->Fill(Proton_greatest,weight);
       h_costhetap_WC->Fill(p_angle,weight);
      }
     

      //std::cout<<"binNscheme1 = "<< binNscheme1<< std::endl;


      
	  //----------------------------------------//

	  // filling in the histo based on the interaction mode	  
	    
	  //----------------------------------------//
	
	} // End of the loop over the events

	//----------------------------------------//	







	//std::cout << "Percetage of events passing the selection cuts = " << 
	//double(CounterEventsPassedSelection)/ double(nentries)*100. << " %" << std::endl; std::cout << std::endl;
//
	//std::cout << "Success percetage in selecting QE events = " << 
	//double(CounterQEEventsPassedSelection)/ double(CounterEventsPassedSelection)*100. << " %" << std::endl; std::cout << std::endl;
//
	//std::cout << "Success percetage in selecting MEC events = " << 
	//double(CounterMECEventsPassedSelection)/ double(CounterEventsPassedSelection)*100. << " %" << std::endl; std::cout << std::endl;
//
	//std::cout << "Success percetage in selecting RES events = " << 
	//double(CounterRESEventsPassedSelection)/ double(CounterEventsPassedSelection)*100. << " %" << std::endl; std::cout << std::endl;
//
	//std::cout << "Success percetage in selecting DIS events = " << 
	//double(CounterDISEventsPassedSelection)/ double(CounterEventsPassedSelection)*100. << " %" << std::endl; std::cout << std::endl;
//
	//std::cout << "Success percetage in selecting COH events = " << 
	//double(CounterCOHEventsPassedSelection)/ double(CounterEventsPassedSelection)*100. << " %" << std::endl; std::cout << std::endl;	
//
	//----------------------------------------//	
	//----------------------------------------//	

	// Division by bin width to get the cross sections	
	// Loop over the interaction processes
std::cout<<"BinWidth Normalizing"<< std::endl;
	/*
	for (int inte = 0; inte < NInte; inte++) {

		//----------------------------------------//
	
		//Reweight(TrueMuonCosThetaPlot[inte]);
        
        
        
		//----------------------------------------//
    

	} // End of the loop over the interaction processes		
*/

//h_pmu_WC->Scale(1,"width");
//h_costheta_WC->Scale(1,"width");
//h_2D_BinN_WC->Scale(1,"width");
//h_p_p_WC->Scale(1,"width");
//h_costhetap_WC->Scale(1,"width");

        std::cout<<"Full 1D Cross Section"<<std::endl;
		
		std::cout<<"Writing Files "<< std::endl; 
		
		
		
	file->cd();
	file->Write();
	fFile->Close();

char Histtitle[1024];
   TCanvas *can = new TCanvas(uniq());
    gStyle->SetOptStat(0); 
    int precision = 2;  // Change this value based on your requirement
	gStyle->SetPaintTextFormat(Form(".%df", precision));
	float textSize = 0.01;  // Adjust this value based on your requirement
	gStyle->SetTextSize(textSize);
    char text_title_pdf1[2024];
  sprintf(Histtitle, "Model: %s binning Scheme 1",fPDF_Name.c_str());
h_UBTH2Poly_WC->SetTitle(Histtitle); 
h_UBTH2Poly_WC->GetXaxis()->SetTitle("Cos(#theta_{#mu})");
h_UBTH2Poly_WC->GetYaxis()->SetTitle("P_{#mu} [Gev/c]");
h_UBTH2Poly_WC->GetZaxis()->SetLabelSize (0.017);
h_UBTH2Poly_WC->Draw("colz text");
 sprintf(text_title_pdf1, "FlatTreePlots_WC.pdf");
  can -> Print(text_title_pdf1);
  
 h_UBTH2Poly_WC->Scale(1, "width");
 h_UBTH2Poly_WC->GetZaxis()->SetTitle("#frac{d#sigma}{dp_{#mu}} [ 10^{-38} #frac{cm^{2}}{GeV/c Ar} ]");
 h_UBTH2Poly_WC->Draw("colz text");
 can -> Print(text_title_pdf1);
  can->Close();   
	//delete can; 
	std::cout << std::endl;
	std::cout << "File " << FileNameAndPath +" has been created created " << std::endl; 
	std::cout << std::endl;
	std::cout << "------------------------------------------------" << std::endl;
	std::cout << "End of Loop " << std::endl;
	//----------------------------------------//		
return; 
} // End of the program
//////////////////////////////////////////////////////////
//	
//////////////////////////////////////////////////////////
void Reweight(TH1D* h) {

  int NBins = h->GetXaxis()->GetNbins();

  for (int i = 0; i < NBins; i++) {

    double CurrentEntry = h->GetBinContent(i+1);
    double NewEntry = CurrentEntry / h->GetBinWidth(i+1);

    double CurrentError = h->GetBinError(i+1);
    double NewError = CurrentError / h->GetBinWidth(i+1);

    h->SetBinContent(i+1,NewEntry); 
    h->SetBinError(i+1,NewError); 
    //h->SetBinError(i+1,0.000001); 

  }

}
//////////////////////////////////////////////////////////
//	
//////////////////////////////////////////////////////////
void run_FlatTreeAnaylizer(std::string rootfileName_input, std::string outputName, std::string PDF_name, std::string target_type){



std::cout<<"Running FlatTreeAnaylizer For: "<<PDF_name.c_str() << std::endl;

TString in(rootfileName_input);
TString out(outputName);



FlatTreeAnalyzer FlatTreeAnalyzer_object(in, out, PDF_name );

FlatTreeAnalyzer_object.Loop(PDF_name);


}


//----------------------------------------//		
//////////////////////////////////////////////////////////
//	
//////////////////////////////////////////////////////////
void Reweight_2DBinning(TH1D &h, UBTH2Poly* UBTH2Poly_input) {

  int NBins = h.GetXaxis()->GetNbins();

  for (int i = 0; i < NBins; i++) {

    double CurrentEntry = h.GetBinContent(i+1);
    double NewEntry = CurrentEntry / UBTH2Poly_input->GetBinWidth(i+1);

    double CurrentError = h.GetBinError(i+1);
    double NewError = CurrentError / UBTH2Poly_input->GetBinWidth(i+1);

    h.SetBinContent(i+1,NewEntry); 
    h.SetBinError(i+1,NewError); 
    //h->SetBinError(i+1,0.000001); 

  }

}
//////////////////////////////////////////////////////////
//	
//////////////////////////////////////////////////////////

int main() {
std::vector<std::string> WhichSample;
std::vector<std::string> WhichName;
std::vector<std::string> TargetType;



   UBTH2Poly instance;
   gInterpreter->Declare("#include \"/exp/uboone/app/users/cnguyen/stv-analysis-new/includes/UBTH2Poly.h\""); 

std::vector<std::string> input{"NuWro.flat","GENIE_v3_0_6_G18_10a_02a.flat","GiBUU.flat","NEUT.flat"};


//WhichSample.push_back("/exp/uboone/data/users/cnguyen/CC0Pi_Selection/modelcomparson/gibuu_2025/GiBUU_2025.flat_250files"); WhichName.push_back("GiBUU_2025_Argon");TargetType.push_back("Argon");
//WhichSample.push_back("/exp/uboone/data/users/cnguyen/CC0Pi_Selection/modelcomparson/G21_11a_00_000/14_1000180400_CC_v3_4_0_G21_11a_00_000.flat");WhichName.push_back("GENIE_G21_11a_00_000_Argon");TargetType.push_back("Argon");
//WhichSample.push_back("/exp/uboone/data/users/cnguyen/CC0Pi_Selection/modelcomparson/G21_11b_00_000/14_1000180400_CC_v3_4_0_G21_11b_00_000.flat");WhichName.push_back("GENIE_G21_11b_00_000_Argon");TargetType.push_back("Argon");
//WhichSample.push_back("/exp/uboone/data/users/cnguyen/CC0Pi_Selection/modelcomparson/AR23_20i_00_000/14_1000180400_CC_v3_4_0_AR23_20i_00_000.flat");WhichName.push_back("AR23_20i_00_000_Argon");TargetType.push_back("Argon");
//WhichSample.push_back("/exp/uboone/data/users/cnguyen/CC0Pi_Selection/modelcomparson/G18_10a_02_11a/14_1000180400_CC_v3_4_0_G18_10a_02_11a.flat");WhichName.push_back("G18_10a_02_11a_Argon");TargetType.push_back("Argon");
//WhichSample.push_back("/exp/uboone/data/users/cnguyen/CC0Pi_Selection/modelcomparson/G18_10b_02_11a/14_1000180400_CC_v3_4_0_G18_10b_02_11a.flat");WhichName.push_back("G18_10b_02_11a_Argon");TargetType.push_back("Argon");
//WhichSample.push_back("/exp/uboone/data/users/cnguyen/CC0Pi_Selection/modelcomparson/G00_00b_00_000/14_1000180400_CC_v3_4_0_G00_00b_00_000.flat");WhichName.push_back("G00_00b_00_000_Argon");TargetType.push_back("Argon");
//WhichSample.push_back("/exp/uboone/data/users/cnguyen/CC0Pi_Selection/modelcomparson/G18_10i_02_11b/14_1000180400_CC_v3_4_0_G18_10i_02_11b.flat");WhichName.push_back("G18_10i_02_11a_Argon");TargetType.push_back("Argon");
//WhichSample.push_back("/exp/uboone/data/users/cnguyen/CC0Pi_Selection/modelcomparson/G18_10i_02_11b/14_1000180400_CC_v3_4_0_G18_10i_02_11b.flat");WhichName.push_back("G18_10i_02_11a_Argon");TargetType.push_back("Argon");
//WhichSample.push_back("/exp/uboone/data/users/cnguyen/CC0Pi_Selection/modelcomparson/G18_10j_02_11b/14_1000180400_CC_v3_4_0_G18_10j_02_11b.flat");WhichName.push_back("G18_10j_02_11a_Argon");TargetType.push_back("Argon");
//WhichSample.push_back("/exp/uboone/data/users/cnguyen/CC0Pi_Selection/modelcomparson/G18_02a_02_11a/14_1000180400_CC_v3_4_0_G18_02a_02_11a.flat");WhichName.push_back("G18_02a_02_11a_Argon");TargetType.push_back("Argon");
//WhichSample.push_back("/exp/uboone/data/users/cnguyen/CC0Pi_Selection/modelcomparson/G18_01a_02_11b/14_1000180400_CC_v3_4_0_G18_01a_02_11b.flat");WhichName.push_back("G18_01a_02_11a_Argon");TargetType.push_back("Argon");
WhichSample.push_back(" /pnfs/uboone/persistent/users/mastbaum/tuning2022/mc/bnb_ub/flat/bnb.ub.num.nuwro_19_02_1.flat");WhichName.push_back("nuwro_19_02_1");TargetType.push_back("Argon");

//WhichSample.push_back("/exp/uboone/data/users/cnguyen/CC0Pi_Selection/modelcomparson_Oxygen/GiBUU_2025/GiBUU_2025_Oxygen_250files"); WhichName.push_back("GiBUU_2025_Oxygen");TargetType.push_back("Oxygen");
//WhichSample.push_back("/exp/uboone/data/users/cnguyen/CC0Pi_Selection/modelcomparson_Oxygen/G21_11a_00_000/14_1000080160_CC_v3_4_0_G21_11a_00_000.flat");WhichName.push_back("GENIE_G21_11a_00_000_Oxygen");TargetType.push_back("Oxygen");
//WhichSample.push_back("/exp/uboone/data/users/cnguyen/CC0Pi_Selection/modelcomparson_Oxygen/G21_11b_00_000/14_1000080160_CC_v3_4_0_G21_11b_00_000.flat");WhichName.push_back("GENIE_G21_11b_00_000_Oxygen");TargetType.push_back("Oxygen");
//WhichSample.push_back("/exp/uboone/data/users/cnguyen/CC0Pi_Selection/modelcomparson_Oxygen/AR23_20i_00_000/14_1000080160_CC_v3_4_0_AR23_20i_00_000.flat");WhichName.push_back("AR23_20i_00_000_Oxygen");TargetType.push_back("Oxygen");
//WhichSample.push_back("/exp/uboone/data/users/cnguyen/CC0Pi_Selection/modelcomparson_Oxygen/G18_10a_02_11a/14_1000080160_CC_v3_4_0_G18_10a_02_11a.flat");WhichName.push_back("G18_10a_02_11a_Oxygen");TargetType.push_back("Oxygen"); // Old FLux
//WhichSample.push_back("/exp/uboone/data/users/cnguyen/CC0Pi_Selection/modelcomparson_Oxygen_newflux/G18_10a_02_11a/14_1000080160_CC_v3_4_0_G18_10a_02_11a.flat");WhichName.push_back("G18_10a_02_11a_Oxygen");TargetType.push_back("Oxygen");

//WhichSample.push_back("/exp/uboone/data/users/cnguyen/CC0Pi_Selection/modelcomparson_Oxygen/G18_10b_02_11a/14_1000080160_CC_v3_4_0_G18_10b_02_11a.flat");WhichName.push_back("G18_10b_02_11a_Oxygen");TargetType.push_back("Oxygen");
//WhichSample.push_back("/exp/uboone/data/users/cnguyen/CC0Pi_Selection/modelcomparson_Oxygen/G00_00b_00_000/14_1000080160_CC_v3_4_0_G00_00b_00_000.flat");WhichName.push_back("G00_00b_00_000_Oxygen");TargetType.push_back("Oxygen");
//WhichSample.push_back("/exp/uboone/data/users/cnguyen/CC0Pi_Selection/modelcomparson_Oxygen/G18_10i_02_11a/14_1000080160_CC_v3_4_0_G18_10i_02_11a.flat");WhichName.push_back("G18_10i_02_11a_Oxygen");TargetType.push_back("Oxygen");
//WhichSample.push_back("/exp/uboone/data/users/cnguyen/CC0Pi_Selection/modelcomparson_Oxygen/G18_10j_02_11a/14_1000080160_CC_v3_4_0_G18_10j_02_11a.flat");WhichName.push_back("G18_10j_02_11a_Oxygen");TargetType.push_back("Oxygen");
//WhichSample.push_back("/exp/uboone/data/users/cnguyen/CC0Pi_Selection/modelcomparson_Oxygen/G18_02a_02_11a/14_1000080160_CC_v3_4_0_G18_02a_02_11a.flat");WhichName.push_back("G18_02a_02_11a_Oxygen");TargetType.push_back("Oxygen");
//WhichSample.push_back("/exp/uboone/data/users/cnguyen/CC0Pi_Selection/modelcomparson_Oxygen/G18_01a_02_11b/14_1000080160_CC_v3_4_0_G18_01a_02_11b.flat");WhichName.push_back("G18_01a_02_11a_Oxygen");TargetType.push_back("Oxygen");
//

//WhichSample.push_back("/exp/uboone/data/users/cnguyen/CC0Pi_Selection/modelcomparson_Oxygen/nuwro_25.11/NuWro.flat_H2O_annie_numu");WhichName.push_back("nuwro_25.11_H2O");TargetType.push_back("Oxygen");  
//WhichSample.push_back("/exp/uboone/data/users/cnguyen/CC0Pi_Selection/modelcomparson_Oxygen/nuwro_25.11/NuWro.flat_Oxg_annie_numu");WhichName.push_back("nuwro_25.11_Oxygen");TargetType.push_back("Oxygen");
/*
for(auto rootname:input ){


std::string inputTfile = "mySamples/" + rootname;
run_FlatTreeAnaylizer(inputTfile, name);
}
*/

 TCanvas *can1 = new TCanvas(uniq());
 
  can1 -> Print("FlatTreePlots_WC.pdf(");

for(int i =0 ; i < WhichSample.size(); ++i ){


std::string name = "/exp/uboone/data/users/cnguyen/CC0Pi_Selection/EventSelection_3_13_26_CrossSectionOpimz/CC0Pi_WC_thresholds_Flat_" + WhichName.at(i);

run_FlatTreeAnaylizer(WhichSample.at(i), name, WhichName.at(i), TargetType.at(i) );
std::cout<<"Finished FlatTreeAnaylizer for ::  "<< WhichName.at(i).c_str()<< std::endl;


}

can1 -> Print("FlatTreePlots_WC.pdf)");
can1->Close(); 
delete can1;
  return 0;
}
