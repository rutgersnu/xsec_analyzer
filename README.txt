Guide for ANNIE Usage of STV-Analysis

--------------
***Workflow***
--------------
1. Enter the ANNIE ToolAnalysis container: `./start_singularity.sh`
2. Set up the environment: `source setup_stv.sh`
4. Clean and compile: `make clean && make`
3. Prepare input files if needed: `bin/annie_stv_prep`
5. Univmake: './univmake files_to_process.txt tutorial_bin_config.txt output.root'
6. Systematic uncertainty: 'root -l tutorial_slice_plots.C'
7. XSec & Event Rate Plots: ./chi_square_cc0pi_christian

--------------
***DVShiftE***
--------------
In place of detector systematics, ANNIE uses a generalized model based on through-going muons and michel electrons to smear reco observables to match the uncertainty of the detector.
-Find way to use multiple universes???
-UniverseMaker.hh lines 803-804 assign model and apply to reco observable.
-Official model does not yet exist, currently being worked on by Luis Mora-Lepin. Temporary model is gaussian function with mean of 0 and sigma of 1.5 added to observable.

------------------------------------------------
***Notable Differences between ANNIE & uBooNE***
------------------------------------------------
-TTree Name - SystematicsCalculator.hh line 854 & univmake.C line 28 & chi_square_cc0pi_christian.cpp line 364 - TTree name for ANNIE files is 'phaseIITriggerTree'
-POT Problem - SystematicsCalculator.hh line 861 - POT is hardcoded
  ANNIE MC ntuples do not include correct POT. POT is calculated externally for original GENIE files. The POT used for ANNIE MC ntuples is approximated based on the number of ntuple files to their corresponding GENIE files. A fix for this would be required in WCSim.
  Ensure POT on line 861 matches POT for all ANNIE MC files.
  POT for ANNIE numuMC, DV, and DVCV must be the same.
  The POT is ultimately an estimate. POT per event currently does not exist with ANNIE.
-Spline weight - UniverseMaker.hh lines 89, 102 - Spline Weight is non-existent and unnecessary for ANNIE so the weight was replaced with 1.
-univmake output - tutorial_slice_plots.C line 33 & chi_square_cc0pi_christian.cpp lines 179,666 - currently, univmake output is saved as 'output.root'
-ANNIE FV & flux - includes/AnnieGeometryTools.* - uses ANNIE FV and integrated flux window

---------------------------
***Most up to date files***
---------------------------
-BG file - (just a regular ntuple with only a single event, will need an actaul BG file in future)
/exp/annie/data/users/jminock/standard_tank_ntuples/PhaseIITree_bg2_ntuple.root
/pnfs/annie/persistent/users/jminock/v1_3_3_stv_ntuples/PhaseIITree_bg_ntuple.root

-MC files - File name, ~# of events, ~POT, been through stvPrep?
/pnfs/annie/persistent/users/jminock/v1_3_3_stv_ntuples/PhaseIITree_40k_stv_ntuple.root    40000  1.508e19 Y
/pnfs/annie/persistent/users/jminock/v1_3_3_stv_ntuples/PhaseIITree_400k_stv_ntuple.root   390000 1.493e20 Y
/pnfs/annie/persistent/users/jminock/v1_3_3_stv_ntuples/PhaseIITree_1mil_stv_ntuple_1.root 966999 3.706e20 Y
/pnfs/annie/persistent/users/jminock/v1_3_3_stv_ntuples/PhaseIITree_1mil_stv_ntuple_2.root 961997 3.694e20 Y
/pnfs/annie/persistent/users/jminock/v1_3_3_stv_ntuples/PhaseIITree_1mil_stv_ntuple_3.root 954998 3.649e20 Y
/pnfs/annie/persistent/users/jminock/v1_3_3_stv_ntuples/PhaseIITree_1mil_stv_ntuple_4.root 959996 3.676e20 Y
/pnfs/annie/persistent/users/jminock/v1_3_3_stv_ntuples/PhaseIITree_400k_BNB_ntuple.root  390000 1.493e20  Y
/pnfs/annie/persistent/users/jminock/v1_3_3_stv_ntuples/PhaseIITree_400k_DVCV_ntuple.root 390000 1.493e20  Y
/pnfs/annie/persistent/users/jminock/v1_3_3_stv_ntuples/PhaseIITree_400k_DV_ntuple.root   390000 1.493e20  Y

-fake data - Ntuples with only CV weight (use until we get real data working)
/pnfs/annie/persistent/users/jminock/fake-data-nuwro/stv_ntuples/ 39742  1.453e19 Y
/pnfs/annie/persistent/users/jminock/fake-data-nuwro/stv_ntuples/ 397478 1.431e20 Y
/pnfs/annie/persistent/users/jminock/fake-data-nuwro/stv_ntuples/ 992843 3.618e20 Y
/exp/annie/data/users/mastbaum/cc0pi/fake_data/fake_data_cv.root       444865 1.7e20 Y
/exp/annie/data/users/mastbaum/cc0pi/fake_data/fake_data_ccqel-up.root 474530 1.7e20 Y
/exp/annie/data/users/mastbaum/cc0pi/fake_data/fake_data_ccqel-dn.root 412453 1.7e20 Y
/exp/annie/data/users/mastbaum/cc0pi/fake_data/fake_data_ccres-up.root 465738 1.7e20 Y
/exp/annie/data/users/mastbaum/cc0pi/fake_data/fake_data_ccres-dn.root 420431 1.7e20 Y
/exp/annie/data/users/mastbaum/cc0pi/fake_data/fake_data_ccmec-up.root 449936 1.7e20 Y
/exp/annie/data/users/mastbaum/cc0pi/fake_data/fake_data_ccmec-dn.root 437033 1.7e20 Y
/exp/annie/data/users/mastbaum/cc0pi/fake_data/fake_data_flux-up.root  514323 1.7e20 Y
/exp/annie/data/users/mastbaum/cc0pi/fake_data/fake_data_flux-dn.root  373855 1.7e20 Y
/exp/annie/data/users/mastbaum/cc0pi/fake_data/fake_data_pmu-up.root   444577 1.7e20 Y
/exp/annie/data/users/mastbaum/cc0pi/fake_data/fake_data_pmu-dn.root   442997 1.7e20 Y

-stv output - Currently ran with corresponding MC as both MC and fake beam data - OUTDATED
/pnfs/annie/persistent/users/jminock/stv-output/stv-40k-output.root
/pnfs/annie/persistent/users/jminock/stv-output/stv-400k-output.root


ATM Refactor 9/2026
===================
1. `bin/annie_stv_prep INPUT.root [OUTPUT.root]` - Add cut branches, prep files
   * Filename defaults `path/INPUT.root` -> `./INPUT.stv.root`
   * Uses MRD efficiency map in `$STV_DATA_DIR` (i.e., `data/`)
   * Added true event category flags (for slice plot stacks)
   * TODO: Add in weight branches currently in other preprocessing scripts
2. `./bin/univmake [FILES_TO_PROCESS] [BIN_CONFIG] [OUTPUT.root] [FILE_PROPS]` - Build universes
   * Files to process: `config/input_files/files_to_process_closure.txt`
   * Bin config: `config/ccinc/bin_config_ccinc.txt`
   * Output: `/exp/annie/data/...` preferred
   * File properties: `config/input_files/file_properties_closure.txt`
   * TODO: Fix hard-coded POT (c. SystematicsCalculator.cc:861)
3. `./bin/slice_plots_ccinc [UNIVMAKE] [SYST_CONFIG] [BIN_CONFIG]` - Plot event rates and fractional uncertanties
   * Now plots background stacks based on event categories
   * TODO: Fix hard-coded POT
4. TODO: `./bin/chi_square_cc0pi_christian` - Unfold and produce cross sections
   * TODO: Review, tidy, simplify
   * TODO: Unify with slice_plots to produce results in one step

Other to-do items:
* Restore the DV shifts in some form
* Check MRD z efficiency
* Final cleanup

