#!/bin/bash

# Number of expected command-line arguments
num_expected=2

STV_PATH=/pnfs/annie/persistent/users/jminock/v1_3_4_world_stv_noWeight_ntuples

echo "Combining to make parts"
hadd ${STV_PATH}/PhaseIITree_400k_stv_ntuple_0.root ${STV_PATH}/PhaseIITree_0.?.0.root ${STV_PATH}/PhaseIITree_0.??.0.root ${STV_PATH}/PhaseIITree_0.???.0.root
hadd ${STV_PATH}/PhaseIITree_400k_stv_ntuple_1.root ${STV_PATH}/PhaseIITree_0.1???.0.root
hadd ${STV_PATH}/PhaseIITree_400k_stv_ntuple_2.root ${STV_PATH}/PhaseIITree_0.2???.0.root
hadd ${STV_PATH}/PhaseIITree_400k_stv_ntuple_3.root ${STV_PATH}/PhaseIITree_0.3???.0.root
hadd ${STV_PATH}/PhaseIITree_400k_stv_ntuple_4.root ${STV_PATH}/PhaseIITree_0.4???.0.root

echo "Combining All"
hadd ${STV_PATH}/PhaseIITree_jupyter_final_stv_ntuple.root ${STV_PATH}/PhaseIITree_400k_stv_ntuple_0.root ${STV_PATH}/PhaseIITree_400k_stv_ntuple_1.root ${STV_PATH}/PhaseIITree_400k_stv_ntuple_2.root ${STV_PATH}/PhaseIITree_400k_stv_ntuple_3.root ${STV_PATH}/PhaseIITree_400k_stv_ntuple_4.root

rm ${STV_PATH}/PhaseIITree_400k_stv_ntuple_0.root ${STV_PATH}/PhaseIITree_400k_stv_ntuple_1.root ${STV_PATH}/PhaseIITree_400k_stv_ntuple_2.root ${STV_PATH}/PhaseIITree_400k_stv_ntuple_3.root ${STV_PATH}/PhaseIITree_400k_stv_ntuple_4.root

echo "Done! Don't forget the NuWro samples and counting POT!"

