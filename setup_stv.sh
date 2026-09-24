#!/bin/bash

# Sets up the local environment for working with the STV analysis scripts

# Sets up local environment in ANNIE ToolAnalysis container
# (not all of these may be necessary, but it works as is)

ulimit -n 4096

export LIBGL_ALWAYS_INDIRECT=1

source /ToolAnalysis/ToolDAQ/root-6.24.06/install/bin/thisroot.sh

export CLHEP_DIR=/ToolAnalysis/ToolDAQ/2.4.0.2/CLHEP_install

export LD_LIBRARY_PATH=/ToolAnalysis/lib:/ToolAnalysis/ToolDAQ/zeromq-4.0.7/lib:/ToolAnalysis/ToolDAQ/boost_1_66_0/install/lib:/ToolAnalysis/ToolDAQ/root/lib:/ToolAnalysis/ToolDAQ/WCSimLib/:/ToolAnalysis/ToolDAQ/MrdTrackLib/src:/ToolAnalysis/ToolDAQ/RATEventLib/lib/:/ToolAnalysis/UserTools/PlotWaveforms:/ToolAnalysis/ToolDAQ/log4cpp/lib:/ToolAnalysis/ToolDAQ/Pythia6Support/v6_424/lib:${CLHEP_DIR}/lib:/ToolAnalysis/ToolDAQ/LHAPDF-6.3.0/install/lib:/ToolAnalysis/ToolDAQ/GENIE-v3-master/lib:/ToolAnalysis/ToolDAQ/Reweight-3_00_04_ub3/lib:$LD_LIBRARY_PATH

export ROOT_INCLUDE_PATH=/ToolAnalysis/ToolDAQ/boost_1_66_0/install/include:/ToolAnalysis/ToolDAQ/WCSimLib/include/:/ToolAnalysis/ToolDAQ/MrdTrackLib/include:/ToolAnalysis/ToolDAQ/RATEventLib/include/:/ToolAnalysis/UserTools/PlotWaveforms:$ROOT_INCLUDE_PATH

export PYTHIA6_DIR=/ToolAnalysis/ToolDAQ/Pythia6Support/v6_424/
export LHAPATH=/ToolAnalysis/ToolDAQ/LHAPDF-6.3.0/install/share/LHAPDF/
export GENIE=/ToolAnalysis/ToolDAQ/GENIE-v3-master
export GENIE_REWEIGHT=/ToolAnalysis/ToolDAQ/Reweight-3_00_04_ub3/

export PATH=/ToolAnalysis/ToolDAQ/LHAPDF-6.3.0/install/bin:$GENIE/bin:$GENIE_REWEIGHT/bin:$PATH

export PATH=/ToolAnalysis/ToolDAQ/fsplit:$PATH
export TF_CPP_MIN_LOG_LEVEL=2

# Finds the directory where this script is located. This method isn't
# foolproof. See https://stackoverflow.com/a/246128/4081973 if you need
# something more robust for edge cases (e.g., you're calling the script using
# symlinks).
THIS_DIRECTORY="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

export STV_ANALYSIS_DIR=${THIS_DIRECTORY}
export PATH=${STV_ANALYSIS_DIR}/bin:${PATH}
export LD_LIBRARY_PATH=${STV_ANALYSIS_DIR}/lib:${LD_LIBRARY_PATH}

export STV_DATA_DIR=${STV_ANALYSIS_DIR}/data

