#!/bin/bash

# Submit multiple wcsim parallel job batches
# James Minock (code based on other people's code as usual)

export SCRIPT_PATH=/exp/annie/app/users/jminock/grid

#jobsub_submit -G annie --dag file://${SCRIPT_PATH}/wcsim.dagnabbit

QUEUE=long

INPUT_FOLDER=/pnfs/annie/persistent/users/jminock/input
RUNS_FOLDER=/pnfs/annie/persistent/processed/CC_analysis_v1.0
#OUTPUT_FOLDER=/pnfs/annie/scratch/users/jminock/ntuple
OUTPUT_FOLDER=/pnfs/annie/persistent/users/jminock/v1_3_4_data

if [ ! -d ${OUTPUT_FOLDER} ]
then
    echo "${OUTPUT_FOLDER} does not exist"
    mkdir -p $OUTPUT_FOLDER
fi

redoruns=(1614 1629 1641 1650 1661 1679 1690 1695 1724 1805 1816 1822 1829 1831 1834 1837 1839 1842 1856 1885 1887 1888 1894 1903 1904 1905 1911 1916 1933 1934)

#there are 20000 events per GENIE and ANNIEDirt file
#only 1000 events per wcsim can be run on the grid in a reasonable amount of time
#loop through batches of 20 parallel jobs that use the same inputs files
while IFS= read -r LINE; do
    if [[ ! -d ${RUNS_FOLDER}/${LINE} ]]
    then
        continue
    fi
    if [[ -f ${OUTPUT_FOLDER}/PhaseIITree_0.${LINE}.root ]]
    then
        continue
    fi
    if ! compgen -G "${OUTPUT_FOLDER}/PhaseIITree_0.ProcessedData_PMTMRDLAPPD_${LINE}*.root" > /dev/null
    then
        continue
    fi
    hadd $OUTPUT_FOLDER/PhaseIITree_0.${LINE}.root $OUTPUT_FOLDER/PhaseIITree_0.ProcessedData_PMTMRDLAPPD_${LINE}*
    rm $OUTPUT_FOLDER/PhaseIITree_0.ProcessedData_PMTMRDLAPPD_${LINE}*
done < runs.txt


