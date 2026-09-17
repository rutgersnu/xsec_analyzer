# Submit multiple wcsim parallel job batches
# James Minock (code based on other people's code as usual)

export SCRIPT_PATH=/exp/annie/app/users/jminock/grid

#jobsub_submit -G annie --dag file://${SCRIPT_PATH}/wcsim.dagnabbit

QUEUE=long

#for RUN in {0..9} #skip 4661
#do
#	echo "Running ${RUN}..."
#	./univmake f${RUN}.txt tutorial_binInc_config.txt output_inc_${RUN}.root fp_${RUN}.txt
#	./univmake f${RUN}.txt tutorial_bin0pi_config.txt output_0pi_${RUN}.root fp_${RUN}.txt
#done

./univmake files_to_process_nuwro.txt tutorial_binInc_config.txt output_nuwro_inc.root file_properties_nuwro.txt
./univmake files_to_process_nuwro.txt tutorial_bin0pi_config.txt output_nuwro_0pi.root file_properties_nuwro.txt


echo "Finished!"

