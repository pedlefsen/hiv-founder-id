### This is for the partitions, for the main results see (first) README.postprocessing.txt.

### First ensure you have done everything in README.preprocessing.txt (this one is not specific to the paritions).

### And of course you need to run the code [TED, can you please make README.processing_partitions.txt like this one with how you created the partitions and the partitions results?]

## Assuming you've already done the prep in README.postprocessing.txt

### First we need to gather the identify_founders results.  We do this using the postProcessIdentifyFounders.sh script, which essetntially just concatenates all of the identify_founders.tab files from the subdirectories, but retains only one header.

## caprisa_002 v3:
./postProcessIdentifyFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/1m/partitions /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/1m/partitions 1 caprisa002 1M V3
./postProcessIdentifyFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/6m/partitions /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/6m/partitions 1 caprisa002 6M V3
./postProcessIdentifyFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/1m6m/partitions /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/1m6m/partitions 1 caprisa002 1M6M V3

## ERE I AM.  TODO: evaluateAllFoundersCap, etc.  FOR NOW JUST LOOKING AT TIMING RESULTS OF THE PARTITIONS.

### To evaluate the ancestral sequence calls, we use evaluateAllFounders.sh and evaluateAllFoundersCap.sh, and see below for Infer results.

## caprisa_002 v3:
./evaluateAllFoundersCap.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/1m/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/1m/ 1
./evaluateAllFoundersCap.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/6m/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/6m/ 1
./evaluateAllFoundersCap.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/1m6m/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/1m6m/ 1

### To evaluate PREAST's ancestral sequence calls, we use evaluateAllFoundersFromInfer.sh and evaluateAllFoundersFromInferCap.sh.

## caprisa_002 v3:
./evaluateAllFoundersFromInferCap.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/1m/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/1m/infer/ 1
./evaluateAllFoundersFromInferCap.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/6m/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/6m/infer/ 1
./evaluateAllFoundersFromInferCap.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/1m6m/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/1m6m/infer/ 1


### Now we need to gather the evaluateFounders.tbl results.  We do this using prepareForPostProcessEvaluateFounders.sh.

## caprisa_002 v3:
./prepareForPostProcessEvaluateFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/1m 
./prepareForPostProcessEvaluateFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/6m 
./prepareForPostProcessEvaluateFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/1m6m 
infer:
./prepareForPostProcessEvaluateFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/1m/infer
./prepareForPostProcessEvaluateFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/6m/infer
./prepareForPostProcessEvaluateFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/1m6m/infer

#####################################################
## Completes analysis of isMultiple for all three: nflg and both v3s:
#./evaluateIsMultiple.bash /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216
./evaluateTimings.sh 
#./postProcessEvaluateFounders.sh 

