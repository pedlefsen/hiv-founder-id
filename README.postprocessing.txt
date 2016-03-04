## Set up the gold standard founder results files (before running evaluateFounders)
./copyAllIntoGoldStandards.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/nflg/
./copyAllIntoGoldStandards.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/v3/
./copyAllIntoGoldStandards.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/

### First we need to gather the identify_founders results.  We do this using the postProcessIdentifyFounders.sh script, which essetntially just concatenates all of the identify_founders.tab files from the subdirectories, but retains only one header.

## nflg
./postProcessIdentifyFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/nflg/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/nflg_copy_20160222/1m/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/nflg_copy_20160222/1m/ 1
./postProcessIdentifyFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/nflg/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/nflg_copy_20160222/6m/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/nflg_copy_20160222/6m/ 1
./postProcessIdentifyFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/nflg/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/nflg_copy_20160222/1m6m/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/nflg_copy_20160222/1m6m/ 1

## caprisa_002 v3:
./postProcessIdentifyFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/1m/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/1m/ 1
./postProcessIdentifyFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/6m/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/6m/ 1
./postProcessIdentifyFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/1m6m/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/1m6m/ 1

## rv217 v3:
./postProcessIdentifyFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/v3/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/1m/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/rv217_v3_edited_20160216/1m/ 1 &
./postProcessIdentifyFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/v3/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/6m/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/rv217_v3_edited_20160216/6m/ 1 &
./postProcessIdentifyFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/v3/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/1m6m/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/rv217_v3_edited_20160216/1m6m/ 1 &

### Next we need to prepare for running evaluateTimings by extracting the true timings-es:

## caprisa_002 v3:
./prepareForEvaluateTimings.sh /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/1m/ 1 > /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/1m/sampleDates.tbl
./prepareForEvaluateTimings.sh /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/6m/ 1 > /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/6m/sampleDates.tbl
./prepareForEvaluateTimings.sh /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/1m6m/ 1 > /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/1m6m/sampleDates.tbl

## rv217 nflg:
./prepareForEvaluateTimings.sh /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/nflg_copy_20160222/1m/ 1 > /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/nflg_copy_20160222/1m/sampleDates.tbl
./prepareForEvaluateTimings.sh /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/nflg_copy_20160222/6m/ 1 > /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/nflg_copy_20160222/6m/sampleDates.tbl
./prepareForEvaluateTimings.sh /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/nflg_copy_20160222/1m6m/ 1 > /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/nflg_copy_20160222/1m6m/sampleDates.tbl

## rv217 v3 (just copy the caprisa v3 ones, which contain both):
cp /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/1m/sampleDates.tbl /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/rv217_v3_edited_20160216/1m/sampleDates.tbl 
cp /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/6m/sampleDates.tbl /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/rv217_v3_edited_20160216/6m/sampleDates.tbl 
cp /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/1m6m/sampleDates.tbl /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/rv217_v3_edited_20160216/1m6m/sampleDates.tbl 

### Next we need to symlink the v3 and nflg dirs because they are not named simply.
ln -s /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/nflg_copy_20160222 /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/nflg
ln -s /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216 /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3 
ln -s /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/rv217_v3_edited_20160216 /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/rv217_v3 

### To evaluate the ancestral sequence calls, we use evaluateAllFounders.sh and evaluateAllFoundersCap.sh, and see below for Infer results.

## nflg:
./evaluateAllFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/nflg/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/nflg_copy_20160222/1m/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/nflg_copy_20160222/1m/ 1
./evaluateAllFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/nflg/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/nflg_copy_20160222/6m/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/nflg_copy_20160222/6m/ 1
./evaluateAllFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/nflg/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/nflg_copy_20160222/1m6m/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/nflg_copy_20160222/1m6m/ 1

## caprisa_002 v3:
./evaluateAllFoundersCap.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/1m/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/1m/ 1
./evaluateAllFoundersCap.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/6m/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/6m/ 1
./evaluateAllFoundersCap.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/1m6m/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/1m6m/ 1

## rv217 v3:
./evaluateAllFoundersCap.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/v3/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/1m/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/rv217_v3_edited_20160216/1m/ 1
./evaluateAllFoundersCap.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/v3/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/6m/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/rv217_v3_edited_20160216/6m/ 1
./evaluateAllFoundersCap.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/v3/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/1m6m/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/rv217_v3_edited_20160216/1m6m/ 1

### To evaluate PREAST's ancestral sequence calls, we use evaluateAllFoundersFromInfer.sh and evaluateAllFoundersFromInferCap.sh.

## nflg:
./evaluateAllFoundersFromInfer.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/nflg/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/nflg_copy_20160222/1m /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/nflg_copy_20160222/1m/infer/ 1
./evaluateAllFoundersFromInfer.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/nflg/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/nflg_copy_20160222/6m/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/nflg_copy_20160222/6m/infer/ 1
./evaluateAllFoundersFromInfer.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/nflg/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/nflg_copy_20160222/1m6m/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/nflg_copy_20160222/1m6m/infer/ 1

## caprisa_002 v3:
./evaluateAllFoundersFromInferCap.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/1m/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/1m/infer/ 1
./evaluateAllFoundersFromInferCap.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/6m/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/6m/infer/ 1
./evaluateAllFoundersFromInferCap.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/1m6m/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/1m6m/infer/ 1

## rv217 v3:
./evaluateAllFoundersFromInferCap.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/v3/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/1m/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/rv217_v3_edited_20160216/1m/infer/ 1
./evaluateAllFoundersFromInferCap.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/v3/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/6m/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/rv217_v3_edited_20160216/6m/infer/ 1
./evaluateAllFoundersFromInferCap.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/v3/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/1m6m/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/rv217_v3_edited_20160216/1m6m/infer/ 1


### Now we need to gather the evaluateFounders.tbl results.  We do this using prepareForPostProcessEvaluateFounders.sh.

## caprisa_002 v3:
./prepareForPostProcessEvaluateFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/1m 
./prepareForPostProcessEvaluateFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/6m 
./prepareForPostProcessEvaluateFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/1m6m 

## rv217 v3:
./prepareForPostProcessEvaluateFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/v3/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/rv217_v3_edited_20160216/1m

#####################################################
## Completes analysis of isMultiple for all three: nflg and both v3s:
./evaluateIsMultiple.bash /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216
./evaluateTimings.sh

