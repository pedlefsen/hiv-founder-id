###### FOR evaluateFounders ONLY: ######################
### First ensure you have done everything else in README.preprocessing.txt first.

### To evaluate the ancestral sequence calls, we use evaluateAllFounders.sh and evaluateAllFoundersCap.sh, and see below for Infer results.

### 0) unfiltered:
## nflg:
./evaluateAllFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/nflg/ /fast/bakeoff_merged_analysis_sequences_unfiltered/results/raw_fixed/nflg/1m/ /fast/bakeoff_merged_analysis_sequences_unfiltered/results/raw_fixed/nflg/1m/ 1
./evaluateAllFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/nflg/ /fast/bakeoff_merged_analysis_sequences_unfiltered/results/raw_fixed/nflg/6m/ /fast/bakeoff_merged_analysis_sequences_unfiltered/results/raw_fixed/nflg/6m/ 1
./evaluateAllFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/nflg/ /fast/bakeoff_merged_analysis_sequences_unfiltered/results/raw_fixed/nflg/1m6m/ /fast/bakeoff_merged_analysis_sequences_unfiltered/results/raw_fixed/nflg/1m6m/ 1

## caprisa_002 v3:
./evaluateAllFoundersCap.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/1m/ /fast/bakeoff_merged_analysis_sequences_unfiltered/results/raw_fixed/v3/1m/ 1
./evaluateAllFoundersCap.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/6m/ /fast/bakeoff_merged_analysis_sequences_unfiltered/results/raw_fixed/v3/6m/ 1
./evaluateAllFoundersCap.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/1m6m/ /fast/bakeoff_merged_analysis_sequences_unfiltered/results/raw_fixed/v3/1m6m/ 1

## rv217 v3:
./evaluateAllFoundersCap.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/v3/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/1m/ /fast/bakeoff_merged_analysis_sequences_unfiltered/results/raw_fixed/v3/1m/ 1
./evaluateAllFoundersCap.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/v3/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/6m/ /fast/bakeoff_merged_analysis_sequences_unfiltered/results/raw_fixed/v3/6m/ 1
./evaluateAllFoundersCap.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/v3/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/1m6m/ /fast/bakeoff_merged_analysis_sequences_unfiltered/results/raw_fixed/v3/1m6m/ 1

### 1) pre-2017 filtered:
## nflg:
./evaluateAllFounders.sh /fast/bakeoff/gold_standard/rv217/nflg/ /fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/raw_fixed/nflg/1m /fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/raw_fixed/nflg/1m
./evaluateAllFounders.sh /fast/bakeoff/gold_standard/rv217/nflg/ /fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/raw_fixed/nflg/6m /fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/raw_fixed/nflg/6m
./evaluateAllFounders.sh /fast/bakeoff/gold_standard/rv217/nflg/ /fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/raw_fixed/nflg/1m6m /fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/raw_fixed/nflg/1m6m

## v3:
./evaluateAllFounders.sh /fast/bakeoff/gold_standard/rv217/v3/ /fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/raw_fixed/v3/1m /fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/raw_fixed/v3/1m
./evaluateAllFounders.sh /fast/bakeoff/gold_standard/rv217/v3/ /fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/raw_fixed/v3/6m /fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/raw_fixed/v3/6m
./evaluateAllFounders.sh /fast/bakeoff/gold_standard/rv217/v3/ /fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/raw_fixed/v3/1m6m /fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/raw_fixed/v3/1m6m


### 2) post-2017 filtered:
## nflg:
./evaluateAllFounders.sh /fast/bakeoff/gold_standard/rv217/nflg/ /fast/bakeoff_merged_analysis_sequences_filtered2019/results/raw_fixed/nflg/1m /fast/bakeoff_merged_analysis_sequences_filtered2019/results/raw_fixed/nflg/1m
./evaluateAllFounders.sh /fast/bakeoff/gold_standard/rv217/nflg/ /fast/bakeoff_merged_analysis_sequences_filtered2019/results/raw_fixed/nflg/6m /fast/bakeoff_merged_analysis_sequences_filtered2019/results/raw_fixed/nflg/6m
./evaluateAllFounders.sh /fast/bakeoff/gold_standard/rv217/nflg/ /fast/bakeoff_merged_analysis_sequences_filtered2019/results/raw_fixed/nflg/1m6m /fast/bakeoff_merged_analysis_sequences_filtered2019/results/raw_fixed/nflg/1m6m

## v3:
./evaluateAllFounders.sh /fast/bakeoff/gold_standard/rv217/v3/ /fast/bakeoff_merged_analysis_sequences_filtered2019/results/raw_fixed/v3/1m /fast/bakeoff_merged_analysis_sequences_filtered2019/results/raw_fixed/v3/1m
./evaluateAllFounders.sh /fast/bakeoff/gold_standard/rv217/v3/ /fast/bakeoff_merged_analysis_sequences_filtered2019/results/raw_fixed/v3/6m /fast/bakeoff_merged_analysis_sequences_filtered2019/results/raw_fixed/v3/6m
./evaluateAllFounders.sh /fast/bakeoff/gold_standard/rv217/v3/ /fast/bakeoff_merged_analysis_sequences_filtered2019/results/raw_fixed/v3/1m6m /fast/bakeoff_merged_analysis_sequences_filtered2019/results/raw_fixed/v3/1m6m

### To evaluate PREAST's ancestral sequence calls, we use evaluateAllFoundersFromInfer.sh and evaluateAllFoundersFromInferCap.sh.

### 0) unfiltered:
## nflg:
./evaluateAllFoundersFromInfer.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/nflg/ /fast/bakeoff_merged_analysis_sequences_unfiltered/results/raw_fixed/nflg/1m /fast/bakeoff_merged_analysis_sequences_unfiltered/results/raw_fixed/nflg/1m/infer/ 1
./evaluateAllFoundersFromInfer.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/nflg/ /fast/bakeoff_merged_analysis_sequences_unfiltered/results/raw_fixed/nflg/6m/ /fast/bakeoff_merged_analysis_sequences_unfiltered/results/raw_fixed/nflg/6m/infer/ 1
./evaluateAllFoundersFromInfer.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/nflg/ /fast/bakeoff_merged_analysis_sequences_unfiltered/results/raw_fixed/nflg/1m6m/ /fast/bakeoff_merged_analysis_sequences_unfiltered/results/raw_fixed/nflg/1m6m/infer/ 1

## caprisa_002 v3:
./evaluateAllFoundersFromInferCap.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/1m/ /fast/bakeoff_merged_analysis_sequences_unfiltered/results/raw_fixed/v3/1m/infer/ 1
./evaluateAllFoundersFromInferCap.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/6m/ /fast/bakeoff_merged_analysis_sequences_unfiltered/results/raw_fixed/v3/6m/infer/ 1
./evaluateAllFoundersFromInferCap.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/1m6m/ /fast/bakeoff_merged_analysis_sequences_unfiltered/results/raw_fixed/v3/1m6m/infer/ 1

## rv217 v3:
./evaluateAllFoundersFromInferCap.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/v3/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/1m/ /fast/bakeoff_merged_analysis_sequences_unfiltered/results/raw_fixed/v3/1m/infer/ 1
./evaluateAllFoundersFromInferCap.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/v3/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/6m/ /fast/bakeoff_merged_analysis_sequences_unfiltered/results/raw_fixed/v3/6m/infer/ 1
./evaluateAllFoundersFromInferCap.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/v3/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/1m6m/ /fast/bakeoff_merged_analysis_sequences_unfiltered/results/raw_fixed/v3/1m6m/infer/ 1

### 1) pre-2017 filtered:
## nflg:
./evaluateAllFoundersFromInfer.sh /fast/bakeoff/gold_standard/rv217/nflg/ /fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/raw_fixed/nflg/1m /fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/raw_fixed/nflg/1m
./evaluateAllFoundersFromInfer.sh /fast/bakeoff/gold_standard/rv217/nflg/ /fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/raw_fixed/nflg/6m /fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/raw_fixed/nflg/6m
./evaluateAllFoundersFromInfer.sh /fast/bakeoff/gold_standard/rv217/nflg/ /fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/raw_fixed/nflg/1m6m /fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/raw_fixed/nflg/1m6m

## v3:
./evaluateAllFoundersFromInfer.sh /fast/bakeoff/gold_standard/rv217/v3/ /fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/raw_fixed/v3/1m /fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/raw_fixed/v3/1m
./evaluateAllFoundersFromInfer.sh /fast/bakeoff/gold_standard/rv217/v3/ /fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/raw_fixed/v3/6m /fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/raw_fixed/v3/6m
./evaluateAllFoundersFromInfer.sh /fast/bakeoff/gold_standard/rv217/v3/ /fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/raw_fixed/v3/1m6m /fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/raw_fixed/v3/1m6m


### 2) post-2017 filtered:
## nflg:
./evaluateAllFoundersFromInfer.sh /fast/bakeoff/gold_standard/rv217/nflg/ /fast/bakeoff_merged_analysis_sequences_filtered2019/results/raw_fixed/nflg/1m /fast/bakeoff_merged_analysis_sequences_filtered2019/results/raw_fixed/nflg/1m
./evaluateAllFoundersFromInfer.sh /fast/bakeoff/gold_standard/rv217/nflg/ /fast/bakeoff_merged_analysis_sequences_filtered2019/results/raw_fixed/nflg/6m /fast/bakeoff_merged_analysis_sequences_filtered2019/results/raw_fixed/nflg/6m
./evaluateAllFoundersFromInfer.sh /fast/bakeoff/gold_standard/rv217/nflg/ /fast/bakeoff_merged_analysis_sequences_filtered2019/results/raw_fixed/nflg/1m6m /fast/bakeoff_merged_analysis_sequences_filtered2019/results/raw_fixed/nflg/1m6m

## v3:
./evaluateAllFoundersFromInfer.sh /fast/bakeoff/gold_standard/rv217/v3/ /fast/bakeoff_merged_analysis_sequences_filtered2019/results/raw_fixed/v3/1m /fast/bakeoff_merged_analysis_sequences_filtered2019/results/raw_fixed/v3/1m
./evaluateAllFoundersFromInfer.sh /fast/bakeoff/gold_standard/rv217/v3/ /fast/bakeoff_merged_analysis_sequences_filtered2019/results/raw_fixed/v3/6m /fast/bakeoff_merged_analysis_sequences_filtered2019/results/raw_fixed/v3/6m
./evaluateAllFoundersFromInfer.sh /fast/bakeoff/gold_standard/rv217/v3/ /fast/bakeoff_merged_analysis_sequences_filtered2019/results/raw_fixed/v3/1m6m /fast/bakeoff_merged_analysis_sequences_filtered2019/results/raw_fixed/v3/1m6m

### Now we need to gather the evaluateFounders.tbl results.  We do this using prepareForPostProcessEvaluateFounders.sh.

### 0) unfiltered:
## nflg:
./prepareForPostProcessEvaluateFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/nflg/ /fast/bakeoff_merged_analysis_sequences_unfiltered/results/raw_fixed/nflg/1m
./prepareForPostProcessEvaluateFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/nflg/ /fast/bakeoff_merged_analysis_sequences_unfiltered/results/raw_fixed/nflg/6m
./prepareForPostProcessEvaluateFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/nflg/ /fast/bakeoff_merged_analysis_sequences_unfiltered/results/raw_fixed/nflg/1m6m
# infer:
./prepareForPostProcessEvaluateFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/nflg/ /fast/bakeoff_merged_analysis_sequences_unfiltered/results/raw_fixed/nflg/1m/infer
./prepareForPostProcessEvaluateFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/nflg/ /fast/bakeoff_merged_analysis_sequences_unfiltered/results/raw_fixed/nflg/6m/infer
./prepareForPostProcessEvaluateFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/nflg/ /fast/bakeoff_merged_analysis_sequences_unfiltered/results/raw_fixed/nflg/1m6m/infer

## caprisa_002 v3:
./prepareForPostProcessEvaluateFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fast/bakeoff_merged_analysis_sequences_unfiltered/results/raw_fixed/v3/1m 
./prepareForPostProcessEvaluateFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fast/bakeoff_merged_analysis_sequences_unfiltered/results/raw_fixed/v3/6m 
./prepareForPostProcessEvaluateFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fast/bakeoff_merged_analysis_sequences_unfiltered/results/raw_fixed/v3/1m6m 
infer:
./prepareForPostProcessEvaluateFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fast/bakeoff_merged_analysis_sequences_unfiltered/results/raw_fixed/v3/1m/infer
./prepareForPostProcessEvaluateFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fast/bakeoff_merged_analysis_sequences_unfiltered/results/raw_fixed/v3/6m/infer
./prepareForPostProcessEvaluateFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fast/bakeoff_merged_analysis_sequences_unfiltered/results/raw_fixed/v3/1m6m/infer

### 1) pre-2017 filtered:
## nflg:
./prepareForPostProcessEvaluateFounders.sh /fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/raw_fixed/nflg/1m /fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/raw_fixed/nflg/1m
./prepareForPostProcessEvaluateFounders.sh /fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/raw_fixed/nflg/6m /fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/raw_fixed/nflg/6m
./prepareForPostProcessEvaluateFounders.sh /fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/raw_fixed/nflg/1m6m /fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/raw_fixed/nflg/1m6m

## v3:
./prepareForPostProcessEvaluateFounders.sh /fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/raw_fixed/v3/1m /fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/raw_fixed/v3/1m
./prepareForPostProcessEvaluateFounders.sh /fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/raw_fixed/v3/6m /fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/raw_fixed/v3/6m
./prepareForPostProcessEvaluateFounders.sh /fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/raw_fixed/v3/1m6m /fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/raw_fixed/v3/1m6m

### 2) post-2017 filtered:
## nflg:
./prepareForPostProcessEvaluateFounders.sh /fast/bakeoff_merged_analysis_sequences_filtered2019/results/raw_fixed/nflg/1m /fast/bakeoff_merged_analysis_sequences_filtered2019/results/raw_fixed/nflg/1m
./prepareForPostProcessEvaluateFounders.sh /fast/bakeoff_merged_analysis_sequences_filtered2019/results/raw_fixed/nflg/6m /fast/bakeoff_merged_analysis_sequences_filtered2019/results/raw_fixed/nflg/6m
./prepareForPostProcessEvaluateFounders.sh /fast/bakeoff_merged_analysis_sequences_filtered2019/results/raw_fixed/nflg/1m6m /fast/bakeoff_merged_analysis_sequences_filtered2019/results/raw_fixed/nflg/1m6m

## v3:
./prepareForPostProcessEvaluateFounders.sh /fast/bakeoff_merged_analysis_sequences_filtered2019/results/raw_fixed/v3/1m /fast/bakeoff_merged_analysis_sequences_filtered2019/results/raw_fixed/v3/1m
./prepareForPostProcessEvaluateFounders.sh /fast/bakeoff_merged_analysis_sequences_filtered2019/results/raw_fixed/v3/6m /fast/bakeoff_merged_analysis_sequences_filtered2019/results/raw_fixed/v3/6m
./prepareForPostProcessEvaluateFounders.sh /fast/bakeoff_merged_analysis_sequences_filtered2019/results/raw_fixed/v3/1m6m /fast/bakeoff_merged_analysis_sequences_filtered2019/results/raw_fixed/v3/1m6m


## Completes analysis of founders, creates */*/postProcessEvaluateFounders.tab results:
./postProcessEvaluateFounders.sh 


