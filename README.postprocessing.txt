## Set up the gold standard founder results files (before running evaluateFounders)
./copyAllIntoGoldStandards.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/nflg/
./copyAllIntoGoldStandards.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/v3/
./copyAllIntoGoldStandards.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/

## Copy the bounds.
tar xzf /fast/bakeoff_merged_analysis_sequences_results/raw_fixed/bounds.tar.gz -C /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/
tar xzf /fast/bakeoff_merged_analysis_sequences_results/raw_fixed/bounds.tar.gz -C /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/

###### REQUIRED FOR EVERYTHING, EVEN IF SKIPPING evaluateFounders (still needed for evaluateIsMultiple and evaluateTimings): ###################
######
### First we need to gather the identify_founders results.  We do this using the postProcessIdentifyFounders.sh script, which essetntially just concatenates all of the identify_founders.tab files from the subdirectories, but retains only one header.

### This is for the main results, at this time called raw_fixed.  For the partitions, see README.postprocessing_partitions.txt. We do this for two ways of hypermutation/recombination:
### 1) in /fast/bakeoff_merged_analysis_sequences_results/results/ we have it using hypermut2.0 and RAP Beta (pre-2017). We used "remove" rather than "fix" for handling hypermutants.
### 2) in /fast/bakeoff_merged_analysis_sequences_results_2019/results/ we have it using hypermutR and RAPR. This time we used the "fix" rather than "remove" hypermutR option.
### 0) Note also the previous results used no filtering, in /fast/bakeoff_analysis_results/raw_edited_20160216/

## Note that for our purposes here in #1 and #2, postProcessIdentifyFounders.sh takes the same argument repeated three times: the processed_lists dir. For #0 the args are different, though.

## SEE BELOW WHERE WE CLEAN THE LISTS TOO

### 0) Results with no filtering.
./postProcessIdentifyFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/nflg/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/nflg_copy_20160222/1m/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/nflg_copy_20160222/1m/ 1
./postProcessIdentifyFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/nflg/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/nflg_copy_20160222/6m/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/nflg_copy_20160222/6m/ 1
./postProcessIdentifyFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/nflg/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/nflg_copy_20160222/1m6m/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/nflg_copy_20160222/1m6m/ 1

### Next we need to symlink the v3 and nflg dirs because they are not named simply.
ln -s /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/nflg_copy_20160222 /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/nflg
ln -s /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216 /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3 
ln -s /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/rv217_v3_edited_20160216 /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/rv217_v3 

### 1) pre-2017 results
## nflg
./postProcessIdentifyFounders.sh /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/nflg/1m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/nflg/1m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/nflg/1m
./postProcessIdentifyFounders.sh /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/nflg/6m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/nflg/6m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/nflg/6m
./postProcessIdentifyFounders.sh /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/nflg/1m6m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/nflg/1m6m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/nflg/1m6m

## v3
./postProcessIdentifyFounders.sh /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/v3/1m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/v3/1m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/v3/1m
./postProcessIdentifyFounders.sh /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/v3/6m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/v3/6m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/v3/6m
./postProcessIdentifyFounders.sh /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/v3/1m6m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/v3/1m6m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/v3/1m6m

## rv217_v3
./postProcessIdentifyFounders.sh /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/rv217_v3/1m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/rv217_v3/1m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/rv217_v3/1m
./postProcessIdentifyFounders.sh /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/rv217_v3/6m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/rv217_v3/6m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/rv217_v3/6m
./postProcessIdentifyFounders.sh /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/rv217_v3/1m6m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/rv217_v3/1m6m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/rv217_v3/1m6m

### 2) post-2017 results
## nflg
./postProcessIdentifyFounders.sh /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/nflg/1m /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/nflg/1m /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/nflg/1m
./postProcessIdentifyFounders.sh /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/nflg/6m /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/nflg/6m /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/nflg/6m
./postProcessIdentifyFounders.sh /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/nflg/1m6m /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/nflg/1m6m /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/nflg/1m6m

## v3
./postProcessIdentifyFounders.sh /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/v3/1m /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/v3/1m /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/v3/1m
./postProcessIdentifyFounders.sh /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/v3/6m /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/v3/6m /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/v3/6m
./postProcessIdentifyFounders.sh /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/v3/1m6m /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/v3/1m6m /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/v3/1m6m

### CLEAN THE LISTS
## !!For now this requires opening the file cleanLists.R and changing the SEQUENCES.DIR and RESULTS.DIRNAME entries.
### 0) Results with no filtering.
#SEQUENCES.DIR <- "/fh/fast/edlefsen_p/bakeoff/analysis_sequences/";
#RESULTS.DIRNAME <- "raw_edited_20160216";
cleanLists.sh
### 1) pre-2017 filtering
#SEQUENCES.DIR <- "/fast/bakeoff_merged_analysis_sequences_results/results/";
#RESULTS.DIRNAME <- "raw_fixed";
cleanLists.sh
### 2) post-2017 (2019) filtering
#SEQUENCES.DIR <- "/fast/bakeoff_merged_analysis_sequences_results_2019/results/";
#RESULTS.DIRNAME <- "raw_fixed";
cleanLists.sh



###### END REQUIRED FOR EVERYTHING, EVEN IF SKIPPING evaluateFounders (still needed for evaluateIsMultiple and evaluateTimings): ###################

###### FOR evaluateFounders ONLY: ######################3
######
### First ensure you have done everything in README.preprocessing.txt above this point.

### To evaluate the ancestral sequence calls, we use evaluateAllFounders.sh and evaluateAllFoundersCap.sh, and see below for Infer results.

### 0) unfiltered:
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

### 1) pre-2017 filtered:
## nflg:
./evaluateAllFounders.sh /fast/bakeoff/gold_standard/rv217/nflg/ /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/nflg/1m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/nflg/1m
./evaluateAllFounders.sh /fast/bakeoff/gold_standard/rv217/nflg/ /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/nflg/6m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/nflg/6m
./evaluateAllFounders.sh /fast/bakeoff/gold_standard/rv217/nflg/ /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/nflg/1m6m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/nflg/1m6m

## v3:
./evaluateAllFounders.sh /fast/bakeoff/gold_standard/rv217/v3/ /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/v3/1m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/v3/1m
./evaluateAllFounders.sh /fast/bakeoff/gold_standard/rv217/v3/ /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/v3/6m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/v3/6m
./evaluateAllFounders.sh /fast/bakeoff/gold_standard/rv217/v3/ /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/v3/1m6m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/v3/1m6m

## rv217_v3:
./evaluateAllFounders.sh /fast/bakeoff/gold_standard/rv217/rv217_v3/ /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/rv217_v3/1m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/rv217_v3/1m
./evaluateAllFounders.sh /fast/bakeoff/gold_standard/rv217/rv217_v3/ /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/rv217_v3/6m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/rv217_v3/6m
./evaluateAllFounders.sh /fast/bakeoff/gold_standard/rv217/rv217_v3/ /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/rv217_v3/1m6m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/rv217_v3/1m6m


### 2) post-2017 filtered:
## nflg:
./evaluateAllFounders.sh /fast/bakeoff/gold_standard/rv217/nflg/ /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/nflg/1m /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/nflg/1m
./evaluateAllFounders.sh /fast/bakeoff/gold_standard/rv217/nflg/ /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/nflg/6m /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/nflg/6m
./evaluateAllFounders.sh /fast/bakeoff/gold_standard/rv217/nflg/ /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/nflg/1m6m /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/nflg/1m6m

## v3:
./evaluateAllFounders.sh /fast/bakeoff/gold_standard/rv217/v3/ /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/v3/1m /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/v3/1m
./evaluateAllFounders.sh /fast/bakeoff/gold_standard/rv217/v3/ /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/v3/6m /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/v3/6m
./evaluateAllFounders.sh /fast/bakeoff/gold_standard/rv217/v3/ /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/v3/1m6m /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/v3/1m6m

### To evaluate PREAST's ancestral sequence calls, we use evaluateAllFoundersFromInfer.sh and evaluateAllFoundersFromInferCap.sh.

### 0) unfiltered:
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

### 1) pre-2017 filtered:
## nflg:
./evaluateAllFoundersFromInfer.sh /fast/bakeoff/gold_standard/rv217/nflg/ /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/nflg/1m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/nflg/1m
./evaluateAllFoundersFromInfer.sh /fast/bakeoff/gold_standard/rv217/nflg/ /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/nflg/6m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/nflg/6m
./evaluateAllFoundersFromInfer.sh /fast/bakeoff/gold_standard/rv217/nflg/ /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/nflg/1m6m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/nflg/1m6m

## v3:
./evaluateAllFoundersFromInfer.sh /fast/bakeoff/gold_standard/rv217/v3/ /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/v3/1m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/v3/1m
./evaluateAllFoundersFromInfer.sh /fast/bakeoff/gold_standard/rv217/v3/ /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/v3/6m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/v3/6m
./evaluateAllFoundersFromInfer.sh /fast/bakeoff/gold_standard/rv217/v3/ /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/v3/1m6m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/v3/1m6m

## rv217_v3:
./evaluateAllFoundersFromInfer.sh /fast/bakeoff/gold_standard/rv217/rv217_v3/ /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/rv217_v3/1m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/rv217_v3/1m
./evaluateAllFoundersFromInfer.sh /fast/bakeoff/gold_standard/rv217/rv217_v3/ /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/rv217_v3/6m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/rv217_v3/6m
./evaluateAllFoundersFromInfer.sh /fast/bakeoff/gold_standard/rv217/rv217_v3/ /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/rv217_v3/1m6m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/rv217_v3/1m6m


### 2) post-2017 filtered:
## nflg:
./evaluateAllFoundersFromInfer.sh /fast/bakeoff/gold_standard/rv217/nflg/ /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/nflg/1m /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/nflg/1m
./evaluateAllFoundersFromInfer.sh /fast/bakeoff/gold_standard/rv217/nflg/ /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/nflg/6m /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/nflg/6m
./evaluateAllFoundersFromInfer.sh /fast/bakeoff/gold_standard/rv217/nflg/ /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/nflg/1m6m /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/nflg/1m6m

## v3:
./evaluateAllFoundersFromInfer.sh /fast/bakeoff/gold_standard/rv217/v3/ /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/v3/1m /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/v3/1m
./evaluateAllFoundersFromInfer.sh /fast/bakeoff/gold_standard/rv217/v3/ /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/v3/6m /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/v3/6m
./evaluateAllFoundersFromInfer.sh /fast/bakeoff/gold_standard/rv217/v3/ /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/v3/1m6m /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/v3/1m6m

### Now we need to gather the evaluateFounders.tbl results.  We do this using prepareForPostProcessEvaluateFounders.sh.

### 0) unfiltered:
## nflg:
./prepareForPostProcessEvaluateFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/nflg/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/nflg_copy_20160222/1m
./prepareForPostProcessEvaluateFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/nflg/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/nflg_copy_20160222/6m
./prepareForPostProcessEvaluateFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/nflg/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/nflg_copy_20160222/1m6m
# infer:
./prepareForPostProcessEvaluateFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/nflg/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/nflg_copy_20160222/1m/infer
./prepareForPostProcessEvaluateFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/nflg/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/nflg_copy_20160222/6m/infer
./prepareForPostProcessEvaluateFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/nflg/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/nflg_copy_20160222/1m6m/infer

## caprisa_002 v3:
./prepareForPostProcessEvaluateFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/1m 
./prepareForPostProcessEvaluateFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/6m 
./prepareForPostProcessEvaluateFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/1m6m 
infer:
./prepareForPostProcessEvaluateFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/1m/infer
./prepareForPostProcessEvaluateFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/6m/infer
./prepareForPostProcessEvaluateFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/1m6m/infer

## rv217 v3:
./prepareForPostProcessEvaluateFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/v3/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/rv217_v3_edited_20160216/1m
./prepareForPostProcessEvaluateFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/v3/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/rv217_v3_edited_20160216/6m
./prepareForPostProcessEvaluateFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/v3/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/rv217_v3_edited_20160216/1m6m
# infer:
./prepareForPostProcessEvaluateFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/v3/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/rv217_v3_edited_20160216/1m/infer
./prepareForPostProcessEvaluateFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/v3/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/rv217_v3_edited_20160216/6m/infer
./prepareForPostProcessEvaluateFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/v3/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/rv217_v3_edited_20160216/1m6m/infer

### 1) pre-2017 filtered:
## nflg:
./prepareForPostProcessEvaluateFounders.sh /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/nflg/1m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/nflg/1m
./prepareForPostProcessEvaluateFounders.sh /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/nflg/6m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/nflg/6m
./prepareForPostProcessEvaluateFounders.sh /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/nflg/1m6m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/nflg/1m6m

## v3:
./prepareForPostProcessEvaluateFounders.sh /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/v3/1m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/v3/1m
./prepareForPostProcessEvaluateFounders.sh /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/v3/6m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/v3/6m
./prepareForPostProcessEvaluateFounders.sh /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/v3/1m6m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/v3/1m6m

## rv217_v3:
./prepareForPostProcessEvaluateFounders.sh /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/rv217_v3/1m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/rv217_v3/1m
./prepareForPostProcessEvaluateFounders.sh /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/rv217_v3/6m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/rv217_v3/6m
./prepareForPostProcessEvaluateFounders.sh /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/rv217_v3/1m6m /fast/bakeoff_merged_analysis_sequences_results/results/raw_fixed/rv217_v3/1m6m

### 2) post-2017 filtered:
## nflg:
./prepareForPostProcessEvaluateFounders.sh /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/nflg/1m /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/nflg/1m
./prepareForPostProcessEvaluateFounders.sh /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/nflg/6m /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/nflg/6m
./prepareForPostProcessEvaluateFounders.sh /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/nflg/1m6m /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/nflg/1m6m

## v3:
./prepareForPostProcessEvaluateFounders.sh /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/v3/1m /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/v3/1m
./prepareForPostProcessEvaluateFounders.sh /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/v3/6m /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/v3/6m
./prepareForPostProcessEvaluateFounders.sh /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/v3/1m6m /fast/bakeoff_merged_analysis_sequences_results_2019/results/raw_fixed/v3/1m6m


###### END FOR evaluateFounders ONLY ######################3

############################################################
## FOR ANALYSIS OF evaluateIsMultiple and evaluateTimings ##
############################################################
## Completes analysis of isMultiple, creates *_isMultiple.tab results:
### NOTE THIS REQUIRES HACKING THE FILE evaluateIsMultiple.R to set RESULTS.DIR
# 0) unfiltered
### !!!!!!!! SET THIS FIRST in evaluateIsMultiple.R !!!!!
### RESULTS.DIR <- "/fh/fast/edlefsen_p/bakeoff_analysis_results/";
### RESULTS.DIRNAME <- "raw_edited_20160216";
./evaluateIsMultiple.sh
# 1) pre-2017 filtered
### !!!!!!!! SET THIS FIRST in evaluateIsMultiple.R !!!!!
### RESULTS.DIR <- "/fast/bakeoff_merged_analysis_sequences_results/results/";
### RESULTS.DIRNAME <- "raw_fixed";
./evaluateIsMultiple.sh
# 2) post-2017 filtered
### !!!!!!!! SET THIS FIRST in evaluateIsMultiple.R !!!!!
### RESULTS.DIR <- "/fast/bakeoff_merged_analysis_sequences_results_2019/results/";
### RESULTS.DIRNAME <- "raw_fixed";
./evaluateIsMultiple.sh

#### Completes analysis of timings, creates */*/evaluateTimings.tab results:
### DO EACH OF THESE FOUR TIMES, FOR EACH COMBO OF HELPFUL.ADDITIONAL.COLS and INCLUDE.INTERCEPT (set in evalateTimings.R:)
#HELPFUL.ADDITIONAL.COLS <- c( "lPVL" );
#HELPFUL.ADDITIONAL.COLS <- c();
#
#INCLUDE.INTERCEPT <- FALSE;
#INCLUDE.INTERCEPT <- TRUE;

# 0) unfiltered
### !!!!!!!! SET THESE FIRST in evaluateTimings.R !!!!!
### SEQUENCES.DIR <- "/fh/fast/edlefsen_p/bakeoff/analysis_sequences/";
### RESULTS.DIR <- "/fh/fast/edlefsen_p/bakeoff_analysis_results/";
### RESULTS.DIRNAME <- "raw_edited_20160216";
./evaluateTimings.sh
### !!!!!!!! SET THIS FIRST in evaluateTimings.R !!!!!
### RESULTS.DIR <- SEQUENCES.DIR <- "/fast/bakeoff_merged_analysis_sequences_results/results/"; RESULTS.DIRNAME <- "raw_fixed";
./evaluateTimings.sh
### !!!!!!!! SET THIS FIRST in evaluateTimings.R !!!!!
### RESULTS.DIR <- SEQUENCES.DIR <- "/fast/bakeoff_merged_analysis_sequences_results_2019/results/"; RESULTS.DIRNAME <- "raw_fixed";
./evaluateTimings.sh

##########################################
## FOR ANALYSIS OF evaluateFounders ONLY #
##########################################
## Completes analysis of founders, creates */*/postProcessEvaluateFounders.tab results:
./postProcessEvaluateFounders.sh 


