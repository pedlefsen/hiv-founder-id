## Set up the gold standard founder results files (before running evaluateFounders)
./copyAllIntoGoldStandards.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/nflg/
./copyAllIntoGoldStandards.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/v3/
./copyAllIntoGoldStandards.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/v3/

### CLEAN THE LISTS
## !!For now this requires opening the file cleanLists.R and changing the SEQUENCES.DIR and RESULTS.DIRNAME entries.
### 0) Results with no filtering.
#SEQUENCES.DIR <- "/fast/bakeoff_merged_analysis_sequences_unfiltered/results/"; 
cleanLists.sh
### 1) pre-2017 filtering
#SEQUENCES.DIR <- "/fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/";  
cleanLists.sh
### 2) post-2017 (2019) filtering
#SEQUENCES.DIR <- "/fast/bakeoff_merged_analysis_sequences_filtered2019/results/";
cleanLists.sh

## Copy the bounds.
tar xzf /fast/bakeoff_merged_analysis_sequences_results/raw_fixed/bounds.tar.gz -C /fast/bakeoff_merged_analysis_sequences_filteredPre2017/results/raw_fixed/
tar xzf /fast/bakeoff_merged_analysis_sequences_results/raw_fixed/bounds.tar.gz -C /fast/bakeoff_merged_analysis_sequences_filtered2019/results/raw_fixed/

## FOR ALL RESULTS EXCEPT evaluateFounders (not in the paper), see notes in files
# /fast/bakeoff_merged_analysis_sequences_unfiltered/archive/NOTES.DOPTE.unfiltered.txt
# /fast/bakeoff_merged_analysis_sequences_filteredPre2017/archive/NOTES.DOPTE.filteredPre2017.txt 
# /fast/bakeoff_merged_analysis_sequences_filtered2019/archive/NOTES.DOPTE.filtered2019.txt 

###### FOR evaluateFounders ONLY: ######################
## See README.postprocessing_evaluateFounders.txt
