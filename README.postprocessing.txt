### First we need to gather the identify_founders results.  We do this manually for now (TODO: Automate it).  We just cat all the results tables together then remove the extraneous headers manually in emacs.

## v3:
cat /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/1m/hiv_founder_id_*/identify_founders.tab > /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/1m/identify_founders.tab
cat /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/6m/hiv_founder_id_*/identify_founders.tab > /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/6m/identify_founders.tab
cat /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/1m6m/hiv_founder_id_*/identify_founders.tab > /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/1m6m/identify_founders.tab

## nflg:
cat /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/nflg_copy_20160222/1m/hiv_founder_id_*/identify_founders.tab > /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/nflg_copy_20160222/1m/identify_founders.tab
cat /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/nflg_copy_20160222/6m/hiv_founder_id_*/identify_founders.tab > /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/nflg_copy_20160222/6m/identify_founders.tab
cat /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/nflg_copy_20160222/1m6m/hiv_founder_id_*/identify_founders.tab > /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/nflg_copy_20160222/1m6m/identify_founders.tab

### Next we need to prepare for running evaluateTimings by extracting the true timings-es:

## v3:
./evaluateTimings.sh /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/1m/ 1 > /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/1m/sampleDates.tbl
./evaluateTimings.sh /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/6m/ 1 > /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/6m/sampleDates.tbl
./evaluateTimings.sh /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/1m6m/ 1 > /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/1m6m/sampleDates.tbl

## nflg:
./evaluateTimings.sh /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/nflg_copy_20160222/1m/ 1 > /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/nflg_copy_20160222/1m/sampleDates.tbl
./evaluateTimings.sh /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/nflg_copy_20160222/6m/ 1 > /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/nflg_copy_20160222/6m/sampleDates.tbl
./evaluateTimings.sh /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/nflg_copy_20160222/1m6m/ 1 > /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/nflg_copy_20160222/1m6m/sampleDates.tbl

### Next we need to symlink the v3 and nflg dirs because they are not named simply.
ln -s /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/nflg_copy_20160222 /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/nflg
ln -s /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216 /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3

### To evaluate the ancestral sequence calls, we use evaluateAllFounders.sh and evaluateAllFoundersCap.sh
## nflg:
./evaluateAllFounders.sh /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/nflg/ /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/nflg_copy_20160222/1m/ /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/nflg_copy_20160222/1m/ 1
## ERE I AM 
./evaluateFounders.bash 
