### First we need to symlink the v3 and nflg dirs because they are not named simply.
ln -s /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/nflg_copy_20160222 /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/nflg
ln -s /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216 /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3 

### Next we need to prepare for createArtificialBoundsOnInfectionDate.R (and also evaluateTimings.R, in the postprocessing phase) by extracting the true timings-es:

## caprisa_002 v3:
./extractSampleDates.sh /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/1m/ 1 > /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/1m/sampleDates.tbl
./extractSampleDates.sh /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/6m/ 1 > /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/6m/sampleDates.tbl
./extractSampleDates.sh /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/1m6m/ 1 > /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/1m6m/sampleDates.tbl

## rv217 nflg:
./extractSampleDates.sh /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/nflg_copy_20160222/1m/ 1 > /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/nflg_copy_20160222/1m/sampleDates.tbl
./extractSampleDates.sh /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/nflg_copy_20160222/6m/ 1 > /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/nflg_copy_20160222/6m/sampleDates.tbl
./extractSampleDates.sh /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/nflg_copy_20160222/1m6m/ 1 > /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/nflg_copy_20160222/1m6m/sampleDates.tbl

## rv217 v3 (just copy the caprisa v3 ones, which contain both):
cp /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/1m/sampleDates.tbl /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/rv217_v3_edited_20160216/1m/sampleDates.tbl 
cp /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/6m/sampleDates.tbl /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/rv217_v3_edited_20160216/6m/sampleDates.tbl 
cp /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/1m6m/sampleDates.tbl /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/rv217_v3_edited_20160216/1m6m/sampleDates.tbl 

