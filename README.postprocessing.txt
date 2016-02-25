# First we need to gather the identify_founders results.  We do this manually for now (TODO: Automate it).  We just cat all the results tables together then remove the extraneous headers manually in emacs.

## v3:
cat /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/1m/hiv_founder_id_*/identify_founders.tab > /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/1m/identify_founders.tab
cat /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/6m/hiv_founder_id_*/identify_founders.tab > /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/6m/identify_founders.tab
cat /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/v3_edited_20160216/1m6m/hiv_founder_id_*/identify_founders.tab > /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/v3_edited_20160216/1m6m/identify_founders.tab

## nflg:
cat /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/nflg_copy_20160222/1m/hiv_founder_id_*/identify_founders.tab > /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/nflg_copy_20160222/1m/identify_founders.tab
cat /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/nflg_copy_20160222/6m/hiv_founder_id_*/identify_founders.tab > /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/nflg_copy_20160222/6m/identify_founders.tab
cat /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/nflg_copy_20160222/1m6m/hiv_founder_id_*/identify_founders.tab > /fh/fast/edlefsen_p/bakeoff_analysis_results/raw_edited_20160216/nflg_copy_20160222/1m6m/identify_founders.tab


