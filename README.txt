Code for preprocessing of data in McClain et al 2019. Raw data that was used can be downloaded from https://buzsakilab.nyumc.org/datasets/TingleyD/ . Exact experimental sessions used in this study can be found in “session_names” file. The following preprocessing steps are carried out in “preprocess”:

Step one: linearize session
  uses files: "linearize_track", utilities
  accomplishes: 
    linearization of position
      computes a skeletonized 2D representation of occupancy map, then represents track as a graph 
    speed (arbitrary units)
      computes difference in original coordinates then smooths and interpolates for first value
    sample trajectory for each trial type
      if map has already been computed, linearizes that. otherwise subsamples original trajectories 
      of particular trial type, then averages position and linearizes samples. Then connects graph 
      to connect points
  produces: linearized behavior file

Step two: convert speed units
  accomplishes:
    changing units of behavior to cm/s
      matches different sessions to have roughly same distribution (some conversions known, other guessed)
  amends: linearized behavior file

Step three: compute tuning of cells 
  accomplishes:
    creates tuning object
    occupancy map (unsmoothed)
      measures number of timestamps spent at each position and scales by sampling rate
    spike count map for each cell (unsmoothed)
      measures number of spikes that occur at each position
  produces: Tuning file

Step four: compute firing rate maps
  uses files: utilities
  accomplishes:
    smooths occupancy and spike counts
      gaussian smoothing along graph
    computes firing rate for each trial
      divides smoothed spike count by smoothed occupancy
    computes mean and standard error of firing rate map
    designates usable trial types
      types with more than 20 trials
  amends: behavior and tuning files
    
Step five: compute place fields
  uses files: "find_place_fields", utilities
  accomplishes:
    identifies place fields based on firing rate and reliability
      identifies groups of spatial bins where firing rate exceeds some fraction of peak firing rate,
      the tests that within each group there are at least 10 spatial bins, peak fr is above some
      value (2 hz), optionally spatial coherence above some value (disabled now), and firing rate
      is above 0 for greater than some fraction of all trials (1/3), reports labels of spatial bins
      and field associated
  amends: tuning file

NOTES: 

Consolidated data for PTP model functions is computed in “create_st_model_object”. Further code and demonstrations of PTP model can be found in https://github.com/kmcclain001/ptpModel

Cell type classification is assumed to have already occurred in the pipeline but corresponding code and parameters can be found in “classify_cell_types”, “classify_sessions” and “CellClassificationParams”

These scripts are to be used with functions and data formatting from the buzcode repository https://buzsakilab.com/wp/resources/buzcode/