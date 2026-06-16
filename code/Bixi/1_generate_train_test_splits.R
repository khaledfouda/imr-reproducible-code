source("./code/Bixi/helper.R")
source("./code/helper.R")

# The following will generate all train/test splits for each of the 5 missing data scenarios.
# Total output files are 300.  50*5=250 training data sets + 50 test sets (the test set is fixed for each split)


for(split_id in 1:NUM_SPLITS){
  bixi_generate_one_split(split_id,
                          file_override = FALSE,
                          create_folder = FALSE)
}
# DONE
#-------------------------------------------------------------------
