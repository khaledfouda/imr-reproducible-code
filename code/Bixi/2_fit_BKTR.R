source("./code/helper.R")
source("./code/Bixi/helper.R")


# two steps. 1. fit all training files (250!) and 2. generate the results tables.
# expected time is 22 minutes * 250 files ~ 4 days!
#-------------------------------------------------------
# 1. fit all 250 files
for (split_id in 1:NUM_SPLITS) {
  for (train_size in TRAIN_SEQ) {
    bktr_fit <- bixi_fit_bktr(
      split_id, train_size,
      burn_in_iter = BKTR_ITER_BURN,
      sampling_iter = BKTR_ITER
    )
  }
}
#-------------------------------------------------------------
# 2. generate the results table and save to /output/

all_results <- data.frame()
for (split_id in 1:NUM_SPLITS) {
  for (train_size in TRAIN_SEQ) {
    all_results <- rbind(
      all_results,
      bixi_fit_bktr_post(split_id, train_size) %>%
        mutate(
          split_id = split_id,
          train_size = train_size
        )
    )
  }
  print(paste0("split #", split_id))
}

# save to disk >> 
rw_a_file(
  "results_bktr.rds",
  data = all_results,
  file_override = TRUE,
  create_folder = TRUE,
  directory = "./output/Bixi/",
  type = "write"
)

# checking the results > 
all_results %>%
  select(-model) %>%
  group_by(train_size) %>%
  summarize_all(mean) %>%
  ungroup() %>%
  as.data.frame()

#--- DONE

