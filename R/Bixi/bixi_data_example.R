# the goal of this file is to generate the bixi example for the package.
source("./R/Bixi/helper.R")

split_id <- 1
train_size <- 65
n <- 100; m <- 150

bixi <- bixi_load_split(split_id, train_size)
kernels <- bixi_load_kernels(split_id, train_size, return_distance=TRUE)
bixi_example <- list(Y=bixi$Y, X=bixi$X,Z=bixi$Z,
                     test = as.matrix(bixi$test),
                     spatial_distance = kernels$distance$spatial,
                     temporal_positions = kernels$distance$temporal)



rind <- sample.int(nrow(bixi$Y),n)
cind <- sample.int(ncol(bixi$Y),m)

Bixi_sample <- list( Y = bixi_example$Y[rind,cind],
                     test = bixi_example$test[rind,cind],
                      X = as.matrix(bixi_example$X[rind,]),
                      Z = as.matrix(bixi_example$Z[cind,]),
                      spatial_distance = bixi_example$spatial_distance[cind,cind],
                      temporal_distance = bixi_example$temporal_positions[rind,rind])


usethis::use_data(Bixi_sample, overwrite = TRUE, compress = "xz")
