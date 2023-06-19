library(data.table)
library(spacexr)
library(Matrix)
library(dplyr)
library(tidyr)
library(tibble)
library(spASE)

c57 <- data.frame(fread('smartseq3/SS3_c57_UMIs_concat.csv'))
rownames(c57) <- c57$V1
c57$V1 <- NULL
cast <- data.frame(fread('smartseq3/SS3_cast_UMIs_concat.csv'))
rownames(cast) <- cast$V1
cast$V1 <- NULL

c57 <- Matrix(as.matrix(c57))
cast <- Matrix(as.matrix(cast))
tot <- c57+cast



myfit <- scase(c57, cast, min.cells = 2^7, cores=8, verbose=T)

saveRDS(myfit, 'betabin_results_smartseq3.rds')

