# benchmarking and validation
# PhD dissertation, Chapter 3
# Jess Kunke, 2025 April

# run choose_bench_networks.R before running this script

library(sf) # for reading and possibly writing shapefiles
library(plyr) # for rbind.fill
library(knitr) # for making tables
library(beepr) # for notifying when a task is done
library(igraph) # for making the adjacency matrix
library(tictoc) # for benchmarking (timing to evaluate efficiency/runtime)
library(Matrix) # for efficient sparse adjacency matrix operations
library(mapview) # for visualizing spatial data
library(tidyverse) # for efficient data manipulation
library(shapefiles) # for read.dbf
library(magrittr) # for the eager pipe
library(SSNbler) # to benchmark SSN preprocessing
library(SSN2) # to benchmark SSN distances and estimation
library(BRISC) # version of BRISC I modified to handle S3N estimation and prediction

source("S3N_functions.R")

# set file paths
data_path = "input_data/"
covs_path = paste0(data_path, "envir_covs/") # environmental covariate data
streams_path = paste0(data_path, "flowlines/")  # flowline layer for Region 5
pred_path = paste0(data_path, "predpoints/") # prediction point layer for Region 5
comid_huc12_path = paste0(data_path, "match_COMID_to_HUC12/") # data matching COMIDs to HUC12s

# default location for where to save preprocessing output from this script
preproc_output = "pwdists/input/"
# where to save data to use as input in computing pairwise distances
pwdist_input_dir = "pwdists/input/"
# where to find the preds-obs pairwise distance results from the cluster
pwdist_predsobs_dir = "pwdists/output_predsobs/"
# where to find the obs-obs pairwise distance results from the cluster
pwdist_obsobs_dir = "pwdists/output_obsobs/"

## preprocess before benchmarking and validation ---------------------------

# before the benchmarking and validation, start a fresh session
# load the libraries, source S3N_functions.R, and set file paths (lines 1-35 of this script)
# then skip to here

# load streams
load(paste0(data_path, "streams_region5.rda"))
# create a backup copy
streams_full = streams
# set lsn file path and output path for benchmark/validation results
lsn.path = "lsn"
bench_res_dir = "bench_res/"

bench_nws = data.frame(
  Network = 1:6,
  nreps_S3N = c(50, 50, 50, 10, 10, 10),
  nreps_SSN = c(50, 50, 50, 10, 2, NA)
)


## Subnetwork 1: HUC10 containing stream outlet (284 reaches) -------------

streams = filter(streams_full, HUC10 == "0514020607")
benchmark_and_validate(streams, pred_path,
                       bench_nws$Network[1], 
                       bench_nws$nreps_S3N[1], 
                       bench_nws$nreps_SSN[1])


## Subnetwork 2: HUC8 containing stream outlet (1273 reaches) -------------

streams = filter(streams_full, HUC8 == "05140206")
benchmark_and_validate(streams, pred_path,
                       bench_nws$Network[2], 
                       bench_nws$nreps_S3N[2], 
                       bench_nws$nreps_SSN[2])


## Subnetwork 3: HUC6 containing stream outlet (7146 reaches) -------------

streams = filter(streams_full, HUC6 == "051402")
benchmark_and_validate(streams, pred_path,
                       bench_nws$Network[3], 
                       bench_nws$nreps_S3N[3], 
                       bench_nws$nreps_SSN[3])

## Subnetwork 4: HUC4 containing stream outlet (11540 reaches)-------------

streams = filter(streams_full, HUC4 == "0514")
benchmark_and_validate(streams, pred_path,
                       bench_nws$Network[4], 
                       bench_nws$nreps_S3N[4], 
                       bench_nws$nreps_SSN[4])


## Subnetwork 5: Two neighboring HUC4s (30748 reaches) --------------------

streams = filter(streams_full, HUC4 %in% c("0514", "0513"))
# benchmark_and_validate(streams, pred_path,
#                        bench_nws$Network[5], 
#                        bench_nws$nreps_S3N[5], 
#                        bench_nws$nreps_SSN[5])
# error in SSN in building LSN:
# Error: vector memory limit of 16.0 Gb reached, see mem.maxVSize()
# 32 GB was also not enough
# so set it to 60 GB with mem.maxVSize(60000) and rerun
# indeed, top (in terminal) indicates max memory is > 48 GB
mem.maxVSize(60000)
benchmark_and_validate(streams, pred_path,
                       bench_nws$Network[5], 
                       bench_nws$nreps_S3N[5], 
                       bench_nws$nreps_SSN[5])


## Subnetwork 6: HUC2 (169092 reaches) ------------------------------------

streams = filter(streams_full, HUC2 == "05")
benchmark_and_validate(streams, pred_path,
                       bench_nws$Network[6],
                       bench_nws$nreps_S3N[6], 
                       bench_nws$nreps_SSN[6], 
                       onlyS3N = TRUE)
beep(3)


## Summarize networks -----------------------------------------------------

for(nw in 1:6){
  load(paste0("bench_results/network", nw, "_initial_data.rda"))
  print(paste("Network", nw))
  print(paste("nreach", nrow(streams), "npreds", nrow(preds), "nobs", nrow(obs)))
  # print(paste("streams range", range(streams$COMID)))
  # print(paste("num of negative streams COMIDs", sum(streams$COMID < 0)))
  # print(paste("preds range", range(preds$COMID)))
  print("")
}

## Combine and summarize benchmarking results -----------------------------

## Combine and summarize validation results -------------------------------






