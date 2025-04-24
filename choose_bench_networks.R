# choosing networks to use for S3N benchmarking and validation
# PhD dissertation, Chapter 3
# Jess Kunke, 2025 April

# run this script before bench_and_validate.R

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


# preprocess streams data before identifying benchmark subnetworks --------

## envir covars preprocessing ----------------------------------------------

# only need to do this once unless the covariates change; uncomment if needed

# # combine environmental covariates and write to csv file
# # note default input: comid_huc12_filename = "HUC12_PU_COMIDs_CONUS.csv"
# combine_and_write_national_covariate_data(covs_path, comid_huc12_path)

## streams data preprocessing ----------------------------------------------

# read in streams flowlines for a given region,
# join the combined envir covariate data,
# reassign duplicate COMIDs as negative COMIDs,
# impute and rename covariates,
# write imputed renamed covariates to file
tic("streams data preprocessing")
streams = read_sf(streams_path, "Flowline_MS05_NSI") %!>% # 169463 obs of 16 vars
  # must join envir covar before reassigning dup COMIDs or else the dup COMIDs won't get HUC codes or envir covars
  join_envir_covar_data(covs_path) %!>%
  reassign_dup_COMIDs() %!>%
  impute_and_rename_covars()
toc()

## build the stream network -----------------------------------------------

tic("configure")
streams_res = configure_stream_network(streams)
toc() # 56.729 sec elapsed (55-65 sec)
streams = streams_res$streams
stream_graphs = streams_res$sg
# now streams has 169463 obs of 36 vars
# outflow:
#     0      1 
#    63 169400
# inflow:
#     0     1     2     3 
# 64095 41337 64030     1
# This network has 63 separate components.

## explore the complex confluence -----------------------------------------
# complex_conf = names(which(table(streams$dnstream_node) == 3))
# # "1151540.2996 2059592.6663"
# streams$complex = (streams$dnstream_node == complex_conf)
# mapview(filter(streams, HUC8 == "05040002"), zcol = c("complex", "COMID"))

## remove smallest branch of complex confluence ---------------------------
# - COMID 167484513 and all reaches upstream of it
complex_conf = 167484513
upnd1 = streams$upstream_node[streams$COMID == complex_conf]
dnnd1 = streams$COMID[streams$dnstream_node == upnd1] # 15411283
upnd2 = streams$upstream_node[streams$COMID == dnnd1]
dnnd2 = streams$COMID[streams$dnstream_node == upnd2] # 15409971 15410029
upnd3 = streams$upstream_node[streams$COMID %in% dnnd2]
dnnd3 = streams$COMID[streams$dnstream_node %in% upnd3] # 15409975 15409973
upnd4 = streams$upstream_node[streams$COMID %in% dnnd3]
dnnd4 = streams$COMID[streams$dnstream_node %in% upnd4] # 15409981 15411355
upnd5 = streams$upstream_node[streams$COMID %in% dnnd4]
dnnd5 = streams$COMID[streams$dnstream_node %in% upnd5] # 15411353 15411309
upnd6 = streams$upstream_node[streams$COMID %in% dnnd5]
dnnd6 = streams$COMID[streams$dnstream_node %in% upnd6] # empty

COMIDs_to_drop = c(complex_conf, 
                   dnnd1, 
                   dnnd2, 
                   dnnd3, 
                   dnnd4, 
                   dnnd5, 
                   dnnd6)

# remove these 10 COMIDs and also any COMIDs outside of the main component
streams2 = streams %>%
  filter(!(COMID %in% COMIDs_to_drop), componentID == 1)
# now streams has 169094 obs (369 fewer COMIDs)
100*nrow(streams2)/nrow(streams) # we are still keeping 99.78% of the network
streams = streams2; rm(streams2)

## build the stream network again -----------------------------------------

tic()
streams_res = configure_stream_network(streams)
toc() # 63.004 sec elapsed (55-65 sec)
streams = streams_res$streams
stream_graphs = streams_res$sg
rm(streams_res)
# still 169094 obs of 36 vars
# outflow:
#     0      1 
#     1 169093 
# inflow:
#     0      1      2 
# 63905  41285  63904 
# This network has 1 separate component.


## save stream network ----------------------------------------------------

# save region 5 up to this point
save(streams, stream_graphs, file = paste0(data_path, "streams_region5_allvars.rda"))

# also save with fewer variables
load("streams_region5_allvars.rda")
streams = select(streams, COMID, HUC8, HUC12, LENGTHKM, TotDASqKM, 
                 Development, Agriculture, Elevation, Water_Area, 
                 Ann_Runoff, Baseflow, Ann_Temp, Hydro_Alter, 
                 Fldplain_Dis, geometry)
save(streams, file = paste0(data_path, "streams_region5_fewvars.rda"))
rm(stream_graphs)


# identify stream subnetworks for benchmarking and validation ----------------

load(paste0(data_path, "streams_region5_fewvars.rda"))

# configure the stream network
streams_res = configure_stream_network(streams)
streams = streams_res$streams
stream_graphs = streams_res$sg
stream_graphs$conn_comp$no # should be 1

# locate the stream outlet (sink)
which(streams$is_sink) # 40332
streams$COMID[which(streams$is_sink)] # 1844789
streams$HUC12[which(streams$is_sink)] # 080101000103
outlet = streams$COMID[streams$is_sink] # 1844789
upnd1 = streams$upstream_node[streams$COMID == outlet]
(dnnd1 = streams$COMID[streams$dnstream_node == upnd1]) # 1842647
outlet = streams$HUC12[streams$COMID == dnnd1] # 051402060704

streams$HUC10 = str_sub(streams$HUC12, end = 10)
streams$HUC6 = str_sub(streams$HUC12, end = 6)
streams$HUC4 = str_sub(streams$HUC12, end = 4)
streams$HUC2 = str_sub(streams$HUC12, end = 2)

# all of region 5 has HUC2 = "05" as expected, except for one source reach and 
# one sink reach
unique(streams$HUC2)
# [1] "05" "08" "04"
streams$is_src[streams$HUC2 == "04"]
# [1] TRUE
streams$is_sink[streams$HUC2 == "08"]
# [1] TRUE

save(streams, file = paste0(data_path, "streams_region5.rda"))
streams_full = streams

# let's identify a nested set of networks for benchmarking and validation

# start with HUC12 containing the stream outlet
streams = filter(streams_full, HUC12 == "051402060704")
(nreach = nrow(streams)) # 42 (too small)
# Subnetwork 1: HUC 10 containing stream outlet
streams = filter(streams_full, HUC10 == "0514020607")
(nreach = nrow(streams)) # 284
# Subnetwork 2: HUC 10 containing stream outlet
streams = filter(streams_full, HUC8 == "05140206")
(nreach = nrow(streams)) # 1273
# Subnetwork 3: HUC 10 containing stream outlet
streams = filter(streams_full, HUC6 == "051402")
(nreach = nrow(streams)) # 7146
# Subnetwork 4: HUC 10 containing stream outlet
streams = filter(streams_full, HUC4 == "0514")
(nreach = nrow(streams)) # 11540
# see what HUC4s are upstream
upnd1 = unique(streams$upstream_node)
(dnnd1 = unique(streams$HUC4[streams$dnstream_node %in% upnd1]))
# "0514" "0511" "0513" "0512" "0509" "0510"
# Subnetwork 5: two neighboring HUC4 regions, one of which contains stream outlet
streams = filter(streams_full, HUC4 %in% c("0514", "0513"))
(nreach = nrow(streams)) # 30748 (size close to original subnetwork 5)
# Subnetwork 6: HUC2 containing stream outlet (essentially Region 5)
streams = filter(streams_full, HUC2 == "05")
(nreach = nrow(streams)) # 169092

