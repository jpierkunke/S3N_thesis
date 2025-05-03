# Region 5 fish population estimates with S3N
# PhD dissertation, Chapter 4
# Jess Kunke, 2025 April

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

# fish survey data
fish_path = "/Users/jessicakunke/Documents/UW/Research/Julian-Olden-SAFS/Data_Region5/fishscales/fishdata_2024-05-02/"

# default location for where to save preprocessing output from this script
preproc_output = "pwdists/input/"
# where to save data to use as input in computing pairwise distances
pwdist_input_dir = "pwdists/input/"
# where to find the preds-obs pairwise distance results from the cluster
pwdist_predsobs_dir = "pwdists/output_predsobs/"
# where to find the obs-obs pairwise distance results from the cluster
pwdist_obsobs_dir = "pwdists/output_obsobs/"

# preprocess streams data before identifying benchmark subnetworks --------

# this code is also part of choose_bench_networks.R but is included here for 
# the sake of completeness

# envir covars preprocessing ----------------------------------------------

# only need to do this once unless the covariates change; uncomment if needed

# # combine environmental covariates and write to csv file
# # note default input: comid_huc12_filename = "HUC12_PU_COMIDs_CONUS.csv"
# combine_and_write_national_covariate_data(covs_path, comid_huc12_path)

# streams data preprocessing ----------------------------------------------

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
toc() # streams data preprocessing: 13.156 sec elapsed

# build the stream network -----------------------------------------------

tic("configure stream network")
streams_res = configure_stream_network(streams)
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
complex_conf = names(which(table(streams$dnstream_node) == 3))
# "1151540.2996 2059592.6663"
streams$complex = (streams$dnstream_node == complex_conf)
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

streams_res = configure_stream_network(streams)
streams = streams_res$streams
stream_graphs = streams_res$sg
rm(streams_res)
toc() # configure stream network: 124.363 sec elapsed
# still 169094 obs of 36 vars
# outflow:
#     0      1 
#     1 169093 
# inflow:
#     0      1      2 
# 63905  41285  63904 
# This network has 1 component.

# Compute stream reach upstream distance and AFV --------------------------

tic("Stream updist and AFV")
streams = compute_stream_updist_vars(streams, stream_graphs)
toc() # Stream updist and AFV: 78.175 sec elapsed
# ~1280 iterations

# save(streams, file = paste0(data_path, "streams_for_Region5_results.rda"))

# Add observation points to LSN, compute site updist and AFV --------------

# load(paste0(data_path, "streams_for_Region5_results.rda"))

tic("Add obs to LSN")
# load preprocessed fish data
# read_fish_data(fish_path, with_cutoff_year = TRUE) # preprocessing took 281.165 sec
load(paste0(fish_path, "combined_fish_data_starting_1990.RData"))
# prep_to_compute_pwdist_region5 accomplished several tasks:
# - reads in and preprocess prediction points layer, add updist and AFV
# - joins geometry, updist and AFV from pred points layer to fish dataset and
#   selects only unique points to create the obs points layer
prep_to_compute_pwdist_region5(streams, fish, pred_path,
                               out_dir = pwdist_input_dir)
toc() # Add obs to LSN: 21.916 sec elapsed

# Compute obs-obs pairwise distances and neighbor variables ---------------

# run this part of the code on the cluster

# obs-obs distances are used in both estimation and prediction
# neighbor variables are just for estimation

# tic("Obs-obs distances")
# system('cd pwdists; ./scripts/all_obsobs.sh')
# toc() # Obs-obs distances: sec elapsed
# 
# # Compute preds-obs pairwise distances and neighbor variables for prediction --
# 
# tic("Preds-obs distances")
# # determine batch number and size based on number of (unobserved) prediction points
# batch_info = get_preds_batch_info( nrow(preds) - nrow(obs) )
# m=10 # use 10 neighbors
# args = sprintf("%i %i %i", batch_info$n, batch_info$size, m)
# system(paste('cd pwdists; ./scripts/run_pwdist_locally.sh', args))
# system(paste('cd pwdists; ./scripts/combine_pred_neighbors.sh', args))
# toc() # Preds-obs distances: sec elapsed

# Prepare for estimation and prediction -----------------------------------

# do this once for Region 5, not once per species

# add covariates to preds and obs and save to file
load(paste0(pwdist_input_dir, "preds_obs_pwdist_input_data_Region5.RData"))
# names(preds); names(obs) # they're the same
# # [1] "ptID"      "COMID"     "networkID" "binaryID"  "up_dist"   "log_AFV"   "geometry" 
# names(streams) # 44 variables

# add covariates to preds
preds = left_join(
  preds,
  streams %>%
    st_drop_geometry() %>%
    data.frame() %>%
    select(ptID, Development, Agriculture, Elevation, Water_Area, 
           Ann_Runoff, Baseflow, Ann_Temp, Hydro_Alter, Fldplain_Dis),
  by="ptID"
)

# add covariates to obs
obs = left_join(
  obs,
  streams %>%
    st_drop_geometry() %>%
    data.frame() %>%
    select(ptID, Development, Agriculture, Elevation, Water_Area, 
           Ann_Runoff, Baseflow, Ann_Temp, Hydro_Alter, Fldplain_Dis),
  by="ptID"
)

# add to fish
load(paste0(fish_path, "combined_fish_data_starting_1990.RData"))
fish_region5 = fish %>%
  filter(COMID %in% preds$COMID) %>%
  left_join(streams %>% st_drop_geometry() %>% select(COMID, LENGTHKM)) %>%
  compute_fish_reach_count_from_density()

# table(fish_region5$Source)
# FishScales GiamOlden16 
# 164645        7533 

save(preds, obs, fish_region5, 
     file = paste0(pwdist_input_dir, "preds_obs_fish_for_estimation.RData"))


# # this is the old code I used for the first Region 5 results
# count_each_species_in_fish(
#   fish, out_dir = fish_path, 
#   filename = "nobs_by_species_Region05.csv"
# )
# # 309 species
# # nobs_by_species %>% filter(nCOMIDs > 15) # 202 species left
# 
# # how about we make predictions for all species with (non-zero) obs in at least 
# # as many COMIDs as the number of neighbors we're using?
# # in a sense, they all have obs at all the COMIDs, since the surveyers count 
# # all the fish they find and fish that are not recorded were not observed at 
# # that time and place, but I mean non-zero observations
# # species_for_prediction = filter(nobs_by_species, nCOMIDs > 15)$Common_Name
# 
# get_obs_layer_for_each_species(
#   fish, 
#   out_dir = paste0(fish_path, "individual_species/"), 
#   region = "Region05"
# )




# Estimation and prediction for each species ------------------------------------

pwdist_obsobs_dir = "Region5_cluster_results/"
pwdist_predsobs_dir = "Region5_cluster_results/"

# needed for estimation
# loads obsobs_dist, obsobs_wt
load(paste0(pwdist_obsobs_dir, "obsobs_dist_wt.rda"))
# loads obs_neighbors
load(paste0(pwdist_obsobs_dir, "obs_neighbors.rda"))
# loads preds, obs, fish_region5
load(paste0(pwdist_input_dir, "preds_obs_pwdist_input_data_Region5.RData"))
load(paste0(pwdist_input_dir, "preds_obs_fish_for_estimation.RData"))

common_name = "Creek Chub"
this_species = fish_region5 %>%
  filter(Common_Name == common_name) %>%
  select(COMID, DensityPer100m)

if(nrow(this_species) != n_distinct(this_species$COMID)){
  stop("COMID is not a unique identifier for this_species")
}

obs_this_species = obs %>%
  left_join(this_species, by = "COMID") %>%
  # set to 0 if NA because this species was not observed there
  mutate(DensityPer100m = replace_na(DensityPer100m, 0))

tic("Estimation"); start_time = Sys.time()
estimation = BRISC_estimation_stream(
  coords = as.matrix(1:nrow(obs_this_species)),
  y = as.matrix(obs_this_species$DensityPer100m),
  x = get_X(obs_this_species),
  neighbor = obs_neighbors,
  cov.model = "exponential",
  nugget_status = 1,
  verbose = TRUE
)
toc(); stop_time = Sys.time()
# runtimes[5] = as.numeric(difftime(stop_time, start_time), units="secs")

# params = c(
#   estimation$Beta, # beta (intercept, elevation)
#   estimation$Theta # sigma.sq, tau.sq, phi
# )
# params[5] = 1/params[5] # convert phi to lambda (range parameter)
# names(params)[5] = "lambda"

# save estimated parameter values, n_neighbors, species, and runtime

# prediction
# load pred_neighbors object
load(paste0(pwdist_predsobs_dir, "pred_neighbors.rda"))
prediction = BRISC_prediction_stream(
  BRISC_Out = estimation,
  coords.0 = as.matrix(1:nrow(preds)),
  X.0 = get_X(preds),
  # X.0 = streams %>%
  #       st_drop_geometry() %>%
  #       select(Elevation) %>%
  #       # select(Elevation, Water_Area,
  #       #        Ann_Runoff, Baseflow, Ann_Temp,
  #       #        Development, Agriculture,
  #       #        Hydro_Alter, Fldplain_Dis) %>%
  #       mutate(Intercept = 1, .before = Elevation) %>%
  #       as.matrix(),
  neighbor = pred_neighbors,
  obsobs_dist = obsobs_dist,
  obsobs_wt = obsobs_wt,
  verbose = TRUE
)
toc(); stop_time = Sys.time()
# runtimes[7] = as.numeric(difftime(stop_time, start_time), units="secs")



