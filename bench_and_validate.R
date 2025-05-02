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
library(latex2exp) # for LaTeX formatting in plot labels
library(lemon) # for repositioning plot legend in empty facets
library(paletteer) # for colors in benchmark network maps
library(PNWColors) # for colors in benchmark network maps
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
bench_res_dir = "bench_results/"

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

bench_nws = data.frame(
  Network = 1:6,
  nreps_S3N = c(50, 50, 50, 10,  2,  2),
  nreps_SSN = c(50, 50,  5,  1,  1, NA) # for SSN network 5, do preprocessing only, no estimation
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
# S3N: R session takes <=2 GB memory throughout, 
#      obs-obs distances takes 4-5 GB memory
#      preds-obs distances takes 1-2 GB memory per batch, 4 batches (1440 pred points each)
# SSN: R session takes 8-9 GB memory for building LSN, stream updist;
#      6-8 GB for AFV, ~1.3 GB to add obs to LSN, 4-6 GB for estimation


## Subnetwork 5: Two neighboring HUC4s (30748 reaches) --------------------

streams = filter(streams_full, HUC4 %in% c("0514", "0513"))
# benchmark_and_validate(streams, pred_path,
#                        bench_nws$Network[5], 
#                        bench_nws$nreps_S3N[5], 
#                        bench_nws$nreps_SSN[5])
# error in SSN in building LSN when nobs is capped at 1000:
# Error: vector memory limit of 16.0 Gb reached, see mem.maxVSize()
# 32 GB was also not enough
# so set it to 60 GB with mem.maxVSize(60000) and rerun
# indeed, top (in terminal) indicates max memory is > 48 GB
# mem.maxVSize(60000)
benchmark_and_validate(streams, pred_path,
                       bench_nws$Network[5], 
                       bench_nws$nreps_S3N[5], 
                       bench_nws$nreps_SSN[5],
                       SSN_preproc_only = TRUE)


## Subnetwork 6: HUC2 (169092 reaches) ------------------------------------

streams = filter(streams_full, HUC2 == "05")
benchmark_and_validate(streams, pred_path,
                       bench_nws$Network[6],
                       bench_nws$nreps_S3N[6], 
                       bench_nws$nreps_SSN[6],
                       onlyS3N = TRUE)
beep(3)




## Summarize networks -----------------------------------------------------

# for(nw in 1:6){
#   load(paste0("bench_results/network", nw, "_initial_data.rda"))
#   print(paste("Network", nw))
#   print(paste("nreach", nrow(streams), "npreds", nrow(preds), "nobs", nrow(obs)))
#   # print(paste("streams range", range(streams$COMID)))
#   # print(paste("num of negative streams COMIDs", sum(streams$COMID < 0)))
#   # print(paste("preds range", range(preds$COMID)))
#   print("")
# }

# bench_nws = data.frame(
#   Network = 1:3,
#   nreps_S3N = c(50, 50, 50),
#   nreps_SSN = c(50, 50, 5),
#   nreach = c(284, 1273, 7146),
#   npreds = c(283, 1267, 7123),
#   nobs = c(142, 634, 3562)
# )

bench_nws = data.frame(
  Network = 1:6,
  nreach = NA,
  npreds = NA,
  nobs = NA,
  nreps_S3N = c(50, 50, 50, 9, 1, NA),
  nreps_SSN = c(50, 50, 5, 1, NA, NA),
  preproc_only = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE)
)

for(nw in 1:6){
  load(paste0("bench_results/network", nw, "_initial_data.rda"))
  bench_nws$nreach[nw] = nrow(streams)
  bench_nws$npreds[nw] = nrow(preds)
  bench_nws$nobs[nw] = nrow(obs)
}
# Network nreach npreds  nobs nreps_S3N nreps_SSN preproc_only
# 1       1    284    283   142        50        50        FALSE
# 2       2   1273   1267   634        50        50        FALSE
# 3       3   7146   7123  3562        50         5        FALSE
# 4       4  11540  11515  5758        10         2         TRUE
# 5       5  30748  30698 10000        10         2         TRUE

## Combine and summarize benchmarking results -----------------------------

bench1 = combine_runtimes_bothmodels(
  bench_res_dir, 
  bench_nws$Network[1], 
  bench_nws$nreps_S3N[1], 
  bench_nws$nreps_SSN[1])
bench1$obs_only
bench1$obs_preds

bench2 = combine_runtimes_bothmodels(
  bench_res_dir, 
  bench_nws$Network[2], 
  bench_nws$nreps_S3N[2], 
  bench_nws$nreps_SSN[2])
bench2$obs_only
bench2$obs_preds

bench3 = combine_runtimes_bothmodels(
  bench_res_dir,
  bench_nws$Network[3],
  bench_nws$nreps_S3N[3],
  bench_nws$nreps_SSN[3])
bench3$obs_only
bench3$obs_preds

bench4 = combine_runtimes_bothmodels(
  bench_res_dir,
  bench_nws$Network[4],
  bench_nws$nreps_S3N[4],
  bench_nws$nreps_SSN[4],
  obs_only_S3N = TRUE)
bench4$obs_only
bench4$obs_preds

bench5 = combine_runtimes_bothmodels(
  bench_res_dir,
  bench_nws$Network[5],
  bench_nws$nreps_S3N[5],
  bench_nws$nreps_SSN[5],
  obs_only_S3N = TRUE,
  obs_only_SSN = TRUE)
bench5$obs_only
bench5$obs_preds

### compare to old runtimes with original benchmark network choices ---------

load("../BRISC_fish/bench_res/S3N_results_network2_rep1.rda")
print("Network 1"); nrow(streams); nrow(obs); nrow(preds)
runtimes
load("../BRISC_fish/bench_res/SSN_results_network2_rep1.rda"); runtimes

load("../BRISC_fish/bench_res/S3N_results_network3_rep1.rda")
print("Network 2"); nrow(streams); nrow(obs); nrow(preds)
runtimes
load("../BRISC_fish/bench_res/SSN_results_network3_rep1.rda"); runtimes

load("../BRISC_fish/bench_res/S3N_results_network4_rep1.rda")
print("Network 3"); nrow(streams); nrow(obs); nrow(preds)
runtimes
load("../BRISC_fish/bench_res/SSN_results_network4_rep1.rda"); runtimes

## Combine and summarize validation results -------------------------------

valid1 = combine_params_bothmodels(bench_res_dir, 
                                   bench_nws$Network[1],
                                   bench_nws$nreps_S3N[1],
                                   bench_nws$nreps_SSN[1])

valid2 = combine_params_bothmodels(bench_res_dir, 
                                   bench_nws$Network[2],
                                   bench_nws$nreps_S3N[2],
                                   bench_nws$nreps_SSN[2])

valid3 = combine_params_bothmodels(bench_res_dir, 
                                   bench_nws$Network[3],
                                   bench_nws$nreps_S3N[3],
                                   bench_nws$nreps_SSN[3])



## Plot benchmark networks ------------------------------------------------

plot_network = function(bench_res_dir, network, plot_points = FALSE){
  
  my_colors = pnw_palette(
    "Bay", 6, type="continuous")[c(6, 2, 4, 5, 3, 1)]
  
  # # version 1
  # my_colors = pnw_palette(
  #   "Bay", 25, type="continuous")[c(25, 5, 15, 20, 10, 1)]
  
  # my_colors = pnw_palette(
  #   "Bay", 25, type="continuous")[c(seq(5,25,5),
  #                                   seq(4,24,5),
  #                                   seq(3,23,5),
  #                                   seq(2,22,5),
  #                                   seq(1,21,5))]
  
  # load streams, preds and obs for this network
  load(paste0(bench_res_dir, "network", network, "_initial_data.rda"))
  obs = filter(preds, COMID %in% obs$COMID)
  # get additional streams info
  
  load(paste0(data_path, "streams_region5.rda"))
  streams$nw1 = (streams$HUC10 == "0514020607")
  streams$nw2 = (streams$HUC8 == "05140206")
  streams$nw3 = (streams$HUC6 == "051402")
  streams$nw4 = (streams$HUC4 == "0514")
  streams$nw5 = (streams$HUC4 %in% c("0514", "0513"))
  streams$nw6 = (streams$HUC2 == "05")
  # streams$nw2not1 = (streams$nw2 & !streams$nw1)
  # streams$nw3not2 = (streams$nw3 & !streams$nw2)
  # streams$nw4not3 = (streams$nw4 & !streams$nw3)
  # streams$nw5not4 = (streams$nw5 & !streams$nw4)
  # streams$nw6not5 = (streams$nw6 & !streams$nw5)
  # # check: indeed, each reach is in exactly one of these categories except for the two reaches not in Region 5
  # streams$total = streams$nw1 + streams$nw2not1 + streams$nw3not2 + streams$nw4not3 + streams$nw5not4 + streams$nw6not5
  # sum(streams$total != 1) #2
  streams$color = as.factor(case_when(
    streams$nw1 ~ 1,
    streams$nw2 & !streams$nw1 ~ 2,
    streams$nw3 & !streams$nw2 ~ 3,
    streams$nw4 & !streams$nw3 ~ 4,
    streams$nw5 & !streams$nw4 ~ 5,
    streams$nw6 & !streams$nw5 ~ 6
  ))
  
  if(network == 1){
    streams = filter(streams, HUC10 == "0514020607")
    nw_name = "HUC10 0514020607"
  }
  if(network == 2){
    streams = filter(streams, HUC8 == "05140206")
    nw_name = "HUC8 05140206"
  }
  if(network == 3){
    streams = filter(streams, HUC6 == "051402")
    nw_name = "HUC6 051402"
  }
  if(network == 4){
    streams = filter(streams, HUC4 == "0514")
    nw_name = "HUC4 0514"
  }
  if(network == 5){
    streams = filter(streams, HUC4 %in% c("0514", "0513"))
    nw_name = "HUC4s 0514, 0513"
  }
  if(network == 6){
    streams = filter(streams, HUC2 == "05")
    nw_name = "HUC2 05"
    
    usa <- st_as_sf(maps::map("state", fill=TRUE, plot=FALSE)) %>%
      st_transform(crs = st_crs(streams))
    
    # plot of Subnetwork 6 with Subnetwork 5 highlighted in colors
    ggplot() +
      geom_sf(data = usa,
              color = "#2b2b2b",
              fill = "white") +
      geom_sf(data = st_simplify(streams,
                                 preserveTopology = FALSE,
                                 dTolerance = 1000),
              aes(color = color)) +
      coord_sf(xlim = c(-92, -75),
               ylim = c(25, 48),
               default_crs = sf::st_crs(4326)) +
      scale_color_manual(
        values = my_colors)  +
      ggthemes::theme_map() +
      theme(legend.position = "none")
    
    ggsave(filename = paste0(bench_res_dir, "benchmark_nw", network, ".png"),
           width = 4, height = 4, units = "in")
    
    # plot just of Region 5 flowlines
    ggplot() +
      geom_sf(data = usa,
              color = "#2b2b2b",
              fill = "white") +
      geom_sf(data = st_simplify(streams,
                                 preserveTopology = FALSE,
                                 dTolerance = 1000),
              color = "#1874cd") +
      coord_sf(xlim = c(-92, -75),
               ylim = c(25, 48),
               default_crs = sf::st_crs(4326)) +
      ggthemes::theme_map() +
      theme(legend.position = "none")
    
    ggsave(filename = paste0(bench_res_dir, "flowlines_Region5.png"),
           width = 4, height = 4, units = "in")
    
  } else { # not all of region 5, so small enough to not need st_simplify or the map
    
    if(plot_points){
      ggplot() +
        geom_sf(data = streams, aes(color = color)) +
        geom_sf(data = preds, color = "purple", size = 0.5) +
        geom_sf(data = obs, color = "gold", size = 0.5) +
        coord_sf(datum = st_crs(streams)) +
        labs(title = paste0("Subnetwork ", network, ": ", nw_name)) +
        scale_color_manual(
          values = my_colors)  +
        theme_bw() +
        theme(legend.position = "none")
    } else {
      ggplot() +
        geom_sf(data = streams, aes(color = color)) +
        labs(title = paste0("Subnetwork ", network, ": ", nw_name)) +
        scale_color_manual(
          values = my_colors)  +
        theme_bw() +
        theme(legend.position = "none")
    }
    
    ggsave(filename = paste0(bench_res_dir, "benchmark_nw", network, ".png"),
           width = 4, height = 4, units = "in")
    
  }
  
}

# Subnetwork 1
plot_network(bench_res_dir, network = 1, plot_points = FALSE)

# Subnetwork 2
plot_network(bench_res_dir, network = 2, plot_points = FALSE)

# Subnetwork 3
plot_network(bench_res_dir, network = 3, plot_points = FALSE)

# Subnetwork 4
plot_network(bench_res_dir, network = 4, plot_points = FALSE)

# Subnetwork 5
plot_network(bench_res_dir, network = 5, plot_points = FALSE)

# Subnetwork 6
plot_network(bench_res_dir, network = 6, plot_points = FALSE)


## Generate more benchmark and validation results for large networks --------------

# for benchmarking only, no estimation/validation
# load streams
load(paste0(data_path, "streams_region5.rda"))
# create a backup copy
streams_full = streams
# set lsn file path and output path for benchmark/validation results
lsn.path = "lsn"

bench_nws = data.frame(
  Network = 1:6,
  nreps_S3N = c(50, 50, 50, 10, 10, 10),
  nreps_SSN = c(50, 50, 50, 10,  2, NA) # for SSN network 5, do preprocessing only, no estimation
)

## network 4
streams = filter(streams_full, HUC4 == "0514")
benchmark_and_validate(streams, pred_path,
                       bench_nws$Network[4], 
                       bench_nws$nreps_S3N[4], 
                       bench_nws$nreps_SSN[4],
                       predist_only = TRUE)

## network 6
streams = filter(streams_full, HUC2 == "05")
benchmark_and_validate(streams, pred_path,
                       bench_nws$Network[6], 
                       bench_nws$nreps_S3N[6], 
                       bench_nws$nreps_SSN[6],
                       predist_only = TRUE,
                       onlyS3N = TRUE)

## network 5
streams = filter(streams_full, HUC4 %in% c("0514", "0513"))
mem.maxVSize(60000)
benchmark_and_validate(streams, pred_path,
                       bench_nws$Network[5], 
                       bench_nws$nreps_S3N[5], 
                       bench_nws$nreps_SSN[5],
                       predist_only = TRUE)
mem.maxVSize(16384)
Sys.sleep(30)
beep(3)

# Generate more SSN validation results for Network 3 (slow, takes about 30 hours)
nreps = 50; network = 3; out_dir = "valid_results/"
lsn.path = "lsn"
ndigits = ceiling(log(nreps, base = 10))
for(rep in 1:nreps){
  message(paste("Network", network, "SSN benchmark rep", rep, "of", nreps))
  SSN_preproc_and_estimation(network, 
                             str_pad(rep, ndigits, side = "left", pad = "0"),
                             # need estimation for validation results
                             out_dir, lsn.path, preproc_only = FALSE, predist_only = FALSE)
}
