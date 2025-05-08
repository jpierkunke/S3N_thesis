# benchmarking and validation
# PhD dissertation, Chapter 3
# Jess Kunke, 2025 May

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
  nreps_S3N = c(50, 50, 50, 10, 10, 10),
  nreps_SSN = c(50, 50, 50, 10,  2, NA),
  preproc_only = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE)
)


## Subnetwork 1: HUC10 containing stream outlet (284 reaches) -------------

streams = filter(streams_full, HUC10 == "0514020607")
benchmark_and_validate(streams, pred_path,
                       bench_nws$Network[1], 
                       bench_nws$nreps_S3N[1], 
                       bench_nws$nreps_SSN[1],
                       bench_nws$preproc_only[1])


## Subnetwork 2: HUC8 containing stream outlet (1273 reaches) -------------

streams = filter(streams_full, HUC8 == "05140206")
benchmark_and_validate(streams, pred_path,
                       bench_nws$Network[2], 
                       bench_nws$nreps_S3N[2], 
                       bench_nws$nreps_SSN[2],
                       bench_nws$preproc_only[2])


## Subnetwork 3: HUC6 containing stream outlet (7146 reaches) -------------

streams = filter(streams_full, HUC6 == "051402")
benchmark_and_validate(streams, pred_path,
                       bench_nws$Network[3], 
                       bench_nws$nreps_S3N[3], 
                       bench_nws$nreps_SSN[3],
                       bench_nws$preproc_only[3])

## Subnetwork 4: HUC4 containing stream outlet (11540 reaches)-------------

streams = filter(streams_full, HUC4 == "0514")
benchmark_and_validate(streams, pred_path,
                       bench_nws$Network[4], 
                       bench_nws$nreps_S3N[4], 
                       bench_nws$nreps_SSN[4],
                       bench_nws$preproc_only[4])
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
                       bench_nws$preproc_only[5])


## Subnetwork 6: HUC2 (169092 reaches) ------------------------------------

streams = filter(streams_full, HUC2 == "05")
benchmark_and_validate(streams, pred_path,
                       bench_nws$Network[6],
                       bench_nws$nreps_S3N[6], 
                       bench_nws$nreps_SSN[6],
                       bench_nws$preproc_only[6],
                       onlyS3N = TRUE)
beep(3)




## Summarize networks -----------------------------------------------------

bench_nws = data.frame(
  Network = 1:6,
  nreach = NA,
  npreds = NA,
  nobs = NA,
  nreps_S3N = c(50, 50, 50, 10, 10, 10),
  nreps_SSN = c(50, 50, 50, 10,  2, NA), # eventually do 50 SSN runs on network 3
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
# 3       3   7146   7123  3562        50        50        FALSE
# 4       4  11540  11515  5758        10        10         TRUE
# 5       5  30748  30698 10000        10         2         TRUE
# 6       6 169092 168875 10000        10        NA         TRUE

## Combine and summarize benchmarking results -----------------------------

# network 1: S3N 6 runtimes + params, SSN 15 runtimes + params (estimation for both, obs + preds both)
# network 2: S3N 6 runtimes + params, SSN 15 runtimes + params (estimation for both, obs + preds both)
# network 3: S3N 6 runtimes + params, SSN 15 runtimes + params (estimation for both, obs + preds both)
# network 4: S3N 5 runtimes (last 2 NA), SSN 5 runtimes (last 1 NA) (predistance only + obs only for both)
# network 5: S3N 5 runtimes (last 2 NA), SSN 5 runtimes (last 1 NA) (predistance only + obs only for both)
# network 6: S3N 5 runtimes (last 2 NA), no SSN (predistance only + obs only for S3N)

get_summary_for_table = function(network, benchres){
  if(network < 4){
    nwres = data.frame(
      Network = network,
      Task = benchres$obs_only$task[1:3],
      S3N = benchres$obs_only$S3N_avg[c(1,2,4)],
      SSN = c(benchres$obs_only$SSN_avg[1:2],
              sum(benchres$obs_only$SSN_avg[3:4]))
    )
  } else {
    nwres = data.frame(
      Network = network,
      Task = benchres$obs_only$task[1:3],
      S3N = benchres$obs_only$S3N_avg[1:3],
      SSN = benchres$obs_only$SSN_avg[1:3]
    )
  }
  return(nwres)
}

bench1 = combine_runtimes_bothmodels(
  bench_res_dir, 
  bench_nws$Network[1], 
  bench_nws$nreps_S3N[1], 
  bench_nws$nreps_SSN[1])
bench1$obs_only
bench1$obs_preds
#                    task    S3N_avg      S3N_sd    SSN_avg    SSN_sd
# 1             Build LSN 0.22744269 0.138381518  1.5055646 0.8733712
# 2 Stream updist and AFV 1.81468081 0.209169285  1.1515321 1.0070120
# 3        Add obs to LSN         NA          NA  6.6005514 0.7398649
# 4    Obs updist and AFV 3.64795209 0.195317421  0.1410471 0.5748062
# 5          Assemble SSN         NA          NA  1.4945106 0.8552582
# 6           Obs distmat 2.74083258 0.117723897  0.5760140 0.5998521
# 7            Estimation 0.02619061 0.004516947  1.6473641 0.9296746
# 8                 Total 8.45709879          NA 13.1165839        NA
#
#                    task    S3N_avg      S3N_sd    SSN_avg    SSN_sd
# 1             Build LSN 0.22744269 0.138381518  1.5055646 0.8733712
# 2 Stream updist and AFV 1.81468081 0.209169285  1.1515321 1.0070120
# 3        Add obs to LSN 3.64795209          NA   6.741598        NA
# 5          Assemble SSN         NA          NA  1.4945106 0.8552582
# 6           Obs distmat 2.74083258 0.117723897  0.5760140 0.5998521
# 7            Estimation 0.02619061 0.004516947  1.6473641 0.9296746
# 8                 Total 8.45709879          NA 13.1165839        NA

bench2 = combine_runtimes_bothmodels(
  bench_res_dir, 
  bench_nws$Network[2], 
  bench_nws$nreps_S3N[2], 
  bench_nws$nreps_SSN[2])
bench2$obs_only
bench2$obs_preds
#                    task    S3N_avg     S3N_sd     SSN_avg      SSN_sd
# 1             Build LSN  0.7285639 0.08947661  2.99213656 0.332678843
# 2 Stream updist and AFV  2.3312726 0.11925206  4.72642903 0.376432657
# 3        Add obs to LSN         NA         NA 25.18034643 0.314592672
# 4    Obs updist and AFV  3.8286596 0.31508041  0.06358808 0.009856377
# 5          Assemble SSN         NA         NA  2.32831927 0.044189002
# 6           Obs distmat  5.2026225 0.15863773  0.57475306 0.023485215
# 7            Estimation  0.1761315 0.06651270  8.47340318 1.239037395
# 8                 Total 12.2672500         NA 44.33897562          NA
# 
#                    task    S3N_avg     S3N_sd     SSN_avg      SSN_sd
# 1             Build LSN  0.7285639 0.08947661  2.99213656 0.332678843
# 2 Stream updist and AFV  2.3312726 0.11925206  4.72642903 0.376432657
# 3        Add obs to LSN  3.8286596         NA 25.24393             NA
# 5          Assemble SSN         NA         NA  2.32831927 0.044189002
# 6           Obs distmat  5.2026225 0.15863773  0.57475306 0.023485215
# 7            Estimation  0.1761315 0.06651270  8.47340318 1.239037395
# 8                 Total 12.2672500         NA 44.33897562          NA

bench3 = combine_runtimes_bothmodels(
  bench_res_dir,
  bench_nws$Network[3],
  bench_nws$nreps_S3N[3],
  bench_nws$nreps_SSN[3],
  obs_only_S3N = FALSE,
  obs_only_SSN = TRUE)
bench3$obs_only
bench3$obs_preds
#                    task     S3N_avg    S3N_sd      SSN_avg       SSN_sd
# 1             Build LSN   2.7910486 0.4181868   26.8652709   0.91765570
# 2 Stream updist and AFV   6.6825090 0.3104472   66.0792157   2.96400747
# 3        Add obs to LSN          NA        NA  142.2515071   0.82120930
# 4    Obs updist and AFV   4.0823088 0.3371444    0.1851079   0.01100262
# 5          Assemble SSN          NA        NA   18.9790649   1.78823961
# 6           Obs distmat 106.6794731 0.9506971   14.4666812   0.31992008
# 7            Estimation   0.8290041 0.4319330  912.4498494 107.34362815
# 8                 Total 121.0643436        NA 1181.2766972           NA
#
#                    task     S3N_avg    S3N_sd      SSN_avg       SSN_sd
# 1             Build LSN   2.7910486 0.4181868   26.8652709   0.91765570
# 2 Stream updist and AFV   6.6825090 0.3104472   66.0792157   2.96400747
# 3        Add obs to LSN   4.0823088        NA  142.4366              NA
# 5          Assemble SSN          NA        NA   18.9790649   1.78823961
# 6           Obs distmat 106.6794731 0.9506971   14.4666812   0.31992008
# 7            Estimation   0.8290041 0.4319330  912.4498494 107.34362815
# 8                 Total 121.0643436        NA 1181.2766972           NA

bench4 = combine_runtimes_predistonly(
  bench_res_dir,
  bench_nws$Network[4],
  bench_nws$nreps_S3N[4],
  bench_nws$nreps_SSN[4],
  obs_only_S3N = TRUE,
  obs_only_SSN = TRUE,
  predist_only = TRUE)
# predist_only_S3N = TRUE,
# predist_only_SSN = TRUE)
bench4$obs_only
#                    task   S3N_avg    S3N_sd   SSN_avg    SSN_sd
# 1             Build LSN  5.938159 0.2277486 157.45195 13.985799
# 2 Stream updist and AFV 12.977719 1.5048399 192.78123 10.051779
# 3        Add obs to LSN  4.235065 0.1040615 239.19207  5.820844
# 4          Assemble SSN        NA        NA  42.80426  4.447670
# 5                 Total 23.150942        NA 632.22950        NA

# load("bench_results_thesis/results_with_extra/SSN_results_network4_rep1.rda")
runtimes_SSN = combine_runtimes_onemodel(bench_res_dir = "bench_results_thesis/results_with_extra/", 
                                         network = 4, nreps = 1, model = "SSN", 
                                         obs_only = FALSE)
#           Build LSN  168.5525270
# Stream updist + AFV  252.2298
#      Add obs to LSN  342.9061 (includes updist + AFV)
#        Assemble SSN   61.5808930
#         Obs distmat   61.5735891
#          Estimation 6095.7276890
#      Total of these 6982.571 (1.9 hours)

bench5 = combine_runtimes_predistonly(
  bench_res_dir,
  bench_nws$Network[5],
  bench_nws$nreps_S3N[5],
  bench_nws$nreps_SSN[5],
  obs_only_S3N = TRUE,
  obs_only_SSN = TRUE,
  predist_only = TRUE)
# predist_only_S3N = TRUE,
# predist_only_SSN = TRUE)
bench5$obs_only
#                    task  S3N_avg    S3N_sd   SSN_avg    SSN_sd
# 1             Build LSN 12.53213 0.3854917 2411.9265 127.96005
# 2 Stream updist and AFV 28.13694 1.2296416 1582.2443  62.75788
# 3        Add obs to LSN  6.34922 0.7430853  429.5282  12.77543
# 4          Assemble SSN       NA        NA  399.8852  57.06011
# 5                 Total 47.01829        NA 4823.5842        NA

# load("bench_results_thesis/results_with_extra/S3N_results_network5_rep1.rda")
# params
# #      beta_1      beta_2    sigma.sq      tau.sq      lambda 
# # -44.0828437   0.5001258   4.9754502   0.1473866   5.2363228 
runtimes_S3N = combine_runtimes_onemodel(bench_res_dir = "bench_results_thesis/results_with_extra/", 
                                         network = 5, nreps = 1, model = "S3N", 
                                         obs_only = TRUE); runtimes_S3N

# 1                 Build LSN   11.918482
# 2     Stream updist and AFV   27.487936
# 3 Prep to compute distances    6.754933
# 4               Obs distmat 3606.981447
# 5                Estimation    1.918258
#                       Total 3655.061 (~ 1 hr)

bench6 = combine_runtimes_predistonly(
  bench_res_dir,
  bench_nws$Network[6],
  bench_nws$nreps_S3N[6],
  bench_nws$nreps_SSN[6],
  obs_only_S3N = TRUE,
  obs_only_SSN = TRUE,
  predist_only = TRUE)
# predist_only_S3N = TRUE,
# predist_only_SSN = TRUE)
bench6$obs_only
#                    task   S3N_avg    S3N_sd SSN_avg SSN_sd
# 1             Build LSN  61.07921 2.4412279      NA     NA
# 2 Stream updist and AFV  70.73305 1.2198065      NA     NA
# 3        Add obs to LSN  18.02988 0.8822904      NA     NA
# 4          Assemble SSN        NA        NA      NA     NA
# 5                 Total 149.84214        NA      NA     NA

nw1_res = get_summary_for_table(1, bench1)
nw2_res = get_summary_for_table(2, bench2)
nw3_res = get_summary_for_table(3, bench3)
nw4_res = get_summary_for_table(4, bench4)
nw5_res = get_summary_for_table(5, bench5)
nw6_res = get_summary_for_table(6, bench6)
nwres = rbind(nw1_res, nw2_res, nw3_res, nw4_res, nw5_res, nw6_res)
nwres$Task = c("Build_LSN", "Stream_updist", "Obs_updist")
nwres_long = nwres
nwres = pivot_wider(nwres, names_from = "Task", values_from = S3N:SSN)
nwres$BuildLSNratio = nwres$SSN_Build_LSN/nwres$S3N_Build_LSN
nwres$Streamratio = nwres$SSN_Stream_updist/nwres$S3N_Stream_updist
nwres$Obsratio = nwres$SSN_Obs_updist/nwres$S3N_Obs_updist
nwres = nwres[,c(1, 2,5,8, 3,6,9, 4,7,10)]
#   Network S3N_Build_LSN SSN_Build_LSN BuildLSNratio S3N_Stream_updist SSN_Stream_updist Streamratio S3N_Obs_updist SSN_Obs_updist Obsratio
#     <dbl>         <dbl>         <dbl>         <dbl>             <dbl>             <dbl>       <dbl>          <dbl>          <dbl>    <dbl>
# 1       1         0.227          1.51          6.62              1.81              1.15       0.635           3.65           6.74     1.85
# 2       2         0.729          2.99          4.11              2.33              4.73       2.03            3.83          25.2      6.59
# 3       3         2.79          26.9           9.63              6.68             66.1        9.89            4.08         142.      34.9 
# 4       4         5.94         157.           26.5              13.0             193.        14.9             4.24         239.      56.5 
# 5       5        12.5         2412.          192.               28.1            1582.        56.2             6.35         430.      67.7 
# 6       6        61.1           NA            NA                70.7              NA         NA              18.0           NA       NA   


# make benchmarking plot

bench_res_plot = nwres_long %>%
  left_join(bench_nws, by = "Network") %>%
  mutate(Task = factor(
    ifelse(
      Task == "Build_LSN", 
      "Build stream network",
      ifelse(
        Task == "Stream_updist",
        "Compute stream updist",
        "Compute site updist"
      )), levels = c("Build stream network", 
                     "Compute stream updist", 
                     "Compute site updist")))

bench_long = bench_res_plot %>%
  pivot_longer(cols = S3N:SSN, names_to = "Software", values_to = "Time")

# ggplot(filter(bench_long, Task != "Compute site updist"), 
ggplot(bench_long, 
       aes(x = log(nreach, base = 10), 
           y = log(Time, base = 10), color = Software)) +
  facet_wrap(Task ~ .) +
  geom_line() + geom_point() +
  labs(x = TeX("log$_{10}$(Number of reaches)"),
       y = TeX("log$_{10}$(Time in seconds)")) +
  theme_bw() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.59, 0.21),
        text = element_text(size = 11)) #,
# legend.text = element_text(size = 12))

ggsave(filename = paste0(bench_res_dir, "benchmark_plot_nreaches.png"),
       width = 8, height = 3, units = "in")

ggplot(filter(bench_long, Task == "Compute site updist"), 
       aes(x = log(nobs, base = 10), 
           y = log(Time, base = 10), color = Software)) +
  facet_wrap(Task ~ .) +
  geom_line() + geom_point() +
  labs(x = TeX("log$_{10}$(Number of obs. points)"),
       y = TeX("log$_{10}$(Time in seconds)")) +
  theme_bw() +
  theme(legend.position = "none",
        # legend.position.inside = c(0.59, 0.21),
        text = element_text(size = 11))

ggsave(filename = paste0(bench_res_dir, "benchmark_plot_nobs.png"),
       width = 2.7, height = 3, units = "in")


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

bench_res_dir = "bench_results_thesis/"

valid1 = combine_params_bothmodels(bench_res_dir, 
                                   bench_nws$Network[1],
                                   bench_nws$nreps_S3N[1],
                                   bench_nws$nreps_SSN[1])
#   |Parameter |   bias_S3N|   bias_SSN|    sd_S3N|    sd_SSN|
#   |:---------|----------:|----------:|---------:|---------:|
#   |beta_1    |  0.0673622|  0.0565686| 1.3474402| 1.3499650|
#   |beta_2    | -0.0012077| -0.0011177| 0.0124137| 0.0123792|
#   |sigma.sq  | -0.2394252| -0.2042920| 0.7769548| 0.8108830|
#   |tau.sq    |  0.0775426|  0.1608641| 0.2443687| 0.3163467|
#   |lambda    |  1.5051939|  1.6518596| 7.7301693| 4.9518021|

valid2 = combine_params_bothmodels(bench_res_dir, 
                                   bench_nws$Network[2],
                                   bench_nws$nreps_S3N[2],
                                   bench_nws$nreps_SSN[2])
#   |Parameter |   bias_S3N|   bias_SSN|    sd_S3N|    sd_SSN|
#   |:---------|----------:|----------:|---------:|---------:|
#   |beta_1    | -0.2351882| -0.2253276| 0.4234030| 0.4251077|
#   |beta_2    |  0.0017473|  0.0016751| 0.0033948| 0.0033762|
#   |sigma.sq  |  0.0089140| -0.0710719| 0.3778211| 0.4634365|
#   |tau.sq    | -0.0253604|  0.0896720| 0.0968090| 0.2328690|
#   |lambda    | -0.1518889|  0.3572424| 0.8020116| 1.6412659|

valid3 = combine_params_bothmodels(bench_res_dir, 
                                   bench_nws$Network[3],
                                   bench_nws$nreps_S3N[3],
                                   bench_nws$nreps_SSN[3])
#   |Parameter |   bias_S3N|   bias_SSN|    sd_S3N|    sd_SSN|
#   |:---------|----------:|----------:|---------:|---------:|
#   |beta_1    |  0.0594826|  0.0651499| 0.2366170| 0.2349273|
#   |beta_2    | -0.0003354| -0.0003733| 0.0016183| 0.0016102|
#   |sigma.sq  | -0.0551622| -0.0727642| 0.1720123| 0.1924396|
#   |tau.sq    |  0.0062762|  0.0797815| 0.0500061| 0.1252996|
#   |lambda    |  0.0169189|  0.2995400| 0.4036738| 0.6629289|



## Plot benchmark networks --------------

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
bench_res_dir = "valid_results/"

# bench_nws = data.frame(
#   Network = 1:6,
#   nreps_S3N = c(50, 50, 50, 10, 10, 10),
#   nreps_SSN = c(50, 50, 50, 10,  2, NA) # for SSN network 5, do preprocessing only, no estimation
# )

bench_nws = data.frame(
  Network = 1:6,
  nreps_S3N = c(1, 50, 50, 10, 10, 10),
  nreps_SSN = c(1, 50, 50, 10,  2, NA) # for SSN network 5, do preprocessing only, no estimation
)

## network 1
streams = filter(streams_full, HUC10 == "0514020607")
benchmark_and_validate(streams, pred_path,
                       bench_nws$Network[1],
                       bench_nws$nreps_S3N[1],
                       bench_nws$nreps_SSN[1],
                       predist_only = TRUE,
                       bench_res_dir = bench_res_dir)

# ## network 4
# streams = filter(streams_full, HUC4 == "0514")
# benchmark_and_validate(streams, pred_path,
#                        bench_nws$Network[4], 
#                        bench_nws$nreps_S3N[4], 
#                        bench_nws$nreps_SSN[4],
#                        predist_only = TRUE,
#                        bench_res_dir = bench_res_dir)
# 
# ## network 6
# streams = filter(streams_full, HUC2 == "05")
# benchmark_and_validate(streams, pred_path,
#                        bench_nws$Network[6], 
#                        bench_nws$nreps_S3N[6], 
#                        bench_nws$nreps_SSN[6],
#                        predist_only = TRUE,
#                        onlyS3N = TRUE,
#                        bench_res_dir = bench_res_dir)
# 
# ## network 5
# streams = filter(streams_full, HUC4 %in% c("0514", "0513"))
# mem.maxVSize(60000)
# benchmark_and_validate(streams, pred_path,
#                        bench_nws$Network[5], 
#                        bench_nws$nreps_S3N[5], 
#                        bench_nws$nreps_SSN[5],
#                        predist_only = TRUE,
#                        bench_res_dir = bench_res_dir)
# mem.maxVSize(16384)
# Sys.sleep(30)
# beep(3)

# Generate more SSN validation results for Network 3 (20 min each rep)
# nreps = 20; network = 3; out_dir = bench_res_dir
# lsn.path = "lsn"
# ndigits = ceiling(log(nreps, base = 10))
# for(rep in 1:nreps){
#   message(paste("Network", network, "SSN benchmark rep", rep, "of", nreps))
#   SSN_preproc_and_estimation(network, 
#                              str_pad(rep, ndigits, side = "left", pad = "0"),
#                              # need estimation for validation results
#                              out_dir, lsn.path, preproc_only = FALSE, predist_only = FALSE)
# }; beep(3)
# 
# network = 3; out_dir = bench_res_dir
# lsn.path = "lsn"
# ndigits = ceiling(log(nreps, base = 10))
# for(rep in 21:35){
#   message(paste("Network", network, "SSN benchmark rep", rep, "of", 35))
#   SSN_preproc_and_estimation(network, 
#                              str_pad(rep, ndigits, side = "left", pad = "0"),
#                              # need estimation for validation results
#                              out_dir, lsn.path, preproc_only = FALSE, predist_only = FALSE)
# }
# beep(3)

# network = 3; out_dir = bench_res_dir; nreps = 50
# lsn.path = "lsn"
# ndigits = ceiling(log(nreps, base = 10))
# for(rep in 36:50){
#   message(paste("Network", network, "SSN benchmark rep", rep, "of", 50))
#   SSN_preproc_and_estimation(network, 
#                              str_pad(rep, ndigits, side = "left", pad = "0"),
#                              # need estimation for validation results
#                              out_dir, lsn.path, preproc_only = FALSE, predist_only = FALSE)
# }
# beep(3)
