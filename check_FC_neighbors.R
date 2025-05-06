# check out whether code to prioritize FC before FU points is working
# answer from this script: 
# - yes for prediction points (all points after the first n points), 
# - no for observation points (the first n points)

library(tidyverse)

pwdist_predsobs_dir = "Region5_cluster_results/"
load(paste0(pwdist_predsobs_dir, "pred_neighbors.rda"))
nnFC = (pred_neighbors$nnWght != 0)
nnFC = matrix(nnFC, ncol=10, byrow=TRUE)

rows_with_FC_after_FU = rep(FALSE, nrow(nnFC))
sum_FC_per_row = apply(nnFC, 1, sum)
table(sum_FC_per_row)
# 0     1     2     3     4     5     6     7     8     9    10 
# 4203  8224  9250  5001  9666  7989 10351 10170  8198  5416 90626 
sum_FC_per_row = apply(nnFC[1:8924,], 1, sum)
table(sum_FC_per_row) # obs only
# 1    2    3    4    5    6    7    8    9   10 
# 605  959 1276 1409 1388 1103  833  672  452  227 
sum_FC_per_row = apply(nnFC[8925:nrow(nnFC),], 1, sum)
table(sum_FC_per_row) # preds only
# 0     1     2     3     4     5     6     7     8     9    10 
# 4203  7619  8291  3725  8257  6601  9248  9337  7526  4964 90399 

for(i in 1:nrow(nnFC)){
  # if FALSE is followed by TRUE anywhere in the row (i.e. FC point follows FU point)
  FC_after_FU = str_detect(str_flatten(nnFC[i,], collapse = " "), pattern = "FALSE TRUE")
  if(FC_after_FU){
    rows_with_FC_after_FU[i] = TRUE
  }
}

sum(rows_with_FC_after_FU)
# [1] 7211
min(which(rows_with_FC_after_FU)) # 1
max(which(rows_with_FC_after_FU)) # 8920 < 8924, number of observation points

# so this only happens for observation points
# apparently the code is working properly for prediction points but not obs

unique(nnFC[8925:nrow(nnFC),])
#        [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10]
# [1,]   TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
# [2,]   TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE
# [3,]   TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE
# [4,]   TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE
# [5,]   TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE
# [6,]   TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE
# [7,]   TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE
# [8,]   TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# [9,]   TRUE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# [10,]  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# [11,] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
