################################################################################
# DESCRIPTION: PROCESS ACTIGRAPH DATA FILES
# AUTHOR: SHI XIN
# CONTACT: STEVEXINSHI@GMAIL.COM
# LAST EDIT: 2018-12-03
################################################################################

# GUIDELINES
#   use 15 sec files to mark wear (vm)
#   good days?
#   aggregate to 60 sec to mark sleep (axis1)
#   map the sleep markers back to 15 sec files
#   use 15 sec files to mark activity levels (vm)

# R-PACKAGES
library(PhysicalActivity)  # read AGD file; mark wear status
library(PhysActBedRest)  # mark sleep status
library(RSQLite)  # utilized by PhysicalActivity

# WORK DIRECTORY
wd <- "L:/Studies/Soslow/scripts"
setwd(wd)
datad <- paste0(wd, "/data/")

# CUTPOINTS
cp <- read.csv("cutpoints.csv", stringsAsFactors = F)

d_wd <- c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday")
d_we <- c("Saturday", "Sunday")

# files
file_paths <- dir(datad, full.names = T)
file_names <- dir(datad)

# gather file info for availability ###########################################
file_info <- data.frame(id=1:49, v1_w=NA, 
                        v1_a=NA, v2_w=NA, v2_a=NA, v3_w=NA, v3_a=NA)

for (u in 1:length(file_names)) {
  id <-  as.numeric(substr(file_names[u], 4, 5))
  if (substr(file_names[u], 7, 7) == "W") {
    if (substr(file_names[u], 10, 10) == "1") {
      col <- 2
    } else if (substr(file_names[u], 10, 10) == "2") {
      col <- 4
    } else {
      col <- 6
    }
    
  } else {
    if (substr(file_names[u], 10, 10) == "1") {
      col <- 3
    } else if (substr(file_names[u], 10, 10) == "2") {
      col <- 5
    } else {
      col <- 7
    }
  }
  
  file_info[file_info$id == id , col] <- 1
}
# write.csv(file_info, "file_availability_2018-12-03.csv", row.names = F)




# process wear and sleep ######################################################
summ_all <- data.frame(id=NA, visit=NA, start=NA, end=NA, days=NA, n_wk=NA, n_wknd=NA,
                       wear_sec=NA, sleep_sec=NA, awake_sec=NA, 
                       
                       wrist_vm_per_wear_min=NA, wrist_vm_per_awake_min=NA, wrist_vm_per_sleep_min=NA, 
                       ankle_vm_per_wear_min=NA, ankle_vm_per_awake_min=NA, ankle_vm_per_sleep_min=NA, 
                       
                       wrist_vm=NA, wrist_vm_awake=NA, wrist_vm_sleep=NA,
                       ankle_vm=NA, ankle_vm_awake=NA, ankle_vm_sleep=NA,
                       
                       sedentary_1_ratio=NA, sedentary_2_ratio=NA, sedentary_3_ratio=NA,
                       low_intensity_1_ratio=NA, low_intensity_2_ratio=NA,
                       moderate_intensity_ratio=NA, vigorous_intensity_ratio=NA,
                       
                       sedentary_1=NA, sedentary_2=NA, sedentary_3=NA, low_intensity_1=NA, low_intensity_2=NA, moderate_intensity=NA, vigorous_intensity=NA,
                       sedentary_1_sec=NA, sedentary_2_sec=NA, sedentary_3_sec=NA, low_intensity_1_sec=NA, low_intensity_2_sec=NA, moderate_intensity_sec=NA, vigorous_intensity_sec=NA,
                       
                       sedentary_1_per_min=NA, sedentary_2_per_min=NA, sedentary_3_per_min=NA, low_intensity_1_per_min=NA, low_intensity_2_per_min=NA, moderate_intensity_per_min=NA, vigorous_intensity_per_min=NA,
                       
                       wrist_vm_per_awake_over_per_wear=NA, wrist_vm_per_awake_over_per_sleep=NA, ankle_vm_per_awake_over_per_wear=NA, ankle_vm_per_awake_over_per_sleep=NA,
                       
                       
                       # week day only
                       wk_wear_sec=NA, wk_sleep_sec=NA, wk_awake_sec=NA, 
                       
                       wk_wrist_vm_per_wear_min=NA, wk_wrist_vm_per_awake_min=NA, wk_wrist_vm_per_sleep_min=NA, 
                       wk_ankle_vm_per_wear_min=NA, wk_ankle_vm_per_awake_min=NA, wk_ankle_vm_per_sleep_min=NA, 
                       
                       wk_wrist_vm=NA, wk_wrist_vm_awake=NA, wk_wrist_vm_sleep=NA,
                       wk_ankle_vm=NA, wk_ankle_vm_awake=NA, wk_ankle_vm_sleep=NA,
                       
                       wk_sedentary_1_ratio=NA, wk_sedentary_2_ratio=NA, wk_sedentary_3_ratio=NA,
                       wk_low_intensity_1_ratio=NA, wk_low_intensity_2_ratio=NA,
                       wk_moderate_intensity_ratio=NA, wk_vigorous_intensity_ratio=NA,
                       
                       wk_sedentary_1=NA, wk_sedentary_2=NA, wk_sedentary_3=NA, wk_low_intensity_1=NA, wk_low_intensity_2=NA, wk_moderate_intensity=NA, wk_vigorous_intensity=NA,
                       wk_sedentary_1_sec=NA, wk_sedentary_2_sec=NA, wk_sedentary_3_sec=NA, wk_low_intensity_1_sec=NA, wk_low_intensity_2_sec=NA, wk_moderate_intensity_sec=NA, wk_vigorous_intensity_sec=NA,
                       
                       
                       # weekend day only
                       wknd_wear_sec=NA, wknd_sleep_sec=NA, wknd_awake_sec=NA, 
                       
                       wknd_wrist_vm_per_wear_min=NA, wknd_wrist_vm_per_awake_min=NA, wknd_wrist_vm_per_sleep_min=NA, 
                       wknd_ankle_vm_per_wear_min=NA, wknd_ankle_vm_per_awake_min=NA, wknd_ankle_vm_per_sleep_min=NA, 
                       
                       wknd_wrist_vm=NA, wknd_wrist_vm_awake=NA, wknd_wrist_vm_sleep=NA,
                       wknd_ankle_vm=NA, wknd_ankle_vm_awake=NA, wknd_ankle_vm_sleep=NA,
                       
                       wknd_sedentary_1_ratio=NA, wknd_sedentary_2_ratio=NA, wknd_sedentary_3_ratio=NA,
                       wknd_low_intensity_1_ratio=NA, wknd_low_intensity_2_ratio=NA,
                       wknd_moderate_intensity_ratio=NA, wknd_vigorous_intensity_ratio=NA,
                       
                       wknd_sedentary_1=NA, wknd_sedentary_2=NA, wknd_sedentary_3=NA, wknd_low_intensity_1=NA, wknd_low_intensity_2=NA, wknd_moderate_intensity=NA, wknd_vigorous_intensity=NA,
                       wknd_sedentary_1_sec=NA, wknd_sedentary_2_sec=NA, wknd_sedentary_3_sec=NA, wknd_low_intensity_1_sec=NA, wknd_low_intensity_2_sec=NA, wknd_moderate_intensity_sec=NA, wknd_vigorous_intensity_sec=NA,
                       
                       
                       # sections 1 of day 0700 - 1200 
                       s1_wear_sec=NA, s1_sleep_sec=NA, s1_awake_sec=NA, 
                       
                       s1_wrist_vm_per_wear_min=NA, s1_wrist_vm_per_awake_min=NA, s1_wrist_vm_per_sleep_min=NA, 
                       s1_ankle_vm_per_wear_min=NA, s1_ankle_vm_per_awake_min=NA, s1_ankle_vm_per_sleep_min=NA, 
                       
                       s1_wrist_vm=NA, s1_wrist_vm_awake=NA, s1_wrist_vm_sleep=NA,
                       s1_ankle_vm=NA, s1_ankle_vm_awake=NA, s1_ankle_vm_sleep=NA,
                       
                       s1_sedentary_1_ratio=NA, s1_sedentary_2_ratio=NA, s1_sedentary_3_ratio=NA,
                       s1_low_intensity_1_ratio=NA, s1_low_intensity_2_ratio=NA,
                       s1_moderate_intensity_ratio=NA, s1_vigorous_intensity_ratio=NA,
                       
                       s1_sedentary_1=NA, s1_sedentary_2=NA, s1_sedentary_3=NA, s1_low_intensity_1=NA, s1_low_intensity_2=NA, s1_moderate_intensity=NA, s1_vigorous_intensity=NA,
                       s1_sedentary_1_sec=NA, s1_sedentary_2_sec=NA, s1_sedentary_3_sec=NA, s1_low_intensity_1_sec=NA, s1_low_intensity_2_sec=NA, s1_moderate_intensity_sec=NA, s1_vigorous_intensity_sec=NA,
                       
                       
                       
                       
                       # sections 2 of day 1200 - 1500
                       s2_wear_sec=NA, s2_sleep_sec=NA, s2_awake_sec=NA, 
                       
                       s2_wrist_vm_per_wear_min=NA, s2_wrist_vm_per_awake_min=NA, s2_wrist_vm_per_sleep_min=NA, 
                       s2_ankle_vm_per_wear_min=NA, s2_ankle_vm_per_awake_min=NA, s2_ankle_vm_per_sleep_min=NA, 
                       
                       s2_wrist_vm=NA, s2_wrist_vm_awake=NA, s2_wrist_vm_sleep=NA,
                       s2_ankle_vm=NA, s2_ankle_vm_awake=NA, s2_ankle_vm_sleep=NA,
                       
                       s2_sedentary_1_ratio=NA, s2_sedentary_2_ratio=NA, s2_sedentary_3_ratio=NA,
                       s2_low_intensity_1_ratio=NA, s2_low_intensity_2_ratio=NA,
                       s2_moderate_intensity_ratio=NA, s2_vigorous_intensity_ratio=NA,
                       
                       s2_sedentary_1=NA, s2_sedentary_2=NA, s2_sedentary_3=NA, s2_low_intensity_1=NA, s2_low_intensity_2=NA, s2_moderate_intensity=NA, s2_vigorous_intensity=NA,
                       s2_sedentary_1_sec=NA, s2_sedentary_2_sec=NA, s2_sedentary_3_sec=NA, s2_low_intensity_1_sec=NA, s2_low_intensity_2_sec=NA, s2_moderate_intensity_sec=NA, s2_vigorous_intensity_sec=NA,
                       
                       
                       
                       # sections 3 of day 1500 - 1900
                       s3_wear_sec=NA, s3_sleep_sec=NA, s3_awake_sec=NA, 
                       
                       s3_wrist_vm_per_wear_min=NA, s3_wrist_vm_per_awake_min=NA, s3_wrist_vm_per_sleep_min=NA, 
                       s3_ankle_vm_per_wear_min=NA, s3_ankle_vm_per_awake_min=NA, s3_ankle_vm_per_sleep_min=NA, 
                       
                       s3_wrist_vm=NA, s3_wrist_vm_awake=NA, s3_wrist_vm_sleep=NA,
                       s3_ankle_vm=NA, s3_ankle_vm_awake=NA, s3_ankle_vm_sleep=NA,
                       
                       s3_sedentary_1_ratio=NA, s3_sedentary_2_ratio=NA, s3_sedentary_3_ratio=NA,
                       s3_low_intensity_1_ratio=NA, s3_low_intensity_2_ratio=NA,
                       s3_moderate_intensity_ratio=NA, s3_vigorous_intensity_ratio=NA,
                       
                       s3_sedentary_1=NA, s3_sedentary_2=NA, s3_sedentary_3=NA, s3_low_intensity_1=NA, s3_low_intensity_2=NA, s3_moderate_intensity=NA, s3_vigorous_intensity=NA,
                       s3_sedentary_1_sec=NA, s3_sedentary_2_sec=NA, s3_sedentary_3_sec=NA, s3_low_intensity_1_sec=NA, s3_low_intensity_2_sec=NA, s3_moderate_intensity_sec=NA, s3_vigorous_intensity_sec=NA,
                       
                       
                       
                       # sections 4 of day 1900 - 2200
                       s4_wear_sec=NA, s4_sleep_sec=NA, s4_awake_sec=NA, 
                       
                       s4_wrist_vm_per_wear_min=NA, s4_wrist_vm_per_awake_min=NA, s4_wrist_vm_per_sleep_min=NA, 
                       s4_ankle_vm_per_wear_min=NA, s4_ankle_vm_per_awake_min=NA, s4_ankle_vm_per_sleep_min=NA, 
                       
                       s4_wrist_vm=NA, s4_wrist_vm_awake=NA, s4_wrist_vm_sleep=NA,
                       s4_ankle_vm=NA, s4_ankle_vm_awake=NA, s4_ankle_vm_sleep=NA,
                       
                       s4_sedentary_1_ratio=NA, s4_sedentary_2_ratio=NA, s4_sedentary_3_ratio=NA,
                       s4_low_intensity_1_ratio=NA, s4_low_intensity_2_ratio=NA,
                       s4_moderate_intensity_ratio=NA, s4_vigorous_intensity_ratio=NA,
                       
                       s4_sedentary_1=NA, s4_sedentary_2=NA, s4_sedentary_3=NA, s4_low_intensity_1=NA, s4_low_intensity_2=NA, s4_moderate_intensity=NA, s4_vigorous_intensity=NA,
                       s4_sedentary_1_sec=NA, s4_sedentary_2_sec=NA, s4_sedentary_3_sec=NA, s4_low_intensity_1_sec=NA, s4_low_intensity_2_sec=NA, s4_moderate_intensity_sec=NA, s4_vigorous_intensity_sec=NA
                       )

# loop
for (i in 1:length(file_paths)) {
  print(paste0("Working on ", file_names[i], " ------------------------------"))
  
  # ankle data?
  if (substr(file_names[i], 7, 7) == "A") {
    print("This is ankle data, jump over.")
    next
  } else {
    print("read file...")
    dframe <- readActigraph(file_paths[i])
  }
  
  # CHECK if no data
  if (nrow(dframe) == 0) {
    print("This file does not have data")
    next
  }
  
  dframe <- dframe[, c("TimeStamp", "axis1", "axis2", "axis3", "vm")]
  colnames(dframe)[1] <- "ts"
  
  # mark wearing and sleeping
  print("marking wearing status...")
  dframe <- wearingMarking(dframe, perMinuteCts = 1, cts = "vm", TS = "ts")
  dframe$wear <- ifelse(dframe$wearing == "w", 1, 0)
  
  # aggregate to 60sec to mark sleep, then map it back
  print("marking sleep status...")
  dframe_60 <- aggregate(dframe[, "axis1"], by=list(cut(dframe$ts, "mins")), FUN = sum)
  colnames(dframe_60) <- c("ts_60", "axis1")
  dframe_60$ts_60 <- as.character(dframe_60$ts_60)
  dframe_60 <- markbedrest(dframe_60, age="youth", loc="wrist", TS="ts_60", cts="axis1", EC=F)
  dframe_60$sleep <- ifelse(dframe_60$bedrest == "br", 1, 0)
  dframe$ts_60 <- as.character(cut(dframe$ts, "mins"))
  dframe <- merge(dframe, dframe_60[, c("ts_60", "sleep")],
                 by = "ts_60",
                 all.x = T)
  dframe$ts_60 <- NULL
  
  # valid day (>=600 mins wearing during 06:00 to 21:00)
  days.all <- unique(cut(dframe$ts, 'day'))
  for (h in 1:length(days.all)) {
    d.day <- dframe[cut(dframe$ts, 'day') == days.all[h], ]
    
    # take 6-21 hr data
    d.day$hr <- as.numeric(format(d.day$ts, "%H"))
    d.day <- d.day[d.day$hr >= 6 & d.day$hr < 21, ]
    
    if (nrow(d.day) == 0) {days.all[h] <- NA}
    
    if (nrow(d.day[d.day$wear == 1,])/4 < 600) {
      days.all[h] <- NA
    }
  }
  dframe <- dframe[cut(dframe$ts, 'day') %in% days.all,]
  if (nrow(dframe) == 0) {
    print(paste(file_names[i], "has no valid days", sep = " "))
    next
  }
  n_days <- sum(ifelse(is.na(days.all), 0, 1)) # how many days left
  start_d <- format(dframe[1, "ts"], "%Y-%m-%d")
  end_d <- format(dframe[nrow(dframe), "ts"], "%Y-%m-%d")
  # mark categories
  print("marking activity levels...")
  for (m in 1:nrow(dframe)) {
    for (mm in 1:nrow(cp)) {
      if (dframe[m, "vm"] >= cp[mm, "start"] / 4 &
          dframe[m, "vm"] <= cp[mm, "end"] / 4) {
        dframe[m, "activity"] <- cp[mm, "activity"]
      }
    }
  }
  
  
  # [ANKLE]
  # CHECK if there is an ankle file for this id and visit
  id <- as.numeric(substr(file_names[i], 4, 5))
  visit <- as.numeric(substr(file_names[i], 10, 10))
  
  if (!is.na(file_info[id, (visit*2 + 1)])) {
    print("Ankle file exists. Integrating...")
    id <- as.character(id)
    if (nchar(id) == 1) {
      id <- paste0("0", id)
    }
    ankle_file_name <- paste0("DMD",id,"_A_","V",visit,"_15.AGD")
    ankle_file_path <- paste0(datad, ankle_file_name)
    dframe_a <- readActigraph(ankle_file_path)
    dframe_a <- dframe_a[, c("TimeStamp", "vm")]
    colnames(dframe_a) <- c("ts", "ankle_vm")
    dframe <- merge(dframe, dframe_a,
                   by = "ts",
                   all.x = T)
  }
  
  # clean up
  dframe$wearing <- NULL
  dframe$days <- NULL
  
  # 15 sec summ file -----------------------------------------------------------
  summ_15sec_file_name <- paste0("./results/DMD",id,"_V",visit, "_", start_d, "_", end_d,"_by15sec.csv")
  # write.csv(dframe, summ_15sec_file_name, row.names = F)
  
  # subject level summ file -----------------------------------------------------------
  dframe <- dframe[dframe$wear == 1, ]
  
  # fill out the mother table
  summ_all[i, "id"] <- id
  summ_all[i, "visit"] <- visit
  
  summ_all[i, "start"] <- start_d
  summ_all[i, "end"] <- end_d
  summ_all[i, "days"] <- n_days
  
  dframe_wd <- dframe[dframe$weekday %in% d_wd,]
  if (nrow(dframe_wd) != 0) {
    day <- unique(cut(dframe_wd$ts, 'day'))
    summ_all[i, "n_wk"] <- length(day)
  } else {
    summ_all[i, "n_wk"] <- 0
  }
  
  dframe_we <- dframe[dframe$weekday %in% d_we,]
  if (nrow(dframe_we) != 0) {
    day <- unique(cut(dframe_we$ts, 'day'))
    summ_all[i, "n_wknd"] <- length(day)
  } else {
    summ_all[i, "n_wknd"] <- 0
  }
  
  
  summ_all[i, "wear_sec"] <- nrow(dframe)*15
  summ_all[i, "sleep_sec"] <- nrow(dframe[dframe$sleep == 1, ])*15
  summ_all[i, "awake_sec"] <- nrow(dframe[dframe$sleep == 0, ])*15
  
  summ_all[i, "wrist_vm"] <- sum(dframe[, "vm"])
  summ_all[i, "wrist_vm_awake"] <- sum(dframe[dframe$sleep == 0, "vm"])
  summ_all[i, "wrist_vm_sleep"] <- sum(dframe[dframe$sleep == 1, "vm"])
  
  summ_all[i, "wrist_vm_per_wear_min"] <- summ_all[i, "wrist_vm"] / (summ_all[i, "wear_sec"] / 60)
  summ_all[i, "wrist_vm_per_awake_min"] <- summ_all[i, "wrist_vm_awake"] / (summ_all[i, "awake_sec"] / 60)
  summ_all[i, "wrist_vm_per_sleep_min"] <- summ_all[i, "wrist_vm_sleep"] / (summ_all[i, "sleep_sec"] / 60)
  
  summ_all[i, "sedentary_1"] <- sum(dframe[dframe$sleep == 0 & dframe$activity == "sedentary_1", "vm"])
  summ_all[i, "sedentary_2"] <- sum(dframe[dframe$sleep == 0 & dframe$activity == "sedentary_2", "vm"])
  summ_all[i, "sedentary_3"] <- sum(dframe[dframe$sleep == 0 & dframe$activity == "sedentary_3", "vm"])
  summ_all[i, "low_intensity_1"] <- sum(dframe[dframe$sleep == 0 & dframe$activity == "low_intensity_1", "vm"])
  summ_all[i, "low_intensity_2"] <- sum(dframe[dframe$sleep == 0 & dframe$activity == "low_intensity_2", "vm"])
  summ_all[i, "moderate_intensity"] <- sum(dframe[dframe$sleep == 0 & dframe$activity == "moderate_intensity", "vm"])
  summ_all[i, "vigorous_intensity"] <- sum(dframe[dframe$sleep == 0 & dframe$activity == "vigorous_intensity", "vm"])
  
  summ_all[i, "sedentary_1_sec"] <- nrow(dframe[dframe$sleep == 0 & dframe$activity == "sedentary_1", ])*15
  summ_all[i, "sedentary_2_sec"] <- nrow(dframe[dframe$sleep == 0 & dframe$activity == "sedentary_2", ])*15
  summ_all[i, "sedentary_3_sec"] <- nrow(dframe[dframe$sleep == 0 & dframe$activity == "sedentary_3", ])*15
  summ_all[i, "low_intensity_1_sec"] <- nrow(dframe[dframe$sleep == 0 & dframe$activity == "low_intensity_1", ])*15
  summ_all[i, "low_intensity_2_sec"] <- nrow(dframe[dframe$sleep == 0 & dframe$activity == "low_intensity_2", ])*15
  summ_all[i, "moderate_intensity_sec"] <- nrow(dframe[dframe$sleep == 0 & dframe$activity == "moderate_intensity", ])*15
  summ_all[i, "vigorous_intensity_sec"] <- nrow(dframe[dframe$sleep == 0 & dframe$activity == "vigorous_intensity", ])*15
  
  summ_all[i, "sedentary_1_ratio"] <- summ_all[i, "sedentary_1_sec"] / summ_all[i, "awake_sec"]
  summ_all[i, "sedentary_2_ratio"] <- summ_all[i, "sedentary_2_sec"] / summ_all[i, "awake_sec"]
  summ_all[i, "sedentary_3_ratio"] <- summ_all[i, "sedentary_3_sec"] / summ_all[i, "awake_sec"]
  summ_all[i, "low_intensity_1_ratio"] <- summ_all[i, "low_intensity_1_sec"] / summ_all[i, "awake_sec"]
  summ_all[i, "low_intensity_2_ratio"] <- summ_all[i, "low_intensity_2_sec"] / summ_all[i, "awake_sec"]
  summ_all[i, "moderate_intensity_ratio"] <- summ_all[i, "moderate_intensity_sec"] / summ_all[i, "awake_sec"]
  summ_all[i, "vigorous_intensity_ratio"] <- summ_all[i, "vigorous_intensity_sec"] / summ_all[i, "awake_sec"]
  
  summ_all[i, "sedentary_1_per_min"] <- summ_all[i, "sedentary_1"] / (summ_all[i, "sedentary_1_sec"] / 60)
  summ_all[i, "sedentary_2_per_min"] <- summ_all[i, "sedentary_2"] / (summ_all[i, "sedentary_2_sec"] / 60)
  summ_all[i, "sedentary_3_per_min"] <- summ_all[i, "sedentary_3"] / (summ_all[i, "sedentary_3_sec"] / 60)
  summ_all[i, "low_intensity_1_per_min"] <- summ_all[i, "low_intensity_1"] / (summ_all[i, "low_intensity_1_sec"] / 60)
  summ_all[i, "low_intensity_2_per_min"] <- summ_all[i, "low_intensity_2"] / (summ_all[i, "low_intensity_2_sec"] / 60)
  summ_all[i, "moderate_intensity_per_min"] <- summ_all[i, "moderate_intensity"] / (summ_all[i, "moderate_intensity_sec"] / 60)
  summ_all[i, "vigorous_intensity_per_min"] <- summ_all[i, "vigorous_intensity"] / (summ_all[i, "vigorous_intensity_sec"] / 60)
  
  summ_all[i, "wrist_vm_per_awake_over_per_wear"] <- summ_all[i, "wrist_vm_per_awake_min"] / summ_all[i, "wrist_vm_per_wear_min"]
  summ_all[i, "wrist_vm_per_awake_over_per_sleep"] <- summ_all[i, "wrist_vm_per_awake_min"] / summ_all[i, "wrist_vm_per_sleep_min"]
  summ_all[i, "ankle_vm_per_awake_over_per_wear"] <- summ_all[i, "ankle_vm_per_awake_min"] / summ_all[i, "ankle_vm_per_wear_min"]
  summ_all[i, "ankle_vm_per_awake_over_per_sleep"] <- summ_all[i, "ankle_vm_per_awake_min"] / summ_all[i, "ankle_vm_per_sleep_min"]
  
  
  
  
  
  

  
  # week days and weekend days -----------------------------------------------
  # week days
  dframe_wk <- dframe[dframe$weekday %in% d_wd,]
  if (nrow(dframe_wk) != 0) {
    summ_all[i, "wk_wear_sec"] <- nrow(dframe_wk)*15
    summ_all[i, "wk_sleep_sec"] <- nrow(dframe_wk[dframe_wk$sleep == 1, ])*15
    summ_all[i, "wk_awake_sec"] <- nrow(dframe_wk[dframe_wk$sleep == 0, ])*15
    
    summ_all[i, "wk_wrist_vm"] <- sum(dframe_wk[, "vm"])
    summ_all[i, "wk_wrist_vm_awake"] <- sum(dframe_wk[dframe_wk$sleep == 0, "vm"])
    summ_all[i, "wk_wrist_vm_sleep"] <- sum(dframe_wk[dframe_wk$sleep == 1, "vm"])
    
    summ_all[i, "wk_wrist_vm_per_wear_min"] <- summ_all[i, "wk_wrist_vm"] / (summ_all[i, "wk_wear_sec"] / 60)
    summ_all[i, "wk_wrist_vm_per_awake_min"] <- summ_all[i, "wk_wrist_vm_awake"] / (summ_all[i, "wk_awake_sec"] / 60)
    summ_all[i, "wk_wrist_vm_per_sleep_min"] <- summ_all[i, "wk_wrist_vm_sleep"] / (summ_all[i, "wk_sleep_sec"] / 60)
    
    summ_all[i, "wk_sedentary_1"] <- sum(dframe_wk[dframe_wk$sleep == 0 & dframe_wk$activity == "sedentary_1", "vm"])
    summ_all[i, "wk_sedentary_2"] <- sum(dframe_wk[dframe_wk$sleep == 0 & dframe_wk$activity == "sedentary_2", "vm"])
    summ_all[i, "wk_sedentary_3"] <- sum(dframe_wk[dframe_wk$sleep == 0 & dframe_wk$activity == "sedentary_3", "vm"])
    summ_all[i, "wk_low_intensity_1"] <- sum(dframe_wk[dframe_wk$sleep == 0 & dframe_wk$activity == "low_intensity_1", "vm"])
    summ_all[i, "wk_low_intensity_2"] <- sum(dframe_wk[dframe_wk$sleep == 0 & dframe_wk$activity == "low_intensity_2", "vm"])
    summ_all[i, "wk_moderate_intensity"] <- sum(dframe_wk[dframe_wk$sleep == 0 & dframe_wk$activity == "moderate_intensity", "vm"])
    summ_all[i, "wk_vigorous_intensity"] <- sum(dframe_wk[dframe_wk$sleep == 0 & dframe_wk$activity == "vigorous_intensity", "vm"])
    
    summ_all[i, "wk_sedentary_1_sec"] <- nrow(dframe_wk[dframe_wk$sleep == 0 & dframe_wk$activity == "sedentary_1", ])*15
    summ_all[i, "wk_sedentary_2_sec"] <- nrow(dframe_wk[dframe_wk$sleep == 0 & dframe_wk$activity == "sedentary_2", ])*15
    summ_all[i, "wk_sedentary_3_sec"] <- nrow(dframe_wk[dframe_wk$sleep == 0 & dframe_wk$activity == "sedentary_3", ])*15
    summ_all[i, "wk_low_intensity_1_sec"] <- nrow(dframe_wk[dframe_wk$sleep == 0 & dframe_wk$activity == "low_intensity_1", ])*15
    summ_all[i, "wk_low_intensity_2_sec"] <- nrow(dframe_wk[dframe_wk$sleep == 0 & dframe_wk$activity == "low_intensity_2", ])*15
    summ_all[i, "wk_moderate_intensity_sec"] <- nrow(dframe_wk[dframe_wk$sleep == 0 & dframe_wk$activity == "moderate_intensity", ])*15
    summ_all[i, "wk_vigorous_intensity_sec"] <- nrow(dframe_wk[dframe_wk$sleep == 0 & dframe_wk$activity == "vigorous_intensity", ])*15
    
    summ_all[i, "wk_sedentary_1_ratio"] <- summ_all[i, "wk_sedentary_1_sec"] / summ_all[i, "wk_awake_sec"]
    summ_all[i, "wk_sedentary_2_ratio"] <- summ_all[i, "wk_sedentary_2_sec"] / summ_all[i, "wk_awake_sec"]
    summ_all[i, "wk_sedentary_3_ratio"] <- summ_all[i, "wk_sedentary_3_sec"] / summ_all[i, "wk_awake_sec"]
    summ_all[i, "wk_low_intensity_1_ratio"] <- summ_all[i, "wk_low_intensity_1_sec"] / summ_all[i, "wk_awake_sec"]
    summ_all[i, "wk_low_intensity_2_ratio"] <- summ_all[i, "wk_low_intensity_2_sec"] / summ_all[i, "wk_awake_sec"]
    summ_all[i, "wk_moderate_intensity_ratio"] <- summ_all[i, "wk_moderate_intensity_sec"] / summ_all[i, "wk_awake_sec"]
    summ_all[i, "wk_vigorous_intensity_ratio"] <- summ_all[i, "wk_vigorous_intensity_sec"] / summ_all[i, "wk_awake_sec"]
  }
  
  # weekend days
  dframe_wknd <- dframe[dframe$weekday %in% d_we,]
  if (nrow(dframe_wknd) != 0) {
    summ_all[i, "wknd_wear_sec"] <- nrow(dframe_wknd)*15
    summ_all[i, "wknd_sleep_sec"] <- nrow(dframe_wknd[dframe_wknd$sleep == 1, ])*15
    summ_all[i, "wknd_awake_sec"] <- nrow(dframe_wknd[dframe_wknd$sleep == 0, ])*15
    
    summ_all[i, "wknd_wrist_vm"] <- sum(dframe_wknd[, "vm"])
    summ_all[i, "wknd_wrist_vm_awake"] <- sum(dframe_wknd[dframe_wknd$sleep == 0, "vm"])
    summ_all[i, "wknd_wrist_vm_sleep"] <- sum(dframe_wknd[dframe_wknd$sleep == 1, "vm"])
    
    summ_all[i, "wknd_wrist_vm_per_wear_min"] <- summ_all[i, "wknd_wrist_vm"] / (summ_all[i, "wknd_wear_sec"] / 60)
    summ_all[i, "wknd_wrist_vm_per_awake_min"] <- summ_all[i, "wknd_wrist_vm_awake"] / (summ_all[i, "wknd_awake_sec"] / 60)
    summ_all[i, "wknd_wrist_vm_per_sleep_min"] <- summ_all[i, "wknd_wrist_vm_sleep"] / (summ_all[i, "wknd_sleep_sec"] / 60)
    
    summ_all[i, "wknd_sedentary_1"] <- sum(dframe_wknd[dframe_wknd$sleep == 0 & dframe_wknd$activity == "sedentary_1", "vm"])
    summ_all[i, "wknd_sedentary_2"] <- sum(dframe_wknd[dframe_wknd$sleep == 0 & dframe_wknd$activity == "sedentary_2", "vm"])
    summ_all[i, "wknd_sedentary_3"] <- sum(dframe_wknd[dframe_wknd$sleep == 0 & dframe_wknd$activity == "sedentary_3", "vm"])
    summ_all[i, "wknd_low_intensity_1"] <- sum(dframe_wknd[dframe_wknd$sleep == 0 & dframe_wknd$activity == "low_intensity_1", "vm"])
    summ_all[i, "wknd_low_intensity_2"] <- sum(dframe_wknd[dframe_wknd$sleep == 0 & dframe_wknd$activity == "low_intensity_2", "vm"])
    summ_all[i, "wknd_moderate_intensity"] <- sum(dframe_wknd[dframe_wknd$sleep == 0 & dframe_wknd$activity == "moderate_intensity", "vm"])
    summ_all[i, "wknd_vigorous_intensity"] <- sum(dframe_wknd[dframe_wknd$sleep == 0 & dframe_wknd$activity == "vigorous_intensity", "vm"])
    
    summ_all[i, "wknd_sedentary_1_sec"] <- nrow(dframe_wknd[dframe_wknd$sleep == 0 & dframe_wknd$activity == "sedentary_1", ])*15
    summ_all[i, "wknd_sedentary_2_sec"] <- nrow(dframe_wknd[dframe_wknd$sleep == 0 & dframe_wknd$activity == "sedentary_2", ])*15
    summ_all[i, "wknd_sedentary_3_sec"] <- nrow(dframe_wknd[dframe_wknd$sleep == 0 & dframe_wknd$activity == "sedentary_3", ])*15
    summ_all[i, "wknd_low_intensity_1_sec"] <- nrow(dframe_wknd[dframe_wknd$sleep == 0 & dframe_wknd$activity == "low_intensity_1", ])*15
    summ_all[i, "wknd_low_intensity_2_sec"] <- nrow(dframe_wknd[dframe_wknd$sleep == 0 & dframe_wknd$activity == "low_intensity_2", ])*15
    summ_all[i, "wknd_moderate_intensity_sec"] <- nrow(dframe_wknd[dframe_wknd$sleep == 0 & dframe_wknd$activity == "moderate_intensity", ])*15
    summ_all[i, "wknd_vigorous_intensity_sec"] <- nrow(dframe_wknd[dframe_wknd$sleep == 0 & dframe_wknd$activity == "vigorous_intensity", ])*15
    
    summ_all[i, "wknd_sedentary_1_ratio"] <- summ_all[i, "wknd_sedentary_1_sec"] / summ_all[i, "wknd_awake_sec"]
    summ_all[i, "wknd_sedentary_2_ratio"] <- summ_all[i, "wknd_sedentary_2_sec"] / summ_all[i, "wknd_awake_sec"]
    summ_all[i, "wknd_sedentary_3_ratio"] <- summ_all[i, "wknd_sedentary_3_sec"] / summ_all[i, "wknd_awake_sec"]
    summ_all[i, "wknd_low_intensity_1_ratio"] <- summ_all[i, "wknd_low_intensity_1_sec"] / summ_all[i, "wknd_awake_sec"]
    summ_all[i, "wknd_low_intensity_2_ratio"] <- summ_all[i, "wknd_low_intensity_2_sec"] / summ_all[i, "wknd_awake_sec"]
    summ_all[i, "wknd_moderate_intensity_ratio"] <- summ_all[i, "wknd_moderate_intensity_sec"] / summ_all[i, "wknd_awake_sec"]
    summ_all[i, "wknd_vigorous_intensity_ratio"] <- summ_all[i, "wknd_vigorous_intensity_sec"] / summ_all[i, "wknd_awake_sec"]
  }
  
  
  # sections -----------------------------------------------------------------
  dframe$hr <- as.numeric(format(dframe$ts, "%H"))
  
  
  # s1: [0700 - 1200)
  dframe_s1 <- dframe[dframe$hr %in% c(7:11),]
  if (nrow(dframe_s1) != 0) {
    summ_all[i, "s1_wear_sec"] <- nrow(dframe_s1)*15
    summ_all[i, "s1_sleep_sec"] <- nrow(dframe_s1[dframe_s1$sleep == 1, ])*15
    summ_all[i, "s1_awake_sec"] <- nrow(dframe_s1[dframe_s1$sleep == 0, ])*15
    
    summ_all[i, "s1_wrist_vm"] <- sum(dframe_s1[, "vm"])
    summ_all[i, "s1_wrist_vm_awake"] <- sum(dframe_s1[dframe_s1$sleep == 0, "vm"])
    summ_all[i, "s1_wrist_vm_sleep"] <- sum(dframe_s1[dframe_s1$sleep == 1, "vm"])
    
    summ_all[i, "s1_wrist_vm_per_wear_min"] <- summ_all[i, "s1_wrist_vm"] / (summ_all[i, "s1_wear_sec"] / 60)
    summ_all[i, "s1_wrist_vm_per_awake_min"] <- summ_all[i, "s1_wrist_vm_awake"] / (summ_all[i, "s1_awake_sec"] / 60)
    summ_all[i, "s1_wrist_vm_per_sleep_min"] <- summ_all[i, "s1_wrist_vm_sleep"] / (summ_all[i, "s1_sleep_sec"] / 60)
    
    summ_all[i, "s1_sedentary_1"] <- sum(dframe_s1[dframe_s1$sleep == 0 & dframe_s1$activity == "sedentary_1", "vm"])
    summ_all[i, "s1_sedentary_2"] <- sum(dframe_s1[dframe_s1$sleep == 0 & dframe_s1$activity == "sedentary_2", "vm"])
    summ_all[i, "s1_sedentary_3"] <- sum(dframe_s1[dframe_s1$sleep == 0 & dframe_s1$activity == "sedentary_3", "vm"])
    summ_all[i, "s1_low_intensity_1"] <- sum(dframe_s1[dframe_s1$sleep == 0 & dframe_s1$activity == "low_intensity_1", "vm"])
    summ_all[i, "s1_low_intensity_2"] <- sum(dframe_s1[dframe_s1$sleep == 0 & dframe_s1$activity == "low_intensity_2", "vm"])
    summ_all[i, "s1_moderate_intensity"] <- sum(dframe_s1[dframe_s1$sleep == 0 & dframe_s1$activity == "moderate_intensity", "vm"])
    summ_all[i, "s1_vigorous_intensity"] <- sum(dframe_s1[dframe_s1$sleep == 0 & dframe_s1$activity == "vigorous_intensity", "vm"])
    
    summ_all[i, "s1_sedentary_1_sec"] <- nrow(dframe_s1[dframe_s1$sleep == 0 & dframe_s1$activity == "sedentary_1", ])*15
    summ_all[i, "s1_sedentary_2_sec"] <- nrow(dframe_s1[dframe_s1$sleep == 0 & dframe_s1$activity == "sedentary_2", ])*15
    summ_all[i, "s1_sedentary_3_sec"] <- nrow(dframe_s1[dframe_s1$sleep == 0 & dframe_s1$activity == "sedentary_3", ])*15
    summ_all[i, "s1_low_intensity_1_sec"] <- nrow(dframe_s1[dframe_s1$sleep == 0 & dframe_s1$activity == "low_intensity_1", ])*15
    summ_all[i, "s1_low_intensity_2_sec"] <- nrow(dframe_s1[dframe_s1$sleep == 0 & dframe_s1$activity == "low_intensity_2", ])*15
    summ_all[i, "s1_moderate_intensity_sec"] <- nrow(dframe_s1[dframe_s1$sleep == 0 & dframe_s1$activity == "moderate_intensity", ])*15
    summ_all[i, "s1_vigorous_intensity_sec"] <- nrow(dframe_s1[dframe_s1$sleep == 0 & dframe_s1$activity == "vigorous_intensity", ])*15
    
    summ_all[i, "s1_sedentary_1_ratio"] <- summ_all[i, "s1_sedentary_1_sec"] / summ_all[i, "s1_awake_sec"]
    summ_all[i, "s1_sedentary_2_ratio"] <- summ_all[i, "s1_sedentary_2_sec"] / summ_all[i, "s1_awake_sec"]
    summ_all[i, "s1_sedentary_3_ratio"] <- summ_all[i, "s1_sedentary_3_sec"] / summ_all[i, "s1_awake_sec"]
    summ_all[i, "s1_low_intensity_1_ratio"] <- summ_all[i, "s1_low_intensity_1_sec"] / summ_all[i, "s1_awake_sec"]
    summ_all[i, "s1_low_intensity_2_ratio"] <- summ_all[i, "s1_low_intensity_2_sec"] / summ_all[i, "s1_awake_sec"]
    summ_all[i, "s1_moderate_intensity_ratio"] <- summ_all[i, "s1_moderate_intensity_sec"] / summ_all[i, "s1_awake_sec"]
    summ_all[i, "s1_vigorous_intensity_ratio"] <- summ_all[i, "s1_vigorous_intensity_sec"] / summ_all[i, "s1_awake_sec"]
  }
  
  
  
  
  
  
  
  
  
  # s2: [1200 - 1500)
  dframe_s2 <- dframe[dframe$hr %in% c(12:14),]
  if (nrow(dframe_s2) != 0) {
    summ_all[i, "s2_wear_sec"] <- nrow(dframe_s2)*15
    summ_all[i, "s2_sleep_sec"] <- nrow(dframe_s2[dframe_s2$sleep == 1, ])*15
    summ_all[i, "s2_awake_sec"] <- nrow(dframe_s2[dframe_s2$sleep == 0, ])*15
    
    summ_all[i, "s2_wrist_vm"] <- sum(dframe_s2[, "vm"])
    summ_all[i, "s2_wrist_vm_awake"] <- sum(dframe_s2[dframe_s2$sleep == 0, "vm"])
    summ_all[i, "s2_wrist_vm_sleep"] <- sum(dframe_s2[dframe_s2$sleep == 1, "vm"])
    
    summ_all[i, "s2_wrist_vm_per_wear_min"] <- summ_all[i, "s2_wrist_vm"] / (summ_all[i, "s2_wear_sec"] / 60)
    summ_all[i, "s2_wrist_vm_per_awake_min"] <- summ_all[i, "s2_wrist_vm_awake"] / (summ_all[i, "s2_awake_sec"] / 60)
    summ_all[i, "s2_wrist_vm_per_sleep_min"] <- summ_all[i, "s2_wrist_vm_sleep"] / (summ_all[i, "s2_sleep_sec"] / 60)
    
    summ_all[i, "s2_sedentary_1"] <- sum(dframe_s2[dframe_s2$sleep == 0 & dframe_s2$activity == "sedentary_1", "vm"])
    summ_all[i, "s2_sedentary_2"] <- sum(dframe_s2[dframe_s2$sleep == 0 & dframe_s2$activity == "sedentary_2", "vm"])
    summ_all[i, "s2_sedentary_3"] <- sum(dframe_s2[dframe_s2$sleep == 0 & dframe_s2$activity == "sedentary_3", "vm"])
    summ_all[i, "s2_low_intensity_1"] <- sum(dframe_s2[dframe_s2$sleep == 0 & dframe_s2$activity == "low_intensity_1", "vm"])
    summ_all[i, "s2_low_intensity_2"] <- sum(dframe_s2[dframe_s2$sleep == 0 & dframe_s2$activity == "low_intensity_2", "vm"])
    summ_all[i, "s2_moderate_intensity"] <- sum(dframe_s2[dframe_s2$sleep == 0 & dframe_s2$activity == "moderate_intensity", "vm"])
    summ_all[i, "s2_vigorous_intensity"] <- sum(dframe_s2[dframe_s2$sleep == 0 & dframe_s2$activity == "vigorous_intensity", "vm"])
    
    summ_all[i, "s2_sedentary_1_sec"] <- nrow(dframe_s2[dframe_s2$sleep == 0 & dframe_s2$activity == "sedentary_1", ])*15
    summ_all[i, "s2_sedentary_2_sec"] <- nrow(dframe_s2[dframe_s2$sleep == 0 & dframe_s2$activity == "sedentary_2", ])*15
    summ_all[i, "s2_sedentary_3_sec"] <- nrow(dframe_s2[dframe_s2$sleep == 0 & dframe_s2$activity == "sedentary_3", ])*15
    summ_all[i, "s2_low_intensity_1_sec"] <- nrow(dframe_s2[dframe_s2$sleep == 0 & dframe_s2$activity == "low_intensity_1", ])*15
    summ_all[i, "s2_low_intensity_2_sec"] <- nrow(dframe_s2[dframe_s2$sleep == 0 & dframe_s2$activity == "low_intensity_2", ])*15
    summ_all[i, "s2_moderate_intensity_sec"] <- nrow(dframe_s2[dframe_s2$sleep == 0 & dframe_s2$activity == "moderate_intensity", ])*15
    summ_all[i, "s2_vigorous_intensity_sec"] <- nrow(dframe_s2[dframe_s2$sleep == 0 & dframe_s2$activity == "vigorous_intensity", ])*15
    
    summ_all[i, "s2_sedentary_1_ratio"] <- summ_all[i, "s2_sedentary_1_sec"] / summ_all[i, "s2_awake_sec"]
    summ_all[i, "s2_sedentary_2_ratio"] <- summ_all[i, "s2_sedentary_2_sec"] / summ_all[i, "s2_awake_sec"]
    summ_all[i, "s2_sedentary_3_ratio"] <- summ_all[i, "s2_sedentary_3_sec"] / summ_all[i, "s2_awake_sec"]
    summ_all[i, "s2_low_intensity_1_ratio"] <- summ_all[i, "s2_low_intensity_1_sec"] / summ_all[i, "s2_awake_sec"]
    summ_all[i, "s2_low_intensity_2_ratio"] <- summ_all[i, "s2_low_intensity_2_sec"] / summ_all[i, "s2_awake_sec"]
    summ_all[i, "s2_moderate_intensity_ratio"] <- summ_all[i, "s2_moderate_intensity_sec"] / summ_all[i, "s2_awake_sec"]
    summ_all[i, "s2_vigorous_intensity_ratio"] <- summ_all[i, "s2_vigorous_intensity_sec"] / summ_all[i, "s2_awake_sec"]
  }
  
  
  
  
  
  
  
  
  
  # s3: [1500 - 1900)
  dframe_s3 <- dframe[dframe$hr %in% c(15:18),]
  if (nrow(dframe_s3) != 0) {
    summ_all[i, "s3_wear_sec"] <- nrow(dframe_s3)*15
    summ_all[i, "s3_sleep_sec"] <- nrow(dframe_s3[dframe_s3$sleep == 1, ])*15
    summ_all[i, "s3_awake_sec"] <- nrow(dframe_s3[dframe_s3$sleep == 0, ])*15
    
    summ_all[i, "s3_wrist_vm"] <- sum(dframe_s3[, "vm"])
    summ_all[i, "s3_wrist_vm_awake"] <- sum(dframe_s3[dframe_s3$sleep == 0, "vm"])
    summ_all[i, "s3_wrist_vm_sleep"] <- sum(dframe_s3[dframe_s3$sleep == 1, "vm"])
    
    summ_all[i, "s3_wrist_vm_per_wear_min"] <- summ_all[i, "s3_wrist_vm"] / (summ_all[i, "s3_wear_sec"] / 60)
    summ_all[i, "s3_wrist_vm_per_awake_min"] <- summ_all[i, "s3_wrist_vm_awake"] / (summ_all[i, "s3_awake_sec"] / 60)
    summ_all[i, "s3_wrist_vm_per_sleep_min"] <- summ_all[i, "s3_wrist_vm_sleep"] / (summ_all[i, "s3_sleep_sec"] / 60)
    
    summ_all[i, "s3_sedentary_1"] <- sum(dframe_s3[dframe_s3$sleep == 0 & dframe_s3$activity == "sedentary_1", "vm"])
    summ_all[i, "s3_sedentary_2"] <- sum(dframe_s3[dframe_s3$sleep == 0 & dframe_s3$activity == "sedentary_2", "vm"])
    summ_all[i, "s3_sedentary_3"] <- sum(dframe_s3[dframe_s3$sleep == 0 & dframe_s3$activity == "sedentary_3", "vm"])
    summ_all[i, "s3_low_intensity_1"] <- sum(dframe_s3[dframe_s3$sleep == 0 & dframe_s3$activity == "low_intensity_1", "vm"])
    summ_all[i, "s3_low_intensity_2"] <- sum(dframe_s3[dframe_s3$sleep == 0 & dframe_s3$activity == "low_intensity_2", "vm"])
    summ_all[i, "s3_moderate_intensity"] <- sum(dframe_s3[dframe_s3$sleep == 0 & dframe_s3$activity == "moderate_intensity", "vm"])
    summ_all[i, "s3_vigorous_intensity"] <- sum(dframe_s3[dframe_s3$sleep == 0 & dframe_s3$activity == "vigorous_intensity", "vm"])
    
    summ_all[i, "s3_sedentary_1_sec"] <- nrow(dframe_s3[dframe_s3$sleep == 0 & dframe_s3$activity == "sedentary_1", ])*15
    summ_all[i, "s3_sedentary_2_sec"] <- nrow(dframe_s3[dframe_s3$sleep == 0 & dframe_s3$activity == "sedentary_2", ])*15
    summ_all[i, "s3_sedentary_3_sec"] <- nrow(dframe_s3[dframe_s3$sleep == 0 & dframe_s3$activity == "sedentary_3", ])*15
    summ_all[i, "s3_low_intensity_1_sec"] <- nrow(dframe_s3[dframe_s3$sleep == 0 & dframe_s3$activity == "low_intensity_1", ])*15
    summ_all[i, "s3_low_intensity_2_sec"] <- nrow(dframe_s3[dframe_s3$sleep == 0 & dframe_s3$activity == "low_intensity_2", ])*15
    summ_all[i, "s3_moderate_intensity_sec"] <- nrow(dframe_s3[dframe_s3$sleep == 0 & dframe_s3$activity == "moderate_intensity", ])*15
    summ_all[i, "s3_vigorous_intensity_sec"] <- nrow(dframe_s3[dframe_s3$sleep == 0 & dframe_s3$activity == "vigorous_intensity", ])*15
    
    summ_all[i, "s3_sedentary_1_ratio"] <- summ_all[i, "s3_sedentary_1_sec"] / summ_all[i, "s3_awake_sec"]
    summ_all[i, "s3_sedentary_2_ratio"] <- summ_all[i, "s3_sedentary_2_sec"] / summ_all[i, "s3_awake_sec"]
    summ_all[i, "s3_sedentary_3_ratio"] <- summ_all[i, "s3_sedentary_3_sec"] / summ_all[i, "s3_awake_sec"]
    summ_all[i, "s3_low_intensity_1_ratio"] <- summ_all[i, "s3_low_intensity_1_sec"] / summ_all[i, "s3_awake_sec"]
    summ_all[i, "s3_low_intensity_2_ratio"] <- summ_all[i, "s3_low_intensity_2_sec"] / summ_all[i, "s3_awake_sec"]
    summ_all[i, "s3_moderate_intensity_ratio"] <- summ_all[i, "s3_moderate_intensity_sec"] / summ_all[i, "s3_awake_sec"]
    summ_all[i, "s3_vigorous_intensity_ratio"] <- summ_all[i, "s3_vigorous_intensity_sec"] / summ_all[i, "s3_awake_sec"]
  }
  
  
  
  
  
  
  
  
  # s4: [1900 - 2200)
  dframe_s4 <- dframe[dframe$hr %in% c(19:21),]
  if (nrow(dframe_s4) != 0) {
    summ_all[i, "s4_wear_sec"] <- nrow(dframe_s4)*15
    summ_all[i, "s4_sleep_sec"] <- nrow(dframe_s4[dframe_s4$sleep == 1, ])*15
    summ_all[i, "s4_awake_sec"] <- nrow(dframe_s4[dframe_s4$sleep == 0, ])*15
    
    summ_all[i, "s4_wrist_vm"] <- sum(dframe_s4[, "vm"])
    summ_all[i, "s4_wrist_vm_awake"] <- sum(dframe_s4[dframe_s4$sleep == 0, "vm"])
    summ_all[i, "s4_wrist_vm_sleep"] <- sum(dframe_s4[dframe_s4$sleep == 1, "vm"])
    
    summ_all[i, "s4_wrist_vm_per_wear_min"] <- summ_all[i, "s4_wrist_vm"] / (summ_all[i, "s4_wear_sec"] / 60)
    summ_all[i, "s4_wrist_vm_per_awake_min"] <- summ_all[i, "s4_wrist_vm_awake"] / (summ_all[i, "s4_awake_sec"] / 60)
    summ_all[i, "s4_wrist_vm_per_sleep_min"] <- summ_all[i, "s4_wrist_vm_sleep"] / (summ_all[i, "s4_sleep_sec"] / 60)
    
    summ_all[i, "s4_sedentary_1"] <- sum(dframe_s4[dframe_s4$sleep == 0 & dframe_s4$activity == "sedentary_1", "vm"])
    summ_all[i, "s4_sedentary_2"] <- sum(dframe_s4[dframe_s4$sleep == 0 & dframe_s4$activity == "sedentary_2", "vm"])
    summ_all[i, "s4_sedentary_3"] <- sum(dframe_s4[dframe_s4$sleep == 0 & dframe_s4$activity == "sedentary_3", "vm"])
    summ_all[i, "s4_low_intensity_1"] <- sum(dframe_s4[dframe_s4$sleep == 0 & dframe_s4$activity == "low_intensity_1", "vm"])
    summ_all[i, "s4_low_intensity_2"] <- sum(dframe_s4[dframe_s4$sleep == 0 & dframe_s4$activity == "low_intensity_2", "vm"])
    summ_all[i, "s4_moderate_intensity"] <- sum(dframe_s4[dframe_s4$sleep == 0 & dframe_s4$activity == "moderate_intensity", "vm"])
    summ_all[i, "s4_vigorous_intensity"] <- sum(dframe_s4[dframe_s4$sleep == 0 & dframe_s4$activity == "vigorous_intensity", "vm"])
    
    summ_all[i, "s4_sedentary_1_sec"] <- nrow(dframe_s4[dframe_s4$sleep == 0 & dframe_s4$activity == "sedentary_1", ])*15
    summ_all[i, "s4_sedentary_2_sec"] <- nrow(dframe_s4[dframe_s4$sleep == 0 & dframe_s4$activity == "sedentary_2", ])*15
    summ_all[i, "s4_sedentary_3_sec"] <- nrow(dframe_s4[dframe_s4$sleep == 0 & dframe_s4$activity == "sedentary_3", ])*15
    summ_all[i, "s4_low_intensity_1_sec"] <- nrow(dframe_s4[dframe_s4$sleep == 0 & dframe_s4$activity == "low_intensity_1", ])*15
    summ_all[i, "s4_low_intensity_2_sec"] <- nrow(dframe_s4[dframe_s4$sleep == 0 & dframe_s4$activity == "low_intensity_2", ])*15
    summ_all[i, "s4_moderate_intensity_sec"] <- nrow(dframe_s4[dframe_s4$sleep == 0 & dframe_s4$activity == "moderate_intensity", ])*15
    summ_all[i, "s4_vigorous_intensity_sec"] <- nrow(dframe_s4[dframe_s4$sleep == 0 & dframe_s4$activity == "vigorous_intensity", ])*15
    
    summ_all[i, "s4_sedentary_1_ratio"] <- summ_all[i, "s4_sedentary_1_sec"] / summ_all[i, "s4_awake_sec"]
    summ_all[i, "s4_sedentary_2_ratio"] <- summ_all[i, "s4_sedentary_2_sec"] / summ_all[i, "s4_awake_sec"]
    summ_all[i, "s4_sedentary_3_ratio"] <- summ_all[i, "s4_sedentary_3_sec"] / summ_all[i, "s4_awake_sec"]
    summ_all[i, "s4_low_intensity_1_ratio"] <- summ_all[i, "s4_low_intensity_1_sec"] / summ_all[i, "s4_awake_sec"]
    summ_all[i, "s4_low_intensity_2_ratio"] <- summ_all[i, "s4_low_intensity_2_sec"] / summ_all[i, "s4_awake_sec"]
    summ_all[i, "s4_moderate_intensity_ratio"] <- summ_all[i, "s4_moderate_intensity_sec"] / summ_all[i, "s4_awake_sec"]
    summ_all[i, "s4_vigorous_intensity_ratio"] <- summ_all[i, "s4_vigorous_intensity_sec"] / summ_all[i, "s4_awake_sec"]
  }
  
  
  
  # ankle data ----------------------------------------------------------------
  id <- as.numeric(substr(file_names[i], 4, 5))
  if (!is.na(file_info[id, (visit*2 + 1)])) {
    summ_all[i, "ankle_vm"] <- sum(dframe[, "ankle_vm"])
    summ_all[i, "ankle_vm_awake"] <- sum(dframe[dframe$sleep == 0, "ankle_vm"])
    summ_all[i, "ankle_vm_sleep"] <- sum(dframe[dframe$sleep == 1, "ankle_vm"])
    
    summ_all[i, "ankle_vm_per_wear_min"] <- summ_all[i, "ankle_vm"] / (summ_all[i, "wear_sec"] / 60)
    summ_all[i, "ankle_vm_per_awake_min"] <- summ_all[i, "ankle_vm_awake"] / (summ_all[i, "awake_sec"] / 60)
    summ_all[i, "ankle_vm_per_sleep_min"] <- summ_all[i, "ankle_vm_sleep"] / (summ_all[i, "sleep_sec"] / 60)
    
    
    # weekday
    if (nrow(dframe_wk) != 0) {
      summ_all[i, "wk_ankle_vm"] <- sum(dframe_wk[, "ankle_vm"])
      summ_all[i, "wk_ankle_vm_awake"] <- sum(dframe_wk[dframe_wk$sleep == 0, "ankle_vm"])
      summ_all[i, "wk_ankle_vm_sleep"] <- sum(dframe_wk[dframe_wk$sleep == 1, "ankle_vm"])
      
      summ_all[i, "wk_ankle_vm_per_wear_min"] <- summ_all[i, "wk_ankle_vm"] / (summ_all[i, "wk_wear_sec"] / 60)
      summ_all[i, "wk_ankle_vm_per_awake_min"] <- summ_all[i, "wk_ankle_vm_awake"] / (summ_all[i, "wk_awake_sec"] / 60)
      summ_all[i, "wk_ankle_vm_per_sleep_min"] <- summ_all[i, "wk_ankle_vm_sleep"] / (summ_all[i, "wk_sleep_sec"] / 60)
    }
    
    
    # weekend
    if (nrow(dframe_wknd) != 0) {
      summ_all[i, "wknd_ankle_vm"] <- sum(dframe_wknd[, "ankle_vm"])
      summ_all[i, "wknd_ankle_vm_awake"] <- sum(dframe_wknd[dframe_wknd$sleep == 0, "ankle_vm"])
      summ_all[i, "wknd_ankle_vm_sleep"] <- sum(dframe_wknd[dframe_wknd$sleep == 1, "ankle_vm"])
      
      summ_all[i, "wknd_ankle_vm_per_wear_min"] <- summ_all[i, "wknd_ankle_vm"] / (summ_all[i, "wknd_wear_sec"] / 60)
      summ_all[i, "wknd_ankle_vm_per_awake_min"] <- summ_all[i, "wknd_ankle_vm_awake"] / (summ_all[i, "wknd_awake_sec"] / 60)
      summ_all[i, "wknd_ankle_vm_per_sleep_min"] <- summ_all[i, "wknd_ankle_vm_sleep"] / (summ_all[i, "wknd_sleep_sec"] / 60)
    }
    
    
    # s1
    if (nrow(dframe_s1) != 0) {
      summ_all[i, "s1_ankle_vm"] <- sum(dframe_s1[, "ankle_vm"])
      summ_all[i, "s1_ankle_vm_awake"] <- sum(dframe_s1[dframe_s1$sleep == 0, "ankle_vm"])
      summ_all[i, "s1_ankle_vm_sleep"] <- sum(dframe_s1[dframe_s1$sleep == 1, "ankle_vm"])
      
      summ_all[i, "s1_ankle_vm_per_wear_min"] <- summ_all[i, "s1_ankle_vm"] / (summ_all[i, "s1_wear_sec"] / 60)
      summ_all[i, "s1_ankle_vm_per_awake_min"] <- summ_all[i, "s1_ankle_vm_awake"] / (summ_all[i, "s1_awake_sec"] / 60)
      summ_all[i, "s1_ankle_vm_per_sleep_min"] <- summ_all[i, "s1_ankle_vm_sleep"] / (summ_all[i, "s1_sleep_sec"] / 60)
    }
    
    
    
    # s2
    if (nrow(dframe_s2) != 0) {
      summ_all[i, "s2_ankle_vm"] <- sum(dframe_s2[, "ankle_vm"])
      summ_all[i, "s2_ankle_vm_awake"] <- sum(dframe_s2[dframe_s2$sleep == 0, "ankle_vm"])
      summ_all[i, "s2_ankle_vm_sleep"] <- sum(dframe_s2[dframe_s2$sleep == 1, "ankle_vm"])
      
      summ_all[i, "s2_ankle_vm_per_wear_min"] <- summ_all[i, "s2_ankle_vm"] / (summ_all[i, "s2_wear_sec"] / 60)
      summ_all[i, "s2_ankle_vm_per_awake_min"] <- summ_all[i, "s2_ankle_vm_awake"] / (summ_all[i, "s2_awake_sec"] / 60)
      summ_all[i, "s2_ankle_vm_per_sleep_min"] <- summ_all[i, "s2_ankle_vm_sleep"] / (summ_all[i, "s2_sleep_sec"] / 60)
    }
    
    
    
    
    # s3
    if (nrow(dframe_s3) != 0) {
      summ_all[i, "s3_ankle_vm"] <- sum(dframe_s3[, "ankle_vm"])
      summ_all[i, "s3_ankle_vm_awake"] <- sum(dframe_s3[dframe_s3$sleep == 0, "ankle_vm"])
      summ_all[i, "s3_ankle_vm_sleep"] <- sum(dframe_s3[dframe_s3$sleep == 1, "ankle_vm"])
      
      summ_all[i, "s3_ankle_vm_per_wear_min"] <- summ_all[i, "s3_ankle_vm"] / (summ_all[i, "s3_wear_sec"] / 60)
      summ_all[i, "s3_ankle_vm_per_awake_min"] <- summ_all[i, "s3_ankle_vm_awake"] / (summ_all[i, "s3_awake_sec"] / 60)
      summ_all[i, "s3_ankle_vm_per_sleep_min"] <- summ_all[i, "s3_ankle_vm_sleep"] / (summ_all[i, "s3_sleep_sec"] / 60)
    }
    
    
    
    
    # s4
    if (nrow(dframe_s4) != 0) {
      summ_all[i, "s4_ankle_vm"] <- sum(dframe_s4[, "ankle_vm"])
      summ_all[i, "s4_ankle_vm_awake"] <- sum(dframe_s4[dframe_s4$sleep == 0, "ankle_vm"])
      summ_all[i, "s4_ankle_vm_sleep"] <- sum(dframe_s4[dframe_s4$sleep == 1, "ankle_vm"])
      
      summ_all[i, "s4_ankle_vm_per_wear_min"] <- summ_all[i, "s4_ankle_vm"] / (summ_all[i, "s4_wear_sec"] / 60)
      summ_all[i, "s4_ankle_vm_per_awake_min"] <- summ_all[i, "s4_ankle_vm_awake"] / (summ_all[i, "s4_awake_sec"] / 60)
      summ_all[i, "s4_ankle_vm_per_sleep_min"] <- summ_all[i, "s4_ankle_vm_sleep"] / (summ_all[i, "s4_sleep_sec"] / 60)
    }
    
    
    
    
    
    
    
    
    
    
    
  }
  
  
}

summ_all_1 <- summ_all[!is.na(summ_all$id),]
write.csv(summ_all_1, "summ_data_2.csv", row.names = F)
