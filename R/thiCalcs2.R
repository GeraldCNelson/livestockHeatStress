# calculate animal and human climate stress
{
  source("R/ISIMIPconstants.R")
  source("R/ISIMIPspatialConstants.R")
  #  source("R/globallyUsed.R")
  library(Rcpp)
  library(data.table)
  sourceCpp("R/cpp/THI_functions_combined.cpp")
 library("crayon")
  locOfClimFiles <- "/Volumes/ExtremeSSD2/ISIMIP/cmip6/"
  
  # constants, THI -----
  speciesChoices <- c("generic") #, "humans", "cattle",  "pigs", "chicken", "sheep") # note goats use same THI formula as cattle.
  speciesChoices <- c("generic")
  breedChoices <- c("generic_generic", "humans_adj", "bos_Taurus", "bos_Indicus",  "goats", "pigs_temp", "pigs_trop", "chicken_temp")
  breedChoices <- c("generic_generic")
  stressLevelChoices <- c("noStress", "moderateStress", "extremeStress")
  hours_worked_humans <- 8
  animalCts_global <- as.data.table(read.csv("data/animals/animalCt.csv"))
  bpList <- data.table::as.data.table(readxl::read_excel("data-raw/animals/AnimalbreakpointslistRaw.xlsx"))
  
  colorList <- (RColorBrewer::brewer.pal(4, "YlOrRd")) # yellow to red
  resultsStorage <- data.table::data.table(breed = character(), scenario = character(), startYear = numeric(), stress_cells= numeric(), no_stress_cells = numeric(), stressValue = numeric())
  
  #test values, THI -----
  speciesChoice <- "humans"
  
  f_getStressValue <- function(breedChoice, stressLevel) {
    if (speciesChoice == "humans") {
      if (stressLevel == "extremeStress") stressValue <- 40
      if (stressLevel == "moderateStress") stressValue <- 60
      if (stressLevel == "noStress") stressValue <- 80
    } else {
      stressValue <- bpList[breeds == breedChoice, get(stressLevel)]
    }
    return(stressValue)
  }
  
  f_readRast_thi_ensemble <- function(modelChoice, breedChoice, k, l) {
    yearSpan <- paste0(l, "_", l + yearRange)
    fileName_in <- paste0(locOfDataFiles_THI, "thi.", breedChoice, "_", modelChoice, "_", k,  "_", yearSpan, ".tif")
    print(paste0("breedChoice: ", breedChoice, ", k: ", k, ", modelChoice: ", modelChoice, ", fileName in: ", fileName_in))
    r <- rast(fileName_in)
    indices <- seq(from = 1, to = nlyr(r), 1)
    indices <- paste0("X", as.character(indices))
    names(r) <- indices
    r
  }
  
  # THI species functions moved to cpp file THI_functions_combined.cpp
  
  f_resultsStorage <- function(r_mean, stressValue, k, l, breedChoice) {
    stress_cells <- sum(r_mean < stressValue, na.rm = FALSE) # count of days at each location where r_out is less than stressValue
    no_stress_cells <- sum(r_mean >= stressValue, na.rm = FALSE) # count of days at each location where r_out is greater than or equal to stressValue
    stress_cells_sum <- global(stress_cells, sum, na.rm = TRUE) # total number of days below stressValue
    no_stress_cells_sum <- global(no_stress_cells, sum, na.rm = TRUE) # total number of days at or stressValue
    ratio <- stress_cells_sum/no_stress_cells_sum
    results <- c(breedChoice, k, l, stress_cells_sum, no_stress_cells_sum, stressValue)
  }
  
  f_THI <- function(k, l, yearRange) {
    # k, l, i, yearRange are set outside the function, in THI scenarios and THI historical
    yearSpan <- paste0(l, "_", l + yearRange)
    startDate <- paste0(l, "-01-01"); endDate <- paste0(l + yearRange, "-12-31")
    indices <- seq(as.Date(startDate), as.Date(endDate), 1)
    indices <- paste0("X", as.character(indices))
    indices_day <- format(as.Date(indices, format = "X%Y-%m-%d"), format = "%j") # %j is day of the year
    indices_day <- as.numeric(indices_day)
    
    for (modelChoice in modelChoices) {
      modelChoice_lower <- tolower(modelChoice)
      fileName_rh_mean <- paste0(locOfClimFiles, "mean_daily/mean_daily", "_hurs_", modelChoice_lower, "_", k,  "_",yearSpan, ".tif")
      fileName_tmax_mean <- paste0(locOfClimFiles,  "mean_daily/mean_daily", "_tasmax_", modelChoice_lower, "_", k, "_", yearSpan, ".tif")
      fileName_tmin_mean <- paste0(locOfClimFiles,  "mean_daily/mean_daily", "_tasmin_", modelChoice_lower, "_", k, "_", yearSpan, ".tif")
      fileName_tas_mean <- paste0(locOfClimFiles,  "mean_daily/mean_daily", "_tas_", modelChoice_lower, "_", k, "_", yearSpan, ".tif")
      
      rh.mean <- rast(fileName_rh_mean)
      cat(red(paste0("rh.mean min: ", round(min(minmax(rh.mean)), 2), ", rh.mean max: ", round(max(minmax(rh.mean)), 2), ", rh file: ", fileName_rh_mean, "\n")))
      tmax.mean <- rast(fileName_tmax_mean)
      cat(red(paste0("tmax.mean min: ", round(min(minmax(tmax.mean)), 2), ", tmax.mean max: ", round(max(minmax(tmax.mean)), 2), ", tmax file: ", fileName_tmax_mean, "\n")))
      tmin.mean <- rast(fileName_tmin_mean)
      cat(red(paste0("tmin.mean min: ", round(min(minmax(tmin.mean)), 2), ", tmin.mean max: ", round(max(minmax(tmin.mean)), 2), ", tmin file: ", fileName_tmin_mean, "\n")))
      tas.mean <- rast(fileName_tas_mean)
      cat(red(paste0("tas.mean min: ", round(min(minmax(tas.mean)), 2), ", tas.mean max: ", round(max(minmax(tas.mean)), 2), ", tas file: ", fileName_tas_mean, "\n")))
      if ((max(minmax(tas.mean)) | max(minmax(tmin.mean))) > max(minmax(tmax.mean))) stop("something is wrong!!")
      comb <- sds(rh.mean, tmax.mean) # the combo used for most animals
      for (breedNum in 1:length(breedChoices)) {
        speciesName <- unlist(strsplit(breedChoices[breedNum], "_"))[1]
        breedName <- breedChoices[breedNum]
        fileName_out <- paste0(locOfDataFiles_THI, "thi.", breedName, "_",  modelChoice, "_", k, "_", yearSpan, ".tif")
        print(paste0("fileName out: ", fileName_out))
        funName <- paste0("f_THI_", breedName) ## identifies the function to be used.
        print(paste0("function name: ", funName))
        if (speciesName == "humans") {
          tas.mean <- rast(fileName_tas_mean)
          comb_humans <- sds(rh.mean, tas.mean)
          print(system.time(r_out <- lapp(comb_humans, fun = funName, hours_worked_humans, filename = fileName_out, overwrite = TRUE, wopt = woptList)))
          cols = rev(colorList)
        } 
        if (speciesName %in% c("generic", "bos", "chicken", "sheep", "pigs", "cattle")) { 
          cols <- colorList
          print(system.time(r_out <- lapp(comb, fun = funName, filename = fileName_out, overwrite = TRUE, wopt = woptList)))
        } 
        gc()
        
        # display results from each breed
        print(paste0("--------", breedName, "--------"))
        print(r_out)
        cat(red(paste0("speciesName: ", breedName, ", ssp: ", k, ", start year: ", l, ", model: ", modelChoice, "\n\n")))
        cat(paste0(red(modelChoice, ", THI min: ", round(min(minmax(r_out)), 2), ",  THImax: ", round(max(minmax(r_out)), 2), ", THI out file: ", fileName_out), "\n\n" ))
        breaks = c(bpList[breeds %in% breedName, zeroLevel], bpList[breeds %in% breedName, noStress], bpList[breeds %in% breedName, moderateStress],bpList[breeds %in% breedName, extremeStress], max(minmax(r_out)))
        if (speciesName == "humans") {
          breaks = c(0, 40, 60, 80, 100)
        }
        breaks <- round(breaks, 0)
        extremeStress <- f_getStressValue(breedChoice, stressLevel= "extremeStress")
        #       if (!breedName == "humans_adj") extremeStress <- bpList[breeds %in% breedName, extremeStress]
        plot(r_out, 1, main = paste0(breedName, ", ", modelChoice, ", ", k, ", ", l, ", Jan 1, extreme stress value: ", extremeStress), ylim = c(-60, 90), range = c(0, 100), col = cols, breaks = breaks, axes = FALSE) #breaks = c(0, 40, 60, 80, 100)
        
        # write out breed-specific THI values
        print(paste0("breed: ", breedName))
        if (speciesName == "bos") {
          # cattle-specific breeds and goats
          for (breedName in c("bos_Indicus", "bos_Taurus", "goats")) {
            fileName_out_breed <- paste0(locOfDataFiles_THI, "thi.", breedName, "_",  modelChoice, "_", k, "_", yearSpan, ".tif")
            file.copy(from = fileName_out, to = fileName_out_breed)
          }
        }
        if (speciesName == "sheep") {
          for (breedName in c("sheep_trop", "sheep_temp")) {
            fileName_out_breed <- paste0(locOfDataFiles_THI, "thi.", breedName, "_",  modelChoice, "_", k, "_", yearSpan, ".tif")
            file.copy(from = fileName_out, to = fileName_out_breed)
          }
        }
        if (speciesName == "chicken") {
          for (breedName in c("chicken_trop", "chicken_temp")) {
            fileName_out_breed <- paste0(locOfDataFiles_THI, "thi.", breedName, "_",  modelChoice, "_", k, "_", yearSpan, ".tif")
            file.copy(from = fileName_out, to = fileName_out_breed)
          }
        }
        if (speciesName == "pigs") {
          for (breedName in c("pigs_trop", "pigs_temp")) {
            fileName_out_breed <- paste0(locOfDataFiles_THI, "thi.", breedName, "_",  modelChoice, "_", k, "_", yearSpan, ".tif")
            print(paste0("fileName_out: ", fileName_out, ", fileName_out_breed: ", fileName_out_breed))
            file.copy(from = fileName_out, to = fileName_out_breed)
          } 
          if (speciesName == "humans") {
            fileName_out_breed <- paste0(locOfDataFiles_THI, "thi.", breedName, "_",  modelChoice, "_", k, "_", yearSpan, ".tif")
            print(paste0("fileName_out: ", fileName_out, ", fileName_out_breed: ", fileName_out_breed))
            file.copy(from = fileName_out, to = fileName_out_breed)
          }
        }
      }
    }
  }
  
  f_extremeTHI <- function(cellVector, stressValue, breedChoice) {
    if (breedChoice == "humans_adj") {
      extremeCt <- sum(cellVector <= stressValue, na.rm = FALSE)
    } else {extremeCt <- sum(cellVector >= stressValue, na.rm = FALSE)
    }
    return(extremeCt) 
  }
  
  f_StressCt_setup <- function(k, l, colorList) {
    yearSpan <- paste0(l, "_", l + yearRange)
    for (breedChoice in breedChoices) {
      speciesChoice <- strsplit(breedChoice,"_")[[1]][1]
      fileName_in <- paste0(locOfDataFiles_THI, "ensemble_thi.", breedChoice, "_", k, "_", yearSpan, ".tif")
      r <- rast(fileName_in)
      maxVal <- round(max(minmax(r)), 3); minVal <- round(min(minmax(r)), 3)
      for (stressLevelChoice in stressLevelChoices) { # eg noStress, moderateStress, etc
        stressValue <- f_getStressValue(breedChoice, stressLevel = stressLevelChoice) 
        print(paste0("breed: ", breedChoice, ", minVal: ", minVal, ", maxVal: ", maxVal, ", stress value: ", stressValue, ", fileName in: ", fileName_in))
        fileName_out <- paste0(locOfDataFiles_THI, "stressCt_", breedChoice, "_", k, "_stressLevel_",  stressValue, "_", yearSpan,  ".tif")
        print(system.time(extremeCt <- app(r, fun = f_extremeTHI, stressValue, breedChoice, filename = fileName_out, overwrite = TRUE, wopt = woptList)))
        print(paste0("fileName in: ", fileName_in, ", fileName out: ", fileName_out))
        titleText <- paste0(breedChoice, ", days above stress value of ", stressValue, ", \nensemble means, ", k, ", ", gsub("_", "-", yearSpan))
         cols = colorList
        if (breedChoice == "humans_adj") {
          cols <- colorList
          titleText <- paste0("Days where PWC is less than ", stressValue, " percent, \nensemble means, ", k, ", ", gsub("_", "-", yearSpan))
        }
        # breaks = seq(from = 1, to = max(minmax(extremeCt)), length.out = 5)
        breaks = seq(from = 1, to = 366, length.out = 4)
        breaks <- round(breaks, 0)
        breaks <- c(0, breaks) # add to give 0 - 1 a color
        print(paste0("breaks: ", breaks))
        print(r)
        print(extremeCt)
        print(paste0("stress value : ", stressValue))
        plot(extremeCt, 1, main = titleText, ylim = c(-60, 90), range = c(0, 100), col = cols, breaks = breaks, axes = FALSE, colNA = "gray")
        print("-----------------------------------")
      }
    }
  }
  
  f_extremeStressGenericCt <- function(k, l) { #, colorList
    stressValueList <- c(80, 81, 84, 85, 86, 89, 90, 91, 92, 93, 94)
    yearSpan <- paste0(l, "_", l + yearRange)
    #for (m in breedChoice) {
    breedChoice <- "generic_generic"
    speciesChoice <- strsplit(breedChoice,"_")[[1]][1]
    
    #      for (s in stressLevelChoices) { # eg noStress, moderateStress, etc
    s <- "extremeStress"
    #       stressValue <- f_getStressValue(m, stressLevel = s)
    fileName_in <- paste0(locOfDataFiles_THI, "ensemble_thi.", breedChoice, "_", k, "_", yearSpan, ".tif")
    r <- rast(fileName_in)
    maxVal <- round(max(minmax(r)), 3)
    minVal <- round(min(minmax(r)), 3)
    for (stressValue in stressValueList) {
      print(paste0("species: ", speciesChoice, ", minVal: ", minVal, ", maxVal: ", maxVal, ", stress value: ", stressValue, ", fileName in: ", fileName_in))
      fileName_out <- paste0(locOfDataFiles_THI, "stressCt_", breedChoice, "_", k, "_stressLevel_",  stressValue, "_", yearSpan,  ".tif")
      print(system.time(extremeCt <- app(r, fun = f_extremeTHI, stressValue, breedChoice, filename = fileName_out, overwrite = TRUE, wopt =woptList)))
      print(paste0("fileName in: ", fileName_in, ", fileName out: ", fileName_out))
    }
  }
  
  f_ensemble_THI <- function(k, l, stressLevelChoices) {
    yearSpan <- paste0(l, "_", l + yearRange)
    storageTemp <- resultsStorage[0] 
    for (breedChoice in breedChoices) {
      print(paste0("------------"))
      x <- lapply(modelChoices, f_readRast_thi_ensemble, breedChoice, k, l)
      r <- rast(x)
      indices_day <- rep(seq(1, nlyr(x[[1]]), 1), 5) # 5 is number of models; if omitted should get the same result
      #     maxVal <- round(max(minmax(r)), 2); minVal <- round(min(minmax(r)), 2)
      fileName_out <- paste0(locOfDataFiles_THI, "ensemble_thi.", breedChoice, "_", k, "_", yearSpan, ".tif")
      cat(paste0(red("breed: ", breedChoice, ", ensemble ssp: ", k, ", start year: ", l , ", fileName out: ", fileName_out)))
      print(system.time(r.mean <- tapp(r, indices_day, fun = "mean", na.rm = TRUE, filename = fileName_out, overwrite = TRUE, wopt = woptList)))
      names(r.mean) <- gsub("X", "Day ", names(r.mean))
      cat(paste0(red("ensemble mean min: ", round(min(minmax(r.mean)), 2), ", ensemble mean max: ", round(max(minmax(r.mean)), 2), ", ensemble mean out file: ", fileName_out), "\n\n"))
      for (stressLevel in stressLevelChoices) {
        stressValue <- f_getStressValue(breedChoice, stressLevel) 
        cat(paste0(red("breed: ", breedChoice, ", ensemble ssp: ", k, ", start year: ", l,", stress value: ", stressValue , ", fileName out: ", fileName_out), "\n\n"))
        storageTemp <- data.table::rbindlist(list(storageTemp, f_resultsStorage(r.mean, stressValue, k, l, breedChoice)), use.names = FALSE)
      }
      return(storageTemp)
    }
  }
  
  f_mean_THI <- function(cellVector) { # get mean over all days in a year for each pixel
    if (is.nan(cellVector[1])) {return(NA)}
    return(mean(cellVector, na.rm = TRUE))
  }
  
  f_prepareR <- function() {
    speciesChoice <- strsplit(breedChoice,"_")[[1]][1]
    yearSpan <- paste0(l, "_", l + yearRange)
    fileName_in <- paste0(locOfDataFiles_THI, "extremeCt.", speciesChoice, "_", k, "_", yearSpan, ".tif")
    print(paste0("fileName in: ", fileName_in))
    r <- rast(fileName_in)
    names(r) <- "value"
    return(r)
  }
  
  # function to fill in dt_stressCts
  f_dt_stressCts <- function(r, dt_stressCts, r_mask, stressValue) {
    r_mask_scen <- r_mask
    r_mask_scen[r_mask_scen < maskMin] <- NA
    r_mask_scen[r_mask_scen > 0] <- 1
    r <- mask(r, r_mask_scen)
    r_ext <- r
    r_ext[r_ext < maxCount] <- NA
    r_ext_mask <- mask(r_mask, r_ext)
    totalNum_extreme <- global(r_ext_mask, fun = "sum", na.rm = TRUE)
    totalNum <- 6.512 * 10^9 # world population in 2005 is 6.512 billion
    if (!speciesChoice == "humans") totalNum <- animalCts_global[species == speciesChoice, ct]
    dt_stressCts <- rbindlist(list(dt_stressCts, list(speciesChoice, k, l, as.numeric(totalNum_extreme), totalNum, maxCount, stressValue)), use.names = FALSE)
    return(dt_stressCts)
  }
  
  f_thi_graphing <- function(m, r, r_mask, maskMin, col, dt_stressCts) {
    print("--------------------------")
    speciesChoice <- strsplit(breedChoice,"_")[[1]][1]
    extremeStress <- dt_stressCts[species == speciesChoice, extremeStress]
    totalNum<- dt_stressCts[species == speciesChoice, totalCts]
    totalNum_extreme<- dt_stressCts[species == speciesChoice, stressCts]
    totalNum_extreme <- totalNum_extreme/1000000 # convert to millions
    totalNum <- totalNum/1000000 # convert to millions
    ratio_extreme <- 100 * totalNum_extreme/totalNum # convert to percent
    totalNum_extreme <- round(totalNum_extreme, 0)
    totalNum <- round(totalNum, 0)
    ratio_extreme <- round (ratio_extreme, 1)
    print(paste0("share of ", speciesChoice, " in regions with extreme stress days greater than ", maxCount, " is ", round(100 * totalNum_extreme/totalNum, 2), " %."))
    print(paste0("species: ", speciesChoice, ", min r: ", min(minmax(r)), ", max r: ", max(minmax(r))))
    if (is.na( max(minmax(r)))) stop("max is na")
    # browser()
    # # crop to eliminate Antarctica and project to Robinson
    # r_mask <- crop(r, extent_noAntarctica)
    # r_maskRob <- project(r_mask, crsRob)
    r <- project(r, crsRob)
    r_mask <- project(r_mask, crsRob)
    r_mask_df <- as.data.frame(r_mask, xy = TRUE)
    names(r_mask_df) <-   c("x", "y", "value")
    
    r_df <- as.data.frame(r, xy = TRUE)
    #   r_df$value[r_df$value > maxVal] <- maxVal #set values > maxVal to maxVal
    r_df_mod <- r_df %>%
      mutate(value_2 = cut(value, breaks = custom_bins)) %>%
      group_by(value_2)
    
    #  print(paste0("min r proj: ", min(minmax(r)), ", max r proj: ", max(minmax(r)), ", rname: ", paste0("r_", k, "_", l)))
    titleText <- paste0(m, ", days above extreme stress value, ", k, ", ", gsub("_", "-", yearSpan))
    if (m == "humans") {
      titleText <- paste0(m, ", days with PWC below ", extremeStress, " percent, " , k, ", ", gsub("_", "-", yearSpan))
    }
    caption <- paste0("The extreme stress value for ", speciesChoice, " is ", extremeStress, ". Stress locations are where the species was raised in the early 21st century. \nOf early century species numbers, ", ratio_extreme, "% are in locations with at least ", maxCount, " days with extreme stress during this period.")
    if (m == "humans") {
      caption <- paste0("The extreme stress value for reduction in physical work capacity (PWC) is ", extremeStress, " %.")
    }
    g <- ggplot() +
      geom_tile(data = r_df_mod, aes(x, y, fill = value_2)) +
      #      scale_fill_discrete(colors = colorList, drop = FALSE, na.value = 'grey95') +
      scale_fill_manual(values = col, drop = FALSE, na.value = 'grey95') + # the na.value doesn't work yet. See https://stackoverflow.com/questions/45144630/scale-fill-manual-define-color-for-na-values/45147172 for a possible solution
      labs(title = titleText, fill = legendTitle, x = "", y = "", caption = caption) + 
      theme_bw()  +
      #      theme(plot.title = element_text(size = 12, hjust = 0.5), plot.caption = element_text(hjust = 0, size = 8)) +
      theme(
        legend.text.align = 1,
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 12, hjust = 0.5),
        plot.caption = element_text(hjust = 0, vjust = 7.0, size = 8)
      ) +
      geom_sf(data = coastline_cropped , 
              color = "black", size = 0.1, stat = "sf", fill = NA,
              position = "identity")
    #    if (! m == "humans") g <- g + geom_tile(data = r_mask_df, aes(x, y, fill = value), alpha = .2, show.legend = FALSE) # +
    
    print(g)
    fileName_out <- paste0(lofOfGraphicsFiles, "THI/THIextremeCt.", speciesChoice, "_", k, "_", yearSpan, ".png")
    
    ggsave(filename = fileName_out, plot = g, device = "png", width = 6, height = 6, units = "in", dpi = 300)
    knitr::plot_crop(fileName_out) # gets rid of margins around the plot
    print(paste0("fileName out: ", fileName_out))
  }
  
  f_thi_runs_graphing <- function(m, r, col, dt_stressCts, run_length) {
    print("--------------------------")
    speciesChoice <- strsplit(breedChoice,"_")[[1]][1]
    print(paste0("species: ", speciesChoice, ", min r: ", min(minmax(r)), ", max r: ", max(minmax(r))))
    if (is.na(max(minmax(r)))) stop("max is na")
    extremeStress <- dt_stressCts[species == speciesChoice, extremeStress]
    # maxCount<- dt_stressCts[species == speciesChoice, countDays]
    
    r <- project(r, crsRob)
    r_runsCt <- r$runsCt
    r_runLength <- r$runLength
    r_runsCt_df <- as.data.frame(r_runsCt, xy = TRUE)
    names(r_runsCt_df) <-   c("x", "y", "value")
    r_runLength_df <- as.data.frame(r_runLength, xy = TRUE)
    names(r_runLength_df) <-   c("x", "y", "value")
    
    # use this code to get binned legend
    # r_df_mod <- r_df %>%
    #   mutate(value_2 = cut(value, breaks = custom_bins)) %>%
    #   group_by(value_2)
    
    # graphing for r_runsCt -----
    legendTitle <- "No. of runs"
    titleText <- paste0(m, ", number of ", run_length, " day runs above extreme stress value, \n", k, ", ", gsub("_", "-", yearSpan))
    if (m == "humans") {
      titleText <- paste0(m, ", number of ", run_length, " day runs with PWC below ", extremeStress, " percent, \n" , k, ", ", gsub("_", "-", yearSpan))
    }
    caption <- paste0("The extreme stress value for ", speciesChoice, " is ", extremeStress, ". Stress locations are where the species was raised in the early 21st century.")
    if (m == "humans") {
      caption <- paste0("The extreme stress value for reduction in physical work capacity (PWC) is ", extremeStress, " %.")
    }
    colorList0 <- c("gray95", colorList)
    g <- ggplot() +
      #      geom_tile(data = r_df_mod, aes(x, y, fill = value_2)) + # use for binned legend
      geom_tile(data = r_runsCt_df, aes(x, y, fill = value)) +
      scale_fill_gradientn(colors = colorList0,  na.value = 'grey50') +
      #scale_fill_manual(values = col, drop = FALSE, na.value = 'grey95') + # the na.value doesn't work yet. See https://stackoverflow.com/questions/45144630/scale-fill-manual-define-color-for-na-values/45147172 for a possible solution
      labs(title = titleText, fill = legendTitle, x = "", y = "", caption = caption) + 
      theme_bw()  +
      #      theme(plot.title = element_text(size = 12, hjust = 0.5), plot.caption = element_text(hjust = 0, size = 8)) +
      theme(
        legend.text.align = 1,
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 12, hjust = 0.5),
        plot.caption = element_text(hjust = 0, vjust = 7.0, size = 8)
      ) +
      geom_sf(data = coastline_cropped , 
              color = "black", size = 0.1, stat = "sf", fill = NA, position = "identity")
    #    if (! m == "humans") g <- g + geom_tile(data = r_mask_df, aes(x, y, fill = value), alpha = .2, show.legend = FALSE) # +
    
    print(g)
    fileName_out <- paste0(lofOfGraphicsFiles, "THI/runsCt.", speciesChoice, "_", k, "_", yearSpan, ".png")
    
    ggsave(filename = fileName_out, plot = g, device = "png", width = 6, height = 6, units = "in", dpi = 300)
    knitr::plot_crop(fileName_out) # gets rid of margins around the plot
    print(paste0("fileName out: ", fileName_out))
    
    # graphing for r_runLength-----
    legendTitle <- "Longest run length"
    titleText <- paste0(m, ", length of longest run with at least ", run_length, " \ndays above extreme stress value, ", k, ", ", gsub("_", "-", yearSpan))
    if (m == "humans_adj") {
      stressValue <- f_getStressValue(m, stressLevel)
      titleText <- paste0(m, ", length of longest run with at least ", run_length, " \ndays where PWC is below ", stressValue, " percent, " , k, ", ", gsub("_", "-", yearSpan))
    }
    caption <- paste0("The extreme stress value for ", speciesChoice, " is ", stressValue, ". Stress locations are where the species was raised in the early 21st century.")
    if (m =="humans_adj") {
      caption <- paste0("The extreme stress value for reduction in physical work capacity (PWC) is ", stressValue, " %.")
    }
    colorList0 <- c("gray95", colorList)
    g <- ggplot() +
      #      geom_tile(data = r_df_mod, aes(x, y, fill = value_2)) + # use for binned legend
      geom_tile(data = r_runLength_df, aes(x, y, fill = value)) +
      scale_fill_gradientn(colors = colorList0,  na.value = 'grey50') +
      #scale_fill_manual(values = col, drop = FALSE, na.value = 'grey95') + # the na.value doesn't work yet. See https://stackoverflow.com/questions/45144630/scale-fill-manual-define-color-for-na-values/45147172 for a possible solution
      labs(title = titleText, fill = legendTitle, x = "", y = "", caption = caption) + 
      theme_bw()  +
      #      theme(plot.title = element_text(size = 12, hjust = 0.5), plot.caption = element_text(hjust = 0, size = 8)) +
      theme(
        legend.text.align = 1,
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 12, hjust = 0.5),
        plot.caption = element_text(hjust = 0, vjust = 7.0, size = 8)
      ) +
      geom_sf(data = coastline_cropped , 
              color = "black", size = 0.1, stat = "sf", fill = NA, position = "identity")
    #    if (! m == "humans") g <- g + geom_tile(data = r_mask_df, aes(x, y, fill = value), alpha = .2, show.legend = FALSE) # +
    
    print(g)
    fileName_out <- paste0(lofOfGraphicsFiles, "THI/runLength.", speciesChoice, "_", k, "_", yearSpan, ".png")
    
    ggsave(filename = fileName_out, plot = g, device = "png", width = 6, height = 6, units = "in", dpi = 300)
    knitr::plot_crop(fileName_out) # gets rid of margins around the plot
    print(paste0("fileName out: ", fileName_out))
  }
  
  f_thi_means_graphing <- function(m, r, col) {
    print("--------------------------")
    print(paste0("species: ", speciesChoice, ", min r: ", min(minmax(r)), ", max r: ", max(minmax(r))))
    if (is.na(max(minmax(r)))) stop("max is na")
    
    r <- project(r, crsRob)
    r_df <- as.data.frame(r, xy = TRUE)
    names(r_df) <-   c("x", "y", "value")
    
    # use this code to get binned legend
    # r_df_mod <- r_df %>%
    #   mutate(value_2 = cut(value, breaks = custom_bins)) %>%
    #   group_by(value_2)
    
    # graphing for r_runsCt -----
    legendTitle <- "Global daily mean"
    titleText <- paste0(m, ", mean of daily THI, \n", k, ", ", gsub("_", "-", yearSpan))
    if (m == "humans") {
      titleText <- paste0(m, ", mean of daily PWC, \n" , k, ", ", gsub("_", "-", yearSpan))
    }
    # caption <- paste0("The extreme stress value for ", speciesChoice, " is ", extremeStress, ". Stress locations are where the species was raised in the early 21st century.")
    # if (m == "humans") {
    #   caption <- paste0("The extreme stress value for reduction in physical work capacity (PWC) is ", extremeStress, " %.")
    # }
    colorList0 <- c("gray95", colorList)
    g <- ggplot() +
      #      geom_tile(data = r_df_mod, aes(x, y, fill = value_2)) + # use for binned legend
      geom_tile(data = r__df, aes(x, y, fill = value)) +
      scale_fill_gradientn(colors = colorList0,  na.value = 'grey50') +
      #scale_fill_manual(values = col, drop = FALSE, na.value = 'grey95') + # the na.value doesn't work yet. See https://stackoverflow.com/questions/45144630/scale-fill-manual-define-color-for-na-values/45147172 for a possible solution
      labs(title = titleText, fill = legendTitle, x = "", y = "") + #, caption = caption) + 
      theme_bw()  +
      #      theme(plot.title = element_text(size = 12, hjust = 0.5), plot.caption = element_text(hjust = 0, size = 8)) +
      theme(
        legend.text.align = 1,
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 12, hjust = 0.5) #,
        #     plot.caption = element_text(hjust = 0, vjust = 7.0, size = 8)
      ) +
      geom_sf(data = coastline_cropped , 
              color = "black", size = 0.1, stat = "sf", fill = NA, position = "identity")
    #    if (! m == "humans") g <- g + geom_tile(data = r_mask_df, aes(x, y, fill = value), alpha = .2, show.legend = FALSE) # +
    
    print(g)
    fileName_out <- paste0(lofOfGraphicsFiles, "THI/mean.", speciesChoice, "_", k, "_", yearSpan, ".png")
    
    ggsave(filename = fileName_out, plot = g, device = "png", width = 6, height = 6, units = "in", dpi = 300)
    knitr::plot_crop(fileName_out) # gets rid of margins around the plot
    print(paste0("fileName out: ", fileName_out))
  }
}

#THI model-specific, scenarios -----
for (k in sspChoices) {
  for (l in startYearChoices) {
    f_THI(k, l, yearRange)
  }
}

#THI model-specific, historical -----
k <- "historical"
l = 1991
f_THI(k, l, yearRange)

# combine all the spatrasters by breed and model for the time period and then take the mean across that combo

# THI ensemble means, historical -----
k <- "historical"
l <- 1991
dt_out <- f_ensemble_THI(k, l, stressLevelChoices)
resultsStorage <- data.table::rbindlist(list(resultsStorage, dt_out), use.names = FALSE)

# THI ensemble means, scenarios -----

for (k in sspChoices) {
  for (l in startYearChoices) {
    dt_out <- f_ensemble_THI(k, l, stressLevelChoices)
    resultsStorage <- data.table::rbindlist(list(resultsStorage, dt_out), use.names = FALSE)
  }
}

fileName_out <- paste0(locOfDataFiles_THI, "results/resultsStorage_all.csv")
write.csv(resultsStorage, fileName_out, row.names = FALSE)

for (k in sspChoices) {
  for (l in startYearChoices) {
    f_extremeStressGenericCt(k, l)
  }
}

k <- "historical"
l <- 1991
f_extremeStressGenericCt(k, l)

# mpegs, ensemble means -----

# for (k in sspChoices) {
#   for (l in startYearChoices) {
#     yearSpan <- paste0(l, "_", l + yearRange)
#     #     yearnumberRange <- seq(l, (l + yearRange), 1)
#     
#     # cattle and goats have the same values, but mpegs are created separately
#     for (speciesChoice in speciesChoices) {
#       cols <- colorList
#       if (speciesChoice == "humans") cols = rev(colorList)
#       #      browser()      
#       extremeStress <- f_getStressValue(m, stressLevel) # just to have something to put display
#       print(paste0("species choice: ", speciesChoice))
#       fileName_in <- paste0(locOfDataFiles_THI, "ensemble_thi.", speciesChoice, "_", k, "_", yearSpan, ".tif")
#       r <- rast(fileName_in)
#       maxVal <- round(max(minmax(r)), 3)
#       minVal <- round(min(minmax(r)), 3)
#       print(paste0("raster r, ", "species: ", speciesChoice, ", minVal: ", minVal, ", maxVal: ", maxVal, ", extreme stress: ", extremeStress, ", fileName in: ", fileName_in))      
#       # set up names with Jan-1. Need to set up the endDate_year with a leap year such as 2040
#       startDate_year <- paste0(2040, "-01-01")
#       endDate_year <- paste0(2040, "-12-31")
#       indices <- seq(as.Date(startDate_year), as.Date(endDate_year), 1)
#       indices_day <- format(indices, "%b %d")
#       indices_date <- format(indices, "%j")
#       #      names(r) <- gsub("X", "Day.", names(r))
#       videoName_out <- paste0(lofOfGraphicsFiles, "THI/ensemble_thi.", speciesChoice, "_", k, "_", yearSpan, ".mp4")
#       print(paste0("video fileName out: ", videoName_out))
#       #     x.cv <- tapp(x, indices_date, fun = cv, na.rm = TRUE)
#       system.time(r.mean <- tapp(r, indices_date, fun = mean, na.rm = TRUE)) #, filename = fileName_out, overwrite = TRUE, wopt = woptList))
#       print(r.mean)
#       names(r.mean) <- gsub("X", "Day ", names(r.mean))
#       title_animate <- paste0("THI value, ", speciesChoice, " ", k,  ", ", gsub("_", "-", yearSpan), ", extreme stress value: ", extremeStress)
#       breaks = c(1, bpList[species %in% speciesChoice, noStress], bpList[species %in% speciesChoice, moderateStress],bpList[species %in% speciesChoice, extremeStress], max(minmax(r.mean)))
#       if (m == "humans") {
#         title_animate <- paste0("PWC value (%), ", k, ", ", gsub("_", "-", yearSpan), ", extreme stress PWC value is less than: ", extremeStress)
#         range = c(0, 100)
#         breaks = c(1, 40, 60, 80, 100)
#       }
#       print(paste0("raster r.mean, ", "species: ", speciesChoice, ", minVal: ", minVal, ", maxVal: ", maxVal, ", extreme stress: ", extremeStress, ", fileName in: ", fileName_in))      
#       breaks <- round(breaks, 0)
#       print(paste0("breaks: ", breaks))
#       names(r.mean) <- paste0(title_animate, ", ", indices_day)
#       rangeMaximum <- max(minmax(r.mean))
#       animation::saveVideo(animate(r.mean, n=1, ylim = c(-60, 90), range = c(0, rangeMaximum,  main = "test",), pause = .001, sub = title_animate, col = cols, breaks = breaks, axes = FALSE), ani.height = 800, ani.width = 1200, video.name = videoName_out)
#     }
#   }
# }

#speciesChoice <- "cattle"
#  stress level counts, scenarios -------
for (k in sspChoices) {
  for (l in startYearChoices) {
    f_StressCt_setup(k, l, colorList)
  }
}

#  stress level counts, historical -------
k <- "historical"
l <- 1991
f_StressCt_setup(k, l, colorList)

{#  do graphics -----
  library(ggplot2)
  library(RColorBrewer)
  #library(rworldmap)
  library(maps)
  #remotes::install_github("ropensci/rnaturalearthhires") need to do once to get the library from github
  #library(rnaturalearthhires)
  library(ggspatial)
  library(readxl)
  library(sf)
  library(dplyr)
  library(magick)
  legendTitle <- "Days"
  # test values
  m = "humans"
}

# speciesChoice <- "humans"
#graphics setup, extreme stress-----

# create blank table to hold animal numbers in high stress areas
dt_stressCts <- data.table(species = character(), ssp  = character(), startYear  = character(), stressCts = numeric(), totalCts = numeric(), countDays = numeric(), extremeStress = numeric())

maxCount <- 30 # used to identify regions with extreme stress value more than maxCount days

for (speciesChoice in speciesChoices) {
  print(paste0("------", speciesChoice, "------"))
  stressValue <- f_getStressValue(speciesChoice, "extremeStress") 
  
  # create mask by number of animals in cell of fileName_r_mask
  fileName_r_mask <- paste0("data/animals/raster_ct_", speciesChoice, ".tif") # already cropped Antarctica off
  r_mask <- rast(fileName_r_mask)
  maskMin <- switch(
    speciesChoice,
    "bos_Taurus" = 1000,
    "bos_Indicus" = 1000,
    "humans" = 0,
    "pigs" = 1000,
    "chicken" = 10000,
    "goat" = 1000,
    "sheep" = 1000
  )
  
  #  maxVal <- max(minmax(r_ssp585_2081))  # assumes the maximum value in all years is found at the end of the century in the scenario ssp585
  #  custom_bins <- round(seq.int(from = 0, to = maxVal, length = 4))
  # alternate approach to custom_bins, choose relevant number of days since should be the same for all THIs
  custom_bins <- c(1, 25, 50, 100, 366)
  cols <- colorList
  
  #  graphics, scenarios  -----
  breedChoice <- "humans_adj"
  stressLevel <- "extremeStress"
  stressValue <- f_getStressValue(breedChoice, stressLevel)
  for (k in sspChoices) {
    for (l in startYearChoices) {
      yearSpan <- paste0(l, "_", l + yearRange)
      fileName_in <- paste0(locOfDataFiles_THI, "stressCt_", breedChoice, "_", k, "_stressLevel_", stressValue, "_", yearSpan, ".tif")
      print(paste0("fileName in: ", fileName_in))
      r <- rast(fileName_in)
      names(r) <- "value"
      dt_stressCts <- f_dt_stressCts(r, dt_stressCts, r_mask, stressValue)
      f_thi_graphing(speciesChoice, r, r_mask, maskMin, col, dt_stressCts)
      print("-------------------------")
    }
  }
  
  # graphics, historical -----
  k <- "historical"
  l <- 1991
  yearSpan <- paste0(l, "_", l + yearRange)
  fileName_in <- paste0(locOfDataFiles_THI, "stressCt_", breedChoice, "_", k, "_stressLevel_", stressValue, "_", yearSpan, ".tif")
  print(paste0("fileName in: ", fileName_in))
  r <- rast(fileName_in)
  names(r) <- "value"
  dt_stressCts <- f_dt_stressCts(r, dt_stressCts, r_mask, stressValue)
  f_thi_graphing(speciesChoice, r, r_mask, maskMin, col, dt_stressCts)
  print("-------------------------")
  
  # delta graphics -----
  l = 2081
  yearSpan <- paste0(l, "_", l + yearRange)
  
  for (k in sspChoices) {
    titleText <- paste0(m, ", change in days above extreme stress value, \n early to end century ", k)
    if (m == "humans") {
      titleText <- paste0(m, ", change in days with PWC below ", extremeStress, " percent, \n early to end century ", k)
    }
    fileName_out <- paste0(lofOfGraphicsFiles, "THI/THIextremeCtDelta.", speciesChoice, "_", k, "_", yearSpan, ".png")
    r_ssp <- rast(paste0(locOfDataFiles_THI, "extremeCt.", speciesChoice, "_", k, "_", yearSpan, ".tif"))
    r_historical <- rast(paste0(locOfDataFiles_THI, "extremeCt.", speciesChoice, "_", "historical", "_", "1991_2010", ".tif"))
    r <- r_ssp - r_historical
    names(r) <- "value"
    r[r$value < 0] <- 0
    print(r)
    # r <- mask(r, r_maskRob)
    # r_df <- as.data.frame(r, xy = TRUE)
    r_mask_scen <- r_mask
    r_mask_scen[r_mask_scen < maskMin] <- NA
    r_mask_scen[r_mask_scen > 0] <- 1
    #   r <- eval(parse(text = paste0("r_", k, "_", speciesChoice, "_", l)))
    r <- mask(r, r_mask_scen)
    
    # get counts in extreme stress regions with more than maxCount days
    r_ext <- r >= maxCount
    r_ext_mask <- mask(r_mask, r_ext)
    totalNum_extreme <- global(r_ext_mask, fun = "sum", na.rm = TRUE)
    totalNum <- global(r_mask, fun = "sum", na.rm = TRUE)
    
    
    r <- crop(r, extent_noAntarctica)
    f_thi_graphing(m, r, r_mask, maskMin, col, dt_stressCts)
  }
}

write.csv(dt_stressCts, paste0(locOfDataFiles_THI, "stressCtsTable.csv"))

{# do ppt for THI extreme ct ensemble means -----
  library(officer)
  source("R/pptFunctions.R")
  f_extremeCtSpeciesForPptx <- function(m) {
    fileNameStart <- paste0("THIextremeCt.", speciesChoice, "_")
    fileName_in <- paste0(lofOfGraphicsFiles, "THI/", fileNameStart, k, "_", yearSpan, ".png")
    print(paste0("fileName in: ", fileName_in))
    extImg_favLocs <- external_img(src = fileName_in, width = defaultWidth, height = defaultHeight)
    my_pres <- add_slide(x = my_pres, layout = 'Title Only', master = 'Office Theme')
    my_pres <- ph_with(x = my_pres, value = extImg_favLocs, location = ph_location(left = defaultLeft, top = defaultTop, width = defaultWidth, height = defaultHeight - 0.5) )
    return(my_pres)
  }
  
  f_deltaExtremeCtSpeciesForPptx <- function() {
    print(paste0("fileName in: ", fileName_in))
    extImg_favLocs <- external_img(src = fileName_in, width = defaultWidth, height = defaultHeight)
    my_pres <- add_slide(x = my_pres, layout = 'Title Only', master = 'Office Theme')
    my_pres <- ph_with(x = my_pres, value = extImg_favLocs, location = ph_location(left = defaultLeft, top = defaultTop, width = defaultWidth, height = defaultHeight - 0.5) )
    return(my_pres)
  }
  
  thiList <- c("thi.humans", "thi.cattle", "thi.sheep", "thi.goat", "thi.yak", "thi.broiler", "thi.layer", "thi.chicken", "thi.pigs")
  thiListReduced <- thiList[!thiList %in% c("thi.yak", "thi.broiler", "thi.layer")]
  
  titleString <- paste0("Effects to 2100 of Temperature and Humidity on Productivity of Humans and ", length(thiListReduced)-1 , " Animal Species")
  contentString <- paste0("Preliminary Results: Monthly ensemble means by species for three time periods to 2100. Powerpoint produced on ", Sys.Date())
  startYearChoices_ensemble <-  c(1991, 2041, 2081) 
  
}
IntroText0 <- "The productivity of humans doing field work and animals in producing meat and milk is affected by exposure to combined high levels of temperature and humidity. For humans, the physical work capacity index (PWC) shows the percentage by which  human work capacity is reduced. The temperature and humidity index (THI) for animals is a species-specific measure of those effects with thresholds for low, medium and high negative productivity effects. "
IntroText0.5 <- "\nThe following figures present global graphics of locations where extreme stress is experienced and how many days this happens in a representative year."
IntroText1 <- "\nFor the animals graphics, extreme stress areas are where the species was raised in the early 21st century."
IntroText4 <- "\nThe THI values are averages for three 20 year periods (1991-2010, 2041-2060, and 2081-2100)."
IntroText5 <- "\nThis powerpoint presents work in progress and should not be circulated without permission."
fp_1 <- fp_text(bold = TRUE, color = "pink", font.size = 0)
fp_2 <- fp_text(bold = FALSE, font.size = 12)
fp_3 <- fp_text(italic = TRUE, color = "black", font.size = 14)

blIntro <- block_list(
  fpar(
    ftext(IntroText0, fp_2),
    ftext(IntroText0.5, fp_2)),
  fpar(ftext(IntroText1, fp_2)),
  fpar(ftext(IntroText4, fp_2)),
  fpar(ftext(IntroText5, fp_2))
)

my_pres <- read_pptx()
my_pres <- add_slide(x = my_pres, layout = 'Title Slide', master = 'Office Theme')
my_pres <- ph_with(x = my_pres, value = titleString, location = ph_location_type(type = "ctrTitle"))
my_pres <- ph_with(x = my_pres, value = contentString, location = ph_location_type(type = "subTitle"))

my_pres <- add_slide(my_pres, layout = "Title and Content", master = "Office Theme")
my_pres <-  ph_with(x = my_pres, value = "Introduction", location = ph_location_type(type = "title"))
my_pres <- ph_with(x = my_pres, value = blIntro, location = ph_location_type(type = "body") )

# do all slides for each breed in breedChoice
for (m in breedChoice) {
  m <- gsub("thi.", "", m)
  ensembleTitle <- m
  my_pres <- add_slide(x = my_pres, layout = 'Section Header', master = 'Office Theme')
  my_pres <- ph_with(x = my_pres, value = ensembleTitle, location = ph_location_type(type = "title"))
  
  # do historical first, then ssps and future periods -----
  k <- "historical"
  l <- 1991
  yearSpan <- paste0(l, "_", l + yearRange)
  my_pres <- f_extremeCtSpeciesForPptx(m)
  
  for (k in sspChoices) {
    for (l in startYearChoices) {
      yearSpan <- paste0(l, "_", l + yearRange)
      my_pres <- f_extremeCtSpeciesForPptx(m)
    }
  }
  
  # delta days, ppt -----
  ensembleTitle <- paste0(m, ", Change between early and end of the 21st century in days with extreme risk")
  if (speciesChoice == "humans")   ensembleTitle <- paste0(speciesChoice, ", Change in days with PWC less than 60% between early and end of the 21st century")
  
  my_pres <- add_slide(x = my_pres, layout = 'Section Header', master = 'Office Theme')
  my_pres <- ph_with(x = my_pres, value = ensembleTitle, location = ph_location_type(type = "title"))
  for (k in sspChoices) {
    fileName_in <- paste0(lofOfGraphicsFiles, "THI/THIextremeCtDelta.", speciesChoice, "_", k, "_", "2081_2100.png")
    my_pres <- f_deltaExtremeCtSpeciesForPptx()
  }
  my_pres <- f_addDataSlide()
}

print(my_pres, target = "presentations/cmip6/THI/damageTemp_Ensemble.pptx") %>% browseURL()


# for (j in 1:length(thiListReduced)) {
#   for (k in sspChoices) {
#     yearSpan <- paste0(l, "_", l + yearRange)
#     print(paste0("ssp choice: ", k, ", start year: ", l))
#     speciesName <- gsub("thi.", "", thiListReduced[j])
#     ensembleTitle <- paste("Ensemble Mean for", speciesName)
#     add_slide(my_pres, layout = 'Section Header', master = 'Office Theme')  %>% 
#       ph_with(value = ensembleTitle, location = ph_location_type(type = "body"))
#     
#     fileNameCts <- paste0(lofOfGraphicsFiles, "THI/THIextremeCt.", speciesName, "_", k, "_", yearSpan, ".png")
#     extImgObs <- external_img(src = fileNameCts, width = 5, height = 8)
#     
#     add_slide(my_pres, layout = 'Title Only', master = 'Office Theme') %>% 
#       ph_with(value = extImgObs, location = ph_location(left = 2, top = 0, width = 5, height = 8) )
#     
#     fileNameObserved <- paste0(lofOfGraphicsFiles, "THI/THIextremeCt.", speciesName, "_historical_",  "1991_2010", ".png")
#     
#     extImgObs <- external_img(src = fileNameObserved, width = 5, height = 8)
#     add_slide(my_pres, layout = 'Title Only', master = 'Office Theme') %>% 
#       ph_with(value = extImgObs, location = ph_location(left = 0, top = 0, width = 5, height = 8) )
#     
#     for (l in startYearChoices_ensemble) {
#       yearSpan <- paste0(l, "_", l + yearRange)
# #      fileNameCV <- paste0(lofOfGraphicsFiles, "THI/THI_ensembleCV_masked_",   speciesName, "_",  yearSpan, "_", k, ".jpg")
#       fileNameMean <- paste0(lofOfGraphicsFiles, "THI/THI_ensembleMean_masked_",  speciesName, "_",  yearSpan, "_", k, ".png")
#       
#       extImgMean <- external_img(src = fileNameMean, width = 5, height = 8)
# #      extImgCV <- external_img(src = fileNameCV, width = 5, height = 8)
#       
#       #   add_slide(my_pres, layout = 'Comparison', master = 'Office Theme') %>% 
#       #     ph_with(value = extImgMean, location = ph_location_left(),  use_loc_size = FALSE ) %>%
#       # #  add_slide(my_pres, layout = 'Comparison', master = 'Office Theme') %>% 
#       #     ph_with(value = extImgCV, location = ph_location_right(),  use_loc_size = FALSE )
#       
#       
#       add_slide(my_pres, layout = 'Title Only', master = 'Office Theme') %>% 
#         ph_with(value = extImgMean, location = ph_location(left = 0, top = 0, width = 5, height = 8) ) %>%
#         ph_with(value = extImgCV, location = ph_location(left = 5, top = 0, width = 5, height = 8) )
#       
#       my_pres <- add_slide(my_pres, layout = "Title and Content", master = "Office Theme")
#       my_pres <-  ph_with(x = my_pres, value = "Data Source", location = ph_location_type(type = "title"))
#       my_pres <- ph_with(x = my_pres, value = blData, location = ph_location_type(type = "body") )
#       
#     #     }
#   }
# }

# library(rgl)
# plot3d(x = r_end_ssp126_df$x, y = r_end_ssp126_df$y, z = r_end_ssp126_df$value, cols <- colorList,
#        xlab="longitude", ylab="latitude", zlab="Count of extreme stress")
# 
# library(plotly)
# plot_ly(x = r_end_ssp126_df$x, y = r_end_ssp126_df$y, z = r_end_ssp126_df$value, colors = colorList, type = "scatter3d", mode = "markers",
#         xlab="longitude", ylab="latitude", zlab="Count of extreme stress")

# runs, scenarios -----
thiMax = 88
test_logic <- paste0("x > ", thiMax)
k <- "ssp585"
l = "2041"

f_runs <- function(x, runlength, test_logic) {
  #  print(paste0("test_logic: ", test_logic))
  runResult <- c(NA, NA) 
  if (is.nan(x[1])) {
    return(runResult)
  }
  seqLengthCode <- paste0("1{", runlength, ",}") #A regular expression  to get the first item of gregexpr. It says look for  run_length times See http://xenon.stanford.edu/~xusch/regexp/
  g <- gregexpr(seqLengthCode, paste(+eval(parse(text = test_logic)), collapse = ""))[[1]] # The + converts TRUE and FALSE to 1 and 0
  #  print(paste0("g1: ", g[1]))
  if ((g[1] == -1)) { # no need to write to growing season if g returns -1, return 0,0
    runResult <- c(0, 0) 
    #    print("no runs")
  } else {
    startDays <- unlist(g)
    runLengths <- sum(as.numeric(attributes(g)$match.length))
    runResult <- c(length(startDays), runLengths)
  }
  return(runResult)
}

f_runs_calculator <- function(k, l, runlengthChoices) {
  yearSpan <- paste0(l, "_", l + yearRange)
  for (runlength in runlengthChoices) {
    for (m in breedChoice) {
      for (s in "extremeStress") {
        stressValue <- f_getStressValue(m, s) -1
        logicDirection <- ">"
        if (speciesChoice == "humans") logicDirection <- "<" # pwc less than thiStressVal is bad
        
        if (logicDirection == ">") ldtext <-"gt"
        if (logicDirection == "<") ldtext <-"lt"
        test_logic <- paste0("x ", logicDirection, " ", stressValue)
        #       print(paste0("test_logicFirst: ", test_logic))
        fileName_in <- paste0("data/cmip6/THI/ensemble_thi.", speciesChoice, "_", k, "_", yearSpan, ".tif")
        r <- rast(fileName_in)
        fileName_out <- paste0("data/cmip6/THI/run_", runlength, "_lim_", logicDirection, stressValue, "_ensemble_thi.", speciesChoice, "_", k, "_", yearSpan, ".tif")
        print(system.time(r_runs <- app(r, f_runs, runlength, test_logic, filename = fileName_out,  overwrite = TRUE, wopt= woptList)))
        #   print(system.time(writeRaster(r_runs, filename = fileName_out,  overwrite = TRUE, wopt= woptList))); flush.console()
        mainrl = paste0(m, ", THI gt ", stressValue, ", \nlongest no. of days in a min. run of ", runlength, " days, ", k, ", ", gsub("_", "-", yearSpan))
        if (m =="humans")  mainrl = paste0(m, ", PWC ", logicDirection, " , ", stressValue, ", \nlongest no of days in a min. run of ", runlength, " days, ", k, ", ", gsub("_", "-", yearSpan))
        mainct = paste0(m, ", THI gt ", stressValue, ", \nrun minimum length is ", runlength, " days, ", k, ", ", yearSpan)
        if (m =="humans")  mainct = paste0(m, ", PWC  ", logicDirection, " , ", stressValue, ", \nrun minimum length is ", runlength, " days, ", k, ", ", gsub("_", "-", yearSpan))
        
        plot(r_runs$lyr.1, main = mainct)
        plot(r_runs$lyr.2, main = mainrl)
      }
    }
  }
}

#run_length <- 40
runlengthChoices <- c(5, 10, 20, 40)
for (k in sspChoices) {
  for (l in startYearChoices) {
    f_runs_calculator(k, l, runlengthChoices)
  }
}

#runs, historical -----
k = "historical"
l = 1991
f_runs_calculator(k, l, runlengthChoices)

# runs graphics -----
dt_stressCts<- as.data.table(read.csv(paste0(locOfDataFiles_THI, "stressCtsTable.csv")))

# runs graphics, scenarios -----
for (k in sspChoices) {
  for (l in startYearChoices) {
    #    m = "cattle"
    for (m in breedChoice) {
      thiMax <- f_getStressValue(m, stressLevel) -1
      logicDirection <- ">"
      if (m == "humans") logicDirection <- "<" # pwc less than thiMax is bad
      
      if (logicDirection == ">") ldtext <-"gt"
      if (logicDirection == "<") ldtext <-"lt"
      test_logic <- paste0("x ", logicDirection, " ", thiMax)
      
      yearSpan <- paste0(l, "_", l + yearRange)
      fileName_in <- paste0("data/cmip6/THI/run_", runlength, "_lim_", logicDirection, thiMax, "_ensemble_thi.", speciesChoice, "_", k, "_", yearSpan, ".tif")
      print(paste0("fileName in: ", fileName_in))
      r <- rast(fileName_in)
      names(r) <- c("runsCt", "runLength")
      f_thi_runs_graphing(m, r, col, dt_stressCts, run_length)
      print("-------------------------")
    }
  }
}

# runs graphics, historical -----
k = "historical"
l = 1991
#    m = "cattle"
for (m in breedChoice) {
  thiMax <- f_getStressValue(m, stressLevel) -1
  logicDirection <- ">"
  if (m == "humans") logicDirection <- "<" # pwc less than thiMax is bad
  
  if (logicDirection == ">") ldtext <-"gt"
  if (logicDirection == "<") ldtext <-"lt"
  test_logic <- paste0("x ", logicDirection, " ", thiMax)
  
  yearSpan <- paste0(l, "_", l + yearRange)
  fileName_in <- paste0("data/cmip6/THI/run_", run_length, "_lim_", logicDirection, thiMax, "_ensemble_thi.", speciesChoice, "_", k, "_", yearSpan, ".tif")
  print(paste0("fileName in: ", fileName_in))
  r <- rast(fileName_in)
  names(r) <- c("runsCt", "runLength")
  f_thi_runs_graphing(m, r, col, dt_stressCts, run_length)
  print("-------------------------")
}
# runs, ppt -----
library(officer)
source("R/pptFunctions.R")

f_runsCtForPptx <- function(m) {
  fileName_in_rl <- paste0(lofOfGraphicsFiles, "THI/runLength.", speciesChoice, "_", k, "_", yearSpan, ".png")
  fileName_in_rCt <- paste0(lofOfGraphicsFiles, "THI/runsCt.", speciesChoice, "_", k, "_", yearSpan, ".png")
  print(paste0("fileName in, rl: ", fileName_in_rl, ", fileName in, cts: ", fileName_in_rCt))
  extImg_rl <- external_img(src = fileName_in_rl, width = defaultWidth, height = defaultHeight)
  extImg_rCt <- external_img(src = fileName_in_rCt, width = defaultWidth, height = defaultHeight)
  my_pres <- add_slide(x = my_pres, layout = 'Title Only', master = 'Office Theme')
  my_pres <- ph_with(x = my_pres, value = extImg_rCt, location = ph_location(left = defaultLeft, top = defaultTop, width = defaultWidth, height = defaultHeight - 0.5) )
  my_pres <- add_slide(x = my_pres, layout = 'Title Only', master = 'Office Theme')
  my_pres <- ph_with(x = my_pres, value = extImg_rl, location = ph_location(left = defaultLeft, top = defaultTop, width = defaultWidth, height = defaultHeight - 0.5) )
  return(my_pres)
}
thiList <- c("thi.humans", "thi.cattle", "thi.sheep", "thi.goat", "thi.yak", "thi.broiler", "thi.layer", "thi.chicken", "thi.pigs")
thiListReduced <- thiList[!thiList %in% c("thi.yak", "thi.broiler", "thi.layer")]

titleString <- paste0("Number of ", run_length, " day runs and longest run of extreme stress days for humans and ", length(thiListReduced)-1 , " animal species.")
contentString <- paste0("Preliminary results for two scenarios and three time periods to 2100. Powerpoint produced on ", Sys.Date())

IntroText0 <- "The negative impacts of extreme stress are amplified if the stress comes for an extended period. "
IntroText0.5 <- "\nThe following figures present global graphics of locations where extreme stress is experienced and how many consequtive days this happens in a representative year."
IntroText1 <- "\nFor the animals graphics, extreme stress areas are where the species was raised in the early 21st century."
IntroText4 <- "\nThe daily stress values are averages for three 20 year periods (1991-2010, 2041-2060, and 2081-2100). The input data are described in the last slide."
IntroText5 <- "\nThis powerpoint presents work in progress and should not be circulated without permission."

blIntro <- block_list(
  fpar(
    ftext(IntroText0, fp_2),
    ftext(IntroText0.5, fp_2)),
  fpar(ftext(IntroText1, fp_2)),
  fpar(ftext(IntroText4, fp_2)),
  fpar(ftext(IntroText5, fp_2))
)

my_pres <- read_pptx()
my_pres <- add_slide(x = my_pres, layout = 'Title Slide', master = 'Office Theme')
my_pres <- ph_with(x = my_pres, value = titleString, location = ph_location_type(type = "ctrTitle"))
my_pres <- ph_with(x = my_pres, value = contentString, location = ph_location_type(type = "subTitle"))

my_pres <- add_slide(my_pres, layout = "Title and Content", master = "Office Theme")
my_pres <-  ph_with(x = my_pres, value = "Introduction", location = ph_location_type(type = "title"))
my_pres <- ph_with(x = my_pres, value = blIntro, location = ph_location_type(type = "body") )

# do all slides for each species in speciesChoice
for (speciesChoice in speciesChoices) {
  speciesChoice <- gsub("thi.", "", speciesChoice)
  ensembleTitle <- speciesChoice
  my_pres <- add_slide(x = my_pres, layout = 'Section Header', master = 'Office Theme')
  my_pres <- ph_with(x = my_pres, value = ensembleTitle, location = ph_location_type(type = "title"))
  
  # do historical first, then ssps and future periods -----
  k <- "historical"
  l <- 1991
  yearSpan <- paste0(l, "_", l + yearRange)
  my_pres <- f_runsCtForPptx(m)
  
  for (k in sspChoices) {
    for (l in startYearChoices) {
      yearSpan <- paste0(l, "_", l + yearRange)
      my_pres <- f_runsCtForPptx(m)
    }
  }
}
my_pres <- f_addDataSlide()
fileName_out <- paste("presentations/cmip6/THI/extremeStressRuns_Ensemble_rl_", run_length, ".pptx")
print(my_pres, target = fileName_out) %>% browseURL()

# calculate mean values over a year for each pixel -----

#means, scenarios -----
ext_north <- ext(-180, 180, 50, 90)

{
  r_global <- data.frame(row.names = seq(1, 366, 1))
  for (k in sspChoices) {
    for (l in startYearChoices) {
      yearSpan <- paste0(l, "_", l + yearRange)
      for (speciesChoice in speciesChoices) {
        fileName_in <- paste0(locOfDataFiles_THI, "ensemble_thi.", speciesChoice, "_", k, "_", yearSpan, ".tif")
        print(paste0("fileName in: ", fileName_in))
        r <- rast(fileName_in)
        r <- crop(r, ext_north)
        print(system.time({r_mean <- global(r, "mean", na.rm = TRUE)
        r_max <- global(r, "max", na.rm = TRUE)
        r_min <- global(r, "min", na.rm = TRUE)}))
        r_combined <- cbind(r_min, r_mean, r_max)
        names(r_combined) <- c(paste0(m, "_", k, "_", l, "_min"), paste0(m, "_", k, "_", l, "_mean"), paste0(m, "_", k, "_", l, "_max"))
        r_global <- cbind(r_global, r_combined)
        #     fileName_out <- paste0(locOfDataFiles_THI, "mean_ensemble_thi.", speciesChoice, "_", k, "_", yearSpan, ".tif")
      }
    }
  }
  
  k = "historical"
  l = 1991
  yearSpan <- paste0(l, "_", l + yearRange)
  for (speciesChoice in speciesChoices) {
    fileName_in <- paste0(locOfDataFiles_THI, "ensemble_thi.", speciesChoice, "_", k, "_", yearSpan, ".tif")
    r <- rast(fileName_in)
    r <- crop(r, ext_north)
    print(paste0("fileName in: ", fileName_in))
    print(system.time({r_mean <- global(r, "mean", na.rm = TRUE)
    r_max <- global(r, "max", na.rm = TRUE)
    r_min <- global(r, "min", na.rm = TRUE)}))
    r_combined <- cbind(r_min, r_mean, r_max)
    names(r_combined) <- c(paste0(m, "_", k, "_", l, "_min"), paste0(m, "_", k, "_", l, "_mean"), paste0(m, "_", k, "_", l, "_max"))
    r_global <- cbind(r_global, r_combined)
    #     fileName_out <- paste0(locOfDataFiles_THI, "mean_ensemble_thi.", speciesChoice, "_", k, "_", yearSpan, ".tif")
  }
  r_global <- r_global[, sort(names(r_global))]
  write.csv(r_global, paste0(locOfDataFiles_THI, "thi_global_nh.csv"))
  
}

# mean graphics, scenarios -----
for (k in sspChoices) {
  for (l in startYearChoices) {
    yearSpan <- paste0(l, "_", l + yearRange)
    for (speciesChoice in speciesChoices) {
      fileName_in <- paste0(locOfDataFiles_THI, "mean_ensemble_thi.", speciesChoice, "_", k, "_", yearSpan, ".tif")
      print(paste0("fileName in: ", fileName_in))
      f_thi_means_graphing(m, r, col) 
    }
  }
}




# # get animal raster masks
# animalsOnly <- speciesChoice[!speciesChoice %in% "humans"]
# for (m in animalsOnly) {
#   if (m %in% "cattle") cutoff = 2000
#   fileName_r_mask <- paste0("data/animals/raster_ct_", speciesChoice, ".tif")
#   r_mask <- rast(fileName_r_mask)
#   values(r_mask)[values(r_mask) < 2000] <- NA
#   r_delta_mask <- mask(r_delta, r_mask)
# }
# r_delta_mask <- mask(r_delta, pop_mask)
# pop_mask
# #plot data with masking
# 
# for (m in speciesChoice) {
#   fileName_century_end <- paste0(locOfDataFiles_THI, "extremeCt.", speciesChoice, "_", k, "_", yearSpan_end, ".tif")
#   fileName_century_begin <- paste0(locOfDataFiles_THI, "extremeCt.", speciesChoice, "_", "historical", "_", yearSpan_begin, ".tif")
#   r_raster <- 
#     r_end <- rast(fileName_century_end)
#   r_begin <- rast(fileName_century_begin)
#   r_delta <- r_end - r_begin
#   
#   r_delta_df <- as.data.frame(r_delta_mask, xy = TRUE)
#   names(r_delta_df)[names(r_delta_df) == "lyr.1"] <- "value"
#   
#   titleText <- paste0(m, ", change in days above extreme stress value, \n\n early to end century ", "SSP585 ")
#   
#   legendTitle <- "Change in \n\nHigh Stress Days"
#   print(titleText)
#   gc()
#   colorList <- (RColorBrewer::brewer.pal(5, "YlOrRd"))
#   maxVal <- 250
#   custom_bins <- round(seq.int(from = 0, to = maxVal, length = 6))
#   r <- r_delta_df
#   r$value[r$value > maxVal] <- maxVal
#   g <- ggplot(data = coastline) +
#     labs(title = titleText, fill = legendTitle) + theme(plot.title = element_text(size = 12, hjust = 0.5)) +
#     labs(x = "", y = "") +
#     
#     geom_raster(data = r, aes(x, y, fill = value)) +
#     scale_fill_gradientn(colours = colorList, na.value = "white",
#                          breaks = custom_bins,labels = custom_bins,
#                          limits = c(0, maxVal)) + 
#     geom_sf(fill = NA, color = "gray")
#   g
#   fileName_out <- paste0(lofOfGraphicsFiles, "THI/THIextremeCtDelta.", speciesChoice, "_", k, ".png")
#   png(filename = fileName_out, width = 6, height = 6, units = "in", res = 300)
#   print(g)
#   dev.off()
# }
# 
# # convert population density raster
# pop_mask <- rast("data-raw/gpw_v4_population_density_rev11_2020_30_min.tif")
# setMinMax(pop_mask)
# fileName_r_mask <- paste0("data/animals/raster_ct_", "humans", ".tif")
# writeRaster(pop_mask, fileName_r_mask, overwrite = TRUE, wopt = woptList)
# 
# maskValUpper <- 3000
# maskValLower <- 10
# pop_mask > maskValUpper
# pop_mask[pop_mask > maskValUpper] <- NA
# pop_mask[pop_mask < maskValLower] <- NA
# pop_mask[pop_mask > 0] <- 1
# 
# fileName_in <- locOfDataFiles_THI, "ensemble_thi.humans_ssp585_2081_2100.tif"
# test <- rast(fileName_in)
# test
# test_brick <- raster::brick(fileName_in)
# test_brick
# test_anim <- raster::animate(test_brick, 1)
# 
# # calculate number of sequential days with stress; minimum of 5
# 
# fun <- function(cellVector) {
#   startend <- c(NA, NA) 
#   # if (is.nan(cellVector[1])) {
#   #   return(startend)
#   # }
#   g <- gregexpr("1{5,}", paste(+(cellVector > 89), collapse = ""))[[1]]
#   print(g)
#   if (!g[[1]] == -1) { # no need to write to growing season if g returns 1
#     startend[1] <- g[[1]]
#     matchLength <- attributes(g)$match.length
#     startend[2] <- startend[1] + matchLength - 1
#   }
#   return(startend) 
# }
# 
# fileName_in <- paste0(locOfDataFiles_THI, "ensemble_thi.", speciesChoice, "_", k, "_", yearSpan, ".tif")
# r <- rast(fileName_in)
# 
# print(system.time(extremeStress <- app(r, fun)))
# I discoverd that the functions below were not being used, as of Jan 11, 2021. I put them here temporarily 
# f_readRast_thi_cattle <- function(yearNumber) {
#   fileName_in <- paste0(locOfDataFiles_THI, "thi.cattle_", i, "_", k,  "_", yearNumber, ".tif") # note yearNumber here
#   r <- rast(fileName_in)
# }
# 
# f_readRast_thi_sheep <- function(yearNumber) {
#   fileName_in <- paste0(locOfDataFiles_THI, "thi.sheep_", i, "_", k,  "_", yearNumber, ".tif") # note yearNumber here
#   print(fileName_in)
#   r <- rast(fileName_in)
# }
# f_readRast_thi_goat <- function(yearNumber) {
#   fileName_in <- paste0(locOfDataFiles_THI, "thi.goat_", i, "_", k,  "_", yearNumber, ".tif") # note yearNumber here
#   print(fileName_in)
#   r <- rast(fileName_in)
# }
# f_readRast_thi_chicken <- function(yearNumber) {
#   fileName_in <- paste0(locOfDataFiles_THI, "thi.chicken_", i, "_", k,  "_", yearNumber, ".tif") # note yearNumber here
#   print(fileName_in)
#   r <- rast(fileName_in)
# }
# 
# f_readRast_thi_pigs <- function(yearNumber) {
#   fileName_in <- paste0(locOfDataFiles_THI, "thi.pigs_", i, "_", k,  "_", yearNumber, ".tif") # note yearNumber here
#   print(fileName_in)
#   r <- rast(fileName_in)
# }

PWCadj <- library(readxl)
PWCadj <- read_excel("/Volumes/ExtremeSSD3/data/animals/PWCadj.xlsx")


