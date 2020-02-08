
#--- CITATION ------------------------------------------------------------------

#*******************************************************************************
# Hawkins, BL. 2019. Exploring the role of temperature and possible alternative 
#   stable states in brook trout (Salvelinus fontinalis) population structure.
#   M.S. Thesis. Univeristy of California San Diego.
#*******************************************************************************

#--- SETUP ---------------------------------------------------------------------

# clear environment
rm(list = ls(all = TRUE))

# libraries
library(lubridate)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(quantmod)
library(runjags)

# create directory for today's output
today <- substr(Sys.time(), start = 1, stop = 16)
today <- gsub(pattern = " ", replacement = "-", x = today)
today <- gsub(pattern = ":", replacement = "-", x = today)
results.directory <- paste0("Results/", today)
if(!dir.exists(results.directory)) {
    dir.create(results.directory)
}
rm(today)


#--- PROCESS DATASETS ----------------------------------------------------------

# process lake and fish data for both datasets into the same format

# input:
# 1) lake.data - contains information about each lake (data.frame)
# 2) length.data - contains information about each fish (data.frame)
# output:
# process.fish.data - contains information about lakes and fish, where each row
#                     is a fish (data.frame)
# notes:
# highly specific to Rob Grasso dataset!
processRobGrassoData <- function(lake.data, length.data) {

    length.data <- length.data %>%
       mutate(Site = as.character(Site.Name),
              Species = as.character(FishSpecies),
              FishLength = as.numeric(as.character(FishLength.CM.)) * 10) %>%
       select(Site, Year, Species, FishLength, Comments)
    
    lake.data <- lake.data %>%
        mutate(Site = as.character(CommonName),
               Elevation = Elevation / 3.28084) %>%
        select(Site, Elevation)

    the.data <- left_join(length.data, lake.data, by = "Site")
    
    the.data <- the.data %>%
        filter(!is.na(FishLength)) %>%
        filter(Species != "") %>%
        filter(Comments != "may be recorded in mm?")
    
    the.data$Species <- gsub(pattern = " ", replacement = "", x = the.data$Species)
    
    the.data <- the.data %>%
        mutate(LakeYear = paste0(Site, "-", Year)) %>%
        select(LakeYear, Site, Year, Elevation, Species, FishLength) %>%
        mutate(Site = as.factor(Site),
               LakeYear = as.factor(LakeYear),
               Species = as.factor(Species)) %>%
        droplevels()
    
    return(the.data)
}

# input:
# 1) lake.data - contains information about each lake (data.frame)
# 2) length.data - contains information about each fish (data.frame)
# output:
# process.fish.data - contains information about lakes and fish, where each row
#                     is a fish (data.frame)
# notes:
# highly specific to Roland Knapp dataset!
processRolandKnappData <- function(lake.data, length.data) {
    
    length.data <- length.data %>%
        select(Lake.ID, Date, Species, Length) %>%
        rename(Site = Lake.ID, FishLength = Length) %>%
        mutate(Date = ymd(Date), tz = "America/Los_Angeles") %>%
        mutate(Year = as.integer(year(Date))) %>%
        mutate(LakeYear = as.factor(paste0(Site, "-", Year))) %>%
        select(LakeYear, Site, Year, Species, FishLength)
    
    lake.data <- lake.data %>%
        rename(Site = Lake.ID)
    
    the.data <- left_join(length.data, lake.data, by = "Site")
    
    # remove 15 observations with NA for fish length
    the.data <- the.data %>%
        filter(!is.na(FishLength))
    
    return(the.data)
}

# load and process Rob Grasso data
lake.data <- read.csv(file = "Data/RobGrassoYosemite/meta.data.csv")
length.data <- read.csv(file = "Data/RobGrassoYosemite/length.data.csv")
RG.process.fish.data <- processRobGrassoData(lake.data = lake.data, length.data = length.data)
rm(lake.data, length.data)

# load and process Roland Knapp data
lake.data <- read.csv(file = "Data/RolandKnappYosemite/Lake ID and elevation.csv", header = TRUE)
length.data <- read.csv(file = "Data/RolandKnappYosemite/Symons_YOSE_FishData_2016.csv", header = TRUE)
RK.process.fish.data <- processRolandKnappData(lake.data = lake.data, length.data = length.data)
rm(lake.data, length.data)


#--- FILTER DATASETS -----------------------------------------------------------

# select trout populations of interest

# input:
# 1) process.fish.data - from above section (data.frame)
# 2) species - a factor from the Species column of process.fish.data (character)
# 3) min.cutoff - minimum observations of species measured in a population (integer)
# 4) max.cutoff - maximum observations of species measured in a populations (integer)
# 5) first.year - would you like to include only the first year that a lake-year
#                 was observed with at least cutoff observations? (boolean)
# output:
# filter.fish.data - filtered version of process.fish.data, with new column of
#                    sample size for species of choice (data.frame)
filterData <- function(process.fish.data, species, min.cutoff = 0, max.cutoff = 1e5, first.year = FALSE) {
    
    # calculate sample size for species for each lake-year
    sample.size.data <- process.fish.data %>%
        filter(Species == species) %>%
        group_by(LakeYear) %>%
        summarize(SampleSize = n()) %>%
        arrange(desc(SampleSize))
    
    # filter data based on sample size
    filter.fish.data <- process.fish.data %>%
        filter(Species == species) %>%
        left_join(sample.size.data, by = "LakeYear") %>%
        filter(SampleSize >= min.cutoff & SampleSize <= max.cutoff) %>%
        select(LakeYear, Site, Year, Elevation, Species, SampleSize, FishLength) %>%
        droplevels()
    
    # filter data based on first
    if (first.year) {
        filter.fish.data <- filter.fish.data %>%
            group_by(Site) %>%
            arrange(Year) %>%
            filter(Year == first(Year)) %>%
            ungroup()  %>%
            droplevels()
    }
    
    # print message if filtered data is empty
    if (nrow(filter.fish.data) == 0) print("Filtered dataframe is empty!")
    
    return(filter.fish.data)
}

# filter data
RG.filter.fish.data <- filterData(process.fish.data = RG.process.fish.data,
                                  species = "BT",
                                  min.cutoff = 25,
                                  first.year = TRUE)
RK.filter.fish.data <- filterData(process.fish.data = RK.process.fish.data,
                                  species = "BK",
                                  min.cutoff = 25,
                                  first.year = FALSE)

# remove process.fish.data
rm(RG.process.fish.data, RK.process.fish.data)


#--- STEP 1: CALCULATE AVERAGE FISH SIZE ---------------------------------------

# calculate average fish size (mm) for each trout population

# input:
# 1) filter.fish.data - from above section (data.frame)
# output:
# average.fish.data - filter.fish.data with new column for mean and median 
#                     fish size for each population
calculateAverageFishSize <- function(filter.fish.data) {
    
    average.size.data <- filter.fish.data %>%
        group_by(LakeYear) %>%
        summarize(MeanSize = mean(FishLength),
                  MedianSize = median(FishLength))
    
    average.fish.data <- filter.fish.data %>%
        left_join(average.size.data, by = "LakeYear") %>%
        select(LakeYear, Site, Year, Elevation, Species, MeanSize, MedianSize,
               SampleSize, FishLength)
    
    return(average.fish.data)
}

# calculate average fish size
RG.average.fish.data <- calculateAverageFishSize(filter.fish.data = RG.filter.fish.data)
RK.average.fish.data <- calculateAverageFishSize(filter.fish.data = RK.filter.fish.data)

# remove filter.fish.data
rm(RG.filter.fish.data, RK.filter.fish.data)


#--- STEP 2: FIND PEAKS --------------------------------------------------------

# find the number and location of size classes for each population,
# assuming that each discrete size class will show up as a 'peak' in a density
# plot of a population's fish lengths

# function to return x and y output from density function, called in returnPeaks
returnDensityXY <- function(data) {
    density.output <- density(data$FishLength)
    density.xy <- data.frame(X = density.output$x, Y = density.output$y)
    return(density.xy)
}

# function to return peak location, called in calculateNumberPeaks
returnPeaks <- function(data) {
    peak.indices <- findPeaks(data$Y) - 1
    peak.locations <- data$X[peak.indices]
    return(peak.locations)
}

# input:
# 1) average.fish.data - from above section (tbl_df)
# output:
# output - contains two data.frames (list)
# output[[1]] - lakeyear.data, contains number of size classes, where each row
#               is a population (data.frame)
# output[[2]] - algebraic.peak.data, contains location of size classes, where
#               each row is a size class(data.frame)
# notes:
# algebraic.peak.data was not used for figures in my thesis, but was created 
# to inform the priors for the Bayesian normal mixture models used in step 3; it
# was also created while developing the methods in order to compare output
# between the peak-finding (step 2) and model-fitting (step 3) approaches.
calculateNumberPeaks <- function(average.fish.data) {
    
    # split average.fish.data by lake year
    split.by.lakeyear <- average.fish.data %>%
        select(LakeYear, FishLength) %>%
        split(average.fish.data$LakeYear)
    
    # find density of each lake year
    density.output <- lapply(split.by.lakeyear, returnDensityXY)
    
    # find peaks of each lake year
    peak.output <- lapply(density.output, returnPeaks)
    
    # reformat peak output into dataframe for lakeyear data
    number.peaks <- lapply(peak.output, length)
    LakeYear <- names(number.peaks)
    NumberPeaks <- unlist(number.peaks)
    lakeyear.data <- data.frame(LakeYear, NumberPeaks)
    rownames(lakeyear.data) <- NULL
    rm(LakeYear, NumberPeaks)
    
    # reformat peak output into dataframe for peak data
    LakeYear <- mapply(rep, x = names(peak.output), times = number.peaks) %>% unlist()
    Peak <- mapply(seq, from = 1, to = number.peaks) %>% unlist()
    PeakLocation <- unlist(peak.output)
    peak.data <- data.frame(LakeYear, Peak, PeakLocation)
    colnames(peak.data)[1] <- "LakeYear"
    rownames(peak.data) <- NULL
    rm(LakeYear, Peak, PeakLocation)
    
    # add additional information from input data to output lakeyear data
    add.to.lakeyear <- average.fish.data %>%
        group_by(LakeYear) %>%
        summarize(Site = first(Site),
                  Year = first(Year),
                  Elevation = first(Elevation),
                  Species = first(Species),
                  MeanSize = first(MeanSize),
                  MedianSize = first(MedianSize),
                  SampleSize = first(SampleSize))
    lakeyear.data <- lakeyear.data %>%
        left_join(add.to.lakeyear, by = "LakeYear") %>%
        select(LakeYear, Site, Year, Elevation, Species, MeanSize, MedianSize,
               SampleSize, NumberPeaks)
    
    return(list(lakeyear.data, peak.data))
}

# find peaks
RG.output <- calculateNumberPeaks(RG.average.fish.data)
RG.lakeyear.data <- RG.output[[1]]
RG.algebraic.peak.data <- RG.output[[2]]
rm(RG.output)

# find peaks
RK.output <- calculateNumberPeaks(RK.average.fish.data)
RK.lakeyear.data <- RK.output[[1]]
RK.algebraic.peak.data <- RK.output[[2]]
rm(RK.output)

# remove algebraic peak data, only used bayesian peak data for analysis
# (see note above calculationNumberPeaks)
rm(RG.algebraic.peak.data, RK.algebraic.peak.data)


#--- STEP 3: FIT BAYESIAN NORMAL MIXTURE MODELS --------------------------------

# find the number, location, variation, and proportion of fish in each size
# class, assuming that each discrete size class can be represented by a normal 
# component in a mixture model

# source models
source(file = "Code/hawkins_thesis_models.R")

# inputs:
# 1) lakeyear.data - from section above (data.frame)
# 2) fish.data - created from average.fish.data from section above (data.frame)
# 3) index - refers to row in lakeyear.data (integer)
# output:
# output - contains integer and data.frame (list)
# output[[1]] - bayesian.model.data, contains number of size classes and 
#               deviance, where each row is one model fit to a population
#               (data.frame)
# output[[2]] - bayesian.peak.data, contains information about each size class,
#               where each row is one normal component from the specified model
#               fit to the specified population (data.frame)
# notes:
# this function is run on each population later on in a for loop (hence the
# 'index' argument).
# this function fits two normal mixture models to each population, and it takes
# some time to do so; this is the most time-consuming portion of this script.
fitModels <- function(lakeyear.data, fish.data, index) {
    
    # find this lake-year
    this.lakeyear <- lakeyear.data$LakeYear[index]
    
    # subset fish data for this lake-year
    this.fish.data <- fish.data %>%
        subset(LakeYear == this.lakeyear) %>%
        select(FishLength) %>%
        arrange(FishLength) %>%
        as.matrix()
    
    # put fish data into JAGS readble format
    dat <- dump.format(list(d = this.fish.data, nf = nrow(this.fish.data)))
    
    # assign model for this lake-year
    algebraic.peaks <- lakeyear.data$NumberPeaks[index]
    
    # decide which models to fit
    if (algebraic.peaks > 1) {
        peaks.to.fit <- (algebraic.peaks-1):algebraic.peaks
    } else if (algebraic.peaks == 1) {
        peaks.to.fit <- algebraic.peaks:(algebraic.peaks+1)
    } else {
        print("Something is weird with the number of algebraic peaks!")
        return(NA)
    }

    # fit models
    for (this.peaks in peaks.to.fit) {
        if (this.peaks == 1) {
            model <- one.peak.model
        } else if (this.peaks == 2) {
            model <- two.peak.model
        } else if (this.peaks == 3) {
            model <- three.peak.model
        } else if (this.peaks == 4) {
            model <- four.peak.model
        } else if (this.peaks == 5) {
            model <- five.peak.model
        } else if (this.peaks == 6) {
            model <- six.peak.model
        } else {
            print("There's no model for the number of peaks to fit!")
            return(NA)
        }
        
        # run model for this lake-year and this number of peaks
        model.output <- run.jags(model,
                                 monitor = c("sdf", "m", "omega"),
                                 data = dat, burnin = 5000, sample = 10000,
                                 n.chains = 3)
        
        # summarize parameters and model performance for this lake-year and this number of peaks
        model.summary <- extract.runjags(model.output, "summary")
        model.dic <- extract.runjags(model.output, "dic")
        model.median <- model.summary[ , 2]
        
        # save peak and model characteristics to dataframes
        LakeYear <- this.lakeyear
        NumberPeaks <- this.peaks
        Deviance <- sum(model.dic$deviance)
        Penalty <- sum(model.dic$penalty)
        PenalizedDeviance <- Deviance + Penalty
        this.bayesian.model.data <- data.frame(LakeYear, NumberPeaks, Deviance, Penalty, PenalizedDeviance)
        row.names(this.bayesian.model.data) <- NULL
        
        SD <- model.median[1:this.peaks]
        Mean <- model.median[(this.peaks+1):(this.peaks*2)]
        ProportionFish <- model.median[(this.peaks*2+1):(this.peaks*3)]
        this.bayesian.peak.data <- data.frame(LakeYear, NumberPeaks, Mean, SD, ProportionFish) %>%
            arrange(Mean) %>%
            cbind(Peak = 1:this.peaks) %>%
            select(LakeYear, NumberPeaks, Peak, Mean, SD, ProportionFish)
        
        if (this.peaks == peaks.to.fit[1]) {
            # initialize dataframe for first model
            bayesian.model.data <- this.bayesian.model.data
            bayesian.peak.data <- this.bayesian.peak.data
        } else {
            # add rows to dataframe for following models
            bayesian.model.data <- rbind(bayesian.model.data, this.bayesian.model.data)
            bayesian.peak.data <- rbind(bayesian.peak.data, this.bayesian.peak.data)
        }
    }
    
    # return list of model fit metrics and peak charateristics
    return(list(bayesian.model.data, bayesian.peak.data))
}

# prep data for function
RG.fish.data <- RG.average.fish.data %>%
    select(LakeYear, FishLength)
RK.fish.data <- RK.average.fish.data %>%
    select(LakeYear, FishLength)

# remove average.fish.data
rm(RG.average.fish.data, RK.average.fish.data)

# fit models to each population
for (this.index in 1:nrow(RG.lakeyear.data)) {
    
    # fit models
    output <- fitModels(lakeyear.data = RG.lakeyear.data,
                        fish.data = RG.fish.data,
                        this.index)
    
    # save output
    if(this.index == 1) {
        # initialize dataframe for first index
        RG.bayesian.model.data <- output[[1]]
        RG.bayesian.peak.data <- output[[2]]
    } else {
        # add rows to dataframe for following indices
        RG.bayesian.model.data <- rbind(RG.bayesian.model.data, output[[1]])
        RG.bayesian.peak.data <- rbind(RG.bayesian.peak.data, output[[2]])
    }
}

# fit models to each population
for (this.index in 1:nrow(RK.lakeyear.data)) {
    
    # fit models
    output <- fitModels(lakeyear.data = RK.lakeyear.data,
                        fish.data = RK.fish.data,
                        this.index)
    
    # save output
    if(this.index == 1) {
        # initialize dataframe for first index
        RK.bayesian.model.data <- output[[1]]
        RK.bayesian.peak.data <- output[[2]]
    } else {
        # add rows to dataframe for following indices
        RK.bayesian.model.data <- rbind(RK.bayesian.model.data, output[[1]])
        RK.bayesian.peak.data <- rbind(RK.bayesian.peak.data, output[[2]])
    }
}

# remove temporary variables from for loop
rm(this.index, output)

# remove fish.data
rm(RG.fish.data, RK.fish.data)


#--- CALCULATE ADDITIONAL SIZE STRUCTURE FEATURES ------------------------------

# calculate separation between size classes and the evenness of distribution of
# size classes

# input:
# 1) mean - column from bayesian.peak.data containing mean of each size class 
#           (vector)
# output:
# output - separation between neighboring peaks (vector)
calculatePeakSeparation <- function(mean) {
    return(c(NA, diff(mean, lag = 1)))
}

# input:
# 1) proportion - column from bayesian.peak.data containing proportion of fish
#                 in each size class (vector)
# output:
# output - Pielou's evenness index for populations with multiple size classes or
#          NA for populations with one size class (integer)
calculatePeakEvenness <- function(proportion) {
    if (length(proportion) == 1) {
        return(NA)
    } else {
        H <- -sum(proportion * log(proportion))
        S <- length(proportion)
        Hmax <- log(S)
        J <- H / Hmax
        return(J)
    }
}

# add peak separation and evenness to bayesian peak data
RG.bayesian.peak.data <- RG.bayesian.peak.data %>%
    group_by(LakeYear, NumberPeaks) %>%
    mutate(Separation = calculatePeakSeparation(Mean),
           Evenness = calculatePeakEvenness(ProportionFish)) %>%
    select(LakeYear, NumberPeaks, Peak, Mean, Separation, SD, ProportionFish, Evenness) %>%
    ungroup()
RK.bayesian.peak.data <- RK.bayesian.peak.data %>%
    group_by(LakeYear, NumberPeaks) %>%
    mutate(Separation = calculatePeakSeparation(Mean),
           Evenness = calculatePeakEvenness(ProportionFish)) %>%
    select(LakeYear, NumberPeaks, Peak, Mean, Separation, SD, ProportionFish, Evenness) %>%
    ungroup()


#--- PROCESS ROB GRASSO SIZE CLASSES -------------------------------------------

# find better fitting Bayesian normal mixture model for each population,
# extract and save size class features from better fitting model

# read in lake data again to add additional environmental features
lake.data <- read.csv(file = "Data/RobGrassoYosemite/meta.data.csv")

# find which model has a lower penalized deviance for each lake year
# end product:
# better.model - contains the number of normal components from the better fitting
#                model for each population, each row is a popultaion (tbl_df)
better.model <- RG.bayesian.model.data %>%
    group_by(LakeYear) %>%
    filter(PenalizedDeviance == min(PenalizedDeviance)) %>%
    select(LakeYear, NumberPeaks) %>%
    ungroup()

# find peak characteristics for better fitting model for each lake year
# end product:
# better.peak - contains information about the normal components from the better
#               fitting model for each population, each row is a size class
#               (tbl_df)
for (this.lakeyear in better.model$LakeYear) {
    this.number.peaks <- better.model %>% 
        filter(LakeYear == this.lakeyear) %>%
        select(NumberPeaks) %>%
        as.numeric()
    this.peaks <- RG.bayesian.peak.data %>%
        filter(LakeYear == this.lakeyear) %>%
        filter(NumberPeaks == this.number.peaks)
    if (this.lakeyear == better.model$LakeYear[1]) {
        better.peak <- this.peaks
    } else {
        better.peak <- rbind(better.peak, this.peaks)
    }
}

# find environment data for each population, used to create population.data
environment.data <- lake.data %>%
    rename(Site = CommonName) %>%
    mutate(Depth = Depth / 3.28084,
           Elevation = Elevation / 3.28084) %>%
    select(-Location, -Descirption) %>%
    right_join(RG.lakeyear.data, by = "Site") %>% 
    rename(Elevation = Elevation.x) %>%
    mutate(Year = as.integer(Year),
           Site = as.factor(Site)) %>%
    select(LakeYear, Site, Year, Elevation, Depth, UTMe, UTMn) %>%
    as.tbl()

# find evenness data for each population, used to create population.data
evenness.data <- better.peak %>%
    group_by(LakeYear) %>%
    filter(Peak == first(Peak)) %>%
    ungroup() %>%
    select(LakeYear, Evenness)

# find average length data for each population, used to create population.data
meansize.data <- RG.lakeyear.data %>% 
    select(LakeYear, MeanSize)

# save size structure and environmental information together
# end product:
# population.data - contains size sturcture features and environmental features 
#                   for each population, each row is a populatoin (tbl_df)
population.data <- better.model %>%
    right_join(environment.data, by = "LakeYear") %>%
    left_join(evenness.data, by = "LakeYear") %>%
    left_join(meansize.data, by = "LakeYear") %>%
    select(LakeYear, Site, Year, Elevation, Depth, UTMe, UTMn, 
           MeanSize, NumberPeaks, Evenness)

# save model, size class, and population data for RG dataset
RG.model <- better.model
RG.peak <- better.peak
RG.population <- population.data

# remove temporary variables from for loops
rm(this.lakeyear, this.number.peaks, this.peaks)

# remove data structures no longer needed
rm(RG.bayesian.model.data, RG.bayesian.peak.data, lake.data, better.model, 
   better.peak, environment.data, population.data, evenness.data, meansize.data)


#--- PROCESS ROLAND KNAPP SIZE CLASSES -----------------------------------------

# find better fitting Bayesian normal mixture model for each population,
# extract and save size class features from better fitting model

# read in lake data again to add additional environmental features
lake.data <- read.csv(file = "Data/RolandKnappYosemite/LakeData_ToShurin_10Apr2019data.csv",
                      header = TRUE)

# find which model has a lower penalized deviance for each lake year
# end product:
# better.model - contains the number of normal components from the better fitting
#                model for each population, each row is a popultaion (tbl_df)
better.model <- RK.bayesian.model.data %>%
    group_by(LakeYear) %>%
    filter(PenalizedDeviance == min(PenalizedDeviance)) %>%
    select(LakeYear, NumberPeaks) %>%
    ungroup()

# find peak characteristics for better fitting model for each lake year
# end product:
# better.peak - contains information about the normal components from the better
#               fitting model for each population, each row is a size class
#               (tbl_df)
for (this.lakeyear in better.model$LakeYear) {
    this.number.peaks <- better.model %>% 
        filter(LakeYear == this.lakeyear) %>%
        select(NumberPeaks) %>%
        as.numeric()
    this.peaks <- RK.bayesian.peak.data %>%
        filter(LakeYear == this.lakeyear) %>%
        filter(NumberPeaks == this.number.peaks)
    if (this.lakeyear == better.model$LakeYear[1]) {
        better.peak <- this.peaks
    } else {
        better.peak <- rbind(better.peak, this.peaks)
    }
}

# find environment data for each population, used to create population.data
environment.data <- lake.data %>%
    transmute(Site = as.factor(lake_id),
              Area = area,
              Perimeter = perimeter,
              Elevation = elevation,
              UTMe = utm_east,
              UTMn = utm_north,
              Depth = max_depth,
              StreamSpawn = stream_spawnhab,
              LakeSpawn = lakestream_spawnhab) %>%
    as.tbl()

# find evenness data for each population, used to create population.data
evenness.data <- better.peak %>%
    group_by(LakeYear) %>%
    filter(Peak == first(Peak)) %>%
    ungroup() %>%
    select(LakeYear, Evenness)

# find average length data for each population, used to create population.data
meansize.data <- RK.lakeyear.data %>% 
    select(LakeYear, MeanSize)

# save size structure and environmental information together
# end product:
# population.data - contains size sturcture features and environmental features 
#                   for each population, each row is a populatoin (tbl_df)
population.data <- better.model %>%
    mutate(Site = substr(LakeYear, start = 1, stop = 5),
           Year = substr(LakeYear, start = 7, stop = 11)) %>%
    left_join(environment.data, by = "Site") %>%
    left_join(evenness.data, by = "LakeYear") %>%
    left_join(meansize.data, by = "LakeYear") %>%
    mutate(Site = as.factor(Site),
           Year = as.factor(Year)) %>%
    select(LakeYear, Site, Year, Elevation, Depth, UTMe, UTMn, 
           MeanSize, NumberPeaks, Evenness)

# save model, size class, and population data for RK dataset
RK.model <- better.model
RK.peak <- better.peak
RK.population <- population.data

# remove temporary variables from for loops
rm(this.lakeyear, this.number.peaks, this.peaks)

# remove data structures no longer needed
rm(RK.bayesian.model.data, RK.bayesian.peak.data, lake.data, better.model, 
   better.peak, environment.data, population.data, evenness.data, meansize.data)


#--- COMBINE ROB GRASSO & ROLAND KNAPP DATASETS --------------------------------

# combine separate datasets

# model datasets
RK.model <- RK.model %>% 
    mutate(LakeYear = as.character(LakeYear),
           Dataset = "RK")
RG.model <- RG.model %>% 
    mutate(LakeYear = as.character(LakeYear),
           Dataset = "RG")
best.model <- bind_rows(RK.model, RG.model) %>%
    mutate(LakeYear = as.factor(LakeYear),
           Dataset = as.factor(Dataset))

# remove model datasets
rm(RK.model, RG.model)
# best.model was not used for figures in my thesis, but was created while
# developing the methods in order to compare how well models fit across populations.

# peak datasets
RK.peak <- RK.peak %>% 
    mutate(LakeYear = as.character(LakeYear),
           Dataset = "RK")
RG.peak <- RG.peak %>% 
    mutate(LakeYear = as.character(LakeYear),
           Dataset = "RG")
best.peak <- bind_rows(RK.peak, RG.peak) %>%
    mutate(LakeYear = as.factor(LakeYear),
           Dataset = as.factor(Dataset))

# remove peak datasets
rm(RK.peak, RG.peak)
# best.peak was not used for figures in my thesis, but was created while developing
# the methods in order to explore the features of each size class individually
# in addition to the size structure of each population.

# population datasets
RK.population <- RK.population %>% 
    mutate(LakeYear = as.character(LakeYear),
           Site = as.character(Site),
           Year = as.integer(as.character(Year)),
           Dataset = "RK")
RG.population <- RG.population %>% 
    mutate(LakeYear = as.character(LakeYear),
           Site = as.character(Site),
           Year = as.integer(as.character(Year)),
           Dataset = "RG")
population.data <- bind_rows(RK.population, RG.population) %>%
    mutate(LakeYear = as.factor(LakeYear),
           Site = as.factor(Site),
           Dataset = as.factor(Dataset))
# population.data was used to create the figures in my thesis.


# remove separate population datasets
rm(RK.population, RG.population)


#--- STEP 4: COMBINE CLOSE SIZE CLASSES ----------------------------------------

# if separation between size classes is small, combine them

# define how small that separation is
separation.cutoff <- 50

# find lake-years that should change
collapse.lakeyears <- best.peak %>% 
    filter(Separation < separation.cutoff) %>% 
    select(LakeYear) %>% 
    pull()

# combine size classes
for (this.lakeyear in collapse.lakeyears) {
    
    # subset peak data for lake year
    this.best.peak <- best.peak %>% 
        filter(LakeYear == this.lakeyear)
    
    # find which peaks need to be collapsed
    which.to.collapse <- which(this.best.peak$Separation < separation.cutoff)
    which.to.collapse <- c(which.to.collapse-1, which.to.collapse)
    
    # subset peak data for lake year and peaks to collapse
    this.best.peak <- best.peak %>% 
        filter(LakeYear == this.lakeyear) %>% 
        filter(Peak %in% which.to.collapse)
    
    # create new row of collapsed peak for peak data
    new.peak.row <- data.frame(LakeYear = this.best.peak$LakeYear[1],
                               Mean = mean(c(this.best.peak$Mean[1], this.best.peak$Mean[2])),
                               Separation = 0,
                               SD = mean(c(this.best.peak$SD[1], this.best.peak$SD[2])),
                               ProportionFish = sum(this.best.peak$ProportionFish[1], this.best.peak$ProportionFish[2]),
                               Evenness = this.best.peak$Evenness[1],
                               Dataset = this.best.peak$Dataset[1])
    
    # create new rows for lake year including collapsed peak 
    new.peak.rows <- best.peak %>% 
        filter(LakeYear == this.lakeyear) %>% 
        filter(Peak != which.to.collapse[1]) %>% 
        filter(Peak != which.to.collapse[2]) %>% 
        select(-Peak, -NumberPeaks) %>% 
        rbind(new.peak.row) %>% 
        arrange(Mean) %>% 
        cbind(NumberPeaks = this.best.peak$NumberPeaks[1]-1,
              Peak = 1:(this.best.peak$NumberPeaks[1]-1)) %>%
        select(LakeYear, NumberPeaks, Peak, Mean, Separation, SD, ProportionFish,
               Evenness, Dataset) %>% 
        as.tbl()
    
    # recalculate separation and evenness
    calculatePeakSeparation <- function(mean) {
        return(c(NA, diff(mean, lag = 1)))
    }
    
    calculatePeakEvenness <- function(proportion) {
        H <- -sum(proportion * log(proportion))
        S <- length(proportion)
        Hmax <- log(S)
        J <- H / Hmax
        return(J)
    }
    
    new.peak.rows <- new.peak.rows %>% 
        mutate(Separation = calculatePeakSeparation(new.peak.rows$Mean),
               Evenness = calculatePeakEvenness(new.peak.rows$ProportionFish))
    
    # replace old peak data with new peak data
    best.peak <- best.peak %>% 
        filter(LakeYear != this.lakeyear) %>% 
        rbind(new.peak.rows)
    
    # create new model row for this lake year
    new.model.row <- data.frame(LakeYear = this.lakeyear,
                                NumberPeaks = new.peak.rows$NumberPeaks[1],
                                Dataset = new.peak.rows$Dataset[1])
    
    # replace old model data with new model data
    best.model <- best.model %>% 
        filter(LakeYear != this.lakeyear) %>% 
        rbind(new.model.row)
    
    # subset population data for lake year
    this.population.data <- population.data %>% 
        filter(LakeYear == this.lakeyear)
    
    # create new population row for this lake year
    new.population.row <- data.frame(LakeYear = this.lakeyear,
                                     Site = this.population.data$Site[1],
                                     Year = this.population.data$Year[1],
                                     Elevation = this.population.data$Elevation[1],
                                     Depth = this.population.data$Depth[1],
                                     UTMn = this.population.data$UTMn[1],
                                     UTMe = this.population.data$UTMe[1],
                                     MeanSize = this.population.data$MeanSize[1],
                                     NumberPeaks = new.peak.rows$NumberPeaks[1],
                                     Evenness = new.peak.rows$Evenness[1],
                                     Dataset = this.population.data$Dataset[1])
    
    # replace old population data with new population data
    population.data <- population.data %>% 
        filter(LakeYear != this.lakeyear) %>% 
        rbind(new.population.row)
}

# remove temporary variables from for loops
rm(this.best.peak, this.lakeyear, which.to.collapse, collapse.lakeyears, 
   new.peak.row, new.peak.rows, new.model.row, new.population.row, 
   this.population.data, separation.cutoff)


#--- RESULTS: Figure 2. Distribution of population structure features. ---------

a <- population.data %>% 
    ggplot() + 
    geom_histogram(aes(x = MeanSize)) + 
    xlab("Mean fish size (mm)") +
    theme(axis.title.y = element_blank())
b <- population.data %>% 
    ggplot() + 
    geom_histogram(aes(x = NumberPeaks)) + 
    xlab("Number of size classes") +
    theme(axis.title.y = element_blank())
c <- population.data %>% 
    ggplot() + 
    geom_histogram(aes(x = Evenness)) + 
    xlab("Evenness of size classes") +
    theme(axis.title.y = element_blank())

d <- grid.arrange(a, b, c, ncol = 3, left = textGrob("Count"))
ggsave(plot = d, filename = paste0(results.directory, "/Figure2.png"),
       width = 11, height = 4, units = "in")
rm(a, b, c, d)


#--- RESULTS: Figure 3. Population structure features plotted against elevation.----

a <- population.data %>% 
    ggplot() +
    geom_point(aes(x = Elevation, y = MeanSize)) +
    xlab("Elevation (m)") +
    ylab("Mean fish size (mm)")
b <- population.data %>% 
    ggplot() +
    geom_point(aes(x = Elevation, y = NumberPeaks)) +
    xlab("Elevation (m)") +
    ylab("Number of size classes")
c <- population.data %>% 
    ggplot() +
    geom_point(aes(x = Elevation, y = Evenness)) +
    xlab("Elevation (m)") +
    ylab("Evenness of size classes")

d <- grid.arrange(a, b, c, ncol = 3)

ggsave(plot = d, filename = paste0(results.directory, "/Figure3.png"),
       width = 11, height = 4, units = "in")
rm(a, b, c, d)


#--- RESULTS: Table 1. Twelve populations in fish-rich dataset.-----------------

table1 <- RG.lakeyear.data %>% 
    select(Site, Year, SampleSize)
write.csv(x = table1, 
          file = paste0(results.directory, "/Table1.csv"))
rm(RG.lakeyear.data, table1)


#--- RESULTS: Table 2. Thirty populations in lake-rich dataset.-----------------

table2 <- RK.lakeyear.data %>% 
    select(Site, Year, SampleSize)
write.csv(x = table2, 
          file = paste0(results.directory, "/Table2.csv"))
rm(RK.lakeyear.data, table2)


#--- RESULTS: Table 4. Spearman rank correlation. ------------------------------

# input:
# 1) population.data - from section above (data.frame)
# 2) response.name - population size structure feature (character)
# 2) predictor.name - environmental feature (character)
# output:
# output - population size structure feature, environmental feature, spearman's
#          rank correlation between the two features, p-value for correlation
#          (vector)
# notes:
# this function is run on each population later on in a for loop.
# although variable names imply variables are predictors and responses,
# that correlation just tests for the existence of a relationship. correlation
# does not indicate cause and effect, nor does it indicate that one variable
# predicts another variable's response.
populationCorrelation <- function(population.data, 
                                  response.name,
                                  predictor.name) {
    
    predictor.index <- which(colnames(population.data) == predictor.name)
    response.index <- which(colnames(population.data) == response.name)
    
    predictor <- pull(population.data[predictor.index])
    response <- pull(population.data[response.index])
    
    test <- cor.test(x = predictor, y = response, method = "spearman")
    
    output <- data.frame(
        Response = response.name,
        Predictor = predictor.name,
        Coefficient = test$estimate,
        PValue = test$p.value
    )
    
    row.names(output) <- NULL
    
    return(output)
}

# intialize output
correlation.output <- NULL
this.correlation.output <- NULL

# which environmental feature are you interested in?
predictor.name <- "Elevation"

# which size sturcture features are you interested in?
for (response.name in c("MeanSize", "NumberPeaks", "Evenness")) {

    # test correlation between environmental feature and size structure feature
    this.correlation.output <- populationCorrelation(
        population.data = population.data,
        response.name = response.name,
        predictor.name = predictor.name
    )
    
    # save output
    if(is.null(correlation.output)) {
        correlation.output <- this.correlation.output
    } else {
        correlation.output <- rbind(correlation.output, this.correlation.output)
    }
}

# remove function arguments
rm(predictor.name, response.name)

# remove temporary variables from for loop
rm(this.correlation.output)

# select and rename certain columns from correlation output
correlation.output <- correlation.output %>% 
    select(Response, Coefficient, PValue) %>% 
    rename(PopulationStructureFeature = Response,
           CorrelationCoefficient = Coefficient)

# save output in csv
write.csv(x = correlation.output, 
          file = paste0(results.directory, "/Table4.csv"),
          row.names = FALSE)

# remove table
rm(correlation.output)


#--- SAVE NEW DATASETS ---------------------------------------------------------

# save population data, in which each row is a population
population.data <- population.data %>% 
    select(Dataset, LakeYear, MeanSize, NumberPeaks, Evenness, Site, Year, Elevation, Depth, UTMn, UTMe) %>% 
    rename(NumberSizeClasses = NumberPeaks)
write.csv(x = population.data, 
          file = paste0(results.directory, "/population.csv"),
          row.names = FALSE)

# save size class data, in which each row is a size class
sizeclass.data <- best.peak %>% 
    select(Dataset, LakeYear, Peak, Mean, Separation, SD, ProportionFish) %>% 
    rename(SizeClass = Peak)
write.csv(x = sizeclass.data, 
          file = paste0(results.directory, "/sizeclass.csv"),
          row.names = FALSE)

# remove data structures no longer needed
rm(best.model, best.peak)

