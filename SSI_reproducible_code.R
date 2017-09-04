#import necessary libraries
library(data.table)
library(plyr)
library(dplyr)
library(gtools)
library(gridExtra)
library(ggplot2)

#set working directory
setwd("/Users/userName/folderName")

#import data
#my file consists of five columns: participant, condition, trial, alternative and attribute
infoSearch <- as.data.table(read.csv("fileName.csv", header = T, sep = ";")) 

#data preparation ####

#identify subsequent fixations to one attribute within one alternative
infoSearch$attributeClean <- ifelse(infoSearch$attribute == shift(infoSearch$attribute, 1L) 
                                    & infoSearch$alternative == shift(infoSearch$alternative, 1L), 1, 0) 

#delete subsequent fixations to one attribute i.e., keep only first fixation
infoSearch <- infoSearch[infoSearch$attributeClean != 1 | is.na(infoSearch$attributeClean)]

#delete unnecessary column
infoSearch[, "attributeClean" := NULL]

#create and count alternative-wise transitions (needed for calculating SI)
infoSearch$transAlt <- ifelse(infoSearch$attribute != lag(infoSearch$attribute, n = 1L) 
                              & infoSearch$alternative == lag(infoSearch$alternative, n = 1L) 
                              & infoSearch$trial == lag(infoSearch$trial, n = 1L) 
                              & infoSearch$participant == lag(infoSearch$participant, n = 1L), 1, 0)
infoSearch[is.na(infoSearch)] <- 0

#new object with alternative-wise transitions for each participant, condition and trial
altTrans <- ddply(infoSearch,.(participant, condition, trial), summarize, transAlt = sum(transAlt)) 

#create and count attribute-wise transitions (needed for calculating SI)
infoSearch$transAtt <- ifelse(infoSearch$attribute == lag(infoSearch$attribute, n = 1L) 
                              & infoSearch$alternative != lag(infoSearch$alternative, n = 1L) 
                              & infoSearch$trial == lag(infoSearch$trial, n = 1L) 
                              & infoSearch$participant == lag(infoSearch$participant, n = 1L), 1, 0)
infoSearch[is.na(infoSearch)] <- 0

#new object with attribute-wise transitions for each participant, condition and trial
attTrans <- ddply(infoSearch,.(participant, condition, trial), summarize, transAtt = sum(transAtt)) 

#combine two data sets by columns 
searchIndex <- as.data.table(cbind(attTrans, altTrans)) 

#calculate SI using apporpiate equation
searchIndex$searchIndex <- (searchIndex$transAlt - searchIndex$transAtt) / (searchIndex$transAlt + searchIndex$transAtt) 

#delete unnecessary columns
searchIndex[, c("participant", "condition", "trial") := NULL] 

#set order of columns
setcolorder(searchIndex, c("participant", "condition", "trial", "transAlt", "transAtt", "searchIndex"))  

#calculate total string length of eye fixations per participant per trial (needed for SSI denominator)
stringLength <- ddply(infoSearch, .(participant, condition, trial), function(infoSearch) length(infoSearch$attribute))
setnames(stringLength, "V1", "N") #change column name

#create counter variable for alternative-wise search (focusing on set of attributes when inspecting different alternatives)
#assigning same number to fixations on attributes within same alternatives 
#e.g., if there are fixations on sugar and fat level within one alternative and sugar and protein level within another alternative, 
#values 1,1,2,2 would be assigned to counter variable 
infoSearch <- setDT(infoSearch)[, counterAltwise:= rleid(condition, trial, alternative)] 

#create counter variable for attribute-wise search (focusing on one attribute when inspecting different alternatives)
#since subsequent fixations have been deleted, when there are fixations on one attribute, they must belong to different alternatives
infoSearch <- setDT(infoSearch)[, counterAttwise:= rleid(condition, trial, attribute)] 

#identify alternative-wise patterns ####

#create alternative-wise substrings (i.e., sequences of letters) based on counter variable
#collapsing all atributes within one alternative into string of letters 
#in example above, two substrings of length of two: 'fs' and 'ps' (i.e., 'fat and sugar' and 'protein and sugar')
altwiseStrings <- infoSearch[,list(string <- paste(attribute, collapse = ""), 
                                   participant = unique(participant), 
                                   condition = unique(condition), 
                                   trial = unique(trial)), by = counterAltwise] 

#delete counter variable 
altwiseStrings[, "counterAltwise" := NULL]

#rename column
setnames(altwiseStrings, "V1", "string")

#define function that keeps unique elements and sorts them alphabetically
relaxedFreqOrder <- function(i){
  paste0(unique(sort(unlist(strsplit(i, "")))), collapse = "")
}

#apply function to column with previosuly created alternative-wise substrings
altwiseStrings$formattedString <- lapply(altwiseStrings$string, relaxedFreqOrder)

#delete variable 
altwiseStrings[, "string" := NULL]

#rename column
setnames(altwiseStrings, "formattedString", "string")

#change class
altwiseStrings$string <- as.character(altwiseStrings$string)

#compare substrings to identify equal subsequent alternative-wise substrings  
altwiseStrings$equalStrings <- ifelse(altwiseStrings$string == lag(altwiseStrings$string, n = 1L) 
                                      | altwiseStrings$string == lead(altwiseStrings$string, n = 1L), 1, 0)

#create counter variable based on string variable within each trial 
#i.e., assign new number for every unique string within each trial
altwiseStrings <- setDT(altwiseStrings)[, counter:= rleid(string, trial)]

#extract equal subsequent substrings (equalStrings = 1)
altwiseStrings <- altwiseStrings[altwiseStrings$equalStrings != 0]

#delete substrings of length one
altwiseStrings <- subset(altwiseStrings, nchar(as.character(string)) >= 2)

#combine substrings into alternative-wise patterns using counter variable 
#i.e., all substrings with equal count should be collapsed into one pattern
altwisePatterns <- altwiseStrings[,list(string <- paste(string, collapse = ""), 
                                        participant = unique(participant), 
                                        condition = unique(condition), 
                                        trial = unique(trial)), by = counter]

#delete counter variable  
altwisePatterns[, "counter" := NULL]

#rename column
setnames(altwisePatterns, "V1", "pattern")

#calculate frequencies for each pattern within each trial, condition and for every participant 
altwisePatternsCount <- as.data.table(with(altwisePatterns, table(pattern, trial, condition, participant)))
altwisePatternsCount <- altwisePatternsCount[altwisePatternsCount$N != 0] 
setnames(altwisePatternsCount, "N", "pattFreq") #rename column

#assess whether obtained patterns occurred by chance by making random data set 
#to which we will compare the patterns from the original data set ####

altwiseSim <- function() { #function which contains random version of data
  
  #import data (same file as in line 14)
  infoSearch <- as.data.table(read.csv("fileName.csv", header = T, sep = ";"))

  #delete unnecessary columns
  infoSearch[, c("alternative", "attribute") := NULL]
  
  #create random data by sampling evenly from alternatives and attributes
  sim <- 154355 #number of rows corresponding to number of fixations in original data set (total string length)
  infoSearch$alternative <- sample(1:4, sim, T) #sample numbers from 1 to 4 154355 times
  infoSearch$attribute <- sample(c("b", "f", "p", "s"), sim, T) #sample letters b, f, p and s 154355 times

  #identify subsequent fixations to one attribute within one alternative
  infoSearch$attributeClean <- ifelse(infoSearch$attribute == shift(infoSearch$attribute, 1L) 
                                      & infoSearch$alternative == shift(infoSearch$alternative, 1L), 1, 0)
  
  #delete subsequent fixations to one attribute i.e., keep only the first fixation
  infoSearch <- infoSearch[infoSearch$attributeClean != 1 | is.na(infoSearch$attributeClean)]
  
  #delete unnecessary column
  infoSearch[, "attributeClean" := NULL]
  
  #create counter variable for alternative-wise search (focusing on set of attributes when inspecting different alternatives)
  infoSearch <- setDT(infoSearch)[, counter:= rleid(condition, trial, alternative)]
  
  #create alternative-wise substrings (i.e., sequences of letters) based on counter variable
  altwiseStrings <- infoSearch[,list(string <- paste(attribute, collapse = ""), 
                                     participant = unique(participant), 
                                     condition = unique(condition), 
                                     trial = unique(trial)), by = counter]
  
  #delete counter variable 
  altwiseStrings[, "counter" := NULL]
  
  #rename column
  setnames(altwiseStrings, "V1", "string")
  
  #apply 'relaxedFreqOrder' function to column with previosuly created alternative-wise substrings
  altwiseStrings$formattedString <- lapply(altwiseStrings$string, relaxedFreqOrder)
  
  #delete string variable  
  altwiseStrings[, "string" := NULL]
  
  #rename column
  setnames(altwiseStrings, "formattedString", "string")
  
  #change class
  altwiseStrings$string <- as.character(altwiseStrings$string)
  
  #compare substrings to identify equal subsequent alternative-wise substrings  
  altwiseStrings$equalStrings <- ifelse(altwiseStrings$string == lag(altwiseStrings$string, n = 1L) 
                                        | altwiseStrings$string == lead(altwiseStrings$string, n = 1L), 1, 0)
  
  #create counter variable based on string variable within each trial 
  #i.e., assign new number for every unique string within each trial
  altwiseStrings <- setDT(altwiseStrings)[, counter:= rleid(string, trial)]
  
  #extract equal subsequent substrings (equalStrings = 1)
  altwiseStrings <- altwiseStrings[altwiseStrings$equalStrings != 0]
  
  #delete substrings of length one
  altwiseStrings <- subset(altwiseStrings, nchar(as.character(string)) >= 2)
  
  #combine substrings into alternative-wise patterns using counter variable 
  #i.e., all substrings with equal count should be collapsed into one pattern
  altwisePatterns <- altwiseStrings[,list(string <- paste(string, collapse = ""), 
                                          participant = unique(participant), 
                                          condition = unique(condition), 
                                          trial = unique(trial)), by = counter]
  
  #delete counter variable  
  altwisePatterns[, "counter" := NULL]
  
  #rename column
  setnames(altwisePatterns, "V1", "pattern")
  
  #calculate frequencies for each pattern within each trial, condition and for every participant 
  altwisePatternsCountRan <- as.data.table(with(altwisePatterns, table(pattern, trial, condition, participant)))
  altwisePatternsCountRan <- altwisePatternsCountRan[altwisePatternsCountRan$N != 0] 
  
  return(altwisePatternsCountRan)
}

#replicate 'altwiseSim' function 10000 times ####
altwiseSimRep <- do.call(rbind, replicate(10000, altwiseSim(), simplify=FALSE)) 

#calculate probabilities and probability complements ####

#write function which compares pattern frequencies in original and simulated data sets for each participant, condition and trial
altwiseProb <- function(i){ 
  sum(altwiseSimRep$pattern == altwisePatternsCount$pattern[i] 
      & altwiseSimRep$participant == altwisePatternsCount$participant[i] 
      & altwiseSimRep$condition == altwisePatternsCount$condition[i] 
      & altwiseSimRep$trial == altwisePatternsCount$trial[i] 
      & altwiseSimRep$N >= altwisePatternsCount$pattFreq[i])
}

#apply 'altwiseProb' function
altwisePatternsCount$pattFreqSim <- sapply(1:nrow(altwisePatternsCount), altwiseProb) 

#calculate probabilities
altwisePatternsCount$probability <- altwisePatternsCount$pattFreqSim / 10000

#calculate probability complements (1 - probability)
altwisePatternsCount$prob_complement <- 1 - altwisePatternsCount$probability

#calculate pattern lengths
altwisePatternsCount$pattLength <- nchar(altwisePatternsCount$pattern)

#save table
#in case we want to perform some data analyses without doing simulation again
write.csv(file="fileName.csv", x=altwisePatternsCount) 

#identify attribute-wise patterns ####

#create attribute-wise substrings (i.e., sequences of letters) based on counter variable
#collapsing all atributes between different alternative into a string of letters 
#e.g., fixations on sugar attribute between four different alternatives => substring 'ssss'
attwiseStrings <- infoSearch[,list(string <- paste(attribute, collapse = ""), 
                                   participant = unique(participant), 
                                   condition = unique(condition), 
                                   trial = unique(trial)), by = counterAttwise] 

#delete substrings of length three or less
attwiseStrings <- subset(attwiseStrings, nchar(as.character(V1)) >= 4)

#delete counter variable  
attwiseStrings[, "counterAttwise" := NULL]

#rename column
setnames(attwiseStrings, "V1", "pattern")

#calculate frequencies for each pattern within each trial, condition and for every participant 
attwisePatternsCount <- as.data.table(with(attwiseStrings, table(pattern, trial, condition, participant)))
attwisePatternsCount <- attwisePatternsCount[attwisePatternsCount$N != 0]
setnames(attwisePatternsCount, "N", "pattFreq") #rename frequency column

#assess whether obtained patterns occurred by chance by making random data set 
#to which we will compare the patterns from the original data set ####

attwiseSim <- function() { #function which contains random version of data
  
  #import data (same file as in line 14)
  infoSearch <- as.data.table(read.csv("fileName.csv", header = T, sep = ";"))

  #delete unnecessary columns
  infoSearch[, c("alternative", "attribute") := NULL]
  
  #create random data by sampling evenly from alternatives and attributes
  sim <- 154355 #number of rows corresponding to number of fixations in original data set (total string length)
  infoSearch$alternative <- sample(1:4, sim, T) #sample numbers from 1 to 4 154355 times
  infoSearch$attribute <- sample(c("b", "f", "p", "s"), sim, T) #sample letters b, f, p and s 154355 times
  
  #identify subsequent fixations to one attribute within one alternative
  infoSearch$attributeClean <- ifelse(infoSearch$attribute == shift(infoSearch$attribute, 1L) 
                                      & infoSearch$alternative == shift(infoSearch$alternative, 1L), 1, 0)
  
  #delete subsequent fixations to one attribute i.e., keep only the first fixation
  infoSearch <- infoSearch[infoSearch$attributeClean != 1 | is.na(infoSearch$attributeClean)]
  
  #delete unnecessary column
  infoSearch[, "attributeClean" := NULL]
  
  #create counter variable for attribute-wise search (focusing on one attribute when inspecting different alternatives)
  #since additional fixations have been deleted, when there are fixations on one attribute, they must belong to different alternatives
  infoSearch <- setDT(infoSearch)[, counter:= rleid(condition, trial, attribute)] 
  
  #create attribute-wise patterns (i.e., sequences of letters) based on counter variable
  attwiseStrings <- infoSearch[,list(string <- paste(attribute, collapse = ""), 
                                     participant = unique(participant), 
                                     condition = unique(condition), 
                                     trial = unique(trial)), by = counter]
  
  #delete patterns of length three or less
  attwiseStrings <- subset(attwiseStrings, nchar(as.character(V1)) >= 4)
  
  #delete counter variable  
  attwiseStrings[, "counter" := NULL]
  
  #rename column
  setnames(attwiseStrings, "V1", "pattern")
  
  #calculate frequencies for each pattern within each trial, condition and for every participant 
  attwisePatternsCountRan <- as.data.table(with(attwiseStrings, table(pattern, trial, condition, participant)))
  attwisePatternsCountRan <- attwisePatternsCountRan[attwisePatternsCountRan$N != 0]
  
  return(attwisePatternsCountRan)
}

#replicate 'attwiseSim' function 10000 times ####

attwiseSimRep <- do.call("rbind", replicate(10000, attwiseSim(), simplify=FALSE )) 

#calculate probabilities and probability complements ####

#write function which compares pattern frequences in original and simulated data sets for each participant, condition and trial
attwiseProb <- function(i){
  sum(attwiseSimRep$pattern == attwisePatternsCount$pattern[i] 
      & attwiseSimRep$participant == attwisePatternsCount$participant[i] 
      & attwiseSimRep$condition == attwisePatternsCount$condition[i] 
      & attwiseSimRep$trial == attwisePatternsCount$trial[i] 
      & attwiseSimRep$N >= attwisePatternsCount$pattFreq[i])
}

#apply 'attwiseProb' function
attwisePatternsCount$pattFreqSim <- sapply(1:nrow(attwisePatternsCount), attwiseProb) 

attwisePatternsCount$probability <- attwisePatternsCount$pattFreqSim / 10000

#calculate probability complements
attwisePatternsCount$prob_complement <- 1 - attwisePatternsCount$probability

#calculate pattern lengths
attwisePatternsCount$pattLength <- nchar(attwisePatternsCount$pattern)

#save table
#in case we want to perform some data analyses without doing simulation again
write.csv(file="fileName.csv", x=attwisePatternsCount) 

#calculate numerator for SSI for alternative-wise patterns 
#numerator = length of each unique pattern * frequency of each unique pattern * probabilty complement
altwisePatternsCount$numerator <- altwisePatternsCount$pattFreq * 
                                  altwisePatternsCount$pattLength * 
                                  altwisePatternsCount$prob_complement
sysAltwise <- ddply(altwisePatternsCount,.(participant, condition, trial), summarize, altwiseSum = sum(numerator))

#format data
sysAltwise <- as.data.table(sysAltwise)
sysAltwise$participant <- as.numeric(sysAltwise$participant)
sysAltwise$trial <- as.numeric(sysAltwise$trial)
sysAltwise <- sysAltwise[order(participant, condition, trial),]  

#merge in total string length per participant per trial #calculated in line 61
sysAltwise <- merge(sysAltwise, stringLength, by = c("participant", "condition", "trial"), all = T)
sysAltwise[is.na(sysAltwise)] <- 0

#calculate numerator for SSI for attribute-wise patterns 
#numerator = length of each unique pattern * frequency of each unique pattern * probabilty complement
attwisePatternsCount$numerator <- attwisePatternsCount$pattFreq * 
                                  attwisePatternsCount$pattLength * 
                                  attwisePatternsCount$prob_complement
sysAttwise <- ddply(attwisePatternsCount,.(participant, condition, trial), summarize, attwiseSum = sum(numerator))

#format data
sysAttwise <- as.data.table(sysAttwise)
sysAttwise$participant <- as.numeric(sysAttwise$participant)
sysAttwise$trial <- as.numeric(sysAttwise$trial)
sysAttwise <- sysAttwise[order(participant, condition, trial),]  

#merge in total string length per participant per trial #calculated in the line 61
sysAttwise <- merge(sysAttwise, stringLength, by = c("participant", "condition", "trial"), all = T)
sysAttwise[is.na(sysAttwise)] <- 0

#calculate SSI ####
sysIndex <- merge(sysAltwise, sysAttwise, by = c("participant", "condition", "trial"), all = T)
sysIndex$N.x <- NULL #delete unnecessary column
setnames(sysIndex, "N.y", "stringLength") #rename column
sysIndex$sysIndex <- (sysIndex$altwiseSum + sysIndex$attwiseSum) / sysIndex$stringLength
