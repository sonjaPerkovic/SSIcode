#import necessary libraries
library(data.table)
library(plyr)
library(dplyr)
library(gtools)
library(gridExtra)
library(ggplot2)

#set working directory
setwd("/Users/userName/folderName")

#read in the data file
infoSearch <- as.data.table(read.csv("fileName.csv", header = T, sep = ";")) #my file consisted of five columns: participant, environment, trial, alternative and attribute column

#preparing the data ####

#identify subsequent fixations to an attribute within the same alternative
infoSearch$attributeClean <- ifelse(infoSearch$attribute == shift(infoSearch$attribute, 1L) & infoSearch$alternative == shift(infoSearch$alternative, 1L), 1, 0) 

#delete subsequent fixations to an attribute i.e., keep only the first fixation
infoSearch <- infoSearch[infoSearch$attributeClean != 1 | is.na(infoSearch$attributeClean)]

#delete unnecessary column
infoSearch[, "attributeClean" := NULL]

#create and count alternative-wise transitions (needed for calculating Search Index)
infoSearch$transAlt <- ifelse(infoSearch$attribute != lag(infoSearch$attribute, n = 1L) & infoSearch$alternative == lag(infoSearch$alternative, n = 1L) & infoSearch$trial == lag(infoSearch$trial, n = 1L) & infoSearch$participant == lag(infoSearch$participant, n = 1L), 1, 0)
infoSearch[is.na(infoSearch)] <- 0
altTrans <- ddply(infoSearch,.(participant, environment, trial), summarize, transAlt = sum(transAlt)) #new file with alternative-wise transitions for each participant within each environment and for each trial

#create and count attribute-wise transitions (needed for calculating Search Index)
infoSearch$transAtt <- ifelse(infoSearch$attribute == lag(infoSearch$attribute, n = 1L) & infoSearch$alternative != lag(infoSearch$alternative, n = 1L) & infoSearch$trial == lag(infoSearch$trial, n = 1L) & infoSearch$participant == lag(infoSearch$participant, n = 1L), 1, 0)
infoSearch[is.na(infoSearch)] <- 0
attTrans <- ddply(infoSearch,.(participant, environment, trial), summarize, transAtt = sum(transAtt)) #new file with attribute-wise transitions for each participant within each environment and for each trial

#combine two data sets by columns and calculate Search Index
searchIndex <- as.data.table(cbind(attTrans, altTrans)) #combine two data sets by columns
searchIndex$searchIndex <- (searchIndex$transAlt - searchIndex$transAtt) / (searchIndex$transAlt + searchIndex$transAtt) #calculate Search Index using the apporpiate equation
searchIndex[, c("participant", "environment", "trial") := NULL] #delete unnecessary columns
setcolorder(searchIndex, c("participant", "environment", "trial", "transAlt", "transAtt", "searchIndex")) #set order of columns 

#calculate the length of total string of eye fixations per participant per trial (needed for calculating the denominator of Systematicity of Search Index)
stringLength <- ddply(infoSearch, .(participant, environment, trial), function(infoSearch) length(infoSearch$attribute))
setnames(stringLength, "V1", "N") #change the name of a column

#create counter variable for alternative-wise search (focusing on a set of attributes when inspecting different alternatives)
infoSearch <- setDT(infoSearch)[, counterAltwise:= rleid(environment, trial, alternative)] #assigning the same number to the fixations to the attributes within the same alternatives (e.g., if a participant first fixated on sugar and fat level within one alternative and then sugar and protein level within another alternative, values 1,1,2,2 would have been assigned to the counter variable) 

#create counter variable for attribute-wise search (focusing on the same attribute when inspecting different alternatives)
infoSearch <- setDT(infoSearch)[, counterAttwise:= rleid(environment, trial, attribute)] #since additional fixations have been deleted, when there is a fixation on the same attribute, it must belong to a different alternative

#identify alternative-wise patterns ####

#create alternative-wise strings (i.e., sequences of letters) based on counter variable
altwiseStrings <- infoSearch[,list(string <- paste(attribute, collapse = ""), participant = unique(participant), environment = unique(environment), trial = unique(trial)), by = counterAltwise] #collapsing all atributes within the same alternative into a string of letters; in the above example we would and up with two strings of length two: 'fs" and 'ps' (i.e., 'fat and sugar' and 'protein and sugar')

#delete counter variable 
altwiseStrings[, "counterAltwise" := NULL]

#rename column
setnames(altwiseStrings, "V1", "string")

#define a function that keeps the unique elements in a string and sorts them alphabetically
relaxedFreqOrder <- function(i){
  paste0(unique(sort(unlist(strsplit(i, "")))), collapse = "")
}

#apply the function to the column with previosuly created alternative-wise strings
altwiseStrings$formattedString <- lapply(altwiseStrings$string, relaxedFreqOrder)

#delete string variable 
altwiseStrings[, "string" := NULL]

#rename column
setnames(altwiseStrings, "formattedString", "string")

#change the type of a column into character
altwiseStrings$string <- as.character(altwiseStrings$string)

#compare strings to identify equal subsequent alternative-wise strings  
altwiseStrings$equalStrings <- ifelse(altwiseStrings$string == lag(altwiseStrings$string, n = 1L) | altwiseStrings$string == lead(altwiseStrings$string, n = 1L), 1, 0)

#create a counter variable based on string variable within each trial (i.e., assign a new number for every unique string within each trial)
altwiseStrings <- setDT(altwiseStrings)[, counter:= rleid(string, trial)]

#extract equal subsequent strings (equalStrings = 1)
altwiseStrings <- altwiseStrings[altwiseStrings$equalStrings != 0]

#delete strings of length one
altwiseStrings <- subset(altwiseStrings, nchar(as.character(string)) >= 2)

#combine strings into alternative-wise patterns using the counter variable (i.e., all the strings with the same count should be collapsed into a pattern)
altwisePatterns <- altwiseStrings[,list(string <- paste(string, collapse = ""), participant = unique(participant), environment = unique(environment), trial = unique(trial)), by = counter]

#delete counter variable  
altwisePatterns[, "counter" := NULL]

#rename column
setnames(altwisePatterns, "V1", "pattern")

#calculate the frequency of occurrence for each pattern within each trial, environment and for every participant 
altwisePatternsCount <- as.data.table(with(altwisePatterns, table(pattern, trial, environment, participant)))
altwisePatternsCount <- altwisePatternsCount[altwisePatternsCount$N != 0] 
setnames(altwisePatternsCount, "N", "pattFreq") #rename the frequency column

#assess whether obtained patterns occurred by chance by making a random data set to which we will compare the patterns from the original data set ####

altwiseSim <- function() { #creating a function which will contain the random version of the data set
  
  #read in the data file (the same file as in the line 12)
  infoSearch <- as.data.table(read.csv("fileName.csv", header = T, sep = ";"))

  #delete unnecessary columns
  infoSearch[, c("alternative", "attribute") := NULL]
  
  #randomize data
  sim <- 154355 #the number of rows corresponding to the number of fixations made in the original data set
  infoSearch$alternative <- sample(1:4, sim, T) #sample the numbers from 1 to 4 154355 times
  infoSearch$attribute <- sample(c("b", "f", "p", "s"), sim, T) #sample the letters b, f, p and s 154355 times

  #identify subsequent fixations to an attribute within the same alternative
  infoSearch$attributeClean <- ifelse(infoSearch$attribute == shift(infoSearch$attribute, 1L) & infoSearch$alternative == shift(infoSearch$alternative, 1L), 1, 0)
  
  #delete subsequent fixations to an attribute i.e., keep only the first fixation
  infoSearch <- infoSearch[infoSearch$attributeClean != 1 | is.na(infoSearch$attributeClean)]
  
  #delete unnecessary column
  infoSearch[, "attributeClean" := NULL]
  
  #create counter variable for alternative-wise search (focusing on a set of attributes when inspecting different alternatives)
  infoSearch <- setDT(infoSearch)[, counter:= rleid(environment, trial, alternative)]
  
  #create alternative-wise strings (i.e., sequences of letters) based on counter variable
  altwiseStrings <- infoSearch[,list(string <- paste(attribute, collapse = ""), participant = unique(participant), environment = unique(environment), trial = unique(trial)), by = counter]
  
  #delete counter variable 
  altwiseStrings[, "counter" := NULL]
  
  #rename column
  setnames(altwiseStrings, "V1", "string")
  
  #apply the 'relaxedFreqOrder' function to the column with previosuly created alternative-wise strings
  altwiseStrings$formattedString <- lapply(altwiseStrings$string, relaxedFreqOrder)
  
  #delete string variable  
  altwiseStrings[, "string" := NULL]
  
  #rename column
  setnames(altwiseStrings, "formattedString", "string")
  
  #change the type of a column into character
  altwiseStrings$string <- as.character(altwiseStrings$string)
  
  #compare strings to identify equal subsequent alternative-wise strings  
  altwiseStrings$equalStrings <- ifelse(altwiseStrings$string == lag(altwiseStrings$string, n = 1L) | altwiseStrings$string == lead(altwiseStrings$string, n = 1L), 1, 0)
  
  #create a counter variable based on string variable within each trial (i.e., assign a new number for every unique string within each trial)
  altwiseStrings <- setDT(altwiseStrings)[, counter:= rleid(string, trial)]
  
  #extract equal subsequent strings (equalStrings = 1)
  altwiseStrings <- altwiseStrings[altwiseStrings$equalStrings != 0]
  
  #delete strings of length one
  altwiseStrings <- subset(altwiseStrings, nchar(as.character(string)) >= 2)
  
  #combine strings into alternative-wise patterns using the counter variable (i.e., all the strings with the same count should be collapsed into a pattern)
  altwisePatterns <- altwiseStrings[,list(string <- paste(string, collapse = ""), participant = unique(participant), environment = unique(environment), trial = unique(trial)), by = counter]
  
  #delete counter variable  
  altwisePatterns[, "counter" := NULL]
  
  #rename column
  setnames(altwisePatterns, "V1", "pattern")
  
  #calculate the frequency of occurrence for each pattern within each trial, environment and for every participant 
  altwisePatternsCountRan <- as.data.table(with(altwisePatterns, table(pattern, trial, environment, participant)))
  altwisePatternsCountRan <- altwisePatternsCountRan[altwisePatternsCountRan$N != 0] 
  
  return(altwisePatternsCountRan)
}

#replicate the 'altwiseSim' function 10000 times ####
altwiseSimRep <- do.call(rbind, replicate(10000, altwiseSim(), simplify=FALSE)) 

#calculate the probabilities and probability complements ####

#write a function which compares the pattern frequences in original and simulated data sets for each participant, environment and trial
altwiseProb <- function(i){ 
  sum(altwiseSimRep$pattern == altwisePatternsCount$pattern[i] & altwiseSimRep$participant == altwisePatternsCount$participant[i] & altwiseSimRep$environment == altwisePatternsCount$environment[i] & altwiseSimRep$trial == altwisePatternsCount$trial[i] & altwiseSimRep$N >= altwisePatternsCount$pattFreq[i])
}

#apply the 'altwiseProb' function
altwisePatternsCount$pattFreqSim <- sapply(1:nrow(altwisePatternsCount), altwiseProb) 

#calculate the probabilities
altwisePatternsCount$probability <- altwisePatternsCount$pattFreqSim / 10000

#calculate the probability complements (1 - probability)
altwisePatternsCount$prob_complement <- 1 - altwisePatternsCount$probability

#calculate the pattern length
altwisePatternsCount$pattLength <- nchar(altwisePatternsCount$pattern)

#save the table
write.csv(file="fileName.csv", x=altwisePatternsCount) #in case we want to perform some data analysis without doing the simulation again

#identify attribute-wise patterns ####

#create attribute-wise strings (i.e., sequences of letters) based on counter variable
attwiseStrings <- infoSearch[,list(string <- paste(attribute, collapse = ""), participant = unique(participant), environment = unique(environment), trial = unique(trial)), by = counterAttwise] #collapsing all atributes between different alternative into a string of letters; for instance, if a participant inspected sugar attribute between four different alternatives, we would end up with a string 'ssss'

#delete strings of length three or less
attwiseStrings <- subset(attwiseStrings, nchar(as.character(V1)) >= 4)

#delete counter variable  
attwiseStrings[, "counterAttwise" := NULL]

#rename column
setnames(attwiseStrings, "V1", "pattern")

#calculate the frequency of occurrence for each pattern within each trial, environment and for every participant 
attwisePatternsCount <- as.data.table(with(attwiseStrings, table(pattern, trial, environment, participant)))
attwisePatternsCount <- attwisePatternsCount[attwisePatternsCount$N != 0]
setnames(attwisePatternsCount, "N", "pattFreq") #rename the frequency column

#assess whether obtained patterns occurred by chance by making a random data set to which we will compare the patterns from the original data set ####

attwiseSim <- function() { #creating a function which will contain the random version of the data set
  
  #reading in the data file (the same file as in the line 12)
  infoSearch <- as.data.table(read.csv("fileName.csv", header = T, sep = ";"))

  #delete unnecessary columns
  infoSearch[, c("alternative", "attribute") := NULL]
  
  #randomizing data
  sim <- 154355 #the number of rows corresponding to the number of fixations made in the original data set
  infoSearch$alternative <- sample(1:4, sim, T) #sample the numbers from 1 to4 154355 times
  infoSearch$attribute <- sample(c("b", "f", "p", "s"), sim, T) #sample the letters b, f, p and s 154355 times
  
  #identify subsequent fixations to an attribute within the same alternative
  infoSearch$attributeClean <- ifelse(infoSearch$attribute == shift(infoSearch$attribute, 1L) & infoSearch$alternative == shift(infoSearch$alternative, 1L), 1, 0)
  
  #delete subsequent fixations to an attribute i.e., keep only the first fixation
  infoSearch <- infoSearch[infoSearch$attributeClean != 1 | is.na(infoSearch$attributeClean)]
  
  #delete unnecessary column
  infoSearch[, "attributeClean" := NULL]
  
  #create counter variable for attribute-wise search (focusing on the same attribute when inspecting different alternatives)
  infoSearch <- setDT(infoSearch)[, counter:= rleid(environment, trial, attribute)] #since additional fixations have been deleted, when there is a fixation on the same attribute, it must belong to a different alternative
  
  #create attribute-wise patterns (i.e., sequences of letters) based on counter variable
  attwiseStrings <- infoSearch[,list(string <- paste(attribute, collapse = ""), participant = unique(participant), environment = unique(environment), trial = unique(trial)), by = counter]
  
  #delete patterns of length three or less
  attwiseStrings <- subset(attwiseStrings, nchar(as.character(V1)) >= 4)
  
  #delete counter variable  
  attwiseStrings[, "counter" := NULL]
  
  #rename column
  setnames(attwiseStrings, "V1", "pattern")
  
  #calculate the frequency of occurrence for each pattern within each trial, environment and for every participant 
  attwisePatternsCountRan <- as.data.table(with(attwiseStrings, table(pattern, trial, environment, participant)))
  attwisePatternsCountRan <- attwisePatternsCountRan[attwisePatternsCountRan$N != 0]
  
  return(attwisePatternsCountRan)
}

#replicate the 'attwiseSim' function 10000 times ####

attwiseSimRep <- do.call("rbind", replicate(10000, attwiseSim(), simplify=FALSE )) 

#calculate the probabilities and probability complements ####

#write a function which compares the pattern frequences in original and simulated data sets for each participant, environment and trial
attwiseProb <- function(i){
  sum(attwiseSimRep$pattern == attwisePatternsCount$pattern[i] & attwiseSimRep$participant == attwisePatternsCount$participant[i] & attwiseSimRep$environment == attwisePatternsCount$environment[i] & attwiseSimRep$trial == attwisePatternsCount$trial[i] & attwiseSimRep$N >= attwisePatternsCount$pattFreq[i])
}

#apply the 'attwiseProb' function
attwisePatternsCount$pattFreqSim <- sapply(1:nrow(attwisePatternsCount), attwiseProb) 

attwisePatternsCount$probability <- attwisePatternsCount$pattFreqSim / 10000

#calculate the probability complement
attwisePatternsCount$prob_complement <- 1 - attwisePatternsCount$probability

#calculate the pattern length
attwisePatternsCount$pattLength <- nchar(attwisePatternsCount$pattern)

#save the table
write.csv(file="fileName.csv", x=attwisePatternsCount) #in case we want to perform some data analysis without doing the simulation again

#calculate numerator for Systematicity of Search Index for alternative-wise patterns (numerator = length of each unique pattern * frequency of each unique pattern * probabilty complement)
altwisePatternsCount$numerator <- altwisePatternsCount$pattFreq * altwisePatternsCount$pattLength * altwisePatternsCount$prob_complement
sysAltwise <- ddply(altwisePatternsCount,.(participant, environment, trial), summarize, altwiseSum = sum(numerator))

#format the data
sysAltwise <- as.data.table(sysAltwise)
sysAltwise$participant <- as.numeric(sysAltwise$participant)
sysAltwise$trial <- as.numeric(sysAltwise$trial)
sysAltwise <- sysAltwise[order(participant, environment, trial),]  

#merge in the string length (eye fixations of the entire sample) #calculated in the line 41
sysAltwise <- merge(sysAltwise, stringLength, by = c("participant", "environment", "trial"), all = T)
sysAltwise[is.na(sysAltwise)] <- 0

#calculate numerator for Systematicity of Search Index for attribute-wise patterns (numerator = length of each unique pattern * frequency of each unique pattern * probabilty complement)
attwisePatternsCount$numerator <- attwisePatternsCount$pattFreq * attwisePatternsCount$pattLength * attwisePatternsCount$prob_complement
sysAttwise <- ddply(attwisePatternsCount,.(participant, environment, trial), summarize, attwiseSum = sum(numerator))

#format the data
sysAttwise <- as.data.table(sysAttwise)
sysAttwise$participant <- as.numeric(sysAttwise$participant)
sysAttwise$trial <- as.numeric(sysAttwise$trial)
sysAttwise <- sysAttwise[order(participant, environment, trial),]  

#merge in the string length (eye fixations of the entire sample) #calculated in the line 41
sysAttwise <- merge(sysAttwise, stringLength, by = c("participant", "environment", "trial"), all = T)
sysAttwise[is.na(sysAttwise)] <- 0

#calculate Systematicity of Search Index ####
sysIndex <- merge(sysAltwise, sysAttwise, by = c("participant", "environment", "trial"), all = T)
sysIndex$N.x <- NULL #delete unnecessary column
setnames(sysIndex, "N.y", "stringLength") #rename column
sysIndex$sysIndex <- (sysIndex$altwiseSum + sysIndex$attwiseSum) / sysIndex$stringLength
