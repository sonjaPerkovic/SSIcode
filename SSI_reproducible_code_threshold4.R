#import libraries
library(data.table)
library(dplyr)
library(parallel)

#set working directory
setwd("/Users/sonjaperkovic/Dropbox/SSIreproducibleCode") #for Mac users
#setwd("C:\\Users\\sonjaperkovic\\Dropbox\\SSIreproducibleCode") #for PC users

#import data
#the file consists of five columns: participant, condition, trial, alternative and attribute
infoSearch = as.data.table(read.csv("gaborClean.csv", header = T, sep = ";"))

#data preparation ####

#identify subsequent fixations (refixations) to one attribute within one alternative
infoSearch$attributeClean = ifelse(infoSearch$trial == shift(infoSearch$trial, 1L)
                                  & infoSearch$attribute == shift(infoSearch$attribute, 1L)
                                  & infoSearch$alternative == shift(infoSearch$alternative, 1L), 1, 0)

#delete refixations, i.e. keep only first fixation
infoSearch = infoSearch[infoSearch$attributeClean != 1 | is.na(infoSearch$attributeClean)]

#delete unnecessary column
infoSearch[, "attributeClean" := NULL]

#identify alternative-wise transitions (needed for calculating SI)
n = nrow(infoSearch)
infoSearch$transAlt = c(0L, infoSearch$trial[-n] == infoSearch$trial[-1] & infoSearch$alternative[-n] == infoSearch$alternative[-1])

#new object with alternative-wise transitions for each participant, condition and trial
setkey(infoSearch, "participant", "condition", "trial")
altTrans = infoSearch[, list(transAlt=sum(transAlt)), by=key(infoSearch)] 

#identify attribute-wise transitions (needed for calculating SI)
infoSearch$transAtt = c(0L, infoSearch$trial[-n] == infoSearch$trial[-1] & infoSearch$attribute[-n] == infoSearch$attribute[-1])

#new object with attribute-wise transitions for each participant, condition and trial
attTrans = infoSearch[, list(transAtt=sum(transAtt)), by=key(infoSearch)] 

#merge two data sets
searchIndex = merge(attTrans, altTrans, by = c("participant", "condition", "trial"))

#calculate SI
searchIndex$searchIndex = (searchIndex$transAlt - searchIndex$transAtt) / (searchIndex$transAlt + searchIndex$transAtt)

#calculate total string length of eye fixations per participant per trial (needed for SSI denominator)
stringLength = infoSearch[, list(N=NROW(attribute)), by=key(infoSearch)] 

#create counter variable which will be used for concatenating fixations into substrings
infoSearch$transDiff = c(0L, infoSearch$transAtt[-n] != infoSearch$transAtt[-1] & infoSearch$transAlt[-n] != infoSearch$transAlt[-1])
infoSearch$diff1 <- c(NA, diff(infoSearch$transDiff)) 

infoSearch$counter = c(1L, infoSearch$trial[-n] != infoSearch$trial[-1] 
                       | c(infoSearch$trial[-n] == infoSearch$trial[-1] & 
                             infoSearch$alternative[-n] != infoSearch$alternative[-1] & 
                             infoSearch$attribute[-n] != infoSearch$attribute[-1]) 
                       
                       | c(infoSearch$trial[-n] == infoSearch$trial[-1] & 
                             infoSearch$alternative[-n] != infoSearch$alternative[-1] & 
                             infoSearch$attribute[-n] == infoSearch$attribute[-1] &
                             infoSearch$transAtt[-n] != infoSearch$transAtt[-1] & 
                             infoSearch$transAlt[-n] != infoSearch$transAlt[-1] &
                             infoSearch$transDiff[-n] != infoSearch$transDiff[-1] & 
                             infoSearch$diff1[-n] != infoSearch$diff1[-1])
                       
                       | c(infoSearch$trial[-n] == infoSearch$trial[-1] & 
                             infoSearch$alternative[-n] != infoSearch$alternative[-1] & 
                             infoSearch$attribute[-n] == infoSearch$attribute[-1] &
                             infoSearch$transAtt[-n] != infoSearch$transAtt[-1] & 
                             infoSearch$transAlt[-n] != infoSearch$transAlt[-1] & 
                             infoSearch$transDiff[-n] == infoSearch$transDiff[-1] & 
                             infoSearch$diff1[-n] == infoSearch$diff1[-1])
                       
                       | c(infoSearch$trial[-n] == infoSearch$trial[-1] & 
                             infoSearch$alternative[-n] == infoSearch$alternative[-1] & 
                             infoSearch$attribute[-n] != infoSearch$attribute[-1] &
                             infoSearch$transAtt[-n] != infoSearch$transAtt[-1] & 
                             infoSearch$transAlt[-n] != infoSearch$transAlt[-1] & 
                             infoSearch$transDiff[-n] != infoSearch$transDiff[-1]))
infoSearch$counter = cumsum(infoSearch$counter)

#identify alternative- and attribute-wise patterns ####

#create alternative- and attribute-wise substrings based on counter variable
substrings = infoSearch[, list(string = paste(attribute, collapse = ""),
                               participant = unique(participant),
                               condition = unique(condition),
                               trial = unique(trial)), by = counter]

#function that identifies rows with identical elements only
identElements = function(i){
  length(unique(unlist(strsplit(i, "")))) == 1
}

substrings$identElements = lapply(substrings$string, identElements)

#subset data with alternative-wise substrings only
altwiseSubstrings = subset(substrings, identElements == "FALSE")

#function that keeps unique elements and sorts them alphabetically
relaxedFreqOrder = function(i){
  paste0(unique(sort(unlist(strsplit(i, "")))), collapse = "")
}

altwiseSubstrings$string = lapply(altwiseSubstrings$string, relaxedFreqOrder)

#change class
altwiseSubstrings$string = as.character(altwiseSubstrings$string)
altwiseSubstrings$identElements = as.numeric(altwiseSubstrings$identElements)
substrings$string = as.character(substrings$string)
substrings$identElements = as.numeric(substrings$identElements)

#merge in the formatted data set with altwise substrings
substrings = anti_join(substrings, altwiseSubstrings, by = c("counter", "participant", "condition", "trial", "identElements")) %>% 
             bind_rows(altwiseSubstrings)
substrings = substrings[order(substrings$counter),]

#create counter variable for alternative-wise substrings based on string variable within each trial
substrings = setDT(substrings)[, countEqualSubstrings := rleid(string, trial)]

#create variable that assigns 1 to subsequent equal counter variable values
substrings$equalCounter <- ifelse(substrings$countEqualSubstrings == lag(substrings$countEqualSubstrings, n = 1L) |
                                  substrings$countEqualSubstrings == lead(substrings$countEqualSubstrings, n = 1L), 1, 0)

#extract strings
substrings = substrings[substrings$equalCounter != 0 | substrings$identElements != 0]

#combine substrings into patterns using counter variable
#i.e., all substrings with equal count should be collapsed into one pattern
patterns = substrings[,list(pattern = paste(string, collapse = ""),
                            participant = unique(participant),
                            condition = unique(condition),
                            trial = unique(trial)), by = countEqualSubstrings]

#keep substrings of minimum length four
patterns = subset(patterns, nchar(as.character(pattern)) >= 4)

#calculate frequencies for each pattern within each trial, condition and for every participant
patternsCount = setDT(patterns)[, list(pattFreq = .N), by = c('pattern', 'trial', 'condition', 'participant')]

#assess whether obtained patterns occurred by chance ####
#create random data set and compare it to the original data set

patternsSim = function() { #function which creates and analyses random data

  #import data again
  infoSearch = as.data.table(read.csv("gaborClean.csv", header = T, sep = ";"))

  #delete unnecessary columns
  infoSearch[, c("alternative", "attribute") := NULL]

  #create random data by sampling evenly from alternatives and attributes
  sim = nrow(infoSearch) #number of rows corresponding to number of fixations in original data set (total string length)
  num_alt = 4 #number of alternatives in the experiment
  infoSearch$alternative = sample(1:num_alt, sim, T) #sample numbers from 1 to 4
  infoSearch$attribute = sample(c("b", "f", "p", "s"), sim, T) #sample letters

  #identify refixations within each alternative
  infoSearch$attributeClean = ifelse(infoSearch$trial == shift(infoSearch$trial, 1L)
                                     & infoSearch$attribute == shift(infoSearch$attribute, 1L)
                                     & infoSearch$alternative == shift(infoSearch$alternative, 1L), 1, 0)
  
  #delete refixations i.e., keep only first fixation
  infoSearch = infoSearch[infoSearch$attributeClean != 1 | is.na(infoSearch$attributeClean)]

  #delete unnecessary column
  infoSearch[, "attributeClean" := NULL]
  
  #identify alternative- and attribute-wise transitions 
  n = nrow(infoSearch)
  infoSearch$transAlt = c(0L, infoSearch$trial[-n] == infoSearch$trial[-1] & infoSearch$alternative[-n] == infoSearch$alternative[-1])
  infoSearch$transAtt = c(0L, infoSearch$trial[-n] == infoSearch$trial[-1] & infoSearch$attribute[-n] == infoSearch$attribute[-1])

  #create counter variable which will be used for concatenating fixations into substrings
  infoSearch$transDiff = c(0L, infoSearch$transAtt[-n] != infoSearch$transAtt[-1] & infoSearch$transAlt[-n] != infoSearch$transAlt[-1])
  infoSearch$diff1 <- c(NA, diff(infoSearch$transDiff)) 
  
  infoSearch$counter = c(1L, infoSearch$trial[-n] != infoSearch$trial[-1] 
                         | c(infoSearch$trial[-n] == infoSearch$trial[-1] & 
                               infoSearch$alternative[-n] != infoSearch$alternative[-1] & 
                               infoSearch$attribute[-n] != infoSearch$attribute[-1]) 
                         
                         | c(infoSearch$trial[-n] == infoSearch$trial[-1] & 
                               infoSearch$alternative[-n] != infoSearch$alternative[-1] & 
                               infoSearch$attribute[-n] == infoSearch$attribute[-1] &
                               infoSearch$transAtt[-n] != infoSearch$transAtt[-1] & 
                               infoSearch$transAlt[-n] != infoSearch$transAlt[-1] &
                               infoSearch$transDiff[-n] != infoSearch$transDiff[-1] & 
                               infoSearch$diff1[-n] != infoSearch$diff1[-1])
                         
                         | c(infoSearch$trial[-n] == infoSearch$trial[-1] & 
                               infoSearch$alternative[-n] != infoSearch$alternative[-1] & 
                               infoSearch$attribute[-n] == infoSearch$attribute[-1] &
                               infoSearch$transAtt[-n] != infoSearch$transAtt[-1] & 
                               infoSearch$transAlt[-n] != infoSearch$transAlt[-1] & 
                               infoSearch$transDiff[-n] == infoSearch$transDiff[-1] & 
                               infoSearch$diff1[-n] == infoSearch$diff1[-1])
                         
                         | c(infoSearch$trial[-n] == infoSearch$trial[-1] & 
                               infoSearch$alternative[-n] == infoSearch$alternative[-1] & 
                               infoSearch$attribute[-n] != infoSearch$attribute[-1] &
                               infoSearch$transAtt[-n] != infoSearch$transAtt[-1] & 
                               infoSearch$transAlt[-n] != infoSearch$transAlt[-1] & 
                               infoSearch$transDiff[-n] != infoSearch$transDiff[-1]))
  infoSearch$counter = cumsum(infoSearch$counter)
  
  #create alternative- and attribute-wise substrings based on counter variable
  substrings = infoSearch[, list(string = paste(attribute, collapse = ""),
                                 participant = unique(participant),
                                 condition = unique(condition),
                                 trial = unique(trial)), by = counter]
  
  #function that identifies rows with identical elements only 
  identElements = function(i){
    length(unique(unlist(strsplit(i, "")))) == 1
  }
  
  substrings$identElements = lapply(substrings$string, identElements)
  
  #subset data with alternative-wise substrings only
  altwiseSubstrings = subset(substrings, identElements == "FALSE")
  
  #function that keeps unique elements and sorts them alphabetically 
  relaxedFreqOrder = function(i){
    paste0(unique(sort(unlist(strsplit(i, "")))), collapse = "")
  }
  
  altwiseSubstrings$string = lapply(altwiseSubstrings$string, relaxedFreqOrder)
  
  #change class
  altwiseSubstrings$string = as.character(altwiseSubstrings$string)
  altwiseSubstrings$identElements = as.numeric(altwiseSubstrings$identElements)
  substrings$string = as.character(substrings$string)
  substrings$identElements = as.numeric(substrings$identElements)
  
  #merge in the formatted data set with altwise substrings
  substrings = anti_join(substrings, altwiseSubstrings, by = c("counter", "participant", "condition", "trial", "identElements")) %>%
               bind_rows(altwiseSubstrings)
  substrings = substrings[order(substrings$counter),]
  
  #create counter variable for alternative-wise substrings based on string variable within each trial
  substrings = setDT(substrings)[, countEqualSubstrings := rleid(string, trial)]
  
  #create variable that assigns 1 to subsequent equal counter variable values
  substrings$equalCounter <- ifelse(substrings$countEqualSubstrings == lag(substrings$countEqualSubstrings, n = 1L) | 
                                      substrings$countEqualSubstrings == lead(substrings$countEqualSubstrings, n = 1L), 1, 0)
  
  #extract strings 
  substrings = substrings[substrings$equalCounter != 0 | substrings$identElements != 0]

  #combine substrings into patterns using counter variable
  #i.e., all substrings with equal count should be collapsed into one pattern
  patterns = substrings[,list(pattern = paste(string, collapse = ""),
                              participant = unique(participant),
                              condition = unique(condition),
                              trial = unique(trial)), by = countEqualSubstrings]
  
  #keep substrings of minimum length four
  patterns = subset(patterns, nchar(as.character(pattern)) >= 4)
  
  #calculate frequencies for each pattern within each trial, condition and for every participant
  patternsCountRan = setDT(patterns)[, list(N = .N), by = c('pattern', 'trial', 'condition', 'participant')]
  
  return(patternsCountRan)
}

#replicate 'patternsSim' function 1000 times ####

patternsSimRep = do.call(rbind, mclapply(1:1000, function(i) patternsSim()))

#calculate probabilities and probability complements ####

#function which compares pattern frequencies in original and simulated data sets for each participant, condition and trial
patternsProb = function(i){
  sum(patternsSimRep$pattern == patternsCount$pattern[i]
      & patternsSimRep$participant == patternsCount$participant[i]
      & patternsSimRep$condition == patternsCount$condition[i]
      & patternsSimRep$trial == patternsCount$trial[i]
      & patternsSimRep$N >= patternsCount$pattFreq[i])
}

#apply 'patternsProb' function
patternsCount$pattFreqSim = sapply(1:nrow(patternsCount), patternsProb)

#calculate probabilities
patternsCount$probability = patternsCount$pattFreqSim / 1000

#calculate probability complements (1 - probability)
patternsCount$prob_complement = 1 - patternsCount$probability

#calculate pattern lengths
patternsCount$pattLength = nchar(patternsCount$pattern)

#save table
#in case we want to perform some data analyses without doing simulation again
write.csv(file="fileName.csv", x = patternsCount)

#calculate numerator for SSI 
#numerator = length of each unique pattern * frequency of each unique pattern * probabilty complement
patternsCount$numerator = patternsCount$pattLength *
                          patternsCount$pattFreq *
                          patternsCount$prob_complement

setkey(patternsCount, "participant", "condition", "trial")
sysIndex = patternsCount[, list(patternSum=sum(numerator)), by=key(patternsCount)] 

#format data
sysIndex = as.data.table(sysIndex)
sysIndex$participant = as.numeric(sysIndex$participant)
sysIndex$trial = as.numeric(sysIndex$trial)
sysIndex = sysIndex[order(participant, condition, trial),]

#merge in total string length per participant per trial #calculated in line 48
sysIndex = merge(sysIndex, stringLength, by = c("participant", "condition", "trial"), all = T)
sysIndex[is.na(sysIndex)] = 0

#calculate SSI ####
sysIndex$sysIndex <- sysIndex$patternSum / sysIndex$N

#merge and format data sets with both indecis
bothIndecis = merge(sysIndex, searchIndex, by=c("participant", "condition", "trial"), all = T)
is.nan.data.frame = function(x) do.call(cbind, lapply(x, is.nan))
bothIndecis[is.nan(bothIndecis)] = 0


