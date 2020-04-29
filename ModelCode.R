# Model Code Demo
# Filename: ModelCode.R
# Author: CA

##            ##
# Source Files #
##            ##
source("SimplifiedSampler.R")
source("PlotGraph.R")
source("Summary.R")

##            ##
# Load in Data #
##            ##
vowels <- levels(data$Target.keywords)

# Load in data
load("Recoded Data.RData")
  
vowels <- levels(data$Target.keywords)
vowel <- "goose"
data <- subset(data,Target.keywords==vowel)
  
  
# Specify terms for selection
model <- F1.mean+F2.mean+F3.mean~ gender + Decade + Age + Preceding.POA + Following.POA 
random <- ~Speaker+Target.transcript

## Simplified sample
model.1 <- Simplified_sampler(data,model,random,interactions =2,param.ex=TRUE,nesting=TRUE,iterations=1000,Graphs=FALSE)

Summarise(model.1)

PlotGraph(model.1,mult.models =FALSE)

library(lme4)
model.lme1 <- lmer(F1.mean~gender + Decade + Age + Preceding.POA + Following.POA +(1|Speaker)+(1|Target.transcript),data=data)
model.lme2 <- lmer(F2.mean~gender + Decade + Age + Preceding.POA + Following.POA +(1|Speaker)+(1|Target.transcript),data=data)
model.lme3 <- lmer(F3.mean~gender + Decade + Age + Preceding.POA + Following.POA +(1|Speaker)+(1|Target.transcript),data=data)
