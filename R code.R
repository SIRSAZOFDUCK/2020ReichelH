## Impact of guidance publication on primary care prescribing rates of simple analgesia: an interrupted time series analysis in England   
## Hannah Reichel & Saran Shantikumar (R1 August 2020)
## Contact: saran.shantikumar@warwick.ac.uk
## Template script

## ATTACH PACKAGES -------------------

rm(list = ls())

list.of.packages <- c("foreign","tsModel","lmtest","Epi","splines","vcd","data.table","dplyr","plyr","RColorBrewer","ggplot2","plotly","Cairo","MASS","sandwich","jtools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# load the packages
library(foreign)
library(tsModel)
library(lmtest)
library(Epi) 
library(splines)
library(vcd)
library(data.table)
library(dplyr)
library(plyr)
library(RColorBrewer)
library(ggplot2)
library(plotly)
library(Cairo)
library(MASS)
library(sandwich)
library(jtools)


## DEFINE DATA AND OUTPUT DIRECTORIES -----------------

# Define analysis folder
dirAnalysis <- "C:/Users/Visitor/Analysis_topical"

# Define data folder
dirData <- "C:/Users/Visitor/Data"

# Set WD to analysis folder
setwd(dirAnalysis)

## FILTER PRESCRIBING DATA FOR BNF CODES OF INTEREST ----------------

# read in BNF codes of interest
keep <- read.csv("bnf_codes_drugclass.csv") # read in codes of drug of interest
keep <- as.vector(keep$bnf_code) # make list a vector
keep <- substr(keep,1,15) # keep first 15 letters of BNF code - this is the full BNF code, but occasionally extra characters are appended

# read in .csv file of GP prescription data for one month
setwd(dirData)
data <- fread("2017_01.csv", header=T,sep = ',')

# filter perscriptions by BNF codes of interest
data1 <- data[data$`BNF CODE` %in% keep, ] # filter by BNF codes in "keep"
setwd(dirAnalysis)
write.csv(data1, "2017_01_filtered.csv", row.names=F) # save

# repeat for all other included monthly datasets, for example:

setwd(dirData)
data <- fread("2017_02.csv", header=T,sep = ',')
data1 <- data[data$`BNF CODE` %in% keep, ] # filter by BNF codes in "keep"
setwd(dirAnalysis)
write.csv(data1, "2017_02_filtered.csv", row.names=F) # save

# read in all filtered monthly data frames 

rm(list=ls())
setwd(dirAnalysis)

data201501 <- read.csv("2015_01_filtered.csv")[,-11] # read in all filtered monthly data, GET RID OF EMPTY COLUMN v11 where this exists (up to and including Nov 2018)
data201502 <- read.csv("2015_02_filtered.csv")[,-11]
data201503 <- read.csv("2015_03_filtered.csv")[,-11]
data201504 <- read.csv("2015_04_filtered.csv")[,-11]

  # ... etc. for all monthly filtered datasets

## LINK OTHER PRACTICE DATA AND AGGREGATE BY MONTH ----------------

# Add monthly list size data to monthly prescribing data

setwd(dirData)
listsize201701 <- read.csv("2017_01_listsize.csv")[,c(1,9)] #only read in practice code and list size (this does not include OOH services)
names(listsize201701) <- c("PRACTICE", "LIST.SIZE") # rename columns so they are the same as the prescribing dataset

# Aggregate number of items and cost for monthly data where practice and period is the same
data201701 <- ddply(data201701,.(PRACTICE, PERIOD),summarise,ITEMS = sum(ITEMS),ACT.COST = sum(ACT.COST))

# Add list size data by practice
data201701 <- merge(data201701,listsize201701,by="PRACTICE") # add list size by practice

# Do the same for each monthly dataset, for example:

listsize201702 <- read.csv("2017_02_listsize.csv")[,c(1,9)]
names(listsize201702) <- c("PRACTICE", "LIST.SIZE")
data201702 <- ddply(data201702,.(PRACTICE, PERIOD),summarise,ITEMS = sum(ITEMS),ACT.COST = sum(ACT.COST))
data201702 <- merge(data201702,listsize201702,by="PRACTICE")

# combine monthly datasets into single data frame and save

setwd(dirAnalysis)
alldata <- rbind(data201501,data201502,data201503,data201504,data201505,data201506,data201507,data201508,data201509,data201510,data201511,data201512,data201601,data201602,data201603,data201604,data201605,data201606,data201607,data201608,data201609,data201610,data201611,data201612,data201701,data201702,data201703,data201704,data201705,data201706,data201707,data201708,data201709,data201710,data201711,data201712,data201801,data201802,data201803,data201804,data201805,data201806,data201807,data201808,data201809,data201810,data201811,data201812,data201901,data201902,data201903)
write.csv(alldata,"allRxdata_filtered.csv",row.names = F)

# aggregate number of items and cost, and number of patients where period is the same (thus summing all prescriptions irrespective of practice)

alldata_sum <- ddply(alldata,.(PERIOD),summarise,ITEMS = sum(ITEMS),ACT.COST = sum(ACT.COST),LIST.SIZE = sum(LIST.SIZE))
alldata_sum$ITEMS.PER.1000 <- 1000*alldata_sum$ITEMS/alldata_sum$LIST.SIZE # calculate prescribing rate per 1000 registered patients
alldata_sum$time <- 1:nrow(alldata_sum) #add "time" starting at 1 from month 1
alldata_sum$year <- substr(alldata_sum$PERIOD,1,4) # add column of year (based on period)
alldata_sum$month <- substr(alldata_sum$PERIOD,5,6) # add column of month (based on period)
alldata_sum$intervention <- 0 #starts column of zeros for intervention dummy variable
alldata_sum$intervention[c(40:51)] <- 1 #add 1 to rows which are after intervention
write.csv(alldata_sum,"allRxData_filtered_aggregate.csv",row.names=F) 

## PLOTS & STATISTICS ------------------

# Read in data

rm(list=ls())
data <- read.csv("allRxData_filtered_aggregate.csv")

# Summary statistics
summary(data)

# View prescribing items / rates before and after intervention
summary(data$ITEMS[data$intervention==0])
summary(data$ITEMS[data$intervention==1])

summary(data$ITEMS.PER.1000[data$intervention==0])
summary(data$ITEMS.PER.1000[data$intervention==1])

## Quasi-Poisson regression model (including term for linear trend, but not seasonality)

# Poisson with the standardised population as an offset
model1 <- glm(ITEMS ~ offset(log(LIST.SIZE)) + intervention + time, family=quasipoisson, data)

summary(model1)
coef(model1) # view coefficients
ci.exp(model1) # view 95% confidence intervals of coefficients
round(ci.lin(model1,Exp=T),3)

# Save model outputs
out1a <- capture.output(summary(model1))
cat("Summary of unadjusted model", out1a, file="Results 1a. Summary of unadjusted model.txt", sep="\n", append=TRUE)
out1b <- capture.output(ci.exp(model1))
cat("CIs of unadjusted model", out1b, file="Results 1b. CIs of unadjusted model.txt", sep="\n", append=TRUE)

# Reform data frame with 0.1 time units to improve plotting
datanew <- data.frame(LIST.SIZE=mean(data$LIST.SIZE),intervention=rep(c(0,1),c(390,210)),
                      time= 1:600/10,month=rep(1:120/10,5))
datanew <- datanew[1:510,]

# Generate predicted values based on the model in order to create a plot
pred1 <- predict(model1,type="response",data)/mean(data$LIST.SIZE)*10^3

# Plot
plot(data$ITEMS.PER.1000,type="n",ylim=c(30,50),xlab="Year",ylab="Items per 1000 registered patients",
     bty="l",xaxt="n")
rect(39,30,51,50,col=grey(0.9),border=F)
points(data$ITEMS.PER.1000,cex=0.7)
axis(1,at=0:5*12,labels=F)
axis(1,at=0:4*12,tick=F,labels=2015:2019)
lines((1:51),pred1,col=2)
title("Prescribing rate, 2015-2019")

# Create and plot the counterfactual
datanew <- data.frame(LIST.SIZE=mean(data$LIST.SIZE),intervention=0,time=1:600/10,
                      month=rep(1:120/10,5))
datanew <- datanew[1:510,]
pred1b <- predict(model1,datanew,type="response")/mean(data$LIST.SIZE)*10^3
lines(datanew$time,pred1b,col=4,lty=2)


# return the data frame to the scenario including the intervention
datanew <- data.frame(LIST.SIZE=mean(data$LIST.SIZE),intervention=rep(c(0,1),c(390,210)), # months before and months after x 10
                      time= 1:600/10,month=rep(1:120/10,5))
datanew <- datanew[1:510,]

###

# Adjusting for seasonality as well as linear trend
# Use harmonic terms specifying the number of sin and cosine pairs to include: 2 pairs with period 12
model2 <- glm(ITEMS ~ offset(log(LIST.SIZE)) + intervention + time +
                harmonic(month,2,12), family=quasipoisson, data)
summary(model2)
coef(model2)
ci.exp(model2)
round(ci.lin(model2,Exp=T),3)

# Save model summary
out2a <- capture.output(summary(model2))
cat("Summary of seasonally adjusted model", out2a, file="Results 2a. Summary of seasonally adjusted model.txt", sep="\n", append=TRUE)
out2b <- capture.output(ci.exp(model2))
cat("CIs of seasonally adjusted model", out2b, file="Results 2b. CIs of seasonally adjusted model.txt", sep="\n", append=TRUE)

# Plot seasonally adjusted model
pred2 <- predict(model2,type="response",datanew)/mean(data$LIST.SIZE)*10^3
plot(data$ITEMS.PER.1000,type="n",ylim=c(30,50),xlab="Year",ylab="Items per 1000 registered patients",
     bty="l",xaxt="n")
rect(39,30,51,50,col=grey(0.9),border=F)
points(data$ITEMS.PER.1000,cex=0.7)
axis(1,at=0:5*12,labels=F)
axis(1,at=0:4*12,tick=F,labels=2015:2019)
lines(1:510/10,pred2,col=2)
title("Prescribing rate, 2015-2019")

# Predict and plot the deseasonalised trend
pred2b <- predict(model2,type="response",transform(datanew,month=6))/
  mean(data$LIST.SIZE)*10^3
lines(1:510/10,pred2b,col=3,lty=2)

# Save plots

Cairo(file="Figure 1. Unadjusted model with counterfactual.png", 
      type="png",
      units="in", 
      width=5, 
      height=4, 
      pointsize=10, 
      dpi=1200)

plot(data$ITEMS.PER.1000,type="n",ylim=c(30,50),xlab="Year",ylab="Items per 1000 registered patients",
     bty="l",xaxt="n")
rect(39,30,51,50,col=grey(0.9),border=F)
points(data$ITEMS.PER.1000,cex=0.3)
axis(1,at=0:5*12,labels=F)
axis(1,at=0:4*12,tick=F,labels=2015:2019)
lines((1:51),pred1,col=2)
title("Prescribing rate, 2015-2019\nUnadjusted model with counterfactual")
datanew <- data.frame(LIST.SIZE=mean(data$LIST.SIZE),intervention=0,time=1:600/10,
                      month=rep(1:120/10,5))
datanew <- datanew[1:510,]
lines(datanew$time,pred1b,col=4,lty=2, lwd=0.6)
dev.off()


Cairo(file="Figure 2. Model adjusted for seasonality.png", 
      type="png",
      units="in", 
      width=5, 
      height=4, 
      pointsize=10, 
      dpi=1200)

plot(data$ITEMS.PER.1000,type="n",ylim=c(30,50),xlab="Year",ylab="Items per 1000 registered patients",
     bty="l",xaxt="n")
rect(39,30,51,50,col=grey(0.9),border=F)
points(data$ITEMS.PER.1000,cex=0.3)
axis(1,at=0:5*12,labels=F)
axis(1,at=0:4*12,tick=F,labels=2015:2019)
lines(1:510/10,pred2,col=2)
title("Prescribing rate, 2015-2019\nSeasonally adjusted model with deseasonalised trend")
lines(1:510/10,pred2b,col=3,lty=2,lwd=0.6)
dev.off()


## REFERENCE ----------------

# Bernal JL et al. Interrupted time series regression for the evaluation of public health interventions: a tutorial. Int J Epidemiol 2017;36:348-355
