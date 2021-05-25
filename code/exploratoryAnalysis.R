setwd("~/Google Drive/UW/Wakefield/WakefieldShared/U5MR/")

# lines represent births
data <- data.frame(read_dta("Kenya2014BirthRecode/KEBR70FL.DTA"))

# get mother's age, sort into bins
birthAge = (data$b3 - data$v011) / 12
interviewAge = (data$v008 - data$v011) / 12
interviewAgeBreaks = seq(14.999, 50, by=5)[-c(5, 7)]
labelsBirthAge = c("5 - 9", "10 - 14", "15 - 19", "20 - 24", "25 - 29", "30 - 39", "40 - 49")
birthAgeBreaks = seq(4.999, 50, by=5)[-c(7, 9)]
labelsInterviewAge = c("15 - 19", "20 - 24", "25 - 29", "30 - 39", "40 - 49")
labelsBirths = c("0 - 4", "5 - 9", "10 - 14", "15 - 19", "20 - 24", "25 - 29", "30 - 34", "35 - 39")
interviewAgeBins = cut(interviewAge, interviewAgeBreaks, labels=labelsInterviewAge)
birthAgeBins = cut(birthAge, birthAgeBreaks, labels=labelsBirthAge)
data$interviewAge = interviewAge
data$birthAge = birthAge
data$interviewAgeFactor = as.character(interviewAgeBins)
data$birthAgeFactor = as.character(birthAgeBins)

dataUrban = data[data$v025==1,]
dataRural = data[data$v025!=1,]

# caveman plots of births through time before survey year
par(mfrow=c(2,2))
hist(2014-data$b2, breaks=seq(-0.5, 38.5, by=1), main="2014 KDHS Births", xlab="Time of Birth (Years Before 2014)")
hist(2014-dataUrban$b2, breaks=seq(-0.5, 38.5, by=1), main="2014 KDHS Urban Births", xlab="Time of Birth (Years Before 2014)")
hist(2014-dataRural$b2, breaks=seq(-0.5, 38.5, by=1), main="2014 KDHS Rural Births", xlab="Time of Birth (Years Before 2014)")

h = hist(birthAge, breaks=seq(-0.5, 49.5, by=1), freq=F)
h$density <- cumsum(h$density) # replace the cell freq.s by cumulative freq.s
plot(h, ylab="Cumulative Density", freq=F, main="Mother's Age at Birth", xlab="Mother's Age at Birth") # plot a cumulative histogram of y
mean(birthAge < 15)

library(ggplot2)
out = hist(2014-data$b2, breaks=seq(-0.5, 40.5, by=5), main="2014 KDHS Births", xlab="Time of Birth (Years Before 2014)")
outUrban = hist(2014-dataUrban$b2, breaks=seq(-0.5, 40.5, by=5), main="2014 KDHS Urban Births", xlab="Time of Birth (Years Before 2014)")
outRural = hist(2014-dataRural$b2, breaks=seq(-0.5, 40.5, by=5), main="2014 KDHS Rural Births", xlab="Time of Birth (Years Before 2014)")
timeLabels = c("0 - 4", "5 - 9", "10 - 14", "15 - 19", "20 - 24", "25 - 29", "30 - 34", "35 - 40")
out$time = timeLabels
outUrban$time = timeLabels
outRural$time = timeLabels
histogramToDataFrame = function(x) {
  class(x) = "list"
  x$xname = NULL
  x$equidist = NULL
  x$breaks = NULL
  x = data.frame(x)
}
out = histogramToDataFrame(out)
outUrban = histogramToDataFrame(outUrban)
outRural = histogramToDataFrame(outRural)
allData = rbind(out, outUrban, outRural)
allData

for(i in 1:length(labelsBirthAge)) {
  thisAgeFactor = labelsBirthAge[i]
  ## calculate histogram for this age range
  
  # first subset by this age range
  thisData = data[data$birthAgeFactor == thisAgeFactor,]
  thisDataUrban = dataUrban[dataUrban$birthAgeFactor == thisAgeFactor,]
  thisDataRural = dataRural[dataRural$birthAgeFactor == thisAgeFactor,]
  
  h = hist(2014-thisData$b2, breaks=seq(-0.5, 40.5, by=5), main="2014 KDHS Births", xlab="Time of Birth (Years Before 2014)")
  hUrban = hist(2014-thisDataUrban$b2, breaks=seq(-0.5, 40.5, by=5), main="2014 KDHS Births", xlab="Time of Birth (Years Before 2014)")
  hRural = hist(2014-thisDataRural$b2, breaks=seq(-0.5, 40.5, by=5), main="2014 KDHS Births", xlab="Time of Birth (Years Before 2014)")
  
  thisFullFrame = histogramToDataFrame(h)
  thisFullFrame$motherAge = thisAgeFactor
  thisFullFrame$yearsSinceBirth = labelsBirths
  thisFullFrameUrban = histogramToDataFrame(hUrban)
  thisFullFrameUrban$motherAge = thisAgeFactor
  thisFullFrameUrban$yearsSinceBirth = labelsBirths
  thisFullFrameRural = histogramToDataFrame(hRural)
  thisFullFrameRural$motherAge = thisAgeFactor
  thisFullFrameRural$yearsSinceBirth = labelsBirths
  if(i == 1) {
    fullFrame = thisFullFrame
    fullFrameUrban = thisFullFrameUrban
    fullFrameRural = thisFullFrameRural
  } else {
    fullFrame = rbind(fullFrame, thisFullFrame)
    fullFrameUrban = rbind(fullFrameUrban, thisFullFrameUrban)
    fullFrameRural = rbind(fullFrameRural, thisFullFrameRural)
  }
}

fullFrame$yearsSinceBirth = factor(fullFrame$yearsSinceBirth, levels=labelsBirths)
fullFrame$motherAge = factor(fullFrame$motherAge, levels=labelsBirthAge)
fullFrameUrban$yearsSinceBirth = factor(fullFrameUrban$yearsSinceBirth, levels=labelsBirths)
fullFrameUrban$motherAge = factor(fullFrameUrban$motherAge, levels=labelsBirthAge)
fullFrameRural$yearsSinceBirth = factor(fullFrameRural$yearsSinceBirth, levels=labelsBirths)
fullFrameRural$motherAge = factor(fullFrameRural$motherAge, levels=labelsBirthAge)

fullFrameAtBirth = fullFrame
fullFrameUrbanAtBirth = fullFrameUrban
fullFrameRuralAtBirth = fullFrameRural

pdf("KDHS2014_birthsVsAgeAtBirthHist.pdf", width=6, height=5)
ggplot(fullFrameAtBirth, aes(fill=motherAge, y=counts, x=yearsSinceBirth)) + 
  geom_bar(position="stack", stat="identity") + ggtitle("2014 KDHS Years Since Birth") + 
  labs(x="Years Since Birth", y="Frequency") + labs(fill = "Mother Age\nat Birth") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        plot.title = element_text(hjust = 0.5))
dev.off()

pdf("KDHS2014_birthsVsAgeAtBirthHistUrban.pdf", width=6, height=5)
ggplot(fullFrameUrbanAtBirth, aes(fill=motherAge, y=counts, x=yearsSinceBirth)) + 
  geom_bar(position="stack", stat="identity") + ggtitle("2014 KDHS Years Since Birth (Urban)") + 
  labs(x="Years Since Birth", y="Frequency") + labs(fill = "Mother Age\nat Birth") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        plot.title = element_text(hjust = 0.5))
dev.off()

pdf("KDHS2014_birthsVsAgeAtBirthHistRural.pdf", width=6, height=5)
ggplot(fullFrameRuralAtBirth, aes(fill=motherAge, y=counts, x=yearsSinceBirth)) + 
  geom_bar(position="stack", stat="identity") + ggtitle("2014 KDHS Years Since Birth (Rural)") + 
  labs(x="Years Since Birth", y="Frequency") + labs(fill = "Mother Age\nat Birth") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        plot.title = element_text(hjust = 0.5))
dev.off()

for(i in 1:length(labelsInterviewAge)) {
  thisAgeFactor = labelsInterviewAge[i]
  ## calculate histogram for this age range
  
  # first subset by this age range
  thisData = data[data$interviewAgeFactor == thisAgeFactor,]
  thisDataUrban = dataUrban[dataUrban$interviewAgeFactor == thisAgeFactor,]
  thisDataRural = dataRural[dataRural$interviewAgeFactor == thisAgeFactor,]
  
  h = hist(2014-thisData$b2, breaks=seq(-0.5, 40.5, by=5), main="2014 KDHS Births", xlab="Time of Birth (Years Before 2014)")
  hUrban = hist(2014-thisDataUrban$b2, breaks=seq(-0.5, 40.5, by=5), main="2014 KDHS Births", xlab="Time of Birth (Years Before 2014)")
  hRural = hist(2014-thisDataRural$b2, breaks=seq(-0.5, 40.5, by=5), main="2014 KDHS Births", xlab="Time of Birth (Years Before 2014)")
  
  thisFullFrame = histogramToDataFrame(h)
  thisFullFrame$motherAge = thisAgeFactor
  thisFullFrame$yearsSinceBirth = labelsBirths
  thisFullFrameUrban = histogramToDataFrame(hUrban)
  thisFullFrameUrban$motherAge = thisAgeFactor
  thisFullFrameUrban$yearsSinceBirth = labelsBirths
  thisFullFrameRural = histogramToDataFrame(hRural)
  thisFullFrameRural$motherAge = thisAgeFactor
  thisFullFrameRural$yearsSinceBirth = labelsBirths
  if(i == 1) {
    fullFrame = thisFullFrame
    fullFrameUrban = thisFullFrameUrban
    fullFrameRural = thisFullFrameRural
  } else {
    fullFrame = rbind(fullFrame, thisFullFrame)
    fullFrameUrban = rbind(fullFrameUrban, thisFullFrameUrban)
    fullFrameRural = rbind(fullFrameRural, thisFullFrameRural)
  }
}

fullFrame$yearsSinceBirth = factor(fullFrame$yearsSinceBirth, levels=labelsBirths)
fullFrame$motherAge = factor(fullFrame$motherAge, levels=labelsInterviewAge)
fullFrameUrban$yearsSinceBirth = factor(fullFrameUrban$yearsSinceBirth, levels=labelsBirths)
fullFrameUrban$motherAge = factor(fullFrameUrban$motherAge, levels=labelsInterviewAge)
fullFrameRural$yearsSinceBirth = factor(fullFrameRural$yearsSinceBirth, levels=labelsBirths)
fullFrameRural$motherAge = factor(fullFrameRural$motherAge, levels=labelsInterviewAge)

fullFrameAtInterview = fullFrame
fullFrameUrbanAtInterview = fullFrameUrban
fullFrameRuralAtInterview = fullFrameRural

pdf("KDHS2014_birthsVsAgeAtInterviewHist.pdf", width=6, height=5)
ggplot(fullFrameAtInterview, aes(fill=motherAge, y=counts, x=yearsSinceBirth)) + 
  geom_bar(position="stack", stat="identity") + ggtitle("2014 KDHS Years Since Birth") + 
  labs(x="Years Since Birth", y="Frequency") + labs(fill = "Mother Age\nat Interview") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        plot.title = element_text(hjust = 0.5))
dev.off()

pdf("KDHS2014_birthsVsAgeAtInterviewHistUrban.pdf", width=6, height=5)
ggplot(fullFrameUrbanAtInterview, aes(fill=motherAge, y=counts, x=yearsSinceBirth)) + 
  geom_bar(position="stack", stat="identity") + ggtitle("2014 KDHS Years Since Birth (Urban)") + 
  labs(x="Years Since Birth", y="Frequency") + labs(fill = "Mother Age\nat Interview") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        plot.title = element_text(hjust = 0.5))
dev.off()

pdf("KDHS2014_birthsVsAgeAtInterviewHistRural.pdf", width=6, height=5)
ggplot(fullFrameRuralAtInterview, aes(fill=motherAge, y=counts, x=yearsSinceBirth)) + 
  geom_bar(position="stack", stat="identity") + ggtitle("2014 KDHS Years Since Birth (Rural)") + 
  labs(x="Years Since Birth", y="Frequency") + labs(fill = "Mother Age\nat Interview") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        plot.title = element_text(hjust = 0.5))
dev.off()




