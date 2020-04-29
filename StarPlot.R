# Load in data
load("Recoded Data.RData")

## Create Data
# Formant mean data
mean.names <- c('F1.mean','F2.mean','F3.mean','gender','Recoded.labels','Target.keywords','Corpus','Preceding.POA','Following.POA','Speaker')
mean.data <- data[,mean.names]
colnames(mean.data) <- c('F1','F2','F3','Gender','Vowel','Keyword','Generation','Preceding','Following','Speaker')
# Normalised results
order.data <- mean.data[order(mean.data$Speaker),]
# Subset data by speaker
levels(order.data$Speaker)
norm.formant <- c(0,0)
for(i in 1: length(levels(order.data$Speaker)))  {
  speaker <- levels(order.data$Speaker)[i]
  speaker.sub <- subset(order.data,Speaker== speaker)
  F1.Lobanov <- (speaker.sub$F1-mean(speaker.sub$F1))/sd(speaker.sub$F1)
  F2.Lobanov <- (speaker.sub$F2-mean(speaker.sub$F2))/sd(speaker.sub$F2)
  #F3.Lobanov <- (speaker.sub$F3-mean(speaker.sub$F3))/sd(speaker.sub$F3)
  lob.results <- cbind(F1.Lobanov,F2.Lobanov)
  norm.formant <- rbind(norm.formant,lob.results) }

norm.formant <- norm.formant[-1,]
# order.data and not mean.data

mean.data <- cbind(order.data,norm.formant)


do.plot <- function(d0,variable="Vowel",all.lims=FALSE, lobanov=FALSE, palette=rainbow) {
  names <- c("F1", "F2")
  if (lobanov)
    names <- paste0(names, ".Lobanov")
  plot(formula(paste0(names[2],"~",names[1])), data=if (all.lims) mean.data else d0, type="n")
  n.speakers <- length(levels(d0$Speaker))
  n.secondary <- length(levels(mean.data[[variable]]))
  cols <- palette(n.secondary)
  cols.transparent <- palette(n.secondary, alpha=0.25)
  for (j in sample(1:n.speakers)) {
    d <- subset(d0, Speaker==levels(d0$Speaker)[j])
    for (i in sample(1:n.secondary)) {
      dd <- eval(parse(text=paste0('subset(d, ',variable,'=="',levels(d[[variable]])[i],'")')))
      if (nrow(dd)==0)
        next
      points(mean(dd[[names[1]]]), mean(dd[[names[2]]]), col=cols[i], pch=16, cex=2)
      points(dd[[names[1]]], dd[[names[2]]], col=cols.transparent[i], pch=16, cex=0.5)
      segments(dd[[names[1]]], dd[[names[2]]], mean(dd[[names[1]]]), mean(dd[[names[2]]]), col=cols.transparent[i])
    }
    
  }
}

par(mfrow=1:2)
do.plot(subset(mean.data, Generation=="00-O"), all.lims=TRUE, lobanov=TRUE) 
n.vowels <- length(levels(mean.data$Vowel))
legend("topright", col=rainbow(n.vowels), pch=15, levels(mean.data$Vowel)) 
title("Old speaker in the 00s") 
do.plot(subset(mean.data, Generation=="00-Y"), all.lims=TRUE, lobanov=TRUE) 
title("Young speakers in the 00s")

par(mfrow=1:2)
# Example: Preceding context matters, but not following do.plot(subset(mean.data, Vowel=="uw"), variable="Preceding") legend("topright", col=rainbow(length(levels(mean.data$Preceding))), pch=16, levels(mean.data$Preceding)) title("Vowel uw: Effect of Preceding Context") do.plot(subset(mean.data, Vowel=="uw"), variable="Following") legend("topright", col=rainbow(length(levels(mean.data$Following))), pch=16, levels(mean.data$Following)) title("Vowel uw: Effect of Following Context")


par(mfrow=1:2)
# Generation (for o seems to have an effect) palette <- function(...,alpha=1) heat.colors(4, alpha=alpha)[c(3,1,4,2)] do.plot(subset(mean.data, Vowel=="o"), variable="Generation", palette=palette) legend("topright", col=heat.colors(4), pch=15, levels(mean.data$Generation)[c(3,1,4,2)])
title("Vowel u: Effect of Generation")

# but not for a
palette <- function(...,alpha=1) heat.colors(4, alpha=alpha)[c(3,1,4,2)] 
do.plot(subset(mean.data, Vowel=="a"), variable="Generation", palette=palette) 
legend("topright", col=heat.colors(4), pch=15, levels(mean.data$Generation)[c(3,1,4,2)])
title("Vowel a: Effect of Generation")
