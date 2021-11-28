###############
## Section 3 ##
###############

library(MASS) 

setwd("C:\\dc_book\\Replication_Files_Element")

ITN.data <- read.csv("ITN_data.csv", header=TRUE)

###############
## Table 3.1 ##
###############

ITN.logit <- glm(useditn ~  malariamessage+education+richest+poorest+pregnant+altitude, family=binomial(link="logit"), data=ITN.data, x=TRUE)
beta <- ITN.logit$coefficients
covmat.beta <- vcov(ITN.logit)
ndraws <- 1000
betadraw <- mvrnorm(ndraws, beta, covmat.beta)

baseline.x <- apply(ITN.logit$x, 2, min)
baseline.x["altitude"] <- median(ITN.data$altitude)
baseline.x["education"] <- median(ITN.data$education)
baseline.x["poorest"] <- 1

hypo.case.A0 <- baseline.x
simprob.A0 <- plogis(betadraw%*%hypo.case.A0)

hypo.case.B0 <- baseline.x
hypo.case.B0["malariamessage"] <- 1
simprob.B0 <- plogis(betadraw%*%hypo.case.B0)

simprob.diff0 <- simprob.B0 - simprob.A0

## Varying Altitude, Rural = 1

hypo.case.A1 <- hypo.case.A0
hypo.case.A1["poorest"] <- 0
hypo.case.A1["richest"] <- 1
simprob.A1 <- plogis(betadraw%*%hypo.case.A1)

hypo.case.B1 <- hypo.case.B0
hypo.case.B1["poorest"] <- 0
hypo.case.B1["richest"] <- 1
simprob.B1 <- plogis(betadraw%*%hypo.case.B1)

simprob.diff1 <- simprob.B1 - simprob.A1

## Difference in Differences

simprob.dd <- simprob.diff1 - simprob.diff0

## summarize results

results.A0 <- append(quantile(simprob.A0, probs=c(0.5, 0.025, 0.975)), sd(simprob.A0))
results.B0 <- append(quantile(simprob.B0, probs=c(0.5, 0.025, 0.975)), sd(simprob.B0))
results.diff0 <- append(quantile(simprob.diff0, probs=c(0.5, 0.025, 0.975)), sd(simprob.diff0))

results.A1 <- append(quantile(simprob.A1, probs=c(0.5, 0.025, 0.975)), sd(simprob.A1))
results.B1 <- append(quantile(simprob.B1, probs=c(0.5, 0.025, 0.975)), sd(simprob.B1))
results.diff1 <- append(quantile(simprob.diff1, probs=c(0.5, 0.025, 0.975)), sd(simprob.diff1))

results.dd <- append(quantile(simprob.dd, probs=c(0.5, 0.025, 0.975)), sd(simprob.dd))

first.differences <- rbind(results.A0, results.B0, results.diff0, results.A1, results.B1, results.diff1, results.dd)
colnames(first.differences) <- c("Median", "2.5%", "97.5%", "(se)")
rownames(first.differences) <- c("No Message:Poor", "Message:Poor", "Difference:Poor", "No Message:Rich", "Message:Rich", "Difference:Rich", "Difference:(Rich-Poor)")
print(round(first.differences, digits=3))

################
## Figure 3.1 ##
################

beta <- ITN.logit$coefficients
covmat.beta <- vcov(ITN.logit)
ndraws <- 1000
betadraw <- mvrnorm(ndraws, beta, covmat.beta)

altitude <- seq(0,2.5,0.01)

mean.x <- colMeans(ITN.logit$x)
simcoef.mean <- matrix(NA, length(altitude), 4)
simprob.all <- matrix(NA, ndraws, length(altitude))

for (i in 1:length(altitude)) {
	mean.x["altitude"] <- altitude[i]
	simprob <- plogis(betadraw%*%mean.x)
	simcoef.mean[i,] <- append(altitude[i], quantile(simprob, probs=c(0.025, 0.5, 0.975)))
	simprob.all[,i] <- simprob
}

ci.bottom <- cbind(simcoef.mean[,1], simcoef.mean[,2])
ci.top <- cbind(simcoef.mean[,1], simcoef.mean[,4])
ci.top <- ci.top[order(ci.top[,1], decreasing=TRUE),]
ci <- rbind(ci.bottom, ci.top)


pdf("fig31a.pdf")
par(mar=c(5.1, 4.1, 4.1, 3.1), mgp=c(2.5, 1, 0), cex.lab=1.5)
plot(simcoef.mean[,3]~simcoef.mean[,1], xlim=c(0.5,2.1), ylim=c(0.5,1), type="n", xlab="Altitude (1000s of meters)", ylab="Probability", xaxt="n", las=1)
axis(1, at=c(0.5, 1, 1.5, 2))
polygon(ci[,1], ci[,2], col="grey75", border=NA)
lines(simcoef.mean[,3]~simcoef.mean[,1], lwd=2, lty=1)
box(lty=1) 
dev.off()

pdf("fig31b.pdf")
par(mar=c(5.1, 4.1, 4.1, 3.1), mgp=c(2.5, 1, 0), cex.lab=1.5)
plot(simcoef.mean[,3]~simcoef.mean[,1], xlim=c(0.5,2.1), ylim=c(0.5,1), type="n", xlab="Altitude (1000s of meters)", ylab="Probability", xaxt="n", las=1)
axis(1, at=c(0.5, 1, 1.5, 2))
for (i in 1:length(altitude)){
	points(simprob.all[i,] ~ altitude, type="l", col=rgb(0,0,0,alpha=0.2), lwd=1)
}
dev.off()


###############
## Table 3.2 ##
###############

## aggregate shares 

one <- 1
Xall <- as.data.frame(na.omit(cbind(one, subset(ITN.data, select=c(malariamessage, education, richest, poorest, pregnant, altitude, sample_weight, useditn)))))
X <- subset(Xall, select=c(one, malariamessage, education, richest, poorest, pregnant, altitude))
y <- Xall$useditn
w <- 1/Xall$sample_weight

X2 <- X
X2$malariamessage <- 1

beta <- ITN.logit$coefficients
covmat.beta <- vcov(ITN.logit)


ndraws <- 1000
betadraw <- mvrnorm(ndraws, beta, covmat.beta)

aggdraw <- matrix(NA,ndraws,1)
aggdraw2 <- matrix(NA,ndraws,1)



for (j in 1:ndraws) { 
  proby1 <- plogis(betadraw[j,]%*%t(X))
  aggshare <- sum(proby1*w)/sum(w)
  aggdraw[j] <- aggshare*100
  proby1.2 <- plogis(betadraw[j,]%*%t(X2))
  aggshare2 <- sum(proby1.2*w)/sum(w)
  aggdraw2[j] <- aggshare2*100
  }

aggdraw.diff <- aggdraw2 - aggdraw

summary.aggshare <- append(quantile(aggdraw, probs=c(0.5, 0.025, 0.975)), sd(aggdraw))   
summary.aggshare2 <- append(quantile(aggdraw2, probs=c(0.5, 0.025, 0.975)), sd(aggdraw2))  
summary.aggshare.diff <- append(quantile(aggdraw.diff, probs=c(0.5, 0.025, 0.975)), sd(aggdraw.diff))  


aggregate.prediction <- rbind(summary.aggshare, summary.aggshare2, summary.aggshare.diff)
colnames(aggregate.prediction) <- c("Median", "2.5%", "97.5%", "(se)")
rownames(aggregate.prediction) <- c("Current Message", "Complete Message", "Difference")
print(round(aggregate.prediction, digits=3))

###############
## Table 3.3 ##
###############

beta <- ITN.logit$coefficients
covmat.beta <- vcov(ITN.logit)
ndraws <- 1000
betadraw <- mvrnorm(ndraws, beta, covmat.beta)

w <- ITN.data$sample_weight  

# AME

predprobs <- plogis(betadraw %*% t(ITN.logit$x))

marginal.ame.matrix <- predprobs * (1 - predprobs) * betadraw[,"altitude"]

marginal.ame <- apply(marginal.ame.matrix, 1, function(x) {sum(x * w)/sum(w)})

elasticity.ame.matrix <- marginal.ame.matrix * t(ITN.logit$x[,"altitude"] / t(predprobs))

elasticity.ame <- apply(elasticity.ame.matrix, 1, function(x) {sum(x * w)/sum(w)})

# MEM

mean.x <- colMeans(ITN.logit$x)
beta <- ITN.logit$coefficients
covmat.beta <- vcov(ITN.logit)

ndraws <- 1000

betadraw <- mvrnorm(ndraws, beta, covmat.beta)

predprob <- plogis(betadraw%*%mean.x)

marginal.mem <- predprob * (1 - predprob) * betadraw[,"altitude"]

elasticity.mem <- marginal.mem * (mean.x["altitude"] / predprob)

## combine results

summary.marginal.ame <- append(quantile(marginal.ame, probs=c(0.5, 0.025, 0.975)), sd(marginal.ame))
summary.marginal.mem <- append(quantile(marginal.mem, probs=c(0.5, 0.025, 0.975)), sd(marginal.mem))
summary.elasticity.ame <- append(quantile(elasticity.ame, probs=c(0.5, 0.025, 0.975)), sd(elasticity.ame))
summary.elasticity.mem <- append(quantile(elasticity.mem, probs=c(0.5, 0.025, 0.975)), sd(elasticity.mem))

marginal.effects.elasticities <- rbind(summary.marginal.ame, summary.marginal.mem, summary.elasticity.ame, summary.elasticity.mem)
colnames(marginal.effects.elasticities) <- c("Median", "2.5%", "97.5%", "(se)")
rownames(marginal.effects.elasticities) <- c("Marginal Effect (AME)", "Marginal Effect (MEM)", "Elasticity (AME)", "Elasticity (MEM)")
print(round(marginal.effects.elasticities, digits=3))


## additional example: MEM by brute force --- note this is actually a first difference

mean.x2 <- mean.x
mean.x2["altitude"] <- mean.x2["altitude"] + 0.0001
probs1 <- plogis(betadraw %*% mean.x)
probs2 <- plogis(betadraw %*% mean.x2)
probs.diff <- (probs2 - probs1) * 10000   # for a 1-unit change
quantile(probs.diff, probs=c(0.5, 0.025, 0.975))


## Elasticity at the mean by brute force
elast.mem <- probs.diff * (mean.x2["altitude"] / probs1)
quantile(elast.mem, probs=c(0.5, 0.025, 0.975))

## coefficient ratios

malariamessage.altitude.ratio <- betadraw[,"malariamessage"]/betadraw[,"altitude"]
ratio.results <- quantile(malariamessage.altitude.ratio, probs=c(0.5, 0.025, 0.975))
print(ratio.results)


###############
## Table 3.4 ##
###############

N <- length(ITN.logit$fitted.values)

beta <- ITN.logit$coefficients
covmat.beta <- vcov(ITN.logit)
ndraws <- 1000

# calculating confidence intervals and standard error via simulation

betadraw <- mvrnorm(ndraws, beta, covmat.beta)
exp.betadraw <- exp(betadraw)
quantiles <- apply(exp.betadraw, 2, function(x) {quantile(x, probs=c(0.5, 0.025, 0.975))})
sds <- apply(exp.betadraw, 2, sd)
odds.ratios <- cbind(t(quantiles), sds)
colnames(odds.ratios) <- c("Median", "2.5%", "97.5%", "(se)")
rownames(odds.ratios) <- names(beta) 
print(round(odds.ratios, digits=3))
print(N)

# alternative calculation exponentiating endpoints of coefficient CI and using delta method for standard errors

beta.OR <- exp(ITN.logit$coefficients)
se.OR <- exp(ITN.logit$coefficients) * sqrt(diag(vcov(ITN.logit)))
CI.OR <- exp(confint.default(ITN.logit, level=0.95))

ITN.logit.results.OR <- cbind(beta.OR, CI.OR, se.OR)
N <- length(ITN.logit$fitted.values)
colnames(ITN.logit.results.OR) <- c("Coeff", "2.5%", "97.5%", "(se)")
print((round(ITN.logit.results.OR, digits=3)))

