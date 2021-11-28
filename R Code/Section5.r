###############
## Section 5 ##
###############

#install.packages("mlogit")
library(mlogit)
library(MASS)

setwd("C:\\dc_book\\Replication_Files_Element")

bes.data.raw <- read.csv("BES2015.csv", header=TRUE)

bes.data <- mlogit.data(bes.data.raw, choice="vote", shape="long", chid.var="respondent_id", alt.var="party_id")

###############
## Table 5.1 ##
###############

bes2015.mnl <- mlogit(vote ~ taxspend | approve_EU + NHS_improved + income + female + age, data=bes.data, reflevel=4)
bes2015.mnl.results <- cbind(bes2015.mnl$coefficients, sqrt(diag(vcov(bes2015.mnl))), confint.default(bes2015.mnl, level=0.95))
N <- length(unique(bes.data$respondent_id))
A <- nrow(bes.data)
colnames(bes2015.mnl.results) <- c("Coeff", "(se)", "2.5%", "97.5%")
print((round(bes2015.mnl.results, digits=3)))
print(N)
print(A)

###############
## Table 5.2 ##
###############

beta <- bes2015.mnl$coefficients
covmat.beta <- vcov(bes2015.mnl)
ndraws <- 1000
betadraw <- mvrnorm(ndraws, beta, covmat.beta)

baseline.x <- with(bes.data, data.frame(
constant = rbind(0, diag(3)),
taxspend = tapply(taxspend, idx(bes2015.mnl, 2), median), 
approve_EU = median(approve_EU) * rbind(0, diag(3)), 
NHS_improved = median(NHS_improved) * rbind(0, diag(3)),
income=median(income) * rbind(0, diag(3)),
female=median(female) * rbind(0, diag(3)),
age=median(age) * rbind(0, diag(3))
))

baseline.x <- baseline.x[order(rownames(baseline.x)),]

approveEUcol <- grep("approve_EU", colnames(baseline.x))
diff.x <- baseline.x
diff.x[,approveEUcol] <- rbind(diag(3), 0) * 2

exb.baseline <- exp(betadraw %*% t(baseline.x))
exb.diff <- exp(betadraw %*% t(diff.x))

probs.baseline <- exb.baseline/rowSums(exb.baseline)
probs.diff <- exb.diff/rowSums(exb.diff)

simprob.diff <- probs.diff - probs.baseline 

all.results <- cbind(probs.baseline, probs.diff, simprob.diff)

quantiles <- apply(all.results, 2, function(x) {quantile(x, probs=c(0.5, 0.025, 0.975))})
sds <- apply(all.results, 2, sd)
first.differences <- cbind(t(quantiles), sds)
colnames(first.differences) <- c("Median", "2.5%", "97.5%", "(se)")
rownames(first.differences) <- c("Neutral:Conservative", "Neutral:Labour", "Neutral:LibDem", "Neutral:UKIP", "Disapprove:Conservative", "Disapprove:Labour", "Disapprove:LibDem", "Disapprove:UKIP", "Diff:Conservative", "Diff:Labour", "Diff:LibDem", "Diff:UKIP")
print(round(first.differences, digits=3))



#################
## Figure 5.1 ##
#################

beta <- bes2015.mnl$coefficients
covmat.beta <- vcov(bes2015.mnl)
ndraws <- 1000
betadraw <- mvrnorm(ndraws, beta, covmat.beta)

baseline.x <- with(bes.data, data.frame(
constant = rbind(0, diag(3)),
taxspend = tapply(taxspend, idx(bes2015.mnl, 2), median), 
approve_EU = median(approve_EU) * rbind(0, diag(3)), 
NHS_improved = median(NHS_improved) * rbind(0, diag(3)),
income = median(income) * rbind(0, diag(3)),
female = median(female) * rbind(0, diag(3)),
age = median(age) * rbind(0, diag(3))
))

## the reference level is always first in idx 
# reorder rows to avoid confusion
baseline.x <- baseline.x[order(rownames(baseline.x)),]

agecols <- grep("age", colnames(baseline.x))
age.range <- c(10:100)

mlogit.probs1 <- matrix(NA, ndraws, length(age.range))
mlogit.probs2 <- matrix(NA, ndraws, length(age.range))
mlogit.probs3 <- matrix(NA, ndraws, length(age.range))
mlogit.probs4 <- matrix(NA, ndraws, length(age.range))

mlogit.mean.probs <- matrix(NA, 4, length(age.range))


for (i in 1:length(age.range)) {

	baseline.x[,agecols] <- rbind(diag(3), 0) * age.range[i]
	exb <- exp(betadraw %*% t(baseline.x))	
	
	mlogit.probs1[,i] <- exb[,1]/rowSums(exb)
	mlogit.probs2[,i] <- exb[,2]/rowSums(exb)
	mlogit.probs3[,i] <- exb[,3]/rowSums(exb)
	mlogit.probs4[,i] <- exb[,4]/rowSums(exb)
	
	mprobs.mean <- exp(beta%*%t(baseline.x))
	mlogit.mean.probs[1,i] <- mprobs.mean[1]/sum(mprobs.mean)
	mlogit.mean.probs[2,i] <- mprobs.mean[2]/sum(mprobs.mean)
	mlogit.mean.probs[3,i] <- mprobs.mean[3]/sum(mprobs.mean)
	mlogit.mean.probs[4,i] <- mprobs.mean[4]/sum(mprobs.mean)
	
}


# Conservative, Labour, Liberal Democrat, UKIP

colors <- c("#47E3FF80", "#FF634780", "#B1D87780","#6347FF80")


pdf("fig51.pdf")
par(mar=c(5.1, 4.1, 4.1, 5.6), mgp=c(2.5, 1, 0), cex.lab=1.5)
plot(mlogit.mean.probs[1,]~age.range, xlim=c(18,90), ylim=c(0,0.8), type="n", xlab="Age", ylab="Probability", xaxt="n", las=1)
axis(1, at=c(18, 30, 40, 50, 60, 70, 80, 90))
for (i in 1:length(age.range)){
	points(mlogit.probs1[i,] ~ age.range, type="l", col=colors[1], lwd=1)
	points(mlogit.probs2[i,] ~ age.range, type="l", col=colors[2], lwd=1)
	points(mlogit.probs3[i,] ~ age.range, type="l", col=colors[3], lwd=1)
	points(mlogit.probs4[i,] ~ age.range, type="l", col=colors[4], lwd=1)
}
lines(mlogit.mean.probs[1,]~age.range, col=rgb(0,0,0,0.5), lwd=2)
lines(mlogit.mean.probs[2,]~age.range, col=rgb(0,0,0,0.5), lwd=2)
lines(mlogit.mean.probs[3,]~age.range, col=rgb(0,0,0,0.5), lwd=2)
lines(mlogit.mean.probs[4,]~age.range, col=rgb(0,0,0,0.5), lwd=2)
mtext("Conservative", side=4, at=mlogit.mean.probs[1,84], cex=1.1, las=1, line=0.3)
mtext("Labour", side=4, at=mlogit.mean.probs[2,84]+0.02, cex=1.1, las=1, line=0.3)
mtext("Liberal", side=4, at=mlogit.mean.probs[3,84]-0.01, cex=1.1, las=1, line=0.3)
mtext("Democrat", side=4, at=mlogit.mean.probs[3,84]-0.04, cex=1.1, las=1, line=0.3)
mtext("UKIP", side=4, at=mlogit.mean.probs[4,84], cex=1.1, las=1, line=0.3)
dev.off()


#################
## Figure 5.2 ##
#################

beta <- bes2015.mnl$coefficients
covmat.beta <- vcov(bes2015.mnl)
ndraws <- 1000
betadraw <- mvrnorm(ndraws, beta, covmat.beta)


baseline.x <- with(bes.data, data.frame(
constant = rbind(0, diag(3)),
taxspend = tapply(taxspend, idx(bes2015.mnl, 2), median), 
approve_EU = median(approve_EU) * rbind(0, diag(3)), 
NHS_improved = median(NHS_improved) * rbind(0, diag(3)),
income = median(income) * rbind(0, diag(3)),
female = median(female) * rbind(0, diag(3)),
age = median(age) * rbind(0, diag(3))
))

baseline.x <- baseline.x[order(rownames(baseline.x)),]

taxspendcol <- grep("taxspend", colnames(baseline.x))
taxspend.range <- seq(-5,15,0.1) 
mlogit.probs1 <- matrix(NA, ndraws, length(taxspend.range))
mlogit.probs2 <- matrix(NA, ndraws, length(taxspend.range))
mlogit.probs3 <- matrix(NA, ndraws, length(taxspend.range))
mlogit.probs4 <- matrix(NA, ndraws, length(taxspend.range))

mlogit.mean.probs <- matrix(NA, 4, length(taxspend.range))

for (i in 1:length(taxspend.range)) {

	baseline.x[rownames(baseline.x)=="4", taxspendcol] <- taxspend.range[i] 
	exb <- exp(betadraw %*% t(baseline.x))
	mlogit.probs1[,i] <- exb[,1]/rowSums(exb)
	mlogit.probs2[,i] <- exb[,2]/rowSums(exb)
	mlogit.probs3[,i] <- exb[,3]/rowSums(exb)
	mlogit.probs4[,i] <- exb[,4]/rowSums(exb)
	
	mprobs.mean <- exp(beta%*%t(baseline.x))
	mlogit.mean.probs[1,i] <- mprobs.mean[1]/sum(mprobs.mean)
	mlogit.mean.probs[2,i] <- mprobs.mean[2]/sum(mprobs.mean)
	mlogit.mean.probs[3,i] <- mprobs.mean[3]/sum(mprobs.mean)
	mlogit.mean.probs[4,i] <- mprobs.mean[4]/sum(mprobs.mean)
}

# Conservative, Labour, Liberal Democrat, UKIP

colors <- c("#47E3FF80", "#FF634780", "#B1D87780","#6347FF80")

pdf("fig52.pdf")
par(mar=c(5.1, 4.1, 4.1, 5.6), mgp=c(2.5, 1, 0), cex.lab=1.5)
plot(mlogit.mean.probs[1,]~taxspend.range, xlim=c(1,10), ylim=c(0,0.6), type="n", xlab="Tax/Spend Scale", ylab="Probability", xaxt="n", las=1)
axis(1, at=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
for (i in 1:length(taxspend.range)){
	points(mlogit.probs1[i,] ~ taxspend.range, type="l", col=colors[1], lwd=1)
	points(mlogit.probs2[i,] ~ taxspend.range, type="l", col=colors[2], lwd=1)
	points(mlogit.probs3[i,] ~ taxspend.range, type="l", col=colors[3], lwd=1)
	points(mlogit.probs4[i,] ~ taxspend.range, type="l", col=colors[4], lwd=1)
}
lines(mlogit.mean.probs[1,]~taxspend.range, col=rgb(0,0,0,0.5), lwd=2)
lines(mlogit.mean.probs[2,]~taxspend.range, col=rgb(0,0,0,0.5), lwd=2)
lines(mlogit.mean.probs[3,]~taxspend.range, col=rgb(0,0,0,0.5), lwd=2)
lines(mlogit.mean.probs[4,]~taxspend.range, col=rgb(0,0,0,0.5), lwd=2)
mtext("Conservative", side=4, at=mlogit.mean.probs[1,156], cex=1.1, las=1, line=0.3)
mtext("Labour", side=4, at=mlogit.mean.probs[2,156], cex=1.1, las=1, line=0.3)
mtext("Liberal", side=4, at=mlogit.mean.probs[3,156]+0.01, cex=1.1, las=1, line=0.3)
mtext("Democrat", side=4, at=mlogit.mean.probs[3,156]-0.01, cex=1.1, las=1, line=0.3)
mtext("UKIP", side=4, at=mlogit.mean.probs[4,156], cex=1.1, las=1, line=0.3)
dev.off()


###############
## Table 5.3 ##
###############

beta <- bes2015.mnl$coefficients
covmat.beta <- vcov(bes2015.mnl)
ndraws <- 1000
betadraw <- mvrnorm(ndraws, beta, covmat.beta)

## observed

observed.x <- with(bes.data, data.frame(
cons1 = as.numeric(idx[2]==1),
cons2 = as.numeric(idx[2]==2),
cons3 = as.numeric(idx[2]==3),
taxspend = taxspend, 
approve_EU1 = approve_EU * as.numeric(idx[2]==1), 
approve_EU2 = approve_EU * as.numeric(idx[2]==2), 
approve_EU3 = approve_EU * as.numeric(idx[2]==3), 
NHS_improved1 = NHS_improved * as.numeric(idx[2]==1),
NHS_improved2 = NHS_improved * as.numeric(idx[2]==2),
NHS_improved3 = NHS_improved * as.numeric(idx[2]==3),
income1 = income * as.numeric(idx[2]==1),
income2 = income * as.numeric(idx[2]==2),
income3 = income * as.numeric(idx[2]==3),
female1 = female * as.numeric(idx[2]==1),
female2 = female * as.numeric(idx[2]==2),
female3 = female * as.numeric(idx[2]==3),
age1 = age * as.numeric(idx[2]==1),
age2 = age * as.numeric(idx[2]==2),
age3 = age * as.numeric(idx[2]==3),
id = idx,
weights = wt_combined_main
))

exb.party1.o <- exp(betadraw %*% t(subset(observed.x, id.alt==1, select=c(1:19))))
exb.party2.o <- exp(betadraw %*% t(subset(observed.x, id.alt==2, select=c(1:19))))
exb.party3.o <- exp(betadraw %*% t(subset(observed.x, id.alt==3, select=c(1:19))))
exb.party4.o <- exp(betadraw %*% t(subset(observed.x, id.alt==4, select=c(1:19))))
exb.all.o <- exb.party1.o + exb.party2.o + exb.party3.o + exb.party4.o

probs1.o <- exb.party1.o/exb.all.o
probs2.o <- exb.party2.o/exb.all.o
probs3.o <- exb.party3.o/exb.all.o
probs4.o <- exb.party4.o/exb.all.o

id.weights <- observed.x$weights[observed.x$id.alt==1]

shares.party1.o <- colSums(apply(probs1.o, 1, function(x){x * id.weights}))/sum(id.weights) * 100
shares.party2.o <- colSums(apply(probs2.o, 1, function(x){x * id.weights}))/sum(id.weights) * 100
shares.party3.o <- colSums(apply(probs3.o, 1, function(x){x * id.weights}))/sum(id.weights) * 100
shares.party4.o <- colSums(apply(probs4.o, 1, function(x){x * id.weights}))/sum(id.weights) * 100

## counterfactual

counterfactual.x <- observed.x
counterfactual.x$rand.draw <- rep(runif(length(unique(counterfactual.x$id.chid))), each=4)
counterfactual.x$approve_EU1[counterfactual.x$approve_EU1 > 1 & counterfactual.x$rand.draw <= 0.25] <- counterfactual.x$approve_EU1[counterfactual.x$approve_EU1 > 1 & counterfactual.x$rand.draw <= 0.25] - 1
counterfactual.x$approve_EU2[counterfactual.x$approve_EU2 > 1 & counterfactual.x$rand.draw <= 0.25] <- counterfactual.x$approve_EU2[counterfactual.x$approve_EU2 > 1 & counterfactual.x$rand.draw <= 0.25] - 1
counterfactual.x$approve_EU3[counterfactual.x$approve_EU3 > 1 & counterfactual.x$rand.draw <= 0.25] <- counterfactual.x$approve_EU3[counterfactual.x$approve_EU3 > 1 & counterfactual.x$rand.draw <= 0.25] - 1

exb.party1.c <- exp(betadraw %*% t(subset(counterfactual.x, id.alt==1, select=c(1:19))))
exb.party2.c <- exp(betadraw %*% t(subset(counterfactual.x, id.alt==2, select=c(1:19))))
exb.party3.c <- exp(betadraw %*% t(subset(counterfactual.x, id.alt==3, select=c(1:19))))
exb.party4.c <- exp(betadraw %*% t(subset(counterfactual.x, id.alt==4, select=c(1:19))))
exb.all.c <- exb.party1.c + exb.party2.c + exb.party3.c + exb.party4.c

probs1.c <- exb.party1.c/exb.all.c
probs2.c <- exb.party2.c/exb.all.c
probs3.c <- exb.party3.c/exb.all.c
probs4.c <- exb.party4.c/exb.all.c

shares.party1.c <- colSums(apply(probs1.c, 1, function(x){x * id.weights}))/sum(id.weights) * 100
shares.party2.c <- colSums(apply(probs2.c, 1, function(x){x * id.weights}))/sum(id.weights) * 100
shares.party3.c <- colSums(apply(probs3.c, 1, function(x){x * id.weights}))/sum(id.weights) * 100
shares.party4.c <- colSums(apply(probs4.c, 1, function(x){x * id.weights}))/sum(id.weights) * 100

## difference

probs1.diff <- probs1.c - probs1.o
probs2.diff <- probs2.c - probs2.o
probs3.diff <- probs3.c - probs3.o
probs4.diff <- probs4.c - probs4.o

shares.party1.diff <- colSums(apply(probs1.diff, 1, function(x){x * id.weights}))/sum(id.weights) * 100
shares.party2.diff <- colSums(apply(probs2.diff, 1, function(x){x * id.weights}))/sum(id.weights) * 100
shares.party3.diff <- colSums(apply(probs3.diff, 1, function(x){x * id.weights}))/sum(id.weights) * 100
shares.party4.diff <- colSums(apply(probs4.diff, 1, function(x){x * id.weights}))/sum(id.weights) * 100

all.results <- cbind(shares.party1.o, shares.party2.o, shares.party3.o, shares.party4.o, shares.party1.c, shares.party2.c, shares.party3.c, shares.party4.c, shares.party1.diff, shares.party2.diff, shares.party3.diff, shares.party4.diff)

quantiles <- apply(all.results, 2, function(x) {quantile(x, probs=c(0.5, 0.025, 0.975))})
sds <- apply(all.results, 2, sd)
first.differences <- cbind(t(quantiles), sds)
colnames(first.differences) <- c("Median", "2.5%", "97.5%", "(se)")
rownames(first.differences) <- c("Neutral:Conservative", "Neutral:Labour", "Neutral:LibDem", "Neutral:UKIP", "Disapprove:Conservative", "Disapprove:Labour", "Disapprove:LibDem", "Disapprove:UKIP", "Diff:Conservative", "Diff:Labour", "Diff:LibDem", "Diff:UKIP")
print(round(first.differences, digits=3))

###############
## Table 5.4 ##
###############

beta <- bes2015.mnl$coefficients
covmat.beta <- vcov(bes2015.mnl)
ndraws <- 1000
betadraw <- mvrnorm(ndraws, beta, covmat.beta)

approveEU.coeffs <- cbind(betadraw[,grep("approve_EU", colnames(betadraw))], 0)

## AME ##

observed.x <- with(bes.data, data.frame(
cons1 = as.numeric(idx[2]==1),
cons2 = as.numeric(idx[2]==2),
cons3 = as.numeric(idx[2]==3),
taxspend = taxspend, 
approve_EU1 = approve_EU * as.numeric(idx[2]==1), 
approve_EU2 = approve_EU * as.numeric(idx[2]==2), 
approve_EU3 = approve_EU * as.numeric(idx[2]==3), 
NHS_improved1 = NHS_improved * as.numeric(idx[2]==1),
NHS_improved2 = NHS_improved * as.numeric(idx[2]==2),
NHS_improved3 = NHS_improved * as.numeric(idx[2]==3),
income1 = income * as.numeric(idx[2]==1),
income2 = income * as.numeric(idx[2]==2),
income3 = income * as.numeric(idx[2]==3),
female1 = female * as.numeric(idx[2]==1),
female2 = female * as.numeric(idx[2]==2),
female3 = female * as.numeric(idx[2]==3),
age1 = age * as.numeric(idx[2]==1),
age2 = age * as.numeric(idx[2]==2),
age3 = age * as.numeric(idx[2]==3),
id = idx,
weights = wt_combined_main
))

approveEU <- observed.x[observed.x$id.alt==1, grep("approve_EU1", colnames(observed.x))]

exb.party1 <- exp(betadraw %*% t(subset(observed.x, id.alt==1, select=c(1:19))))
exb.party2 <- exp(betadraw %*% t(subset(observed.x, id.alt==2, select=c(1:19))))
exb.party3 <- exp(betadraw %*% t(subset(observed.x, id.alt==3, select=c(1:19))))
exb.party4 <- exp(betadraw %*% t(subset(observed.x, id.alt==4, select=c(1:19))))
exb.all <- exb.party1 + exb.party2 + exb.party3 + exb.party4

probs1 <- exb.party1/exb.all
probs2 <- exb.party2/exb.all
probs3 <- exb.party3/exb.all
probs4 <- exb.party4/exb.all

id.weights <- observed.x$weights[observed.x$id.alt==1]/observed.x$weights[observed.x$id.alt==1]

## marginal effects

me.approveEU.1 <- probs1 * (approveEU.coeffs[,1] - ((probs1 * approveEU.coeffs[,1]) + (probs2 * approveEU.coeffs[,2]) + (probs3 * approveEU.coeffs[,3]) + (probs4 * approveEU.coeffs[,4])))
me.approveEU.2 <- probs2 * (approveEU.coeffs[,2] - ((probs1 * approveEU.coeffs[,1]) + (probs2 * approveEU.coeffs[,2]) + (probs3 * approveEU.coeffs[,3]) + (probs4 * approveEU.coeffs[,4])))
me.approveEU.3 <- probs3 * (approveEU.coeffs[,3] - ((probs1 * approveEU.coeffs[,1]) + (probs2 * approveEU.coeffs[,2]) + (probs3 * approveEU.coeffs[,3]) + (probs4 * approveEU.coeffs[,4])))
me.approveEU.4 <- probs4 * (approveEU.coeffs[,4] - ((probs1 * approveEU.coeffs[,1]) + (probs2 * approveEU.coeffs[,2]) + (probs3 * approveEU.coeffs[,3]) + (probs4 * approveEU.coeffs[,4])))

me.AME.approveEU.1 <- colSums(apply(me.approveEU.1, 1, function(x){x * id.weights}))/sum(id.weights)
me.AME.approveEU.2 <- colSums(apply(me.approveEU.2, 1, function(x){x * id.weights}))/sum(id.weights)
me.AME.approveEU.3 <- colSums(apply(me.approveEU.3, 1, function(x){x * id.weights}))/sum(id.weights)
me.AME.approveEU.4 <- colSums(apply(me.approveEU.4, 1, function(x){x * id.weights}))/sum(id.weights)


## elasticities

e.approveEU.1 <- approveEU * (approveEU.coeffs[,1] - ((probs1 * approveEU.coeffs[,1]) + (probs2 * approveEU.coeffs[,2]) + (probs3 * approveEU.coeffs[,3]) + (probs4 * approveEU.coeffs[,4])))
e.approveEU.2 <- approveEU * (approveEU.coeffs[,2] - ((probs1 * approveEU.coeffs[,1]) + (probs2 * approveEU.coeffs[,2]) + (probs3 * approveEU.coeffs[,3]) + (probs4 * approveEU.coeffs[,4])))
e.approveEU.3 <- approveEU * (approveEU.coeffs[,3] - ((probs1 * approveEU.coeffs[,1]) + (probs2 * approveEU.coeffs[,2]) + (probs3 * approveEU.coeffs[,3]) + (probs4 * approveEU.coeffs[,4])))
e.approveEU.4 <- approveEU * (approveEU.coeffs[,4] - ((probs1 * approveEU.coeffs[,1]) + (probs2 * approveEU.coeffs[,2]) + (probs3 * approveEU.coeffs[,3]) + (probs4 * approveEU.coeffs[,4])))

e.AME.approveEU.1 <- colSums(apply(e.approveEU.1, 1, function(x){x * id.weights}))/sum(id.weights)
e.AME.approveEU.2 <- colSums(apply(e.approveEU.2, 1, function(x){x * id.weights}))/sum(id.weights)
e.AME.approveEU.3 <- colSums(apply(e.approveEU.3, 1, function(x){x * id.weights}))/sum(id.weights)
e.AME.approveEU.4 <- colSums(apply(e.approveEU.4, 1, function(x){x * id.weights}))/sum(id.weights)

## MEM ##

mean.x <- with(bes.data, data.frame(
constant = rbind(0, diag(3)),
taxspend = tapply(taxspend, idx(bes2015.mnl, 2), mean), 
approve_EU = mean(approve_EU) * rbind(0, diag(3)), 
NHS_improved = mean(NHS_improved) * rbind(0, diag(3)),
income=mean(income) * rbind(0, diag(3)),
female=mean(female) * rbind(0, diag(3)),
age=mean(age) * rbind(0, diag(3))
))

mean.x <- mean.x[order(rownames(mean.x)),]
mean.approveEU <- mean(bes.data$approve_EU)

exb.mean <- exp(betadraw %*% t(mean.x))

probs.mean <- exb.mean/rowSums(exb.mean)

## marginal effects

me.MEM.approveEU.1 <- probs.mean[,1] * (approveEU.coeffs[,1] - ((probs.mean[,1] * approveEU.coeffs[,1]) + (probs.mean[,2] * approveEU.coeffs[,2]) + (probs.mean[,3] * approveEU.coeffs[,3]) + (probs.mean[,4] * approveEU.coeffs[,4])))
me.MEM.approveEU.2 <- probs.mean[,2] * (approveEU.coeffs[,2] - ((probs.mean[,1] * approveEU.coeffs[,1]) + (probs.mean[,2] * approveEU.coeffs[,2]) + (probs.mean[,3] * approveEU.coeffs[,3]) + (probs.mean[,4] * approveEU.coeffs[,4])))
me.MEM.approveEU.3 <- probs.mean[,3] * (approveEU.coeffs[,3] - ((probs.mean[,1] * approveEU.coeffs[,1]) + (probs.mean[,2] * approveEU.coeffs[,2]) + (probs.mean[,3] * approveEU.coeffs[,3]) + (probs.mean[,4] * approveEU.coeffs[,4])))
me.MEM.approveEU.4 <- probs.mean[,4] * (approveEU.coeffs[,4] - ((probs.mean[,1] * approveEU.coeffs[,1]) + (probs.mean[,2] * approveEU.coeffs[,2]) + (probs.mean[,3] * approveEU.coeffs[,3]) + (probs.mean[,4] * approveEU.coeffs[,4])))

## elasticities 

e.MEM.approveEU.1 <- mean.approveEU * (approveEU.coeffs[,1] - ((probs.mean[,1] * approveEU.coeffs[,1]) + (probs.mean[,2] * approveEU.coeffs[,2]) + (probs.mean[,3] * approveEU.coeffs[,3]) + (probs.mean[,4] * approveEU.coeffs[,4])))
e.MEM.approveEU.2 <- mean.approveEU * (approveEU.coeffs[,2] - ((probs.mean[,1] * approveEU.coeffs[,1]) + (probs.mean[,2] * approveEU.coeffs[,2]) + (probs.mean[,3] * approveEU.coeffs[,3]) + (probs.mean[,4] * approveEU.coeffs[,4])))
e.MEM.approveEU.3 <- mean.approveEU * (approveEU.coeffs[,3] - ((probs.mean[,1] * approveEU.coeffs[,1]) + (probs.mean[,2] * approveEU.coeffs[,2]) + (probs.mean[,3] * approveEU.coeffs[,3]) + (probs.mean[,4] * approveEU.coeffs[,4])))
e.MEM.approveEU.4 <- mean.approveEU * (approveEU.coeffs[,4] - ((probs.mean[,1] * approveEU.coeffs[,1]) + (probs.mean[,2] * approveEU.coeffs[,2]) + (probs.mean[,3] * approveEU.coeffs[,3]) + (probs.mean[,4] * approveEU.coeffs[,4])))

all.results <- cbind(me.AME.approveEU.1, me.AME.approveEU.2, me.AME.approveEU.3, me.AME.approveEU.4, me.MEM.approveEU.1, me.MEM.approveEU.2, me.MEM.approveEU.3, me.MEM.approveEU.4, e.AME.approveEU.1, e.AME.approveEU.2, e.AME.approveEU.3, e.AME.approveEU.4, e.MEM.approveEU.1, e.MEM.approveEU.2, e.MEM.approveEU.3, e.MEM.approveEU.4)

quantiles <- apply(all.results, 2, function(x) {quantile(x, probs=c(0.5, 0.025, 0.975))})
sds <- apply(all.results, 2, sd)
marginal.effects <- cbind(t(quantiles), sds)
colnames(marginal.effects) <- c("Median", "2.5%", "97.5%", "(se)")
rownames(marginal.effects) <- c("Conservative:AME", "Labour:AME", "LibDem:AME", "UKIP:AME", "Conservative:MEM", "Labour:MEM", "LibDem:MEM", "UKIP:MEM", "Conservative: Average Elasticity", "Labour: Average Elasticity", "LibDem: Average Elasticity", "UKIP: Average Elasticity", "Conservative: Elasticity at the Mean", "Labour: Elasticity at the Mean", "LibDem: Elasticity at the Mean", "UKIP: Elasticity at the Mean")
print(round(marginal.effects, digits=3))

###############
## Table 5.5 ##
###############

beta <- bes2015.mnl$coefficients
covmat.beta <- vcov(bes2015.mnl)
ndraws <- 1000
betadraw <- mvrnorm(ndraws, beta, covmat.beta)

taxspend.coeffs <- betadraw[,grep("taxspend", colnames(betadraw))]

## AME ##

observed.x <- with(bes.data, data.frame(
cons1 = as.numeric(idx[2]==1),
cons2 = as.numeric(idx[2]==2),
cons3 = as.numeric(idx[2]==3),
taxspend = taxspend, 
approve_EU1 = approve_EU * as.numeric(idx[2]==1), 
approve_EU2 = approve_EU * as.numeric(idx[2]==2), 
approve_EU3 = approve_EU * as.numeric(idx[2]==3), 
NHS_improved1 = NHS_improved * as.numeric(idx[2]==1),
NHS_improved2 = NHS_improved * as.numeric(idx[2]==2),
NHS_improved3 = NHS_improved * as.numeric(idx[2]==3),
income1 = income * as.numeric(idx[2]==1),
income2 = income * as.numeric(idx[2]==2),
income3 = income * as.numeric(idx[2]==3),
female1 = female * as.numeric(idx[2]==1),
female2 = female * as.numeric(idx[2]==2),
female3 = female * as.numeric(idx[2]==3),
age1 = age * as.numeric(idx[2]==1),
age2 = age * as.numeric(idx[2]==2),
age3 = age * as.numeric(idx[2]==3),
id = idx,
weights = wt_combined_main
))

party4.taxspend <- observed.x[observed.x$id.alt==4, grep("taxspend", colnames(observed.x))]

exb.party1 <- exp(betadraw %*% t(subset(observed.x, id.alt==1, select=c(1:19))))
exb.party2 <- exp(betadraw %*% t(subset(observed.x, id.alt==2, select=c(1:19))))
exb.party3 <- exp(betadraw %*% t(subset(observed.x, id.alt==3, select=c(1:19))))
exb.party4 <- exp(betadraw %*% t(subset(observed.x, id.alt==4, select=c(1:19))))
exb.all <- exb.party1 + exb.party2 + exb.party3 + exb.party4

probs1 <- exb.party1/exb.all
probs2 <- exb.party2/exb.all
probs3 <- exb.party3/exb.all
probs4 <- exb.party4/exb.all

id.weights <- observed.x$weights[observed.x$id.alt==1]/observed.x$weights[observed.x$id.alt==1]

## marginal effects for UKIP
me.taxspend <- taxspend.coeffs * (probs4 * (1 - probs4))
me.AME.taxspend <- colSums(apply(me.taxspend, 1, function(x){x * id.weights}))/sum(id.weights)

## cross-marginal effects
cme.taxspend41 <- taxspend.coeffs * (probs1 * probs4)
cme.taxspend42 <- taxspend.coeffs * (probs2 * probs4)
cme.taxspend43 <- taxspend.coeffs * (probs3 * probs4)

cme.AME.taxspend41 <- -1 * colSums(apply(cme.taxspend41, 1, function(x){x * id.weights}))/sum(id.weights)
cme.AME.taxspend42 <- -1 * colSums(apply(cme.taxspend42, 1, function(x){x * id.weights}))/sum(id.weights)
cme.AME.taxspend43 <- -1 * colSums(apply(cme.taxspend43, 1, function(x){x * id.weights}))/sum(id.weights)

cme.AME.taxspend <- cbind(cme.AME.taxspend41, cme.AME.taxspend42, cme.AME.taxspend43)

## elasticities for UKIP
e.taxspend <- taxspend.coeffs * t(party4.taxspend  * t(1 - probs4))
e.AME.taxspend <- colSums(apply(e.taxspend, 1, function(x){x * id.weights}))/sum(id.weights)

## cross-elasticities
ce.taxspend41 <- taxspend.coeffs * t(party4.taxspend * t(probs4))
ce.taxspend42 <- taxspend.coeffs * t(party4.taxspend * t(probs4))
ce.taxspend43 <- taxspend.coeffs * t(party4.taxspend * t(probs4))

ce.AME.taxspend41 <- -1 * colSums(apply(ce.taxspend41, 1, function(x){x * id.weights}))/sum(id.weights)
ce.AME.taxspend42 <- -1 * colSums(apply(ce.taxspend42, 1, function(x){x * id.weights}))/sum(id.weights)
ce.AME.taxspend43 <- -1 * colSums(apply(ce.taxspend43, 1, function(x){x * id.weights}))/sum(id.weights)

ce.AME.taxspend <- cbind(ce.AME.taxspend41, ce.AME.taxspend42, ce.AME.taxspend43)

## MEM ##

mean.x <- with(bes.data, data.frame(
constant = rbind(0, diag(3)),
taxspend = tapply(taxspend, idx(bes2015.mnl, 2), mean), 
approve_EU = mean(approve_EU) * rbind(0, diag(3)), 
NHS_improved = mean(NHS_improved) * rbind(0, diag(3)),
income=mean(income) * rbind(0, diag(3)),
female=mean(female) * rbind(0, diag(3)),
age=mean(age) * rbind(0, diag(3))
))

mean.x <- mean.x[order(rownames(mean.x)),]

taxspendcol <- grep("taxspend", colnames(mean.x))

exb.mean <- exp(betadraw %*% t(mean.x))

probs.mean <- exb.mean/rowSums(exb.mean)

## marginal effects for UKIP
me.MEM.taxspend <- taxspend.coeffs * probs.mean[,colnames(probs.mean)=="4"] * (1 - probs.mean[,colnames(probs.mean)=="4"])

## cross-marginal effects
cme.MEM.taxspend41 <- -1 * taxspend.coeffs * probs.mean[,colnames(probs.mean)=="4"] * probs.mean[,colnames(probs.mean)=="1"]
cme.MEM.taxspend42 <- -1 * taxspend.coeffs * probs.mean[,colnames(probs.mean)=="4"] * probs.mean[,colnames(probs.mean)=="2"]
cme.MEM.taxspend43 <- -1 * taxspend.coeffs * probs.mean[,colnames(probs.mean)=="4"] * probs.mean[,colnames(probs.mean)=="3"]
cme.MEM.taxspend <- cbind(cme.MEM.taxspend41, cme.MEM.taxspend42, cme.MEM.taxspend43)

## elasticities for UKIP
e.MEM.taxspend <- taxspend.coeffs * mean.x[rownames(mean.x)=="4",taxspendcol] * (1 - probs.mean[,colnames(probs.mean)=="4"])

## cross-elasticities
ce.MEM.taxspend41 <- -1 * taxspend.coeffs * mean.x[rownames(mean.x)=="4",taxspendcol] * probs.mean[,colnames(probs.mean)=="4"]
ce.MEM.taxspend42 <- -1 * taxspend.coeffs * mean.x[rownames(mean.x)=="4",taxspendcol] * probs.mean[,colnames(probs.mean)=="4"]
ce.MEM.taxspend43 <- -1 * taxspend.coeffs * mean.x[rownames(mean.x)=="4",taxspendcol] * probs.mean[,colnames(probs.mean)=="4"]
ce.MEM.taxspend <- cbind(ce.MEM.taxspend41, ce.MEM.taxspend42, ce.MEM.taxspend43)

all.results <- cbind(me.AME.taxspend, me.MEM.taxspend, e.AME.taxspend, e.MEM.taxspend, cme.AME.taxspend, cme.MEM.taxspend, ce.AME.taxspend, ce.MEM.taxspend)

quantiles <- apply(all.results, 2, function(x) {quantile(x, probs=c(0.5, 0.025, 0.975))})
sds <- apply(all.results, 2, sd)
marginal.effects <- cbind(t(quantiles), sds)
colnames(marginal.effects) <- c("Median", "2.5%", "97.5%", "(se)")
rownames(marginal.effects) <- c("UKIP: AME", "UKIP: MEM", "UKIP: Average Elasticity", "UKIP: Elasticity at the Mean", "Conservative: Average CME", "Labour: Average CME", "LibDem: Average CME", "Conservative: CME at Mean", "Labour: CME at Mean", "LibDem: CME at Mean", "Conservative: Average CE", "Labour: Average CE", "LibDem: Average CE", "Conservative: CE at Mean", "Labour: CE at Mean", "LibDem: CE at Mean")
print(round(marginal.effects, digits=3))


###############
## Table 5.6 ##
###############

N <- length(unique(bes.data$respondent_id))
A <- nrow(bes.data)
beta <- bes2015.mnl$coefficients
covmat.beta <- vcov(bes2015.mnl)
ndraws <- 1000

betadraw <- mvrnorm(ndraws, beta, covmat.beta)
exp.betadraw <- exp(betadraw)
quantiles <- apply(exp.betadraw, 2, function(x) {quantile(x, probs=c(0.5, 0.025, 0.975))})
sds <- apply(exp.betadraw, 2, sd)
odds.ratios <- cbind(t(quantiles), sds)
colnames(odds.ratios) <- c("Median", "2.5%", "97.5%", "(se)")
rownames(odds.ratios) <- names(beta) 
print(round(odds.ratios, digits=3))
print(N)
print(A)

