###############
## Section 4 ##
###############

#install.packages("ordinal")
library(ordinal)
library(MASS)

setwd("C:\\dc_book\\Replication_Files_Element")

financedata <- read.csv("ANES_2016.csv", header=T) 


###############
## Table 4.1 ##
###############

finance.ologit <- clm(as.ordered(prospective_persfin) ~ pred_favored_win+health+mortgage+gender+age, data = financedata, link="logit")
finance.ologit.results <- cbind(finance.ologit$coefficients, sqrt(diag(vcov(finance.ologit))), confint.default(finance.ologit, level=0.95))
N <- length(finance.ologit$fitted.values)
colnames(finance.ologit.results) <- c("Coeff", "(se)", "2.5%", "97.5%")
print((round(finance.ologit.results, digits=3)))
print(N)

## Note: Warning about formula can be ignored -- clm mistakenly throws an error for long formulas

###############
## Table 4.2 ##
###############

beta <- finance.ologit$coefficients
covmat.beta <- vcov(finance.ologit)
ndraws <- 1000
betadraw <- mvrnorm(ndraws, beta, covmat.beta)

baseline.x <- apply(finance.ologit$model[,-1], 2, median)

hypo.case.A <- baseline.x

hypo.case.B <- baseline.x
hypo.case.B["pred_favored_win"] <- 0

# multiply by -1 because the formula is tau - XB
hypo.matrix.A <- matrix(hypo.case.A, 4, 5, byrow=TRUE) * -1
hypo.matrix.tau.A <- cbind(diag(4), hypo.matrix.A)

hypo.matrix.B <- matrix(hypo.case.B, 4, 5, byrow=TRUE) * -1
hypo.matrix.tau.B <- cbind(diag(4), hypo.matrix.B)


simprob.A <- plogis(betadraw%*%t(hypo.matrix.tau.A))
simprob.A.1 <- simprob.A[,1]
simprob.A.2 <- simprob.A[,2] - simprob.A[,1]
simprob.A.3 <- simprob.A[,3] - simprob.A[,2]
simprob.A.4 <- simprob.A[,4] - simprob.A[,3]
simprob.A.5 <- 1 - simprob.A[,4]

simprob.B <- plogis(betadraw%*%t(hypo.matrix.tau.B))
simprob.B.1 <- simprob.B[,1]
simprob.B.2 <- simprob.B[,2] - simprob.B[,1]
simprob.B.3 <- simprob.B[,3] - simprob.B[,2]
simprob.B.4 <- simprob.B[,4] - simprob.B[,3]
simprob.B.5 <- 1 - simprob.B[,4]


simprob.diff.1 <- simprob.B.1 - simprob.A.1
simprob.diff.2 <- simprob.B.2 - simprob.A.2
simprob.diff.3 <- simprob.B.3 - simprob.A.3
simprob.diff.4 <- simprob.B.4 - simprob.A.4
simprob.diff.5 <- simprob.B.5 - simprob.A.5


all.results <- cbind(
simprob.A.1, simprob.B.1, simprob.diff.1, 
simprob.A.2, simprob.B.2, simprob.diff.2, 
simprob.A.3, simprob.B.3, simprob.diff.3, 
simprob.A.4, simprob.B.4, simprob.diff.4, 
simprob.A.5, simprob.B.5, simprob.diff.5)


quantiles <- apply(all.results, 2, function(x) {quantile(x, probs=c(0.5, 0.025, 0.975))})
sds <- apply(all.results, 2, sd)
first.differences <- cbind(t(quantiles), sds)
colnames(first.differences) <- c("Median", "2.5%", "97.5%", "(se)")
rownames(first.differences) <- c("Much Worse:Yes", "Much Worse:No", "Much Worse:Diff", "Somewhat Worse:Yes", "Somewhat Worse:No", "Somewhat Worse:Diff", "The Same:Yes", "The Same:No", "The Same:Diff", "Somewhat Better:Yes", "Somewhat Better:No", "Somewhat Better:Diff", "Much Better:Yes", "Much Better:No", "Much Better:Diff")
print(round(first.differences, digits=3))

################
## Figure 4.1 ##
################

beta <- finance.ologit$coefficients
covmat.beta <- vcov(finance.ologit)
ndraws <- 1000
betadraw <- mvrnorm(ndraws, beta, covmat.beta)

baseline.x <- apply(finance.ologit$model[,-1], 2, median)
baseline.matrix <- matrix(baseline.x, 4, 5, byrow=TRUE) * -1  # multiply by -1 because the formula is t - XB
baseline.matrix.tau <- cbind(diag(4), baseline.matrix)
colnames(baseline.matrix.tau) <- names(finance.ologit$coefficients)


colors <- c("#004D7380", "#1C82A380", "#17A1CC80", "#0ABA9E80", "#00DBA680")

age.range <- c(10:100)


ologit.probs1 <- matrix(NA, ndraws, length(age.range))
ologit.probs2 <- matrix(NA, ndraws, length(age.range))
ologit.probs3 <- matrix(NA, ndraws, length(age.range))
ologit.probs4 <- matrix(NA, ndraws, length(age.range))
ologit.probs5 <- matrix(NA, ndraws, length(age.range))

ologit.mean.probs <- matrix(NA, 5, length(age.range))

for (i in 1:length(age.range)) {

	baseline.matrix.tau[,"age"] <- (age.range[i] * -1)
	oprobs <- plogis(betadraw%*%t(baseline.matrix.tau))
	ologit.probs1[,i] <- oprobs[,1]
	ologit.probs2[,i] <- oprobs[,2] - oprobs[,1]
	ologit.probs3[,i] <- oprobs[,3] - oprobs[,2]
	ologit.probs4[,i] <- oprobs[,4] - oprobs[,3]
	ologit.probs5[,i] <- 1 - oprobs[,4]
	
	oprobs.mean <- plogis(beta%*%t(baseline.matrix.tau))
	ologit.mean.probs[1,i] <- oprobs.mean[1]
	ologit.mean.probs[2,i] <- oprobs.mean[2] - oprobs.mean[1]
	ologit.mean.probs[3,i] <- oprobs.mean[3] - oprobs.mean[2]
	ologit.mean.probs[4,i] <- oprobs.mean[4] - oprobs.mean[3]
	ologit.mean.probs[5,i] <- 1 - oprobs.mean[4]
}

pdf("fig41.pdf")
par(mar=c(5.1, 4.1, 4.1, 5.6), mgp=c(2.5, 1, 0), cex.lab=1.5)
plot(ologit.mean.probs[1,]~age.range, xlim=c(18,90), ylim=c(0,0.6), type="n", xlab="Age", ylab="Probability", xaxt="n", las=1, lwd=2)
axis(1, at=c(18, 30, 40, 50, 60, 70, 80, 90))
for (i in 1:length(age.range)){
	points(ologit.probs1[i,] ~ age.range, type="l", col=colors[1], lwd=1)
	points(ologit.probs2[i,] ~ age.range, type="l", col=colors[2], lwd=1)
	points(ologit.probs3[i,] ~ age.range, type="l", col=colors[3], lwd=1)
	points(ologit.probs4[i,] ~ age.range, type="l", col=colors[4], lwd=1)
	points(ologit.probs5[i,] ~ age.range, type="l", col=colors[5], lwd=1)
}
lines(ologit.mean.probs[1,]~age.range, col=rgb(0,0,0,0.5), lwd=2)
lines(ologit.mean.probs[2,]~age.range, col=rgb(0,0,0,0.5), lwd=2)
lines(ologit.mean.probs[3,]~age.range, col=rgb(0,0,0,0.5), lwd=2)
lines(ologit.mean.probs[4,]~age.range, col=rgb(0,0,0,0.5), lwd=2)
lines(ologit.mean.probs[5,]~age.range, col=rgb(0,0,0,0.5), lwd=2)
mtext("Much Worse", side=4, at=ologit.mean.probs[1,82], cex=1.1, las=1, line=0.5)
mtext("Somewhat", side=4, at=ologit.mean.probs[2,82]-0.01, cex=1.1, las=1, line=0.5)
mtext("Worse", side=4, at=ologit.mean.probs[2,82]-0.03, cex=1.1, las=1, line=0.5)
mtext("The Same", side=4, at=ologit.mean.probs[3,82], cex=1.1, las=1, line=0.5)
mtext("Somewhat", side=4, at=ologit.mean.probs[4,82]+0.03, cex=1.1, las=1, line=0.5)
mtext("Better", side=4, at=ologit.mean.probs[4,82]+0.01, cex=1.1, las=1, line=0.5)
mtext("Much Better", side=4, at=ologit.mean.probs[5,82], cex=1.1, las=1, line=0.5)
dev.off()


###############
## Table 4.3 ##
###############

beta <- finance.ologit$coefficients
covmat.beta <- vcov(finance.ologit)
ndraws <- 1000
betadraw <- mvrnorm(ndraws, beta, covmat.beta)

# AME

ologit.X <- cbind(0,0,0,0, finance.ologit$model[,-1]) * -1  # multiply by -1 because the formula is t - XB
w <- financedata$sample_weight


for (i in 1:4) {

	ologit.X.temp <- ologit.X
	ologit.X.temp[,i] <- 1
	ologit.p.temp <- plogis(betadraw %*% t(ologit.X.temp))
	assign(paste("ologit.p", i, sep="."), ologit.p.temp)
	assign(paste("ologit.d", i, sep="."), ologit.p.temp * (1 - ologit.p.temp))
}

# multiply each row by draws of health coefficients
ame.1 <- diag(betadraw[,"health"]) %*% (0 - ologit.d.1) 
ame.2 <- diag(betadraw[,"health"]) %*% (ologit.d.1 - ologit.d.2) 
ame.3 <- diag(betadraw[,"health"]) %*% (ologit.d.2 - ologit.d.3) 
ame.4 <- diag(betadraw[,"health"]) %*% (ologit.d.3 - ologit.d.4) 
ame.5 <- diag(betadraw[,"health"]) %*% (ologit.d.4) 

ologit.ame.1 <- apply(ame.1, 1, function(x){sum(x * w)/sum(w)})
ologit.ame.2 <- apply(ame.2, 1, function(x){sum(x * w)/sum(w)})
ologit.ame.3 <- apply(ame.3, 1, function(x){sum(x * w)/sum(w)})
ologit.ame.4 <- apply(ame.4, 1, function(x){sum(x * w)/sum(w)})
ologit.ame.5 <- apply(ame.5, 1, function(x){sum(x * w)/sum(w)})

ologit.ame <- cbind(ologit.ame.1, ologit.ame.2, ologit.ame.3, ologit.ame.4, ologit.ame.5)

ae.1 <- ame.1 * ((1/ologit.p.1) %*% diag(finance.ologit$model[,"health"]))
ae.2 <- ame.2 * ((1/(ologit.p.2 - ologit.p.1)) %*% diag(finance.ologit$model[,"health"]))
ae.3 <- ame.3 * ((1/(ologit.p.3 - ologit.p.2)) %*% diag(finance.ologit$model[,"health"]))
ae.4 <- ame.4 * ((1/(ologit.p.4 - ologit.p.3)) %*% diag(finance.ologit$model[,"health"]))
ae.5 <- ame.5 * ((1/(1 - ologit.p.4)) %*% diag(finance.ologit$model[,"health"]))

ologit.ae.1 <- apply(ae.1, 1, function(x){sum(x * w)/sum(w)})
ologit.ae.2 <- apply(ae.2, 1, function(x){sum(x * w)/sum(w)})
ologit.ae.3 <- apply(ae.3, 1, function(x){sum(x * w)/sum(w)})
ologit.ae.4 <- apply(ae.4, 1, function(x){sum(x * w)/sum(w)})
ologit.ae.5 <- apply(ae.5, 1, function(x){sum(x * w)/sum(w)})

ologit.ae <- cbind(ologit.ae.1, ologit.ae.2, ologit.ae.3, ologit.ae.4, ologit.ae.5)

quantiles.ame <- apply(ologit.ame, 2, function(x) {quantile(x, probs=c(0.5, 0.025, 0.975))})
sd.ame <- apply(ologit.ame, 2, sd)
summary.marginal.ame <- cbind(t(quantiles.ame), sd.ame)

quantiles.ae <- apply(ologit.ae, 2, function(x) {quantile(x, probs=c(0.5, 0.025, 0.975))})
sd.ae <- apply(ologit.ae, 2, sd)
summary.elasticity.ame <- cbind(t(quantiles.ae), sd.ae)

# MEM 

mean.x <- apply(finance.ologit$model[,-1], 2, mean)
mean.x.matrix <- matrix(mean.x, 4, 5, byrow=TRUE) * -1
mean.x.all <- cbind(diag(4), mean.x.matrix)

mem.probs <- plogis(betadraw%*%t(mean.x.all))

mem.1 <- -1 * (mem.probs[,1] * (1 - mem.probs[,1])) * betadraw[,"health"]
mem.2 <- ((mem.probs[,1] * (1 - mem.probs[,1])) - (mem.probs[,2] * (1 - mem.probs[,2]))) * betadraw[,"health"]
mem.3 <- ((mem.probs[,2] * (1 - mem.probs[,2])) - (mem.probs[,3] * (1 - mem.probs[,3]))) * betadraw[,"health"]
mem.4 <- ((mem.probs[,3] * (1 - mem.probs[,3])) - (mem.probs[,4] * (1 - mem.probs[,4]))) * betadraw[,"health"]
mem.5 <- (mem.probs[,4] * (1 - mem.probs[,4]))  * betadraw[,"health"]

ologit.mem <- cbind(mem.1, mem.2, mem.3, mem.4, mem.5)

me.1 <- mem.1 * (mean.x["health"]/(mem.probs[,1]))
me.2 <- mem.2 * (mean.x["health"]/(mem.probs[,2] - mem.probs[,1]))
me.3 <- mem.3 * (mean.x["health"]/(mem.probs[,3] - mem.probs[,2]))
me.4 <- mem.4 * (mean.x["health"]/(mem.probs[,4] - mem.probs[,3]))
me.5 <- mem.5 * (mean.x["health"]/(1 - mem.probs[,4]))

ologit.me <- cbind(me.1, me.2, me.3, me.4, me.5)

quantiles.mem <- apply(ologit.mem, 2, function(x) {quantile(x, probs=c(0.5, 0.025, 0.975))})
sd.mem <- apply(ologit.mem, 2, sd)
summary.marginal.mem <- cbind(t(quantiles.mem), sd.mem)

quantiles.me <- apply(ologit.me, 2, function(x) {quantile(x, probs=c(0.5, 0.025, 0.975))})
sd.me <- apply(ologit.me, 2, sd)
summary.elasticity.mem <- cbind(t(quantiles.me), sd.me)

marginal.effects.elasticities <- rbind(summary.marginal.ame, summary.marginal.mem, summary.elasticity.ame, summary.elasticity.mem)
colnames(marginal.effects.elasticities) <- c("Median", "2.5%", "97.5%", "(se)")
category.names <- c("MW", "SW", "S", "SB", "MB")
rownames(marginal.effects.elasticities) <- c(paste("Marginal Effect (AME)", category.names, sep="."), paste("Marginal Effect (MEM)", category.names, sep="."), paste("Elasticity (AME)", category.names, sep="."), paste("Elasticity (MEM)", category.names, sep="."))
print(round(marginal.effects.elasticities, digits=3))


## coefficient ratios

mortgage.age.ratio <- betadraw[,"mortgage"]/betadraw[,"age"]
ratio.results <- quantile(mortgage.age.ratio, probs=c(0.5, 0.025, 0.975))
print(ratio.results)

###############
## Table 4.4 ##
###############

N <- length(finance.ologit$fitted.values)

beta <- finance.ologit$coefficients
covmat.beta <- vcov(finance.ologit)
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

