#######################
## Section 6, Part 1 ##
#######################

#install.packages("pglm")
#install.packages("mlogit")
#install.packages("statmod")

library(pglm)
library(mlogit)
library(statmod)
library(MASS)

setwd("C:\\dc_book\\Replication_Files_Element")

Scot.panel <- read.csv("Scottish_Ref.csv")

# scotReferendumIntention 1 = leave, 0 = stay 
# euRefVote 1 = stay, 0 = leave (reverse coded from original survey question) 

###############
## Table 6.1 ##
###############


# quadrature
re.logit.quad <- pglm(scotReferendumIntention ~ euRefVote+postEUref+postEUVote+satDemScot+riskUnemployment+Scottish, model=("random"), effect=("individual"), R=64, index=c("id", "wave"), family=binomial(link="logit"), data=Scot.panel)
re.logit.quad.results <- cbind(re.logit.quad$estimate, sqrt(diag(vcov(re.logit.quad))), confint.default(re.logit.quad, level=0.95))
N <- length(unique(Scot.panel$id))
O <- nrow(Scot.panel)
colnames(re.logit.quad.results) <- c("Coeff", "(se)", "2.5%", "97.5%")
print((round(re.logit.quad.results, digits=3)))
print(N)
print(O)

# simulation

Scot.panel.long <- mlogit.data(Scot.panel, choice="scotReferendumIntention", shape="wide", id="id")
Scot.panel.long$constant <- as.numeric(as.character(Scot.panel.long$alt))
re.logit.sim <- mlogit(scotReferendumIntention ~ constant | 0+euRefVote+postEUref+postEUVote+satDemScot+riskUnemployment+Scottish, rpar=c(constant='n'), panel=TRUE, index="id", R=500, data=Scot.panel.long, reflevel="0")
re.logit.sim.results <- cbind(re.logit.sim$coefficients, sqrt(diag(vcov(re.logit.sim))), confint.default(re.logit.sim, level=0.95))
N <- length(unique(Scot.panel$id))
O <- nrow(Scot.panel)
colnames(re.logit.sim.results) <- c("Coeff", "(se)", "2.5%", "97.5%")
print((round(re.logit.sim.results, digits=3)))
print(N)
print(O)

# logit 

pooled.logit <- glm(scotReferendumIntention ~ euRefVote+postEUref+postEUVote+satDemScot+riskUnemployment+Scottish, family=binomial(link="logit"), data=Scot.panel)
pooled.logit.results <- cbind(pooled.logit$coefficients, sqrt(diag(vcov(pooled.logit))), confint.default(pooled.logit, level=0.95))
N <- length(unique(Scot.panel$id))
O <- nrow(Scot.panel)
colnames(pooled.logit.results) <- c("Coeff", "(se)", "2.5%", "97.5%")
print((round(pooled.logit.results, digits=3)))
print(N)
print(O)


################
## Figure 6.1 ##
################

beta.re <- re.logit.quad$estimate
beta.re0 <- re.logit.quad$estimate[-grep("sigma", names(re.logit.quad$estimate))]
beta.pooled <- pooled.logit$coefficients

ndraws <- 1000
beta <- re.logit.quad$estimate
covmat.beta <- vcov(re.logit.quad)
betadraw <- mvrnorm(ndraws, beta, covmat.beta)

betadraw.pooled <- mvrnorm(ndraws, beta.pooled, vcov(pooled.logit))
betadraw.re0 <- as.matrix(betadraw[,-grep("sigma", names(re.logit.quad$estimate))])

baseline.x <- with(Scot.panel, data.frame(
constant = 1,
euRefVote = median(euRefVote),
postEUref = median(postEUref),
postEUVote = median(postEUVote),
satDemScot = median(satDemScot),
riskUnemployment = median(riskUnemployment),
Scottish = median(Scottish)
))

satdemcol <- grep("satDemScot", colnames(baseline.x))
satdem.range <- seq(0, 5, 0.1)

probs.re <- matrix (0, ndraws, length(satdem.range))
probs.re0 <- matrix (0, ndraws, length(satdem.range))
probs.pooled <- matrix (0, ndraws, length(satdem.range))

probs.mean <- matrix(0, length(satdem.range), 3)

R <- 64
# the command "gauss.quad.prob" gives the "probabilist" version of the Hermite polynomial (the command "gauss.quad" gives the "physicist" version)
quadpoints <- gauss.quad.prob(R, dist="normal")
intpoints <- quadpoints$nodes
intweights <- quadpoints$weights



for (i in 1:length(satdem.range)) {

	baseline.x[satdemcol] <- satdem.range[i]
	
	for (j in 1:length(intpoints)) {
	
		x.re <- cbind(baseline.x, intpoints[j])
		p.re <- plogis(beta.re %*% t(x.re)) * intweights[j]
		probs.mean[i,1] <- probs.mean[i,1] + p.re
		
		p.draw.re <- plogis(betadraw %*% t(x.re)) * intweights[j]
		probs.re[,i] <- probs.re[,i] + p.draw.re
	}
	
	probs.mean[i,2] <- plogis(beta.pooled %*% t(baseline.x))
	probs.mean[i,3] <- plogis(beta.re0 %*% t(baseline.x))
	
	probs.pooled[,i] <- plogis(betadraw.pooled %*% t(baseline.x))
	probs.re0[,i] <- plogis(betadraw.re0 %*% t(baseline.x))
 
}

# random effects, pooled logit, random effect (v=0)
colors <- c("#FC660080", "#004D7380", "#86AC4180")

pdf("fig61.pdf")
par(mar=c(5.1, 4.1, 4.1, 5.5), mgp=c(2.5, 1, 0), cex.lab=1.5)
plot(probs.mean[,1] ~ satdem.range, xlim=c(1,4), ylim=c(0, 1), type="n", xaxt="n", xlab="Satisfaction with Democracy", ylab="Probability", xaxt="n", las=1)
axis(1, at=c(1, 2, 3, 4))
for (i in 1:length(satdem.range)){
	points(probs.re[i,] ~ satdem.range, type="l", col=colors[1], lwd=1)
	points(probs.pooled[i,] ~ satdem.range, type="l", col=colors[2], lwd=1)
	points(probs.re0[i,] ~ satdem.range, type="l", col=colors[3], lwd=1)
}
lines(probs.mean[,1]~satdem.range, col=rgb(0,0,0,0.5), lwd=2)
lines(probs.mean[,2]~satdem.range, col=rgb(0,0,0,0.5), lwd=2)
lines(probs.mean[,3]~satdem.range, col=rgb(0,0,0,0.5), lwd=2)
mtext("Random", side=4, at=probs.mean[42,1] + 0.035, cex=1.1, las=1, line=0.5)
mtext("Effects", side=4, at=probs.mean[42,1], cex=1.1, las=1, line=0.5)
mtext("Logit", side=4, at=(probs.mean[42,1] - 0.035), cex=1.1, las=1, line=0.5)
mtext("Pooled Logit", side=4, at=probs.mean[42,2], cex=1.1, las=1, line=0.5)
mtext("Random", side=4, at=probs.mean[42,3]+0.035, cex=1.1, las=1, line=0.5)
mtext("Effects", side=4, at=probs.mean[42,3], cex=1.1, las=1, line=0.5)
mtext(expression(paste("Logit (", alpha," = 0)", sep="")), side=4, at=(probs.mean[42,3] - 0.035), cex=1.1, las=1, line=0.5)
dev.off()




###############
## Table 6.2 ##
###############

ndraws <- 1000
beta <- re.logit.quad$estimate
covmat.beta <- vcov(re.logit.quad)
betadraw <- mvrnorm(ndraws, beta, covmat.beta)

baseline.x <- with(Scot.panel, data.frame(
constant = 1,
euRefVote = median(euRefVote),
postEUref = median(postEUref),
postEUVote = median(postEUVote),
satDemScot = median(satDemScot),
riskUnemployment = median(riskUnemployment),
Scottish = median(Scottish)
))

x.remain.pre <- baseline.x 
x.remain.pre[names(x.remain.pre)=="euRefVote"] <- 1
x.remain.pre[names(x.remain.pre)=="postEUref"] <- 0
x.remain.pre[names(x.remain.pre)=="postEUVote"] <- 0

x.remain.post <- baseline.x 
x.remain.post[names(x.remain.post)=="euRefVote"] <- 1
x.remain.post[names(x.remain.post)=="postEUref"] <- 1
x.remain.post[names(x.remain.post)=="postEUVote"] <- 1

x.leave.pre <- baseline.x
x.leave.pre[names(x.leave.pre)=="euRefVote"] <- 0
x.leave.pre[names(x.leave.pre)=="postEUref"] <- 0
x.leave.pre[names(x.leave.pre)=="postEUVote"] <- 0

x.leave.post <- baseline.x
x.leave.post[names(x.leave.post)=="euRefVote"] <- 0
x.leave.post[names(x.leave.post)=="postEUref"] <- 1
x.leave.post[names(x.leave.post)=="postEUVote"] <- 0

probs <- matrix(0, ndraws, 7)

## quadrature
R <- 64

quadpoints <- gauss.quad.prob(64, dist="normal")
intpoints <- quadpoints$nodes
intweights <- quadpoints$weights

## could use simulation instead
#install.packages("randtoolbox")
#library(randtoolbox)
#R <- 250
#halton.seq <- halton(R)
#intpoints <- qnorm(halton.seq)
#intweights <- rep(1/R, R)

for (i in 1:R) {

	x.remain.pre.q <- cbind(x.remain.pre, intpoints[i])
	x.remain.post.q <- cbind(x.remain.post, intpoints[i])
	x.leave.pre.q <- cbind(x.leave.pre, intpoints[i])
	x.leave.post.q <- cbind(x.leave.post, intpoints[i])

	p.remain.pre <- plogis(betadraw %*% t(x.remain.pre.q)) * intweights[i]
	p.remain.post <- plogis(betadraw %*% t(x.remain.post.q)) * intweights[i]
	p.leave.pre <- plogis(betadraw %*% t(x.leave.pre.q)) * intweights[i]
	p.leave.post <- plogis(betadraw %*% t(x.leave.post.q)) * intweights[i]
	
	probs[,1] <- probs[,1] + p.remain.pre
	probs[,2] <- probs[,2] + p.remain.post
	probs[,3] <- probs[,3] + p.leave.pre
	probs[,4] <- probs[,4] + p.leave.post	
 
}

probs[,5] <- probs[,2] - probs[,1]
probs[,6] <- probs[,4] - probs[,3]
probs[,7] <- probs[,5] - probs[,6]

all.results <- probs	
quantiles <- apply(all.results, 2, function(x) {quantile(x, probs=c(0.5, 0.025, 0.975))})
sds <- apply(all.results, 2, sd)
first.differences <- cbind(t(quantiles), sds)
colnames(first.differences) <- c("Median", "2.5%", "97.5%", "(se)")
rownames(first.differences) <- c("Remain: Pre Brexit", "Remain: Post Brexit", "Remain: Difference", "Leave: Pre Brexit", "Leave: Post Brexit", "Leave: Difference", "Diff in Diff")
print(round(first.differences, digits=3))


###############
## Table 6.3 ##
###############

ndraws <- 1000
beta <- re.logit.quad$estimate
covmat.beta <- vcov(re.logit.quad)
betadraw <- mvrnorm(ndraws, beta, covmat.beta)

R <- 64
quadpoints <- gauss.quad.prob(R, dist="normal")
intpoints <- quadpoints$nodes
intweights <- quadpoints$weights

observed.x <- with(Scot.panel, data.frame(
constant = 1, 
euRefVote = euRefVote,
postEUref = postEUref,
postEUVote = postEUVote,
satDemScot = satDemScot,
riskUnemployment = riskUnemployment,
Scottish = Scottish
))

mean.x <- colMeans(observed.x)

id.weights <- Scot.panel$sample_weight

mem <- 0
ame.all <- matrix(0, nrow(betadraw), nrow(observed.x))

mem.predprob <- matrix(0, nrow(betadraw), 1)
ame.predprob <- matrix(0, nrow(betadraw), nrow(observed.x))

for (i in 1:R) {

	# MEM
	x.mean.all <- cbind(t(mean.x), intpoints[i])
	predprob.mean <- plogis(betadraw %*% t(x.mean.all))
	predprob.mean.q <- predprob.mean * intweights[i]
	mem.predprob <- mem.predprob + (predprob.mean.q)
	mem.temp <- predprob.mean * (1 - predprob.mean) * betadraw[,"satDemScot"] 
	mem <- mem + (mem.temp * intweights[i])

	# AME
	x.vars.all <- cbind(observed.x, intpoints[i])
	predprobs <- plogis(betadraw %*% t(x.vars.all)) 
	ame.predprob <- ame.predprob + (predprobs * intweights[i])
	ame.temp <- predprobs * (1 - predprobs) * betadraw[,"satDemScot"]  
	ame.all <- ame.all + (ame.temp * intweights[i])
 
}


ame <- colSums(apply(ame.all, 1, function(x){x * id.weights}))/sum(id.weights)
elas <- t(ame.all) * observed.x[,"satDemScot"] * t(1/ame.predprob)
avg.elas <- colSums(apply(elas, 2, function(x){x * id.weights}))/sum(id.weights)
mem.elas <- mem * (mean.x["satDemScot"]/mem.predprob)

all.results <- cbind(ame, mem, avg.elas, mem.elas)
quantiles <- apply(all.results, 2, function(x) {quantile(x, probs=c(0.5, 0.025, 0.975))})
sds <- apply(all.results, 2, sd)
marginal.effects <- cbind(t(quantiles), sds)
colnames(marginal.effects) <- c("Median", "2.5%", "97.5%", "(se)")
rownames(marginal.effects) <- c("AME", "MEM", "Average Elasticity", "Elasticity at the Mean")
print(round(marginal.effects, digits=3))


###############
## Table 6.4 ##
###############

## Population Odds Ratios

ndraws <- 1000
beta <- re.logit.quad$estimate
covmat.beta <- vcov(re.logit.quad)
betadraw <- mvrnorm(ndraws, beta, covmat.beta)

R <- 64
quadpoints <- gauss.quad.prob(R, dist="normal")
intpoints <- quadpoints$nodes
intweights <- quadpoints$weights

observed.x <- with(Scot.panel, data.frame(
Intercept = 1,
euRefVote = euRefVote,
postEUref = postEUref,
postEUVote = postEUVote,
satDemScot = satDemScot,
riskUnemployment = riskUnemployment,
Scottish = Scottish
))

id.weights <- Scot.panel$sample_weight 

OR.pop <- matrix(NA, ncol(observed.x), 4)

p.0 <- matrix(0, nrow(betadraw), nrow(observed.x))

# calculate baseline
for (j in 1:R) {

	x.vars.all <- cbind(observed.x, intpoints[j])
	p0 <- plogis(betadraw %*% t(x.vars.all)) 
	p.0 <- p.0 + (p0 * intweights[j])
	
}

# increment each variable by 1, calculate odds ratios
for (i in 1:ncol(observed.x)) {

	p.1 <- matrix(0, nrow(betadraw), nrow(observed.x))
	observed.x1 <- observed.x
	observed.x1[,i] <- observed.x1[,i] + 1
	
	if (i==4) { 
	observed.x1[,2] <- observed.x1[,2] + 1
	observed.x1[,3] <- observed.x1[,3] + 1
	}

	for (j in 1:length(intpoints)) {

		x1.vars.all <- cbind(observed.x1, intpoints[j])
		p1 <- plogis(betadraw %*% t(x1.vars.all)) 
		p.1 <- p.1 + (p1 * intweights[j])
		
	}

		OR.temp <- (p.1 / (1 - p.1)) / (p.0 / (1 - p.0))
		OR.temp.w <- colSums(apply(OR.temp, 1, function(x){x * id.weights}))/sum(id.weights)
		OR.pop[i,1:3] <- quantile(OR.temp.w, probs=c(0.5, 0.025, 0.975))
		OR.pop[i,4] <- sd(OR.temp)

}

N <- length(unique(Scot.panel$id))
O <- nrow(Scot.panel)
colnames(OR.pop) <- c("Median", "2.5%", "97.5%", "(se)")
rownames(OR.pop) <- names(observed.x)
print((round(OR.pop, digits=3)))
print(N)
print(O)



## Subject-Specific Odds ratios

betadraw.or <- betadraw[,-grep("sigma|(Intercept)|postEUVote", colnames(betadraw))]
postEUVote <- betadraw[,grep("euRefVote", colnames(betadraw))] + betadraw[,grep("postEUref", colnames(betadraw))] + betadraw[,grep("postEUVote", colnames(betadraw))]
betadraw.or <- cbind(betadraw.or, postEUVote)

exp.betadraw <- exp(betadraw.or)
N <- length(unique(Scot.panel$id))
O <- nrow(Scot.panel)
quantiles <- apply(exp.betadraw, 2, function(x) {quantile(x, probs=c(0.5, 0.025, 0.975))})
sds <- apply(exp.betadraw, 2, sd)
odds.ratios <- cbind(t(quantiles), sds)
colnames(odds.ratios) <- c("Median", "2.5%", "97.5%", "(se)")
rownames(odds.ratios) <- colnames(exp.betadraw) 
print(round(odds.ratios, digits=3))
print(N)
print(O)


## approximation from Zeger et al. 1988, p. 1054 ##

c <- ((16*sqrt(3))/(15*pi))
scale.pa <- sqrt(1 + c^2 * re.logit.quad$estimate[grep("sigma", names(re.logit.quad$estimate))]^2)
beta.Zeger <- apply(betadraw.or, 2, function(x){x/scale.pa})
exp.beta.Zeger <- exp(beta.Zeger)
quantiles.Zeger <- apply(exp.beta.Zeger, 2, function(x) {quantile(x, probs=c(0.5, 0.025, 0.975))})
sds.Zeger <- apply(exp.beta.Zeger, 2, sd)
Zeger.adj <- cbind(t(quantiles.Zeger), sds.Zeger)
colnames(Zeger.adj) <- c("Median", "2.5%", "97.5%", "(se)")
rownames(Zeger.adj) <- colnames(exp.beta.Zeger) 
print(round(Zeger.adj, digits=3))



## approximation from Hedeker et al. 2018 ## 

betadraw.H <- betadraw[,-grep("sigma|postEUVote", colnames(betadraw))]
postEUVote <- betadraw[,grep("euRefVote", colnames(betadraw))] + betadraw[,grep("postEUref", colnames(betadraw))] + betadraw[,grep("postEUVote", colnames(betadraw))]
sigma <- betadraw[,grep("sigma", colnames(betadraw))]
betadraw.H <- cbind(betadraw.H, postEUVote, sigma)

observed.H <- with(Scot.panel, data.frame(
Intercept = 1,
euRefVote = euRefVote,
postEUref = postEUref,
satDemScot = satDemScot,
riskUnemployment = riskUnemployment,
Scottish = Scottish,
postEUVote = postEUVote
))

predprobs <- matrix(0, nrow(betadraw.H), nrow(observed.H))

for (i in 1:R) {

	x.vars.all <- cbind(observed.H, intpoints[i])
	predprobs.temp <- plogis(betadraw.H %*% t(x.vars.all)) 
	predprobs <- predprobs + (predprobs.temp * intweights[i])
 
}

X <- as.matrix(observed.H)
lambda.matrix <- log(predprobs/(1 - predprobs))
beta.Hedeker <- apply(lambda.matrix, 1, function(x){solve(t(X) %*% X) %*% (t(X) %*% x)}) 
exp.beta.Hedeker <- exp(beta.Hedeker)
quantiles.Hedeker <- apply(exp.beta.Hedeker, 1, function(x) {quantile(x, probs=c(0.5, 0.025, 0.975))})
sds.Hedeker <- apply(exp.beta.Hedeker, 1, sd)
Hedeker.adj <- cbind(t(quantiles.Hedeker), sds.Hedeker)
colnames(Hedeker.adj) <- c("Median", "2.5%", "97.5%", "(se)")
rownames(Hedeker.adj) <- colnames(observed.H)
Hedeker.adj <- Hedeker.adj[-1,]  # remove constant term
print(round(Hedeker.adj, digits=3))

