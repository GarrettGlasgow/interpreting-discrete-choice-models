#######################
## Section 6, Part 2 ##
#######################

#install.packages("mlogit")
#install.packages("randtoolbox")

library(mlogit)
library(randtoolbox)
library(MASS)

setwd("C:\\dc_book\\Replication_Files_Element")

cb.data <- read.csv("cb_data.csv")

cb.data <- mlogit.data(cb.data, choice="choice", shape="long", chid.var="choiceid", alt.var="region", id.var="id")

###############
## Table 6.5 ##
###############

cb.data$n.travelcost <- cb.data$travelcost * -1

cosco.busan.mxl <- mlogit(choice ~ n.travelcost + bridge + marin + eastbay + bayshore + peninsula | 0, data=cb.data, rpar=c(n.travelcost="ln", bridge="n", marin="n", eastbay="n", bayshore="n", peninsula="n"), panel=TRUE, R=500, halton=NA)
cosco.busan.mxl.results <- cbind(cosco.busan.mxl$coefficients, sqrt(diag(vcov(cosco.busan.mxl))), confint.default(cosco.busan.mxl, level=0.95))
N <- length(unique(cb.data$id))
C <- length(unique(cb.data$choiceid))
A <- nrow(cb.data)
colnames(cosco.busan.mxl.results) <- c("Coeff", "(se)", "2.5%", "97.5%")
print((round(cosco.busan.mxl.results, digits=3)))  # interpret negative standard deviations as if they were positive
print(N)
print(C)
print(A)


###############
## Table 6.6 ## 
###############

coeffs <- cosco.busan.mxl$coefficients
covmat <- vcov(cosco.busan.mxl)
ndraws <- 1000

betadraw <- mvrnorm(ndraws, coeffs, covmat)
R <- 500
halton.seq <- qnorm(halton(R, dim=6))

## observed

observed.x <- with(cb.data, data.frame(
travelcost = travelcost,
bridge = bridge,
marin = marin,
eastbay = eastbay, 
bayshore = bayshore,
peninsula = peninsula,
id = idx, 
weights = wts
))

id.weights <- observed.x$weights[observed.x$id.alt==1]

shares <- matrix(0, ndraws, 15)

marin.o <- matrix(0, ndraws, nrow(subset(observed.x, id.alt==1)))
eastbay.o <- matrix(0, ndraws, nrow(subset(observed.x, id.alt==2)))
bayshore.o <- matrix(0, ndraws, nrow(subset(observed.x, id.alt==3)))
peninsula.o <- matrix(0, ndraws, nrow(subset(observed.x, id.alt==4)))
notrip.o <- matrix(0, ndraws, nrow(subset(observed.x, id.alt==5)))
marin.c <- matrix(0, ndraws, nrow(subset(observed.x, id.alt==1)))
eastbay.c <- matrix(0, ndraws, nrow(subset(observed.x, id.alt==2)))
bayshore.c <- matrix(0, ndraws, nrow(subset(observed.x, id.alt==3)))
peninsula.c <- matrix(0, ndraws, nrow(subset(observed.x, id.alt==4)))
notrip.c <- matrix(0, ndraws, nrow(subset(observed.x, id.alt==5)))

for (i in 1:R) {

	travelcost <- -1 * exp(betadraw[,"n.travelcost"] + (halton.seq[i,1] * abs(betadraw[,"sd.n.travelcost"])))
	bridge <- betadraw[,"bridge"] + (halton.seq[i,2] * abs(betadraw[,"sd.bridge"]))
	marin <- betadraw[,"marin"] + (halton.seq[i,3] * abs(betadraw[,"sd.marin"]))
	eastbay <- betadraw[,"eastbay"] + (halton.seq[i,4] * abs(betadraw[,"sd.eastbay"]))
	bayshore <- betadraw[,"bayshore"] + (halton.seq[i,5] * abs(betadraw[,"sd.bayshore"]))
	peninsula <- betadraw[,"peninsula"] + (halton.seq[i,6] * abs(betadraw[,"sd.peninsula"]))

	beta.R <- cbind(travelcost, bridge, marin, eastbay, bayshore, peninsula)
	
	exb.region1.o <- exp(beta.R %*% t(subset(observed.x, id.alt==1, select=c(1:6))))
	exb.region2.o <- exp(beta.R %*% t(subset(observed.x, id.alt==2, select=c(1:6))))
	exb.region3.o <- exp(beta.R %*% t(subset(observed.x, id.alt==3, select=c(1:6))))
	exb.region4.o <- exp(beta.R %*% t(subset(observed.x, id.alt==4, select=c(1:6))))
	exb.region5.o <- exp(beta.R %*% t(subset(observed.x, id.alt==5, select=c(1:6))))
	exb.all.o <- exb.region1.o + exb.region2.o + exb.region3.o + exb.region4.o + exb.region5.o
	
	probs.region1.o <- exb.region1.o/exb.all.o
	probs.region2.o <- exb.region2.o/exb.all.o
	probs.region3.o <- exb.region3.o/exb.all.o
	probs.region4.o <- exb.region4.o/exb.all.o
	probs.region5.o <- exb.region5.o/exb.all.o

	marin.o <- marin.o + (probs.region1.o * (1/R))
	eastbay.o <- eastbay.o + (probs.region2.o * (1/R)) 
	bayshore.o <- bayshore.o + (probs.region3.o * (1/R)) 
	peninsula.o <- peninsula.o + (probs.region4.o * (1/R)) 
	notrip.o <- notrip.o + (probs.region5.o * (1/R)) 
	
	# counterfactual -- set alternative specific constant to drive probability of East Bay to zero
	
	beta.R.c <- beta.R
	beta.R.c[,"eastbay"] <- -25

	exb.region1.c <- exp(beta.R.c %*% t(subset(observed.x, id.alt==1, select=c(1:6))))
	exb.region2.c <- exp(beta.R.c %*% t(subset(observed.x, id.alt==2, select=c(1:6))))
	exb.region3.c <- exp(beta.R.c %*% t(subset(observed.x, id.alt==3, select=c(1:6))))
	exb.region4.c <- exp(beta.R.c %*% t(subset(observed.x, id.alt==4, select=c(1:6))))
	exb.region5.c <- exp(beta.R.c %*% t(subset(observed.x, id.alt==5, select=c(1:6))))
	exb.all.c <- exb.region1.c + exb.region2.c + exb.region3.c + exb.region4.c + exb.region5.c
	
	probs.region1.c <- exb.region1.c/exb.all.c
	probs.region2.c <- exb.region2.c/exb.all.c
	probs.region3.c <- exb.region3.c/exb.all.c
	probs.region4.c <- exb.region4.c/exb.all.c
	probs.region5.c <- exb.region5.c/exb.all.c
	
	marin.c <- marin.c + (probs.region1.c * (1/R))
	eastbay.c <- eastbay.c + (probs.region2.c * (1/R)) 
	bayshore.c <- bayshore.c + (probs.region3.c * (1/R)) 
	peninsula.c <- peninsula.c + (probs.region4.c * (1/R)) 
	notrip.c <- notrip.c + (probs.region5.c * (1/R)) 
	

}

# apply sample enumeration

marin.so <- (colSums(apply(marin.o, 1, function(x) {x * id.weights}))/sum(id.weights)) * 100
eastbay.so <- (colSums(apply(eastbay.o, 1, function(x) {x * id.weights}))/sum(id.weights)) * 100
bayshore.so <- (colSums(apply(bayshore.o, 1, function(x) {x * id.weights}))/sum(id.weights)) * 100
peninsula.so <- (colSums(apply(peninsula.o, 1, function(x) {x * id.weights}))/sum(id.weights)) * 100
notrip.so <- (colSums(apply(notrip.o, 1, function(x) {x * id.weights}))/sum(id.weights)) * 100

marin.sc <- (colSums(apply(marin.c, 1, function(x) {x * id.weights}))/sum(id.weights)) * 100
eastbay.sc <- (colSums(apply(eastbay.c, 1, function(x) {x * id.weights}))/sum(id.weights)) * 100
bayshore.sc <- (colSums(apply(bayshore.c, 1, function(x) {x * id.weights}))/sum(id.weights)) * 100
peninsula.sc <- (colSums(apply(peninsula.c, 1, function(x) {x * id.weights}))/sum(id.weights)) * 100
notrip.sc <- (colSums(apply(notrip.c, 1, function(x) {x * id.weights}))/sum(id.weights)) * 100

marin.s.diff <- marin.sc - marin.so
eastbay.s.diff <- eastbay.sc - eastbay.so
bayshore.s.diff <- bayshore.sc - bayshore.so
peninsula.s.diff <- peninsula.sc - peninsula.so
notrip.s.diff <- notrip.sc - notrip.so

shares[,1] <- marin.so
shares[,2] <- eastbay.so
shares[,3] <- bayshore.so
shares[,4] <- peninsula.so
shares[,5] <- notrip.so
shares[,6] <- marin.sc
shares[,7] <- eastbay.sc
shares[,8] <- bayshore.sc
shares[,9] <- peninsula.sc
shares[,10] <- notrip.sc
shares[,11] <- marin.s.diff 
shares[,12] <- eastbay.s.diff 
shares[,13] <- bayshore.s.diff
shares[,14] <- peninsula.s.diff
shares[,15] <- notrip.s.diff
	
quantiles <- apply(shares, 2, function(x) {quantile(x, probs=c(0.5, 0.025, 0.975))})
sds <- apply(shares, 2, sd)
first.differences <- cbind(t(quantiles), sds)
colnames(first.differences) <- c("Median", "2.5%", "97.5%", "(se)")
rownames(first.differences) <- c("Baseline: Marin", "Baseline: East Bay", "Baseline: Bayshore", "Baseline: Peninsula", "Baseline: No Trip", "Counterfactual: Marin", "Counterfactual: East Bay", "Counterfactual: Bayshore", "Counterfactual: Peninsula", "Counterfactual: No Trip", "Difference: Marin", "Difference: East Bay", "Difference: Bayshore", "Difference: Peninsula", "Difference: No Trip")
print(round(first.differences, digits=3))


###############
## Table 6.7 ##
###############

coeffs <- cosco.busan.mxl$coefficients
covmat <- vcov(cosco.busan.mxl)
ndraws <- 1000

betadraw <- mvrnorm(ndraws, coeffs, covmat)
R <- 500
halton.seq <- qnorm(halton(R, dim=6))

## observed

observed.x <- with(cb.data, data.frame(
travelcost = travelcost,
bridge = bridge,
marin = marin,
eastbay = eastbay, 
bayshore = bayshore,
peninsula = peninsula,
id = idx,
weights = wts
))

id.weights <- observed.x$weights[observed.x$id.alt==1]

me.MEM.travelcost <- matrix(0, nrow(betadraw), 1)
cme.MEM.travelcost21 <- matrix(0, nrow(betadraw), 1)
cme.MEM.travelcost23 <- matrix(0, nrow(betadraw), 1)
cme.MEM.travelcost24 <- matrix(0, nrow(betadraw), 1)
cme.MEM.travelcost25 <- matrix(0, nrow(betadraw), 1)

me.AME.travelcost.all <- matrix(0, nrow(betadraw), nrow(observed.x[observed.x$id.alt==4,]))
cme.AME.travelcost21.all <- matrix(0, nrow(betadraw), nrow(observed.x[observed.x$id.alt==1,]))
cme.AME.travelcost23.all <- matrix(0, nrow(betadraw), nrow(observed.x[observed.x$id.alt==2,]))
cme.AME.travelcost24.all <- matrix(0, nrow(betadraw), nrow(observed.x[observed.x$id.alt==3,]))
cme.AME.travelcost25.all <- matrix(0, nrow(betadraw), nrow(observed.x[observed.x$id.alt==5,]))

probs.region1.all <- matrix(0, nrow(betadraw), nrow(observed.x[observed.x$id.alt==1,]))
probs.region2.all <- matrix(0, nrow(betadraw), nrow(observed.x[observed.x$id.alt==2,]))
probs.region3.all <- matrix(0, nrow(betadraw), nrow(observed.x[observed.x$id.alt==3,]))
probs.region4.all <- matrix(0, nrow(betadraw), nrow(observed.x[observed.x$id.alt==4,]))
probs.region5.all <- matrix(0, nrow(betadraw), nrow(observed.x[observed.x$id.alt==5,]))

probs.region1.all.m <- matrix(0, nrow(betadraw), 1)
probs.region2.all.m <- matrix(0, nrow(betadraw), 1)
probs.region3.all.m <- matrix(0, nrow(betadraw), 1)
probs.region4.all.m <- matrix(0, nrow(betadraw), 1)
probs.region5.all.m <- matrix(0, nrow(betadraw), 1)

for (i in 1:R) {

	travelcost <- -1 * exp(betadraw[,"n.travelcost"] + (halton.seq[i,1] * abs(betadraw[,"sd.n.travelcost"])))
	bridge <- betadraw[,"bridge"] + (halton.seq[i,2] * abs(betadraw[,"sd.bridge"]))
	marin <- betadraw[,"marin"] + (halton.seq[i,3] * abs(betadraw[,"sd.marin"]))
	eastbay <- betadraw[,"eastbay"] + (halton.seq[i,4] * abs(betadraw[,"sd.eastbay"]))
	bayshore <- betadraw[,"bayshore"] + (halton.seq[i,5] * abs(betadraw[,"sd.bayshore"]))
	peninsula <- betadraw[,"peninsula"] + (halton.seq[i,6] * abs(betadraw[,"sd.peninsula"]))

	beta.R <- cbind(travelcost, bridge, marin, eastbay, bayshore, peninsula)
	
	# MEM
	
	exb.region1.m <- exp(beta.R %*% colMeans(subset(observed.x, id.alt==1, select=c(1:6))))
	exb.region2.m <- exp(beta.R %*% colMeans(subset(observed.x, id.alt==2, select=c(1:6))))
	exb.region3.m <- exp(beta.R %*% colMeans(subset(observed.x, id.alt==3, select=c(1:6))))
	exb.region4.m <- exp(beta.R %*% colMeans(subset(observed.x, id.alt==4, select=c(1:6))))
	exb.region5.m <- exp(beta.R %*% colMeans(subset(observed.x, id.alt==5, select=c(1:6))))
	exb.all.m <- exb.region1.m + exb.region2.m + exb.region3.m + exb.region4.m + exb.region5.m
	
	probs.region1.m <- exb.region1.m/exb.all.m  
	probs.region2.m <- exb.region2.m/exb.all.m
	probs.region3.m <- exb.region3.m/exb.all.m
	probs.region4.m <- exb.region4.m/exb.all.m
	probs.region5.m <- exb.region5.m/exb.all.m
	
	probs.region1.all.m <- probs.region1.all.m + (probs.region1.m * (1/R))
	probs.region2.all.m <- probs.region2.all.m + (probs.region2.m * (1/R))
	probs.region3.all.m <- probs.region3.all.m + (probs.region3.m * (1/R))
	probs.region4.all.m <- probs.region4.all.m + (probs.region4.m * (1/R))
	probs.region5.all.m <- probs.region5.all.m + (probs.region5.m * (1/R))
	
	## marginal effects for East Bay
	me.MEM.travelcost <- me.MEM.travelcost + (travelcost * (probs.region2.m * (1 - probs.region2.m)) * (1/R)) 

	## cross-marginal effects
	cme.MEM.travelcost21 <- cme.MEM.travelcost21 + ((-1 * travelcost * (probs.region1.m * probs.region2.m)) * (1/R)) 
	cme.MEM.travelcost23 <- cme.MEM.travelcost23 + ((-1 * travelcost * (probs.region3.m * probs.region2.m)) * (1/R))
	cme.MEM.travelcost24 <- cme.MEM.travelcost24 + ((-1 * travelcost * (probs.region4.m * probs.region2.m)) * (1/R))
	cme.MEM.travelcost25 <- cme.MEM.travelcost25 + ((-1 * travelcost * (probs.region5.m * probs.region2.m)) * (1/R))
	
	
	# AME
	
	exb.region1 <- exp(beta.R %*% t(subset(observed.x, id.alt==1, select=c(1:6))))
	exb.region2 <- exp(beta.R %*% t(subset(observed.x, id.alt==2, select=c(1:6))))
	exb.region3 <- exp(beta.R %*% t(subset(observed.x, id.alt==3, select=c(1:6))))
	exb.region4 <- exp(beta.R %*% t(subset(observed.x, id.alt==4, select=c(1:6))))
	exb.region5 <- exp(beta.R %*% t(subset(observed.x, id.alt==5, select=c(1:6))))
	exb.all <- exb.region1 + exb.region2 + exb.region3 + exb.region4 + exb.region5
	
	probs.region1 <- exb.region1/exb.all  
	probs.region2 <- exb.region2/exb.all
	probs.region3 <- exb.region3/exb.all
	probs.region4 <- exb.region4/exb.all
	probs.region5 <- exb.region5/exb.all
	
	probs.region1.all <- probs.region1.all + (probs.region1 * (1/R))
	probs.region2.all <- probs.region2.all + (probs.region2 * (1/R))
	probs.region3.all <- probs.region3.all + (probs.region3 * (1/R))
	probs.region4.all <- probs.region4.all + (probs.region4 * (1/R))
	probs.region5.all <- probs.region5.all + (probs.region5 * (1/R))
	
	## marginal effects for Bayshore
	me.AME.travelcost.all <- me.AME.travelcost.all + ((travelcost * (probs.region2 * (1 - probs.region2))) * (1/R)) 

	## cross-marginal effects
	cme.AME.travelcost21.all <- cme.AME.travelcost21.all + ((-1 * travelcost * (probs.region1 * probs.region2)) * (1/R)) 
	cme.AME.travelcost23.all <- cme.AME.travelcost23.all + ((-1 * travelcost * (probs.region3 * probs.region2)) * (1/R))
	cme.AME.travelcost24.all <- cme.AME.travelcost24.all + ((-1 * travelcost * (probs.region4 * probs.region2)) * (1/R))
	cme.AME.travelcost25.all <- cme.AME.travelcost25.all + ((-1 * travelcost * (probs.region5 * probs.region2)) * (1/R))
	
}

me.AME.travelcost <- colSums(apply(me.AME.travelcost.all, 1, function(x) {x * id.weights}))/sum(id.weights)
cme.AME.travelcost21 <- colSums(apply(cme.AME.travelcost21.all, 1, function(x) {x * id.weights}))/sum(id.weights) 
cme.AME.travelcost23 <- colSums(apply(cme.AME.travelcost23.all, 1, function(x) {x * id.weights}))/sum(id.weights) 
cme.AME.travelcost24 <- colSums(apply(cme.AME.travelcost24.all, 1, function(x) {x * id.weights}))/sum(id.weights) 
cme.AME.travelcost25 <- colSums(apply(cme.AME.travelcost25.all, 1, function(x) {x * id.weights}))/sum(id.weights) 

## Elasticities

travelcost.2 <- observed.x[observed.x$id.alt==2, grep("travelcost", colnames(observed.x))]
travelcost.2.m <- mean(observed.x[observed.x$id.alt==2, grep("travelcost", colnames(observed.x))])

## elasticities for East Bay -- MEM
e.MEM.travelcost <- me.MEM.travelcost * (travelcost.2.m / probs.region2.all.m)


## elasticities for East Bay -- AME
e.travelcost <- me.AME.travelcost.all * t(travelcost.2 / t(probs.region2.all))
e.AME.travelcost <- colSums(apply(e.travelcost, 1, function(x) {x * id.weights}))/sum(id.weights)

## cross-elasticities -- MEM
ce.MEM.travelcost21 <- cme.MEM.travelcost21 * (travelcost.2.m/probs.region1.all.m)
ce.MEM.travelcost23 <- cme.MEM.travelcost23 * (travelcost.2.m/probs.region3.all.m)
ce.MEM.travelcost24 <- cme.MEM.travelcost24 * (travelcost.2.m/probs.region4.all.m)
ce.MEM.travelcost25 <- cme.MEM.travelcost25 * (travelcost.2.m/probs.region5.all.m)

## cross-elasticities -- AME
ce.travelcost21 <- cme.AME.travelcost21.all * t(apply(probs.region1.all, 1, function(x){travelcost.2/x}))
ce.travelcost23 <- cme.AME.travelcost23.all * t(apply(probs.region3.all, 1, function(x){travelcost.2/x}))
ce.travelcost24 <- cme.AME.travelcost24.all * t(apply(probs.region4.all, 1, function(x){travelcost.2/x}))
ce.travelcost25 <- cme.AME.travelcost25.all * t(apply(probs.region5.all, 1, function(x){travelcost.2/x}))

ce.AME.travelcost21 <- colSums(apply(ce.travelcost21, 1, function(x) {x * id.weights}))/sum(id.weights)  
ce.AME.travelcost23 <- colSums(apply(ce.travelcost23, 1, function(x) {x * id.weights}))/sum(id.weights) 
ce.AME.travelcost24 <- colSums(apply(ce.travelcost24, 1, function(x) {x * id.weights}))/sum(id.weights) 
ce.AME.travelcost25 <- colSums(apply(ce.travelcost25, 1, function(x) {x * id.weights}))/sum(id.weights) 

# note marginal effects are multiplied by 100 for ease of interpretation
margins <- cbind(
me.AME.travelcost, me.MEM.travelcost, 
cme.AME.travelcost21, cme.AME.travelcost23, cme.AME.travelcost24, cme.AME.travelcost25, 
cme.MEM.travelcost21, cme.MEM.travelcost23, cme.MEM.travelcost24, cme.MEM.travelcost25
) * 100

elasticities <- cbind(
e.AME.travelcost, e.MEM.travelcost,
ce.AME.travelcost21, ce.AME.travelcost23, ce.AME.travelcost24, ce.AME.travelcost25,
ce.MEM.travelcost21, ce.MEM.travelcost23, ce.MEM.travelcost24, ce.MEM.travelcost25
)

all.results <- cbind(margins, elasticities)

quantiles <- apply(all.results, 2, function(x) {quantile(x, probs=c(0.5, 0.025, 0.975))})
sds <- apply(all.results, 2, sd)
marginal.effects <- cbind(t(quantiles), sds) 
colnames(marginal.effects) <- c("Median", "2.5%", "97.5%", "(se)")
rownames(marginal.effects) <- c("East Bay: AME", "East Bay: MEM", "Marin: Average CME", "Bayshore: Average CME", "Peninsula: Average CME", "No Trip: Average CME", "Marin: CME at Mean", "Bayshore: CME at Mean", "Peninsula: CME at Mean", "No Trip: CME at Mean", "East Bay: Average Elasticity", "East Bay: Elasticity at the Mean", "Marin: Average CE", "Bayshore: Average CE", "Peninsula: Average CE", "No Trip: Average CE", "Marin: CE at Mean", "Bayshore: CE at Mean", "Peninsula: CE at Mean", "No Trip: CE at Mean")
print(round(marginal.effects, digits=3))


## WTP calculation ##

wtp.summary <- matrix(NA, ndraws, 2)

for (i in 1:ndraws) {

	travelcost <- -1 * exp(betadraw[i,"n.travelcost"] + (halton.seq[,1] * abs(betadraw[i,"sd.n.travelcost"])))
	bridge <- betadraw[i,"bridge"] + (halton.seq[,2] * abs(betadraw[i,"sd.bridge"]))
	
	all.bridge <- expand.grid(bridge, travelcost)
	wtp.bridge <- all.bridge[,1]/(-1 * all.bridge[,2])
	wtp.summary[i,1] <- median(wtp.bridge)
	wtp.summary[i,2] <- mean(wtp.bridge)
	
}

quantiles <- apply(wtp.summary, 2, function(x) {quantile(x, probs=c(0.5, 0.025, 0.975))})
sds <- apply(wtp.summary, 2, sd)
wtp.all <- cbind(t(quantiles), sds) 
colnames(wtp.all ) <- c("Median", "2.5%", "97.5%", "(se)")
rownames(wtp.all ) <- c("Median WTP", "Mean WTP")
print(round(wtp.all , digits=3))



###############
## Table 6.8 ##
###############

coeffs <- cosco.busan.mxl$coefficients
covmat <- vcov(cosco.busan.mxl)
ndraws <- 1000

betadraw <- mvrnorm(ndraws, coeffs, covmat)
intdraws <- 500
halton.seq <- qnorm(halton(intdraws, dim=6))

## observed

observed.x <- with(cb.data, data.frame(
travelcost = travelcost,
bridge = bridge,
marin = marin,
eastbay = eastbay, 
bayshore = bayshore,
peninsula = peninsula,
id = idx,
weights = wts
))

p.01 <- matrix(0, nrow(betadraw), nrow(observed.x))
p.02 <- matrix(0, nrow(betadraw), nrow(observed.x))
p.03 <- matrix(0, nrow(betadraw), nrow(observed.x))
p.04 <- matrix(0, nrow(betadraw), nrow(observed.x))
p.05 <- matrix(0, nrow(betadraw), nrow(observed.x))

p.tr1 <- matrix(0, nrow(betadraw), nrow(observed.x))
p.tr2 <- matrix(0, nrow(betadraw), nrow(observed.x))
p.tr3 <- matrix(0, nrow(betadraw), nrow(observed.x))
p.tr4 <- matrix(0, nrow(betadraw), nrow(observed.x))
p.tr5 <- matrix(0, nrow(betadraw), nrow(observed.x))

p.b1 <- matrix(0, nrow(betadraw), nrow(observed.x))
p.b2 <- matrix(0, nrow(betadraw), nrow(observed.x))
p.b3 <- matrix(0, nrow(betadraw), nrow(observed.x))
p.b4 <- matrix(0, nrow(betadraw), nrow(observed.x))
p.b5 <- matrix(0, nrow(betadraw), nrow(observed.x))

# travelcost
observed.x.tr <- observed.x
observed.x.tr[,"travelcost"][observed.x.tr[,"id.alt"] != 5] <- observed.x.tr[,"travelcost"][observed.x.tr[,"id.alt"] != 5] + 1

# bridge
observed.x.b <- observed.x
observed.x.b[,"bridge"][observed.x.b["id.alt"] != 5] <- observed.x.b["bridge"][observed.x.b["id.alt"] != 5]  + 1

for (i in 1:ndraws){

	travelcost <- -1 * exp(betadraw[i,"n.travelcost"] + (halton.seq[,1] * abs(betadraw[i,"sd.n.travelcost"])))
	bridge <- betadraw[i,"bridge"] + (halton.seq[,2] * abs(betadraw[i,"sd.bridge"]))
	marin <- betadraw[i,"marin"] + (halton.seq[,3] * abs(betadraw[i,"sd.marin"]))
	eastbay <- betadraw[i,"eastbay"] + (halton.seq[,4] * abs(betadraw[i,"sd.eastbay"]))
	bayshore <- betadraw[i,"bayshore"] + (halton.seq[,5] * abs(betadraw[i,"sd.bayshore"]))
	peninsula <- betadraw[i,"peninsula"] + (halton.seq[,6] * abs(betadraw[i,"sd.peninsula"]))

	beta.R <- cbind(travelcost, bridge, marin, eastbay, bayshore, peninsula)
	
	exb.region1.0 <- exp(beta.R %*% t(subset(observed.x, id.alt==1, select=c(1:6))))
	exb.region2.0 <- exp(beta.R %*% t(subset(observed.x, id.alt==2, select=c(1:6))))
	exb.region3.0 <- exp(beta.R %*% t(subset(observed.x, id.alt==3, select=c(1:6))))
	exb.region4.0 <- exp(beta.R %*% t(subset(observed.x, id.alt==4, select=c(1:6))))
	exb.region5.0 <- exp(beta.R %*% t(subset(observed.x, id.alt==5, select=c(1:6))))
	exb.all.0 <- exb.region1.0 + exb.region2.0 + exb.region3.0 + exb.region4.0 + exb.region5.0
	
	p.01[i,] <- colMeans(exb.region1.0/exb.all.0)
	p.02[i,] <- colMeans(exb.region2.0/exb.all.0)
	p.03[i,] <- colMeans(exb.region3.0/exb.all.0)
	p.04[i,] <- colMeans(exb.region4.0/exb.all.0)
	p.05[i,] <- colMeans(exb.region5.0/exb.all.0)

	exb.region1.tr <- exp(beta.R %*% t(subset(observed.x.tr, id.alt==1, select=c(1:6))))
	exb.region2.tr <- exp(beta.R %*% t(subset(observed.x.tr, id.alt==2, select=c(1:6))))
	exb.region3.tr <- exp(beta.R %*% t(subset(observed.x.tr, id.alt==3, select=c(1:6))))
	exb.region4.tr <- exp(beta.R %*% t(subset(observed.x.tr, id.alt==4, select=c(1:6))))
	exb.region5.tr <- exp(beta.R %*% t(subset(observed.x.tr, id.alt==5, select=c(1:6))))
	exb.all.tr <- exb.region1.tr + exb.region2.tr + exb.region3.tr + exb.region4.tr + exb.region5.tr
	
	p.tr1[i,] <- colMeans(exb.region1.tr/exb.all.tr)
	p.tr2[i,] <- colMeans(exb.region2.tr/exb.all.tr)
	p.tr3[i,] <- colMeans(exb.region3.tr/exb.all.tr)
	p.tr4[i,] <- colMeans(exb.region4.tr/exb.all.tr)
	p.tr5[i,] <- colMeans(exb.region5.tr/exb.all.tr)
	
	exb.region1.b <- exp(beta.R %*% t(subset(observed.x.b, id.alt==1, select=c(1:6))))
	exb.region2.b <- exp(beta.R %*% t(subset(observed.x.b, id.alt==2, select=c(1:6))))
	exb.region3.b <- exp(beta.R %*% t(subset(observed.x.b, id.alt==3, select=c(1:6))))
	exb.region4.b <- exp(beta.R %*% t(subset(observed.x.b, id.alt==4, select=c(1:6))))
	exb.region5.b <- exp(beta.R %*% t(subset(observed.x.b, id.alt==5, select=c(1:6))))
	exb.all.b <- exb.region1.b + exb.region2.b + exb.region3.b + exb.region4.b + exb.region5.b
	
	p.b1[i,] <- colMeans(exb.region1.b/exb.all.b)
	p.b2[i,] <- colMeans(exb.region2.b/exb.all.b)
	p.b3[i,] <- colMeans(exb.region3.b/exb.all.b)
	p.b4[i,] <- colMeans(exb.region4.b/exb.all.b)
	p.b5[i,] <- colMeans(exb.region5.b/exb.all.b)
	
	
}

id.weights <- cb.data$wts

t1 <- (p.tr1 / p.tr5) / (p.01 / p.05)
t2 <- (p.tr2 / p.tr5) / (p.02 / p.05)
t3 <- (p.tr3 / p.tr5) / (p.03 / p.05)
t4 <- (p.tr4 / p.tr5) / (p.04 / p.05)

b1 <- (p.b1 / p.b5) / (p.01 / p.05)
b2 <- (p.b2 / p.b5) / (p.02 / p.05)
b3 <- (p.b3 / p.b5) / (p.03 / p.05)
b4 <- (p.b4 / p.b5) / (p.04 / p.05)

OR.tr.1.w <- colSums(apply(t1, 1, function(x) {x * id.weights}))/sum(id.weights)
OR.tr.2.w <- colSums(apply(t2, 1, function(x) {x * id.weights}))/sum(id.weights)
OR.tr.3.w <- colSums(apply(t3, 1, function(x) {x * id.weights}))/sum(id.weights)
OR.tr.4.w <- colSums(apply(t4, 1, function(x) {x * id.weights}))/sum(id.weights)

OR.b.1.w <- colSums(apply(b1, 1, function(x) {x * id.weights}))/sum(id.weights)
OR.b.2.w <- colSums(apply(b2, 1, function(x) {x * id.weights}))/sum(id.weights)
OR.b.3.w <- colSums(apply(b3, 1, function(x) {x * id.weights}))/sum(id.weights)
OR.b.4.w <- colSums(apply(b4, 1, function(x) {x * id.weights}))/sum(id.weights)

OR <- cbind(OR.tr.1.w, OR.tr.2.w, OR.tr.3.w, OR.tr.4.w, OR.b.1.w, OR.b.2.w, OR.b.3.w, OR.b.4.w)
quantiles <- apply(OR, 2, function(x) {quantile(x, probs=c(0.5, 0.025, 0.975))})
sds <- apply(OR, 2, sd)
pop.odds <- cbind(t(quantiles), sds)
colnames(pop.odds) <- c("Median", "2.5%", "97.5%", "(se)")
rownames(pop.odds) <- c("Travel Cost: Marin", "Travel Cost: East Bay", "Travel Cost: Bayshore", "Travel Cost: Peninsula", "Bridge Crossing: Marin", "Bridge Crossing: East Bay", "Bridge Crossing: Bayshore", "Bridge Crossing: Peninsula")
print(round(pop.odds, digits=3))


