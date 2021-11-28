###############
## Section 2 ##
###############

library(MASS)  

setwd("C:\\dc_book\\Replication_Files_Element")

ITN.data <- read.csv("ITN_data.csv", header=TRUE)

###############
## Table 2.1 ##
###############

## Logit model

ITN.logit <- glm(useditn ~  malariamessage+education+richest+poorest+pregnant+altitude, family=binomial(link="logit"), data=ITN.data, x=TRUE)
ITN.logit.results <- cbind(ITN.logit$coefficients, sqrt(diag(vcov(ITN.logit))), confint.default(ITN.logit, level=0.95))
N <- length(ITN.logit$fitted.values)
colnames(ITN.logit.results) <- c("Coeff", "(se)", "2.5%", "97.5%")
print((round(ITN.logit.results, digits=3)))
print(N)

###############
## Table 2.2 ##
###############

## Simulated Coefficients

beta <- ITN.logit$coefficients
covmat.beta <- vcov(ITN.logit)
ndraws <- 1000
betadraw <- mvrnorm(ndraws, beta, covmat.beta)

mean.x <- colMeans(ITN.logit$x)
simprobs <- plogis(betadraw%*%mean.x)

simcoef.results <- append(quantile(simprobs, probs=c(0.5, 0.025, 0.975)), sd(simprobs))

## Bootstrapping

nsamples <- 1000
betas.boot <- matrix(NA,nsamples,length(beta))

for (i in 1:1000) {

	ITN.data.boot <- ITN.data[sample(seq(1,nrow(ITN.data),1), replace=TRUE),]
	logit.boot <- glm(useditn ~  malariamessage+education+richest+poorest+pregnant+altitude, family=binomial(link="logit"), data=ITN.data.boot)
	betas.boot[i,] <- logit.boot$coefficients	
}

probs.boot <- plogis(betas.boot%*%mean.x)

boot.results <- append(quantile(probs.boot, probs=c(0.5, 0.025, 0.975)), sd(probs.boot))

## Delta Method

xb <- mean.x %*% beta
predicted.prob <- plogis(xb)

deltamethod.var <- (dlogis(xb)^2) * (t(mean.x) %*% covmat.beta %*% mean.x)
deltamethod.se <- sqrt(deltamethod.var)

delta.95ci.u <- predicted.prob + 1.96*deltamethod.se
delta.95ci.l <- predicted.prob - 1.96*deltamethod.se

delta.results <- cbind(predicted.prob, delta.95ci.l, delta.95ci.u, deltamethod.se)

## Summarize Results

statistical.uncertainty <- rbind(simcoef.results, boot.results, delta.results)
colnames(statistical.uncertainty) <- c("Median", "2.5%", "97.5%", "(se)")
rownames(statistical.uncertainty) <- c("Simulated Coefficients", "Bootstrap", "Delta Method")

print(round(statistical.uncertainty, digits=3))

