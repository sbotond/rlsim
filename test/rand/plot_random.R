#!/usr/bin/env Rscript

n <- 2 * 10^5

pdf("Poisson.pdf")
# Poisson with mean 10:
data <- read.table("poisson_10.tab")$V1
true <-rpois(n,10)
qqplot(data, true)

# Poisson with mean 100:
data <- read.table("poisson_100.tab")$V1
true <-rpois(n,100)
qqplot(data, true)

pdf("Binomial.pdf")
# Binomial n=20, p = 0.3
data <- read.table("binomial_20_0.3.tab")$V1
true <-rbinom(n,20,0.3)
qqplot(data, true)

# Binomial n=1000, p = 0.75
data <- read.table("binomial_1000_0.75.tab")$V1
true <-rbinom(n,1000,0.75)
qqplot(data, true)
