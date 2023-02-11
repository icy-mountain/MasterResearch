rm(list = ls(all.names = TRUE))
# install.packages('Rsolnp')
# install.packages('mosaicCalc')
library(Rsolnp)
library(mosaicCalc)
source("./tables.R")
source("./utilities.R")
source("./ASkfModel.R")
# params$f <- "(1 - t)^2" # χ2-divergence(Pearsonian distance).
# params$f <- "-log(t)" # reverse KL divergence.
# params$f <- "t * log(t)" # KLdivergence.
# params$f <- "(t - 1) * log(t)" # power-divergence??
# params$f <- "(3*(3 + 1))^(-1) * (t^(4) - t)" # power-divergence lambda = 3

# for screen shot ####
freq <-
  c(98,   150,  135,   53,
    37,   131,  133,   43,
    9,    16,   33,   15,
    4,    1,    4,    21)
ASkf_result <- DisplayASkfResult(freq=freq, f="t * log(t)", name="t", score=1:4, k=1)
ASkf_result <- DisplayASkfResult(freq=freq, f="(1 - t)^2", name="t", score=1:4, k=2)
ASkf_result <- DisplayASkfResult(freq=freq, f="-log(t)", name="t", score=1:4, k=3)
ASkf_result <- DisplayASkfResult(freq=freq, f="(1 - sqrt(x))^2", name="x", score=1:4, k=3)
ASkf_result <- DisplayASkfResult(freq=freq, f="1/x - 1", name="x", score=1:4, k=3)
ASkf_result <- DisplayASkfResult(freq=freq, f="-(x+1) * log((x+1)/2) + x*log(x)", name="x", score=1:4, k=3)

ASkf_result <- DisplayASkfResult(freq=freq, f="x^3", name="x", score=1:4, k=3)
ASkf_result <- DisplayASkfResult(freq=freq, f="(1 - x)^4", name="x", score=1:4, k=3)

# time speed ####
freq <- time1_freq
r <- CountRow(freq)
system.time(DisplayASkfResult(freq=freq, f="x * log(x)", name="x", score=1:r, k=1))


ASkf_result <- DisplayASkfResult(freq=freq, f="x * log(x)", name="x", score=1:4, k=1)
ASkf_result$pars
ASkf_result <- DisplayASkfResult(freq=freq, f="(1 - x)^2", name="x", score=1:4, k=1)
ASkf_result <- DisplayASkfResult(freq=freq, f="-log(x)", name="x", score=1:4, k=2)
ASkf_result <- DisplayASkfResult(freq=freq, f="(1 - sqrt(x))^2", name="x", score=1:4, k=3)
ASkf_result <- DisplayASkfResult(freq=freq, f="1/x - 1", name="x", score=1:4, k=3)
ASkf_result <- DisplayASkfResult(freq=freq, f="-(x+1) * log((x+1)/2) + x*log(x)", name="x", score=1:4, k=3)

ASkf_result <- DisplayASkfResult(freq=freq, f="x^3", name="x", score=1:4, k=3)
ASkf_result <- DisplayASkfResult(freq=freq, f="(1 - x)^4", name="x", score=1:4, k=3)