# chr1	91239355	0.223904	120.452070 100001
# chr1	92088305	0.615182	121.275592 101000
p <- 2500
sample_info <- read.delim("sample_info")
sample_lis <- list()
for (pop in c("AFR", "AMR", "EAS", "EUR", "SAS")) {
  sample_lis[[pop]] = sample_info$sample[which(sample_info$super_pop==pop)]
}
# summary(sample_lis)
# Length Class  Mode
# AFR 661    -none- character
# AMR 347    -none- character
# EAS 504    -none- character
# EUR 503    -none- character
# SAS 489    -none- character
vcf_file <- file.path("test_SMME_phase3_chr1.vcf.gz")
vcf<- sim1000G::readVCF(vcf_file, maxNumberOfVariants = p, 
                        min_maf = 0.02 ,max_maf = NA)
# [##......] Chromosome:   1  Mbp:  91.23936  Region Size:  848.95 kb  Num of individuals: 2504
# [##......] Before filtering  Num of variants: 23501 Num of individuals: 2504
# [###.....] After filtering  Num of variants: 2500 Num of individuals: 2504
idx = list()
for (spl in vcf$individual_ids) {
  for (pop in c("AFR", "AMR", "EAS", "EUR", "SAS")) {
    if (spl %in% sample_lis[[pop]]) {
      idx[[pop]] = c(idx[[pop]], which(vcf$individual_ids==spl))
    }
  }
}
geno <- list()
for (pop in c("AFR", "AMR", "EAS", "EUR", "SAS")) {
  geno[[pop]] = vcf$gt1[, idx[[pop]]] + vcf$gt2[, idx[[pop]]]
}
# dim(geno)
# [1] 322 267
remove(vcf)
ssize = c()
for (pop in c("AFR", "AMR", "EAS", "EUR", "SAS")) {
  cat(pop, length(colnames(geno[[pop]])), "\n")
  ssize = c(ssize, length(colnames(geno[[pop]])))
}

# install.packages("SMME")
library("SMME")
##size of example
set.seed(42)
G <- length(ssize); n <- ssize
x <- y <- list()

##group design matrices
# for(g in 1:G){x[[g]] <- matrix(rnorm(n[g] * p), n[g], p)}
for(g in 1:G){
  x[[g]] <- t(geno[[c("AFR", "AMR", "EAS", "EUR", "SAS")[g]]])
}
remove(geno)

##common features and effects
common_features <- rbinom(p, 1, 0.1)
common_effects <- rnorm(p) * common_features

##group response
for(g in 1:G){
  bg <- rnorm(p, 0, 0.5) * (1 - common_features) + common_effects
  mu <- x[[g]] %*% bg
  y[[g]] <- rnorm(n[g]) + mu
}

##fit model for range of lambda and zeta
system.time(fit <- softmaximin(x, y, zeta = c(0.1, 1), penalty = "lasso", alg = "npg"))
betahat <- fit$coef
# minimum number of backtraking steps reached for model no. 17  for zeta(s) 0.1
# Multithreading enabled using 4 threads
# When p = 5000, n = 2500, the running time is as:
# user   system   elapse
# 677.53  617.48 1021.45

##estimated common effects for specific lambda and zeta
modelno <- 6; zetano <- 2
m <- min(betahat[[zetano]][ , modelno], common_effects)
M <- max(betahat[[zetano]][ , modelno], common_effects)
plot(common_effects, type = "p", ylim = c(m, M), col = "red")
lines(betahat[[zetano]][ , modelno], type = "h")

system.time(fit <- softmaximin(x, y, zeta = c(0.1, 1), penalty = "lasso", alg = "fista"))
# maximum number of backtraking steps reached for model no. 17  for zeta(s) 0.1
# Multithreading enabled using 4 threads
# When p = 5000, n = 2500, the running time is as:
# user  system  elapse
# 661.32 593.00 983.05
## evaluation
# input: pred = Xb+intercept, ytest, family = 'binomial' / 'gaussian'
eval <- function(pred, ytest, family){
  if (family == 'binomial') {
    ## type.measure="class"
    y = as.factor(ytest)
    ntab = table(y)
    nc = as.integer(length(ntab))
    y = diag(nc)[as.numeric(y), , drop=FALSE]
    predmat=1 / (1+exp(-pred))
    class = y[, 1] * (predmat > 0.5) + y[, 2] * (predmat <= 0.5)
    return(sum(class)/length(class))
  }
  if (family == 'gaussian') {
    ## type.measure="mse"
    mse = mean((ytest - pred)^2)
    rsq = (cor(ytest, pred))^2
    return(list(mse = mse, nrsq = -rsq))
  }
}

cross_validation <- function(x, y, betahat){
  metric = list()
  for (z_idx in 1:length(betahat)) {
    temp = c()
    met_all = c()
    lam_idx = c()
    for (lam in 1:dim(betahat[[z_idx]])[2]) {
      vg = c()
      for (i in 1:length(y)) {
        beta <- betahat[[z_idx]][, lam]
        pred <- x[[i]] %*% beta
        vg = c(vg, eval(pred, y[[i]], family = "gaussian"))$nrsq
      }
      vg = ifelse(is.na(vg), 0, vg)
      vg = min(vg, na.rm = T)
      temp = c(temp, vg)
    }
    met_all = c(met_all, min(temp, na.rm = T))
    lam_idx = c(lam_idx, which(temp==met_all)[1])
    metric[[z_idx]] = temp
  }
  z_idx =  which(met_all == min(met_all, na.rm = T))
  mse = metric[[z_idx]][lam_idx]
  return(list(z_idx = z_idx, lam_idx = lam_idx, metric = metric, mse = mse))
}
res = cross_validation(x, y, betahat)
z = fit$zeta[res$z_idx]
lamda_cv = fit$lambda[[res$z_idx]][res$lam_idx]
beta = fit$coef[[res$z_idx]][, res$lam_idx]
vg = c()
for (i in 1:length(y)) {
  pred <- x[[i]] %*% beta
  vg = c(vg, (cor(pred, y[[i]])^2))
  vg = ifelse(is.na(vg), 0, vg)
}
vg
for (dx in 3:8) {
  beta = fit$coef[[res$z_idx]][, dx]
  vg = c()
  for (i in 1:length(y)) {
    pred <- x[[i]] %*% beta
    vg = c(vg, (cor(pred, y[[i]])^2))
    vg = ifelse(is.na(vg), 0, vg)
  }
  print(vg)
}

