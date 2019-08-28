test_cor <- function(trait.cor, cor.nperm = 10000) {
  tmp <- svd(trait.cor)

  eigenvalue <- tmp$d
  eigenvalue[eigenvalue < max(eigenvalue) / 1] <- 0
  eigenvalue[eigenvalue != 0] <- 1 / eigenvalue[eigenvalue != 0]
  multi.df1 <- sum(eigenvalue != 0)
  D <- diag(eigenvalue)
  t.mat1 <- tmp$u %*% D %*% t(tmp$v)


  eigenvalue <- tmp$d
  eigenvalue[eigenvalue < max(eigenvalue) / 10] <- 0
  eigenvalue[eigenvalue != 0] <- 1 / eigenvalue[eigenvalue != 0]
  multi.df10 <- sum(eigenvalue != 0)
  D <- diag(eigenvalue)
  t.mat10 <- tmp$u %*% D %*% t(tmp$v)


  eigenvalue <- tmp$d
  eigenvalue[eigenvalue < max(eigenvalue) / 30] <- 0
  eigenvalue[eigenvalue != 0] <- 1 / eigenvalue[eigenvalue != 0]
  multi.df30 <- sum(eigenvalue != 0)
  D <- diag(eigenvalue)
  t.mat30 <- tmp$u %*% D %*% t(tmp$v)


  eigenvalue <- tmp$d
  eigenvalue[eigenvalue < max(eigenvalue) / 50] <- 0
  eigenvalue[eigenvalue != 0] <- 1 / eigenvalue[eigenvalue != 0]
  multi.df50 <- sum(eigenvalue != 0)
  D <- diag(eigenvalue)
  t.mat50 <- tmp$u %*% D %*% t(tmp$v)

  tmp <- svd(trait.cor)
  eigenvalue <- tmp$d
  # eigenvalue[eigenvalue < max(eigenvalue) / 10000] <- 0
  eigenvalue[eigenvalue < 0] <- 0
  D <- diag(eigenvalue)
  CovSsqrt <- tmp$u %*% sqrt(D) %*% t(tmp$v)

  T1s <- matrix(NA, cor.nperm, 4)
  start.time <- proc.time()[3]

  for (i in 1:cor.nperm) {
    tmp.z <- rnorm(dim(CovSsqrt)[1], 0, 1) %*% CovSsqrt
    tmp.z <- as.vector(tmp.z)

    pAT1 <- MAT_fast(tmp.z, t.mat1, multi.df1)$multi.p
    pAT10 <- MAT_fast(tmp.z, t.mat10, multi.df10)$multi.p
    pAT30 <- MAT_fast(tmp.z, t.mat30, multi.df30)$multi.p
    pAT50 <- MAT_fast(tmp.z, t.mat50, multi.df50)$multi.p

    T1s[i,] <- c(pAT1, pAT10, pAT30, pAT50)
  }

  colSums(T1s < 0.01) / cor.nperm
  T1s.save <- T1s
  T1s[T1s > 0.99] <- 0.99
  T1s <- qnorm(1 - T1s)
  setbased_corEst1 <- cor(T1s[1:cor.nperm,])

  return(setbased_corEst1)
}



MAT <- function(Z.score, trait.cor, cutoff = 30) {
  Z.score <- as.numeric(Z.score)

  tmp <- svd(trait.cor)
  eigenvalue <- tmp$d
  eigenvalue[eigenvalue < max(eigenvalue) / cutoff] <- 0
  eigenvalue[eigenvalue != 0] <- 1 / eigenvalue[eigenvalue != 0]

  multi.df <- sum(eigenvalue != 0)
  D <- diag(eigenvalue)
  t.mat2 <- tmp$u %*% D %*% t(tmp$v)

  multi.stat <- Z.score %*% t.mat2 %*% Z.score
  multi.p <- pchisq(multi.stat, df = multi.df, lower.tail = F)


  return(list(multi.stat = multi.stat, multi.p = multi.p, multi.df = multi.df))
}

aMAT_fast <- function(z.tmp, t.mat1, t.mat10, t.mat30, t.mat50, multi.df1, multi.df10, multi.df30, multi.df50, setbased_corEst1) {
  pAT1 <- MAT_fast(tmp.z, t.mat1, multi.df1)$multi.p
  pAT10 <- MAT_fast(tmp.z, t.mat10, multi.df10)$multi.p
  pAT30 <- MAT_fast(tmp.z, t.mat30, multi.df30)$multi.p
  pAT50 <- MAT_fast(tmp.z, t.mat50, multi.df50)$multi.p

  omni_stat <- min(pAT1, pAT10, pAT30, pAT50)

  omni_p <- 1 - mvtnorm::pmvnorm(lower = -Inf, upper = rep(qnorm(1 - omni_stat), 4), sigma = setbased_corEst1)[1]
  return(list(MAT1 = pAT1, MAT10 = pAT10, MAT30 = pAT30, MAT50 = pAT50, aMAT = omni_p))
}


MAT_fast <- function(Z.score, t.mat2, multi.df) {
  # tmp is  svd(trait.cor)
  multi.stat <- Z.score %*% t.mat2 %*% Z.score
  multi.p <- pchisq(multi.stat, df = multi.df, lower.tail = F)

  return(list(multi.stat = multi.stat, multi.p = multi.p, multi.df = multi.df))
}