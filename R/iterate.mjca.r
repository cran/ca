################################################################################
iterate.mjca <- function(B, lev.n, nd = 2, maxit = 50, epsilon = 0.0001) {
  if (is.na(maxit) & is.na(epsilon)){
    maxit   <- 50
    epsilon <- 0.0001
    }

  coord <- NA
  lev   <- lev.n
  n     <- sum(B)
  J     <- sum(lev)
  Q     <- length(lev)
  foo0  <- rep(1:Q, lev.n)
  foo1  <- (foo0) %*% t(rep(1, sum(lev.n))) - t((foo0) %*% t(rep(1, sum(lev.n))))
  dummy <- ifelse(foo1 == 0, 1, 0) 
 # li    <- as.vector(c(0,cumsum(lev)))
 # dummy <- matrix(0, J, J)
 # for (i in 1:(length(li)-1)) {
 #   ind.lo <- li[i]+1
 #   ind.up <- li[i+1]
 #   ind.to <- diff(li)[i]
 #   dummy[rep(ind.lo:ind.up, ind.to) +
 #        (rep(ind.lo:ind.up, each = ind.to)-1) * J] <- 1
 #   }
  iterate <- function(obj, dummy, nd, adj = FALSE) {
    Bp     <- obj/n
    cm     <- apply(Bp, 2, sum)
    eP     <- cm %*% t(cm)
    cm.mat <- diag(cm^(-0.5))
    S      <- cm.mat %*% (Bp - eP) %*% cm.mat
    dec    <- eigen(S)
    lam    <- dec$values
    u      <- dec$vectors
    phi    <- u[, 1:nd] / matrix(rep(sqrt(cm), nd), ncol = nd)
    if (adj){
      lam <- (Q / (Q - 1))^2 * (lam[lam >= 1 / Q] - 1 / Q) ^ 2
      }
    for (s in 1:nd) {
#      if (exists("coord")) {
      if (!is.na(coord[1])) {
        coord <- coord + lam[s] * (phi[,s] %*% t(phi[,s]))
        } else {
        coord <- lam[s] * (phi[,s] %*% t(phi[,s]))
        }
      }
    return(obj * (1 - dummy) + n * eP * dummy * (1 + coord))
    }
 # first iteration (adjusted lambda)
  B.star <- iterate(B, dummy, nd, adj = TRUE)
 # subsequent iterations
  k  <- 1
  it <- TRUE
  while (it) {
    temp    <- iterate(B.star, dummy, nd)
    delta.B <- max(abs(B.star - temp))
    B.star  <- temp
    if (delta.B <= epsilon | k >= maxit) it <- FALSE
    k       <- k + 1
    }
  return(list(B.star, crit=c(k-1, delta.B)))
  }
