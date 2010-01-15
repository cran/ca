mjca <- function(obj, nd = 2, lambda = "adjusted", supcol = NA, subsetcol = NA, ps = "", maxit = 50, epsilon = 0.0001) {

 ### SOME CHECKS:
 ## subsets only for Burt-case:
  if ((!(length(unlist(subsetcol)) == 1 & is.na(unlist(subsetcol)[1]))) & lambda != "Burt"){
    lambda <- "Burt"
    warning("Subset analysis is only supported for lambda=\"Burt\" \n lambda is set to \"Burt\"")
    }
 ## valid argument for 'lambda'
  lam.v0 <- c("indicator", "Burt", "adjusted", "JCA")
  lam.v1 <- c("i", "B", "a", "J")
 # is lambda in lam.v0?
  if (length(grep(tolower(lambda), tolower(lam.v0), fixed = TRUE)) != 0){
    lambda <- lam.v0[grep(tolower(lambda), tolower(lam.v0), fixed = TRUE)[1]]
    } else {
   # only first letter specified?
    if (length(grep(tolower(lambda), lam.v1)) == 1){
      lambda <- lam.v0[grep(tolower(lambda), lam.v1)]
      } else {
     # "fuzzy matching" unique?
      if (length(agrep(tolower(lambda), tolower(lam.v0))) == 1 ) {
        lambda <- lam.v0[agrep(tolower(lambda), tolower(lam.v0))]
        } else {
        stop(paste("\nInvalid 'lambda' specification. Valid values are:\n", 
             paste("\"", lam.v0, "\" ", c(",", ",", "and ", "."), sep = "", collapse = ""), 
                   collapse = "", sep = ""))
        }
      }
    }


#  if(!is.data.frame(obj)){
    obj <- data.frame(lapply(data.frame(obj), as.factor)) 
#    }
  I        <- dim(obj)[1]
  cases    <- dim(obj)[2]
  supcol0  <- supcol

# Subset check:
  subsetcol.in <- NA
  if (!(is.na(subsetcol)[1] & length(subsetcol) == 1)){
   # check if given as vector
    if (mode(subsetcol) != "list"){
### bcn, 2009_11: negative indexes:
      foo01 <- unlist(lapply(obj, nlevels))
      if (sum(subsetcol < 0) == length(subsetcol)){
        subsetcol <- (1:sum(foo01))[subsetcol]
        }
      lut  <- cumsum(foo01) - unlist(foo01)
      s0   <- (1:sum(foo01))[subsetcol]
### bcn, 2009_11: ??? not neccessary below:
    #  s0 <- list()
    #  for (i in 2:length(obj)){
    #    foo00     <- subsetcol[subsetcol > lut[i-1] & subsetcol <= lut[i]] - lut[i-1]
    #    if (length(foo00) > 0) { 
    #      s0[[i-1]] <- foo00 
    #      } else {
    #      s0[[i-1]] <-NA
    #      }
    #    }
      } # end subset-vector

    if (mode(subsetcol) == "list"){
      s0 <- list()
      if (length(subsetcol) < length(obj)){
        for (i in (length(subsetcol)+1):length(obj)){
          subsetcol[[i]] <- NA
          }
        }
      for (i in 1:length(obj)){
        if (is.na(subsetcol[[i]])[1]){
          s0[[i]] <- NA
          } else {
          s0[[i]] <- (1:nlevels(obj[[i]]))[subsetcol[[i]]]
          }
        }
      } 
    subsetcol <- s0
    } 

# supplementary columns (variables)
  if (!is.na(supcol)[1]){
    obj.supcol   <- data.frame(obj[,supcol])
    colnames(obj.supcol) <- colnames(obj)[supcol]
    obj          <- data.frame(obj[,-supcol])
    Q.star       <- dim(obj.supcol)[2]
    levels.n.sup <- unlist(lapply(obj.supcol, nlevels))
    n.star       <- cumsum(levels.n.sup)
    J.star       <- sum(levels.n.sup)
    supcol       <- 1:J.star + sum(unlist(lapply(obj, nlevels)))
    Z.star       <- matrix(0, nrow = I, ncol = n.star[length(n.star)])
    newdat.star  <- lapply(obj.supcol, as.numeric)
    offset.star  <- (c(0, n.star[-length(n.star)]))
    for (i in 1:Q.star) 
      Z.star[1:I + (I * (offset.star[i] + newdat.star[[i]] - 1))] <- 1
    fn.star      <- rep(names(obj.supcol), unlist(lapply(obj.supcol, nlevels)))
    ln.star      <- unlist(lapply(obj.supcol,levels))
    cn.star      <- dimnames(obj.supcol)[2]
   # check: subset and sup?
    if (!(is.na(subsetcol)[1] & length(subsetcol) == 1)){
      for (i in sort(supcol0, decr = TRUE)){
        subsetcol <- subsetcol[-i]
        }
      }
    }

  levels.n <- unlist(lapply(obj, nlevels))
  rn       <- dimnames(obj)[[1]]
  cn       <- dimnames(obj)[[2]]


 # prepare data
  n        <- cumsum(levels.n)
  Q        <- dim(obj)[2]
 # get column indexes for subsets:
### BCN 2009_11: ??? below:
 # if (!(is.na(subsetcol)[1] & length(subsetcol) == 1)){
 #   subsetcol.in <- numeric(0)
 #   lut  <- cumsum(unlist(lapply(obj, nlevels))) - unlist(lapply(obj, nlevels))
 #   for (i in 1:length(subsetcol)){
 #     if (!is.na(subsetcol[[i]][1])){
 #       subsetcol.in <- c(subsetcol.in, subsetcol[[i]] + lut[i])
 #       }
 #     }
 #   }
  subsetcol.in <- subsetcol
  
 # indicator and burt matrix:
  Z        <- matrix(0, nrow = I, ncol = n[length(n)])
  newdat   <- lapply(obj, as.numeric)
  offset   <- (c(0, n[-length(n)]))
  for (i in 1:Q) 
    Z[1:I + (I * (offset[i] + newdat[[i]] - 1))] <- 1
  fn       <- rep(names(obj), unlist(lapply(obj, nlevels)))
  ln       <- unlist(lapply(obj,levels))
  B        <- t(Z)%*%Z
  J        <- dim(B)[1]

  col.names        <- paste(fn, ln, sep = ps)
  dimnames(Z)[[2]] <- col.names
  dimnames(Z)[[1]] <- as.character(1:I)
  dimnames(B)[[2]] <- col.names
  dimnames(B)[[1]] <- col.names

 # some placeholders for adjusted/JCA
  B.star       <- NA
  lambda.adj   <- NA
  JCA.it       <- list(NA, c(NA, NA))
  subin        <- subinr(B, unlist(lapply(obj, nlevels)))

 # INDICATOR APPROACH:
  nd.max     <- min(J-Q, I-1)
  if (is.na(nd) | nd > nd.max)
    nd <- nd.max
  P          <- Z / sum(Z)
  rm         <- apply(P, 1, sum)
  cm         <- apply(P, 2, sum)
  eP         <- rm %*% t(cm)
  S          <- (P-eP)/sqrt(eP)
  rowinertia <- apply(S^2, 1, sum)
  rowdist    <- sqrt(rowinertia / rm)
  rowmass    <- rep(1/I, I)
  dec        <- svd(S)
  lambda0    <- dec$d[1:nd.max]^2
 ## lambda0    <- dec$d[1:nd.max]
  rowcoord   <- as.matrix(dec$u[,1:nd.max]) / sqrt(rm)
  colcoord   <- as.matrix(dec$v[,1:nd.max]) / sqrt(cm)
  colinertia <- apply(S^2, 2, sum)
  coldist    <- sqrt(colinertia / cm)
  lambda.t   <- sum(lambda0)
  lambda.e   <- lambda0 / lambda.t
  lambda.et  <- 1
 # SIGNS FROM SVD AND EV ARE NOT THE SAME, ADJUSTMENT BELOW:
  coord.sign <- sign(colcoord[1,])
 # temporary fix
  colinertia.Burt <- NA
  lambda.Burt     <- NA
 # SUPPLEMENTARY COLUMNS FOR INDICATOR BELOW:
  if (!is.na(supcol[1])){
    cs.star       <- apply(Z.star, 2, sum)
    base.star     <- Z.star / matrix(rep(cs.star, I), nrow = I, byrow=T)
    svgam.star    <- matrix(sqrt(lambda0), nrow = sum(levels.n.sup), ncol = nd.max, byrow = TRUE)
    colcoord.star <- (t(base.star) %*% rowcoord) / svgam.star
    colcoord      <- rbind(colcoord, colcoord.star)
    colinertia    <- c(colinertia, rep(NA, J.star))
    coldist       <- c(coldist, rep(NA, J.star))
    cm            <- c(cm, rep(NA, J.star))
    cn            <- c(cn, dimnames(obj.supcol)[[2]])
    col.names.sup <- paste(fn.star, ln.star, sep = ".")
    if ((is.na(subsetcol)[1] & length(subsetcol) == 1)){
      col.names     <- c(col.names, col.names.sup)
      } else {
      col.names     <- c(col.names[subsetcol.in], col.names.sup)
      supcol <- supcol - (J - sum(!is.na(unlist(subsetcol))))
      }
   # Burt thingie for suppl
    B.star2       <- t(Z)%*%Z.star
    }

 # NON-INDICATOR CASES BELOW:
  if (lambda != "indicator"){
   # BURT, ADJUSTED AND JCA BELOW:
    nd.max     <- min(J-Q, I-1)

    cm         <- apply(P, 2, sum)
    P          <- B/sum(B)
    eP         <- cm %*% t(cm)
    S          <- (P - eP) / sqrt(eP)
### bcn, 2009_11:
   ##### PLACEHOLDER FOR SUBSET ! 
    if (!(is.na(subsetcol)[1] & length(subsetcol) == 1)){
      B         <- B[subsetcol.in,subsetcol.in]
      S         <- S[subsetcol.in,subsetcol.in]
      nd.max    <- dim(B)[1]
     # col.names <- col.names[subsetcol.in]
      cm        <- cm[subsetcol.in]
     # Z.star    <- Z.star[,subsetcol.in]
     # Z         <- Z[,subsetcol.in]
     # B.star2   <- t(Z)%*%Z.star
     # nd.max    <- min(length(subsetcol), I-1)
      }
   ##### EO PLACEHOLDER, SUBSET.
    if (is.na(nd) | nd > nd.max)
      nd <- nd.max

    dec        <- eigen(S)
   # lambda0    <- (dec$values[1:nd.max])
    lambda0    <- (dec$values[1:nd.max])^2
    colcoord   <- as.matrix(dec$vectors[,1:nd.max]) / sqrt(cm)
   # MATCH THE SIGNS FROM INDICATOR ANALYSIS (FOR ROWS)
    c.sign     <- sign(colcoord[1,])
    col.tr     <- rep(1, nd.max)
    col.tr[coord.sign[1:nd.max] != c.sign] <- -1
    rowcoord   <- rowcoord[,1:nd.max]%*%diag(col.tr)
    colinertia <- apply(S^2, 2, sum)
      colinertia.Burt <- colinertia  # temporary fix for COR/CTR in adjusted case
      lambda.Burt     <- lambda0     #  " "
    coldist    <- sqrt(colinertia / cm)
    lambda.t   <- sum(lambda0)
    lambda.e   <- lambda0 / lambda.t
    lambda.et  <- 1
   # ADJUSTED CASE BELOW:
    if (lambda == "adjusted"){
      nd.max      <- sum(sqrt(lambda0) >= 1/Q)
      if (is.na(nd) | nd > nd.max)
        nd <- nd.max
      lambda.adj  <- ((Q/(Q-1))^2 * (sqrt(lambda0)[1:nd.max] - 1/Q)^2)
      lambda.t    <- (Q/(Q-1)) * (sum(lambda0) - ((J - Q) / Q^2))
      lambda.e    <- lambda.adj / lambda.t
      lambda.et   <- NA
      lambda0     <- lambda.adj
      colinertia  <- (Q/(Q-1))^2 * (colinertia - (1/Q)*((1/Q)-cm))
      colcoord    <- as.matrix(dec$vectors[,1:nd.max]) / sqrt(cm)
      rowcoord    <- rowcoord[,1:nd.max]
      } else {
       # JCA CASE BELOW:
      if (lambda == "JCA"){
        nd.max     <- sum(sqrt(lambda0) >= 1/Q)
        if (is.na(nd) | nd > nd.max)
          nd <- nd.max
        B.it       <- iterate.mjca(B, lev.n = levels.n, nd = nd, maxit = maxit, epsilon = epsilon)
        B.star     <- B.it[[1]]
        JCA.it     <- B.it[[2]]
        subin      <- subinr(B.star, unlist(lapply(obj, nlevels)))
        P          <- B.star / sum(B.star)
        cm         <- apply(P, 2, sum)
        eP         <- cm %*% t(cm)
        S          <- (P - eP) / sqrt(eP)
        dec        <- eigen(S)
        lambda0    <- (dec$values[1:nd.max])^2
       ## lambda0    <- (dec$values[1:nd.max])
        colcoord   <- as.matrix(dec$vectors[,1:nd.max]) / sqrt(cm)
        rowcoord   <- rowcoord[,1:nd.max]
        colinertia <- apply(S^2, 2, sum)
        coldist    <- sqrt(colinertia / cm)
        lambda.e   <- rep(NA, nd.max)
        lambda.t   <- sum(subin)
        lambda.et  <- (sum(lambda0[1:nd]) - sum(diag(subin))) / (sum(subin)-sum(diag(subin)))
        } # end JCA
      }# end else
   # SUPPLEMENTARY CATEGORIES FOR BURT CASES:
    if (!is.na(supcol[1])){
      if (!(is.na(subsetcol)[1] & length(subsetcol) == 1)) {
        B.star2 <- B.star2[subsetcol.in,]
       # adjust supplementary indexes:
        } # else {
     #   B.star22 <- B.star2
     #   }
      rs.star      <- apply(B.star2, 2, sum)
      P.star       <- B.star2%*%diag(1/rs.star)
      cm           <- c(cm, apply(B.star2/sum(Z), 2, sum))
      colcoord.sup <- t(P.star)%*%colcoord
      colcoord.sup <- colcoord.sup%*%diag(1/sqrt(lambda0))
      colcoord     <- rbind(colcoord, colcoord.sup)
      colinertia   <- c(colinertia, apply(P.star^2, 2, sum))
      coldist      <- c(coldist, sqrt(apply(P.star^2, 2, sum)/apply(B.star2/sum(Z), 2, sum)))
      levels.n     <- c(levels.n, levels.n.sup)
      }
    }# end !indicator

  mjca.output <- list(sv         = sqrt(lambda0), 
                      lambda     = lambda,
                      inertia.e  = lambda.e,
                      inertia.t  = lambda.t,
                      inertia.et = lambda.et,
                      levelnames = col.names,
                      levels.n   = levels.n,
                      nd         = nd,
                      nd.max     = nd.max,
                      rownames   = rn, 
                      rowmass    = rowmass,
                      rowdist    = rowdist,
                      rowinertia = rowinertia,
                      rowcoord   = rowcoord,
                     # rowsup     = suprow, 
                      colnames   = cn, 
                      colmass    = cm, 
                      coldist    = coldist,
                      colinertia = colinertia, 
                      colcoord   = colcoord, 
                      colsup     = supcol,
                      subsetcol  = subsetcol.in,
                      Burt       = B,
                      Burt.upd   = B.star,
                      subinertia = subin,
                      JCA.iter   = JCA.it,
                     # bcn 09
                     # temp0      = colinertia.Burt, # temporary fix for COR/CTR in adjusted case
                     # temp1      = sqrt(lambda.Burt),     #  " "
                      call       = match.call())
  class(mjca.output) <- "mjca"
  return(mjca.output)

  }
