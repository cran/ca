
################################################################################
# 
#     ca:  Computation of Simple CA
#
#  Input:  obj       Two-way freq. table
#          nd        Number of dimensions to be included (default: NA)
#          suprow    Indexes of supplementary columns (default: NA)
#          supcol    Indexes of supplementary rows (default: NA)
#          subsetrow Column indexes of subset
#          subsetcol Row indexes of subset
#
# Output:  An object of type 'ca'
# 
################################################################################

ca <- function(obj, 
               nd        = NA, 
               suprow    = NA, 
               supcol    = NA,
               subsetrow = NA,
               subsetcol = NA)
{
  nd0 <- nd
  I  <- dim(obj)[1] ; J <- dim(obj)[2]
  rn <- dimnames(obj)[[1]]
  cn <- dimnames(obj)[[2]]
  N  <- matrix(as.matrix(obj), nrow = I, ncol = J)

 # Temporary remove supplementray rows/columns:
  Ntemp <- N ; NtempC <- NtempR <- N
  suprow <- sort(suprow) ; supcol <- sort(supcol)

 # back to supplementary profiles:
  if (!is.na(supcol[1]) & !is.na(suprow[1])) { 
    NtempC <- Ntemp[-suprow,] 
    NtempR <- Ntemp[,-supcol] 
    }
  if (!is.na(supcol[1])) { 
    SC <- as.matrix(NtempC[,supcol])
    Ntemp <- Ntemp[,-supcol]
   }
  if (!is.na(suprow[1])) {
    SR <- matrix(as.matrix(NtempR[suprow,]), nrow = length(suprow))
    Ntemp <- Ntemp[-suprow,]
    }
  N <- matrix(as.matrix(Ntemp), nrow = dim(Ntemp)[1], ncol = dim(Ntemp)[2])

 # Adjustment for subset CA
  subsetrowt <- subsetrow
  if (!is.na(subsetrow[1]) & !is.na(suprow[1])) {
    subsetrowi <- subsetrow
    subsetrowt <- sort(c(subsetrow, suprow))
    subsetrowt <- subsetrowt[!duplicated(subsetrowt)]
    I <- length(subsetrowt)
   # adjust subset index
    for (q in length(suprow):1) {
      subsetrow <- subsetrow[subsetrow != suprow[q]]
      subsetrow <- subsetrow - as.numeric(suprow[q] < subsetrow)
      }
   # adjust supplementary row indexes:
    for (q in 1:length(suprow)) suprow[q] <- (1:length(subsetrowt))[subsetrowt == suprow[q]]
    }
  subsetcolt <- subsetcol
  if (!is.na(subsetcol[1]) & !is.na(supcol[1])) {
    subsetcoli <- subsetcol
    subsetcolt <- sort(c(subsetcol, supcol))
    subsetcolt <- subsetcolt[!duplicated(subsetcolt)]
    J <- length(subsetcolt)
   # adjust subset index
    for (q in length(supcol):1) {
      subsetcol <- subsetcol[subsetcol != supcol[q]]
      subsetcol <- subsetcol - as.numeric(supcol[q] < subsetcol)
      }
   # adjust supplementary column indexes:
    for (q in 1:length(supcol)) supcol[q] <- (1:length(subsetcolt))[subsetcolt == supcol[q]]
    }


 # check for subset CA
  dim.N <- dim(N)
  if (!is.na(subsetrow[1])) {
    dim.N[1] <- min(c(length(subsetrow) + 1, dim.N[1]))
    if (!is.na(supcol[1])) SC <- as.matrix(SC[subsetrow,])
    }
  if (!is.na(subsetcol[1])) {
    dim.N[2] <- min(c(length(subsetcol) + 1, dim.N[2]))
    if (!is.na(suprow[1])) SR <- matrix(as.matrix(SR[,subsetcol]), nrow = length(suprow))
    }
 # end subset CA
  nd.max <- min(dim.N) - 1
  if (is.na(nd) | nd > nd.max ) nd <- nd.max

 # Init:
  n <- sum(N) ; P <- N/n
  rm <- apply(P, 1, sum) ; cm <- apply(P, 2, sum)

 # SVD:
  eP     <- rm %*% t(cm)
  eN     <- eP * n
  S      <- (P - eP) / sqrt(eP)
 # subset CA
  if (!is.na(subsetcol[1])) { 
    S  <- S[,subsetcol] 
    cm <- cm[subsetcol]
    cn <- cn[subsetcolt] }
  if (!is.na(subsetrow[1])) { 
    S  <- S[subsetrow,]
    rm <- rm[subsetrow]
    rn <- rn[subsetrowt] }
 # end sCA
  chimat <- S^2 * n
  dec    <- svd(S)
  #sv     <- dec$d[1:nd]
  sv     <- dec$d[1:(min(dim(S))-1)]
  u      <- dec$u
  v      <- dec$v
  ev     <- sv^2
  cumev  <- cumsum(ev)

 # Inertia:
  totin <- sum(ev)
  rin   <- apply(S^2, 1, sum)
  cin   <- apply(S^2, 2, sum)

 # chidist
  rachidist <- sqrt(rin / rm)
  cachidist <- sqrt(cin / cm)
  rchidist  <- rep(NA, I) 
  cchidist  <- rep(NA, J) 

  if (!is.na(subsetrow[1])) {
    obj <- obj[subsetrowt,]
    }
  if (!is.na(subsetcol[1])) {
    obj <- obj[,subsetcolt]
    }

 # supplementary rows/columns:
  if (!is.na(suprow[1])) {
    if (is.na(supcol[1])) {
      P.stemp <- matrix(as.matrix(obj[suprow,]), nrow = length(suprow))
      } else P.stemp <- matrix(as.matrix(obj[suprow, -supcol]), nrow = length(suprow))
    P.stemp           <- P.stemp / apply(P.stemp, 1, sum)
    P.stemp           <- t((t(P.stemp) - cm) / sqrt(cm))
    rschidist         <- sqrt(apply(P.stemp^2, 1, sum))
    rchidist[-suprow] <- rachidist
    rchidist[suprow]  <- rschidist
    } else rchidist <- rachidist

  if (!is.na(supcol[1])) {
    if (is.na(suprow[1])) {
      P.stemp <- as.matrix(obj[, supcol])
      } else P.stemp <- as.matrix(obj[-suprow, supcol])
    P.stemp           <- t(t(P.stemp) / apply(P.stemp, 2, sum))
    P.stemp           <- (P.stemp - rm) / sqrt(rm)
    cschidist         <- sqrt(apply(P.stemp^2, 2, sum))
    cchidist[-supcol] <- cachidist
    cchidist[supcol]  <- cschidist
    } else cchidist <- cachidist

 # Standard coordinates:
  phi <- as.matrix(u[,1:nd]) / sqrt(rm)
  gam <- as.matrix(v[,1:nd]) / sqrt(cm)

 # Standard coordinates for supplementary rows/columns
  if (!is.na(suprow[1])) {
    rs             <- apply(SR, 1, sum)
    base2          <- SR / matrix(rep(rs, dim(gam)[1]), ncol = dim(gam)[1])
    svphi          <- matrix(sv[1:nd], nrow = length(suprow), ncol = nd, 
                             byrow = TRUE)
    phi2           <- (as.matrix(base2) %*% gam) / svphi
    phi3           <- matrix(NA, ncol = nd, nrow = I)
    phi3[suprow,]  <- phi2
    phi3[-suprow,] <- phi
    rm0            <- rep(NA, I)
    rm0[-suprow]   <- rm
    rm             <- rm0
    P.star        <- SR / n
    rm[suprow]    <- apply(P.star, 1, sum)
    rin0           <- rep(NA, I)
    rin0[-suprow]  <- rin
    rin            <- rin0
    }
  if (!is.na(supcol[1])) {
    cs             <- apply(SC, 2, sum)
    base2          <- SC / matrix(rep(cs, dim(phi)[1]), nrow = dim(phi)[1], 
                                  byrow = TRUE)
    svgam          <- matrix(sv[1:nd], nrow = length(supcol), ncol = nd, 
                             byrow = TRUE)
    gam2           <- (as.matrix(t(base2)) %*% phi) / svgam
    gam3           <- matrix(NA, ncol = nd, nrow = J)
    gam3[supcol,]  <- gam2
    gam3[-supcol,] <- gam
    cm0            <- rep(NA, J)
    cm0[-supcol]   <- cm
    cm             <- cm0
    P.star        <- SC / n
    cm[supcol]    <- apply(P.star, 2, sum)
    cin0           <- rep(NA, J)
    cin0[-supcol]  <- cin
    cin            <- cin0
    }
  if (exists("phi3")) phi <- phi3
  if (exists("gam3")) gam <- gam3

ca.output <- 
list(sv         = sv, 
     nd         = nd0,
     rownames   = rn, 
     rowmass    = rm, 
     rowdist    = rchidist,
     rowinertia = rin, 
     rowcoord   = phi, 
     rowsup     = suprow, 
     colnames   = cn, 
     colmass    = cm, 
     coldist    = cchidist,
     colinertia = cin, 
     colcoord   = gam, 
     colsup     = supcol,
     call       = match.call())
class(ca.output) <- "ca"
return(ca.output)  
}

################################################################################
