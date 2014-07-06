# MF: changed default ps=":" to allow working with the output coordinates
# MF: cleaned up code for checking lambda
# MF: mjca now accepts a table object, using expand.dft

mjca <- function(obj, 
                 nd = 2, 
                 lambda = c("adjusted", "indicator", "Burt", "JCA"), 
                 supcol = NA, 
                 subsetcol = NA, 
                 ps = ":", 
                 maxit = 50, 
                 epsilon = 0.0001){

################################################################################
 ##### Part 1: Input checks:
################################################################################
 ### Check for valid argument in 'lambda'
  lambda <- match.arg(lambda)
#  lam.v0 <- c("indicator", "Burt", "adjusted", "JCA")
#  lam.v1 <- c("i", "B", "a", "J")
#  if (length(grep(tolower(lambda), tolower(lam.v0), fixed = TRUE)) != 0){
#    lambda <- lam.v0[grep(tolower(lambda), tolower(lam.v0), fixed = TRUE)[1]]
#    } else {
# # Only first letter specified?
#    if (length(grep(tolower(lambda), lam.v1)) == 1){
#      lambda <- lam.v0[grep(tolower(lambda), lam.v1)]
#      } else {
# # "Fuzzy matching" unique?
#	  if (length(agrep(tolower(lambda), tolower(lam.v0))) == 1 ) {
#        lambda <- lam.v0[agrep(tolower(lambda), tolower(lam.v0))]
#        } else {
#        stop(paste("\nInvalid 'lambda' specification. Valid values are:\n", 
#                   paste("\"", lam.v0, "\" ", c(",", ",", "and ", "."), 
#                         sep = "", collapse = ""), collapse = "", sep = ""))
#        }
#      }
#    }
 ### End check 'lambda'

 ### check input data

 # allow for table input
 if (is.table(obj)) {
	 obj <- expand.dft(obj)
 }
 ### BELOW: edit from GR (2011-09):
 # if(!is.data.frame(obj)){
    obj <- data.frame(lapply(data.frame(obj), factor)) 
 #	}
 ### End check input data


################################################################################
 ##### Part 2: Data preparation
################################################################################
 # Indicator and Burt matrix:
  levels.n.0 <- unlist(lapply(obj, nlevels))
  rn.0       <- dimnames(obj)[[1]]
  cn.0       <- dimnames(obj)[[2]]
  n.0        <- cumsum(levels.n.0)
  I.0        <- dim(obj)[1]
  Q.0        <- dim(obj)[2]
  Q.sup      <- NA
  J.0        <- sum(levels.n.0)
  J.sup      <- NA
  Qind.0     <- 1:Q.0
  Qind       <- Qind.0
  Qind.sup   <- NA
  ind.0      <- 1:J.0
  ind.sub    <- NA # subset index, initially set to NA
  ind.sup.foo<- NA # subset index, initially set to NA
  ind.sup    <- NA # passive index, initially set to NA

  cn               <- dimnames(obj)[[2]]
  fn               <- rep(names(obj), unlist(lapply(obj, nlevels)))
  ln               <- unlist(lapply(obj,levels))
  col.names        <- paste(fn, ln, sep = ps)
 # dimnames(Z)[[2]] <- col.names
 # dimnames(Z)[[1]] <- as.character(1:I)
 # dimnames(B)[[2]] <- col.names
 # dimnames(B)[[1]] <- col.names


 # Subsets, supplementary check:
  if (!(is.na(subsetcol)[1] & length(subsetcol) == 1)){
   # check if given as vector
    if (mode(subsetcol) != "list"){
      if (sum(subsetcol < 0) == length(subsetcol)){ # check for negative indexes
        subsetcol <- (1:sum(levels.n))[subsetcol]
        }
      lut  <- cumsum(levels.n.0) - unlist(levels.n.0)
      s0   <- (1:sum(levels.n.0))[subsetcol]
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
      } # end subset-list
    subsetcol <- s0
    } # end subset 

 # Supplementary points:
  if (!is.na(supcol)[1]){
	Qind     <- Qind.0[-supcol]
	Qind.sup <- supcol
   # get indices for Burt matrix:
    for (k in 1:length(supcol)){
	  ind.sup <- c(ind.sup, (c(0,n.0)[supcol[k]] + 1):(c(0,n.0)[supcol[k]+1]))
	  }
    ind.sup  <- ind.sup[-1]
	ind      <- ind.0[-ind.sup]
	Q.sup    <- length(supcol)
	Q        <- Q.0 - Q.sup
	J.sup    <- sum(levels.n.0[Qind.sup])
	J        <- sum(levels.n.0[Qind])
   # check: subset and supplementary?
    ind.sup.foo <- ind.sup
	if (!(is.na(subsetcol)[1] & length(subsetcol) == 1)){
      ind.sup.foo <- ind.sup.foo - (length(ind) - length(subsetcol))
	  }
   #   for (i in sort(supcol, decr = TRUE)){
   #	    subsetcol <- subsetcol[-i]
   #     }
   #   } # end check subset&sup
    } else {
	  ind.sup <- NA
	  ind     <- ind.0
  	  Q.sup   <- NA
	  Q       <- Q.0
	  J.sup   <- NA
	  J       <- J.0
	} # end supplementary 

  levels.n <- levels.n.0[Qind]

 # Subset indexes:
  if (!(is.na(subsetcol)[1] & length(subsetcol) == 1)){
    ind.sub <- subsetcol
   # Levels in Subset:
#    levels.n.sub        <- table(rep(1:Q, each = levels.n)[subsetcol])
### CHECK THAT ONE BELOW: (GR; 2011)
    levels.n.sub        <- table(rep(1:Q, levels.n)[subsetcol])
    names(levels.n.sub) <- names(levels.n)
    Q.sub               <- Q - sum(levels.n.sub == 0)
    levels.n.sub        <- levels.n.sub[levels.n.sub != 0]
    }


################################################################################
 ##### Part 3: 'Core' Computation
################################################################################

 ### Set up data
 # Indicator and Burt matrix:
  Z.0      <- matrix(0, nrow = I.0, ncol = J.0)
  newdat   <- lapply(obj, as.numeric)
  offset.b <- c(0, n.0)
  offset   <- c(0, n.0[-length(n.0)])
  for (i in 1:Q.0) 
    Z.0[1:I.0 + (I.0 * (offset[i] + newdat[[i]] - 1))] <- 1
  fn         <- rep(names(obj), unlist(lapply(obj, nlevels)))
  ln         <- unlist(lapply(obj,levels))
  B.0        <- t(Z.0) %*% Z.0
  B          <- B.0[ind, ind]
  Z          <- Z.0[,ind]
  P          <- B / sum(B)
  cm         <- apply(P, 2, sum)
  rm         <- apply(Z.0/sum(Z.0), 1, sum)
  S          <- diag(sqrt(1/cm)) %*% (P - cm %*% t(cm)) %*% diag(sqrt(1/cm))
  evd.S      <- eigen(S)
### FIX (2011-09):
 #  obj.num    <- apply(obj[,Qind], 2, as.numeric)
  obj.num    <- as.matrix(data.frame(lapply(obj[,Qind], as.numeric)))
  rowmass    <- NA # rep(1/I.0, I.0)
  rowinertia <- NA # apply(S^2, 1, sum)
  rowdist    <- NA # sqrt(rowinertia / rm)
  colinertia <- apply(S^2, 2, sum)
  coldist    <- sqrt(colinertia / cm)

 # Burt bits for supplementary variables:
  if (!is.na(ind.sup[1])){
    B.sup <- B.0[ind.sup, ind]
    }

### FIX THIS BELOW (2011-09):

  if(!is.na(subsetcol[1])){

# added this here (2011-09):
    if(!is.na(ind.sup[1])){
      ind.sub <- c(ind.sup,ind.sub)
      ind.sub <- (ind.sub[!duplicated(ind.sub)])[-(1:length(ind.sup))]
      subsetcol <- ind.sub
      ind.sup.foo <- ind.sup - (length(ind.0)-length(c(ind.sub,ind.sup)))
     # ind.sup <- ind.sup.foo
      }
###
    B.sub <- B[ind.sub,ind.sub]
	}

 # some placeholders for adjusted/JCA
  B.star       <- NA
  lambda.adj   <- NA
  JCA.it       <- list(NA, c(NA, NA))
  subin        <- subinr(B, levels.n)

 # placxeholders for rows:
  row.ctr <- NA
  row.cor <- NA

 # coldist <- NA
 # colinertia <- NA

 ##### 3.1: 'lambda' = "indicator"
 
  if (lambda == "indicator"){
    nd.max  <- J - Q
    col.sc  <- diag(1 / sqrt(cm)) %*% evd.S$vectors[,1:(J-Q)] 
    col.pc  <- col.sc %*% diag(sqrt(evd.S$values[1:(J-Q)]))

   # Computations for rows:
    indices <- t(t(obj.num) + offset[Qind])
    row.pc  <- matrix(0, nrow = nrow(obj.num), ncol = J-Q)
    for(i in 1:nrow(obj.num)) 
      row.pc[i,] <- apply(col.sc[indices[i,],], 2, sum) / Q
   # (maybe you can figure out a way to do this without a loop!)
    row.sc  <- row.pc %*% diag(1/sqrt(evd.S$values[1:(J-Q)]))

    col.ctr <- evd.S$vectors[,1:(J-Q)]^2
    row.ctr <- (1/nrow(obj.num)) * row.pc^2 %*% diag(1/evd.S$values[1:(J-Q)])
    col.cor <- col.pc^2 / apply(col.pc^2, 1, sum)
    row.cor <- row.pc^2 / apply(row.pc^2, 1, sum)

    lambda0    <- evd.S$values[1:nd.max]
    lambda.t   <- sum(lambda0)
    lambda.e   <- lambda0 / lambda.t
    lambda.et  <- 1

   # EDIT (2011-09):
  #  subsetcol <- ind.sub 
  
   # Subset analysis:
    if (!is.na(subsetcol)[1]){
      nd.max  <- min(length(ind.sub), J-Q)
      evd.S   <- eigen(S[subsetcol,subsetcol])
      col.sc  <- diag(1 / sqrt(cm[subsetcol])) %*% evd.S$vectors[,1:nd.max] 
      col.pc  <- col.sc %*% diag(sqrt(evd.S$values[1:nd.max]))
      col.ctr <- evd.S$vectors[,1:nd.max]^2
      col.cor <- col.pc^2 / apply(col.pc^2, 1, sum)

	  lookup  <- offset[1:Q]
      rpm     <- as.matrix(data.frame(newdat))[,1:Q]
      indices <- t(t(rpm) + lookup)
	  row.pc  <- matrix(0, nrow = I.0, ncol = nd.max)
      for(i in 1:(I.0)) {
        profile              <- -cm
        profile[indices[i,]] <- profile[indices[i,]]+1/Q
        profile              <- profile[subsetcol]
        row.pc[i,]           <- t(profile) %*% col.sc
        }
      row.sc  <- row.pc %*% diag(1/sqrt(evd.S$values[1:nd.max]))
	  row.ctr <- (1/I.0) * row.pc^2 %*% diag(1/evd.S$values[1:nd.max])
      row.cor <- row.pc^2 / apply(row.pc^2, 1, sum)
     # Subset & Supplementary variables:
	  if(!is.na(supcol)[1]){
	    cols.pc  <- sweep((B.sup / apply(B.sup, 1, sum))[,subsetcol], 2, 
	                      cm[subsetcol]) %*% col.sc
        cols.sc  <- cols.pc %*% diag(1/sqrt(evd.S$values[1:nd.max]))
        cols.cor <- cols.pc^2 / apply(cols.pc^2,1,sum)
		}
	  } # END Subset

   # Supplementary points:
    if (!is.na(supcol)[1] & is.na(subsetcol)[1]){  
     # B.sup    <- B.0[ind.sup, ind]
      cols.pc  <- (B.sup / apply(B.sup, 1, sum)) %*% col.sc
      cols.sc  <- cols.pc %*% diag(1 / evd.S$values[1:nd.max])
      cols.sqd <- apply((sweep(sweep((B.sup / apply(B.sup,1,sum)), 2, cm), 2, 
                               sqrt(cm), FUN = "/"))^2, 1, sum)
      cols.cor <- cols.pc^2 / apply(cols.pc^2, 1, sum)
      } # End supplementary

	  
    } else # END if "indicator"
	{

 ##### 3.2: 'lambda' = "Burt"
    col.sc  <- diag(1/sqrt(cm)) %*% evd.S$vectors[,1:(J-Q)] 
    col.pc  <- col.sc %*% diag(evd.S$values[1:(J-Q)])
    row.sc  <- col.sc
    row.pc  <- col.pc

    col.ctr <- evd.S$vectors[,1:(J-Q)]^2
    col.cor <- col.pc^2 / apply(col.pc^2, 1, sum)

   # Subset analysis:
    if (!is.na(subsetcol)[1]){
      nd.max  <- min(length(ind.sub), J-Q)
      evd.S   <- eigen(S[subsetcol,subsetcol])
      col.sc  <- diag(1/sqrt(cm[subsetcol])) %*% evd.S$vectors[,1:nd.max] 
      col.pc  <- col.sc %*% diag(evd.S$values[1:nd.max])
      row.sc  <- col.sc
      row.pc  <- col.pc
      col.ctr <- evd.S$vectors[,1:nd.max]^2
      col.cor <- col.pc^2 / apply(col.pc^2, 1, sum)
     # Subset & Supplementary variables:
	  if(!is.na(supcol)[1]){
	    cols.pc  <- sweep((B.sup / apply(B.sup, 1, sum))[,subsetcol], 2, 
	                      cm[subsetcol]) %*% col.sc
        cols.sc  <- cols.pc %*% diag(1 / evd.S$values[1:nd.max])
        cols.cor <- cols.pc^2 / apply(cols.pc^2, 1, sum)
        }
	  } # End Subset

   # Supplementary points:
    if (!is.na(supcol)[1] & is.na(subsetcol)[1]){
      B.sup    <- B.0[ind.sup, 1:J]
      cols.pc  <- (B.sup / apply(B.sup, 1, sum)) %*% col.sc
      cols.sc  <- cols.pc %*% diag(1/evd.S$values[1:(J-Q)])
      cols.cor <- cols.pc^2 / apply(cols.pc^2, 1, sum)
      } # End supplementary

   # if(!is.na(supcol)[1]){
	#  colinertia <- c(colinertia, rep(NA, J.sup))
	#  coldist    <- c(coldist, rep(NA, J.sup))
	#  }

    nd.max     <- J - Q
    lambda0    <- (evd.S$values[1:nd.max])^2
    lambda.t   <- sum(lambda0)
    lambda.e   <- lambda0 / lambda.t
    lambda.et  <- 1

 ##### 3.3: 'lambda' = "adjusted"
    if (lambda != "Burt"){
      nd.max <- sum(sqrt(lambda0) >= 1/Q)
      B.null <- B - diag(diag(B))
      P.null <- B.null / sum(B.null)
      S.null <- diag(sqrt(1/cm)) %*% (P.null - cm %*% t(cm)) %*% 
                  diag(sqrt(1/cm))
  
      evd.S.null <- eigen(S.null)
      K0         <- length(which(evd.S.null$values > 1e-8))
 
      Pe <- P
      for(q in 1:Q) 
        Pe[(offset.b[q]+1):offset.b[q+1],(offset.b[q]+1):offset.b[q+1]] <- 
          cm[(offset.b[q]+1):offset.b[q+1]] %*% 
          t(cm[(offset.b[q]+1):offset.b[q+1]])
      Se          <- diag(sqrt(1/cm)) %*% (Pe-cm%*%t(cm)) %*% diag(sqrt(1/cm))
      inertia.adj <- sum(Se^2) * Q / (Q-1)

      col.sc  <- diag(1/sqrt(cm)) %*% evd.S.null$vectors[,1:K0] 
      col.pc  <- col.sc %*% diag(evd.S.null$values[1:K0])
      row.sc  <- col.sc
      row.pc  <- col.pc

      col.ctr <- evd.S.null$vectors[,1:K0]^2

      col.inr.adj <- apply(Se^2,2,sum) * Q/(Q-1)
      col.cor     <- diag(cm) %*% col.pc^2 / col.inr.adj

      lambda.adj  <- ((Q/(Q-1))^2 * (sqrt(lambda0)[1:nd.max] - 1/Q)^2)
      lambda.t    <- (Q/(Q-1)) * (sum(lambda0) - ((J - Q) / Q^2))
      lambda.e    <- lambda.adj / lambda.t
      lambda.et   <- NA
      lambda0     <- lambda.adj

     # Subset analysis:
      if (!is.na(subsetcol)[1]){
		evd.S0 <- eigen(S.null[subsetcol,subsetcol])
        K0     <- length(which(evd.S0$values>1e-8))
		nd.max <- K0
		lookup <- offset.b[1:(Q+1)]
		Pe <- P
        for(q in 1:Q) 
		  Pe[(lookup[q]+1):lookup[q+1],(lookup[q]+1):lookup[q+1]] <- 
		    cm[(lookup[q]+1):lookup[q+1]]%*%t(cm[(lookup[q]+1):lookup[q+1]])
        Se <- diag(sqrt(1/cm))%*%(Pe-cm%*%t(cm))%*%diag(sqrt(1/cm))
       # inertia.adj.subset <- sum(Se[subsetcol,subsetcol]^2)*Q / (Q-1)
        lambda.adj  <- evd.S0$values[1:K0]^2
        lambda0     <- lambda.adj
		lambda.t    <- sum(Se[subsetcol,subsetcol]^2)*Q / (Q-1)
        lambda.e    <- lambda.adj / lambda.t
		col.sc  <- diag(1/sqrt(cm[subsetcol])) %*% evd.S0$vectors[,1:K0] 
        col.pc  <- col.sc %*% diag(evd.S0$values[1:K0])
        row.sc  <- col.sc
        row.pc  <- col.pc
        col.ctr <- evd.S0$vectors[,1:K0]^2
        col.inr.adj.subset <- apply(Se[subsetcol,subsetcol]^2,2,sum)*Q/(Q-1)
        col.cor <- diag(cm[subsetcol]) %*% col.pc^2 / col.inr.adj.subset
       # Subset & Supplementary variables:
        if(!is.na(supcol)[1]){
		  cols.pc  <- sweep((B.sup / apply(B.sup, 1, sum))[,subsetcol], 2, 
		                    cm[subsetcol]) %*% col.sc
          cols.sc  <- cols.pc %*% diag(1/evd.S0$values[1:K0])
          cols.sqd <- apply((sweep(sweep((B.sup/apply(B.sup, 1, 
                         sum))[,subsetcol], 2, cm[subsetcol]), 2, 
                         sqrt(cm[subsetcol]), FUN="/"))^2, 1, sum)
          cols.cor <- cols.pc^2 / cols.sqd
          }
        } # End Subset

     # Supplementary points:
      if (!is.na(supcol)[1] & is.na(subsetcol)[1]){
        B.sup    <- B.0[ind.sup, 1:J]
        cols.pc  <- (B.sup / apply(B.sup, 1, sum)) %*% col.sc
        cols.sc  <- cols.pc %*% diag(1/evd.S.null$values[1:K0])
        cols.cor <- cols.pc^2 / apply(cols.pc^2, 1, sum)
        }

 ##### 3.4: 'lambda' = "JCA"
      if (lambda == "JCA"){
        if (is.na(nd) | nd > nd.max)
          nd <- nd.max
        B.it        <- iterate.mjca(B, lev.n = levels.n, nd = nd, 
                                    maxit = maxit, epsilon = epsilon)
        B.star      <- B.it[[1]]
        JCA.it      <- B.it[[2]]
        subin       <- subinr(B.star, levels.n)
        P           <- B.star / sum(B.star)
        cm          <- apply(P, 2, sum)
        eP          <- cm %*% t(cm)
        S           <- (P - eP) / sqrt(eP)
        dec         <- eigen(S)
        lambda0     <- (dec$values[1:nd.max])^2

        col.sc      <- as.matrix(dec$vectors[,1:nd.max]) / sqrt(cm)
        col.pc      <- col.sc %*% diag(sqrt(lambda0))
        row.sc      <- col.sc
        row.pc      <- col.pc

        inertia.mod      <- sum(subin - diag(diag(subin)))
        inertia.discount <- sum(diag(subin))
        inertia.expl     <- (sum(lambda0[1:nd]) - inertia.discount) / 
                             inertia.adj

        lambda.e   <- rep(NA, nd.max)
        lambda.t   <- sum(subin)
        lambda.et  <- (sum(lambda0[1:nd]) - sum(diag(subin))) / 
                      (sum(subin)-sum(diag(subin)))

        Pm <- B.star / sum(B.star)
        Sm <- diag(sqrt(1/cm))%*%(Pm-cm%*%t(cm))%*%diag(sqrt(1/cm))

        inertia.col.discount <- rep(0, J)
        for(q in 1:Q) 
          inertia.col.discount[(offset.b[q]+1):(offset.b[q+1])] <- 
            apply(Sm[(offset.b[q]+1):(offset.b[q+1]),
                     (offset.b[q]+1):(offset.b[q+1])]^2, 2, sum)
        inertia.col.adj <- apply(Sm^2, 2, sum) - inertia.col.discount

        col.ctr <- (apply(cm*col.pc[,1:nd]^2, 1, sum) - inertia.col.discount)/
                     (sum(lambda0[1:nd])-inertia.discount)
        col.cor <- (apply(cm*col.pc[,1:nd]^2, 1, sum) - inertia.col.discount) / 
                   inertia.col.adj

       # Subset analysis:
        if (!is.na(subsetcol)[1]){
         # template matrix:
          foo0 <- rep(1:Q.sub, each = levels.n.sub)
          foo1 <- (foo0) %*% t(rep(1, sum(levels.n.sub))) - t((foo0) %*% 
                  t(rep(1, sum(levels.n.sub))))
          upd.template <- ifelse(foo1 == 0, TRUE, FALSE)
          cat.template <- rep(FALSE, J)
		  cat.template[subsetcol] <- TRUE
		  
		  Bsub.margin     <- apply(B.star, 1, sum) / sum(B.star)
		  Bsub.red.margin <- Bsub.margin[cat.template]
		  Bsub.red        <- B.star[cat.template,cat.template]
		  Bsub.red.P      <- Bsub.red / sum(B.star)
		  Bsub.red.S      <- diag(1/sqrt(Bsub.red.margin)) %*% (Bsub.red.P - 
		                       Bsub.red.margin %*% t(Bsub.red.margin)) %*% 
		                       diag(1/sqrt(Bsub.red.margin))
		  Bsub.red.SVD    <- svd(Bsub.red.S)		  
		  Bsub.red.est    <- Bsub.red.margin %*% t(Bsub.red.margin) + 
		                       diag(sqrt(Bsub.red.margin)) %*% 
		                       (Bsub.red.SVD$u[,1:nd] %*% 
		                       diag(Bsub.red.SVD$d[1:nd]) %*% 
		                       t(Bsub.red.SVD$v[,1:nd])) %*% 
		                       diag(sqrt(Bsub.red.margin))
          Bsub.red.P.mod  <- (1-upd.template) * Bsub.red.P + upd.template * 
                               Bsub.red.est
          Bsub.red.S      <- diag(1/sqrt(Bsub.red.margin)) %*% 
                               (Bsub.red.P.mod - Bsub.red.margin %*% 
                               t(Bsub.red.margin)) %*% 
                               diag(1/sqrt(Bsub.red.margin))
          Bsub.red.SVD    <- svd(Bsub.red.S)
		  
		 # iterations:
		  it <- TRUE
          k  <- 0 
		  while (it){
		    Bsub.red.P.previous <- Bsub.red.P.mod
			Bsub.red.est        <- Bsub.red.margin %*% t(Bsub.red.margin) + 
			                         diag(sqrt(Bsub.red.margin)) %*% 
			                         (Bsub.red.SVD$u[,1:nd] %*% 
			                         diag(Bsub.red.SVD$d[1:nd]) %*% 
			                         t(Bsub.red.SVD$v[,1:nd])) %*% 
			                         diag(sqrt(Bsub.red.margin))
			Bsub.red.P.mod      <- (1 - upd.template) * Bsub.red.P + 
			                         upd.template * Bsub.red.est
			Bsub.red.S          <- diag(1/sqrt(Bsub.red.margin)) %*% 
			                         (Bsub.red.P.mod - Bsub.red.margin %*% 
			                         t(Bsub.red.margin)) %*% 
			                         diag(1/sqrt(Bsub.red.margin))
			Bsub.red.SVD        <- svd(Bsub.red.S)
			if (max(abs(Bsub.red.P.mod - Bsub.red.P.previous)) < epsilon | 
			    k >= maxit){
			  it <- FALSE
			  }
            k <- k + 1
		    }
		  
		  inertia.adj.red.discount <- sum((upd.template * Bsub.red.S)^2)
          inertia.adj.red.subset   <- sum(Bsub.red.S^2) - 
                                      inertia.adj.red.discount
		  col.sc  <- sqrt(1/cm[cat.template]) * Bsub.red.SVD$v[,1:nd]
		  col.pc  <- col.sc %*% diag(Bsub.red.SVD$d[1:nd])
		  row.sc  <- col.sc
		  row.pc  <- col.pc

		  Sm <- Bsub.red.S
		  inertia.col.red.discount <- apply((upd.template * Sm)^2, 2, sum )
		  inertia.col.red.adj      <- apply(Sm^2, 2, sum) - 
		                              inertia.col.red.discount #[subsetcol]
		  col.ctr <- (apply(cm[cat.template]*col.pc[,1:nd]^2, 1, sum) - 
		              inertia.col.red.discount)/(sum(Bsub.red.SVD$d[1:nd]^2) - 
		              inertia.adj.red.discount)
		  col.cor <- (apply(cm[cat.template]*col.pc[,1:nd]^2, 1, sum) - 
		              inertia.col.red.discount) / inertia.col.red.adj

         # Subset & Supplementary variables:
          if(!is.na(supcol)[1]){
		    cols.pc  <- sweep((B.sup / apply(B.sup, 1, sum))[,cat.template], 2, 
		                      cm[cat.template]) %*% col.sc
			cols.sc  <- cols.pc %*% diag(1 / Bsub.red.SVD$d[1:nd])
			cols.sqd <- apply((sweep(sweep((B.sup / 
			              apply(B.sup, 1, sum))[,cat.template], 2, 
			              cm[cat.template]), 2, sqrt(cm[cat.template]), 
			              FUN = "/"))^2, 1, sum)
			cols.cor <- apply(cols.pc[,1:nd]^2, 1, sum) / cols.sqd
		    }
	 	  } # End Subset

       # Supplementary points:
        if (!is.na(supcol)[1] & is.na(subsetcol)[1]){
          B.sup    <- B.0[ind.sup, 1:J]
          cols.pc  <- (B.sup / apply(B.sup, 1, sum)) %*% col.sc
          cols.sc  <- cols.pc %*% diag(1 / lambda0)
          cols.sqd <- apply((sweep(sweep((B.sup/apply(B.sup,1,sum)), 2, cm), 2, 
                                   sqrt(cm), FUN="/"))^2, 1, sum)
          cols.cor <- apply(cols.pc[,1:nd]^2, 1, sum) / cols.sqd
          }


        } # END if "JCA"
      } # END if !"Burt"
    } # END else if "indicator"

  if (!is.na(supcol)[1]){
    colcoord  <- rbind(col.sc, cols.sc)
    colpcoord <- rbind(col.pc, cols.pc)
    if (lambda != "JCA"){
      col.ctr   <- rbind(col.ctr, matrix(NA, nrow = length(ind.sup), 
                         ncol = ncol(col.ctr)))
	  colcor    <- rbind(col.cor, cols.cor)
	  } else {
      col.ctr <- c(col.ctr, rep(NA, length(ind.sup)))
	  colcor  <- c(col.cor, cols.cor)
	  }
    } else {
    colcoord  <- col.sc
    colpcoord <- col.pc
    colcor    <- col.cor
    }
  col.names0 <- col.names
  if (!is.na(subsetcol[1])){
    cm         <- cm[ind.sub]
	coldist    <- coldist[ind.sub]
	colinertia <- colinertia[ind.sub]
    col.names0 <- col.names[ind.sub]
    if(!is.na(supcol)[1]){
      col.names0 <- c(col.names0,col.names[ind.sup])
	  }
	B.out <- B.sub
    } else {
	B.out <- B.0
	}
  colctr    <- col.ctr
  rowcoord  <- row.sc
  rowpcoord <- row.pc

  if(!is.na(supcol)[1]){
    colinertia <- c(colinertia, rep(NA, J.sup))
    coldist    <- c(coldist, rep(NA, J.sup))
    cm         <- c(cm, rep(NA, J.sup))
    }
  col.names <- col.names0

# balkan solution for colcor > 1 (2011-09):
# adjusted analysis only!
if (lambda == "adjusted"){
  foo0 <- apply(colcor,1,max) > 1
  if (sum(foo0) > 0){
    insert <- c(1, rep(0, dim(colcor)[2]-1))
    colcor[foo0,] <- matrix(insert, nrow = sum(foo0), ncol = dim(colcor)[2], 
                            byrow = TRUE)
    }
  }
 # wrap up results
  mjca.output <- list(sv         = sqrt(lambda0), 
                      lambda     = lambda,
                      inertia.e  = lambda.e,
                      inertia.t  = lambda.t,
                      inertia.et = lambda.et,
                      levelnames = col.names,
                      levels.n   = levels.n.0,
                      nd         = nd,
                      nd.max     = nd.max,

                      rownames   = rn.0, 
                      rowmass    = rowmass,
                      rowdist    = rowdist,
                      rowinertia = rowinertia,
                      rowcoord   = rowcoord,
					  
					  rowpcoord  = rowpcoord,
					  rowctr     = row.ctr,
					  rowcor     = row.cor,

                      colnames   = cn.0, 
                      colmass    = cm, 
                      coldist    = coldist,
                      colinertia = colinertia, 
                      colcoord   = colcoord, 
					  
					  colpcoord  = colpcoord,
					  colctr     = col.ctr,
					  colcor     = colcor,

                      colsup     = ind.sup.foo,
                      subsetcol  = subsetcol,
                      Burt       = B.out,
                      Burt.upd   = B.star,
                      subinertia = subin,
                      JCA.iter   = JCA.it,
                      call       = match.call())
  class(mjca.output) <- "mjca"
  return(mjca.output)
  }


