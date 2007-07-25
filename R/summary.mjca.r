################################################################################
# 
# summary.mjca: Summarizing 'mjca'-objects
#
#      Input: object The 'mjca' object which should be summarized
#             scree  Logical indicating if a scree-plot should be included
#                    (default: TRUE)
#             nd     Maximum dimension to include in scree-plot (default: 0=all)
#
#     Output: summary for 'mjca' object
# 
################################################################################


summary.mjca <- function(object, scree = TRUE, rows = FALSE, ...)
  {
  obj <- object
  nd  <- obj$nd
  if (is.na(nd)){
    nd <- 2
    } else {
    if (nd > length(obj$sv)) nd <- length(obj$sv)
    }

 # principal coordinates:
 # K   <- dim(obj$rowcoord)[2]
  K   <- obj$nd.max
  I   <- length(obj$rownames) # dim(obj$rowcoord)[1] 
  J   <- dim(obj$colcoord)[1]
  evF <- matrix(rep(sqrt(obj$sv[1:K]), I), I, K, byrow = TRUE)
  evG <- matrix(rep(sqrt(obj$sv[1:K]), J), J, K, byrow = TRUE)
  rpc <- obj$rowcoord[,1:K] * evF
  cpc <- obj$colcoord[,1:K] * evG
 # row profiles:
  r.names <- abbreviate(obj$rownames, 3)
  sr      <- obj$rowsup
 # if (!is.na(sr[1])) r.names[sr] <- paste("(*)", r.names[sr], sep = "")
  r.mass  <- obj$rowmass
  r.inr   <- obj$rowinertia / sum(obj$rowinertia, na.rm = TRUE)
  r.ccc   <- matrix(NA, nrow = length(r.names), ncol = nd * 3)
  for (i in 1:nd)  {
    r.ccc[,3 * (i - 1) + 1] <- rpc[,i]
    r.ccc[,3 * (i - 1) + 2] <- obj$rowmass * rpc[,i]^2 / 
                               obj$rowinertia
    r.ccc[,3 * (i - 1) + 3] <- obj$rowmass * rpc[,i]^2 /
                               obj$sv[i]
    if (obj$lambda == "indicator"){
      r.ccc[,3 * (i - 1) + 3] <- obj$rowmass * rpc[,i]^2 /
                                 sqrt(obj$sv[i])
      }
    }
  if (nd > 1) {
    r.qlt <- apply(r.ccc[,((1:nd-1) * 3 + 2)], 1, sum) 
    } else {
    r.qlt <- r.ccc[,((1:nd-1) * 3 + 2)] 
    }

  r1              <- paste(" k=", 1:nd, sep = "")
  r2              <- rep("cor", nd)
  r3              <- rep("ctr", nd)
  rcclab          <- as.vector(rbind(r1, r2, r3))
  dimnames(r.ccc) <- list(r.names, rcclab)
  r.out           <- data.frame(r.names, 
                                round(1000 * r.mass, 0), 
                                round(1000 * r.qlt, 0),
                                round(1000 * r.inr, 0), 
                                round(1000 * r.ccc, 0))
  dimnames(r.out) <- list(as.character(1:length(r.names)),
                          c("name", "mass", " qlt", " inr", rcclab))

 # column profiles:
  getfirst <- function(input) input[1]
  getlast  <- function(input) input[length(input)]
  c.part1  <- unlist(lapply(strsplit(obj$levelnames, "\\."), getfirst))
  c.part2  <- unlist(lapply(strsplit(obj$levelnames, "\\."), getlast))
  c.names  <- paste(abbreviate(c.part1, 3), c.part2, sep = ".") 
  sc       <- obj$colsup
  if (!is.na(sc[1])) c.names[sc] <- paste("(*)", c.names[sc], sep = "")
  c.mass   <- obj$colmass
  c.inr    <- obj$colinertia / sum(obj$colinertia, na.rm = TRUE)
  c.ccc    <- matrix(NA, nrow = length(c.names), ncol = nd * 3)
  for (i in 1:nd){
    c.ccc[,3 * (i - 1) + 1] <- cpc[,i]
    c.ccc[,3 * (i - 1) + 2] <- obj$colmass * cpc[,i]^2 /
                               obj$colinertia
    c.ccc[,3 * (i - 1) + 3] <- obj$colmass * cpc[,i]^2 / 
                               obj$sv[i]
    if (obj$lambda == "indicator"){
      c.ccc[,3 * (i - 1) + 3] <- obj$colmass * cpc[,i]^2 / 
                                 sqrt(obj$sv[i])
      }
    }
  if (obj$lambda == "adjusted"){
    temp.pc  <- obj$colcoord[,1:length(obj$sv)] %*% diag(sqrt(obj$sv))
    temp.rcc <- diag(obj$colmass)%*%temp.pc^2
    cpc.new  <- diag(1/obj$colinertia)%*%temp.rcc
    c.ccc[,3 * ((1:nd) - 1) + 2] <- cpc.new[,1:nd]
    }
 # cor and quality for supplementary columns
  if (!is.na(obj$colsup[1])){
    i0 <- obj$colsup
    for (i in 1:nd){
      c.ccc[i0,3 * (i - 1) + 2] <- obj$colmass[i0] * cpc[i0,i]^2 / obj$colinertia[i0]
      c.ccc[i0,3 * (i - 1) + 3] <- NA
      }
    }

  if (nd > 1) { 
    c.qlt <- apply(c.ccc[,((1:nd - 1) * 3 + 2)], 1, sum) 
    } else {
    c.qlt <- c.ccc[,((1:nd - 1) * 3 + 2)] 
    }

  c1              <- paste(" k=", 1:nd, sep = "")
  c2              <- rep("cor", nd)
  c3              <- rep("ctr", nd)
  ccclab          <- as.vector(rbind(c1, c2, c3))
  dimnames(c.ccc) <- list(c.names, ccclab)
  c.out           <- data.frame(c.names, 
                                round(1000 * c.mass, 0), 
                                round(1000 * c.qlt, 0),
                                round(1000 * c.inr, 0), 
                                round(1000 * c.ccc, 0))
  dimnames(c.out) <- list(as.character(1:length(c.names)),
                          c("name", "mass", " qlt", " inr", ccclab))

 # scree plot:
  sev.0 <- round(100*obj$inertia.et, 1)
  if (scree) {
   # if (obj$lambda=="indicator"){
   #   values     <- obj$sv^2
   #   } else {
      values     <- obj$sv
   #   }
    values2    <- round(100*values/sum(values), 1)
    scree.out  <- cbind(1:length(obj$sv), round(values, 6), values2, round(cumsum(100*values/sum(values)), 1))
   # sev.0      <- round(100*sum(values/sum(values)), 1)
    if (obj$lambda == "adjusted"){
      values     <- round(obj$sv, 6)
      values2    <- round(100*obj$inertia.e, 1)
      values3    <- rep(NA, length(values))
      scree.out  <- cbind(1:length(obj$sv), round(values, 6), values2, values3)
      }
    if (obj$lambda == "JCA"){
      values     <- round(obj$sv, 6)
      values2    <- rep(NA, length(values))
      values3    <- rep(NA, length(values))
      scree.out  <- cbind(1:length(obj$sv), round(values, 6), values2, values3)
      }
   #   scree.out  <- cbind(1:length(obj$sv), round(values, 6), values2, round(cumsum(100*values/sum(values)), 1))
    } else {
    scree.out <- NA
    }

 # output:
  out <- list(scree   = scree.out,
              rows    = r.out,
              columns = c.out, 
              sev     = sev.0, 
              JCA     = obj$JCA.iter,
              tin     = obj$inertia.t, 
              JCA.nd  = obj$nd, 
              JCA.ind = sum(diag(obj$subinertia)),
              JCA.nit = obj$JCA.iter[1],
              JCA.eps = obj$JCA.iter[2])
  class(out) <- "summary.mjca"
  return(out)
  }

################################################################################

