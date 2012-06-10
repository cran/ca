
################################################################################
# 
# print.mjca: (Generic) function for printing 'ca'-objects
#
#      Input: x   The 'mjca' object which should be printed
#
#     Output:     'mjca' object neatly formatted
# 
################################################################################


print.mjca <- function(x, ...)
  {
  obj <- x
  nd  <- obj$nd
  if (is.na(nd)){
    nd <- 2
    } else {
    if (nd > length(obj$sv)) nd <- length(obj$sv)
    }

 # Eigenvalues:
  Dimension  <- 1:length(obj$sv)
  Value      <- round(obj$sv^2, 6)
  Percentage <- paste(as.character(round(100*obj$inertia.e, 2)), "%", sep = "")

  tmp <- rbind(Value = as.character(Value), 
               Percentage = as.character(Percentage))
  dimnames(tmp)[[2]] <- Dimension  
  Eigenvalues <- tmp

 # Row Profiles:
 # tmp <- rbind(obj$rowmass, obj$rowdist, obj$rowinertia, t(obj$rowcoord[,1:nd]))
 # tmpnames <- obj$rownames
 # if (!is.na(obj$rowsup[1]))
 #   {
 #   tmpnames[obj$rowsup] <- paste(tmpnames[obj$rowsup],"(*)")
 #   }
 # dimnames(tmp)[[2]] <- tmpnames
 # dn <- paste("Dim.", 1:nd)
 # dimnames(tmp)[[1]] <- c("Mass", "ChiDist", "Inertia", dn)
 # Row.profiles <- tmp

 # Column Profiles:
  tmp <- rbind(obj$colmass, obj$coldist, obj$colinertia, t(obj$colcoord[,1:nd]))
  tmpnames <- obj$levelnames
  if (!is.na(obj$colsup[1]))
    {
    tmpnames[obj$colsup] <- paste(tmpnames[obj$colsup],"(*)",sep="")
    }
# BCN 2009_11:
 # if (!is.na(obj$subsetcol[1])){
 #   dimnames(tmp)[[2]] <- tmpnames[obj$subsetcol]
 #   } else {
    dimnames(tmp)[[2]] <- tmpnames
 #   }
  dn <- paste("Dim.", 1:nd)
  dimnames(tmp)[[1]] <- c("Mass", "ChiDist", "Inertia", dn)
  Column.profiles <- tmp

cat("\n Eigenvalues:\n")
print.table(Eigenvalues, width = 4)

#cat("\n\n Rows:\n")
#print(round(Row.profiles, 6))

cat("\n\n Columns:\n")
print(round(Column.profiles, 6))

}


################################################################################

