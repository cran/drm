"kronecker.drm" <-
  function (y, nclass = length(levels(factor(y))), dep = "I") 
{
  nrep <- nre <- length(y)
  mis <- which(is.na(y) == TRUE)
  if (length(grep("M", dep)) > 0) {
    if (length(mis) > 0) {
      if (any(!is.na(y[min(mis):nrep]))) 
        stop("Only monotone missing patterns for M-structure currently allowed")
      y <- y[seq(min(mis) - 1)]
      nre <- length(y)
    }
    if (length(grep("M2", dep)) > 0) {
      ma <- matrix(0, ncol = nrep - max(nre, 3), nrow = nclass^3)
      ma[1, ] <- 1
      if (nre == 1) 
        w <- cbind(matrix(kronecker.drm(c(y, NA, NA), 
                                       nclass = nclass), ncol = 1), ma)
      else {
        if (nre == 2) 
          w <- cbind(matrix(kronecker.drm(c(y, NA), nclass = nclass), 
                            ncol = 1), ma)
        else {
          w <- cbind(sapply(1:(nre - 2), function(j, 
                                                  y, nclass) kronecker.drm(y[j:(j + 2)], nclass = nclass), 
                            nclass = nclass, y = y), ma)
        }
      }
      if (nrep > 3) {
        if (nre < 4) 
          w <- cbind(w, ma)
        else {
          w <- cbind(w, rbind(sapply(2:(nre - 2), function(j, 
                                                           y, nclass) kronecker.drm(y[j:(j + 1)], nclass = nclass), 
                                     nclass = nclass, y = y), matrix(0, ncol = nre - 
                                                        3, nrow = (nclass^3 - nclass^2))), ma)
        }
      }
      w
    }
    else {
      ma <- matrix(0, ncol = nrep - max(nre, 2), nrow = nclass^2)
      ma[1, ] <- 1
      if (nre == 1) 
        w <- cbind(matrix(kronecker.drm(c(y, NA), nclass = nclass), 
                          ncol = 1), ma)
      else {
        w <- cbind(sapply(1:(nre - 1), function(j, y, 
                                                nclass) kronecker.drm(y[j:(j + 1)], nclass = nclass), 
                          nclass = nclass, y = y), ma)
      }
      if (nrep > 2) {
        if (nre < 3) 
          w <- cbind(w, ma)
        else {
          w <- cbind(w, rbind(sapply(2:(nre - 1), function(j, 
                                                           y, nclass) kronecker.drm(y[j], nclass = nclass), 
                                     nclass = nclass, y = y), matrix(0, ncol = nre - 
                                                        2, nrow = (nclass^2 - nclass))), ma)
        }
      }
      w
    }
  }
  else {
    if (length(mis) > 0) {
      y <- y[seq(min(mis) - 1)]
      w <- kronecker.drm(y, nclass = nclass)
      c(w, rep(0, (nclass^nrep) - (nclass^length(y))))
    }
    else {
      ymat <- matrix(0, ncol = nrep, nrow = nclass)
      ymat[, y == 1] <- c(1, rep(-1, (nclass - 1)))
      for (i in 2:nclass) {
        ymat[, y == i] <- c(rep(0, (i - 1)), 1, rep(0, 
                 (nclass - i)))
      }
      if (nrep == 1) 
        w <- ymat
      else (w <- eval(parse(text = paste(paste(rep("c(", 
                              nrep - 1), collapse = ""),
                              paste("ymat[,1]",
                                    paste("%*% t(ymat[,", 2:nrep, "]))", collapse = ""))))))
      w
    }
  }
}

