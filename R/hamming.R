hamming <- function(Y){
   
   d <- ncol(Y)
#       N <- nrow(Y)
#       
#       H <- matrix(0,N,N)
#       
#       for (i in seq(1,N)){
#          for (j in seq(1,N)){
#             H[i,j] <- sum(Y[i,] != Y[j,], na.rm = TRUE)
#          }
#       }
   
   B <- Y
   B[is.na(Y)] <- 0
   B1 <- (B==1)
   
   C <- Y
   C[is.na(Y)] <- 1
   C0 <- (C==0)
   
   H <- (B1%*%t(C0) + C0%*%t(B1))/d
   
   return(H)
   
}
