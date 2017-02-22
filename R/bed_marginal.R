#' Marginal epiallele distribution
#'
#' Marginalises over the W parameter in a fitted epiallele model
#' @param Y an N by d numeric matrix of observed DNAm reads
#' @param X a Q by d numeric matrix of inferred epialleles
#' @param W an N by 1 numeric vector assigning reads to epiallels
#' @return An N by d numerc matrix where where element [i,mu] is the probability
#' that read i corresponds to epiallele mu.
#' @examples
#' set.seed(1314)
#'
#' # Generate 50 reads that are fully methylated
#' Y1 <- matrix(rep(c(1,1,1,1,1,1,1,1),50), nrow=50, ncol=8)
#'
#' # Generate 30 reads that are half methylated
#' Y2 <- matrix(rep(c(1,1,1,1,0,0,0,0),30),nrow=30, ncol=8, byrow=TRUE)
#'
#' # Generate 20 reads that are fully unmethylated
#' Y3 <- matrix(rep(c(0,0,0,0,0,0,0,0),20),nrow=20,ncol=8)
#'
#' # Combine them
#' Y <- rbind(Y1,Y2,Y3)
#'
#' # Add some noise by randomly regenerating 5 percent of values
#' Y[sample(1:800,40,replace=FALSE)] <- round(runif(40))
#'
#' # Make a few data missing
#' Y[sample(1:800,5,replace=FALSE)] <- NA
#'
#' # Infer which epialleles are present
#' fit <- bed_fit(Y)
#'
#' # Plot AIC score
#' plot(fit$AIC[1:10], type='b', xlab='Q = number of epialleles',
#'     ylab='AIC score', ylim=c(200,300))
#'
#' # The fit$rho field contains the optimal number of epialleles
#' cat(' Optimal number of epialleles =', fit$rho,'\n')
#'
#' #Marginalise over W posterior
#' w.margin <- bed_marginal(Y, fit$model[[fit$rho]]$X, fit$model[[fit$rho]]$W)
#'
#' # Summing over the columns gives the total proportion of reads that are
#' # attributed to each epiallele. Just over 4 percent of reads are attributed to
#' # a spurious epiallele (due to noise)
#' p <- apply(w.margin, 2, FUN='sum')/100
#' cat(' Inferred Epiallele [1] =',paste(fit$model[[4]]$X[1,],collapse=""),
#'     ', Proportion =', p[1],'\n',
#'     'Inferred Epiallele [2] =',paste(fit$model[[4]]$X[2,],collapse=""),
#'     ', Proportion =', p[2],'\n',
#'     'Inferred Epiallele [3] =',paste(fit$model[[4]]$X[3,],collapse=""),
#'     ', Proportion =', p[3],'\n',
#'     'Inferred Epiallele [4] =',paste(fit$model[[4]]$X[4,],collapse="")
#'     ,', Proportion =', p[4])
#'
#' @seealso \code{\link{bed_fit}}

bed_marginal <- function(Y,X,W){

   Q <- nrow(X)
   N <- nrow(Y)

   N.plus <- numeric(Q)
   N.minus <- numeric(Q)

   for (q in seq(1,Q)){
      N.q <- sum(W==q)
      term1 <- colSums(Y[W==q,,drop=FALSE], na.rm = TRUE)
      term3 <- colSums(!is.na(Y[W==q,,drop=FALSE]), na.rm = TRUE)
      term2 <- term3 - term1
      N.plus[q] <- sum(term1[X[q,]==1]) + sum(term2[X[q,]==0])
      N.minus[q] <- sum(term3) - N.plus[q]
   }

   matches <- sum(N.plus)
   mismatches <- sum(N.minus)

   # change w
   w.margin <- matrix(NA,N,Q)

   for (i in 1:N){
      w <- numeric(Q)
      for (mu in 1:Q){

         # subtract MAP matches and add current matches
         w.matches <- matches - sum(Y[i,] == X[W[i],], na.rm = TRUE) +
            sum(Y[i,] == X[mu,], na.rm = TRUE)
         w.mismatches <- mismatches - sum(Y[i,] != X[W[i],], na.rm = TRUE) +
            sum(Y[i,] != X[mu,], na.rm = TRUE)

         beta <- max(mismatches,1)/(matches+mismatches)

         w[mu] <- w.matches*log(1-beta) + w.mismatches*log(beta)
      }
      w.margin[i,] <- exp(w - mean(w))/sum(exp(w - mean(w)))
   }

   return(w.margin)
}
