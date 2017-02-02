#' Fit a Bayesian epiallele model
#'
#' Fits an epiallele model to DNAm sequencing data.
#' @param Y is an N by d numeric matrix of 0's and 1's where each row is a
#' sequencing read and each column corresponds to a CpG site. Missing data can
#' be represented as NA.
#' @usage fit = bed_fit(Y)
#' @return A list structure where rho is the optimal number of epialleles
#' and AIC is an N-dimensional vector where element q contains the Akaike
#' information criterion score for a model with q epialleles. The model
#' field contains fitted models for different values of q. Each model has a
#' vector W which indicates which epiallele each observed read corresponds to,
#' a matrix X where each row is an inferred epiallale, matches and
#' mismatches are simply the total number of mis/matches, and beta is
#' the inferred noise hyperparameter.
#'
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
#' # Add some noise by randomly regenerating 20 percent of values
#' Y[sample(1:800,80,replace=FALSE)] <- round(runif(80))
#'
#' # Make a few data missing
#' Y[sample(1:800,5,replace=FALSE)] <- NA
#'
#' # Infer which epialleles are present
#' fit <- bed_fit(Y)
#'
#' # Plot AIC score
#' plot(fit$AIC[1:10], type='b', xlab='Q = number of epialleles',
#'     ylab='AIC score')
#'
#' # The fit$rho field contains the optimal number of epialleles
#' cat(' Optimal number of epialleles =', fit$rho,'\n')
#'
#' #Marginalise over W posterior
#' w.margin <- bed_marginal(Y, fit$model[[fit$rho]]$X, fit$model[[fit$rho]]$W)
#'
#' # Summing over the columns gives the total proportion of reads that are
#' # attributed to each epiallele
#' p <- apply(w.margin, 2, FUN='sum')/100
#' cat(' Inferred Epiallele [1] =',paste(fit$model[[3]]$X[1,],collapse=""),
#'     ', Proportion =', p[1],'\n',
#'     'Inferred Epiallele [2] =',paste(fit$model[[3]]$X[2,],collapse=""),
#'     ', Proportion =', p[2],'\n',
#'     'Inferred Epiallele [3] =',paste(fit$model[[3]]$X[3,],collapse="")
#'     ,', Proportion =', p[3])
#'
#' @seealso \code{\link{bed_marginal}}

bed_fit <- function(Y){

   # Total number of reads
   N <- nrow(Y)
   # Preallocate a list for the output
   results <- list(model=vector('list',N), AIC=rep(NA,N))

   H <- hamming(Y)
   fit <- hclust(as.dist(H), method = 'ward.D')

   # ---------- Begin loop over N -------- #
   for (q in 1:N){

      # BREAK if there's only one barcode pair observed
      if (N==1){
         results$rho <- 1
         break
      }

      # Fit model using LL function
      results$model[[q]] <- MAPxw(Y, q, fit)
      results$AIC[q] <- results$model[[q]]$AIC
      results$rho <- which.min(results$AIC) # Current best estimate for corrected read count

      # BREAK if there are no mismatches and q=1. This means the model has achieved a
      # perfect fit so we can stop.
      # BREAK if there are no mismatches.
      if (results$model[[q]]$mismatches==0){
         results$rho <- which.min(results$AIC[1:q])
         break
      }

      # Search at least until q=10, then check to see if we've hit the minimum AIC score.
      # The AIC is non monotonic so don't search in the last five values of q.
      if ((q>=10) & (which.min(results$AIC[1:q])<(q-5))) break
   }
   # ---------- End loop over N -------- #

   return(results)
}
