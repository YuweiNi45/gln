#'
#' Provide basic summary statistics
#'
#' Generic function to provide basic summary statistics analysis for numeric vectors
#'
#' @param y An object for which a method has been defined, or a numeric vector containing the values whose length, mean, median, quartiles are to be computed.
#' @return The method returns \code{y} containing the length, mean, median, 25% and 75% quartiles of the object. Missing values are removed by default.
#' @seealso \code{\link[base]{length}}, \code{\link[base]{mean}},\code{\link[stats]{quantile}}
#' @references Becker, R.A., Chambers, J.M., and Wilks, A.R.(1988) The New S Language. Wadsworth & Brooks/Cole.
#' @examples
#' x<-c(0:10,50)
#' data_sum(x)
#'@export

data_sum<- function(y){
  out <- data.frame(n = length(y),
                    n.complete = length(y[!is.na(y)]),
                    mean = mean(y, na.rm=T),
                    sd = sd(y, na.rm=T),
                    median = median(y, na.rm=T),
                    lower_quartile = quantile(y, 0.25, na.rm=T),
                    upper_quartile = quantile(y, 0.75, na.rm=T))
  return(out)
}
