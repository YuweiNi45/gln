#' create tables providing for post hoc analysis
#'
#' To detect between which groups the difference is, the package produces different tables for post hoc analysis depending on one factor or two factors in the data.
#'
#' @param x  a vector that contains different levels whose summary results are wanted for comparison.
#' @param y  a numeric vector of data values
#' @param dat  a data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables to be summarized.
#' @param method	 the method to be used; for General linear hypotheses and multiple comparisons for parametric models; must be one of " tukey " (default), "dunnett", "general", "holm" or "none".
#' @param K	a matrix that contains the differences between different levels that are interested in
#' @param digits   the digits of the data values (default is 3)
#' @param ...	further arguments to be passed to or from methods.
#' @return The posthoc_test returns a table containing p-values for selected method, must be one of " tukey " (default), "dunnett", "general", "holm" or "none".  Missing values are removed by default.
#' @seealso \code{\link[multcomp]{glht}}
#' @references Becker, R.A., Chambers, J.M., and Wilks, A.R.(1988) The New S Language. Wadsworth & Brooks/Cole.
#' @examples
#'x <- c("A","B","C","D","A","B","C","D","A","B","C","D")
#'y <- c(1,3,2,3, 1,3,2,3, 1,3,2,3)
#'data1 <- as.data.frame(cbind(x,as.numeric(y)))
#'data1$V2 <- as.numeric(data1$V2)
#'posthoc_test(data1$x, data1$V2, data1)
#'
#'
#'@export
posthoc_test <- function(x, y, dat, method, K, digits=3, ...){

  x <- as.factor(x)

  test<-function(x){
    if(x < 0.001) x.txt<- "<0.001"
    else if(x<0.05) x.txt<-"<0.05"
    else x.txt<- format(x,digits = 2)
    x.txt
  }

  if(!missing(K)){

    out <- lm(y ~ x, data = dat, ...)

    ##parametric

    dunnett <- summary(glht(out, linfct = mcp(x = "Dunnett")))

    tukey <- summary(glht(out, linfct = mcp(x = "Tukey")))

    general <- summary(glht(out, linfct=mcp(x=K)))

    holm <- summary(glht(out, linfct=mcp(x=K)), test=adjusted(type="holm"))

    non_adjusted <- summary(glht(out, linfct=mcp(x=K)), test=adjusted(type="none"))

    ##dunnett

    dunnett_p <- dunnett$test$pvalues
    dunnett_p1 <- unname(dunnett_p, force = FALSE)

    np1 <- length(dunnett_p1)

    p.dunnett<-numeric(np1)

    for (i in 1:np1) {
      p.dunnett[i] <- test(dunnett_p1[i])
    }

    ##tukey

    tukey_p <- tukey$test$pvalues
    tukey_p1 <- unname(tukey_p, force = FALSE)

    np2 <- length(tukey_p1)

    p.tukey<-numeric(np2)

    for (i in 1:np2) {
      p.tukey[i] <- test(tukey_p1[i])
    }

    ##general

    general_p <- general$test$pvalues
    general_p1 <- unname(general_p, force = FALSE)

    np3 <- length(general_p1)

    p.general<-numeric(np3)

    for (i in 1:np3) {
      p.general[i] <- test(general_p1[i])
    }

    ##holm
    holm_p <- holm$test$pvalues
    holm_p1 <- unname(holm_p, force = FALSE)

    np4 <- length(holm_p1)

    p.holm<-numeric(np4)

    for (i in 1:np4) {
      p.holm[i] <- test(holm_p1[i])
    }

    #non adjusted
    non_adjusted_p <- non_adjusted$test$pvalues
    non_adjusted_p1 <- unname(non_adjusted_p, force = FALSE)

    np5 <- length(non_adjusted_p1)

    p.nonad<-numeric(np5)

    for (i in 1:np5) {
      p.nonad[i] <- test(non_adjusted_p1[i])
    }

    if(!missing(method)){
      if (method == "dunnett"){
        out <- data.frame(pm.p.dunnett=p.dunnett,
                          pm.p=p.nonad
        )
      } else if (method == "tukey"){
        out <- data.frame(pm.p.tukey=p.tukey,
                          pm.p=p.nonad
        )
      } else if (method == "general"){
        out <- data.frame(pm.p.general=p.general,
                          pm.p=p.nonad
        )
      } else if (method == "holm"){
        out <- data.frame(pm.p.holm=p.holm,
                          pm.p=p.nonad
        )
      }
    } else {
      out <- data.frame(pm.p.tukey=p.tukey,
                        pm.p=p.nonad
      )
    }

    rownames(out) <- rownames(K)
  }
  else {
    #don't have k matrix

    out <- lm(y ~ x, data = dat, ...)

    dunnett <- summary(glht(out, linfct = mcp(x = "Dunnett")))

    tukey <- summary(glht(out, linfct = mcp(x = "Tukey")))

    ##dunnett

    dunnett_p <- dunnett$test$pvalues
    dunnett_p1 <- unname(dunnett_p, force = FALSE)

    np1 <- length(dunnett_p1)

    p.dunnett<-numeric(np1)

    for (i in 1:np1) {
      p.dunnett[i] <- test(dunnett_p1[i])
    }

    ##tukey

    tukey_p <- tukey$test$pvalues
    tukey_p1 <- unname(tukey_p, force = FALSE)

    np2 <- length(tukey_p1)

    p.tukey<-numeric(np2)

    for (i in 1:np2) {
      p.tukey[i] <- test(tukey_p1[i])
    }

    if(!missing(method)){
      if (method == "dunnett"){
        out <- data.frame(pm.p.dunnett=p.dunnett
        )
        rownames(out) <- rownames(dunnett$linfct)
      } else if (method == "tukey"){
        out <- data.frame(pm.p.tukey=p.tukey
        )
        rownames(out) <- rownames(tukey$linfct)
      }
    }else {
      out <- data.frame(pm.p.tukey=p.tukey
      )
      rownames(out) <- rownames(tukey$linfct)
    }
  }

  return(out)

}
