#' create tables providing for post hoc analysis
#'
#' To detect between which groups the difference is, the package produces different tables for post hoc analysis depending on one factor or two factors in the data.
#'
#'@param data	a data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables to be summarized
#'@param y	a numeric vector of data value
#'@param x	a vector that contains different levels whose summary results are wanted for comparison
#'@param x1        a vector that contains different levels whose summary results are wanted for comparison. Only necessary for two-factor analysis
#'@param K	a matrix that contains the differences between different levels that are interested in
#'@param method	         the method to be used; for General linear hypotheses and multiple comparisons for parametric models; must be one of " tukey " (default), "dunnett", "general", "holm" or "none"
#'@param model        the parametric models to be used; for linear regression and multiple comparisons; must be one of "cellmean"(default) or "maineffect"
#'@param digits	the digits of the data values (default is 3)
#'@param ...	further arguments to be passed to or from methods
#'
#' @return The posthoc_test returns a table containing p-values for selected method, must be one of " tukey " (default), "dunnett", "general", "holm" or "none", and selected model, must be one of "cellmean" or "maineffect". Missing values are removed by default
#' @seealso \code{\link[multcomp]{glht}}
#' @references Becker, R.A., Chambers, J.M., and Wilks, A.R.(1988) The New S Language. Wadsworth & Brooks/Cole.
#' @examples
#'x <- c("A","B","C","D","A","B","C","D","A","B","C","D")
#'y <- c(6,4,2,7, 3,10,9,3, 1,2,8,1)
#'data1 <- as.data.frame(cbind(x,as.numeric(y)))
#'data1$V2 <- as.numeric(data1$V2)
#'posthoc_test(data1, data1$V2, data1$x)
#'
#'
#'@export
posthoc_test <-
  function(data,
           y,
           x,
           x1,
           K,
           method = c("tukey", "dunnett", "general", "holm", "notadjusted"),
           model = c("cellmean", "maineffect"),
           digits = 3,
           ...) {
    x <- as.factor(x)

    method <- match.arg(method)

    model <- match.arg(model)

    test <- function(x) {
      if (x < 0.001)
        x.txt <- "<0.001"
      else if (x < 0.05)
        x.txt <- "<0.05"
      else
        x.txt <- format(x, digits = 2)
      x.txt
    }

    ##main effect model
    if (!missing(x1) | model == "maineffect") {
      x1 <- factor(x1)

      ##linear regression
      out1 <- lm(y ~ x * x1 - 1, data = data, ...)


      if (!missing(K)) {
        ##with K matrix

        general <- summary(glht(out1, linfct = K))

        holm <-
          summary(glht(out1, linfct = K), test = adjusted(type = "holm"))

        non_adjusted <-
          summary(glht(out1, linfct = K), test = adjusted(type = "none"))

        ##getting the p-values for different methods

        ##general

        general_p <- general$test$pvalues
        general_p1 <- unname(general_p, force = FALSE)

        np3 <- length(general_p1)

        p.general <- numeric(np3)

        for (i in 1:np3) {
          p.general[i] <- test(general_p1[i])
        }

        ##holm
        holm_p <- holm$test$pvalues
        holm_p1 <- unname(holm_p, force = FALSE)

        np4 <- length(holm_p1)

        p.holm <- numeric(np4)

        for (i in 1:np4) {
          p.holm[i] <- test(holm_p1[i])
        }

        #non adjusted
        non_adjusted_p <- non_adjusted$test$pvalues
        non_adjusted_p1 <- unname(non_adjusted_p, force = FALSE)

        np5 <- length(non_adjusted_p1)

        p.nonad <- numeric(np5)

        for (i in 1:np5) {
          p.nonad[i] <- test(non_adjusted_p1[i])
        }

        #non adjusted
        non_adjusted_p <- non_adjusted$test$pvalues
        non_adjusted_p1 <- unname(non_adjusted_p, force = FALSE)

        np5 <- length(non_adjusted_p1)

        p.nonad <- numeric(np5)

        for (i in 1:np5) {
          p.nonad[i] <- test(non_adjusted_p1[i])
        }

        ##choosing the method

        if (!missing(method)) {
          if (method == "dunnett") {
            ##Dunnett

            K1 <- glht(out1, linfct = mcp(x = "Dunnett"))$linfct
            K2 <- glht(out1, linfct = mcp(x1 = "Dunnett"))$linfct

            dunnett <- summary(glht(out1, linfct = rbind(K1, K2)))

            ##preparing matrixes

            K5 <- as.matrix(dunnett$linfct)

            ##dunnett's nonadjusted
            non_adjusted_1 <-
              summary(glht(out1, linfct = K5), test = adjusted(type = "none"))

            ##dunnett

            dunnett_p <- dunnett$test$pvalues
            dunnett_p1 <- unname(dunnett_p, force = FALSE)

            np1 <- length(dunnett_p1)

            p.dunnett <- numeric(np1)

            for (i in 1:np1) {
              p.dunnett[i] <- test(dunnett_p1[i])
            }

            #non adjusted_dunnett
            non_adjusted_1_p <- non_adjusted_1$test$pvalues
            non_adjusted_1_p1 <-
              unname(non_adjusted_1_p, force = FALSE)

            p.nonad_1 <- numeric(np1)

            for (i in 1:(np1)) {
              p.nonad_1[i] <- test(non_adjusted_1_p1[i])
            }

            out <- data.frame(pm.p.dunnett = p.dunnett,
                              pm.p = p.nonad_1)
            rownames(out) <- rownames(dunnett$linfct)
          } else if (method == "general") {
            out <- data.frame(pm.p.general = p.general,
                              pm.p = p.nonad)
            rownames(out) <- rownames(K)
          } else if (method == "holm") {
            out <- data.frame(pm.p.holm = p.holm,
                              pm.p = p.nonad)
            rownames(out) <- rownames(K)
          } else if (method == "notadjusted") {
            out <- data.frame(pm.p = p.nonad)
            rownames(out) <- rownames(K)
          } else {
            ##Tukey

            K3 <- glht(out1, linfct = mcp(x = "Tukey"))$linfct
            K4 <- glht(out1, linfct = mcp(x1 = "Tukey"))$linfct

            tukey <- summary(glht(out1, linfct = rbind(K3, K4)))

            ##preparing matrixes

            K6 <- as.matrix(tukey$linfct)

            ##tukey's nonadjusted
            non_adjusted_2 <-
              summary(glht(out1, linfct = K6), test = adjusted(type = "none"))


            ##tukey

            tukey_p <- tukey$test$pvalues
            tukey_p1 <- unname(tukey_p, force = FALSE)

            np2 <- length(tukey_p1)

            p.tukey <- numeric(np2)

            for (i in 1:np2) {
              p.tukey[i] <- test(tukey_p1[i])
            }


            #non adjusted_tukey
            non_adjusted_2_p <- non_adjusted_2$test$pvalues
            non_adjusted_2_p1 <-
              unname(non_adjusted_2_p, force = FALSE)

            p.nonad_2 <- numeric(np2)

            for (i in 1:np2) {
              p.nonad_2[i] <- test(non_adjusted_2_p1[i])
            }


            out <- data.frame(pm.p.tukey = p.tukey,
                              pm.p = p.nonad_2)
            rownames(out) <- rownames(tukey$linfct)
          }
        } else {
          ##Tukey

          K3 <- glht(out1, linfct = mcp(x = "Tukey"))$linfct
          K4 <- glht(out1, linfct = mcp(x1 = "Tukey"))$linfct

          tukey <- summary(glht(out1, linfct = rbind(K3, K4)))

          ##preparing matrixes

          K6 <- as.matrix(tukey$linfct)

          ##tukey's nonadjusted
          non_adjusted_2 <-
            summary(glht(out1, linfct = K6), test = adjusted(type = "none"))


          ##tukey

          tukey_p <- tukey$test$pvalues
          tukey_p1 <- unname(tukey_p, force = FALSE)

          np2 <- length(tukey_p1)

          p.tukey <- numeric(np2)

          for (i in 1:np2) {
            p.tukey[i] <- test(tukey_p1[i])
          }


          #non adjusted_tukey
          non_adjusted_2_p <- non_adjusted_2$test$pvalues
          non_adjusted_2_p1 <-
            unname(non_adjusted_2_p, force = FALSE)

          p.nonad_2 <- numeric(np2)

          for (i in 1:np2) {
            p.nonad_2[i] <- test(non_adjusted_2_p1[i])
          }

          out <- data.frame(pm.p.tukey = p.tukey,
                            pm.p = p.nonad_2)
          rownames(out) <- rownames(tukey$linfct)
        }



      } else {
        ## missing k matrix

        if (!missing(method)) {
          if (method == "dunnett") {
            ##Dunnett

            K1 <- glht(out1, linfct = mcp(x = "Dunnett"))$linfct
            K2 <- glht(out1, linfct = mcp(x1 = "Dunnett"))$linfct

            dunnett <- summary(glht(out1, linfct = rbind(K1, K2)))

            ##preparing matrixes

            K5 <- as.matrix(dunnett$linfct)

            ##dunnett's nonadjusted
            non_adjusted_1 <-
              summary(glht(out1, linfct = K5), test = adjusted(type = "none"))

            ##dunnett

            dunnett_p <- dunnett$test$pvalues
            dunnett_p1 <- unname(dunnett_p, force = FALSE)

            np1 <- length(dunnett_p1)

            p.dunnett <- numeric(np1)

            for (i in 1:np1) {
              p.dunnett[i] <- test(dunnett_p1[i])
            }

            #non adjusted_dunnett
            non_adjusted_1_p <- non_adjusted_1$test$pvalues
            non_adjusted_1_p1 <-
              unname(non_adjusted_1_p, force = FALSE)

            p.nonad_1 <- numeric(np1)

            for (i in 1:(np1)) {
              p.nonad_1[i] <- test(non_adjusted_1_p1[i])
            }

            out <- data.frame(pm.p.dunnett = p.dunnett,
                              pm.p = p.nonad_1)
            rownames(out) <- rownames(dunnett$linfct)
          } else if (method == "tukey") {
            ##Tukey

            K3 <- glht(out1, linfct = mcp(x = "Tukey"))$linfct
            K4 <- glht(out1, linfct = mcp(x1 = "Tukey"))$linfct

            tukey <- summary(glht(out1, linfct = rbind(K3, K4)))

            ##preparing matrixes

            K6 <- as.matrix(tukey$linfct)

            ##tukey's nonadjusted
            non_adjusted_2 <-
              summary(glht(out1, linfct = K6), test = adjusted(type = "none"))


            ##tukey

            tukey_p <- tukey$test$pvalues
            tukey_p1 <- unname(tukey_p, force = FALSE)

            np2 <- length(tukey_p1)

            p.tukey <- numeric(np2)

            for (i in 1:np2) {
              p.tukey[i] <- test(tukey_p1[i])
            }


            #non adjusted_tukey
            non_adjusted_2_p <- non_adjusted_2$test$pvalues
            non_adjusted_2_p1 <-
              unname(non_adjusted_2_p, force = FALSE)

            p.nonad_2 <- numeric(np2)

            for (i in 1:np2) {
              p.nonad_2[i] <- test(non_adjusted_2_p1[i])
            }

            out <- data.frame(pm.p.tukey = p.tukey,
                              pm.p = p.nonad_2)
            rownames(out) <- rownames(tukey$linfct)
          } else {
            out <- "Missing K matrix"
          }
        } else {
          ##Tukey

          K3 <- glht(out1, linfct = mcp(x = "Tukey"))$linfct
          K4 <- glht(out1, linfct = mcp(x1 = "Tukey"))$linfct

          tukey <- summary(glht(out1, linfct = rbind(K3, K4)))

          ##preparing matrixes

          K6 <- as.matrix(tukey$linfct)

          ##tukey's nonadjusted
          non_adjusted_2 <-
            summary(glht(out1, linfct = K6), test = adjusted(type = "none"))


          ##tukey

          tukey_p <- tukey$test$pvalues
          tukey_p1 <- unname(tukey_p, force = FALSE)

          np2 <- length(tukey_p1)

          p.tukey <- numeric(np2)

          for (i in 1:np2) {
            p.tukey[i] <- test(tukey_p1[i])
          }


          #non adjusted_tukey
          non_adjusted_2_p <- non_adjusted_2$test$pvalues
          non_adjusted_2_p1 <-
            unname(non_adjusted_2_p, force = FALSE)

          p.nonad_2 <- numeric(np2)

          for (i in 1:np2) {
            p.nonad_2[i] <- test(non_adjusted_2_p1[i])
          }

          out <- data.frame(pm.p.tukey = p.tukey,
                            pm.p = p.nonad_2)
          rownames(out) <- rownames(tukey$linfct)
        }

      }


    } else {
      ## cell mean model

      out1 <- lm(y ~ x - 1, data = data, ...)

      dunnett <- summary(glht(out1, linfct = mcp(x = "Dunnett")))

      tukey <- summary(glht(out1, linfct = mcp(x = "Tukey")))

      ##preparing matrixes

      K5 <- as.matrix(dunnett$linfct)

      K6 <- as.matrix(tukey$linfct)

      ##dunnett's nonadjusted
      non_adjusted_1 <-
        summary(glht(out1, linfct = K5), test = adjusted(type = "none"))

      ##tukey's nonadjusted
      non_adjusted_2 <-
        summary(glht(out1, linfct = K6), test = adjusted(type = "none"))


      ##dunnett

      dunnett_p <- dunnett$test$pvalues
      dunnett_p1 <- unname(dunnett_p, force = FALSE)

      np1 <- length(dunnett_p1)

      p.dunnett <- numeric(np1)

      for (i in 1:np1) {
        p.dunnett[i] <- test(dunnett_p1[i])
      }

      ##tukey

      tukey_p <- tukey$test$pvalues
      tukey_p1 <- unname(tukey_p, force = FALSE)

      np2 <- length(tukey_p1)

      p.tukey <- numeric(np2)

      for (i in 1:np2) {
        p.tukey[i] <- test(tukey_p1[i])
      }

      #non adjusted_dunnett
      non_adjusted_1_p <- non_adjusted_1$test$pvalues
      non_adjusted_1_p1 <- unname(non_adjusted_1_p, force = FALSE)

      np6 <- length(non_adjusted_1_p1)

      p.nonad_1 <- numeric(np6)

      for (i in 1:np6) {
        p.nonad_1[i] <- test(non_adjusted_1_p1[i])
      }

      #non adjusted_tukey
      non_adjusted_2_p <- non_adjusted_2$test$pvalues
      non_adjusted_2_p1 <- unname(non_adjusted_2_p, force = FALSE)

      np7 <- length(non_adjusted_2_p1)

      p.nonad_2 <- numeric(np7)

      for (i in 1:np7) {
        p.nonad_2[i] <- test(non_adjusted_2_p1[i])
      }

      if (!missing(K)) {
        ##performing different methods

        general <- summary(glht(out1, linfct = mcp(x = K)))

        holm <-
          summary(glht(out1, linfct = mcp(x = K)), test = adjusted(type = "holm"))

        non_adjusted <-
          summary(glht(out1, linfct = mcp(x = K)), test = adjusted(type = "none"))


        ##getting the pvalues

        ##general

        general_p <- general$test$pvalues
        general_p1 <- unname(general_p, force = FALSE)

        np3 <- length(general_p1)

        p.general <- numeric(np3)

        for (i in 1:np3) {
          p.general[i] <- test(general_p1[i])
        }

        ##holm
        holm_p <- holm$test$pvalues
        holm_p1 <- unname(holm_p, force = FALSE)

        np4 <- length(holm_p1)

        p.holm <- numeric(np4)

        for (i in 1:np4) {
          p.holm[i] <- test(holm_p1[i])
        }

        #non adjusted
        non_adjusted_p <- non_adjusted$test$pvalues
        non_adjusted_p1 <- unname(non_adjusted_p, force = FALSE)

        np5 <- length(non_adjusted_p1)

        p.nonad <- numeric(np5)

        for (i in 1:np5) {
          p.nonad[i] <- test(non_adjusted_p1[i])
        }


        ##choosing the method
        if (!missing(method)) {
          if (method == "dunnett") {
            out <- data.frame(pm.p.dunnett = p.dunnett,
                              pm.p = p.nonad_1)
            rownames(out) <- rownames(dunnett$linfct)
          } else if (method == "general") {
            out <- data.frame(pm.p.general = p.general,
                              pm.p = p.nonad)
            rownames(out) <- rownames(K)
          } else if (method == "holm") {
            out <- data.frame(pm.p.holm = p.holm,
                              pm.p = p.nonad)
            rownames(out) <- rownames(K)
          } else if (method == "notadjusted") {
            out <- data.frame(pm.p = p.nonad)
            rownames(out) <- rownames(K)
          } else {
            out <- data.frame(pm.p.tukey = p.tukey,
                              pm.p = p.nonad_2)
            rownames(out) <- rownames(tukey$linfct)
          }
        } else {
          out <- data.frame(pm.p.tukey = p.tukey,
                            pm.p = p.nonad_2)
          rownames(out) <- rownames(tukey$linfct)
        }


      }
      else {
        #don't have k matrix

        if (!missing(method)) {
          if (method == "dunnett") {
            out <- data.frame(pm.p.dunnett = p.dunnett,
                              pm.p = p.nonad_1)
            rownames(out) <- rownames(dunnett$linfct)
          } else if (method == "tukey") {
            out <- data.frame(pm.p.tukey = p.tukey,
                              pm.p = p.nonad_2)
            rownames(out) <- rownames(tukey$linfct)
          } else {
            out <- "Missing K matrix"
          }
        } else {
          out <- data.frame(pm.p.tukey = p.tukey,
                            pm.p = p.nonad_2)
          rownames(out) <- rownames(tukey$linfct)
        }
      }
    }


    return(out)

  }
