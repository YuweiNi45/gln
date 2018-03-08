#' create tables and graphs providing overall multiple comparison
#'
#' Generic function to create tables and plots providing basic overall numerical and graphical summary in different groups for comparison, based on parametric and/or non-parametric approach.
#'
#' @param x   a vector that contains different levels whose summary results are wanted for comparison.
#' @param y   a numeric vector of data values
#' @param data  a data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables to be summarized.
#' @param alternative	  a character string specifying the alternative approaches, must be one of "overall" (default), "parametric" or " non-parametric".
#' @param graph	  logical. If TRUE the corresponding plots (bar plot, box plot) are returned.
#' @param digits     the digits of the data values (default is 3)
#' @param ...	  further arguments to be passed to or from methods.
#'
#' @return  For parametric approach, it returns a table containing sample size, number of observations without missing value, mean, standard deviation for each level and F-test results, and a bar plot showing mean, standard error for each level and F-test results. For non-parametric approach, it returns a table containing sample size, number of observations without missing value, median, IQR for each level and Kruskal-Wallis test results, and a boxplot showing median, IQR for each level and Kruskal-Wallis test results. For overall, it returns a table containing both the summary results in parametric and non-parametric approach, a bar plot showing mean, standard error for each level and F-test results, and a boxplot showing median, IQR for each level and Kruskal-Wallis test results.
#'
#' @references Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988) The New S Language. Wadsworth & Brooks/Cole.
#'
#' @seealso \code{\link[stats]{anova}}
#'
#' \code{\link[stats]{kruskal.test}}
#'
#' @examples
#' x <- c(1,1,2,2,3,3,4,4)
#'
#'y <- c(1,3,2,3,3,3,2,4)
#'
#'data <- as.data.frame(cbind(x,y))
#'
#'all_sum_test("x", "y", data, alternative = "parametric")
#'
#'all_sum_test(x,y,data)
#'
#'@import ggplot2
#'@import multcomp
#'
#'@export

assumption_check <-
  function(x,
           y,
           data,
           digits = 3,
           table = c(TRUE, FALSE),
           graph = c(TRUE, FALSE),
           ...) {
    out <- lm(y ~ x - 1, data = data, ...)

    qqplot <-
      ggplot(data, aes(sample = rstudent(out))) +
      stat_qq(shape = 1, size = 3) +
      geom_abline(slope = 1, intercept = 0) +
      theme(panel.background = element_rect(fill = "white", colour = "black")) +
      ggtitle("Normal Q-Q Plot") +
      xlab("Theoreticla Quantiles") +
      ylab("Sample Quantiles")

    f <- length(x)

    data1 <- data.frame(1:f, rstudent(out))

    variance <- ggplot(data1, aes(1:f, rstudent(out))) +
      geom_point(shape = 3, size = 3) +
      theme(panel.background = element_rect(fill = "white", colour = "black")) +
      xlab("Index") +
      ylab("Studentized Residual")

    shap.text <- shapiro.test(y)
    bar.test <- bartlett.test(y, x)

    ##preparing for the y lab

    y_length <- length(y)

    data_names <- colnames(data)

    name_length <- length(data_names)

    result <- numeric(name_length)

    for (i in 1:name_length) {
      result[i] <- sum(y %in% factor(unlist(data[data_names[i]])))
    }

    y_graph <- data_names[which(result == y_length)]

    graph2 <-
      ggplot(data, aes(x, y, fill = x)) +
      geom_boxplot() +
      labs(y = y_graph, fontface = "bold") +
      theme(
        panel.background = element_rect(fill = "white", colour = "black"),
        axis.title.x = element_blank()
      ) +
      theme(legend.position = "none")

    hist <-
      ggplot(data, aes(y)) +
      geom_histogram(binwidth = (max(y) - min(y)) / length(y)) +
      theme(
        panel.background = element_rect(fill = "white", colour = "black"),
        axis.title.x = element_blank()
      ) +
      labs(title = paste("Histogram of ", y_graph, sep = ""),
           fontface = "bold")

    table1 <-
      data.frame(
        p.shap.text = format(shap.text$p.value, digits = digits),
        p.bar.test = format(bar.test$p.value, digits = digits)
      )

    if (!missing(table) && !missing(graph)) {
      if (table == TRUE && graph == TRUE) {
        lst <- list(table1, qqplot, variance, graph2, hist)
        names(lst) <-
          c(
            "Shapiro-Wilk Normality Test & Bartlett test",
            "QQ-plot",
            "Assessing equal variance assumption",
            "Boxplot",
            "Histogram"
          )
      } else if (table == TRUE && graph == FALSE) {
        lst <- list(table1)
        names(lst) <-
          c("Shapiro-Wilk Normality Test & Bartlett test")
      } else if (table == FALSE && graph == TRUE) {
        lst <- list(qqplot, variance, graph2, hist)
        names(lst) <-
          c("QQ-plot",
            "Assessing equal variance assumption",
            "Boxplot",
            "Histogram")
      } else {
        lst <- "Null"
      }
    } else if (!missing(table) && missing(graph)) {
      if (table == TRUE) {
        lst <- list(table1, qqplot, variance, graph2, hist)
        names(lst) <-
          c(
            "Shapiro-Wilk Normality Test & Bartlett test",
            "QQ-plot",
            "Assessing equal variance assumption",
            "Boxplot",
            "Histogram"
          )
      } else {
        lst <- list(qqplot, variance, graph2, hist)
        names(lst) <-
          c("QQ-plot",
            "Assessing equal variance assumption",
            "Boxplot",
            "Histogram")
      }
    } else if (missing(table) && !missing(graph)) {
      if (graph == TRUE) {
        lst <- list(table1, qqplot, variance, graph2, hist)
        names(lst) <-
          c(
            "Shapiro-Wilk Normality Test & Bartlett test",
            "QQ-plot",
            "Assessing equal variance assumption",
            "Boxplot",
            "Histogram"
          )
      } else {
        lst <- list(table1)
        names(lst) <-
          c("Shapiro-Wilk Normality Test & Bartlett test")
      }
    } else {
      lst <- list(table1, qqplot, variance, graph2, hist)
      names(lst) <-
        c(
          "Shapiro-Wilk Normality Test & Bartlett test",
          "QQ-plot",
          "Assessing equal variance assumption",
          "Boxplot",
          "Histogram"
        )
    }

    return(lst)

  }
