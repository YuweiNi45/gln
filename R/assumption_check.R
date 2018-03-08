#' Create tables and graphs providing cell means model assumption checking
#'
#' Based on cell means model, model assumption checking was performed by displaying tables containing test results and graphs. Normality and heteroscedasticity are examined in this function. Shapiro-Wilk normality test, qq-plot and histogram are provided to test normality. Bartlett test, residual plot and boxplot are for heteroscedasticity test.
#'
#' @param x   a vector that contains different levels whose variances are to be compared
#' @param y   a numeric vector of data values
#' @param data  a data frame, list or environment (or object coercible by as.data.frame to a data frame)
#' @param digits   the digits of the data values (default is 3)
#' @param table    a boolean value specifying providing table or not, must be TRUE (default) or FALSE
#' @param graph   a boolean value specifying providing graph or not, must be TRUE (default) or FALSE
#' @param ...	further arguments to be passed.
#'
#'
#' @return By default, it returns both tables and graphs. Tables contain results of Shapiro-Wilk normality test and Bartlett test. Graphs are qq-plot, boxplot, residual plot and histogram. If only asking for graphs, specify table = False. If only asking for tables, specify graph = False.
#'
#' @references Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988) The New S Language. Wadsworth & Brooks/Cole.
#'
#' @seealso \code{\link[stats]{shapiro.test}}
#'
#' \code{\link[stats]{bartlett.test}}
#'
#' @examples
#' x <- c (1,1,1,1,2,2,2,2,3,3,3,3)
#'
#'y <- c (1,3,3,2,4,3,2,3,2,3,3,2)
#'
#'data <- data.frame(cbind(x,y))
#'
#'assumption_check(x,y,data)
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
