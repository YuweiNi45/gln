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

all_sum_test <- function(x, y, data, alternative, graph, digits = 3, ...){

  data_sum <- function(y){
    out <- data.frame(n = length(y),
                      n.complete = length(y[!is.na(y)]),
                      mean = mean(y, na.rm=T),
                      sd = sd(y, na.rm=T),
                      median = median(y, na.rm=T),
                      lower_quartile = quantile(y, 0.25, na.rm=T),
                      upper_quartile = quantile(y, 0.75, na.rm=T))
    return(out)
  }

  smry <- tapply(y, x, data_sum)
  smry <- do.call(rbind, smry)

  lm <- lm (y ~ x, data = data, ...)

  anova <- anova(lm (y ~ 1, data = data, ...),lm)

  kruskal <- kruskal.test(y~x)

  test<-function(x){
    if(x < 0.001) x.txt<- "<0.001"
    else if(x<0.05) x.txt<-"<0.05"
    else x.txt<- format(x,digits = 2)
    x.txt
  }

  numrow <- nrow(smry)-1

  p_anova <- test(anova$`Pr(>F)`[2])

  p_kw <- test(kruskal$p.value)

  out <- data.frame(smry[,1:2],
                    mean.sd = paste0(format(smry$mean, digits=digits),
                                     " +/- ",
                                     format(smry$sd, digits=digits)),
                    median.IQR = paste0(format(smry$median, digits=digits),
                                        " (",
                                        format(smry$lower_quartile, digits=digits), ",", format(smry$upper_quartile, digits = digits), ")"),
                    p.anova=paste0(c(rep("", numrow), p_anova)),
                    p.kw=paste0(c(rep("", numrow), p_kw))
  )
  names(out)[3] <- "mean +/- sd"
  names(out)[4] <- "median(IQR)"

  graph_test <- function(x){
    if(x < 0.001) x.txt<- "<0.001"
    else if(x<0.05) x.txt<-"<0.05"
    else x.txt<- paste("=", format(x,digits = 2))
    x.txt
  }

  p_anova_graph <- graph_test(anova$`Pr(>F)`[2])

  p_kw_graph <- graph_test(kruskal$p.value)

  ##making the plot
  x_cat <- levels(x)

  x_num <- length(x_cat)

  x_mean <- numeric(x_num)

  x_se <- numeric(x_num)

  for (i in 1:x_num) {
    x_mean[i] <- mean(y[which(x == x_cat[i])])
  }

  for (i in 1:x_num) {
    x_se[i] <- sd(y[which(x == x_cat[i])])/sqrt(length(y[which(x == x_cat[i])]))
  }

  data1 <- data.frame(x_mean, x_cat, x_se)

  dodge <- position_dodge(width = 0.9)
  limits <- aes(ymax = data1$x_mean + data1$x_se,
                ymin = data1$x_mean - data1$x_se)

  ##preparing for the y lab name

  y_length <- length(y)

  data_names <- colnames(data)

  name_length <- length(data_names)

  result <- numeric(name_length)

  for (i in 1:name_length) {
    result[i] <- sum(y %in% factor(unlist(data[data_names[i]])))
  }

  y_graph <- data_names[which(result == y_length)]

  ##plotting
  graph1 <-
    ggplot(data1, aes(x_cat, x_mean, fill = x_cat)) +
    geom_bar(stat="identity") +
    geom_errorbar(limits, position = dodge, width = 0.25) +
    annotate("text", label = paste("P", p_anova_graph, sep = ""), x = 1, y = ceiling(min(x_mean + x_se) + max(x_mean + x_se)), size = 3) +
    theme(panel.background = element_rect(fill = "white", colour = "black"), axis.title.x=element_blank()) +
    labs(y = y_graph, fontface="bold") +
    theme(legend.position="none")

  graph2 <-
    ggplot(data, aes(x, y, fill = x)) +
    geom_boxplot() +
    annotate("text", label = paste("P", p_kw_graph, sep = ""), x = 1, y = ceiling(min(x_mean + x_se) + max(x_mean + x_se)), size = 3) +
    theme(panel.background = element_rect(fill = "white", colour = "black"), axis.title.x=element_blank()) +
    labs(y = y_graph, fontface="bold") +
    theme(legend.position="none")

  if (!missing(alternative)) {
    if (alternative == "overall"){
      if (!missing(graph)){
        if (graph == T){
          lst <- list(out, graph1, graph2)
          names(lst) <- c("Numerical Analysis", "Barplot", "Boxplot")
          return(lst)
        }
        else {
          return(out)
        }
      }else {
        lst <- list(out, graph1, graph2)
        names(lst) <- c("Numerical Analysis", "Barplot", "Boxplot")
        return(lst)
      }
    }
    else if (alternative == "parametric"){
      if (!missing(graph)){
        if (graph == F){
          return(out[,c(1:3, 5)])
        }
        else {
          lst <- list(out[,c(1:3, 5)], graph1)
          names(lst) <- c("Numerical Analysis", "Barplot")
          return(lst)
        }
      } else {
        lst <- list(out[,c(1:3, 5)], graph1)
        names(lst) <- c("Numerical Analysis", "Barplot")
        return(lst)
      }
    }
    else if (alternative == "nonparametric") {
      if (!missing(graph)){
        if (graph == F) {
          return(out[,c(1:2, 4, 6)])
        }
        else{
          lst <- list(out[,c(1:2, 4, 6)], graph2)
          names(lst) <- c("Numerical Analysis", "Boxplot")
          return(lst)
        }
      } else {
        lst <- list(out[,c(1:2, 4, 6)], graph2)
        names(lst) <- c("Numerical Analysis", "Boxplot")
        return(lst)
      }
    }
  }
  else {
    if (!missing(graph)){
      if (graph == F){
        return(out)
      }
      else {
        lst <- list(out, graph1, graph2)
        names(lst) <- c("Numerical Analysis", "Barplot", "Boxplot")
        return(lst)
      }
    }else {
      lst <- list(out, graph1, graph2)
      names(lst) <- c("Numerical Analysis", "Barplot", "Boxplot")
      return(lst)
    }
  }
}
