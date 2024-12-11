library(ggpubr)
library(ggplot2)

#' Plotting Univariate Mark-Specific Proportional Hazards Model Fits Using ggplot
#'
#' \code{ggplot} method for class \code{summary.sievePH} for univariate marks. Point and interval estimates of the mark-specific treatment effect parameter specified by component  \code{contrast} in \code{\link{summary.sievePH}} are plotted, together with scatter/box plots of the observed mark values by treatment.
#' @param x an object returned by \code{\link{summary.sievePH}}
#' @param mark a numeric vector specifying a univariate continuous mark. For subjects with a right-censored time-to-event, the value(s) in \code{mark} should be set to \code{NA}.
#' @param tx a numeric vector indicating the treatment group (1 if treatment, 0 if placebo)
#' @param xlim a numeric vector of length 2 specifying the x-axis range (\code{NULL} by default)
#' @param ylim a numeric vector of length 2 specifying the y-axis range (\code{NULL} by default)
#' @param xtickAt a numeric vector specifing the position of x-axis tickmarks (\code{NULL} by default)
#' @param xtickLab a numeric vector specifying labels for tickmarks listed in \code{xtickAt}. If \code{NULL} (default), the labels are determined by \code{xtickAt}.
#' @param ytickAt a numeric vector specifing the position of y-axis tickmarks (\code{NULL} by default)
#' @param ytickLab a numeric vector specifying labels for tickmarks listed in \code{ytickAt}. If \code{NULL} (default), the labels are determined by \code{ytickAt}.
#' @param tickLabSize a numeric value specifying the font size of tickmark labels along both axes in the bottom panel (\code{14} by default)
#' @param xlab a character string specifying the x-axis label (\code{NULL} by default)
#' @param ylab a character string specifying the y-axis label (\code{NULL} by default)
#' @param axisLabSize a numeric value specifying the font size of both axis labels in the bottom panel (\code{15} by default)
#' @param title a character string specifying the plot title (\code{NULL} by default)
#' @param titleSize a numeric value specifying the font size of the plot title (\code{16} by default)
#' @param subtitle a character string specifying the plot subtitle (\code{NULL} by default)
#' @param subtitleSize a numeric value specifying the font ize of the plot subtitle (\code{10} by default)
#' @param txLab a character vector of length 2 specifying the placebo and treatment labels (in this order). The default labels are \code{placebo} and \code{treatment}.
#' @param txLabSize a numeric value specifying the font size of labels \code{txLab} (\code{5} by default)
#' @param legendLabSize a numeric value specifying the font size of legend labels in the bottom panel (\code{11} by default)
#' @param estLineSize a numeric value specifying the line width for the point estimate of the mark-specific treatment effect (\code{1.6} by default)
#' @param ciLineSize a numeric value specifying the line width for the confidence limits for the mark-specific treatment effect (\code{1.2} by default)
#' @param boxplotWidth a numeric value specifying the width of each box in the box plot (\code{0.8}) by default
#' @param jitterFactor a numeric value specifying the amount of vertical jitter (\code{0.1} by default)
#' @param jitterSeed a numeric value setting the seed of R's random number generator for jitter in the scatter plot (\code{0} by default)
#' @param pointSize a numeric value specifying the size of data points in the scatter plot (\code{1.7} by default)
#' @param bottomPlotMargin a numeric vector, using cm as the unit, passed on to argument \code{plot.margin} in \code{theme()} for the bottom panel (\code{c(-0.5, 0.3, 0, 0)} by default)
#' @param topPlotMargin a numeric vector, using \code{"lines"} as the unit, passed on to argument \code{plot.margin} in \code{theme()} for the top panel (\code{c(0, 0.3, -0.12, 1.83)} by default)
#' @param plotHeights a numeric vector specifying relative heights of the top and bottom panels (\code{c(0.33, 0.67)} by default) passed on to argument \code{heights} in \code{ggpubr::ggarrange()}
#' @param ... other arguments to be passed to plotting functions
#'
#' @return A ggplot object.
#'
#' @examples
#' n <- 500
#' tx <- rep(0:1, each=n/2)
#' tm <- c(rexp(n/2, 0.2), rexp(n/2, 0.2 * exp(-0.4)))
#' cens <- runif(n, 0, 15)
#' eventTime <- pmin(tm, cens, 3)
#' eventInd <- as.numeric(tm <= pmin(cens, 3))
#' mark <- ifelse(eventInd==1, c(rbeta(n/2, 2, 5), rbeta(n/2, 2, 2)), NA)
#' markRng <- range(mark, na.rm=TRUE)
#'
#' # fit a model with a univariate mark
#' fit <- sievePH(eventTime, eventInd, mark, tx)
#' sfit <- summary(fit, markGrid=seq(markRng[1], markRng[2], length.out=10))
#' print(ggplot(sfit, mark, tx))
#'
#' @seealso \code{\link{plot.summary.sievePH}}, \code{\link{sievePH}} and \code{\link{summary.sievePH}}
#'
#' @export
ggplot.summary.sievePH <- function(x, mark=NULL, tx=NULL, xlim=NULL, ylim=NULL, xtickAt=NULL, xtickLab=NULL, ytickAt=NULL,
                                   ytickLab=NULL, tickLabSize = 14, xlab=NULL, ylab=NULL, axisLabSize = 15,
                                   title=NULL, titleSize = 16, subtitle=NULL, subtitleSize=10,
                                   txLab=c("Placebo", "Treatment"), txLabSize = 5, legendLabSize = 12,
                                   estLineSize=1.6, ciLineSize=1.2, boxplotWidth = 0.8,
                                   jitterFactor = 0.1, jitterSeed = 0, pointSize=1.7,
                                   bottomPlotMargin=c(-0.5, 0.5, 0.1, 0), topPlotMargin=c(1, 0.5, -0.12, 1.83),
                                   plotHeights=c(0.33, 0.67), ...){

  contrast <- names(x)[length(names(x))]

  # a 2-dimensional plot only when the mark is univariate
  if (NCOL(x[[contrast]])==4){
    if (is.null(xlab)){ xlab <- colnames(x[[contrast]])[1] }
    if (is.null(ylab)){ ylab <- switch(colnames(x[[contrast]])[2], TE="Treatment Efficacy", HR="Hazard Ratio", LogHR="Log Hazard Ratio") }

    fit <- x[[contrast]]
    colnames(fit) <- c("mark", "EST", "LB", "UB")
    if (is.null(xlim)){ xlim <- range(fit[, 1]) }
    if (is.null(ylim)){ ylim <- range(fit[, -1], na.rm=TRUE) }

    p1 <- ggplot(fit) +
      geom_hline(yintercept=ifelse(colnames(x[[contrast]])[2]=="HR", 1, 0), color="gray70", size=1) +
      geom_line(aes(x = mark, y = EST), linetype = "solid", size = estLineSize, na.rm = TRUE) +
      geom_line(aes(x = mark, y = LB, linetype = '95% Pointwise CI', size = '95% Pointwise CI'), na.rm = TRUE) +
      geom_line(aes(x = mark, y = UB, linetype = '95% Pointwise CI', size = '95% Pointwise CI'), na.rm = TRUE) +
      scale_linetype_manual(name="", labels = c('95% Pointwise CI'), values = c('95% Pointwise CI'= "dashed")) +
      scale_size_manual(name="", labels = c('95% Pointwise CI'), values = c('95% Pointwise CI'= ciLineSize)) +
      xlab(xlab) +
      ylab(ylab) +
      theme_bw() +
      theme(legend.key.size = unit(0.65, "cm"),
            legend.margin=margin(grid::unit(0,"cm")),
            legend.text=element_text(size=legendLabSize),
            legend.position = c(0.96, 1.08),
            legend.justification = c(1, 1),
            legend.key = element_blank(),
            legend.key.width = unit(1.4,"cm"),
            legend.background = element_blank(),
            axis.title.x = element_text(size = axisLabSize, margin = margin(t = 10)),
            axis.title.y = element_text(size = axisLabSize, margin = margin(r = 10)),
            axis.text.x = element_text(size = tickLabSize, colour = "black"),
            axis.text.y = element_text(size = tickLabSize, colour = "black"),
            plot.margin=unit(bottomPlotMargin, "cm"))

    # Add x-axis limits, ticks and labels
    if (is.null(xtickAt) && is.null(xtickLab)){
      p1 <- p1 + scale_x_continuous(limits = xlim, breaks = scales::pretty_breaks(n = 5))

    } else if (!is.null(xtickAt) && is.null(xtickLab)){
      p1 <- p1 + scale_x_continuous(limits = xlim, breaks = xtickAt, labels = xtickAt)

    } else if (is.null(xtickAt) && !is.null(xtickLab)){
      p1 <- p1 + scale_x_continuous(limits = xlim, breaks = xtickLab, labels = xtickLab)

    } else {
      p1 <- p1 + scale_x_continuous(limits = xlim, breaks = xtickAt, labels = xtickLab)
    }

    # Add y-axis limits, ticks and labels
    if (is.null(ytickAt) && is.null(ytickLab)){
      p1 <- p1 + scale_y_continuous(breaks = scales::pretty_breaks(n = 5))

    } else if (!is.null(ytickAt) && is.null(ytickLab)){
      p1 <- p1 + scale_y_continuous(breaks = ytickAt, labels = ytickAt)

    } else if (is.null(ytickAt) && !is.null(ytickLab)){
      p1 <- p1 + scale_y_continuous( breaks = ytickLab, labels = ytickLab)

    } else {
      p1 <- p1 + scale_y_continuous(breaks = ytickAt, labels = ytickLab)
    }
    p1 <- p1 + coord_cartesian(ylim = ylim)
    data <- data.frame("tx" = as.factor(tx), "mark" = mark)
    set.seed(jitterSeed)
    p2 <- ggplot(data) +
      geom_boxplot(aes(x = mark, y = tx), width = boxplotWidth, fill = "gray80", color = "black", lwd = 0.8, outlier.shape = NA, na.rm = TRUE, data = data)
    if (is.null(jitterFactor)){
      p2 <- p2 + geom_jitter(aes(x = mark, y = tx, shape = as.factor(tx), color = as.factor(tx)), alpha = 1, size = pointSize,
                             width = 0, fill = "white", stroke = 1, na.rm = TRUE, data = data)
    } else {
      p2 <- p2 + geom_jitter(aes(x = mark, y = tx, shape = as.factor(tx), color = as.factor(tx)), alpha = 1,size = pointSize, height = jitterFactor,
                             width = 0, fill = "white", stroke = 1, na.rm = TRUE, data = data)
    }
    p2 <- p2 + scale_shape_manual(name = "", labels = c("0" = txLab[1], "1" = txLab[2]), values = c(21, 24)) +
      scale_color_manual(name="", values = c("red3", "blue"), breaks = c("1", "0"))
    p2 <- p2 + geom_text(x = min(xlim) - 0.07 * (xlim[2] - xlim[1]), y = "0", label = txLab[1], size = txLabSize, hjust = 1, check_overlap = TRUE)
    p2 <- p2 + geom_text(x = min(xlim) - 0.07 * (xlim[2] - xlim[1]), y = "1", label = txLab[2], size = txLabSize, hjust = 1, check_overlap = TRUE)

    # Add the same x-axis limits and ticks as for 'p1'
    if (is.null(xtickAt)){
      p2 <- p2 + scale_x_continuous(limits = xlim, breaks = scales::pretty_breaks(n = 5), labels = NULL)
    } else {
      p2 <- p2 + scale_x_continuous(limits = xlim, breaks = xtickAt, labels = NULL)
    }

    if (!is.null(subtitle)){ p2 <- p2 + labs(subtitle = subtitle) }
    if (!is.null(title)){ p2 <- p2 + labs(title = title) }

    p2 <- p2 +
      xlab("") +
      ylab("") +
      scale_y_discrete(labels = c("1" = txLab[2], "0" = txLab[1]), breaks = c("1", "0")) +
      theme_bw() +
      guides(shape = "none") +
      guides(color = "none") +
      theme(plot.title = element_text(size = titleSize, hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(size = subtitleSize, hjust = 0.5),
            axis.text.y = element_blank(),
            axis.ticks.x = element_blank(),
            plot.margin=unit(topPlotMargin, "lines")) +
      coord_cartesian(clip = "off")

    p <- ggpubr::ggarrange(p2, p1, heights = plotHeights, ncol=1, nrow=2, align = "v")+theme(plot.margin = margin(0.1,0.1,0.1,0.1, "cm"))
    # p <- p2 + p1 + patchwork::plot_layout(ncol=1, height=plotHeights)

  } else {
    stop("Plotting of results is available for univariate marks only using this function.\n
         Consider plot.summary.sievePH() for bivariat marks.")
  }

  return(p)
}

