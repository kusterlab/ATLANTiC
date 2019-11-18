plotFACI_FIX <- function (effList, indAxis = c("ED", "EF"), caRef = TRUE, showPoints = FALSE, 
          add = FALSE, ylim, ...) 
{
  indAxis <- match.arg(indAxis)
  faValues <- effList[["EDvec"]]
  minfa <- min(faValues)
  faValues[faValues < 0] <- -(100 - abs(faValues[faValues < 
                                                   0]))
  plotMat <- switch(indAxis, ED = effList[["CAx"]], EF = effList[["CAy"]])
  xVec <- as.numeric(rownames(plotMat))
  xVec[faValues < 0] <- rev(xVec[faValues < 0])
  xVec[faValues > 0] <- rev(xVec[faValues > 0])
  yVec <- plotMat[, 1]
  seVec <- plotMat[, 2]
  if (caRef) {
    yLimits <- c(0, max(yVec))
  }
  else {
    yLimits <- range(yVec)
  }
  if (!missing(ylim)) {
    yLimits <- ylim
  }
  if (!add) {
    plot(xVec, yVec, type = "l", xlim = c(minfa - 1, max(faValues) + 
                                            1), ylim = yLimits, #xlab = "Fraction affected", ylab = "Combination index", 
         ...)
    # abline(h = 1, lty = 2, lwd = 2)
  }
  else {
    lines(xVec, yVec, ...)
  }
  dispersion(xVec, yVec, ulim = 1.96 * seVec, arrow.gap = 0.15, 
             ...)
  if (showPoints) {
    points(xVec, yVec, pch = 1)
  }
  invisible(plotMat)
}
environment(plotFACI_FIX) <- asNamespace('drc')