multiSetMEsFIX <- function (exprData, colors, universalColors = NULL, useSets = NULL, 
          useGenes = NULL, impute = TRUE, nPC = 1, align = "along average", scale = T,
          excludeGrey = FALSE, grey = if (is.null(universalColors)) {
            if (is.numeric(colors)) 
              0
            else "grey"
          } else if (is.numeric(universalColors)) 0 else "grey", subHubs = TRUE, 
          trapErrors = FALSE, returnValidOnly = trapErrors, softPower = 6, 
          verbose = 1, indent = 0) 
{
  spaces = indentSpaces(indent)
  nSets = length(exprData)
  setsize = checkSets(exprData, useSets = useSets)
  nGenes = setsize$nGenes
  nSamples = setsize$nSamples
  if (verbose > 0) 
    printFlush(paste(spaces, "multiSetMEs: Calculating module MEs."))
  MEs = vector(mode = "list", length = nSets)
  consValidMEs = NULL
  if (!is.null(universalColors)) 
    consValidColors = universalColors
  if (is.null(useSets)) 
    useSets = c(1:nSets)
  if (is.null(useGenes)) {
    for (set in useSets) {
      if (verbose > 0) 
        printFlush(paste(spaces, "  Working on set", 
                         as.character(set), "..."))
      if (is.null(universalColors)) {
        setColors = colors[, set]
      }
      else {
        setColors = universalColors
      }
      setMEs = moduleEigengenes(expr = exprData[[set]]$data, 
                                colors = setColors, impute = impute, nPC = nPC, 
                                align = align, excludeGrey = excludeGrey, grey = grey, 
                                trapErrors = trapErrors, subHubs = subHubs, returnValidOnly = FALSE, 
                                softPower = softPower, verbose = verbose - 1, 
                                indent = indent + 1, scale = scale)
      if (!is.null(universalColors) && (!setMEs$allOK)) {
        if (is.null(consValidMEs)) {
          consValidMEs = setMEs$validMEs
        }
        else {
          consValidMEs = consValidMEs * setMEs$validMEs
        }
        consValidColors[setMEs$validColors != universalColors] = setMEs$validColors[setMEs$validColors != 
                                                                                      universalColors]
      }
      MEs[[set]] = setMEs
      names(MEs[[set]])[names(setMEs) == "eigengenes"] = "data"
    }
  }
  else {
    for (set in useSets) {
      if (verbose > 0) 
        printFlush(paste(spaces, "  Working on set", 
                         as.character(set), "..."))
      if (is.null(universalColors)) {
        setColors = colors[useGenes, set]
      }
      else {
        setColors = universalColors[useGenes]
      }
      setMEs = moduleEigengenes(expr = exprData[[set]]$data[, 
                                                            useGenes], colors = setColors, impute = impute, 
                                nPC = nPC, align = align, excludeGrey = excludeGrey, 
                                grey = grey, trapErrors = trapErrors, subHubs = subHubs, 
                                returnValidOnly = FALSE, softPower = softPower, 
                                verbose = verbose - 1, indent = indent + 1, scale = scale)
      if (!is.null(universalColors) && (!setMEs$allOK)) {
        if (is.null(consValidMEs)) {
          consValidMEs = setMEs$validMEs
        }
        else {
          consValidMEs = consValidMEs * setMEs$validMEs
        }
        consValidColors[setMEs$validColors != universalColors[useGenes]] = setMEs$validColors[setMEs$validColors != 
                                                                                                universalColors[useGenes]]
      }
      MEs[[set]] = setMEs
      names(MEs[[set]])[names(setMEs) == "eigengenes"] = "data"
    }
  }
  if (!is.null(universalColors)) {
    for (set in 1:nSets) {
      if (!is.null(consValidMEs)) 
        MEs[[set]]$validMEs = consValidMEs
      MEs[[set]]$validColors = consValidColors
    }
  }
  for (set in 1:nSets) {
    MEs[[set]]$allOK = (sum(!MEs[[set]]$validMEs) == 0)
    if (returnValidOnly) {
      valid = (MEs[[set]]$validMEs > 0)
      MEs[[set]]$data = MEs[[set]]$data[, valid]
      MEs[[set]]$averageExpr = MEs[[set]]$averageExpr[, 
                                                      valid]
      MEs[[set]]$varExplained = MEs[[set]]$varExplained[, 
                                                        valid]
      MEs[[set]]$isPC = MEs[[set]]$isPC[valid]
      MEs[[set]]$allPC = (sum(!MEs[[set]]$isPC) == 0)
      MEs[[set]]$isHub = MEs[[set]]$isHub[valid]
      MEs[[set]]$validAEs = MEs[[set]]$validAEs[valid]
      MEs[[set]]$allAEOK = (sum(!MEs[[set]]$validAEs) == 
                              0)
      MEs[[set]]$validMEs = rep(TRUE, times = ncol(MEs[[set]]$data))
    }
  }
  names(MEs) = names(exprData)
  MEs
}
environment(multiSetMEsFIX) <- asNamespace('WGCNA')