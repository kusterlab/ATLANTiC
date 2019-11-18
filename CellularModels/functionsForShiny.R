require("reshape2")

#myPlot is part of that function, this function adjust the information in order to plot the plot
adjustments2 <- function(total_result,drug, expression, total_normalized_data, desired_drug, desired_cellline, desired_protein, normalizing_method){
  total <- drug[,c("Gene.names", desired_drug)]
  total <- total[!is.na(total[,desired_drug]),]
  
  total[,2] <- -log10(total[,2]*10e-9)
  
  
  total <- merge(total, expression[, c("Gene.names", desired_cellline)], by = "Gene.names")
  total["scaleup"] <- NA
  total["scaleup_2"] <- NA


  dgg <- dcast(total_normalized_data[total_normalized_data$Cellline.names == desired_cellline], Gene.names ~ Dataset, value.var = "Expression")
  
  total <- cbind(total, dgg[dgg$Gene.names %in% total$Gene.names,-1])
  total <- total[order(total[,desired_drug], decreasing = T),]
  rownames(total) <- total$Gene.names
  
  for(i in c(6:ncol(total))){
    colnames(total)[6] <- paste0("normalized_",i-5)
  }
  
  y1<- ifelse(-max(total_normalized_data$Expression) < min(total_normalized_data$Expression), 1.2*-max(total_normalized_data$Expression), 1.2*min(total_normalized_data$Expression))
  y2 <- 1.2*max(total_normalized_data$Expression)
  
  myPlot2(total_result,expression, total, desired_cellline, desired_drug, desired_protein, y_lim=c(y1,y2), normalizing_method = normalizing_method)
}


#plots the plots when selecting the selectivity score
myPlot2 <- function(total_result,expression, total, desired_cellline, desired_drug, desired_protein, y_lim, normalizing_method){
  
  #create empty plot
  par(mar = c(5,4,4,5), xpd = F)
  
  plot <- barplot(height = as.matrix(t(total[,grep("normalized", colnames(total))])), 
                  names.arg = total$Gene.names,
                  las = 2 ,
                  cex.names = ifelse(nrow(total) < 50, 1, 0.7),
                  ylim = y_lim,
                  cex.axis = 0.7,
                  beside = T)
  
  #creating rectangles
  pKd <- min(total[total$Gene.names %in% desired_protein,desired_drug])
  pKdx10 <- pKd - 1
  pKdx100 <- pKd - 2
  
  #darkgrey_box <- total[total$AK_1pKd >= pKdx10,]$AK_1 
  darkgrey_box <- sum(total[, desired_drug] >= pKdx10)
  lightgrey_box <- sum(total[, desired_drug] < pKdx10 & total[, desired_drug] >= pKdx100) 
  
  
  rect(xleft=0, ybottom = y_lim[1], xright = plot[nrow(plot), darkgrey_box]+1, ytop = y_lim[2]  ,col = "darkgrey")
  mtext1 <- plot[nrow(plot), darkgrey_box]/2
  if(darkgrey_box != 0){
    mtext("x10", side = 3, line = -2, at = mtext1, cex = 0.8)
  }
  
  if(lightgrey_box != 0){
    mtext2<- plot[nrow(plot), darkgrey_box]+plot[nrow(plot), lightgrey_box]/2
    rect(xleft = (plot[nrow(plot), darkgrey_box]) + 1 , ybottom = y_lim[1], xright = (plot[nrow(plot), darkgrey_box] + plot[nrow(plot), lightgrey_box])+1.5, ytop = y_lim[2], col = "lightgrey") 
    mtext("x100", side = 3, line = -2, at = mtext2, cex = 0.8)
  }
  
  #plotting
  plot <- barplot(height = as.matrix(t(total[,grep("normalized", colnames(total))])), 
                  names.arg = total$Gene.names,
                  las = 2 ,
                  cex.names = ifelse(nrow(total) < 50, 1, 0.7),
                  ylim = y_lim,
                  ylab = paste0( normalizing_method, " of ",desired_cellline),
                  cex.axis = 0.7,
                  add = T, 
                  main = paste("Interactions in ", desired_cellline, "with", desired_drug), beside = T,
                  col="royalblue")
  
  #create new y-axis
  axis(side = 4, las = 2, at = c( 2*y_lim[1]*1.1/3,y_lim[1]*1.1/3, y_lim[1]*1.1-y_lim[1]), 
       labels = (c(round(max(total[,desired_drug])/3, digits = 1),
                   round(2*max(total[,desired_drug])/3, digits = 1),
                   round(max(total[,desired_drug]), digits = 1) 
       )),
       cex.axis = 0.8)
  
  xmargin <- diff(grconvertX(0:1, "inches","user"))*par("cin")[2]*par("cex")*par("lheight")
  text(par("usr")[2]+2.5*xmargin,mean(par("usr")[3:4]), expression("        "~pK[d]^app), srt = -90, xpd = TRUE, pos = 4, adj = 0.25)
  
  #create scale-up
  var1 <- abs(y_lim[1]/2)*1.1
  var2 <- (max(total[,desired_drug])/2)*1.1
  scaleup <-  (var1/var2)
  total$scaleup <- total[,desired_drug]*scaleup + (y_lim[1]*1.1)
  
  #create points
  #ii <- sapply(1:nrow(total), function(i){ifelse(i ==1, 1.5, i*ifelse(is.null(ncol(total[,grep("normalized", colnames(total))])),2,ncol(total[,grep("normalized", colnames(total))])+1)-0.5)})
  points(plot[1,], y = total$scaleup,
         pch = 19, cex = 0.6)
  if(is.null(ncol(total[,grep("normalized", colnames(total))]))){
    
    text(par("usr")[2]+2.75*xmargin,mean(par("usr")[3:4]), expression("log2(raw intensities)       "), srt = -90, xpd = TRUE, pos = 2, adj = 0.25)
    
    axis(side = 4, las = 2, at = c(0, (y_lim[2]/2), y_lim[2]), 
         labels = c("0", round(max(expression[,-1])/2, digits = 0), round(max(expression[,-1]), digits = 0)),
         cex.axis = 0.7)
    
    #absolut values of iBAQ
    var3 <- y_lim[2]/2
    var4 <- max(expression[,-1])/2
    scaleup_2 <-  (var3/var4)
    total$scaleup_2 <- total[,desired_cellline]*scaleup_2
    points(plot[1,], y = total$scaleup_2,
           pch = 24, col = "black", bg = "red", cex = 0.6)
    
  }
  
  #add line depending on protein
  for(i in c(1:length(desired_protein))){
    vertical_1 <- which(total$Gene.names == desired_protein[i])*ifelse(which(total$Gene.names == desired_protein[1])==1,2, ifelse(is.null(ncol(total[,grep("normalized", colnames(total))])), 2,ncol(total[,grep("normalized", colnames(total))])+1))-0.5
    abline(v = vertical_1, lty = 3)
  }
  
  #text with targetselectivity
  # mtext(round(total_result[total_result$Gene.Names == desired_protein & total_result$Drugs == desired_drug & total_result$Cellline.names == desired_cellline, "Selectivity_Score"], digits = 4), side = 3)
  mtext(signif(total_result[total_result$Gene.Names == desired_protein & total_result$Drugs == desired_drug & total_result$Cellline.names == desired_cellline, "Selectivity_Score"], digits = 3), side = 3)
  
  
  abline(h = 0, v = 0, col = "black")  
}

