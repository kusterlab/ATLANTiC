library(DT)
library(shiny)
source("functionsForShiny.R")

server <- function(input, output) {
  
  dat <- readRDS("dat.RDS")
  
  output$tab <- DT::renderDataTable(
    dat$selectivityScore,
    rownames = FALSE,
    colnames = c("Data panel", "Drug", "Gene name", "Cell line", "Selectivity score"),
    filter = 'top',
    selection = "multiple",
    options = list(scrollX = TRUE, pageLength = 20)
  )
  
  output$plots <- renderUI({
    req(input$tab_rows_selected)
    np <- min(3, length(input$tab_rows_selected))
    comp <- rev(rev(input$tab_rows_selected)[1:np])
    
    lapply(comp, function(i) {
      isolate({
        output[[paste0("plot", i)]] <- renderPlot({ 
          dataset <- as.character(dat$selectivityScore$Dataset[i])
          adjustments2(total_result = dat$selectivityScore[dat$selectivityScore$Dataset %in% dataset, ], 
                       drug = dat$cpd, 
                       expression = dat[[dataset]] , 
                       total_normalized_data = dat$zscore[[dataset]], 
                       desired_drug =  as.character(dat$selectivityScore$Drugs[i]), 
                       desired_cellline = as.character(dat$selectivityScore$Cellline.names[i]), 
                       desired_protein = as.character(dat$selectivityScore$Gene.Names[i]), 
                       normalizing_method = "z-score")
        }, height = 300)
      })
    })
  })
}
