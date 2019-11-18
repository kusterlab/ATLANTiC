library(DT)
library(shiny)

ui <- fluidPage(
  fluidRow(
    column(width = 6, DT::dataTableOutput('tab')),
    column(width = 6, uiOutput("plots"))
  )
)
