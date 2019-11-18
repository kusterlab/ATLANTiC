
#Colors
TUMmauve <-"#69085A" 
TUMviolet <-"#0F1B5F" 
TUMblackblue <- "#003359"
TUMdarkblue <- "#005293"
TUMblue <- "#0073CF"
TUMlightblue <- "#64A0C8"
TUMciel <- "#98C6EA"
TUMturk <- "#00778A"
TUMgreen <- "#007C30"
TUMlightgreen <- "#679A1D"
TUMolive <- "#A2AD00"
TUMbeige <- "#DAD7CB"
TUMyellow <- "#FFDC00"
TUMgold <- "#F9BA00"
TUMorange <- "#E37222"
TUMbrick <- "#D64C13"
TUMred <- "#C4071B"
TUMblood <- "#9C0D16"

require(gplots)

colorp <- c(TUMmauve,
            TUMblue,
            TUMorange,
            TUMolive,
            TUMgold,
            TUMturk, 
            TUMbeige,
            col2hex("dodgerblue4"),
            TUMciel, 
            col2hex("mediumorchid"),
            col2hex("dimgrey"),
            col2hex("darkgreen"),
            TUMred, 
            TUMviolet, 
            TUMlightgreen, 
            TUMblood,
            TUMyellow, 
            TUMbrick,
            TUMdarkblue,
            TUMgreen
)

plcols <- c(TUMblue,
            TUMolive,
            TUMorange,
            TUMviolet,
            TUMgold,
            TUMblood,
            TUMturk, 
            TUMbeige,
            col2hex("dodgerblue4"),
            TUMciel, 
            col2hex("mediumorchid"),
            col2hex("dimgrey"),
            col2hex("darkgreen"),
            TUMred, 
            TUMmauve, 
            TUMlightgreen, 
            TUMyellow, 
            TUMbrick,
            TUMdarkblue,
            TUMgreen
)

plcols2 <- c(TUMblue,
             TUMolive,
             TUMorange,
             TUMviolet,
             TUMgold,
             TUMblood,
             TUMturk, 
             TUMbeige,
             col2hex("dodgerblue4"),
             TUMciel, 
             col2hex("mediumorchid"),
             col2hex("dimgrey"),
             col2hex("darkgreen"),
             TUMred, 
             TUMmauve, 
             TUMlightgreen, 
             TUMyellow, 
             TUMbrick,
             TUMdarkblue,
             TUMgreen
)

#Define ggplot2 theme
my_theme <- function (base_size = 10) # , base_family="Arial"
{
  require(grid)
  theme(text = element_text(size = base_size, face
                            = "plain",
                            colour = "black", hjust = 0.5, vjust = 0.5,
                            angle = 0, lineheight = 0.9), # , family = base_family
        axis.text  = element_text(size = rel(0.8), colour="black"),
        axis.ticks = element_line(colour = "black"),
        axis.title.y = element_text(size = rel(1), vjust=0.5,angle=90),
        axis.title.x = element_text(size = rel(1), hjust=0.5, vjust=0.5),
        legend.key = element_rect(colour = "white", fill= "white"),
        legend.title = element_text(size = rel(1), face="plain"),
        legend.text = element_text(size = rel(0.8), face="plain"),
        legend.margin = unit(0,"cm"),
        legend.background = element_rect(fill= NA),
        plot.title = element_text(size = rel(1.2), colour = "black",face="plain", vjust = 1.5),
        panel.background = element_rect(fill = "white",colour = NA),
        panel.border = element_rect(fill = NA, colour = "black"),
        panel.grid.major = element_line(colour = "grey70", size = 0.33,linetype="dotted"),
        panel.grid.minor = element_line(colour = "grey70", size = 0.33,linetype="dotted"),
        strip.background = element_rect(fill = "grey90", colour = "grey50", size=0.3)
  )
}
