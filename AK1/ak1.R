source('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/src/colors.R')
require(data.table)
require(ggplot2)
library(ggsignif)
library(extrafont)
# font_import()
# fonts()
# fonttable()
loadfonts(device="postscript")

auc <- fread('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/followups/ak1/nucleotides.txt', integer64 = 'double')
setnames(auc, c('Nucleotide', '-_1', '-_2', '-_3', '+_1', '+_2', '+_3'))
dt <- melt.data.table(auc, id.vars = 'Nucleotide', value.name =  'Abundance [AUC]', variable.name = 'AK1', variable.factor = F, value.factor = F, verbose = F)
dt[,Replicate:=as.integer(unlist(lapply(strsplit(AK1, '_'), function(i) tail(i, 1))))]
dt[,AK1:=factor(unlist(lapply(strsplit(AK1, '_'), function(i) head(i, 1))))]
dt[,Nucleotide:=factor(Nucleotide, levels = c('ATP', sort(unique(Nucleotide)[!unique(Nucleotide)%in%'ATP'])), ordered = T)]

dt <- dt[,list(`Abundance [AUC]`=mean(`Abundance [AUC]`, na.rm = T),
               SD=sd(`Abundance [AUC]`, na.rm = T)),by = list(Nucleotide, AK1)]

pal <- c(tumblue, tumred)
names(pal) <- c('-', '+')

my_theme <- function (base_size = 12) # , base_family="Arial"
{
  require(grid)
  theme(text = element_text(size = base_size, face
                            = "plain",
                            colour = "black", hjust = 0.5, vjust = 0.5,
                            angle = 0, lineheight = 1, family = 'TUM Neue Helvetica 55 Regular'), # , family = base_family
        axis.text  = element_text(size = rel(1), colour="black"),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
        axis.ticks = element_line(colour = "black", lineend = 'round'),
        axis.ticks.length = unit(2, units = 'pt'),
        axis.title.y = element_text(margin = margin(0,2,0,0, unit = 'pt'), size = rel(1), vjust=0.5,angle=90),
        axis.title.x = element_text(margin = margin(2,0,0,0,unit = 'pt'), size = rel(1), hjust=0.5, vjust=0.5),
        legend.key = element_blank(),
        legend.key.size = unit(0.7, units = 'lines'),
        legend.title = element_text(size = rel(1), face="plain"),
        legend.text = element_text(size = rel(1), face="plain"),
        legend.spacing = unit(0,"cm"),
        legend.background = element_blank(),
        legend.position = c(0.11,0.82),
        legend.direction = 'vertical',
        legend.box.background = element_blank(),
        legend.box.margin = margin(0,0,0,0, unit = 'pt'),
        legend.box.spacing = unit(0, units = 'pt'),
        legend.margin = margin(0,0,0,2),
        plot.title = element_text(size = rel(1), colour = "black",face="plain", vjust = 0.5, hjust = 0.5, margin = margin(0,0,4,0, unit = 'pt')),
        plot.margin = margin(1,14,5,1,unit = 'pt'),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", size = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey90", colour = "grey50", size=0.3)
  )
}

setEPS()
postscript('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/res/20180419_AK1/nucleotides.eps', width = 5.7 /2.539998, height = 5.7/2.539998, pointsize = 10, family = 'TUM Neue Helvetica 55 Regular')
ggplot(dt, aes(x = Nucleotide, y = `Abundance [AUC]`, fill = `AK1`)) +
  geom_bar(position="dodge", stat="identity", colour = 'black', size = 0.4) +
  geom_errorbar(aes(ymin=`Abundance [AUC]`-SD, ymax=`Abundance [AUC]`+SD), width=.2,
                position=position_dodge(.9)) +
  geom_signif(annotations = c(rep('***', 5)), y_position = c(6600000, 6800000, 8600000, 10100000, 8000000), xmin= 0.77+0:4, xmax=1.23+0:4) +
  my_theme(base_size = 10) +
  labs(title = 'Dephophorylation of\nantimetabolites by AK1') +
  scale_fill_manual(values = pal) +
  ylim(c(0, 1.5e7)) +
  guides(colour = guide_legend(ncol = 1, title.position = "top", title.hjust = 0, legend.position = 'right',keyheight = 1))
dev.off()
