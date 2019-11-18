source('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/src/colors.R')
al2hex <- function(alphacol, bg = 'white') {
  compcol <- col2rgb(alphacol, alpha = T)
  compbg <- col2rgb(bg, alpha = F)
  compcol <- rbind(compcol[4,]*compcol[1:3,,drop=F]/255, 255-compcol[4,,drop = F])
  rescol <- compcol[1:3,,drop=F] + (compcol[4,]*compbg)/255
  return(rgb(t(rescol), maxColorValue = 255))
}
library(extrafont)
# font_import()
# fonts()
# fonttable()
loadfonts(device="postscript")
loadfonts()

# the aim is to use the 'drc' package to fit models to data and then extract the data and use for plotting in ggplot
# the data could be growth curves, `Dose [uM]`-response, Michaelis-Menten etc.
# here, the S.alba data from drc is used for `Dose [uM]`-response

require(tidyverse)
require(broom)
require(drc)
require(modelr)
# attach(S.alba) # the data used in this gist
# library(egg) 

basedir <- '/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/res/20180413_5fu/'

df <- read_tsv(file.path(basedir, 'dosedata.tsv'))
# df <- df %>% filter(Replicate!=3)
df <- df %>% mutate(Prediction = factor(Prediction, levels = c('Sensitive', 'Resistant'), ordered = T))
df3 <- df %>%
  group_by(Prediction, `Dose [uM]`, `Cell line`) %>%
  summarise(mean(Viability), min(Viability), max(Viability))

# define drm function to use with map
drm.func <- function(x) {
  drm(Viability ~ `Dose [uM]`, 
      fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")), 
      data = x)
}

predict.fun <- function(x) {
  add_predictions(data.frame(`Dose [uM]` = 10^seq(log10(min(df$`Dose [uM]`)), log10(max(df$`Dose [uM]`)), length.out = 1000), check.names = F), x)
}

coefs.fun <- function(x) {coef(x) %>% tidy}

df2 <- df %>% group_by(Prediction, `Cell line`) %>% nest() %>%
  mutate(drmod = map(data, drm.func), 
         pred = map(drmod, predict.fun),
         coefs = map(drmod, coefs.fun))

pal <- c(tumred, tumblue)
names(pal) <- c('Sensitive', 'Resistant')
lty <- c('solid', 'dotted', 'dashed', 'dotdash')
names(lty) <- c('CCRFCEM', 'HCC2998', 'RPMI8226', 'KM12')


my_theme <- function (base_size = 12, base_family = 'TUM Neue Helvetica 55 Regular') # , base_family="Arial"
{
  require(grid)
  theme(text = element_text(size = base_size, face
                            = "plain",
                            colour = "black", hjust = 0.5, vjust = 0.5,
                            angle = 0, lineheight = 1), # , family = base_family
        line = element_line(),
        axis.text  = element_text(size = rel(1), colour="black"),
        axis.text.x = element_text(margin = margin(t = 3, r = 0, b = 1,l = 0, unit = 'pt'), vjust = 0),
        axis.ticks = element_line(colour = "black", size = 0.4, lineend = 'round'),
        axis.ticks.length = unit(2, units = 'pt'),
        axis.title.y = element_text(margin = margin(0,2,0,0, unit = 'pt'), size = rel(1), vjust=0.5,angle=90),
        axis.title.x = element_text(margin = margin(2,0,0,0,unit = 'pt'), size = rel(1), hjust=0.5, vjust=0.5),
        # axis.line = element_line(lineend = 'round', size = 0.4),
        legend.key = element_blank(),
        legend.key.size = unit(1.5, units = 'lines'),
        legend.title = element_text(size = rel(1), face="plain"),
        legend.text = element_text(size = rel(0.7), face="plain"),
        legend.spacing = unit(5,"pt"),
        legend.background = element_blank(),
        legend.position = c(0.4,0.33), # 
        legend.direction = 'vertical',
        legend.box.background = element_blank(),
        legend.box.margin = margin(0.2,0,0,-3, unit = 'lines'),
        legend.box.spacing = unit(3, units = 'pt'),
        legend.box = 'vertical',
        legend.margin = margin(0,0,0,2),
        plot.title = element_text(size = rel(1), colour = "black",face="plain", vjust = 0.5, hjust = 0.5, margin = margin(0,0,4,0, unit = 'pt')),
        panel.background = element_rect(fill = "white",colour = NA),
        panel.border = element_rect(fill = NA, colour = "black", size = 0.6),
        # panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(3,14,5,1,unit = 'pt'),
        strip.background = element_rect(fill = "grey90", colour = "grey50", size=0.3)
  )
}

setEPS()
postscript(file.path(basedir, 'drm_new.eps'), width = 5.8 /2.539998, height = 5.8/2.539998, pointsize = 10, family = 'TUM Neue Helvetica 55 Regular', encoding = 'Greek')
# plot raw data, model and ED50 line
df2 %>% unnest(data) %>% 
  ggplot() + 
  # geom_vline(aes(xintercept = log10(x), color = Prediction, linetype = `Cell line`),
  #            
  #            data = df2 %>% unnest(coefs) %>% filter(names == "ED50:(Intercept)"), show.legend = F) +
  geom_line(aes(log10(`Dose [uM]`*1000), pred, color = Prediction, linetype = `Cell line`), data =
              df2 %>% unnest(pred)) +
  geom_errorbar(data = df3, aes(log10(`Dose [uM]`*1000), ymin = `min(Viability)`, ymax = `max(Viability)`, color = Prediction), linetype = 'solid', width = .2) +
  scale_colour_manual(values = pal) +
  scale_linetype_manual(values = lty) +
  labs(title = 'Validation of\n5-FU Predictions', x = bquote('log'[10]*'(Dose [nM])'), y = 'Relative viability') +
  guides(colour = guide_legend(ncol = 1, title.position = "top", title.hjust = 0, keyheight = 0.5, default.unit = 'lines')) +
  guides(linetype = guide_legend(ncol = 1, title.position = "top", title.hjust = 0, keyheight = 0.5, default.unit = 'lines')) +
  my_theme(base_size = 10)
dev.off()

ED50s <- sapply(df2$drmod, function(i) ED(i, respLev = 50)[1])
names(ED50s) <- df2$`Cell line`

# pdf(file.path(basedir, 'drm.pdf'), width = 7.3509 /2.539998, height = 5.0474/2.539998, pointsize = 12, family = 'TUM Neue Helvetica 55 Regular', useDingbats = F)
# df2 %>% unnest(data) %>% 
#   ggplot() + 
#   geom_vline(aes(xintercept = log10(x), color = Prediction),
#              linetype = 5,
#              data = df2 %>% unnest(coefs) %>% filter(names == "ED50:(Intercept)"), show.legend = F) +
#   geom_line(aes(log10(`Dose [uM]`), pred, color = Prediction), data =
#               df2 %>% unnest(pred)) +
#   geom_errorbar(data = df3, aes(log10(`Dose [uM]`), ymin = `min(Viability)`, ymax = `max(Viability)`, color = Prediction), width = .2) +
#   scale_colour_manual(values = pal) +
#   labs(title = 'Combination Prediction', x = bquote('log10(Dose ['*mu*'M])'), y = 'Relative viability') +
#   guides(colour = guide_legend(nrow = 2, title.position = "left", title.hjust = 0, keyheight = 1, default.unit = 'lines')) +
#   my_theme(base_size = 12)
# dev.off()
# embed_fonts(file = file.path(basedir, 'drm.pdf'), outfile = file.path(basedir, 'drm_embed.pdf'))

# summary of coefficients
df2 %>% unnest(coefs) %>% spread(names, x)
