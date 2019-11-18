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
my_theme <- function (base_size = 12, base_family = 'TUM Neue Helvetica 55 Regular') # , base_family="Arial"
{
  require(grid)
  theme(text = element_text(size = base_size, face
                            = "plain",
                            colour = "black", hjust = 0.5, vjust = 0.5,
                            angle = 0, lineheight = 1), # , family = base_family
        axis.text  = element_text(size = rel(1), colour="black"),
        axis.ticks = element_line(colour = "black", lineend = 'round'),
        axis.ticks.length = unit(2, units = 'pt'),
        axis.title.y = element_text(margin = margin(0,2,0,0, unit = 'pt'), size = rel(1), vjust=0.5,angle=90),
        axis.title.x = element_text(margin = margin(2,0,0,0,unit = 'pt'), size = rel(1), hjust=0.5, vjust=0.5),
        legend.key = element_blank(),
        legend.key.size = unit(0.5, units = 'lines'),
        legend.title = element_text(size = rel(1), face="plain"),
        legend.text = element_text(size = rel(0.8), face="plain"),
        legend.spacing = unit(0,"cm"),
        legend.background = element_blank(),
        legend.position = 'bottom', # c(0.16,0.33)
        legend.direction = 'horizontal',
        legend.box.background = element_blank(),
        legend.box.margin = margin(0.5,0,0,-1.6, unit = 'lines'),
        legend.box.spacing = unit(0, units = 'pt'),
        legend.margin = margin(0,0,0,2),
        plot.title = element_text(size = rel(1), colour = "black",face="plain", vjust = 0.5, hjust = 0.5, margin = margin(0,0,4,0, unit = 'pt')),
        plot.margin = margin(3,14,5,1,unit = 'pt'),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", size = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey90", colour = "grey50", size=0.3)
  )
}
# the aim is to use the 'drc' package to fit models to data and then extract the data and use for plotting in ggplot
# the data could be growth curves, `Dose [uM]`-response, Michaelis-Menten etc.
# here, the S.alba data from drc is used for `Dose [uM]`-response

library(tidyverse)
library(broom)
library(drc)
library(modelr)
# attach(S.alba) # the data used in this gist
# library(egg) 

basedir <- '/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/res/20180413_tepotinib_mk/'

df <- read_tsv(file.path(basedir, 'dosedata.tsv'))
df <- df %>% filter(Replicate!=3)
df <- df %>% mutate(Treatment = factor(Treatment, levels = c('Tepotinib', 'MK-2461', 'Combination (Tepotinib-scale)', 'Combination (MK-2461-scale)'), labels = c('Tepotinib', 'MK-2461', 'Combination (relative to Tepotinib)', 'Combination (relative to MK-2461)'), ordered = T))
df3 <- df %>%
  group_by(Treatment, `Dose [uM]`, `Cell line`) %>%
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

df2 <- df %>% group_by(Treatment) %>% nest() %>%
  mutate(drmod = map(data, drm.func), 
         pred = map(drmod, predict.fun),
         coefs = map(drmod, coefs.fun))

pal <- c(al2hex(alpha(tumblue, alpha = 0.5)), tumblue, al2hex(alpha(tumred, alpha = 0.5)), tumred)
names(pal) <- c('Tepotinib', 'Combination (relative to Tepotinib)', 'MK-2461', 'Combination (relative to MK-2461)')

setEPS()
postscript(file.path(basedir, 'drm.eps'), width = 8.7 /2.539998, height = (6.5433)/2.539998, pointsize = 10, family = 'TUM Neue Helvetica 55 Regular', encoding = 'Greek')
# plot raw data, model and ED50 line
df2 %>% unnest(data) %>% 
  ggplot() + 
  # geom_vline(aes(xintercept = log10(x), color = Treatment),
  #            linetype = 5,
  #            data = df2 %>% unnest(coefs) %>% filter(names == "ED50:(Intercept)"), show.legend = F) +
  geom_line(aes(log10(`Dose [uM]`*1000), pred, color = Treatment), data =
              df2 %>% unnest(pred)) +
  geom_errorbar(data = df3, aes(log10(`Dose [uM]`*1000), ymin = `min(Viability)`, ymax = `max(Viability)`, color = Treatment), width = .2) +
  scale_colour_manual(values = pal) +
  labs(title = 'Combination treatment of HDC-8 cells', x = bquote('log'[10]*'(Dose [nM])'), y = 'Relative\nviability') + # \ntargeting MET and MST1R
  guides(colour = guide_legend(nrow = 2, title.position = "top", title.hjust = 0.5, keyheight = 1, default.unit = 'lines')) +
  my_theme(base_size = 10)
dev.off()

# pdf(file.path(basedir, 'drm.pdf'), width = 7.3509 /2.539998, height = 5.0474/2.539998, pointsize = 12, family = 'TUM Neue Helvetica 55 Regular', useDingbats = F)
# df2 %>% unnest(data) %>% 
#   ggplot() + 
#   geom_vline(aes(xintercept = log10(x), color = Treatment),
#              linetype = 5,
#              data = df2 %>% unnest(coefs) %>% filter(names == "ED50:(Intercept)"), show.legend = F) +
#   geom_line(aes(log10(`Dose [uM]`), pred, color = Treatment), data =
#               df2 %>% unnest(pred)) +
#   geom_errorbar(data = df3, aes(log10(`Dose [uM]`), ymin = `min(Viability)`, ymax = `max(Viability)`, color = Treatment), width = .2) +
#   scale_colour_manual(values = pal) +
#   labs(title = 'Combination treatment', x = bquote('log10(Dose ['*mu*'M])'), y = 'Relative viability') +
#   guides(colour = guide_legend(nrow = 2, title.position = "left", title.hjust = 0, keyheight = 1, default.unit = 'lines')) +
#   my_theme(base_size = 12)
# dev.off()
# embed_fonts(file = file.path(basedir, 'drm.pdf'), outfile = file.path(basedir, 'drm_embed.pdf'))

# summary of coefficients
df2 %>% unnest(coefs) %>% spread(names, x)
