setwd('J:/Timelapses_data/Timelapses_glass_2021/MatlabAnalysis_210209_210215_210222_210301/Table folder')

library(gridExtra)
library(ggplot2)
library(lifecycle)
library(dplyr)
library(grid)
library(cowplot)


CT = read.csv("CT.csv",sep=",",header = TRUE)
IL2 = read.csv("IL2.csv",sep=",",header = TRUE)
Ttz = read.csv("Ttz.csv",sep=",",header = TRUE)
Il2Ttz = read.csv("Il2Ttz.csv",sep=",",header = TRUE)

CT$Area = CT$Area * 0.3225^2 ; IL2$Area = IL2$Area * 0.3225^2; Ttz$Area = Ttz$Area * 0.3225^2; Il2Ttz$Area = Il2Ttz$Area * 0.3225^2
CT$Bin = factor(CT$Bin,levels = c("IN", "EDGE", "OUT", 'NaN'))
IL2$Bin = factor(IL2$Bin,levels = c("IN", "EDGE", "OUT", 'NaN'))
Ttz$Bin = factor(Ttz$Bin,levels = c("IN", "EDGE", "OUT", 'NaN'))
Il2Ttz$Bin = factor(Il2Ttz$Bin,levels = c("IN", "EDGE", "OUT", 'NaN'))


MERGE = rbind(CT,IL2,Ttz,Il2Ttz)
MERGE$Exp = c(rep('CT',nrow(CT)),rep('IL2',nrow(IL2)),rep('Ttz',nrow(Ttz)),rep('Il2Ttz',nrow(Il2Ttz)))
MERGE$Exp = as.factor(MERGE$Exp)
MERGE$Bin = as.factor(MERGE$Bin)
# Process mean Areas
MERGE$Mean_Area = rep(0, nrow(MERGE))
for(i in unique(MERGE$File)){
  for (j in unique(MERGE$ID)){
    MERGE$Mean_Area[MERGE$File == i & MERGE$ID == j] = mean(MERGE$Area[MERGE$File == i & MERGE$ID == j])
  }
}
MERGEsorted = with(MERGE, MERGE[order(Exp, Bin, File, ID),])  #ordered by experiment, region, file, cell ID
write.csv(MERGEsorted,'J:/Timelapses_data/Timelapses_glass_2021/MatlabAnalysis_210209_210215_210222_210301/Table folder/Data_sorted.csv')
MERGEu = distinct(MERGEsorted, File, ID,.keep_all = TRUE)    #keep only the mean values = keep 1 row per cell
write.csv(MERGEu,'J:/Timelapses_data/Timelapses_glass_2021/MatlabAnalysis_210209_210215_210222_210301/Table folder/MeanData_sorted.csv')

nrow(distinct(MERGE, File))
nrow(MERGEu[MERGEu$ID == 1,])

# Polar histogram
ggplot(CT, aes(x = Theta)) +   coord_polar(theta = "x", start = 0, direction = +1) +
  geom_histogram(binwidth = 10) + 
  scale_x_continuous(breaks = seq(0, 90, 10), limits = c(0,360)) 

ggplot(IL2, aes(x = Theta)) +   coord_polar(theta = "x", start = 0, direction = +1) +
  geom_histogram(binwidth = 10) + 
  scale_x_continuous(breaks = seq(0, 90, 10), limits = c(0,360))

ggplot(Ttz, aes(x = Theta)) +   coord_polar(theta = "x", start = 0, direction = +1) +
  geom_histogram(binwidth = 10) + 
  scale_x_continuous(breaks = seq(0, 90, 10), limits = c(0,360))

ggplot(Il2Ttz, aes(x = Theta)) +   coord_polar(theta = "x", start = 0, direction = +1) +
  geom_histogram(binwidth = 10) + 
  scale_x_continuous(breaks = seq(0, 90, 10), limits = c(0,360))




# Violins
Vio1 = ggplot(MERGE[MERGE$Bin != 'NaN',], aes(x=Bin[Bin != 'NaN'], y=Vi[Bin != 'NaN'], fill = Exp[Bin != 'NaN'])) + 
  geom_violin() + xlab('Region') + ylab('Instant Velocities') + ylim(0,20)
Vio2 = ggplot(MERGEu[MERGEu$Bin != 'NaN',], aes(x=Bin[Bin != 'NaN'], y=Vm[Bin != 'NaN'], fill = Exp[Bin != 'NaN'])) + 
  geom_violin() + xlab('Region')  + ylab('Mean Velocities') + ylim(0,15)
Vio3 = ggplot(MERGE[MERGE$Bin != 'NaN',], aes(x=Bin[Bin != 'NaN'], y=Area[Bin != 'NaN'], fill = Exp[Bin != 'NaN'])) + 
  geom_violin() + xlab('Region')+ ylab('Instant Areas') + ylim(0,150)
Vio4 = ggplot(MERGEu[MERGEu$Bin != 'NaN',], aes(x=Bin[Bin != 'NaN'], y=Mean_Area[Bin != 'NaN'], fill = Exp[Bin != 'NaN'])) + 
  geom_violin() + xlab('Region') + ylab('Mean Areas') + ylim(0,75)
lgd = get_legend(Vio1 + labs(fill = 'Conditions'))
VIOLINS = plot_grid(
          plot_grid(Vio1 + theme(legend.position="none"), Vio2 + theme(legend.position="none"), 
          Vio3 + theme(legend.position="none"), Vio4 + theme(legend.position="none"), 
          ncol = 2, nrow = 2, labels = c("A.1","A.2","B.1","B.2"),label_size = 12), lgd, ncol = 2, nrow = 1, rel_widths = c(.8, .1))
VIOLINS = plot_grid( #only the means
          plot_grid(Vio2 + theme(legend.position="none"), Vio4 + theme(legend.position="none"), 
            ncol = 1, nrow = 2, labels = c("A.1","A.2","B.1","B.2"),label_size = 12), lgd, ncol = 2, nrow = 1, rel_widths = c(.8, .1))

# Histograms - Vi/Vm - CT/DRUGS
plot_grid(
ggplot(MERGE[MERGE$Bin != 'NaN',], aes(x=Vi[Bin != 'NaN'], fill = Exp[Bin != 'NaN'])) + 
  geom_histogram(aes(y=..density..), bins = 100, position = "identity", alpha =.5) + geom_density(aes(colour = Exp[Bin != 'NaN']), alpha = .1) + facet_wrap(vars(Bin[Bin != 'NaN'])),

ggplot(MERGEu[MERGEu$Bin != 'NaN',], aes(x=Vm[Bin != 'NaN'], fill = Exp[Bin != 'NaN'])) + 
  geom_histogram(aes(y = ..density..), bins = 25, position = "identity", alpha =.5) + geom_density(aes(colour = Exp[Bin != 'NaN']), alpha = .1) + facet_wrap(vars(Bin[Bin != 'NaN'])),
ncol = 1, nrow = 2) 



# 2D histogram, density, binned density - Vi=f(Area)
a = ggplot(CT, aes( y = Vi, x = Area))
b = ggplot(IL2, aes( y = Vi, x = Area))
c = ggplot(Ttz, aes( y = Vi, x = Area))
d = ggplot(Il2Ttz, aes( y = Vi, x = Area))
grid.arrange(
  a + geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + xlim(0, 200),
  a +   stat_density_2d(aes(fill = ..density..), geom = "raster", na.rm = TRUE, contour = FALSE) + 
    scale_fill_continuous(type = "viridis") + xlim(0, 50) + ylim(0,10),
  a + geom_density2d_filled(na.rm = TRUE) + xlim(0, 50) + ylim(0,10),
  
  b + geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + xlim(0, 200),
  b +   stat_density_2d(aes(fill = ..density..), geom = "raster", na.rm = TRUE, contour = FALSE) + 
    scale_fill_continuous(type = "viridis") + xlim(0, 50) + ylim(0,10),
  b +   geom_density2d_filled(na.rm = TRUE) + xlim(0, 50) + ylim(0,10),

  c + geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + xlim(0, 200),
  c +   stat_density_2d(aes(fill = ..density..), geom = "raster", na.rm = TRUE, contour = FALSE) + 
    scale_fill_continuous(type = "viridis") + xlim(0, 50) + ylim(0,10),
  c +   geom_density2d_filled(na.rm = TRUE) + xlim(0, 50) + ylim(0,10),

  d + geom_bin2d(bins = 100) + scale_fill_continuous(type = "viridis") + xlim(0, 200),
  d +   stat_density_2d(aes(fill = ..density..), geom = "raster", na.rm = TRUE, contour = FALSE) + 
    scale_fill_continuous(type = "viridis") + xlim(0, 50) + ylim(0,10),
  d +   geom_density2d_filled(na.rm = TRUE) + xlim(0, 50) + ylim(0,10),
   ncol = 3, nrow = 4) 

# 2D density Vi=f(Area) - Each region - CT/DRUGS
a = ggplot(MERGE[MERGE$Bin != 'NaN',], aes(x = Area[Bin != 'NaN'], y = Vi[Bin != 'NaN'])) + geom_density2d_filled(contour_var = "ndensity", na.rm = TRUE) + xlim(0, 50) + ylim(0,10)+
  xlim(0, 50) + ylim(0,10) + labs( x = expression(Instant~Cell~Area~(µm^2)), y = expression(Instant~Cell~Velocity~(µm.min^{-1}))) +
  facet_grid(Exp[Bin != 'NaN'] ~ Bin[Bin != 'NaN']) +
  theme_bw(base_size = 11, base_line_size = 1, base_rect_size = 1) +
  theme(panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines"), plot.title = element_text(face = "bold"))   
  #+ ggsave('MigMorpho.png',plot = last_plot())

# 2D density Vm=f(Aream) - Each region - CT/DRUGS
b = ggplot(MERGEu[MERGEu$Bin != 'NaN',], aes(x = Mean_Area[Bin != 'NaN'], y = Vm[Bin != 'NaN'])) + geom_density2d_filled(contour_var = "ndensity", na.rm = TRUE) + 
  xlim(0, 50) + ylim(0,10) + labs( x = expression(Mean~Cell~Area~(µm^2)), y = expression(Mean~Cell~Velocity~(µm.min^{-1}))) +
  facet_grid(Exp[Bin != 'NaN'] ~ Bin[Bin != 'NaN']) + 
  theme_bw(base_size = 11, base_line_size = 1, base_rect_size = 1) +
  theme(panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines"), plot.title = element_text(face = "bold"))  
  #+ ggsave('MigMorpho.png',plot = last_plot(), heigth = 2, width = 3, units = scale(""))

lgd = get_legend(a + labs(fill = 'Density'))
CORPLOTS = plot_grid(plot_grid(a + theme(legend.position="none"),b+ theme(legend.position="none"), ncol = 1, nrow= 2, labels = c("C.1", "C.2"), label_size = 12),
          lgd, ncol = 2, nrow = 1, rel_widths = c(.8, .2))


plot_grid(VIOLINS, CORPLOTS, nrow = 1, ncol = 2)


# Cells roundess, aspect ratio, circularity - histograms
plot_grid(
ggplot(CT[CT$Bin != 'NaN',], aes(x = Circularity[Bin != 'NaN'])) + 
  geom_histogram(aes(fill = Bin[Bin != 'NaN']), bins = 50) + labs(title='CT'),
ggplot(DRUGS[DRUGS$Bin != 'NaN',], aes(x = Circularity[Bin != 'NaN'])) + 
  geom_histogram(aes(fill = Bin[Bin != 'NaN']), bins = 50) + labs(title='DRUGS'),
ggplot(CT[CT$Bin != 'NaN',], aes(x = Ri[Bin != 'NaN'])) + 
  geom_histogram(aes(fill = Bin[Bin != 'NaN']), bins = 50) + labs('CT'),
ggplot(DRUGS[DRUGS$Bin != 'NaN',], aes(x = Ri[Bin != 'NaN'])) + 
  geom_histogram(aes(fill = Bin[Bin != 'NaN']), bins = 50) + labs('DRUGS'),
ggplot(CT[CT$Bin != 'NaN',], aes(x = AFi[Bin != 'NaN'])) + 
  geom_histogram(aes(fill = Bin[Bin != 'NaN']), bins = 50) + labs('CT'),
ggplot(DRUGS[DRUGS$Bin != 'NaN',], aes(x = AFi[Bin != 'NaN'])) + 
  geom_histogram(aes(fill = Bin[Bin != 'NaN']), bins = 50) + labs('DRUGS'),
ncol = 2, nrow = 3)

# Aspect ratio, roundess, each bin - maps
plot_grid(
  ggplot(MERGE[MERGE$Bin != 'NaN',], aes(x = AFi[Bin != 'NaN'], y = Vi[Bin != 'NaN'])) + geom_density2d_filled(contour_var = "ndensity", na.rm = TRUE) + 
    facet_grid(Exp[Bin !='NaN'] ~ Bin[Bin != 'NaN']) + ylim(0,10),
  ggplot(MERGE[MERGE$Bin != 'NaN',], aes(x = Ri[Bin != 'NaN'], y = Vi[Bin != 'NaN'])) + geom_density2d_filled(contour_var = "ndensity", na.rm = TRUE) + 
    facet_grid(Exp[Bin !='NaN'] ~ Bin[Bin != 'NaN']) + ylim(0,10),
  ncol = 1, nrow = 2)

# Same with mean values 
plot_grid(
  ggplot(MERGEu[MERGEu$Bin != 'NaN',], aes(x = AFm[Bin != 'NaN'], y = Vm[Bin != 'NaN'])) + geom_density2d_filled(contour_var = "ndensity", na.rm = TRUE) + 
    facet_grid(Exp[Bin !='NaN'] ~ Bin[Bin != 'NaN']) + ylim(0,10),
  ggplot(MERGEu[MERGEu$Bin != 'NaN',], aes(x = Rm[Bin != 'NaN'], y = Vm[Bin != 'NaN'])) + geom_density2d_filled(contour_var = "ndensity", na.rm = TRUE) + 
    facet_grid(Exp[Bin !='NaN'] ~ Bin[Bin != 'NaN']) + ylim(0,10),
  ncol = 1, nrow = 2)


# Orientation - histogram
plot_grid(
  ggplot(CT[CT$Bin != 'NaN',], aes(x = Orientation[Bin != 'NaN'])) + geom_histogram(aes(fill = Bin[Bin != 'NaN'])),
  ggplot(DRUGS[DRUGS$Bin != 'NaN',], aes(x = Orientation[Bin != 'NaN'])) + geom_histogram(aes(fill = Bin[Bin != 'NaN'])), 
  ncol = 1, nrow = 2)
# angle are toward the x-axis -> need to convert toward the axis of each spheroid ...


library(stringr)
# Morpho params of cells entering the spheroid
CROSS = MERGE[FALSE,]
for(i in unique(MERGE$File)){
  for (j in unique(MERGE$ID)){
    if ('IN' %in% unique(MERGE$Bin[ MERGE$File == i & MERGE$ID == j])){
      pos = which(MERGE$Bin[MERGE$ID == j & MERGE$File == i] == 'EDGE')
      CROSS = rbind(CROSS, 
                    MERGE[MERGE$ID == j & MERGE$File == i,][which(MERGE$Bin[MERGE$ID == j & MERGE$File == i][pos+1]=='IN'),])
    }
 }
}

unique(MERGE$Bin[MERGE$ID == 70 & MERGE$File == i])

MERGE$Bin[MERGE$ID == 48 & MERGE$File == i]
pos = unique(MERGE$Bin[MERGE$ID == 48 & MERGE$File == i])
MERGE$Bin[MERGE$ID == 48 & MERGE$File == i] == pos[1]


MERGE$Bin[MERGE$ID == 39 & MERGE$File == i]
unique(MERGE$Bin[MERGE$ID == 48 & MERGE$File == i])
unique(MERGE$Bin[MERGE$ID == 39 & MERGE$File == i])







# Density plots for Alice
grid.arrange(
ggplot(CT, aes(x = Vm)) + geom_histogram(bins = 150) + xlim(0,15),
ggplot(CT, aes(x = Vm)) + geom_histogram(bins = 500, aes( y = ..density..)) + geom_density(fill = "red", alpha = .3) + xlim(0,15) + ylim(0,0.3),
ggplot(CT, aes(x = Vm)) + geom_histogram(bins = 250, aes( y = ..density..)) + geom_density(fill = "red", alpha = .3) + xlim(0,15)+ ylim(0,0.3),
ggplot(CT, aes(x = Vm)) + geom_histogram(bins = 20, aes( y = ..density..)) + geom_density(fill = "red", alpha = .3) + xlim(0,15)+ ylim(0,0.3),
ncol = 1, nrow = 4)



