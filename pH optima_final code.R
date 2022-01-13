library(dplyr)
library(fANCOVA)
library(ggplot2)
library(nlme)
library(pastecs)
library(interactions)
library(leaps)
library(patchwork)
library(ggpubr)
library(purrr)

#### LOADING IN DATA FILES ####
### Raw data ###
setwd("~/Dropbox")
dat.all <- read.csv("Wade Margenot/pH optimization/Data and Quantitative Analyses/Processed data for R.csv")
dat.all$Soil.id <- factor(dat.all$Soil.id, levels=c("Soil 1", "Soil 2", "Soil 3", "Soil 4", "Soil 5", "Soil 6", "Soil 7", "Soil 8", "Soil 9", "Soil 10", "Soil 11", "Soil 12", "Soil 13", "Soil 14", "Soil 15", "Soil 16", "Soil 17", "Soil 18", "Soil 19", "Soil 20", "Soil 21", "Soil 22", "Soil 23", "Soil 24", "Soil 25", "Soil 26", "Prunus dulcis", "Aspergillus niger", "Escherichia coli", "Wheat germ"))
dat.all$Rel.Act <- dat.all$Activity.std*100
dat.phos <- subset(dat.all, Enzyme=="Phos")
dat.glu <- subset(dat.all, Enzyme=="B.Glu")

### Loess models and extraction of pH maxes ###
## Standardized values by pH: Fit loess curve and extract max values ###
#NB: there has to be a better way to do this than copy/paste, but couldn't get the damn for loop to work
levels(dat.phos$Soil.id)
DF <- subset(dat.phos, Soil.id=="Wheat germ")  # Change data first!!!

std.span.model <- loess.as(DF$Buffer.pH, DF$Activity.std, criterion="aicc", plot=TRUE)
std.span.pred <- std.span.model[[11]][[1]]

std.model <- loess(Activity.std ~ Buffer.pH, data=DF, span=std.span.pred)
x.range <- range(DF$Buffer.pH)
x.seq <- seq(from=x.range[1], to=x.range[2], length=200)
std.pred <- predict(std.model, newdata = data.frame(Buffer.pH = x.seq), se=TRUE)
Y <- std.pred$fit
CI <- std.pred$se.fit * qt(0.95 / 2 + .5, std.pred$df)
Y.min <- Y - CI
Y.max <- Y + CI
loess.DF.std <- round(data.frame(x = x.seq, Y, Y.min, Y.max),4)

clipr::write_clip(loess.DF.std) #copy/paste to excel file cuz I suck at coding - becomes the 'modeled.data.all' file

### Load and wrangle the modeled data ###
modeled.data.all <- read.csv("/Users/TheDankness/Dropbox/Wade Margenot/pH optimization/Data and Quantitative Analyses/Standardized data for R.csv")
modeled.data.all$Rel.Act_est <- modeled.data.all$Activity.std_est*100
modeled.data.all$Rel.Act_min <- modeled.data.all$Activity.std_min*100
modeled.data.all$Rel.Act_max <- modeled.data.all$Activity.std_max*100

modeled.data.all$Soil.id <- factor(modeled.data.all$Soil.id, levels=c("Soil 1", "Soil 2", "Soil 3", "Soil 4", "Soil 5", "Soil 6", "Soil 7", "Soil 8", "Soil 9", "Soil 10", "Soil 11", "Soil 12", "Soil 13", "Soil 14", "Soil 15", "Soil 16", "Soil 17", "Soil 18", "Soil 19", "Soil 20", "Soil 21", "Soil 22", "Soil 23", "Soil 24", "Soil 25", "Soil 26", "Prunus dulcis", "Aspergillus niger", "Escherichia coli", "Wheat germ"))

modeled.data <- subset(modeled.data.all, Type=="Soil")
modeled.pure.data <- subset(modeled.data.all, Type=="Pure")

### Summarize soil characteristics ###
soil.sum <- as.data.frame(dat.all[,c(2:9)] %>% dplyr::group_by(Enzyme, Soil.id, Soil.ID, pH.category, Clay.category) %>% dplyr::summarize(Clay.pct=mean(Clay.pct), SOC.pct=mean(SOC.pct), soil.pH=mean(soil.pH)))
soil.sum <- soil.sum[-c(27,28,55,56),]
NRCS.clim.sum <- read.csv("/Users/TheDankness/Dropbox/Wade Margenot/pH optimization/Data and Quantitative Analyses/Site data final (long).csv")
NRCS.clim.dat <- NRCS.clim.sum[,-1]
soil.summary <- cbind(soil.sum, NRCS.clim.dat)

pH.vec <- as.data.frame(soil.summary %>% subset(Enzyme=="B.Glu") %>% dplyr::select(Soil.id, soil.pH))


### Figure 1: activity profiles ###
Fig1 <- ggplot(data=modeled.data.all, aes(x=pH_est, y=Rel.Act_est, group=Enzyme, color=Enzyme)) + geom_ribbon(aes(ymin=Rel.Act_min, ymax=Rel.Act_max, color=NULL), alpha=0.3) + geom_point(size=0.18) + facet_wrap(.~Soil.id, nrow=6, ncol=5) + scale_color_manual(values=c("#005AB5", "#DC3220")) + theme_bw() + theme(legend.position="none") + labs(y="Relative Activity (%)", x="Buffer pH") + scale_y_continuous(breaks=seq(0,100, by=25)) + scale_x_continuous(breaks=seq(3,13,2)) + geom_vline(aes(xintercept=soil.pH), data=pH.vec, color="black", lty="dashed") + coord_cartesian(ylim=c(0,100), xlim=c(3,12))
Fig1
ggsave("Figure 1 - activity profiles.pdf", width=7, height=7.75)

avg.loess.plot <- ggplot(data=modeled.data.all, aes(x=pH_est, y=Activity.std_est, group=Soil.id)) + scale_x_continuous(limits=c(3,13), breaks=seq(3,13,2)) + geom_line(aes(group=Soil.id), colour="#c9c9c9", alpha=0.3, size=2) + theme_bw() + facet_wrap(~Enzyme, nrow=2) + theme(legend.position="none") + stat_summary(aes(group=Enzyme, color=Enzyme), fun.y=mean, geom="line", size=1.5) + scale_color_manual(values=c("#005AB5", "#DC3220"))
avg.loess.plot # export at 500 x 527

avg.loess.plot + (glu.pH.optima.dp / phos.pH.optima.dp)


data.phos <- subset(data.all, Enzyme=="Phos")
data.glu <- subset(data.all, Enzyme=="B.Glu")



span.model <- loess.as(data.all$Buffer.pH, data.all$Activity, criterion="aicc", plot=FALSE)
summary(span.model)
span.model$pars$span
span.pred <- span.model[[11]][[1]]


#### PREDEFINED PALETTES ####
cbPalette_grey <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#F0E442")
cbPalette_JustColors <- c("#E69F00", "#56B4E9", "#009E73", "#D55E00")
cbPalette_black <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")