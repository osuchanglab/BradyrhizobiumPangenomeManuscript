######################
# Necessary Packages #
######################

library(ggplot2)
library(reshape2)
library(gridExtra)
library(scales)

#####################
# Generic functions #
#####################

# Not %in%
`%!in%` = Negate(`%in%`)

# Standard error with NA removal
se <- function(x) sd(x, na.rm=T)/sqrt(length(x))

## Regus 2015 equation for scaled RG
# x.rg is sample biomass
# c.rg is control biomass
calculate.RGD <-  function(x.rg, c.rg){
  return (((x.rg - c.rg) / c.rg ) * 100)
}

## Wendlandt 2019 for RG ratio
# x.rg is sample biomass
# c.rg is control biomass
calculate.RG <-  function(x.rg, c.rg){
  return (x.rg / c.rg)
}

# Canonical Function For Generating Barplot with SE and datapoints
symbio.plot <- function(dataset, response, factor1, factor2=NULL){
  # Defining necessary function inside of this R function  
  se <- function(x) sd(x, na.rm=T)/sqrt(length(x))
  
  # Generation of aggregated mean and se dataset and merge them
  if (is.null(factor2)){
    temp.dataset <- dataset[, c(response, factor1)]
    colnames(temp.dataset) <- c('y', 'x1')
    agg.mean <- aggregate(y~x1, data=temp.dataset, mean, na.rm=T) 
    agg.se <- aggregate(y~x1, data=temp.dataset, se) 
    colnames(agg.se) [2] <- 'se'
    agg <- merge(agg.mean, agg.se, by='x1')
    # Plotting the aggregated dataset
    ggplot(agg, aes(x=x1, y=y, fill=x1)) + 
      geom_bar(stat='summary',position='dodge', fun.y='mean', na.rm=T) + 
      geom_errorbar(aes( ymin=y-se, ymax=y+se), position = position_dodge(width = 0.90), width=0.2) +
      geom_point(data=temp.dataset, aes(x=x1, color=x2), colour='black', pch=21, position=position_dodge(width=0.9)) +
      #ggtitle('Median Nodules by Nod+ Treatments') +
      xlab(factor1) + ylab(response) +
      facet_grid(~x1, space='free_x', scales='free_x')
    
    
  } else {
    temp.dataset <- dataset[, c(response, factor1, factor2)]
    colnames(temp.dataset) <- c('y', 'x1', 'x2')
    agg.mean <- aggregate(y~x1+x2, data=temp.dataset, mean, na.rm=T)
    agg.se <- aggregate(y~x1+x2, data=temp.dataset, se) 
    colnames(agg.se) [3] <- 'se'
    agg <- merge(agg.mean, agg.se, by=c('x1', 'x2'))
    
    # Plotting the aggregated dataset
    ggplot(agg, aes(x=x1, y=y, fill=x2)) + 
      geom_bar(stat='summary',position='dodge', fun.y='mean', na.rm=T) + 
      geom_errorbar(aes( ymin=y-se, ymax=y+se), position = position_dodge(width = 0.90), width=0.2) +
      geom_point(data=temp.dataset, aes(x=x1, color=x2), colour='black', pch=21, position=position_dodge(width=0.9)) +
      #ggtitle('Median Nodules by Nod+ Treatments') +
      xlab(factor1) + ylab(response) +
      facet_grid(~x1, space='free_x', scales='free_x')
  }
  
}

#########

# Function to generate diagnostic plot of lm model and summary
library(ggfortify)
mod.diagnostics <-  function(model.lm) {
  diagnostics <- ggplot2::autoplot(model.lm) + theme_light()
  grid.arrange(grobs=diagnostics@plots, top=format(model.lm$call$formula))
  summary(model.lm)
}

# Function to generate interaction plots between two variables on response
interaction.plot <- function(dataset, response, x_group1, group2, fun.y=function (x) {x}){
  temp.dataset <- dataset[, c(response, x_group1, group2)]
  colnames(temp.dataset) <- c('y', 'x1', 'x2')
  #temp.dataset
  
  ggplot(temp.dataset, aes(x=factor(x1), y=fun.y(y), group=x2, color=x2)) + 
    stat_summary(fun.y=mean, 
                 fun.ymin = function(x) mean(x) - sd(x),
                 fun.ymax = function(x) mean(x) + sd(x),
                 geom='linerange', position=position_dodge(width=0.3)) + 
    stat_summary(fun.y=mean, geom='line', size=2) +
    ggtitle('Interaction Plot') + xlab(as.character(x_group1)) + ylab(as.character(response)) +
    labs(color=as.character(group2))
}

####################
# Inoculation Data #
####################
innoc <- data.frame('Treatment'=c('2', '131', '156', '187', '186', '4', '200'), 
                    'Mean_Count'=c(4.9E7, 4.3E7, 5.9E7, 2.1E6, 3.3E7, 3.4E7, 6.7E7),
                    'SE'=c(2.5E7, 2E7, 2.9E7, 1.3E6, 1.3E7, 1.4E7, 3.2E7))

head(innoc)

ggplot(innoc, aes(x=Treatment, y=log10(Mean_Count))) + 
  geom_bar(stat='identity', fill='steelblue') +
  geom_errorbar(aes(ymin=log10(Mean_Count - SE), ymax=log10(Mean_Count + SE))) +
  scale_y_continuous(limits=c(4, 9), oob=rescale_none) +
  ggtitle('Mean CFU count from two dilutions, two replications each')

######################
## Dataset formation #
######################

# Loading harvest nodule count and root-shoot biomass datasheet
data.load <- read.csv('Harvest/Host_Generalism_Harvest_final.csv', header = TRUE)
data <- data.load
head(data.load)
biomass.load <- read.csv('Harvest/Host_Generalism_Biomass_Datasheet.csv')
data <- merge(data, biomass.load, by='Plant_ID')


## Checking data equality after merging
data[which(data$Species != data$Host),]
data[which(as.character(data$Treatment.x) != as.character(data$Treatment.y)),]
dim(data.load)
dim(biomass.load)
dim(data)
data <- within(data, rm(Treatment.y))
colnames(data)[4] <- 'Treatment'
head(data)
str(data)
table(data[,c('Treatment', 'Host')])


# Calculate day-post-innoculation or dpi 
data$dpi <- as.Date(data$Date, format='%m/%d/%y') - as.Date('10/12/2019', format='%m/%d/%y')

## Loading Traits
trait <- read.csv('Harvest/trait.csv', header=T)
trait
data <- merge(data, trait, by='Treatment')

##############################
# Calculating relative growth #
##############################

## Biomass transformation to mg 
data$ShootMass_mg <- as.numeric(as.character(data$ShootMass_g))*1000 # g to mg conversion
data$RootMass_mg <- as.numeric(as.character(data$RootMass_g))*1000 # g to mg conversion
data$Biomass_mg <- data$ShootMass_mg + data$RootMass_mg

## RG for Relative Growth Ratio, Wendlandt et al. 2019
data$ShootRG <- NA
data$RootRG <- NA
data$RG <- NA

## RGD for Scaled Relative Growth Difference, Regus et al. 2015?
data$ShootRGD <- NA
data$RootRGD <- NA
data$RGD <- NA


## Calculating control plant mean biomass
control.RG.AcS <- apply(data[which(data$Host == 'AcS' & data$Treatment =='H2O'), ][,c('ShootMass_mg', 'RootMass_mg', 'Biomass_mg')], 2, mean, na.rm=T)
control.RG.AcW <- apply(data[which(data$Host == 'AcW' & data$Treatment =='H2O'), ][,c('ShootMass_mg', 'RootMass_mg', 'Biomass_mg')], 2, mean, na.rm=T)
control.RG.Bic <- apply(data[which(data$Host == 'Bic' & data$Treatment =='H2O'), ][,c('ShootMass_mg', 'RootMass_mg', 'Biomass_mg')], 2, mean, na.rm=T)
control.RG.Nan <- apply(data[which(data$Host == 'Nan' & data$Treatment =='H2O'), ][,c('ShootMass_mg', 'RootMass_mg', 'Biomass_mg')], 2, mean, na.rm=T)
control.RG.Suc <- apply(data[which(data$Host == 'Suc' & data$Treatment =='H2O'), ][,c('ShootMass_mg', 'RootMass_mg', 'Biomass_mg')], 2, mean, na.rm=T)

# Relative Growth Calculation  (Wendlandt 2019)
data[which(data$Host == 'AcS'), ][, c('ShootRG', 'RootRG', 'RG')] <- sweep(data[which(data$Host == 'AcS'), ]
                                                                           [, c('ShootMass_mg', 'RootMass_mg', 'Biomass_mg')], 2, control.RG.AcS, calculate.RG)
data[which(data$Host == 'AcW'), ][, c('ShootRG', 'RootRG', 'RG')] <- sweep(data[which(data$Host == 'AcW'), ]
                                                                           [, c('ShootMass_mg', 'RootMass_mg', 'Biomass_mg')], 2, control.RG.AcW, calculate.RG)
data[which(data$Host == 'Bic'), ][, c('ShootRG', 'RootRG', 'RG')] <- sweep(data[which(data$Host == 'Bic'), ]
                                                                           [, c('ShootMass_mg', 'RootMass_mg', 'Biomass_mg')], 2, control.RG.Bic, calculate.RG)
data[which(data$Host == 'Nan'), ][, c('ShootRG', 'RootRG', 'RG')] <- sweep(data[which(data$Host == 'Nan'), ]
                                                                           [, c('ShootMass_mg', 'RootMass_mg', 'Biomass_mg')], 2, control.RG.Nan, calculate.RG)
data[which(data$Host == 'Suc'), ][, c('ShootRG', 'RootRG', 'RG')] <- sweep(data[which(data$Host == 'Suc'), ]
                                                                           [, c('ShootMass_mg', 'RootMass_mg', 'Biomass_mg')], 2, control.RG.Suc, calculate.RG)

p1 <- ggplot(data, aes(x=Species, y=Biomass_mg)) + 
  geom_boxplot() + 
  geom_jitter(aes(color=Block)) +
  ggtitle('Biomass Distribution by Speceies')
p2 <- ggplot(data, aes(x=Species, y=RG)) + 
  geom_boxplot() + 
  geom_jitter(aes(color=Block)) +
  ggtitle('RG Distribution by Species (Wendlandt 2019)')

# Relative Growth Calculation (Regus 2015)
data[which(data$Host == 'AcS'), ][, c('ShootRGD', 'RootRGD', 'RGD')] <- sweep(data[which(data$Host == 'AcS'), ]
                                                                           [, c('ShootMass_mg', 'RootMass_mg', 'Biomass_mg')], 2, control.RG.AcS, calculate.RGD)
data[which(data$Host == 'AcW'), ][, c('ShootRGD', 'RootRGD', 'RGD')] <- sweep(data[which(data$Host == 'AcW'), ]
                                                                           [, c('ShootMass_mg', 'RootMass_mg', 'Biomass_mg')], 2, control.RG.AcW, calculate.RGD)
data[which(data$Host == 'Bic'), ][, c('ShootRGD', 'RootRGD', 'RGD')] <- sweep(data[which(data$Host == 'Bic'), ]
                                                                           [, c('ShootMass_mg', 'RootMass_mg', 'Biomass_mg')], 2, control.RG.Bic, calculate.RGD)
data[which(data$Host == 'Nan'), ][, c('ShootRGD', 'RootRGD', 'RGD')] <- sweep(data[which(data$Host == 'Nan'), ]
                                                                           [, c('ShootMass_mg', 'RootMass_mg', 'Biomass_mg')], 2, control.RG.Nan, calculate.RGD)
data[which(data$Host == 'Suc'), ][, c('ShootRGD', 'RootRGD', 'RGD')] <- sweep(data[which(data$Host == 'Suc'), ]
                                                                           [, c('ShootMass_mg', 'RootMass_mg', 'Biomass_mg')], 2, control.RG.Suc, calculate.RGD)

p3 <- ggplot(data, aes(x=Treatment, y=RGD)) + 
  geom_boxplot() +
  geom_jitter(aes(color=Block)) + 
  ggtitle('Scaled RG Distribution by Species (Regus 2015 equation)')

grid.arrange(p1, p2, p3, nrow=1)


############################
# Process Nodule Size Data #
############################
nod <- read.csv('Harvest/nodule_area.csv', header = TRUE, sep=',')

nod$plant_id <- as.numeric(nod$plant_id)
nod$nodule_area <- (nod$nodule_area)*100
nod$dryMass <- ((nod$nodule_area)*100 - 0.9097853)/5.525844 # Based on Wendlandt 2019 equation
nodule_biomass.sum <- aggregate(nod[, c('dryMass')], by=list(nod$plant_id), FUN=sum)
nodule_biomass.mean <- aggregate(nod[, c('dryMass')], by=list(nod$plant_id), FUN=mean)
nodule_biomass.area <- aggregate(nod[, c('nodule_area')], by=list(nod$plant_id), FUN=sum)
nodule_biomass.mean_area <- aggregate(nod[, c('nodule_area')], by=list(nod$plant_id), FUN=mean)

colnames(nodule_biomass.sum)[2] <- 'total_nodule_biomass'
colnames(nodule_biomass.mean)[2] <- 'mean_nodule_biomass'
colnames(nodule_biomass.area)[2] <- 'total_nodule_area'
colnames(nodule_biomass.mean_area)[2] <- 'mean_nodule_area'

data <- merge(data, nodule_biomass.sum, by.y='Group.1', by.x='Plant_ID', all.x = T)
data <- merge(data, nodule_biomass.mean, by.y='Group.1', by.x='Plant_ID', all.x = T)
data <- merge(data, nodule_biomass.area, by.y='Group.1', by.x='Plant_ID', all.x = T)
data <- merge(data, nodule_biomass.mean_area, by.y='Group.1', by.x='Plant_ID', all.x = T)
#data <- merge(data, nod, by.y='Group.1', by.x='Plant_ID', all.x = T)

  
###################################
# Process Dry Nodule Biomass Data #
###################################
raw_nodule_biomass <- read.csv('Harvest/dry_nodule_biomass.csv', header=T)
data <- merge(data, raw_nodule_biomass, by=c('Plant_ID', 'Treatment', 'Species'), all.x=T)
dim(data)

## Eyeballing dry nodule biomass with nodule numbers
#investment_dataset <- data[which(data$Nodules != 0), c('Host', 'Treatment', 'Nodules', 'dry_nodule_biomass', 'RG','RGD')]
#investment_summary <- aggregate(investment_dataset, by=list(investment_dataset$Host, investment_dataset$Treatment), 
#                                mean, na.rm=T)
#investment_summary <- melt(investment_summary[, c('Group.1', 'Group.2', 'Nodules', 'dry_nodule_biomass', 'RG', 'RGD')], id.vars=c('Group.1', 'Group.2'))
#head(investment_summary)
#
#ggplot(investment_summary, aes(x=variable, y=value, fill=Group.1)) + geom_bar(stat='identity', position='dodge') + facet_wrap(~Group.2, scales='free')
  


# Calculating Investment
data$investment  <- as.numeric(data$dry_nodule_biomass) / as.numeric(data$Biomass_mg)

p1 <- ggplot(data[which(data$Treatment != 'H2O'),], aes(x=Treatment, y=investment, fill=Host)) + 
  geom_bar(stat='summary', fun.y=mean, na.rm=T, position='dodge') + 
  ggtitle('Investment (Dry Nodule Biomass / Dry Total Biomass)') + ylab('Investment') +
  facet_grid(~Treatment, scales='free_x')
  
p2 <- ggplot(data[which(data$Treatment != 'H2O'),], aes(x=Treatment, y=Nodules, color=Host, size=dry_nodule_biomass, group=Host)) + 
  geom_point(alpha=0.7, position=position_dodge(width=0.9)) +
  scale_size(range=c(.1, 10), name='Nodule Biomass (mg)') +
  ggtitle('Total nodules and dry-biomass') +  
  facet_grid(~Treatment, space='free_x', scales='free_x')

grid.arrange(p1, p2, nrow=2)


###################
# Leaf count data #
###################

leaf.count <-  read.csv(file='~/Sachs/host_generalism/2019GH16/Inoculation/leaf_count_block_data.csv', header=T)
leaf.count <- merge(leaf.count, data, by='Plant_ID')
head(leaf.count)

# The following plot shows that effect of  #true_leaves during inoculation do not exist, probably except for AcS
ggplot(leaf.count, aes(x=True_leaves_count, y=Nodules, color=Species.x)) + geom_point()

# Investigate AcS more detail on #Nodules and #True_leaves
leaf.count[which(leaf.count$Species.x=='AcS'), c('Treatment.x', 'True_leaves_count', 'Nodules', 'Comment')]
  # Based on the data, #200 failed to infect all AcS plant which are small, and smaller plants tend to get few nodules in Fix-

# Investigate other species on #Nodules and #True_leaves
leaf.count[which(leaf.count$Species.x=='Bic'), c('Treatment.x', 'True_leaves_count', 'Nodules', 'Comment')]
  # 200, almost all plants dies, 2 plant lives, one of them has a single nodule, not pink.
  # 131, 3 dies, nodules are not pink.
  # 2, 3 dies, one plant has nodule
  # 187, single plant, pink nodule
  # 186, many nodule, pink
  # 4, pink nodule


#leaf.count[which(leaf.count$Species.x=='Suc'), c('Treatment.x', 'True_leaves_count', 'Nodules', 'Comment')]



  

#############################
# Heat-plot showing summary #
#############################

## Total nodules
# aggregate data
agg.total_nod <- aggregate(Nodules~Species+Treatment, data=data[, c('Treatment', 'Species', 'Nodules')], sum, na.rm=T)
agg.total_nod <- merge(agg.total_nod, trait, by='Treatment')
ggplot(agg.total_nod, aes(x=Species, y=Trait)) + geom_tile(aes(fill=Nodules), colour='white') + 
  scale_fill_gradient(low='white', high='steelblue') +
  geom_text(aes(label=signif(Nodules, digits=3))) + ggtitle('Total Nodules') 

## Mean nodules
# aggregate data
agg.mean_nod <- aggregate(Nodules~Species+Treatment, data=data[which(data$Nodules>=1), c('Treatment', 'Species', 'Nodules')], mean, na.rm=T)
agg.mean_nod <- merge(agg.mean_nod, trait, by='Treatment')
ggplot(agg.mean_nod, aes(x=Species, y=Trait)) + geom_tile(aes(fill=Nodules),  colour='white') + 
  scale_fill_gradient(low='white', high='steelblue') +
  geom_text(aes(label=signif(Nodules, digits=3))) + ggtitle('Mean Nodule')


## Infectivity
# How many plants in a treatment formed nodule
agg.infectivity <- aggregate(Nodules~Species+Treatment, data=data[which(data$Comment!='dead'), c('Treatment', 'Species', 'Nodules')], FUN= function(x) mean(x != 0, na.rm=T))
agg.infectivity <- merge(agg.infectivity, trait, by='Treatment')
ggplot(agg.infectivity, aes(Species, Trait)) + geom_tile(aes(fill=Nodules), colour='white') + 
  scale_fill_gradient(low='white', high='steelblue') +
  geom_text(aes(label=signif(Nodules, digits=3))) + ggtitle('Infectivity') 


## RG
agg.RG <- aggregate(RG~Species+Treatment, data=data, sum, na.rm=T)
agg.RG <- merge(agg.RG, trait, by='Treatment')
agg.RG$RG <- log(agg.RG$RG/100)

ggplot(agg.RG, aes(x=Species, y=Trait)) + geom_tile(aes(fill=log(RG)), colour='white') + 
  scale_fill_gradient(low='white', high='steelblue') +
  geom_text(aes(label=signif(RG, digits=3))) + ggtitle('ShootRG')
 


#############################
# Detecting the of outliers #
#############################

# outliers.ID <- c(101, 34, 235, 6, 72, 131, 195, 118)
# 
# nodulated.plants.wo.outliers <- data[which(data$Nodules >= 1 & !(data$Plant_ID %in% outliers.ID)),]
# 
# ggplot(nodulated.plants.wo.outliers, aes(x=Treatment, y=ShootRG/100, fill=Species)) + 
#   geom_bar(stat='summary',position='dodge', fun.y='mean', na.rm=T)
# ggplot(nodulated.plants.wo.outliers, aes(x=Treatment, y=log(ShootRG/100), fill=Species)) + 
#   geom_bar(stat='summary',position='dodge', fun.y='mean', na.rm=T)
# ggplot(nodulated.plants.wo.outliers, aes(x=Treatment, y=RG/100, fill=Species)) + geom_boxplot()



################################################################
# Eliminating single occurance from treatment-host combination #
################################################################


# Nodulated plants
nodulated.plants <- data[which(data$Nodules >= 1), ]


# Aggregate nodulated plants nodule counts by host and treatment
nodulated.plants$Count <- 1
agg.nod.plants <- aggregate(Count~Species+Treatment, data=nodulated.plants, table)
agg.nod.plants[which(agg.nod.plants$Count <= 1), ]

# Getting the plant ID of single occurance plants
nodulated.plants[which(nodulated.plants$Species == 'AcW' & nodulated.plants$Treatment == 186), c('Plant_ID', 'Species', 'Treatment')]
nodulated.plants[which(nodulated.plants$Species == 'Bic' & nodulated.plants$Treatment == 187), c('Plant_ID', 'Species', 'Treatment')]
nodulated.plants[which(nodulated.plants$Species == 'Bic' & nodulated.plants$Treatment == 2), c('Plant_ID', 'Species', 'Treatment')]
nodulated.plants[which(nodulated.plants$Species == 'AcS' & nodulated.plants$Treatment == 200), c('Plant_ID', 'Species', 'Treatment')]
nodulated.plants[which(nodulated.plants$Species == 'Bic' & nodulated.plants$Treatment == 200), c('Plant_ID', 'Species', 'Treatment')]

single.occurance <- c(33, 94, 131, 162, 252)


# Dataset containing atleast two plants nodulated in a treatment-host combination
nodulated.plants.wo.single <- nodulated.plants[which(!(nodulated.plants$Plant_ID %in% single.occurance)), ]
aggregate(Count~Species+Treatment, data=nodulated.plants.wo.single, table)


###################################
# Final Visualization and Testing #
###################################

# Other possible parameters
# - Nodule quality (size)
# - Infection Rate (Number of nodule by root size)
# - Infection rate vs. RG
# - Nodule area to size ratio

nodulated.plants.wo.single$infection_rate <- nodulated.plants.wo.single$Nodules / nodulated.plants.wo.single$RootMass_mg
nodulated.plants.wo.single$mean_dry_nodule_biomass <- nodulated.plants.wo.single$dry_nodule_biomass / nodulated.plants.wo.single$Nodules
nodulated.plants.wo.single$RG_per_nodule <- nodulated.plants.wo.single$RG / nodulated.plants.wo.single$Nodules


p.nod <- symbio.plot(nodulated.plants.wo.single, 'Nodules', 'Treatment', 'Host')
p.infection_rate <- symbio.plot(nodulated.plants.wo.single, 'infection_rate', 'Treatment', 'Host')
grid.arrange(p.nod, p.infection_rate, nrow=2, top='#Nodules and #Nodules per mg of root')

p.rg <- symbio.plot(nodulated.plants.wo.single, 'RG', 'Treatment', 'Host')
p.rgd <- symbio.plot(nodulated.plants.wo.single, 'RGD', 'Treatment', 'Host')

p.shoot_rg <- symbio.plot(nodulated.plants.wo.single, 'ShootRG', 'Treatment', 'Host')
p.shoot_rgd <- symbio.plot(nodulated.plants.wo.single, 'ShootRGD', 'Treatment', 'Host')

grid.arrange(p.rg, p.shoot_rg, nrow=2, top='RG vs ShootRG on Nodulated Plants w/ >=2 occurance')
grid.arrange(p.shoot_rgd, nrow=1, top='ShootRGD on Nodulated Plants w/ >=2 Occurance')


p.investment <- symbio.plot(nodulated.plants.wo.single, 'investment', 'Treatment', 'Host')
p.nod_biomass <- symbio.plot(nodulated.plants.wo.single, 'dry_nodule_biomass', 'Treatment', 'Host')
p.mean_nod_biomass <- symbio.plot(nodulated.plants.wo.single, 'mean_dry_nodule_biomass', 'Treatment', 'Host')

grid.arrange(p.nod_biomass, p.investment, nrow=2, top='Dry Nodule Biomass and Investment into Nodulation')
grid.arrange(p.mean_nod_biomass, p.investment, nrow=2, top='Mean Dry Nodule Biomass and Investment into Nodulation')


p.efficiencty <- symbio.plot(nodulated.plants.wo.single, 'RG_per_nodule', 'Treatment', 'Host')

grid.arrange(p.investment, p.efficiencty, nrow=2, top='Investment into Nodulation & RG per Nodules')

p.efficiencty.wo.single <- ggplot(nodulated.plants.wo.single, aes(x=Treatment, y=Nodules, color=Host, size=dry_nodule_biomass, group=Host)) + 
  geom_point(alpha=0.7, position=position_dodge(width=0.9)) +
  scale_size(range=c(.1, 10), name='Nodule Biomass (mg)') +
  ggtitle('Total nodules and dry-biomass') +  
  facet_grid(~Treatment, space='free_x', scales='free_x')


grid.arrange(p.investment, p.efficiencty.wo.single, nrow=2, top='Investment into Nodulation & Nodule count-size distribution')


##########################
# Based on Lorena's code #
##########################

    ## Trait: Relative Growth Gain
    ## Formula: Regus 2015

# Main model describing all factors
lm.logRG.1 <- lm(log(ShootRG) ~ Treatment + dpi + Host, data=nodulated.plants.wo.single)
summary(lm.logRG.1)
Anova(lm.logRG.1, type='III')

  #Testing model assumptions:
  hist(resid(lm.logRG.1, type='deviance'))
  qqnorm(resid(lm.logRG.1, type='deviance'))
  qqline(resid(lm.logRG.1, type='deviance'))
  shapiro.test(resid(lm.logRG.1, type='deviance'))
  plot(lm.logRG.1)
  
  #Testing random factor
  lm.logRG.2 <- lm(log(ShootRG) ~ Treatment + dpi, data=nodulated.plants.wo.single)
  summary(lm.logRG.2)
  anova(lm.logRG.1, lm.logRG.2)
  
  #Post-hoc test
  #1.
  lsmean.logRG1.1 <- lsmeans(lm.logRG.1, adjust='tukey')
  
  


###########################
# Publication Ready Plots #
###########################
# Fix + = Green
# Fix - = Orange
# Non-nodulating = Blue
library(ggsci)
  
## Creating dataset that shows all treatment-host combination, but counts traits from the combination which has at least two plants nodulated
bic.187 <- data[which(data$Treatment == 187 & data$Host == 'Bic'), ]
bic.2 <- data[which(data$Treatment == 2 & data$Host == 'Bic'), ]
nan.2 <- data[which(data$Treatment == 2 & data$Host == 'Nan'), ]
acs.200 <- data[which(data$Treatment == 200 & data$Host == 'AcS'), ]
bic.200 <- data[which(data$Treatment == 200 & data$Host == 'Bic'), ]

no.stat <-  rbind(bic.187, bic.2, nan.2, acs.200, bic.200)  
no.stat$Nodules <- NA  
no.stat$investment <- NA
no.stat$dry_nodule_biomass <- NA
no.stat$ShootRGD <- NA
no.stat$Count <- NA
no.stat$infection_rate <- NA
no.stat$mean_dry_nodule_biomass <- NA
no.stat$RG_per_nodule <- NA

dim(no.stat)
dim(nodulated.plants.wo.single)

nodulated.plants.wo.single <- rbind(nodulated.plants.wo.single, no.stat)
nodulated.plants.wo.single <- nodulated.plants.wo.single[order(nodulated.plants.wo.single$Fix, decreasing=T),]
fix_order <- c('4', '131', '2', '156', '186', '187', '200')

## Plot total nodules and dry-bioass
p.nod <- ggplot(nodulated.plants.wo.single, aes(x=factor(Treatment, level=fix_order), y=Nodules, fill=Fix, size=dry_nodule_biomass, group=Host)) + 
  geom_point(alpha=0.9, position=position_dodge(width=0.9), shape=21) +
  scale_size(range=c(.1, 10), name='Nodule Biomass (mg)') +
  ggtitle('B. Nodules') +  
  facet_grid(~Host, space='free_y', scales='free_y') +
  scale_fill_manual(values=c('#FFA500', "#0d98ba"), labels=c('Ineffective', 'Effective')) + 
  theme(legend.title = element_blank(), axis.text.x=element_text(angle=0)) + xlab('') +
  coord_cartesian(ylim=c(0, 27))


## Plot investment
investment.data <- nodulated.plants.wo.single[, c('Treatment', 'Host', 'investment', 'Fix')]
agg.mean.in <- aggregate(investment~Treatment+Host+Fix, data=investment.data, mean, na.rm=T)
agg.se.in <- aggregate(investment~Treatment+Host+Fix, data=investment.data, se)
colnames(agg.se.in)[4] <- 'se'
agg.in <- merge(agg.mean.in, agg.se.in, by=c('Treatment', 'Host', 'Fix'))
#colnames(agg.in)[4] <- 'investment'

Host <- c('AcS', 'Bic', 'Bic', 'Bic', 'Nan')
Treatment <- c(200, 187, 2, 200, 2)
Fix <- c('minus', 'minus', 'minus', 'minus', 'minus')
investment <- c(NA, NA, NA, NA, NA)
se <- c(NA, NA, NA, NA, NA)
new.agg.in <- data.frame(Treatment, Host, Fix, investment, se)
agg.in <- rbind(agg.in, new.agg.in)

head(agg.in)

p.inv <- ggplot(agg.in, aes(x=factor(Treatment, level=fix_order), y=investment, fill=Fix)) +
  geom_bar(stat='summary', position='dodge', fun.y='mean', na.rm=T) +
  geom_errorbar(aes(ymin=investment-se, ymax=investment+se), position=position_dodge(width=0.90), width=0.2)+
  facet_grid(~Host, space='free_x', scales='free_x') +
  xlab('') + ylab('Investment') + ggtitle('C. Investment') +
  scale_fill_manual(values=c('#FFA500', "#0d98ba") , labels=c('Ineffective', 'Effective'))  + 
  theme(legend.title = element_blank(), axis.text.x=element_text(angle=0))
  

## Plot relative growth
RGD.data <- nodulated.plants.wo.single[, c('Treatment', 'Host', 'ShootRGD', 'Fix')]
agg.mean.rgd <- aggregate(ShootRGD~Treatment+Host+Fix, data=RGD.data, mean, na.rm=T)
agg.se.rgd <- aggregate(ShootRGD~Treatment+Host+Fix, data=RGD.data, se)
colnames(agg.se.rgd)[4] <- 'se'
agg.rgd <- merge(agg.mean.rgd, agg.se.rgd, by=c('Treatment', 'Host', 'Fix'))
head(agg.rgd)

Host <- c('AcS', 'Bic', 'Bic', 'Bic', 'Nan')
Treatment <- c(200, 187, 2, 200, 2)
Fix <- c('minus', 'minus', 'minus', 'minus', 'minus')
ShootRGD <- c(NA, NA, NA, NA, NA)
se <- c(NA, NA, NA, NA, NA)
new.agg.rgd <- data.frame(Treatment, Host, Fix, ShootRGD, se)
agg.rgd <- rbind(agg.rgd, new.agg.rgd)

p.rgb <- ggplot(agg.rgd, aes(x=Host, y=ShootRGD, fill=Fix)) +
  geom_bar(stat='summary', position='dodge', fun.y='mean', na.rm=T) +
  geom_errorbar(aes(ymin=ShootRGD-se, ymax=ShootRGD+se), position=position_dodge(width=0.90), width=0.2)+
  facet_grid(~factor(Treatment, level=fix_order), space='free_x', scales='free_x') +
  xlab('') + ylab('RGB') + ggtitle('A. Relative Growth Benefit') +
  scale_fill_manual(values=c('#FFA500', "#0d98ba"), labels=c('Ineffective', 'Effective')) +
  theme(legend.title = element_blank(), axis.text.x=element_text(angle=0))


grid.arrange(p.rgb, p.nod, p.inv, nrow=3, top='Host Generalism Experiment Outcome')


########################
# Statistical Testing  #
########################
library(ggfortify)
library(lme4)
library(car)
library(agricolae)

nodulated.plants.wo.single$Block_code <- 0 
nodulated.plants.wo.single[which(nodulated.plants.wo.single$Block == 'Large'), ]$Block_code <- 1

nod.2 <- nodulated.plants.wo.single[which(nodulated.plants.wo.single$Treatment == 2), ]
mod.nod2 <- lm(Nodules~Species, data=nod.2)
autoplot(mod.nod2)

nod.4 <- nodulated.plants.wo.single[which(nodulated.plants.wo.single$Treatment == 4), ]
mod.nod4 <- lm(Nodules~Species, data=nod.4)
autoplot(mod.nod2)

nod.186 <- nodulated.plants.wo.single[which(nodulated.plants.wo.single$Treatment == 186), ]
mod.nod186 <- lm(Nodules~Species, data=nod.186)
autoplot(mod.nod2)
Anova(nod.186)



nod.data <- nodulated.plants.wo.single
t.test(nod.data[which(nod.data$Host == 'AcS' & nod.data$Treatment == '156'), ]$ShootRGD, 
       data[which(data$Host == 'AcS' & data$Treatment =='H2O'), ]$ShootRGD, alternative='greater')

# Compared to 0
t.test(nod.data[which(nod.data$Host == 'AcS' & nod.data$Treatment == '131'), ]$ShootRGD, alternative='greater')


###########################################################
# Linear model for host investment in terms of nodulation #
###########################################################
lm.nod <- lm(Nodules ~ Host*Treatment + as.numeric(dpi), data=nodulated.plants.wo.single)
mod.diagnostics(lm.nod)

shapiro.test(residuals(lm.nod))

summary(aov(lm.nod))

plot(TukeyHSD(aov(lm.nod, type='II'), which='Host', conf.level=0.9))
(HSD.test(aov(lm.nod, type='II'), 'Host', group=T, alpha=0.9))


####################################
# Linear model for host investment #
####################################
lm.invest <- lm(log(investment) ~ Treatment*Host + dpi, data=nodulated.plants.wo.single)
mod.diagnostics(lm.invest)
shapiro.test(residuals(lm.invest))

summary(aov(lm.invest))
plot(TukeyHSD(aov(lm.invest), which='Host', conf.level=0.95))
(HSD.test(aov(lm.invest), 'Host', group=T, alpha=0.95))



####################################
# Linear model for RG #
####################################
lm.RG <- lm(log(mean_nodule_biomass) ~ Treatment+Host, data=nodulated.plants.wo.single)
mod.diagnostics(lm.RG)
shapiro.test(residuals(lm.RG))

summary(aov(lm.RG))
plot(TukeyHSD(aov(lm.RG), which='Host', conf.level=0.95))
(HSD.test(aov(lm.RG), 'Host', group=T, alpha=0.95))



#################
# Two way ANOVA #
#################
library(car)
library(agricolae)

# Modeling nodules 
mod.nodule <- lm(log10(Nodules)~Treatment*Species, data = nodulated.plants.wo.single)

mod.diagnostics(mod.nodule)
interaction.plot(nodulated.plants.wo.single, 'Nodules', 'Species', 'Treatment')

Anova(mod.nodule, type = 'II')

posthoc <- HSD.test(mod.nodule, c('Treatment', 'Species'), group=T)
posthoc

# Modeling RGD
mod.RG <- lm(log(ShootRG)~Treatment*Species + Nodules + total_nodule_area + dry_nodule_biomass + (1|Block_code), 
             data=nodulated.plants.wo.single)
mod.diagnostics(mod.RG)

interaction.plot(nodulated.plants.wo.single, 'ShootRG', 'Species', 'Treatment')

Anova(mod.RG, type='II')

posthoc <- HSD.test(mod.RG, c('Species', 'Treatment'))
posthoc

# Modeling Investment
mod.investment <- lm(log(investment)~Treatment*Species + Nodules + log(ShootRG) + total_nodule_area + dry_nodule_biomass + (1|Block_code), 
                     data=nodulated.plants.wo.single)
mod.diagnostics(mod.investment)
interaction.plot(nodulated.plants.wo.single, 'investment', 'Species', 'Treatment')
Anova(mod.investment, type='II')

posthoc <- HSD.test(mod.diagnostics, c('Species', 'Treatment'), group=T)
posthoc

# visualize nodulated plant dataset 
# Which response variable resolves the dataset better?
p1 <- symbio.plot(nodulated.plants, 'RG', 'Treatment', 'Host')
p2 <- symbio.plot(nodulated.plants, 'ShootRG', 'Treatment', 'Host')
grid.arrange(p1, p2, nrow=2, top='Comparing RG vs ShootRG')


##################################################################################################################
##################################################################################################################

#' #' Do simple ANOVA using aov
#' #' Present ANOVA summary
#' #' Do a TUKEY HSD posthoc test using given variable
#' #' Plot TUKEY HSD output
#' anova.report <- function(model, posthoc.var){
#'   res <-aov(model)
#'   summary(res)
#'   posthoc <- HSD.test(model, posthoc.var)
#'   posthoc
#'   plot(posthoc)
#' }

################
# By Host line #
################
signf.bar.plot <- function (dataset, response, factor, trait, signf=NULL){
  
  # Defining custom function
  se <- function(x) sd(x, na.rm=T)/sqrt(length(x))
  
  # Dataset aggregation
  temp.dataset <- dataset[, c(response, factor, trait)]
  colnames(temp.dataset) <- c('y', 'x', 't')
  agg.mean <-aggregate(y~x+t, data=temp.dataset, mean, na.rm=T)
  agg.se <- aggregate(y~x+t, data=temp.dataset, se)
  colnames(agg.se)[3] <- 'se'
  agg <- merge(agg.mean, agg.se, by=c('x', 't'))
  
  # Plotting Barplot with errorbar and significance
  if(is.null(signf) == T){
  p <- ggplot(data=agg, aes(x=x, y=y, fill=t)) +
    geom_bar(stat='identity') +
    geom_errorbar(aes( ymin=y-se, ymax=y+se), width=0.2)
  } else {
  p + geom_text(aes(y=y, label=signf, color=signf), vjust=-0.25)
  }
}

## Seperate dataset by host
nodulated.plants.wo.single$Block_code <- 1
nodulated.plants.wo.single[which(nodulated.plants.wo.single$Block == 'Large'), ]$Block_code <- 2

AcS <- nodulated.plants.wo.single[which(nodulated.plants.wo.single$Species == 'AcS'), ]
Bic <- nodulated.plants.wo.single[which(nodulated.plants.wo.single$Species == 'Bic'), ]
Suc <- nodulated.plants.wo.single[which(nodulated.plants.wo.single$Species == 'Suc'), ]
Nan <- nodulated.plants.wo.single[which(nodulated.plants.wo.single$Species == 'Nan'), ]



#############################################
## Host Benefit in each Hosts by Treatments #
#############################################

p <- interaction.plot(nodulated.plants.wo.single, 'RGD', 'Treatment', 'Host', log)
p

## Relative growth AcS
model.RGD <- lm( log(RG)~Treatment + dpi, data=AcS)
mod.diagnostics(model.RGD)
Anova(model.RGD, type='III')
res <- aov(model.RGD)
summary(res)
posthoc <- HSD.test(model.RGD, 'Treatment', group=T)
plot(posthoc)

AcS$ShootRGD <- log(AcS$ShootRGD + 1)
p <- signf.bar.plot(AcS, 'ShootRGD', 'Treatment', 'Fix', c('a', 'ab', 'ab', 'b','ab', 'a'))
p + ggtitle('AcS Host Response') + xlab('Treatments') + ylab('log(ShootRGD + 1)')

## Relative growth BIC
model.RGD <- lm(log(ShootRGD+1)~Treatment+log10(Nodules), data=Bic)
mod.diagnostics(model.RGD)

res <- aov(model.RGD)
summary(res)
posthoc <- HSD.test(model.RGD, 'Treatment', group=T)
posthoc$groups
plot(posthoc)

Bic$ShootRGD <- log(Bic$ShootRGD + 1)
p <- signf.bar.plot(Bic, 'ShootRGD', 'Treatment', 'Fix', c('ab', 'b', 'a', 'ab'))
p + ggtitle('Bic Host Response') + xlab('Treatments') + ylab('log(ShootRGD + 1)')

## Relative growth Suc
model.RGD <- lm(log(ShootRGD+1)~Treatment+log10(Nodules), data=Suc)
mod.diagnostics(model.RGD)

res <- aov(model.RGD)
summary(res)
posthoc <- HSD.test(model.RGD, 'Treatment')
plot(posthoc)

Suc$ShootRGD <- log(Suc$ShootRGD + 1)
p <- signf.bar.plot(Suc, 'ShootRGD', 'Treatment', 'Fix')
p + ggtitle('Suc Host Response') + xlab('Treatments') + ylab('log(ShootRGD + 1)')

## Relative growth Nan
model.RGD <- lm(log(ShootRGD+1)~Treatment+log10(Nodules), data=Nan)
mod.diagnostics(model.RGD)

res <- aov(model.RGD)
summary(res)
posthoc <- HSD.test(model.RGD, 'Treatment')
plot(posthoc)

Nan$ShootRGD <- log(Nan$ShootRGD + 1)
p <- signf.bar.plot(Nan, 'ShootRGD', 'Treatment', 'Fix')
p + ggtitle('Nan Host Response') + xlab('Treatments') + ylab('log(ShootRGD + 1)')


# t-test for % relative growth grain 
t.test(Suc[which(Suc$Treatment == '200'), ]$ShootMass_mg, 
       data[which(data$Host =='Suc' & data$Treatment == 'H2O'), ]$ShootMass_mg, alternative='greater')

# mixed linear model
library(lme4)
library(car)
library(stargazer)

lmm.rgd <- lmer(log10(ShootRGD) ~ Treatment + (1|Host), data=nodulated.plants.wo.single)
summary(lmm.rgd)
plot(lmm.rgd)
qqnorm(resid(lmm.rgd))
qqline(resid(lmm.rgd))

mm_plot <- ggplot(nodulated.plants.wo.single, aes(x=Treatment, y=log10(ShootRGD))) +
  facet_wrap(~Host, nrow=2) +
  geom_point(alpha=0.5) +
  geom_line(data=cbind(nodulated.plants.wo.single, pred=predict(lmm.rgd)), aes(y=pred), size=1 )

Anova(lmm.rgd, type='III')

stargazer(lmm.rgd, type='text', digits=3, star.cutoffs=c(0.05, 0.01, 0.001), digit.separator = "")

#########################################
# Nodulation in each hosts by treatment #
#########################################

interaction.plot(nodulated.plants.wo.single, 'Nodules', 'Treatment', 'Host')

## Nodulation in AcS
model.Nod <- lm(log10(Nodules)~Treatment + dpi, data=AcS)
mod.diagnostics(model.Nod)

res <- aov(model.Nod)
summary(res)
posthoc <- HSD.test(model.Nod, 'Treatment')
plot(posthoc)

AcS$Nodules <- log10(AcS$Nodules)
p <- signf.bar.plot(AcS, 'Nodules', 'Treatment', 'Fix')
p + ggtitle('AcS Host Nodulation') + xlab('Treatments') + ylab('log10(Nodules)')

## Nodulation in Bic
model.Nod <- lm(log10(Nodules)~Treatment+log(RootMass_mg), data=Bic)
mod.diagnostics(model.Nod)

res <- aov(model.Nod)
summary(res)
posthoc <- HSD.test(model.Nod, 'Treatment')
plot(posthoc)
posthoc$groups


Bic$Nodules <- log10(Bic$Nodules)
p <- signf.bar.plot(Bic, 'Nodules', 'Treatment', 'Fix',  c('a','a', 'b', 'b'))
p + ggtitle('Bic Host Nodulation') + xlab('Treatments') + ylab('log10(Nodules)')

## Nodulation in Suc
model.Nod <- lm(log10(Nodules)~Treatment+log(RootMass_mg), data=Suc)
mod.diagnostics(model.Nod)

res <- aov(model.Nod)
summary(res)
posthoc <- HSD.test(model.Nod, 'Treatment')
plot(posthoc)

Suc$Nodules <- log10(Suc$Nodules)
p <- signf.bar.plot(Suc, 'Nodules', 'Treatment', 'Fix')
p + ggtitle('Suc Host Nodulation') + xlab('Treatments') + ylab('log10(Nodules)')

## Nodulation in Nan
model.Nod <- lm(log10(Nodules)~Treatment+log(RootMass_mg), data=Nan)
mod.diagnostics(model.Nod)

res <- aov(model.Nod)
summary(res)
posthoc <- HSD.test(model.Nod, 'Treatment')
plot(posthoc)

Nan$Nodules <- log10(Nan$Nodules)
p <- signf.bar.plot(Nan, 'Nodules', 'Treatment', 'Fix')
p + ggtitle('Suc Host Nodulation') + xlab('Treatments') + ylab('log10(Nodules)')


#########################################
# Investment in each hosts by treatment #
#########################################
interaction.plot(nodulated.plants.wo.single, 'investment', 'Treatment', 'Host')

## Investment in AcS
model.Inv <- lm(investment~Treatment+log10(Nodules)+log(ShootRGD)+log(dry_nodule_biomass), data=AcS)
mod.diagnostics(model.Inv)

Anova(model.Inv)

res <- aov(model.Inv)
summary(res)
posthoc <- HSD.test(model.Inv, 'Treatment')
plot(posthoc)
posthoc$groups

p <- signf.bar.plot(AcS, 'investment', 'Treatment', 'Fix', c('ab', 'a', 'ab', 'b', 'ab', 'b'))
p + ggtitle('AcS Host Investment') + xlab('Treatments') + ylab('Investment')

## Investment in Bic
model.Inv <- lm(investment~Treatment+log10(Nodules)+log(ShootRGD)+log(dry_nodule_biomass), data=Bic)
mod.diagnostics(model.Inv)

res <- aov(model.Inv)
summary(res)
posthoc <- HSD.test(model.Inv, 'Treatment')
plot(posthoc)
posthoc$groups

p <- signf.bar.plot(Bic, 'investment', 'Treatment', 'Fix', c('bc', 'c', 'a', 'b'))
p + ggtitle('Bic Host Investment') + xlab('Treatments') + ylab('Investment')

## Investment in Suc
model.Inv <- lm(investment~Treatment+log10(Nodules)+log(ShootRGD)+log(dry_nodule_biomass), data=Suc)
mod.diagnostics(model.Inv)

res <- aov(model.Inv)
summary(res)
posthoc <- HSD.test(model.Inv, 'Treatment')
plot(posthoc)

p <- signf.bar.plot(Suc, 'investment', 'Treatment', 'Fix')
p + ggtitle('Suc Host Investment') + xlab('Treatments') + ylab('Investment')

## Investment in Bic
model.Inv <- lm(investment~Treatment+log10(Nodules)+log(ShootRGD)+log(dry_nodule_biomass), data=Nan)
mod.diagnostics(model.Nod)

res <- aov(model.Inv)
summary(res)
posthoc <- HSD.test(model.Inv, 'Treatment')
plot(posthoc)

p <- signf.bar.plot(Bic, 'investment', 'Treatment', 'Fix')
p + ggtitle('Bic Host Investment') + xlab('Treatments') + ylab('Investment')
###############################################################################################
###############################################################################################

# Nod+ plants (filter out the marginally nodulating plants)

# This does filter out unnodulated plant
single_occurance <- c(237, 131, 94, 33) #plant id
#nod.plus.plants <-data[which(!data$Plant_ID %in% single_occurance), ] 

nod.plus.plants <- data[which(data$Nodules >= 1 & !data$Plant_ID %in% single_occurance), ]


nod.part <- aggregate(Nodules~Species+Treatment, data=nod.plus.plants, mean)
nod_se.part <- aggregate(Nodules~Species+Treatment, data=nod.plus.plants, se)
colnames(nod_se.part)[3] <- 'Nodules_se'
nod <- merge(nod.part, nod_se.part, by=c('Species', 'Treatment'))

rg.part <- aggregate(RG~Species+Treatment, data=nod.plus.plants, mean, na.rm=T)
rg_se.part <- aggregate(RG~Species+Treatment, data=nod.plus.plants, se)
colnames(rg_se.part)[3] <- 'RG_se'
rg <- merge(rg.part, rg_se.part, by=c('Species', 'Treatment'))

agg <- merge(nod, rg, by=c('Species', 'Treatment'))            

ggplot(agg[which(agg$Treatment != 'H2O'),], aes(x=Treatment, y=RG, fill=Species)) + 
  geom_bar(stat='identity',position='dodge', na.rm=T) + 
  geom_point(data=nod.plus.plants[which(nod.plus.plants$Treatment != 'H2O'),], aes(x=Treatment, color=Species), colour='black', pch=21, position=position_dodge(width=0.9)) +
  geom_errorbar(aes( ymin=RG-RG_se, ymax=RG+RG_se), position = position_dodge(width = 0.90), width=0.2) +
  ggtitle('Relative Growth of Nod+ Treatments (n >= 2)') + geom_hline(yintercept=1, linetype='dashed') +
  facet_grid(~Treatment, space='free_x', scales='free_x')

ggplot(agg[which(agg$Treatment != 'H2O'),], aes(x=Treatment, y=Nodules, fill=Species)) + 
  geom_bar(stat='summary',position='dodge', fun.y='mean', na.rm=T) + 
  geom_errorbar(aes( ymin=Nodules-Nodules_se, ymax=Nodules+Nodules_se), position = position_dodge(width = 0.90), width=0.2) +
  ggtitle('Median Nodules by Nod+ Treatments') +
  facet_grid(~Treatment, space='free_x', scales='free_x')



symbio.plot(data, 'RG', 'Treatment')

