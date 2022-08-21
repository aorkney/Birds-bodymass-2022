# The purpose of this script is to compute the scheme of allometric dependencies of skeletal element size
# upon body mass across birds, and to discern whether the organisation is modular. 

library( geomorph )
library( ape )
library( nlme )
# Load required analytical R packages

setwd('D:/Documents/Alex_birds/Ideas_2022/SICBY')
# Set work directory as appropriate

load('Csize.29.07.2022.RData')
# Load centroid sizes for study birds. 
# In our study, we defined this object as an array called 'GPA.Csize'
# Be aware the names do not yet match the names in the phylogeny we will use.

load('tree.29.07.2022.RData')
# Load phylogenetic tree of birds
# (Pruned from Prum et al., 2015)
# The object is a phylogeny called 'pruned.tree'.

load('tree.names.29.07.2022.RData')
# Load the tree names required to match birds to the closest genus on the tree.
# The object is a vector called 'tree_names'.

load('masses.29.07.2022.Rdata')
# Load bird masses
# The object is a named numeric vector called 'masses'.
# The names will not yet match those on the phylogeny.

# Let's now ensure names match the phylogeny. 

names(masses) <- tree_names

for( i in 1:length(GPA.Csize) ){
	names(GPA.Csize[[i]]) <- tree_names 
}

# For 30 random subsamples of 140 taxa,  
# compute the allometric dependency of log10(centroid size) upon log10(body mass), 
# across all skeletal elements. 
# This will be undertaken in a phylogenetic context, using generalised least squares regression. 
# Phylogenetic variance-covariance is represented by the corPagel correlation structure of 
# the Prum et al., 2015 phylogeny. 

ordered_bones<-c('skull','mandible','carpometacarpus','radius','ulna','humerus',
'sternum','coracoid','scapula','synsacrum','femur','tibiotarsus','tarsometatarsus')
# Bone names, arranged by anatomical module.

dataframe<-matrix(NA,30,length(ordered_bones))
colnames(dataframe)<-ordered_bones
# Prepare a data frame.

for(i in 1:length(ordered_bones)){ # For each bone
	for(j in 1:30){ # For 30 replicates
		taxa<-sample(names(masses),140) # Produce a subsample of 140 taxa
		newphy <- drop.tip( pruned.tree , pruned.tree$tip.label[ !pruned.tree$tip.label %in% taxa ] ) # Prune the phylogeny
		df <- data.frame( mass= log10(masses[taxa]), Csize = log10( GPA.Csize[[ ordered_bones[i] ]][taxa] ) ) # Arrange a dataframe
 		fit<- gls( Csize ~ mass, correlation = corPagel( 0, phy=newphy ), data=df ) # Perform a GLS fit
		dataframe[j,i] <- fit$coef[[2]] # Extract the coefficient of the relationship
	}
}

library(ggplot2)
library(reshape2)
# These packages are required to plot and organise the output.

molten.df<-melt(dataframe)
molten.df$Var1[which(molten.df$Var2=='skull' | molten.df$Var2=='mandible')]<-'head'
molten.df$Var1[which(molten.df$Var2=='coracoid' | molten.df$Var2=='scapula' | molten.df$Var2=='sternum'| molten.df$Var2=='synsacrum')]<-'trunk'
molten.df$Var1[which(molten.df$Var2=='humerus' | molten.df$Var2=='radius' | molten.df$Var2=='ulna'| molten.df$Var2=='carpometacarpus')]<-'wing'
molten.df$Var1[which(molten.df$Var2=='femur' | molten.df$Var2=='tibiotarsus'| molten.df$Var2=='tarsometatarsus')]<-'leg'
# Organise the output

dev.new(height=12,width=16.9,unit='cm')
# Prepare a new plot

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# This is a colour blind friendly palette. 


ggplot(data=molten.df)+
geom_violin( aes(x=Var2,y=value,group=Var2, col=Var1),lwd=2)+
labs(y='coefficient',col='module')+
scale_colour_manual( values=c(rep(cbbPalette[7],1),rep(cbbPalette[8],1),rep(cbbPalette[6],1),rep(cbbPalette[4],1)) )+
  theme_bw() + theme(axis.text.x=element_text(size=15, angle=90, vjust=0.3,colour=
	c(rep(cbbPalette[7],2),rep(cbbPalette[4],4),rep(cbbPalette[6],4),rep(cbbPalette[8],3)) ,face = "bold"),
	axis.title.x=element_blank(),
	axis.title.y=element_text(size=15,face='bold'),
	axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.y=element_text(size=15,colour='black'),
      axis.ticks=element_line(size=2),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())

# Produce a violin plot to visualise the ouput. 


setwd('D:/Documents/Alex_birds/Ideas_2022/SICBY')
# ggsave(filename='allometric_scaling_20_08_2022.pdf')
# Save the output.



