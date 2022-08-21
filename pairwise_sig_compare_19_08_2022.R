# The purpose of this script is to load 
# centroid size information for 13 skeletal elements across a sample
# of 149 birds from across crown Aves. 
# The landmark constellations from which these data were extracted 
# are available at: https://doi.org/10.18563/journal.m3.125
# (Bjarnason & Benson, 2021)
# Thereafter, the birds' masses will be used to perform
# a phylogenetically gnostic Generalised Least Squares regression, 
# and the residuals will be extracted.
# This will remove allometric scaling in centroid size. 
# The phylogeny of Prum et al., 2015 was used:
# https://doi.org/10.1038/nature15697
# Residual centroid sizes of skeletal elements across birds are known 
# to exhibit a modular organisation of evolutionary covariances. 
# https://doi.org/10.1038/s41559-021-01509-w
# (Orkney, Bjarnason, et al., 2021)

# Having prepared the dataset, 
# phylogenetic two block Partial Least Squares analyses will be performed,
# in order to assess the degree of integration between different pairwise 
# combinations of bones (carpometacarpus-humerus, sternum-scapula,
# sternum-carpometacarpus). These integration statistics will be computed
# in different cohorts of 30 birds of varying mass distributions.
# The significance of difference between the integration statistics will 
# then be computed between the cohorts, using the compare.pls function
# from the R Geomorph package. 
# This will clarify whether integration within the wing, within the trunk
# and between the wing and trunk, varies significantly as a function of body mass,
# across birds. 


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


# We may wish to investigate a specific taxonomic subset of 
# the available birds. 

setwd("D:/Documents/Alex_birds/")
# Set work directory as appropriate

Telluraves<-read.csv('Telluraves.csv')
# This is a pre-prepared list of birds belonging to Telluraves.
# This will not actually be used in analysis, but may be helpful
# for other researchers who wish to perform custom analyses. 

Apodiformes<- c('Archilochus_colubris',
'Topaza_pella',
'Chaetura_brachyura',
'Streptoprocne_zonaris',
'Hemiprocne_comata')
# This is a pre-prepared list of birds belonging to Apodiformes.

setwd('D:/Documents/Alex_birds/Ideas_2022/SICBY')
# Set work directory as appropriate

get.item <- function( X , item ) { X[[ item ]] }
# This is a simple indexing function


get.residual.Csize <- function( array, masses, phylogeny, taxa ){
	allometry.Csize <- list()
	newphy <- drop.tip( phylogeny , phylogeny$tip.label[ !phylogeny$tip.label %in% taxa ] )
	for(i in 1:length(array) ){
		df <- data.frame( mass= log10(masses[taxa]), Csize = log10( array[[ i ]][taxa] ) )
		allometry.Csize[[ i ]] <- gls( Csize ~ mass, correlation = corPagel( 0, phy=newphy ), data=df )
	}
	names(allometry.Csize) <- names(array)
	residual_Csize <- lapply( lapply( allometry.Csize , get.item , item = "residuals" ) , c )
	return(residual_Csize)
}

# This function uses phylogenetic Generalised Least Squares to compute the allometrically adjusted centroid sizes of skeletal elements.
# The user must provide the array of original bird skeletal element centroid sizes, bodymasses, phylogeny and taxa of interest. 

residual.Csize <- get.residual.Csize( array = GPA.Csize, masses = masses, phylogeny = pruned.tree, taxa=names(masses) )
# Which taxa do we wish to investigate? 
# In this example, we will investigate all birds. 

residual.Csize.nApods <- get.residual.Csize( array = GPA.Csize, masses = masses, phylogeny = pruned.tree, taxa=names(masses)[ !names(masses) %in% Apodiformes] )
residual.Csize.T <- get.residual.Csize( array = GPA.Csize, masses = masses, phylogeny = pruned.tree, taxa= names(masses)[ names(masses) %in% as.character(Telluraves$x)]  )
# Residual Centroid sizes, computed within pre-defined taxonomic subsets. 

masses.all<-masses
masses.nApods<-masses[names(masses)[ !names(masses) %in% Apodiformes]]
masses.T<-masses[ names(masses) %in% as.character(Telluraves$x)]
# Mass distributions of the pre-defined taxonomic subsets. 


width <- 30
# The number of birds to be included in any cohort. 

# The following is a function that will compute significance of difference of pairwise integration
# between two bones among various combinations of lighter and heavier cohorts of birds, 
# within a specific taxonomic subset. 

pair.sig <- function( mass, proportions, phylogeny, bone1, bone2){
	lights <- c(1:width) # Begin by defining the lightest cohort of birds within the birds of interest. 
	bin.mass<-matrix(NA,(length( mass )-2*(length(lights))),(length( mass )-2*(length(lights))))  
	p.dff<-matrix(NA,(length( mass )-2*(length(lights))),(length( mass )-2*(length(lights))))
	# Define matrices to receive data. 

	for(k in 1: dim(bin.mass)[1]){ # For all possible combinations of cohorts of birds. 
		lighter<-which( mass[order(mass)] < mass[order(mass)][min(lights)] ) # Find all birds less massive than the 'light' cohort.
		heavies<-c(1:length(mass))[-c(lighter,lights)] # Define a cohort of birds heavier than the 'light' cohort. 
		for(i in 1:(length(heavies)-width) ){ # For all possible cohorts of birds within the heavier group. 
			start <- heavies[i] 
			end <- heavies[i+(width-1)]
			bin.mass[k,i+length(lighter)] <- mean(mass[order(mass)][start:end]) # Compute the mean mass of the heavy cohort. 
			taxa.heavies <- names(mass[order(mass)][start:end]) # Name the birds in the heavy cohort. 
			newphy.heavies <- drop.tip( phylogeny, phylogeny$tip.label[ !phylogeny$tip.label  %in% taxa.heavies] ) # Prune the phylogeny. 
			wing.int.heavies <- phylo.integration(as.matrix(proportions[[bone1]][taxa.heavies]),as.matrix(proportions[[bone2]][taxa.heavies]),phy=newphy.heavies,print.progress=F)
			# Compute integration within the heavier cohort. 
			taxa.lights <- names(mass[order(mass)][lights ]) # Name the birds in the light cohort. 
			newphy.lights  <- drop.tip( phylogeny, phylogeny$tip.label[ !phylogeny$tip.label  %in% taxa.lights ] ) # Prune the phylogeny. 
			wing.int.lights  <- phylo.integration(as.matrix(proportions[[bone1]][taxa.lights ]),as.matrix(proportions[[bone2]][taxa.lights ]),phy=newphy.lights,print.progress=F )
			# Compute integration within the light cohort. 
			p.dff[k,i+length(lighter)] <- compare.pls(wing.int.lights,wing.int.heavies)$pairwise.P[2]
			# Compute the significance of difference between the two integration analyses. 
		}
	print(round(100*(k/dim(bin.mass)[1]))) # Update on progress. 
	lights<-lights+1 # Shift the light cohort to exclude the lightest bird within, and include the lightest bird of the remaining heavier birds. 
	}

	lights<-c(1:width)
	light.masses<-list()
	for(k in 1: dim(bin.mass)[1]){
		light.masses[[k]]<-mean(mass[order(mass)][lights])
		lights<-lights+1
	}
	# Compute the masses of the light cohorts. 

	colnames(p.dff)<-bin.mass[1,]
	rownames(p.dff)<-unlist(light.masses)
	# Tidy the data output. 

	output <- list(p.dff,bin.mass)
	return(output)
	# Produce output.
}

wing.all <- pair.sig(masses, residual.Csize, pruned.tree, 'humerus', 'carpometacarpus')
wing.nApods <- pair.sig(masses.nApods, residual.Csize.nApods, pruned.tree, 'humerus', 'carpometacarpus')
# Significance of integration change across body mass, for all birds and for birds excluding Apodiformes. 

library(reshape2)
library(MBA)
# These packages are necessary to reshape and interpolate data. 

df<-melt(wing.all[[1]])
df<-df[complete.cases(df),]
orig.df<-df
orig.df<-log10(orig.df)
colnames(orig.df)<-c('x','y','z')
grid_mba<-mba.surf(orig.df, no.X = 100, no.Y = 100, extend = T)
dimnames(grid_mba$xyz.est$z) <- list(grid_mba$xyz.est$x, grid_mba$xyz.est$y)
grid_mba <- melt(grid_mba$xyz.est$z, varnames = c('x', 'y'), value.name = 'z')
grid_mba[grid_mba$x>grid_mba$y,3]<-NA
# Interpolate the significance statistics onto a dense grid. 

df2<-melt(wing.nApods[[1]])
df2<-df2[complete.cases(df2),]
orig.df2<-df2
orig.df2<-log10(orig.df2)
colnames(orig.df2)<-c('x','y','z')
grid_mba2<-mba.surf(orig.df2, no.X = 100, no.Y = 100, extend = T)
dimnames(grid_mba2$xyz.est$z) <- list(grid_mba2$xyz.est$x, grid_mba2$xyz.est$y)
grid_mba2 <- melt(grid_mba2$xyz.est$z, varnames = c('x', 'y'), value.name = 'z')
grid_mba2[grid_mba2$x>grid_mba2$y,3]<-NA

wing.lims.x<-c(min(c(orig.df$y,orig.df2$y)),max(c(orig.df$y,orig.df2$y)))
wing.lims.y<-c(min(c(orig.df$x,orig.df2$x)),max(c(orig.df$x,orig.df2$x)))
wing.lims.p<-c(min(c(grid_mba$z,grid_mba2$z),na.rm=T),max(c(grid_mba$z,grid_mba2$z),na.rm=T))
# Define plotting limits. 

mp.wing.all <- ggplot(grid_mba, aes(y, x, fill = z, z=z)) +
lims(x=wing.lims.x,y=wing.lims.y)+
geom_raster( )+
geom_contour(colour='black',breaks=c(log10(0.05)),lwd=2) +
labs(y='mass light cohort', x='mass heavy cohort',fill=expression(paste(log[10],'(p)')))+
scale_fill_gradient2(high='red',low='blue',mid='white',na.value='white',midpoint=log10(0.05), limits=wing.lims.p)+
  theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=10,colour='black'),
      axis.text.y=element_text(size=10,colour='black'),
	axis.title.x.bottom=element_text(size=15,colour="black"),
	axis.title.y.left=element_text(size=15,colour="black"),
      axis.ticks=element_line(size=2),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())
# Produce a plot of the field of interpolated significance statistics across all combinations of lighter and heavier cohorts of birds, 
# illustrating significant shifts of integration within the wing across the study birds. 

# From hereon the analyses are repeated for different taxonomic subsets and combinations of bones; comments will be desisted if they are repetitive. 

mp.wing.nApods <- ggplot(grid_mba2, aes(y, x, fill = z, z=z)) +
lims(x=wing.lims.x,y=wing.lims.y)+
geom_raster( )+
geom_contour(colour='black',breaks=c(log10(0.05)),lwd=2) +
labs(y='mass light cohort', x='mass heavy cohort',fill=expression(paste(log[10],'(p)')))+
scale_fill_gradient2(high='red',low='blue',mid='white',na.value='white',midpoint=log10(0.05), limits=wing.lims.p)+
  theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=10,colour='black'),
      axis.text.y=element_text(size=10,colour='black'),
	axis.title.x.bottom=element_text(size=15,colour="black"),
	axis.title.y.left=element_text(size=15,colour="black"),
      axis.ticks=element_line(size=2),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())



# It is probably best not to try to assess patterns within Telluraves 
# The mass distribution range is simply not large enough.


trunk.all <- pair.sig(masses, residual.Csize, pruned.tree, 'sternum', 'scapula')
trunk.nApods <- pair.sig(masses.nApods, residual.Csize.nApods, pruned.tree, 'sternum', 'scapula')

df3<-melt(trunk.all[[1]])
df3<-df3[complete.cases(df3),]
orig.df3<-df3
orig.df3<-log10(orig.df3)
colnames(orig.df3)<-c('x','y','z')
grid_mba3<-mba.surf(orig.df3, no.X = 100, no.Y = 100, extend = T)
dimnames(grid_mba3$xyz.est$z) <- list(grid_mba3$xyz.est$x, grid_mba3$xyz.est$y)
grid_mba3 <- melt(grid_mba3$xyz.est$z, varnames = c('x', 'y'), value.name = 'z')
grid_mba3[grid_mba3$x>grid_mba3$y,3]<-NA

df4<-melt(trunk.nApods[[1]])
df4<-df4[complete.cases(df4),]
orig.df4<-df4
orig.df4<-log10(orig.df4)
colnames(orig.df4)<-c('x','y','z')
grid_mba4<-mba.surf(orig.df4, no.X = 100, no.Y = 100, extend = T)
dimnames(grid_mba4$xyz.est$z) <- list(grid_mba4$xyz.est$x, grid_mba4$xyz.est$y)
grid_mba4 <- melt(grid_mba4$xyz.est$z, varnames = c('x', 'y'), value.name = 'z')
grid_mba4[grid_mba4$x>grid_mba4$y,3]<-NA

trunk.lims.x<-c(min(c(orig.df3$y,orig.df4$y)),max(c(orig.df3$y,orig.df4$y)))
trunk.lims.y<-c(min(c(orig.df3$x,orig.df4$x)),max(c(orig.df3$x,orig.df4$x)))
trunk.lims.p<-c(min(c(grid_mba3$z,grid_mba4$z),na.rm=T),max(c(grid_mba3$z,grid_mba4$z),na.rm=T))

mp.trunk.all <- ggplot(grid_mba3, aes(y, x, fill = z, z=z)) +
lims(x=trunk.lims.x,y=trunk.lims.y)+
geom_raster( )+
geom_contour(colour='black',breaks=c(log10(0.05)),lwd=2) +
labs(y='mass light cohort', x='mass heavy cohort',fill=expression(paste(log[10],'(p)')))+
scale_fill_gradient2(high='red',low='blue',mid='white',na.value='white',midpoint=log10(0.05), limits=trunk.lims.p)+
  theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=10,colour='black'),
      axis.text.y=element_text(size=10,colour='black'),
	axis.title.x.bottom=element_text(size=15,colour="black"),
	axis.title.y.left=element_text(size=15,colour="black"),
      axis.ticks=element_line(size=2),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())

mp.trunk.nApods <- ggplot(grid_mba4, aes(y, x, fill = z, z=z)) +
lims(x=trunk.lims.x,y=trunk.lims.y)+
geom_raster( )+
geom_contour(colour='black',breaks=c(log10(0.05)),lwd=2) +
labs(y='mass light cohort', x='mass heavy cohort',fill=expression(paste(log[10],'(p)')))+
scale_fill_gradient2(high='red',low='blue',mid='white',na.value='white',midpoint=log10(0.05), limits=trunk.lims.p)+
  theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=10,colour='black'),
      axis.text.y=element_text(size=10,colour='black'),
	axis.title.x.bottom=element_text(size=15,colour="black"),
	axis.title.y.left=element_text(size=15,colour="black"),
      axis.ticks=element_line(size=2),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())


cross.all <- pair.sig(masses, residual.Csize, pruned.tree, 'sternum', 'carpometacarpus')
cross.nApods <- pair.sig(masses.nApods, residual.Csize.nApods, pruned.tree, 'sternum', 'carpometacarpus')

df5<-melt(cross.all[[1]])
df5<-df5[complete.cases(df5),]
orig.df5<-df5
orig.df5<-log10(orig.df5)
colnames(orig.df5)<-c('x','y','z')
grid_mba5<-mba.surf(orig.df5, no.X = 100, no.Y = 100, extend = T)
dimnames(grid_mba5$xyz.est$z) <- list(grid_mba5$xyz.est$x, grid_mba5$xyz.est$y)
grid_mba5 <- melt(grid_mba5$xyz.est$z, varnames = c('x', 'y'), value.name = 'z')
grid_mba5[grid_mba5$x>grid_mba5$y,3]<-NA

df6<-melt(cross.nApods[[1]])
df6<-df6[complete.cases(df6),]
orig.df6<-df6
orig.df6<-log10(orig.df6)
colnames(orig.df6)<-c('x','y','z')
grid_mba6<-mba.surf(orig.df6, no.X = 100, no.Y = 100, extend = T)
dimnames(grid_mba6$xyz.est$z) <- list(grid_mba6$xyz.est$x, grid_mba6$xyz.est$y)
grid_mba6 <- melt(grid_mba6$xyz.est$z, varnames = c('x', 'y'), value.name = 'z')
grid_mba6[grid_mba6$x>grid_mba6$y,3]<-NA

cross.lims.x<-c(min(c(orig.df5$y,orig.df6$y)),max(c(orig.df5$y,orig.df6$y)))
cross.lims.y<-c(min(c(orig.df5$x,orig.df6$x)),max(c(orig.df5$x,orig.df6$x)))
cross.lims.p<-c(min(c(grid_mba5$z,grid_mba6$z),na.rm=T),max(c(grid_mba5$z,grid_mba6$z),na.rm=T))


mp.cross.all <- ggplot(grid_mba5, aes(y, x, fill = z, z=z)) +
lims(x=cross.lims.x,y=cross.lims.y)+
geom_raster( )+
geom_contour(colour='black',breaks=c(log10(0.05)),lwd=2) +
labs(y='mass light cohort', x='mass heavy cohort',fill=expression(paste(log[10],'(p)')))+
scale_fill_gradient2(high='red',low='blue',mid='white',na.value='white',midpoint=log10(0.05), limits=cross.lims.p)+
  theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=10,colour='black'),
      axis.text.y=element_text(size=10,colour='black'),
	axis.title.x.bottom=element_text(size=15,colour="black"),
	axis.title.y.left=element_text(size=15,colour="black"),
      axis.ticks=element_line(size=2),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())

mp.cross.nApods <- ggplot(grid_mba6, aes(y, x, fill = z, z=z)) +
lims(x=cross.lims.x,y=cross.lims.y)+
geom_raster( )+
geom_contour(colour='black',breaks=c(log10(0.05)),lwd=2) +
labs(y='mass light cohort', x='mass heavy cohort',fill=expression(paste(log[10],'(p)')))+
scale_fill_gradient2(high='red',low='blue',mid='white',na.value='white',midpoint=log10(0.05), limits=cross.lims.p)+
  theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=10,colour='black'),
      axis.text.y=element_text(size=10,colour='black'),
	axis.title.x.bottom=element_text(size=15,colour="black"),
	axis.title.y.left=element_text(size=15,colour="black"),
      axis.ticks=element_line(size=2),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())

dev.new(width=16.9,height=16.9,unit='cm')
# Prepare a new plot.
library(ggpubr)
# Load package for arranging plots. 
ggarrange(mp.wing.all,mp.wing.nApods,mp.trunk.all,mp.trunk.nApods,mp.cross.all,mp.cross.nApods,ncol=2,nrow=3,align='hv',
labels=c('a','b','c','d','e','f'),label.x=-0.025,label.y=1.035,font.label=list(size=30))
# Arrange subplots. 

getwd() # Check working directory
# ggsave(filename='pairwise_sig_compare_plot_19_08_2022b.pdf')
# save the output. 

write.csv(x=df,file='wing_all_19_08_2022.csv')
write.csv(x=df2,file='wing_nApods_19_08_2022.csv')
write.csv(x=df3,file='trunk_all_19_08_2022.csv')
write.csv(x=df4,file='trunk_nApods_19_08_2022.csv')
write.csv(x=df5,file='cross_all_19_08_2022.csv')
write.csv(x=df6,file='cross_nApods_19_08_2022.csv')
# The above analyses took a long time to run; the raw data output can be saved so that they can be plotted again 
# later without having to re-run the analysis from scratch. 


