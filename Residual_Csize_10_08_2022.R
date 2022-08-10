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
# The residual centroid sizes may then be investigated by other scripts. 
# The phylogeny of Prum et al., 2015 was used:
# https://doi.org/10.1038/nature15697
# Residual centroid sizes of skeletal elements across birds are known 
# to exhibit a modular organisation of evolutionary covariances. 
# https://doi.org/10.1038/s41559-021-01509-w
# (Orkney, Bjarnason, et al., 2021)


library( geomorph )
library( ape )
library( nlme )
# Load required analytical R packages

setwd('D:/Documents/Alex_birds/Ideas_2022/SICBY')
# Set work directory as appropriate

load('Csize.29.07.2022.RData')
# Load centroid sizes for study birds. 
# In our study, we defined an this object as an array called 'GPA.Csize'
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

# An array of residual centroid sizes is now available for investigation. 


