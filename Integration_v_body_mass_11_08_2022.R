# This script will produce csv files for 10 different bins of body mass, 
# showing how integration among pairwise combinations of bones changes with mass.
# For the studied birds (149 species), these bins will be overlapping, 
# which will generate statistical non-independence between some bins. 
# This can be regarded as being similar to a rolling-average.
# Because the birds are not uniformly distributed in frequency over log10(bodymass), 
# some bins will have wider or narrower ranges of mass represented than others. 

# You will need to run 'Residual_Csize_10_08_2022.R' 
# to first obtain the array of residual centroid sizes. 

library(abind)
# The 'array bind' function will be needed. 


bones<-names(residual.Csize)
# Retrieve the names of the skeletal elements

# The following specifies a function 
# that computes pairwise Z-scores of integration between 
# different combinations of bones, for a stated group of taxa and a provided phylogeny relating the taxa.

pair.int.phy <- function( input, phylogeny, taxa ){

	newphy<-drop.tip( phylogeny, phylogeny$tip.label[ !phylogeny$tip.label  %in% taxa] )
	
	bones <- names(input)
	output <- matrix(NA,length(bones),length(bones))
	output.p <- matrix(NA,length(bones),length(bones))
	rownames(output)<-bones
	colnames(output)<-rownames(output)
	colnames(output.p)<-colnames(output)
	rownames(output.p)<-rownames(output)
	bone_combinations<- combn(bones,2,simplify=T)
	
	int_list<-list()
	for(i in 1:dim(bone_combinations)[2] ){	
		x<-paste(bone_combinations[,i][1])
		y<-paste(bone_combinations[,i][2])
		int_list[[i]] <- phylo.integration(as.matrix(input[[which(bones==x)]][taxa]),as.matrix(input[[which(bones==y)]][taxa]),phy=newphy,iter=999,print.progress=F)
	}

	comparison<-compare.pls(int_list)
	
	z.score<-unlist(comparison$sample.z)

	p.val<-unlist( lapply(int_list,function(x) x$P.value ) )

	for(i in 1:dim(bone_combinations)[2] ){
		x<-paste(bone_combinations[,i][1])
		y<-paste(bone_combinations[,i][2])
		output[which(rownames(output)==x),which(colnames(output)==y)]<-z.score[i]
		output[which(rownames(output)==y),which(colnames(output)==x)]<-z.score[i]
	}

	for(i in 1:dim(bone_combinations)[2] ){
		x<-paste(bone_combinations[,i][1])
		y<-paste(bone_combinations[,i][2])
		output.p[which(rownames(output.p)==x),which(colnames(output.p)==y)]<-p.val[i]
		output.p[which(rownames(output.p)==y),which(colnames(output.p)==x)]<-p.val[i]
	}
		
	return(abind( output, output.p,along=3 ))
}


breaks<-round(seq(40,length(masses),length.out=10))
# This defines the upper-bounds of the overlapping bins for body mass. 

bin.width<-40
# This will be the number of taxa in any individual bin. 

n<-30 
# This will be the number of replicates for the test statistic. 

subsample.size<-30
# Each replicate will be computed by subsampling 30 taxa within each bin.  

k <- dim(combn(names(residual.Csize),2,simplify=2))[2]
# This is the number of pairwise combinations of skeletal elements.

size.scale.integration<-matrix( NA, length(breaks)*n, k+4 ) 
size.scale.integration.p<-matrix( NA, length(breaks)*n, k+4 ) 
# These lines define a matrix to receive analysis output.

break.vector<-rep(breaks,each=n)
# This defines a vector of indices for the overlapping bins of body masses.

for( i in 1 : dim(size.scale.integration)[1] ){
# For each replicate, within each bin...

	end<-break.vector[i]
	start<-end-bin.width
	# Define the bin of body mass.

	bin.mass <- mean(masses[order(masses)][start:end])
	bin.birds <- names(masses[order(masses)][start:end])
	# Compute the mean mass and name the birds within the bin. 

	selection<-sample(bin.birds,subsample.size)
	# Extract a subsample of taxa.

	subsample.mass <- mean(masses[order(masses)][start:end][selection])
	min.mass <- min(masses[order(masses)][start:end][selection])
	max.mass <- max(masses[order(masses)][start:end][selection])
	# Compute the mean, min and max body masses within the subsample of taxa. 

	size.scale.integration[i,1]<-bin.mass
	size.scale.integration[i,2]<-subsample.mass
	size.scale.integration.p[i,1]<-bin.mass
	size.scale.integration.p[i,2]<-subsample.mass
	size.scale.integration[i,3]<-min.mass
	size.scale.integration[i,4]<-max.mass
	size.scale.integration.p[i,3]<-min.mass
	size.scale.integration.p[i,4]<-max.mass
	# Add these information to the output matrix. 

	result<-pair.int.phy( residual.Csize, pruned.tree, selection)
	# Compute pairwise integration between all combinations of skeletal elements. 

	size.scale.integration[i,-c(1:4)] <- result[,,1][upper.tri(result[,,1])]
	size.scale.integration.p[i,-c(1:4)] <- result[,,2][upper.tri(result[,,2])]
	# Add the resultant test statistics to the data output matrix. 

	print(round(i/dim(size.scale.integration)[1]*100,digit=2))
	# Print the progress to the screen (this code will take some time to run).
}

# This has produced a matrix of Z-scores and p-values as body mass changes

ind<- which( upper.tri(result[,,1],diag=F) , arr.ind = TRUE )
# This is the index for the upper triangle of pairwise integration scores.
# The lower triangle is identical, so we do not need it. 

colpaste<-paste(bones[ind[,1]],bones[ind[,2]],sep='.')
# Define column names for the data output matrix. 

colnames(size.scale.integration) <- c('bin.mass','subsample.mass','min.mass','max.mass',colpaste)
colnames(size.scale.integration.p) <- c('bin.mass','subsample.mass','min.mass','max.mass',colpaste)
# Set the column names. 


getwd()
# Check the working directory is as expected.

write.csv( x=size.scale.integration, file='size_scale_integration_11_08_2022.csv')
write.csv( x=size.scale.integration.p, file='size_scale_integration_p_11_08_2022.csv')
# Save the output for subsequent analyses. 
