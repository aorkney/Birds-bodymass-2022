# This script will produce csv files for 10 different levels of sample size, 
# and 10 different levels of body mass,
# showing how integration among pairwise combinations of bones evolves over these variables. 

# You will need to run 'Residual_Csize_10_08_2022.R' 
# to first obtain the array of residual centroid sizes. 

library(abind)
# The 'array bind' function will be needed. 

bones<-names(residual.Csize)


# The following specifies a function 
# that computes pairwise Z-scores of integration between 
# different combinations of bones, for a stated group of taxa. 

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


# The logarithm of the bird masses is a normal distribution, 
# therefore the number of birds with very extreme masses is very small
# and the number with 'average' masses is quite large.
# We need to have confidence that any changes in intergration strength across the different size bins we select
# are not driven by small sampling numbers. 
# Therefore we need to assess what the minimum sample size we can select is, without severely inflating uncertainty.

# We shall do this by taking iterative samples of decreasing size, beginning at 10, and increasing to 100. 
# We should draw 30 random samples of each size
# and we should compute the resultant matrices of integration for each one
# and store them in a matrix
# with each row represented by 'result[upper.tri(result)]'
# and a corresponding value for the sample size.

sample.size <- seq(10,100,by=10)
n <- 30 # number of replicates
k <- dim(combn(names(residual.Csize),2,simplify=2))[2]

size.scale.uncertainty<-matrix( NA, length(sample.size)*n, k+1 ) 
size.scale.uncertainty.p<-matrix( NA, length(sample.size)*n, k+1 )
# where 'k' is the number of pairwise combinations of bones

sample.size.vector<-rep(sample.size,each=n)

for( i in 1 : dim(size.scale.uncertainty)[1] ){
	size.scale.uncertainty[i,1] <- sample.size.vector[i]
	selection<-sample(pruned.tree$tip,sample.size.vector[i])
	result<-pair.int.phy( residual.Csize, pruned.tree, selection)
	size.scale.uncertainty[i,-1] <- result[,,1][upper.tri(result[,,1])]
	size.scale.uncertainty.p[i,-1] <- result[,,2][upper.tri(result[,,2])]
	print(round(i/dim(size.scale.uncertainty)[1]*100,digit=2))
}

ind<- which( upper.tri(result[,,1],diag=F) , arr.ind = TRUE )

colpaste<-paste(bones[ind[,1]],bones[ind[,2]],sep='.')
colnames(size.scale.uncertainty) <- c('sample.size',colpaste)
colnames(size.scale.uncertainty.p) <- c('sample.size',colpaste)

# This has made a pair of matrices representing the Z-score of pairwise integration of
# allometrically-adjusted bone sizes, varying as a function of sample size. 



write.csv( x=size.scale.uncertainty, file='size_scale_uncertainty_10_08_2022.csv')
write.csv( x=size.scale.uncertainty.p, file='size_scale_uncertainty_p_10_08_2022.csv')
# Save the output for subsequent analysis. 




