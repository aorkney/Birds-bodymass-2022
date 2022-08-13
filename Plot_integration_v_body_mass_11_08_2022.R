# This script will produce plots demonstrating how integration with 
# avian body modules changes with increasing body mass. 

# This script requires that previous analysis output have been produced and saved.
# This can be achieved by...
# Running 'Residual_Csize_10_08_2022.R' 
# to first obtain the array of residual centroid sizes,
# for the desired selection of taxa. 
# Then running 'Integration_v_body_mass_11_08_2022.R'
# to compute the pairwise integration across the skeleton 
# over several intervals of body mass. 
# These individual scripts should be consulted for furhter details
# or customisation. 

setwd('D:/Documents/Alex_birds/Ideas_2022/SICBY')
# Set the working directory.

library(ggplot2)
library(reshape2)
library(ggpubr)
library(ggnewscale)
# Load packages to organise and plot data. 

size.int<-read.csv('size_scale_integration_11_08_2022.csv')
size.int<-size.int[,-1]
size.int[size.int < 0 ]<-0
size.int.p<-read.csv('size_scale_integration_p_11_08_2022.csv')
size.int.p<-size.int.p[,-1]
size.int.p[size.int.p < 0 ]<-0
# Load and prepare data. 

head.bones<-
combn(c('skull','mandible'),2) 
within.head<-matrix(NA,1,dim(head.bones)[2])
for(i in 1:dim(head.bones)[2] ){
	within.head[,i]<-intersect(grep(head.bones[1,i],colnames(size.int)),
	grep(head.bones[2,i],colnames(size.int)) )
}
within.head<-as.vector(within.head)
# Define the indices for the bones in the head module.

wing.bones<-
combn(c('humerus','radius','ulna','carpometacarpus'),2) 
within.wing<-matrix(NA,1,dim(wing.bones)[2])
for(i in 1:dim(wing.bones)[2] ){
	within.wing[,i]<-intersect(grep(wing.bones[1,i],colnames(size.int)),
	grep(wing.bones[2,i],colnames(size.int)) )
}
within.wing<-as.vector(within.wing)
# Define the indices for the bones in the wing module. 

trunk.bones<-
combn(c('scapula','coracoid','sternum','synsacrum'),2) 
within.trunk<-matrix(NA,1,dim(trunk.bones)[2])
for(i in 1:dim(trunk.bones)[2] ){
	within.trunk[,i]<-intersect(grep(trunk.bones[1,i],colnames(size.int)),
	grep(trunk.bones[2,i],colnames(size.int)) )
}
within.trunk<-as.vector(within.trunk)
# Define the indices for the bones in the trunk module.

leg.bones<-
combn(c('femur','tibiotarsus','tarsometatarsus'),2) 
within.leg<-matrix(NA,1,dim(leg.bones)[2])
for(i in 1:dim(leg.bones)[2] ){
	within.leg[,i]<-intersect(grep(leg.bones[1,i],colnames(size.int)),
	grep(leg.bones[2,i],colnames(size.int)) )
}
within.leg<-as.vector(within.leg)
# Define the indices for the bones in the leg modules. 


inc<-unique(size.int$bin.mass)
mean.size.int<-matrix(NA,length(inc),21)
mean.size.int[,1]<-log10(inc)
colnames(mean.size.int)<-c('bin.mass',
'head','headmin','headmax','headminmin','headmaxmax',
'wing','wingmin','wingmax','wingminmin','wingmaxmax',
'trunk','trunkmin','trunkmax','trunkminmin','trunkmaxmax',
'leg','legmin','legmax','legminmin','legmaxmax')
# Define a matrix to recieve binned values

for(i in 1: length(inc) ){
	ind<-which(size.int$bin.mass==inc[i])
	mean.size.int[i,2]<-mean( size.int[ind,within.head] )

	boot<-list()
	for(j in 1:1000){
		 boot[[j]]<-sample(size.int[ind,within.head],replace=T)
	}
	boot<-unlist(boot)
	mean.size.int[i,3]<-quantile(boot,probs=0.32)
	mean.size.int[i,4]<-quantile(boot,probs=0.68)
	mean.size.int[i,5]<-quantile(boot,probs=0.05)
	mean.size.int[i,6]<-quantile(boot,probs=0.95)

	mean.size.int[i,7]<-mean( rowMeans( size.int[ind,within.wing] ))

	boot<-list()
	for(j in 1:1000){
		 boot[[j]]<-sample(rowMeans(size.int[ind,within.wing]),replace=T)
	}
	boot<-unlist(boot)
	mean.size.int[i,8]<-quantile(boot,probs=0.32)
	mean.size.int[i,9]<-quantile(boot,probs=0.68)
	mean.size.int[i,10]<-quantile(boot,probs=0.05)
	mean.size.int[i,11]<-quantile(boot,probs=0.95)

	mean.size.int[i,12]<-mean( rowMeans( size.int[ind,within.trunk] ))

	boot<-list()
	for(j in 1:1000){
		 boot[[j]]<-sample(rowMeans(size.int[ind,within.trunk]),replace=T)
	}
	boot<-unlist(boot)
	mean.size.int[i,13]<-quantile(boot,probs=0.32)
	mean.size.int[i,14]<-quantile(boot,probs=0.68)
	mean.size.int[i,15]<-quantile(boot,probs=0.05)
	mean.size.int[i,16]<-quantile(boot,probs=0.95)

	mean.size.int[i,17]<-mean( rowMeans( size.int[ind,within.leg] ))

	boot<-list()
	for(j in 1:1000){
		 boot[[j]]<-sample( rowMeans(size.int[ind,within.leg]),replace=T)
	}
	boot<-unlist(boot)
	mean.size.int[i,18]<-quantile(boot,probs=0.32)
	mean.size.int[i,19]<-quantile(boot,probs=0.68)
	mean.size.int[i,20]<-quantile(boot,probs=0.05)
	mean.size.int[i,21]<-quantile(boot,probs=0.95)

}
# Populate the matrix with mean within-module integration statistics, 
# as well as 2 and 1 sigma confidence envelope statistics. 

df<-as.data.frame(mean.size.int)
molten.df<-melt(df,id='bin.mass')
# Re-organise the data.

ribbons<-c(grep('min',molten.df$variable),grep('max',molten.df$variable))
# Define the confidence envelopes.

melt.df<-melt(df[-ribbons,],id='bin.mass')
# Re-organise a version of the data with confidence envelopes.

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# This is a colour blind friendly palette. 

mp<-
ggplot()+
geom_line(aes(x=bin.mass,y=value,group=variable,linetype=variable,col=variable),size=1,data=molten.df[-ribbons,])+
scale_colour_manual(values=c('head'=cbbPalette[7],'wing'=cbbPalette[4],
'trunk'=cbbPalette[6],'leg'=cbbPalette[8]))+
labs(colour='module',linetype='module')+
theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=15,colour='black'),
      axis.text.y=element_text(size=15,colour='black'),
	axis.title.x.bottom=element_text(size=20,colour="black"),
	axis.title.y.left=element_text(size=20,colour="black"),
      axis.ticks=element_line(size=2),
      legend.position="right",
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank(),
	legend.key.width = unit(2, "cm"))
# Produce a line plot (the purpose of this plot is just to extract its legend).

mp.legend<-cowplot::get_legend(mp)
# Make sure the cowplot package is installed. 



wing.mp<-
ggplot()+
geom_ribbon(aes(x=bin.mass,ymin=wingminmin,ymax=wingmaxmax),alpha=1,data=df,bg="#80FEF3")+
geom_ribbon(aes(x=bin.mass,ymin=wingmin,ymax=wingmax),alpha=1,data=df,bg="#40CEB3")+
geom_line(aes(x=bin.mass,y=wing),size=1,col=cbbPalette[4], data=df,linetype='22')+
geom_vline(xintercept=log10(300),size=3/2)+
#lims(y=c(2.5,5.5),x=c(1,4))+ 
#lims(y=c(3,5))+ for telluraves
#lims(y=c(3,5.5),x=c(1,4))+
lims(y=c(2.5,5.5))+ 
labs(colour='module',linetype='module')+
labs(x=expression(paste(log[10],'(',italic(mass),')',' [g]')), y=expression(italic(Z)))+
  theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=15,colour='black'),
      axis.text.y=element_text(size=15,colour='black'),
	axis.title.x.bottom=element_text(size=20,colour="black"),
	axis.title.y.left=element_text(size=20,colour="black"),
      axis.ticks=element_line(size=2),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())

# Produce a plot illustrating how the mean pairwise Z-score of integration
# within the wing varies with increasing body mass across the studied birds. 

yvals<-c(2.5,2.9,3.3,3.5,3.6,3.7)
#yvals<-c(0,3.2,0,3.5,3.6,3.7) # for Telluraves
#yvals<-c(3.2,0,0,3.2,0,3.2) # for non-Telluraves
#yvals<-c(3.1,0,3.3,3.5,3.6,3.7) # for non-Apodiformes
kseq<-inc[seq(1,10,by=2)]
kseq<-inc[c(1,2,3,5,7,9)]
mp.comp<-wing.mp
for(k in 1:length(kseq)){
	recompile <- matrix(NA,13,13)
	rownames(recompile)<-c('skull','mandible','carpometacarpus','radius','ulna','humerus',
	'sternum','coracoid','scapula','synsacrum','femur','tibiotarsus','tarsometatarsus')
	colnames(recompile)<-rownames(recompile)

	for(i in 1:length(size.int[1,-c(1:2)])  ){
		names<-colnames(size.int[,-c(1:2)])
		x<-  sub("\\..*","", names)[i]
		y<-  sub(".*\\.","", names)[i]
		ind<-which(size.int$bin.mass==kseq[k])
		recompile[which(rownames(recompile)==x),which(colnames(recompile)==y)]<-mean( size.int[ind,-c(1:2)][,i] )
		recompile[which(rownames(recompile)==y),which(colnames(recompile)==x)]<-mean( size.int[ind,-c(1:2)][,i] )
	}

	result.rev<-apply(t(recompile),2,rev)
	longData<-melt(result.rev)
	longData$value[which(longData$value<0)]<-0
	longData$Var1<-(as.numeric(longData$Var1)/(100/(5.5-2.5)))+yvals[k]
	longData$Var2<-(as.numeric(longData$Var2)/(100/2))+log10(kseq[k])

	mp.comp<- mp.comp+
	coord_fixed(ratio=2/3)+
  	geom_raster(aes(fill=value,x = Var2, y = Var1),data=longData)+
   	#scale_fill_gradientn(colours = c('darkblue','red','yellow') ,na.value="white")+
	scale_fill_gradientn(colours = c('yellow','red','black') ,na.value="white",lim=c(0,6))+
	new_scale("fill") +
	theme(legend.position='none')
}

# Add overlaid microplots of matrices of pairwise integration at different
# representative ranges of body mass. 


head.mp<-
ggplot()+
geom_ribbon(aes(x=bin.mass,ymin=headminmin,ymax=headmaxmax),alpha=1,data=df,bg="#F5CE80")+
geom_ribbon(aes(x=bin.mass,ymin=headmin,ymax=headmax),alpha=1,data=df,bg="#F58E40")+
geom_line(aes(x=bin.mass,y=head),size=1,col=cbbPalette[7], data=df,linetype='solid')+
geom_vline(xintercept=log10(300),size=3/2)+
#lims(y=c(2.5,5.5),x=c(1,4))+ 
#lims(y=c(3,5))+ for telluraves
#lims(y=c(3,5.5),x=c(1,4))+
#lims(y=c(2.5,5.5))+ 
labs(colour='module',linetype='module')+
labs(x=expression(paste(log[10],'(',italic(mass),')',' [g]')), y=expression(italic(Z)))+
  theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=15,colour='black'),
      axis.text.y=element_text(size=15,colour='black'),
	axis.title.x.bottom=element_text(size=20,colour="black"),
	axis.title.y.left=element_text(size=20,colour="black"),
      axis.ticks=element_line(size=2),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())

# Produce a plot illustrating how the mean pairwise Z-score of integration
# within the head varies with increasing body mass across the studied birds. 

trunk.mp<-
ggplot()+
geom_ribbon(aes(x=bin.mass,ymin=trunkminmin,ymax=trunkmaxmax),alpha=1,data=df,bg="#80E2F2")+
geom_ribbon(aes(x=bin.mass,ymin=trunkmin,ymax=trunkmax),alpha=1,data=df,bg="#40A2F2")+
geom_line(aes(x=bin.mass,y=trunk),size=1,col=cbbPalette[6], data=df,linetype='31')+
geom_vline(xintercept=log10(300),size=3/2)+
#lims(y=c(2.5,5.5),x=c(1,4))+ 
#lims(y=c(3,5))+ for telluraves
#lims(y=c(3,5.5),x=c(1,4))+
labs(colour='module',linetype='module')+
labs(x=expression(paste(log[10],'(',italic(mass),')',' [g]')), y=expression(italic(Z)))+
  theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=15,colour='black'),
      axis.text.y=element_text(size=15,colour='black'),
	axis.title.x.bottom=element_text(size=20,colour="black"),
	axis.title.y.left=element_text(size=20,colour="black"),
      axis.ticks=element_line(size=2),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())

# Produce a plot illustrating how the mean pairwise Z-score of integration
# within the trunk varies with increasing body mass across the studied birds. 

leg.mp<-
ggplot()+
geom_ribbon(aes(x=bin.mass,ymin=legminmin,ymax=legmaxmax),alpha=1,data=df,bg="#FCE9F7")+
geom_ribbon(aes(x=bin.mass,ymin=legmin,ymax=legmax),alpha=1,data=df,bg="#FCA9E7")+
geom_line(aes(x=bin.mass,y=leg),size=1,col=cbbPalette[8], data=df,linetype='33')+
geom_vline(xintercept=log10(300),size=3/2)+
#lims(y=c(2.5,5.5),x=c(1,4))+ 
#lims(y=c(3,5))+ for telluraves
#lims(y=c(3,5.5),x=c(1,4))+
labs(colour='module',linetype='module')+
labs(x=expression(paste(log[10],'(',italic(mass),')',' [g]')), y=expression(italic(Z)))+
  theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_text(size=15,colour='black'),
      axis.text.y=element_text(size=15,colour='black'),
	axis.title.x.bottom=element_text(size=20,colour="black"),
	axis.title.y.left=element_text(size=20,colour="black"),
      axis.ticks=element_line(size=2),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())

# Produce a plot illustrating how the mean pairwise Z-score of integration
# within the leg varies with increasing body mass across the studied birds. 

# The following code will compute Ordinary Least Squares fits for change in pairwise-Z scores
# of integration across all combinations of skeletal elements as a function of mean body mass
# within each bin. 
# This is a plot produced for illustrative purposes; be aware that the bins of body mass
# may overlap, causing statistical non-independence among the data. 

recompile <- matrix(NA,13,13)
rownames(recompile)<-c('skull','mandible','carpometacarpus','radius','ulna','humerus',
'sternum','coracoid','scapula','synsacrum','femur','tibiotarsus','tarsometatarsus')
colnames(recompile)<-rownames(recompile)

inc<-unique(size.int[,1])
for(i in 1:length(size.int[1,-c(1:2)])  ){
	names<-colnames(size.int[,-c(1:2)])
	x<-  sub("\\..*","", names)[i]
	y<-  sub(".*\\.","", names)[i]
	ind<-which(size.int$bin.mass==kseq[k])

	temp<-list()
	for(j in 1:length(inc)){
		ind<-which(size.int[,1]==inc[j])
		temp[[j]]<-mean(size.int[,-c(1:2)][,i][ind])
		
	}
	temp<-unlist(temp)

	coef<-lm(temp~log10(inc))

	if(summary(coef)$coef[2,4] <1){ # You might also use a threshold of 0.05
		recompile[which(rownames(recompile)==x),which(colnames(recompile)==y)]<-coef$coef[2]
		recompile[which(rownames(recompile)==y),which(colnames(recompile)==x)]<-coef$coef[2]
	} 
}

result.rev<-apply(t(recompile),2,rev)
longData<-melt(result.rev)
# Re-organise the output. 

bone_colours<-rev(c(rep(cbbPalette[7],2),rep(cbbPalette[4],4),rep(cbbPalette[6],4),rep(cbbPalette[8],3)))
# Colours for different modules. 

mp.change<-
ggplot(longData, aes(x = Var2, y = Var1)) + 
	coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")+
  geom_tile(aes(fill=value)) + 
   scale_fill_gradient2(low='blue',mid='white',high='red' ,na.value="white",midpoint=0) +
	labs(fill = "slope")+
  labs(x="", y="") +
  theme_bw() + theme(axis.text.x=element_text(size=7, angle=90, vjust=0.3,colour=rev(bone_colours),face = "bold"),
                     axis.text.y=element_text(size=7,colour=(bone_colours),face = "bold"),
				legend.key.width=unit(1/2,'cm'),
				legend.title=element_text(size=20),
				plot.margin=unit(c(-0.5,-1,0.5,-1),'cm'),
                     plot.title=element_text(size=11,hjust=.5), legend.position="right",panel.background = element_blank())+
	labs(fill = expression(Delta))+
geom_hline(yintercept = 3.5,size=1,colour="black")+
geom_hline(yintercept = 11.5,size=1,colour="black")+
geom_hline(yintercept = 7.5,size=1,colour="black")+
geom_vline(xintercept = 10.5,size=1,colour="black")+
geom_line(data=as.data.frame(cbind(c(2.5,2.5),c(0.5,13.5))),aes(x=V1,y=V2),size=1,colour='black')+
geom_line(data=as.data.frame(cbind(c(6.5,6.5),c(0.5,13.5))),aes(x=V1,y=V2),size=1,colour='black')+
geom_line(data=as.data.frame(cbind(c(10.5,10.5),c(0.5,13.5))),aes(x=V1,y=V2),size=1,colour='black')+
scale_x_discrete(position = "top") 

# Produce a plot that shows general trends in pairwise-Z of integration across the skeleton as a function
# of increasing body mass. 

# The following lines will assemble all the subplots into one combined plot. 

mp.modules<-ggarrange(head.mp,trunk.mp,leg.mp,ncol=1,nrow=3,align='hv',
labels=c('b','d','e'),label.x=-0.025,label.y=1,font.label=list(size=30))
mp.left<-ggarrange(mp.comp,mp.change,ncol=1,
labels=c('a','c'),label.x=0.15,label.y=1,font.label=list(size=30))
dev.new(height=12,width=16.9,unit='cm')
ggarrange(mp.left,mp.modules,cowplot::plot_grid(mp.legend),ncol=3,widths=c(1,0.8,0.25))

# ggsave(filename='mass_v_integration_11_08_2022b.pdf')
# You may choose to save the plot. 

