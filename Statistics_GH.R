#!/usr/bin/env Rscript

library(RColorBrewer)
library(dunn.test)
library(ggplot2) 
library(ggsignif) 
library(gridExtra)

dataset <- function(path){
mito=read.csv(path, header=TRUE, sep=",", quote = "\"", na.strings=c("","NA") )
return(mito)
}

boxplot <- function(a,mito,vet_groups,kruskal,ret)
	{
	stats=colnames(mito)[a]
	#The first column (data) contains the static values and the second column (clade) contains the groups
	df=data.frame(data=c(),clade=c())

	
	vet_clades=c()
	unlist_g=unlist(vet_groups)
	#For each clade I obtain the rank
	for (b in unlist_g){
		vet_clades=append(vet_clades,colnames(mito)[apply(mito, 2, function(x) any(x==b & !is.na(x)))  ])
		}
		
	#Retriving the static for each clade and adding them to df
	for (b in seq(1,length(unlist_g))) {
		data=mito[,a][ mito[,vet_clades[b]] == unlist_g[b] & !is.na(mito[,vet_clades[b]]) ]
		clade=rep(unlist_g[b],length(data))
		df = data.frame( data=c(df$data,data), clade=c(as.character(df$clade),clade) )
		}
	
	#If some clades are merged in a single group, I change the name of the clades with the name of the group they belong to
	for (b in seq(1,length(vet_groups))) {
		df$clade[df$clade %in% vet_groups[[b]]] = names(vet_groups)[b]
		}
	
	df = data.frame( data=df$data, clade=factor(df$clade,levels=names(vet_groups) ))

	sp=split(df$data,df$clade)
	
	#Calculating the max and min limits of the whiskers
	max=max(unlist(lapply(sp, function(x) quantile(x,0.75)+(1.5*(quantile(x,0.75)-quantile(x,0.25)) )    )))
	min=min(unlist(lapply(sp, function(x) quantile(x,0.25)-(1.5*(quantile(x,0.75)-quantile(x,0.25)) )    )))
	diff=max-min
	
	#I get rid of the outliers in order to avoid errors with ggplot
	df$data[df$data > max+diff*0.1] = max+diff*0.1
	df$data[df$data < min-diff*0.1] = min-diff*0.1


	#Calculating the numbers of samples for each group and creating the labels for the plot
	labels=c()
	summ=summary(factor(df$clade))
	for (b in seq(1,length(summ) ) ) {
		labels=append(labels,paste(names(summ[b]),"\n(N:",unname(summ[b]),")",sep="") )
		}

	vet_color=rep(c(brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(12,"Set3")),3)[1:length(levels(df$clade))]
	
	if(kruskal){
	#Function to calculate the significance between groups
	dunn_t=box_test(df$data,df$clade)	
	
	#Removing the non-significance pairs, which won't be displayed in the plot
	dunn_t$comparisons=dunn_t$comparisons[dunn_t$signif != ""]
	dunn_t$signif=dunn_t$signif[dunn_t$signif != ""]
	pos=c()
	
	#Creating a vector with the position for each significance bar. 
	for (b in seq(1,length(dunn_t$comparisons)) ){
		pos=append(pos,max+diff*0.2)
		max=max+diff*0.2
		}
	#Obtaining the plot with significance bars
	tot_plot=ggplot(df, aes_string(x="clade", y="data") ) + geom_boxplot(fill=vet_color, colour="black", outlier.shape= NA) +  geom_signif(comparisons = dunn_t$comparisons, annotations= dunn_t$signif, y_position=pos ) + scale_x_discrete(labels = labels )  + theme(axis.text=element_text(size=10), axis.title=element_text(size=25,face="bold") )  + scale_y_continuous(limits=c(min-diff*0.1,pos[length(pos)]+diff*0.1) )+ ylab(stats) + xlab("") 
	}
	else{
	
	#Obtainig the plot without significance bars
	tot_plot=ggplot(df, aes_string(x="clade", y="data") ) + geom_boxplot(fill=vet_color, colour="black", outlier.shape= NA ) + scale_x_discrete(labels = labels )  + theme(axis.text=element_text(size=10), axis.title=element_text(size=25,face="bold") )+ ylab(stats) + xlab("")   + scale_y_continuous(limits=c(min-diff*0.1,max+diff*0.1) ) 
	
		}	
	if(ret){
		return(tot_plot)
		}
	
	plot(tot_plot)
	cat("Do you want to save it? (F o T)") 
	salva=NA
	#A loop is needed since with Rstudio the first input is skipped most of the times
	while(is.na(salva)) {salva <- as.logical(readline())}
	if (salva) {
		cat("Name it: ")
		nome <-	readline()
		ggsave(paste(nome,".pdf",sep=""), plot=tot_plot)
		}
		 
	}


box_test<-function(data,groups){
	#Calculating the Kruskal-Wallis test
	kruskal=kruskal.test(data~groups)$p.value
	dun=dunn.test(data,g=groups, method="bonferroni")
	signif=c()
	#Taking notes of the pairwise significant groups
	if (kruskal <= 0.05) {
		for (a in dun$P.adjusted) {
			if (a > 0.05){
				signif=append(signif,"")
			} else if (a <= 0.05 & a > 0.01){
				signif=append(signif,"*")
			} else if (a <=0.01 & a > 0.001){
				signif=append(signif,"**")
			} else  {
				signif=append(signif,"***")
				}
		}
	}
	else {
		signif=rep("",length(dun$comparisons))
	}
	#I also built a table, when we have a lot of groups it is more easy to display the table instead of the significance bars
	vet_groups=levels(factor(groups))	
	tab=matrix(NA,length(vet_groups),length(vet_groups))
	colnames(tab)=vet_groups
	rownames(tab)=vet_groups
	for (a in seq(1,length(dun$comparisons)) ){
		comp=strsplit(dun$comparisons[a]," - ")[[1]]
		tab[comp[1],comp[2]]=paste(round(dun$P.adjusted[a],3),signif[a],paste="")
		tab[comp[2],comp[1]]=paste(round(dun$P.adjusted[a],3),signif[a],paste="")
		}
		
	#Creating the list dunn_t, which contains list of comaprisons, the significance for each group and the table
	dunn_t=list(comparisons=strsplit(dun$comparisons," - "), signif=signif, tab=tab)

	return(dunn_t)
			
}


corr <- function(mito,sub) {
	#Obtaining the rank of the clade
	rank=colnames(mito)[apply(mito, 2, function(x) any(x==sub & !is.na(x)))  ]	
	#Extracting the sub-dataset
	sub_data=mito[mito[,rank]==sub & !is.na(mito[,rank]),]
	library(Hmisc)
  	library(corrplot)
	mito_corr=rcorr(as.matrix(sub_data[,c(13,3:10,12)]),type="spearman")
	#Bonferroni correction
	mito_corr$P=mito_corr$P*length(sub_data$index)
	mito_corr$sub=sub		
	
	corrplot(mito_corr$r,type="upper",p.mat=mito_corr$P, sig.level=0.05, insig="blank",title=mito_corr$sub,method="pie",diag=F, mar=c(0,0,1,0))
	
	
	cat("Do you want to save it? (F o T)")
	salva=NA
	#A loop is needed since with Rstudio the first input is skipped most of the times
	while(is.na(salva)) {salva <- as.logical(readline())}
	
	if (salva) {
		cat("Name it: ")
		nome <-	readline()
		pdf(file=paste(nome,".pdf",sep=""))
		corrplot(mito_corr$r,type="upper",p.mat=mito_corr$P, sig.level=0.001, insig="blank", title=sub,method="pie",diag=F, mar=c(0,0,1,0))

		dev.off()
		}

		 

}

pca <- function(mito,sub,clade_groups,num_ellissi)
	{

	pca1=mito[,c(3:10,12:13)]

	clade=colnames(mito)[apply(mito, 2, function(x) any(x==sub & !is.na(x)))  ]
	col=which(colnames(mito) == clade)
	col_groups=which(colnames(mito) == clade_groups)
	sub_data=pca1[mito[,col] == sub & !is.na(mito[,col]),]
	
	
	library(ade4)
	#Just the first 4 PCs can be displayed, increase nf if you want to display the others
	pca=dudi.pca(sub_data,scale=T, center=T, scannf= F, nf=4)

	barplot(pca$eig)

	#In this vector is reported at which sub_clade belongs each sample
	groups=mito[,col_groups][mito[,col] == sub & !is.na(mito[,col]) ]


	#NA are called Others
	if ( any(is.na(groups)) ) {
	groups=factor(groups, levels=c(names(sort(summary(factor(groups)),decreasing=T)),"Others"))
	groups[is.na(groups)]="Others"
	groups=factor(groups)
	}
	else
	{
	groups=factor(groups, levels=names(sort(summary(factor(groups)),decreasing=T)) )
	groups=factor(groups)
	}



	

	vet_color=rep(c(brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(12,"Set3")),4)[1:length(levels(groups))]
	labels=levels(groups)
	
	#Defining for which sub_clades will be calculated the ellipse 
	if (length(levels(groups)) > num_ellissi) {

		five_clades=levels(groups)[1:num_ellissi]		
		count=unname(summary(groups)[names(summary(groups)) %in% five_clades])
		labels[labels %in% five_clades]=paste(labels[labels %in% five_clades]," N=",count,sep="")		
			
		pca_gr=pca$li[groups %in% five_clades,]
		group_b=groups[groups %in% five_clades]
		colour_b=vet_color[1:num_ellissi]
	}
	else {
		group_b=groups
		pca_gr=pca$li
		colour_b=vet_color[1:length(levels(groups))]
		count=unname(summary(groups))
		labels=paste(labels," N=",count,sep="")
		
		}
	plot_abomba=T
	while(plot_abomba){

		cat("Choose which PC do you want to display (Insert arabic numbers between 1 and 4)\nx:")
		pcax_n = ""
		while (pcax_n==""){	pcax_n <- readline() }
		pcax=paste("Axis",pcax_n,sep="")
		pcax_s=paste("PC ",pcax_n,sep="")
		pcx_cs=paste("CS",pcax_n,sep="")		
		
		cat("y:")
		pcay_n <- readline()
		pcay=paste("Axis",pcay_n,sep="")
		pcay_s=paste("PC ",pcay_n,sep="")
		pcy_cs=paste("CS",pcay_n,sep="")
		
		#The explained variability
		tot=round((sum(pca$eig[as.integer(pcax_n)],pca$eig[as.integer(pcay_n)])/sum(pca$eig))*100,2)
		print(paste("The variability explained by these two PC is ", tot, sep=""))		
		
		pca_plot=ggplot(pca$li, aes_string(x=pcax, y=pcay, colour=groups)) + geom_point()  + stat_ellipse(data=pca_gr, geom="polygon", aes_string(x=pcax, y=pcay, color=group_b, fill=group_b),level= 0.95, inherit.aes=F, show.legend=F, alpha = 0.2) + ylab(pcay_s) + xlab(pcax_s) + scale_color_manual(values=vet_color, labels=labels) + scale_fill_manual(values=colour_b) + labs(color=clade_groups, title=sub, caption=paste("variance explained: ",tot,"%",sep="") ) 
		
		#Function for the circle plot
		circ <- circleFun(c(0,0),2,npoints = 500)
		#Setting the labels
		label_string=paste(rownames(pca$c1),collapse='","')
		label_geom_text=paste('c("',label_string,'")')

		plot(pca_plot)

		#I ask if the user want to display the circle plot 
		cat("Do you want to plot also the correlation circle? (T o F)")
		like <- as.logical(readline())
		if (like) {	
		cat("Choose whete to place the corrplot \nx0: ")
		x0 <- as.numeric(readline())	
		cat("x_end:")
		x_end <- as.numeric(readline())
		cat("y0:")
		y0 <- as.numeric(readline())
		cat("y_end:")
		y_end <- as.numeric(readline())
		
		c1_plot=ggplot() + geom_path(data = circ,aes(x,y), lty = 2, color = "grey", alpha = 0.7) + geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) + geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9)+geom_segment( data= pca$c1, aes_string(x=0, xend=pcx_cs, y=0, yend=pcy_cs), color="gray43", alpha=0.7, arrow = arrow(length = unit(0.025, "npc"), type = "open"), lwd =1) +   ylab(pcay_s) + theme_minimal() + theme(panel.grid = element_blank(), panel.border = element_rect(fill= "transparent")) + geom_text(data = pca$c1, aes_string( x=pcx_cs, y=pcy_cs, label=label_geom_text), check_overlap = F, size = 3) +xlab(pcax_s) + theme(axis.title=element_text(size=5))

		tot_plot=pca_plot + annotation_custom(ggplotGrob(c1_plot), xmin = x0, xmax = x_end, ymin = y0, ymax = y_end)
		}
		else { tot_plot=pca_plot}
		plot(tot_plot)
		
		cat("Do you want to save it (T o F)")
		like <- as.logical(readline())
		if (like) {
			cat("Name it:")
			name <- readline()
			ggsave(paste(name,".pdf",sep=""))
			}
			

		cat("Do you want to try with other PCs? (T o F)")
		plot_abomba <- as.logical(readline())
	}
} 

circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}


