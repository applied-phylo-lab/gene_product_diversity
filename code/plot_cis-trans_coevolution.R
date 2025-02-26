# Analysis results from simulations of cis-trans coevolution (produced by sim_cis-trans_coevolution.R) and make plots

library(ggplot2)

n1=100 # Number of deleterious modification events
opt.trans=2 # Optimum for trans-genotypic value Q for overall modification activity
l.all=c(2,5,10) # All values to consider for the number of cis-loci
sig.all=c(2,20) # All values to consider for width of fitness function for selection on trans-genotypic value Q
Ne.all=10^c(2,2.5,3,3.5,4,4.5,5) # Effective population sizes to consider

# Go through data files, extract information, and generate a summary data file
out=matrix(0,nrow=length(l.all)*length(sig.all)*length(Ne.all),ncol=10)
out[,1]=n1;out[,4]=opt.trans
row=1
for(l in l.all){
	dir1=paste("./l=",l,sep="") # Go to directory corresponding to the right cis-loci number
	setwd(dir1)
	for(sig in sig.all){
		dir2=paste("./sig=",sig,sep="") # Go to directory corresponding to the right selection on Q
		setwd(dir2)
		for(Ne in Ne.all){ # Go through values of Ne
			out[row,2]=l;out[row,3]=as.integer(l/2)
			out[row,5]=sig
			out[row,6]=log10(Ne)
			fn1=paste("lv_out_",n1,"_",log10(Ne),".txt",sep="") # Data file for modification levels (each column for a modification event)
			fn2=paste("trans_out_",n1,"_",log10(Ne),".txt",sep="") # Data file for end-point Q (1 column)
			fn3=paste("fitness_out_",n1,"_",log10(Ne),".txt",sep="") # Data file for end-point fitness (1 column)
			d1<-read.table(fn1,sep="\t");v=rep(0,100);for(i in 1:100){v[i]=median(d1[,i])};out[row,7]=median(v)
			d2<-read.table(fn2,sep="\t");out[row,8]=median(d2[,1]);out[row,9]=mean(d2[,1])
			d3<-read.table(fn3,sep="\t");out[row,10]=median(d3[,1])
			row=row+1
		}
		setwd("..")
	}
	setwd("..")
}
colnames(out)=c("ngene","l","v0","opt_Q","sig_Q","log10Ne","median_median_lv","median_Q","mean_Q","median_fitness")
write.table(out,file="sum_gaussian.txt",sep="\t")

d<-read.table("sum_gaussian.txt",sep="\t",header=TRUE)
d$sig=factor(d$sig,levels=sort(unique(d$sig)),ordered=TRUE)

# Plot Q against Ne
for(l.test in l.all){
	dsub=d[which(d$l==l.test),]
	g<-ggplot(dsub,aes(x=log10Ne,y=mean_Q,colour=sig))
	g=g+geom_point()+geom_line(linewidth=2)+scale_color_manual(values=c("red","blue"))
	g=g+theme_classic()
	g=g+xlab(expression(paste("Lo",g[10],N[e])))+ylab("Mean Q")
	g=g+theme(axis.text=element_text(size=12),axis.title=element_text(size=15),legend.text=element_text(size=15),legend.title=element_text(size=15),legend.key.size=unit(0.7,'cm')) # Adjust font size of axis labels, axis text, and legend
	fn.out=paste("plot_trans_gaussian_l=",l.test,".pdf",sep="")
	ggsave(fn.out,plot=g,width=6,height=4)
}

# Function to calculate a summary statistic for degree of sharing
# Input: data matrix containing modification level in each lineage, an identification cutoff
# All genes have the same alpha & gamma parameters, so the same modification level cutoff is applied to all genes
share.frac <- function(d,cutoff){
	share=rep(0,ncol(d))
	for(i in 1:ncol(d)){
		share[i]=length(which(d[,i]>cutoff))/nrow(d)
	}
	return(median(share))
}

# Examine pattern of sharing among lineages
cutoff=0.005 # Detection cutoff; a modification is considered as shared by 2 lineages if both have modification levels above the cutoff
T.all=c(1,2,4,6,8) # Timescales to consider
l.test=2;sig=2 # Examine a subset with a given number of cis-loci and strength of selection on Q
dir=paste("./l=",l.test,"/sig=",sig,sep="") # Go to the directory for the right parameter combination
setwd(dir)
out=matrix(0,nrow=length(Ne.all),ncol=length(T.all)) # Data matrix to write output (each row for a value of Ne, each column for a timescale)
for(i in 1:length(T.all)){
	T=T.all[i]
	dir=paste("./T=",T,sep="") # Go to the directory for the right timescale
	setwd(dir)
	for(j in 1:length(Ne.all)){ # Go through Ne values
		Ne=Ne.all[j]
		fn=paste("lv_out_",n1,"_",log10(Ne),".txt",sep="") # Data file for modification levels
		d<-read.table(fn,sep="\t");d=data.matrix(d)
		out[j,i]=share.frac(d,cutoff)
	}
	setwd("..")
}
# Add results from maximum timespan
out=cbind(out,rep(0,length(Ne.all)))
for(j in 1:length(Ne.all)){
	Ne=Ne.all[j]
	fn=paste("lv_out_",n1,"_",log10(Ne),".txt",sep="")
	d<-read.table(fn,sep="\t");d=data.matrix(d)
	out[j,ncol(out)]=share.frac(d,cutoff)
}
# Add starting point
out=cbind(rep(1,nrow(out)),out)
colnames(out)=c(0,T.all,10)
# Rearrange into a format that is better suited for plotted
out.new=matrix(0,nrow=ncol(out)*length(Ne.all),ncol=3)
row=1
for(i in 1:ncol(out)){
	T=as.numeric(colnames(out)[i])
	for(j in 1:length(Ne.all)){
		Ne=Ne.all[j]
		out.new[row,1]=T*1e7
		out.new[row,2]=log10(Ne)
		out.new[row,3]=out[j,i]
		row=row+1
	}
}
colnames(out.new)=c("T","log10Ne","share")

# Plot degree of sharing against timescale of evolution
d=data.frame(out.new)
d$log10Ne=factor(d$log10Ne,levels=sort(log10(Ne.all)),ordered=TRUE) # Factorize Ne for plotting
g=ggplot(d,aes(x=T,y=share,colour=log10Ne))
g=g+geom_point()+geom_line(linewidth=2,alpha=0.8)
g=g+theme_classic()
g=g+xlab("Time")+ylab("Sharing")+scale_x_continuous(breaks=seq(0,1e8,2e7))
g=g+theme(axis.text=element_text(size=12),axis.title=element_text(size=15),legend.text=element_text(size=15),legend.title=element_text(size=15),legend.key.size=unit(0.7,'cm')) # Adjust font size of axis labels, axis text, and legend

fn.out=paste("plot_share_l=",l.test,"_sig=",sig,".pdf",sep="")
setwd("..");setwd("..")
ggsave(g,file=fn.out,width=6,height=4)


