setwd("/Users/daohanji/Desktop/gene_product_diversity/out/trans/")

n1=100;opt.trans=2
l.all=c(2,5,10)
sig.all=c(2,20)
Ne.all=10^c(2,2.5,3,3.5,4,4.5,5)

# Extract information and generate a summary data file
out=matrix(0,nrow=length(l.all)*length(sig.all)*length(Ne.all),ncol=10)
out[,1]=n1;out[,4]=opt.trans
row=1
for(l in l.all){
	dir1=paste("./l=",l,sep="")
	setwd(dir1)
	for(sig in sig.all){
		dir2=paste("./sig=",sig,sep="")
		setwd(dir2)
		for(Ne in Ne.all){
			out[row,2]=l;out[row,3]=as.integer(l/2)
			out[row,5]=sig
			out[row,6]=log10(Ne)
			fn1=paste("lv_out_",n1,"_",log10(Ne),".txt",sep="")
			fn2=paste("trans_out_",n1,"_",log10(Ne),".txt",sep="")
			fn3=paste("fitness_out_",n1,"_",log10(Ne),".txt",sep="")
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

library(ggplot2)

d<-read.table("sum_gaussian.txt",sep="\t",header=TRUE)
d$sig=factor(d$sig,levels=sort(unique(d$sig)),ordered=TRUE)

# Plot trans- gv against Ne
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

# Plot conservation against time
# A summary statistic for degree of sharing
# Input: data matrix containing modification level in each lineage, identification cutoff
# All genes have the same alpha & gamma parameters, so the same modification level cutoff is applied to all genes
share.frac <- function(d,cutoff){
	share=rep(0,ncol(d))
	for(i in 1:ncol(d)){
		share[i]=length(which(d[,i]>cutoff))/nrow(d)
	}
	return(median(share))
}

cutoff=0.005
T.all=c(1,2,4,6,8)
l.test=2;sig=2
dir=paste("./l=",l.test,"/sig=",sig,sep="")
setwd(dir)
out=matrix(0,nrow=length(Ne.all),ncol=length(T.all))
for(i in 1:length(T.all)){
	T=T.all[i]
	dir=paste("./T=",T,sep="")
	setwd(dir)
	for(j in 1:length(Ne.all)){
		Ne=Ne.all[j]
		fn=paste("lv_out_",n1,"_",log10(Ne),".txt",sep="")
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

d=data.frame(out.new)
d$log10Ne=factor(d$log10Ne,levels=sort(log10(Ne.all)),ordered=TRUE)
g=ggplot(d,aes(x=T,y=share,colour=log10Ne))
g=g+geom_point()+geom_line(linewidth=2,alpha=0.8)#+scale_color_manual(values=c("red","blue"))
g=g+theme_classic()
g=g+xlab("Time")+ylab("Sharing")+scale_x_continuous(breaks=seq(0,1e8,2e7))
g=g+theme(axis.text=element_text(size=12),axis.title=element_text(size=15),legend.text=element_text(size=15),legend.title=element_text(size=15),legend.key.size=unit(0.7,'cm')) # Adjust font size of axis labels, axis text, and legend

fn.out=paste("plot_share_l=",l.test,"_sig=",sig,".pdf",sep="")
setwd("..");setwd("..")
ggsave(g,file=fn.out,width=6,height=4)


