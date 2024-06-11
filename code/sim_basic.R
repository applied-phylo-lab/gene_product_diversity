# Basic setting (editing-type modification)
# 2 isoforms (unmidified I_0 is functional, modified I_1 is toxic), all loci have equal effects

library(expm)
library(ggplot2)

# Calculate normalized cis-genotypic value
# Input: genotype (a vector of 0 and 1), effect sizes of effector allele at all loci (a vector)
cis.gv <- function(gt,effect){
	v=sum(gt*effect)
	vmax=sum(effect)
	return(v/vmax)
}

# Calculate modification rate beta
# beta as a linear function of cis-genotypic value
# Input: normalized cis-genotypic value (a number, calculated using cis.gv()), trans- genotypic value, shape parameter characterizing relationship between cis- loci and enzyme-substrate affinity (a number), background/non-specific modification activity (a number)
beta.calc <- function(cis,trans,C,epsilon=0){
	beta=trans*(C*cis+epsilon)
	return(beta)
}

# Genotype-phenotype mapping (2 isoforms)
# Input: genotypic value of overall expression level, rate of modification (calculated using beta.calc()), decay rates (numbers)
g2p <- function(alpha,beta,gamma_0,gamma_1){
	P_0=alpha/(beta+gamma_0)
	P_1=alpha*beta/(gamma_1*(beta+gamma_0))
	return(c(P_0,P_1))
}

# Calculate fitness w/ respect to an isoform under stabilizing selection
# Input: isoform abundance (a number), optimum (a number), SD of fitness function
# Input values are not log-transformed and should be log-transformed by the function
fitness.stab <- function(z,opt,sig){
	d=log(z/opt) # Convert to log scale
	if(sig>0){
		w=dnorm(d,mean=0,sd=sig)/dnorm(0,mean=0,sd=sig)
	}else{
		w=1
	}
	return(w)
}

# Calculate fitness w/ respect to a deleterious isoform
# Input: isoform abundance, shape parameter (a positive number)
fitness.del <- function(z,par){
	w=exp(-par*z)
	return(w)
}

# Overall fitness as a product of fitness w/ respect to each isoform
# Input: isoform abundances (a vector), selection parameter (a list; vector-->stabilizing, number-->deleterious, NULL=neutral/equivalent to I_0), whether each modified isoform is functionally equivalent to I_0 (a vector, 1 means equivalent to I_0)
fitness.prod <- function(z,selection,eq=NULL){
	w=rep(0,length(z))
	type=lengths(selection) # Type of selection, which corresponds to number of selection parameters
	if(length(eq)>0){
		z[1]=z[1]+sum(z[which(eq==1)])
	}
	for(i in 1:length(z)){
		if(type[i]==2){
			w[i]=fitness.stab(z[i],selection[[i]][1],selection[[i]][2])
		}else{
			if(type[i]==1){
				w[i]=fitness.del(z[i],selection[[i]])
			}else{
				w[i]=1
			}
		}
	}
	return(prod(w))
}

# Fixation probability given ancestral and mutant fitness
fix.prob <- function(wa,wm,Ne){
	if(wa>0){
		s=(wm/wa)-1 # Coefficient of selection
		if(s==0){
			p=1/(2*Ne)
		}else{
			p=(1-exp(-2*s))/(1-exp(-4*Ne*s))
		}
	}else{ # The ancestral fitness is close to 0 (close enough to be recognized as zero by R)
		if(wm==0){ # Both ancestral and mutant fitness are close to 0
			p=1/(2*Ne) # The mutation is considered neutral
		}else{ # Ancestral fitness is close to 0 while mutant fitness isn't
			p=1 # The mutation is considered strongly beneficial
		}
	}
	return(p)
}

# Obtain transition matrix given other parameters
# Simplest version: 2 isoforms, constant transition matrix (no trans- substitutions, no environmental change), all cis- loci have equal effect size (i.e., no. of states = no. of loci+1) and mutational spectrum
# Input: number of cis- loci, mutation rates (a vector, elements being u01 and u10, respectively; 2*Ne*u should not exceed 1), selection parameter (a list), effective population size, other parameters (alpha,trans,C,gamma_0,gamma_1,epsilon,eq)
tm <- function(nloci,u,selection,Ne,alpha,trans,C,epsilon=0,gamma_0,gamma_1,eq=NULL){
	mat=matrix(0,nrow=nloci+1,ncol=nloci+1)
	for(i in 1:nrow(mat)){
		v1=(i-1)/nloci
		beta1=beta.calc(v1,trans,C,epsilon)
		z1=g2p(alpha,beta1,gamma_0,gamma_1)
		w1=fitness.prod(z1,selection,eq)
		for(j in 1:ncol(mat)){
			if(abs(i-j)==1){ # Consider adjescent states only
				v2=(j-1)/nloci
				beta2=beta.calc(v2,trans,C,epsilon)
				z2=g2p(alpha,beta2,gamma_0,gamma_1)
				w2=fitness.prod(z2,selection,eq)
				fp=fix.prob(w1,w2,Ne)
				if(i<j){
					uij=2*Ne*(nloci-(i-1))*u[1]
				}else{
					uij=2*Ne*(i-1)*u[2]
				}
				mat[i,j]=uij*fp
			}
		}
	}
	for(i in 1:nrow(mat)){
		mat[i,i]=1-sum(mat[i,])
	}
	return(mat)
}


Ne.all=c(1e2,1e3,1e4,1e5) # Effective population sizes to consider
alpha.all=exp(0:5) # Optimum to be set as alpha/gamma_0
nloci.all=1:10 # Number of cis-loci 
u.all=1e-9*rbind(c(1,1),c(1.5,0.5),c(0.5,1.5)) # Per-locus mutation rates
trans.all=(1:10)/5 # Trans-genotypic value
C=1 # Scaling between cis-genotypic value and regulatory efffect
gamma_0=1;gamma_1=1 # Decay rates
T=1e8 # Time (may not be enough to reach stationary distribution, final dependent on initial state)
sig=10 # Width of fitness function for P_0
lambda=1e-3 # Strength of selection on P_1

# Get a data matrix containing parameter combinations to consider
comb.all=rep(0,6)
for(Ne in Ne.all){
	for(alpha in alpha.all){
		for(nloci in nloci.all){
			for(trans in trans.all){
				for(i in 1:nrow(u.all)){
					row=c(Ne,alpha,nloci,trans,u.all[i,])
					comb.all=rbind(comb.all,row)
				}
			}
		}
	}
}
comb.all=comb.all[2:nrow(comb.all),]
rownames(comb.all)=NULL

# Get expected mean phenotype when I_0 is functional and I_1 is deleterious
# Output data matrix (columns: mean normalized cis- genotypic value, isoform abundances and modification frequency based on the mean gv)
out=matrix(0,nrow=nrow(comb.all),ncol=4)
for(i in 1:nrow(comb.all)){
	# Read parameter values of the row
	Ne=comb.all[i,1]
	alpha=comb.all[i,2]
	nloci=comb.all[i,3]
	trans=comb.all[i,4]
	u=comb.all[i,5:6]
	opt=alpha/gamma_0
	par=list(c(opt,sig),lambda)
	mat=tm(nloci,u,par,Ne,alpha,trans,C,gamma_0,gamma_1) # Obtain transition matrix
	start=rep(0,nloci+1);start[1]=1 # Start from zero cis-genotypic value
	distr=start%*%(mat %^% T) # Distribution of cis-genotypic value
	v_mean=sum((distr*(0:nloci))/nloci) # Mean cis-genotypic value
	beta_mean=beta.calc(v_mean,trans,C) # Calculate the corresponding beta
	phe=g2p(alpha,beta_mean,gamma_0,gamma_1) # Calculate the corresponding isoform abundances
	out[i,1]=v_mean;out[i,2:3]=phe;out[i,4]=phe[2]/sum(phe) # Data to write: mean cis-genotypic value, isoform abundances, modification level
	if(is.complex(out[i,1])==TRUE|is.complex(out[i,2])==TRUE|is.complex(out[i,3])==TRUE|is.complex(out[i,4])==TRUE){ # If complex number is produced, there is something wrong
		break
	}
}

# Write data file
out=data.frame(comb.all,out)
colnames(out)=c("Ne","alpha","l","Q","u01","u10","v","P_0","P_1","frac")
out$alpha=log(out$alpha)
write.table(out,file="out_basic.txt",sep="\t")

# Plot a subset
d<-read.table("out_basic.txt",sep="\t",header=TRUE)

# Check effect of Ne, nloci, and mutational bias
sub=which((d$alpha==0)&(d$Q==1)&(d$u01==d$u10)) # Rows with mutational bias
#sub=which((d$alpha==0)&(d$Q==1)&(d$u01>d$u10)) # Rows with mutational bias towards the effector allele
#sub=which((d$alpha==0)&(d$Q==1)&(d$u01<d$u10)) # Rows with mutational bias towards the null allele
dsub=d[sub,]
dsub$l=factor(dsub$l,levels=sort(unique(dsub$l)),ordered=TRUE)
g<-ggplot(dsub,aes(x=log10(Ne),y=frac,colour=l))
g=g+geom_point()+geom_line()+scale_color_brewer(palette="Paired")
g=g+theme_classic()
g=g+xlab(expression(paste("Lo",g[10],N[e])))+ylab("Modification level")
lt=expression(paste("Number of ",italic(cis),"- loci"))
g=g+labs(color=lt)
g=g+ylim(c(0,0.15))
g=g+theme(axis.text=element_text(size=12),axis.title=element_text(size=15),legend.text=element_text(size=15),legend.title=element_text(size=15),legend.key.size=unit(0.7,'cm')) # Adjust font size of axis labels, axis text, and legend
ggsave("plot_basic.pdf",plot=g,width=6,height=5)
#ggsave("plot_basic_mb1.pdf",plot=g,width=6,height=5)
#ggsave("plot_basic_mb2.pdf",plot=g,width=6,height=5)

# Check effect of expression level
l.test=2 # Choose a subset of rows with the same number of cis-loci 
sub=which((d$l==l.test)&(d$Q==1)&(d$u01==d$u10))
dsub=d[sub,]
dsub$alpha=dsub$alpha-log(gamma_0) # convert to optimum of P_0
dsub$alpha=factor(dsub$alpha,levels=sort(unique(dsub$alpha)),ordered=TRUE)
g<-ggplot(dsub,aes(x=log10(Ne),y=frac,colour=alpha))
g=g+geom_point()+geom_line(linewidth=1.5)+scale_color_brewer(palette="Paired")
g=g+theme_classic()
g=g+xlab(expression(paste("Lo",g[10],N[e])))+ylab("Modification level")
lt=expression(paste(italic(ln),P["0.opt"]))
g=g+labs(color=lt)
g=g+ylim(c(0,0.1))
g=g+theme(axis.text=element_text(size=12),axis.title=element_text(size=15),legend.text=element_text(size=15),legend.title=element_text(size=15),legend.key.size=unit(0.7,'cm')) # Adjust font size of axis labels, axis text, and legend
fn=paste("plot_basic_l",l.test,".pdf",sep="")
ggsave(fn,plot=g,width=6,height=5)

