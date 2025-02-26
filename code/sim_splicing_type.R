# Generate model predictions regarding splicing-type modification
# 3 isoforms, I_0, I_1, and I_2; I_1 is functional while I_2 is deleterious

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
# Two effector alleles, one increasing beta_1, the other increasing beta_2 (e.g., alternative splicing)
# Input: normalized cis-genotypic value for beta_1 (that for beta_2 is 1-cis), trans- genotypic value, shape parameter (a vector), background/non-specific modification activity (a vector; epsilon[1] should be high enough to ensure low mis-splicing rate)
beta.calc.type.2 <- function(cis,trans,C,epsilon=c(0,0)){
	beta1=trans*(C[1]*cis+epsilon[1])
	beta2=trans*(C[2]*(1-cis)+epsilon[2])
	return(c(beta1,beta2))
}

# Genotype-phenotype mapping (>2 isoforms)
# Input: genotypic value of overall expression rate, rate of modification (a vector, each element calculated using beta.calc()), decay rate of unmodified isoform (a number), decay rates of modified isoforms (a vector)
g2p <- function(alpha,beta,gamma_0,gamma){
	ni=length(beta) # Number of modified isoforms
	A=matrix(0,nrow=ni+1,ncol=ni+1)
	A[1,1]=sum(beta)+gamma_0
	A[2:nrow(A),1]=beta
	diag(A)[2:nrow(A)]=-gamma
	b=c(alpha,rep(0,ni))
	return(solve(A,b))
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

# Obtain transition matrix given other parameters (type-2 modification)
# Input: number of cis- loci, mutation rates (a vector, elements being u01 and u10, respectively; 2*Ne*u should not exceed 1), selection parameter (a list), effective population size, other parameters (alpha,trans,C,gamma_0,gamma_1,epsilon,eq)
tm.type.2 <- function(nloci,u,selection,Ne,alpha,trans,C,epsilon=c(0,0),gamma_0,gamma,eq=NULL){
	mat=matrix(0,nrow=nloci+1,ncol=nloci+1)
	for(i in 1:nrow(mat)){
		v1=(i-1)/nloci # cis- genotypic value for beta_1
		beta1=beta.calc.type.2(v1,trans,C,epsilon)
		z1=g2p(alpha,beta1,gamma_0,gamma) # Remember to use the generic version of g2p()
		w1=fitness.prod(z1,selection,eq)
		for(j in 1:ncol(mat)){
			if(abs(i-j)==1){ # Consider adjescent states only
				v2=(j-1)/nloci
				beta2=beta.calc.type.2(v2,trans,C,epsilon)
				z2=g2p(alpha,beta2,gamma_0,gamma)
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

# Parameter range
Ne.all=c(1e2,1e3,1e4,1e5) # Effective population sizes to consider
alpha.all=exp(0:5) # Rate of expression; optimum to be set as alpha/gamma_1
nloci.all=10*(1:5) # Number of cis-loci (should be greater than editing-type)
u.all=1e-9*rbind(c(1,1),c(1.5,0.5),c(0.5,1.5)) # Per-locus mutation rates (each vector contains u01 and u10, respectively)
trans=100 # trans-genotypic value characterizing modification activity (should be high enough such that P_0 is much smaller than P_0+P_1+P_2)
C=c(1,1) # Scaling between cis- genotypic value and regulatory efffect
epsilon=c(1,1)*1e-3 # Non-spcific interaction strength
gamma_0=0 # Decay rate of unmodified isoform (set to be 0 for this study)
gamma_1=1 # Decay rate of functional isoform
gamma_2.all=(1:5)*20 # Decay rate of mis-spliced isoform, reflecting efficiency of quality-control mechanisms like stop-mediated degradation
T=1e8 # Time (may not be enough to reach stationary distribution, final dependent on initial state)
sig=10 # Width of fitness function, characterizing strength of stabilizing on abundance of the functional isoform (P_1)
lambda=1e-3 # Strength of selection on abundance of the deleterious isoform (P_2)

# Get a data matrix containing parameter combinations to consider
comb.all=rep(0,4)
for(Ne in Ne.all){
	for(alpha in alpha.all){
		for(nloci in nloci.all){
			for(gamma_2 in gamma_2.all){
				for(i in 1:nrow(u.all)){				
					row=c(Ne,alpha,nloci,gamma_2,u.all[i,])
					comb.all=rbind(comb.all,row)
				}			
			}
		}
	}
}
comb.all=comb.all[2:nrow(comb.all),]
rownames(comb.all)=NULL

# Get expected mean phenotype when I_1 is functional and I_2 is deleterious
# Output data matrix (columns: mean normalized cis- genotypic value, modified isoform abundances corresponding to the mean cis-genotypic value, relative abundance of I_2 corresponding to the mean cis-genotypic value)
out=matrix(0,nrow=nrow(comb.all),ncol=4)
for(i in 1:nrow(comb.all)){
	# Read parameter values of the row
	Ne=comb.all[i,1]
	alpha=comb.all[i,2]
	nloci=comb.all[i,3]
	gamma_2=comb.all[i,4];gamma=c(gamma_1,gamma_2)
	opt=alpha/gamma_1 # Make optimal P_0 the value reached in the absence of mis-splicing
	par=list(NULL,c(opt,sig),lambda)
	mat=tm.type.2(nloci,u,par,Ne,alpha,trans,C,epsilon,gamma_0,gamma)
	start=rep(0,nloci+1);start[nloci+1]=1
	distr=start%*%(mat %^% T) # Stationary distribution of cis- genotypic value
	v_mean=sum((distr*(0:nloci))/nloci) # Mean cis- genotypic value
	if(v_mean>1){
		v_mean=1
	}
	beta_mean=beta.calc.type.2(v_mean,trans,C) # Calculate the corresponding beta
	phe=g2p(alpha,beta_mean,gamma_0,gamma) # Calculate the corresponding isoform abundances
	out[i,1]=v_mean;out[i,2:3]=phe[2:3];out[i,4]=phe[3]/(phe[2]+phe[3]) # Write results
	if(is.complex(out[i,1])==TRUE|is.complex(out[i,2])==TRUE|is.complex(out[i,3])==TRUE|is.complex(out[i,4])==TRUE){ # If complex number is produced, there is something wrong
		break
	}
}

# Write data file
out=data.frame(comb.all,out)
colnames(out)=c("Ne","alpha","l","gamma_2","v","P_1","P_2","frac")
out$alpha=log(out$alpha)
write.table(out,file="out_as.txt",sep="\t")

d<-read.table("out_as.txt",sep="\t")

# Check effect of Ne, nloci, and degradation efficiency on minor isoform fraction
alpha.test=1;gamma_2.test=100
sub=which((d$alpha==alpha.test)&(d$gamma_2==gamma_2.test)) # Choose a subset of data with a given combination alpha and gamma_2
dsub=d[sub,]
dsub$l=factor(dsub$l,levels=sort(unique(dsub$l)),ordered=TRUE) # Factorize the number of cis-loci for plotting
g<-ggplot(dsub,aes(x=log10(Ne),y=frac,colour=l))
g=g+geom_point()+geom_line()+scale_color_brewer(palette="Paired")
g=g+theme_classic()
g=g+xlab("")+ylab("")
lt=expression(paste("Number of ",italic(cis),"- loci"))
g=g+labs(color=lt)
g=g+ylim(c(0,0.005))
g=g+theme(axis.text=element_text(size=12),axis.title=element_text(size=15),legend.text=element_text(size=15),legend.title=element_text(size=15),legend.key.size=unit(0.7,'cm')) # Adjust font size of axis labels, axis text, and legend
fn=paste("plot_as_exp_a",alpha.test,"r",gamma_2.test,".pdf",sep="")
ggsave(fn,plot=g,width=6,height=5)

# Check expression level's effect on minor isoform fraction
l.test=50;gamma_2.test=60
sub=which((d$l==l.test)&(d$gamma_2==gamma_2.test)) # Choose a subset of data with a given combination cis-loci number and gamma_2
dsub=d[sub,];gamma_1=rep(1,nrow(dsub))
dsub$alpha=dsub$alpha-log(gamma_1) # convert to optimum of P_1
dsub$alpha=factor(dsub$alpha,levels=sort(unique(dsub$alpha)),ordered=TRUE)
g<-ggplot(dsub,aes(x=log10(Ne),y=frac,colour=alpha))
g=g+geom_point()+geom_line(linewidth=1.5)+scale_color_brewer(palette="Paired")
g=g+theme_classic()
g=g+xlab("")+ylab("")
lt=expression(paste(italic(ln),P["1.opt"]))
g=g+labs(color=lt)
g=g+ylim(c(0,0.005))
g=g+theme(axis.text=element_text(size=12),axis.title=element_text(size=15),legend.text=element_text(size=15),legend.title=element_text(size=15),legend.key.size=unit(0.7,'cm')) # Adjust font size of axis labels, axis text, and legend
fn=paste("plot_as_l",l.test,"r",gamma_2.test,".pdf",sep="")
ggsave(fn,plot=g,width=6,height=5)

# Plot normalized cis-genotypic value instead
alpha.test=1;gamma_2.test=100
sub=which((d$alpha==alpha.test)&(d$gamma_2==gamma_2.test)) # Choose a subset of data with a given combination alpha and gamma_2
dsub=d[sub,]
dsub=d[sub,]
dsub$l=factor(dsub$l,levels=sort(unique(dsub$l)),ordered=TRUE)
g<-ggplot(dsub,aes(x=log10(Ne),y=1-v,colour=l))
g=g+geom_point()+geom_line()+scale_color_brewer(palette="Paired")
g=g+theme_classic()
g=g+xlab("")+ylab("")
lt=expression(paste("Number of ",italic(cis),"- loci"))
g=g+labs(color=lt)
g=g+ylim(c(0,0.1))
g=g+theme(axis.text=element_text(size=12),axis.title=element_text(size=15),legend.text=element_text(size=15),legend.title=element_text(size=15),legend.key.size=unit(0.7,'cm')) # Adjust font size of axis labels, axis text, and legend
fn=paste("plot_as_v_exp_a",alpha.test,"r",gamma_2.test,".pdf",sep="")
ggsave(fn,plot=g,width=6,height=5)

l.test=50;gamma_2.test=60
sub=which((d$l==l.test)&(d$gamma_2==gamma_2.test)) # Choose a subset of data with a given combination cis-loci number and gamma_2
dsub=d[sub,];gamma_1=rep(1,nrow(dsub))
dsub$alpha=dsub$alpha-log(gamma_1) # convert to optimum of P_1
dsub$alpha=factor(dsub$alpha,levels=sort(unique(dsub$alpha)),ordered=TRUE)
g<-ggplot(dsub,aes(x=log10(Ne),y=1-v,colour=alpha))
g=g+geom_point()+geom_line(linewidth=1.5)+scale_color_brewer(palette="Paired")
g=g+theme_classic()
g=g+xlab("")+ylab("")
lt=expression(paste(italic(ln),P["1.opt"]))
g=g+labs(color=lt)
g=g+ylim(c(0,0.1))
g=g+theme(axis.text=element_text(size=12),axis.title=element_text(size=15),legend.text=element_text(size=15),legend.title=element_text(size=15),legend.key.size=unit(0.7,'cm')) # Adjust font size of axis labels, axis text, and legend
fn=paste("plot_as_v_l",l.test,"r",gamma_2.test,".pdf",sep="")
ggsave(fn,plot=g,width=6,height=5)


