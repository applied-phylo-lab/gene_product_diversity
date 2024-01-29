# 3 isoforms, I_0, I_1, and I_2; I_1 is functional while I_2 is deleterious
# Representing alternative splicing
# Transcriptome-wide distribution
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
fitness.prod <- function(z,par,eq=NULL){
	w=rep(0,length(z))
	type=lengths(par) # Type of selection, which corresponds to number of selection parameters
	if(length(eq)>0){
		z[1]=z[1]+sum(z[which(eq==1)])
	}
	for(i in 1:length(z)){
		if(type[i]==2){
			w[i]=fitness.stab(z[i],par[[i]][1],par[[i]][2])
		}else{
			if(type[i]==1){
				w[i]=fitness.del(z[i],par[[i]])
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
tm.type.2 <- function(nloci,u,par,Ne,alpha,trans,C,epsilon=c(0,0),gamma_0,gamma,eq=NULL){
	mat=matrix(0,nrow=nloci+1,ncol=nloci+1)
	for(i in 1:nrow(mat)){
		v1=(i-1)/nloci # cis- genotypic value for beta_1
		beta1=beta.calc.type.2(v1,trans,C,epsilon)
		z1=g2p(alpha,beta1,gamma_0,gamma) # Remember to use the generic version of g2p()
		w1=fitness.prod(z1,par,eq)
		for(j in 1:ncol(mat)){
			if(abs(i-j)==1){ # Consider adjescent states only
				v2=(j-1)/nloci
				beta2=beta.calc.type.2(v2,trans,C,epsilon)
				z2=g2p(alpha,beta2,gamma_0,gamma)
				w2=fitness.prod(z2,par,eq)
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

# Parameter values
# Gene-specific parameters sampled from pre-chosen distributions
T=1e8 # Number of time steps to run for each gene/modification event
u=c(1e-9,1e-9)
gamma_0=0;trans=5;epsilon=c(1e-3,1e-3);C=c(1,1)

ngene=2e4
alpha.all=exp(rnorm(ngene,mean=0,sd=1)) # Expression rate, sampled from log-normal distribution
nloci.all=as.integer(rgamma(ngene,shape=3,scale=5))+1 # Number of cis- loci, sampled from a gamma distribution (mean~15, peak near mean)
gamma_1.all=rep(1,ngene) # Decay rate of I_0, distribution type TBD
gamma_2.all=rep(20,ngene) # Decay rate of I_1, distribution type TBD
sig.all=1/(10^(rnorm(ngene,mean=-1,sd=0.1))) # Width of Gaussian fitness function (for P_0 or P_1); sample "importance" first, then take reciprocal to be width
lambda.all=rep(1e-3,ngene) # Strength of selection on P_1; distribution type TBD

Ne=1e3
comb.all=data.frame(alpha.all,nloci.all,C.all,gamma_1.all,gamma_2.all,sig.all,lambda.all)
out=matrix(0,nrow=nrow(comb.all),ncol=4)
for(i in 1:ngene){
	alpha=comb.all[i,1]
	nloci=comb.all[i,2]
	gamma_1=comb.all[i,4]
	gamma_2=comb.all[i,5]
	gamma=c(gamma_1,gamma_2)
	sig=comb.all[i,6]
	lambda=comb.all[i,7]
	par=list(NULL,c(alpha/gamma_1,sig),lambda)
	mat=tm.type.2(nloci,u,par,Ne,alpha,trans,C,epsilon,gamma_0,gamma)
	start=rep(0,nloci+1);start[nloci+1]=1
	distr=start%*%(mat %^% T)
	v_sample=sample(0:nloci,1,prob=distr)/nloci # Sample a cis- genotypic value from the distribution
	beta_sample=beta.calc.type.2(v_sample,trans,C,epsilon)
	phe=g2p(alpha,beta_sample,gamma_0,gamma) # Get a final phenotype
	out[i,1]=v_sample;out[i,2:3]=phe[2:3];out[i,4]=phe[3]/(phe[2]+phe[3])
}
colnames(out)=c("gv","P_1","P_2","lv")
d=data.frame(comb.all,out)
write.table(d,file="distr_mix_as.txt",sep="\t")

id.cutoff=0 # Detectability cutoff
dsub=d[which(d[,1]*d[,11]>id.cutoff),]
g=ggplot(dsub,aes(x=lv))+geom_histogram(binwidth=0.01)
g=g+theme_classic()
g=g+xlab("AS level")+ylab("Gene number")
g=g+theme(axis.text=element_text(size=12),axis.title=element_text(size=15))
fn.out=paste("hist_mix_as_",id.cutoff,"_",Ne,".pdf",sep="")
ggsave(fn.out,plot=g,width=5.5,height=5)




