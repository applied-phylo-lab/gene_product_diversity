# Generate transcriptome-wide distributions
library(expm) # Required to do matrix exponential (in A %^% T form)
library(ggplot2)
#library(MASS)

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

# Obtain transition matrix given other parameters
# Simplest version: 2 isoforms, constant transition matrix (no trans- substitutions, no environmental change), all cis- loci have equal effect size (i.e., no. of states = no. of loci+1) and mutational spectrum
# Input: number of cis- loci, mutation rates, mutation rates (a vector, elements being u01 and u10, respectively; 2*Ne*u should not exceed 1), selection parameter (a list), effective population size, "background" parameters (alpha,trans,C,gamma_0,gamma_1)
tm <- function(nloci,u,par,Ne,alpha,trans,C,epsilon=0,gamma_0,gamma_1,eq=NULL){
	mat=matrix(0,nrow=nloci+1,ncol=nloci+1)
	for(i in 1:nrow(mat)){
		v1=(i-1)/nloci
		beta1=beta.calc(v1,trans,C,epsilon)
		z1=g2p(alpha,beta1,gamma_0,gamma_1)
		w1=fitness.prod(z1,par,eq)
		for(j in 1:ncol(mat)){
			if(abs(i-j)==1){ # Consider adjescent states only
				v2=(j-1)/nloci
				beta2=beta.calc(v2,trans,C,epsilon)
				z2=g2p(alpha,beta2,gamma_0,gamma_1)
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
Ne=1e3
u=c(1e-9,1e-9)
trans=1;epsilon=1e-3

n1=2e4 # I_0 is functional
n2=0 # I_1 is functional
n3=0 # I_0 and I_1 are functionally equivalent
ngene=n1+n2+n3
type.all=c(rep(1,n1),rep(2,n2),rep(3,n3))
alpha.all=exp(rnorm(ngene,mean=0,sd=1)) # Expression rate, sampled from log-normal distribution
nloci.all=sample(1:10,ngene,replace=TRUE) # Number of cis- loci, sampled uniformly from a pre-chosen set
C.all=exp(rnorm(ngene,mean=0,sd=1)) # Parameter characterizing relationship between cis- loci and enzyme-substrate affinity, sampled from log-normal distribution
gamma_0.all=rep(1,ngene) # Decay rate of I_0, distribution type TBD
gamma_1.all=rep(1,ngene) # Decay rate of I_1, distribution type TBD
sig.all=1/(10^(rnorm(ngene,mean=-1,sd=0.1))) # Width of Gaussian fitness function (for P_0 or P_1); sample "importance" first, then take reciprocal to be width
lambda.all=rep(1e-3,ngene) # Strength of selection on P_1; distribution type TBD

comb.all=data.frame(alpha.all,nloci.all,C.all,gamma_0.all,gamma_1.all,type.all,sig.all,lambda.all)

out=matrix(0,nrow=nrow(comb.all),ncol=4)
for(i in 1:ngene){
	alpha=comb.all[i,1]
	nloci=comb.all[i,2]
	C=comb.all[i,3]
	gamma_0=comb.all[i,4]
	gamma_1=comb.all[i,5]
	type=comb.all[i,6]
	if(type==3){ # Modification is neutral
		eq=c(0,1)
		par=list(c(alpha/gamma_0,sig),NULL)
	}else{ # Modification is non-neutral
		eq=c(0,0)
		sig=comb.all[i,7]
		lambda=comb.all[i,8]
		if(type==1){ # I_0 is functional
			par=list(c(alpha/gamma_0,sig),lambda)
		}else{ # I_1 is functional
			par=list(lambda,c(alpha/gamma_1,sig))
		}
	}
	mat=tm(nloci,u,par,Ne,alpha,trans,C,epsilon,gamma_0,gamma_1,eq)
	start=rep(0,nloci+1);start[1]=1
	distr=start%*%(mat %^% T)
	v_sample=sample(0:nloci,1,prob=distr)/nloci # Sample a cis- genotypic value from the distribution
	beta_sample=beta.calc(v_sample,trans,C,epsilon)
	phe=g2p(alpha,beta_sample,gamma_0,gamma_1) # Get a final phenotype
	out[i,1]=v_sample;out[i,2:3]=phe;out[i,4]=phe[2]/sum(phe)
}

d=data.frame(comb.all,out)

# Impose a detectability cutoff and get a subset of "detected" sites/genes
id.cutoff=1e-3
d.id=d[which(d[,1]*d[,12]>id.cutoff),]

