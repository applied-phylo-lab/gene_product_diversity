# Simulate evolution of isoform abundances along a branch
# Parameters other than cis- genotypic value can also evolve

library(expm)
library(ggplot2)

# Calculate normalized cis-genotypic value
# Input: genotype (a vector of 0 and 1), effect sizes of effector allele at all loci (a vector)
cis.gv <- function(gt,effect=1){
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

# Function to simulate evolution under selection on isoform abundances (2 isoforms, cis- loci have equal effect)
# For one gene only, trans- factor assumed unchangeable
# Parameters: values of constant parameters (a vector), starting state (a vector), mutation parameters (a list), selection parameters (a list), Ne, time
# Constants: nloci, trans, C, epsilon, gamma_0, gamma_1
# Starting state: cis- genotypic value (a number), alpha (a number)
# Mutation parameters: cis- loci (a vector containing u01 and u10), alpha (a vector containing rate and effect SD)
# Selection parameters in the same format as used for fitness.prod
sim.evo <- function(constant,start,U,selection,Ne,T){
	nmut=rpois(1,lambda=U[[2]][1]*2*Ne*T) # Number of mutations that would occur (should be <T, otherwise the model doesn't work)
	if(nmut>=T){
		return("error")
	}else{
		if(nmut>0){
			nloci=constant[[1]];trans=constant[[2]];C=constant[[3]];epsilon=constant[[4]];gamma_0=constant[[5]];gamma_1=constant[[6]]
			tmut=sort(sample(1:T,nmut)) # Time points at which mutations happen
			v=start[1];alpha=start[2] # Set initial values of cis- gv and alpha
			distr=rep(0,nloci+1);distr[v+1]=1 # Set initial distribution of cis- gv
			v.all=rep(0,nmut);alpha.all=rep(0,nmut) # vectors storing all values of cis- gv an alpha
			mat=tm(nloci,U[[1]],selection,Ne,alpha,trans,C,epsilon,gamma_0,gamma_1) # Set initial transition matrix
			t.last=0
			for(i in 1:nmut){
				t=tmut[i]
				distr=distr%*%(mat %^% (t-t.last)) # Calculate distribution at the time point under consideration
				v=sample(0:nloci,1,prob=distr);v.all[i]=v # Sample a value from the distribution
				distr=rep(0,nloci+1);distr[v+1]=1 # Reset the distribution
				beta=beta.calc(v/nloci,trans,C,epsilon) # Calculate the corresponding beta
				z=g2p(alpha,beta,gamma_0,gamma_1) # Calculate the corresponding phenotype
				alpha.mut=exp(log(alpha)+rnorm(1,mean=0,sd=U[[2]][2]))
				z.mut=g2p(alpha.mut,beta,gamma_0,gamma_1)
				fp=fix.prob(fitness.prod(z,selection),fitness.prod(z.mut,selection),Ne)
				if.fix=rbinom(n=1,size=1,prob=fp)
				if(if.fix==1){
					alpha=alpha.mut
					mat=tm(nloci,U[[1]],selection,Ne,alpha,trans,C,epsilon,gamma_0,gamma_1) # Re-calculate transition matrix
				}
				alpha.all[i]=alpha
			}
			return(list(v.all,alpha.all))
		}else{
			return("no.mut")
		}
	}
}

Nrep=100
T=1e8
Ne.all=c(1e2,1e3,1e4,1e5)
alpha.all=exp(0:5)
nloci.all=1:10
comb.all=rep(0,3)
for(Ne in Ne.all){
	for(alpha in alpha.all){
		for(nloci in nloci.all){
			row=c(Ne,alpha,nloci)
			comb.all=rbind(comb.all,row)
		}
	}
}
comb.all=comb.all[2:nrow(comb.all),]
rownames(comb.all)=NULL

U=list(c(1e-9,1e-9),c(1e-8,0.1))
trans=1;C=1;epsilon=1e-2;gamma_0=1;gamma_1=1
sig=10 # Width of fitness function for P_0
lambda=1e-3 # Strength of selection on P_1

v.out=matrix(0,nrow=nrow(comb.all),ncol=Nrep);alpha.out=matrix(0,nrow=nrow(comb.all),ncol=Nrep)
for(c in 1:nrow(comb.all)){
	Ne=comb.all[c,1]
	start=c(0,comb.all[c,2])
	constant=c(comb.all[c,3],trans,C,epsilon,gamma_0,gamma_1)
	selection=list(c(comb.all[c,2]/gamma_0,sig),lambda)
	for(n in 1:Nrep){
		out=sim.evo(constant,start,U,selection,Ne,T)
		if(length(out)==2){
			v.out[c,n]=out[[1]][length(out[[1]])]
			alpha.out[c,n]=out[[1]][length(out[[2]])]
		}
	}
}




