# Simulate cis-trans coevolution
# Selective advantage of trans- factor mediated by something other than focal modifications 

library(expm)
library(ggplot2)

setwd("/Users/daohanji/Desktop/gene_product_diversity/out/")

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

# Simulate cis-trans coevolution
# Two isoforms per gene, one function and one deleterious
sim.evo.2<-function(type.all,nloci.all,alpha.all,gamma_0.all,gamma_1.all,selection.all,T,Ne,u.cis,u.trans,C,epsilon,start){
	nmut=rpois(1,lambda=u.trans[1]*2*Ne*T)
	if(nmut>=T){
		return("error")
	}else{
		ngene=length(type.all)
		tm.all=list()
		distr.all=list()
		v.all=start[[1]]
		trans=start[[2]]
		beta.all=rep(0,ngene)
		z.all=matrix(0,nrow=ngene,ncol=2)
		for(i in 1:ngene){
			tm.all[[i]]=tm(nloci.all[i],u.cis,selection.all[[i]],Ne,alpha.all[i],trans,C,epsilon,gamma_0.all[i],gamma_1.all[i])
			distr.all[[i]]=rep(0,nloci.all[i]+1);distr.all[[i]][v.all[i]+1]=1
			beta.all[i]=beta.calc(v.all[i]/nloci.all[i],trans,C,epsilon)
			z.all[i,]=g2p(alpha.all[i],beta.all[i],gamma_0.all[i],gamma_1.all[i])
		}
		tmut=sort(sample(1:T,nmut))
		t.last=0
		for(j in 1:nmut){
			t=tmut[j]
			w0=1
			for(i in 1:ngene){
				distr.all[[i]]=distr.all[[i]]%*%(tm.all[[i]] %^% (t-t.last)) # Calculate distribution at the time point under consideration
				v.all[i]=sample(0:nloci.all[i],1,prob=distr.all[[i]])
				w0=w0*fitness.prod(z.all[[i]],selection.all[[i]])
			}
			w0.trans=fitness.stab(trans,opt.trans,sig.trans)
			w0=w0*w0.trans
			trans.mut=exp(log(trans)+rnorm(1,mean=0,sd=u.trans[2]))
			w1=1
			beta.mut=rep(0,ngene)
			z.mut=matrix(0,nrow=ngene,ncol=2)
			for(i in 1:ngene){
				beta.mut[i]=beta.calc(v.all[i]/nloci.all[i],trans.mut,C,epsilon)
				z.mut[i,]=g2p(alpha.all[i],beta.mut[i],gamma_0.all[i],gamma_1.all[i])
				w1=w1*fitness.prod(z.mut[[i]],selection.all[[i]])
			}
			w1.trans=fitness.stab(trans.mut,opt.trans,sig.trans)
			w1=w1*w1.trans
			fp=fix.prob(w0,w1,Ne)
			if.fix=rbinom(n=1,size=1,prob=fp)
			if(if.fix==1){
				trans=trans.mut
				beta.all=beta.mut
				z.all=z.mut
				for(i in 1:ngene){
					tm.all[[i]]=tm(nloci.all[i],u.cis,selection.all[[i]],Ne,alpha.all[i],trans.mut,C,epsilon,gamma_0.all[i],gamma_1.all[i])
				}
			}
			t.last=t
		}
		return(list(z.all,trans))
	}
}

T=1e8 # Number of time steps to run for each gene/modification event
Ne=1e3
u.cis=c(1e-9,1e-9) # Mutation rate at cis- loci, elements representing mutation to effector ellele and to null allele, respectively
u.trans=c(1e-8,1e-1) # Rate and SD of effect (log scale) of mutations affecting trans genotypic value
C=1 # Affinity parameter, assumed as constant across genes
epsilon=1e-3 # Non-specific modification intensity, assumed as constant across genes
# Selection on Q independent of focal modifications
opt.trans=1;sig.trans=20

# Assign gene-specific parameters
n1=100;n2=0;ngene=n1+n2
type.all=c(rep(1,n1),rep(2,n2))
nloci.all=sample(1:10,ngene,replace=TRUE)
alpha.all=exp(rnorm(ngene,mean=0,sd=1));alpha.all=sort(alpha.all)
gamma_0.all=rep(1,ngene) # Decay rate of I_0, distribution type TBD
gamma_1.all=rep(1,ngene) # Decay rate of I_1, distribution type TBD
sig.all=1/(10^(rnorm(ngene,mean=-1,sd=0.1))) # Strength of stabilizing selection on expression level
lambda.all=rep(1e-3,ngene)
selection.all=list()
for(i in 1:ngene){
	if(type.all[i]==1){ # Deleterious
		selection.all[[i]]=list(c(alpha.all[i]/gamma_0.all[i],sig.all[i]),lambda.all[i])
	}else{ # Required modification
		selection.all[[i]]=list(lambda.all[i],c(alpha.all[i]/gamma_1.all[i],sig.all[i]))
	}
}

# Null distribution of cis- genotypic value
pgv=list()
for(l in 1:10){
	pgv[[l]]=rep(0,l+1)
	for(i in 0:l){
		pgv[[l]][i+1]=choose(l,i)
	}
	pgv[[l]]=pgv[[l]]/(sum(pgv[[l]]))
}

# Starting cis- and trans- genotypic values
start=list(rep(0,ngene),1)
for(i in 1:ngene){
	if(type.all[i]==1){
		start[[1]][i]=sample(0:(nloci.all[i]),1,prob=pgv[[nloci.all[i]]])
	}else{
		start[[1]][i]=nloci.all[i]
	}
}

# Test run
x=sim.evo.2(type.all,nloci.all,alpha.all,gamma_0.all,gamma_1.all,selection.all,T,Ne,u.cis,u.trans,C,epsilon,start)
x[[2]]
hist(x[[1]][,2]/(x[[1]][,2]+x[[1]][,1]),breaks=20)
# End-point fitness
w=1
for(i in 1:ngene){
	w=w*fitness.prod(x[[1]][i,],selection.all[[i]])
}
w=w*fitness.stab(x[[2]],opt.trans,sig.trans)
w

# The same set of gene, varying Ne, toxicity, and selection on Q
Ne.all=c(1e2,1e3,1e4,1e5)
lambda.scale.all=c(.1,1,10)
opt.trans.all=exp(c(-2,-1,0,1,2))
par.all=rep(0,3)
for(i in 1:length(Ne.all)){
	for(j in 1:length(lambda.scale.all)){
		for(k in 1:length(opt.trans.all)){
			par.all=rbind(par.all,c(Ne.all[i],lambda.scale.all[j],opt.trans.all[k]))
		}
	}
}
par.all=par.all[2:nrow(par.all),]

Nrep=100
lv.out=matrix(0,nrow=nrow(par.all),ncol=Nrep) # Median modification level
trans.out=matrix(0,nrow=nrow(par.all),ncol=Nrep) # End-point Q
fitness.out=matrix(0,nrow=nrow(par.all),ncol=Nrep) # End-point finess
for(c in 1:nrow(par.all)){
	Ne=par.all[c,1]
	lambda.all=rep(1e-3,ngene)*par.all[c,2]
	opt.trans=par.all[c,3]
	for(n in 1:Nrep){
		x=sim.evo.2(type.all,nloci.all,alpha.all,gamma_0.all,gamma_1.all,selection.all,T,Ne,u.cis,u.trans,C,epsilon,start)
		lv=x[[1]][,2]/(x[[1]][,2]+x[[1]][,1])
		lv.out[c,n]=median(lv)
		trans.out[c,n]=x[[2]]
		w=1
		for(i in 1:ngene){
			w=w*fitness.prod(x[[1]][i,],selection.all[[i]])
		}
		w=w*fitness.stab(x[[2]],opt.trans,sig.trans)
		fitness.out[c,n]=w
	}
}

